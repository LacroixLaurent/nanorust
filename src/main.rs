mod bam;
mod bgzf;
mod budget;
mod coverage;
mod rds;
mod rstats;
mod signal;
mod x87;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser)]
#[command(version, about = "Fast nanoT BrdU parsing of dorado mod-call BAMs")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Parse a dorado BAM into <prefix>_nanoT_alldata.rds and a bamCoverage-like <prefix>.bw
    Run(RunArgs),
    /// Compare two bigWig files base by base
    BwCompare { a: PathBuf, b: PathBuf },
}

#[derive(Parser)]
struct RunArgs {
    /// Input BAM (dorado basecaller output, unsorted is fine; "-" for stdin)
    #[arg(short, long)]
    bam: PathBuf,
    /// Output prefix: writes <prefix>_nanoT_alldata.rds and <prefix>.bw
    #[arg(short, long)]
    out_prefix: Option<String>,
    /// Output alldata RDS (overrides the prefix-derived name)
    #[arg(long)]
    rds: Option<PathBuf>,
    /// Output coverage bigWig (overrides the prefix-derived name)
    #[arg(long)]
    bw: Option<PathBuf>,
    /// Total CPU cores to use (default: all cores). Decompression, processing and
    /// output compression share this budget, e.g. set it to HTCondor/SLURM's allocated CPUs
    #[arg(short, long)]
    threads: Option<usize>,
    /// Modification to extract, as in the MM tag: <base><strand><code>, e.g. T+B (BrdU),
    /// T+e (EdU), C+m (5mC), C+h (5hmC), A+a (6mA), or a ChEBI code like C+17802
    #[arg(long = "mod", default_value = "T+B")]
    modification: String,
    /// Binarise probabilities for signalbin: prob < THR -> 0, else 1 (R's binarise/bin_thr).
    /// med_signal is still computed on raw probabilities, as in R
    #[arg(long, value_name = "THR")]
    binarise: Option<f64>,
    /// Drop supplementary mappings (default keeps them like parsing_DoradoRemora_v18_Br.r)
    #[arg(long)]
    no_supplementary: bool,
    /// supp_filter max distance to the primary mapping
    #[arg(long, default_value_t = 15000)]
    max_dist: i64,
    /// Mappings need (end - start) > min_len
    #[arg(long, default_value_t = 1)]
    min_len: i64,
    /// Bin size for signalbin
    #[arg(long, default_value_t = 1000)]
    bin_size: i64,
    /// Bin size of the coverage bigWig (bamCoverage --binSize)
    #[arg(long, default_value_t = 50)]
    cov_bin_size: u32,
    /// Chromosome name prefixes excluded from the signal (not from coverage)
    #[arg(long, default_value = "chrM", value_delimiter = ',')]
    exclude_prefix: Vec<String>,
    /// gzip level of the RDS (R's saveRDS uses 6)
    #[arg(long, default_value_t = 6)]
    rds_level: u32,
    /// Floating-point flavour of R's mean(): "x87" reproduces R on x86-64 Linux
    /// (80-bit long double), "f64" reproduces R on arm64 macOS
    #[arg(long, value_enum, default_value_t = FloatMode::X87)]
    float_mode: FloatMode,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum FloatMode {
    X87,
    F64,
}

fn main() -> Result<()> {
    match Cli::parse().cmd {
        Cmd::Run(a) => run(a),
        Cmd::BwCompare { a, b } => bw_compare(&a, &b),
    }
}

struct Batch {
    data: Vec<u8>,
    ends: Vec<usize>,
}

const BATCH_BYTES: usize = 8 << 20;

fn run(mut a: RunArgs) -> Result<()> {
    let t0 = Instant::now();
    rstats::set_x87(matches!(a.float_mode, FloatMode::X87));
    if let Some(p) = &a.out_prefix {
        a.rds.get_or_insert_with(|| format!("{p}_nanoT_alldata.rds").into());
        a.bw.get_or_insert_with(|| format!("{p}.bw").into());
    }
    if a.rds.is_none() && a.bw.is_none() {
        bail!("nothing to do: give --out-prefix, --rds and/or --bw");
    }
    let threads = a.threads.unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get())).max(1);
    let input: Box<dyn std::io::Read + Send> = if a.bam.as_os_str() == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(File::open(&a.bam).with_context(|| format!("opening {}", a.bam.display()))?)
    };
    let input = std::io::BufReader::with_capacity(1 << 20, input);
    let params = signal::Params {
        min_len: a.min_len,
        bin_size: a.bin_size,
        keep_supplementary: !a.no_supplementary,
        exclude_prefixes: a.exclude_prefix.clone(),
        modspec: signal::ModSpec::parse(&a.modification)?,
        bin_values: signal::Params::bin_values(a.binarise),
    };
    let cov_bin = a.cov_bin_size;

    let (chrom_names, cov, maps, counters) = if threads == 1 {
        // everything inline in this thread: exactly one core
        let mut reader = bgzf::InlineBgzfReader::new(input);
        let (chrom_names, cov) = read_header(&mut reader, cov_bin)?;
        let mut w = Worker::new(&cov, &params, &chrom_names);
        let mut buf = Vec::new();
        loop {
            buf.clear();
            if !bam::read_raw_record(&mut reader, &mut buf)? {
                break;
            }
            w.record(&buf)?;
        }
        let (maps, cnt) = (w.out, w.cnt);
        (chrom_names, cov, maps, cnt)
    } else {
        // -t 2: one thread reads+inflates inline, one processes records.
        // -t N>2: compute stages (inflate, process) share N-1 permits; the
        // remaining core goes to the light block-reader and record-splitter threads.
        let compute = threads - 1;
        let budget = budget::Budget::new(compute);
        let mut reader: Box<dyn std::io::Read + Send> = if threads == 2 {
            Box::new(bgzf::InlineBgzfReader::new(input))
        } else {
            Box::new(bgzf::ParallelBgzfReader::new(input, compute, budget.clone()))
        };
        let (chrom_names, cov) = read_header(&mut reader, cov_bin)?;
        let (tx, rx) = crossbeam_channel::bounded::<Batch>(compute * 2);
        let (maps, cnt) = std::thread::scope(|s| -> Result<_> {
            let producer = s.spawn(move || -> Result<()> {
                let mut batch = Batch { data: Vec::with_capacity(BATCH_BYTES), ends: Vec::new() };
                while bam::read_raw_record(&mut reader, &mut batch.data)? {
                    batch.ends.push(batch.data.len());
                    if batch.data.len() >= BATCH_BYTES {
                        let full = std::mem::replace(
                            &mut batch,
                            Batch { data: Vec::with_capacity(BATCH_BYTES), ends: Vec::new() },
                        );
                        if tx.send(full).is_err() {
                            break;
                        }
                    }
                }
                if !batch.ends.is_empty() {
                    let _ = tx.send(batch);
                }
                Ok(())
            });
            let workers: Vec<_> = (0..compute)
                .map(|_| {
                    let rx = rx.clone();
                    let (cov, params, chrom_names, budget) = (&cov, &params, &chrom_names, &budget);
                    s.spawn(move || -> Result<_> {
                        let mut w = Worker::new(cov, params, chrom_names);
                        for batch in rx {
                            let _permit = budget.acquire();
                            let mut from = 0;
                            for &to in &batch.ends {
                                w.record(&batch.data[from..to])?;
                                from = to;
                            }
                        }
                        Ok((w.out, w.cnt))
                    })
                })
                .collect();
            drop(rx);
            let mut maps = Vec::new();
            let mut cnt = signal::Counters::default();
            for w in workers {
                let (m, c) = w.join().expect("worker panicked")?;
                maps.extend(m);
                cnt += c;
            }
            producer.join().expect("reader panicked")?;
            Ok((maps, cnt))
        })?;
        (chrom_names, cov, maps, cnt)
    };
    eprintln!(
        "[{:.1}s] read {} records: {} candidate mappings, {} with signal{}",
        t0.elapsed().as_secs_f64(),
        counters.records,
        counters.candidates,
        counters.mappings,
        if counters.mm_overflow > 0 { format!(", {} skipped (more MM calls than target bases)", counters.mm_overflow) } else { String::new() }
    );

    let (maps, missing_mq) = signal::supp_filter(maps, a.max_dist);
    if missing_mq > 0 {
        eprintln!("warning: {missing_mq} multi-mapping reads without SA tag kept without overlap check");
    }
    eprintln!("[{:.1}s] {} mappings after supp_filter", t0.elapsed().as_secs_f64(), maps.len());

    // Outputs run one after the other so they never exceed the core budget.
    if let Some(p) = &a.rds {
        write_rds(p, &maps, &chrom_names, a.rds_level, threads)?;
        eprintln!("[{:.1}s] wrote {}", t0.elapsed().as_secs_f64(), p.display());
    }
    if let Some(p) = &a.bw {
        cov.write_bigwig(p, 1)?;
        eprintln!("[{:.1}s] wrote {}", t0.elapsed().as_secs_f64(), p.display());
    }
    Ok(())
}

fn read_header<R: std::io::Read>(r: &mut R, cov_bin: u32) -> Result<(Vec<String>, coverage::Coverage)> {
    let header = bam::read_header(r)?;
    let names = header.references.iter().map(|r| r.name.clone()).collect();
    Ok((names, coverage::Coverage::new(&header, cov_bin)))
}

/// Per-thread record processing: coverage + signal extraction.
struct Worker<'a> {
    cov: &'a coverage::Coverage,
    params: &'a signal::Params,
    chrom_names: &'a [String],
    scratch: signal::Scratch,
    out: Vec<signal::Mapping>,
    cnt: signal::Counters,
}

impl<'a> Worker<'a> {
    fn new(cov: &'a coverage::Coverage, params: &'a signal::Params, chrom_names: &'a [String]) -> Self {
        Worker { cov, params, chrom_names, scratch: Default::default(), out: Vec::new(), cnt: Default::default() }
    }

    fn record(&mut self, data: &[u8]) -> Result<()> {
        let rec = bam::Record::parse(data)?;
        self.cnt.records += 1;
        if rec.flag & (bam::FLAG_UNMAPPED | bam::FLAG_SECONDARY) != 0 || rec.ref_id < 0 {
            return Ok(());
        }
        self.cov.add(&rec);
        let name = &self.chrom_names[rec.ref_id as usize];
        if let Some(m) = signal::extract(&rec, name, self.params, &mut self.scratch, &mut self.cnt) {
            self.out.push(m);
        }
        Ok(())
    }
}

/// Serializes in memory, then gzip-compresses chunks in parallel as
/// concatenated gzip members (read transparently by R's gzfile/readRDS).
fn write_rds(path: &Path, maps: &[signal::Mapping], levels: &[String], level: u32, threads: usize) -> Result<()> {
    let mut w = rds::RdsWriter::new(Vec::with_capacity(64 << 20))?;
    w.write_alldata(maps, levels)?;
    let raw = w.finish();
    const CHUNK: usize = 16 << 20;
    let chunks: Vec<&[u8]> = raw.chunks(CHUNK).collect();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut parts: Vec<Vec<u8>> = vec![Vec::new(); chunks.len()];
    let slots: Vec<std::sync::Mutex<&mut Vec<u8>>> = parts.iter_mut().map(std::sync::Mutex::new).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.min(chunks.len()).max(1) {
            s.spawn(|| -> std::io::Result<()> {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= chunks.len() {
                        return Ok(());
                    }
                    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
                    enc.write_all(chunks[i])?;
                    **slots[i].lock().unwrap() = enc.finish()?;
                }
            });
        }
    });
    drop(slots);
    let mut f = BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    for p in &parts {
        if p.is_empty() && !raw.is_empty() {
            bail!("RDS compression failed");
        }
        f.write_all(p)?;
    }
    f.flush()?;
    Ok(())
}

/// Base-level comparison of two bigWigs (interval layout may differ).
fn bw_compare(a: &Path, b: &Path) -> Result<()> {
    let ra = coverage::read_bigwig(a)?;
    let rb = coverage::read_bigwig(b)?;
    let mut total_diff = 0u64;
    for (name, len, ia) in &ra {
        let Some((_, lenb, ib)) = rb.iter().find(|(n, _, _)| n == name) else {
            println!("{name}: missing in {}", b.display());
            continue;
        };
        if len != lenb {
            println!("{name}: length {len} vs {lenb}");
        }
        let expand = |iv: &[(u32, u32, f32)]| {
            let mut v = vec![f32::NAN; *len as usize];
            for &(s, e, x) in iv {
                for p in s..e.min(*len) {
                    v[p as usize] = x;
                }
            }
            v
        };
        let (va, vb) = (expand(ia), expand(ib));
        let mut ndiff = 0u64;
        let mut first = None;
        for (p, (x, y)) in va.iter().zip(&vb).enumerate() {
            let same = (x.is_nan() && y.is_nan()) || x == y;
            if !same {
                ndiff += 1;
                first.get_or_insert((p, *x, *y));
            }
        }
        total_diff += ndiff;
        println!(
            "{name}\tlen={len}\tintervals={}/{}\tdiff_bases={ndiff}{}",
            ia.len(),
            ib.len(),
            first.map_or(String::new(), |(p, x, y)| format!("\tfirst_diff@{p}: {x} vs {y}"))
        );
    }
    for (name, _, _) in &rb {
        if !ra.iter().any(|(n, _, _)| n == name) {
            println!("{name}: missing in {}", a.display());
        }
    }
    println!("TOTAL differing bases: {total_diff}");
    Ok(())
}
