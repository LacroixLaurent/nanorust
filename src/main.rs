mod bam;
mod bgzf;
mod coverage;
mod rds;
mod rstats;
mod signal;

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
    /// Worker threads (default: all cores)
    #[arg(short, long)]
    threads: Option<usize>,
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
    let mut reader = bgzf::ParallelBgzfReader::new(std::io::BufReader::with_capacity(1 << 20, input), threads);
    let header = bam::read_header(&mut reader)?;
    let chrom_names: Vec<String> = header.references.iter().map(|r| r.name.clone()).collect();
    let cov = coverage::Coverage::new(&header, a.cov_bin_size);
    let params = signal::Params {
        min_len: a.min_len,
        bin_size: a.bin_size,
        keep_supplementary: !a.no_supplementary,
        exclude_prefixes: a.exclude_prefix.clone(),
    };

    let (tx, rx) = crossbeam_channel::bounded::<Batch>(threads * 2);
    let (maps, counters) = std::thread::scope(|s| -> Result<_> {
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
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let rx = rx.clone();
                let (cov, params, chrom_names) = (&cov, &params, &chrom_names);
                s.spawn(move || -> Result<_> {
                    let mut out = Vec::new();
                    let mut cnt = signal::Counters::default();
                    let mut scratch = signal::Scratch::default();
                    for batch in rx {
                        let mut from = 0;
                        for &to in &batch.ends {
                            let rec = bam::Record::parse(&batch.data[from..to])?;
                            from = to;
                            cnt.records += 1;
                            if rec.flag & (bam::FLAG_UNMAPPED | bam::FLAG_SECONDARY) != 0 || rec.ref_id < 0 {
                                continue;
                            }
                            cov.add(&rec);
                            let name = &chrom_names[rec.ref_id as usize];
                            if let Some(m) = signal::extract(&rec, name, params, &mut scratch, &mut cnt) {
                                out.push(m);
                            }
                        }
                    }
                    Ok((out, cnt))
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
    eprintln!(
        "[{:.1}s] read {} records: {} candidate mappings, {} with signal{}",
        t0.elapsed().as_secs_f64(),
        counters.records,
        counters.candidates,
        counters.mappings,
        if counters.mm_overflow > 0 { format!(", {} skipped (MM longer than T count)", counters.mm_overflow) } else { String::new() }
    );

    let (maps, missing_mq) = signal::supp_filter(maps, a.max_dist);
    if missing_mq > 0 {
        eprintln!("warning: {missing_mq} multi-mapping reads without SA tag kept without overlap check");
    }
    eprintln!("[{:.1}s] {} mappings after supp_filter", t0.elapsed().as_secs_f64(), maps.len());

    std::thread::scope(|s| -> Result<()> {
        let bw_job = a.bw.as_ref().map(|p| s.spawn(|| cov.write_bigwig(p, 2)));
        if let Some(p) = &a.rds {
            write_rds(p, &maps, &chrom_names, a.rds_level, threads)?;
            eprintln!("[{:.1}s] wrote {}", t0.elapsed().as_secs_f64(), p.display());
        }
        if let Some(j) = bw_job {
            j.join().expect("bigWig writer panicked")?;
            eprintln!("[{:.1}s] wrote {}", t0.elapsed().as_secs_f64(), a.bw.as_ref().unwrap().display());
        }
        Ok(())
    })
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
