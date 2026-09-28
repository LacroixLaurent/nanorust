//! Replica of deeptools 3.5.4 `bamCoverage -b in.bam -o out.bw` with default
//! settings (binSize 50, no normalisation, no extension, all mapped reads).
//!
//! Per read, aligned blocks come from pysam `get_blocks()` (M/=/X blocks, split
//! on D/N) and each read increments every 50 bp bin touched by a block at most
//! once (deeptools `last_eIdx` logic).

use crate::bam::{Header, Record};
use anyhow::{Result, anyhow};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigRead, BigWigWrite, Value};
use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

pub struct Coverage {
    pub bin: u32,
    pub chroms: Vec<(String, u32)>,
    counts: Vec<Vec<AtomicU32>>,
}

impl Coverage {
    pub fn new(header: &Header, bin: u32) -> Self {
        let chroms: Vec<(String, u32)> =
            header.references.iter().map(|r| (r.name.clone(), r.len)).collect();
        let counts = chroms
            .iter()
            .map(|(_, len)| (0..len.div_ceil(bin)).map(|_| AtomicU32::new(0)).collect())
            .collect();
        Coverage { bin, chroms, counts }
    }

    /// Adds one mapped record (caller filters flags).
    pub fn add(&self, rec: &Record) {
        let Some(bins) = self.counts.get(rec.ref_id as usize) else { return };
        let nbins = bins.len() as i64;
        let bs = self.bin as i64;
        let mut pos = rec.pos as i64;
        let mut last_e: Option<i64> = None;
        for &(op, l) in &rec.cigar {
            let l = l as i64;
            match op {
                // M, =, X: aligned block
                0 | 7 | 8 => {
                    let (fs, fe) = (pos, pos + l);
                    pos += l;
                    if fe - fs == 0 {
                        continue;
                    }
                    let mut s = fs.div_euclid(bs).max(0);
                    let e = ((fe + bs - 1).div_euclid(bs)).min(nbins);
                    if let Some(le) = last_e {
                        s = s.max(le);
                        if s >= e {
                            continue;
                        }
                    }
                    for b in &bins[s as usize..e.max(s) as usize] {
                        b.fetch_add(1, Ordering::Relaxed);
                    }
                    last_e = Some(e);
                }
                // D, N: advance reference
                2 | 3 => pos += l,
                _ => {}
            }
        }
    }

    /// bedGraph-like runs of equal values, per chromosome in header order.
    fn intervals(&self) -> Vec<(String, Value)> {
        let mut out = Vec::new();
        for ((name, len), bins) in self.chroms.iter().zip(&self.counts) {
            let mut i = 0usize;
            while i < bins.len() {
                let v = bins[i].load(Ordering::Relaxed);
                let mut j = i + 1;
                while j < bins.len() && bins[j].load(Ordering::Relaxed) == v {
                    j += 1;
                }
                let start = i as u32 * self.bin;
                let end = (j as u32 * self.bin).min(*len);
                out.push((name.clone(), Value { start, end, value: v as f32 }));
                i = j;
            }
        }
        out
    }

    pub fn write_bigwig(&self, path: &Path, threads: usize) -> Result<()> {
        let sizes: HashMap<String, u32> = self.chroms.iter().cloned().collect();
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(threads.max(1))
            .build()?;
        let vals = BedParserStreamingIterator::wrap_infallible_iter(self.intervals().into_iter(), true);
        let out = BigWigWrite::create_file(path, sizes)?;
        out.write(vals, runtime).map_err(|e| anyhow!("bigWig write failed: {e}"))?;
        Ok(())
    }
}

/// Reads all intervals of a bigWig: chrom -> (len, [(start, end, value)]).
pub fn read_bigwig(path: &Path) -> Result<Vec<(String, u32, Vec<(u32, u32, f32)>)>> {
    let mut bw = BigWigRead::open_file(path).map_err(|e| anyhow!("{e}"))?;
    let chroms: Vec<(String, u32)> = bw.chroms().iter().map(|c| (c.name.clone(), c.length)).collect();
    let mut out = Vec::new();
    for (name, len) in chroms {
        let iv = bw
            .get_interval(&name, 0, len)
            .map_err(|e| anyhow!("{e}"))?
            .map(|v| v.map(|v| (v.start, v.end, v.value)))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow!("{e}"))?;
        out.push((name, len, iv));
    }
    Ok(out)
}
