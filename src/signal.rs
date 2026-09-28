//! Per-mapping BrdU signal extraction, mirroring `parsing_DoradoRemora_v18_Br.r`
//! (read.bam.batch -> process.mod.tag -> process.signal -> signalbin).

use crate::bam::{self, BASE_A, BASE_T, OP_D, OP_I, OP_M, OP_N, OP_S, Record, Tag};
use crate::rstats::{r_mean, r_median, r_median_ml};

pub struct Params {
    pub min_len: i64,
    pub bin_size: i64,
    pub keep_supplementary: bool,
    /// chromosomes whose name starts with any of these are skipped (step 01 `$3 !~ /^chrM/`)
    pub exclude_prefixes: Vec<String>,
}

pub struct Mapping {
    pub read_id: Box<str>,
    pub flag: u16,
    pub chrom: u32,
    pub minus: bool,
    pub start: i64,
    pub end: i64,
    /// (bin start position, mean BrdU probability)
    pub bins: Vec<(f64, f64)>,
    pub med_signal: f64,
    pub med_signalbin: f64,
    /// Query intervals [from, to) (1-based) covered by M ops, i.e. `na.omit(mapping_table)$read_pos`.
    /// Only kept for reads with an SA tag, as it is only needed by the supplementary filter.
    pub m_query: Option<Vec<(u32, u32)>>,
}

#[derive(Default, Clone, Copy)]
pub struct Counters {
    pub records: u64,
    pub candidates: u64,
    pub mappings: u64,
    pub mm_overflow: u64,
}

impl std::ops::AddAssign for Counters {
    fn add_assign(&mut self, o: Self) {
        self.records += o.records;
        self.candidates += o.candidates;
        self.mappings += o.mappings;
        self.mm_overflow += o.mm_overflow;
    }
}

/// An M block from R's `parseCigar`: query start (1-based), length, reference start (1-based).
struct Block {
    q: i64,
    len: i64,
    r: i64,
}

pub fn extract(
    rec: &Record,
    chrom_name: &str,
    p: &Params,
    scratch: &mut Scratch,
    cnt: &mut Counters,
) -> Option<Mapping> {
    let flag = rec.flag;
    if p.exclude_prefixes.iter().any(|pre| chrom_name.starts_with(pre.as_str())) {
        return None;
    }
    let minus = flag & bam::FLAG_REVERSE != 0;
    let is_primary = flag == 0 || flag == 16;
    let is_supp = flag == 2048 || flag == 2064;
    if !is_primary && !(is_supp && p.keep_supplementary) {
        return None;
    }

    let (mm, ml) = match (bam::find_tag(rec.aux, b"MM"), bam::find_tag(rec.aux, b"ML")) {
        (Some(Tag::Str(mm)), Some(Tag::Array(b'C', ml))) => (mm, ml),
        _ => match (bam::find_tag(rec.aux, b"Mm"), bam::find_tag(rec.aux, b"Ml")) {
            (Some(Tag::Str(mm)), Some(Tag::Array(b'C', ml))) => (mm, ml),
            _ => return None,
        },
    };
    if ml.is_empty() {
        return None;
    }
    let sa = match bam::find_tag(rec.aux, b"SA") {
        Some(Tag::Str(s)) => Some(s),
        _ => None,
    };
    if is_supp {
        // keep only if the first SA entry is on the same chrom and strand
        let sa = sa?;
        let mut f = sa.split(|&c| c == b',');
        let sa_chr = f.next()?;
        let _pos = f.next();
        let sa_str = f.next()?;
        if sa_chr != chrom_name.as_bytes() || sa_str != (if minus { b"-" } else { b"+" }) {
            return None;
        }
    }

    // rlen = sum of M and D ops
    let rlen: i64 = rec
        .cigar
        .iter()
        .filter(|(op, _)| *op == OP_M || *op == OP_D)
        .map(|&(_, l)| l as i64)
        .sum();
    let start = rec.pos as i64 + 1;
    let end = start + rlen - 1;
    if end - start <= p.min_len {
        return None;
    }
    cnt.candidates += 1;

    // --- parseCigar ---
    let blocks = &mut scratch.blocks;
    blocks.clear();
    let (mut rpos, mut qidx) = (start, 1i64);
    let mut last_op = u8::MAX;
    let mut last_len = 0i64;
    for &(op, l) in &rec.cigar {
        let l = l as i64;
        match op {
            OP_M => {
                blocks.push(Block { q: qidx, len: l, r: rpos });
                rpos += l;
                qidx += l;
            }
            OP_I | OP_S => qidx += l,
            OP_D | OP_N => rpos += l,
            _ => {}
        }
        last_op = op;
        last_len = l;
    }
    // `max(mapping_table$read_pos)` used to flip minus-strand positions
    let read_length = if last_op == OP_S && last_len > 0 {
        qidx - 1
    } else {
        match blocks.last() {
            Some(b) => b.q + b.len - 1,
            None => return None,
        }
    };

    // --- process.mod.tag: probabilities for each T in read orientation ---
    // In read orientation, T's are BAM-seq T's (+) or BAM-seq A's read backwards (-).
    let target = if minus { BASE_A } else { BASE_T };
    let tpos = &mut scratch.tpos; // BAM-orientation 0-based indices, in read-orientation order
    tpos.clear();
    if minus {
        for i in (0..rec.l_seq).rev() {
            if rec.base_code(i) == target {
                tpos.push(i as u32);
            }
        }
    } else {
        for i in 0..rec.l_seq {
            if rec.base_code(i) == target {
                tpos.push(i as u32);
            }
        }
    }
    // prob code per T: 0..=255, or NONE (unreported with '?')
    const NONE: u16 = u16::MAX;
    let probs = &mut scratch.probs;
    probs.clear();
    let mut ml_off = 0usize;
    let mut found_b = false;
    let mm = mm.strip_suffix(b";").unwrap_or(mm);
    for entry in mm.split(|&c| c == b';') {
        let mut fields = entry.split(|&c| c == b',');
        let code = fields.next().unwrap_or(b"");
        let n_rel = fields.clone().count();
        let (code_core, unreported) = match code.last() {
            Some(b'?') => (&code[..code.len() - 1], NONE),
            Some(b'.') => (&code[..code.len() - 1], 0u16),
            _ => (code, 0u16),
        };
        if code_core == b"T+B" && !found_b {
            found_b = true;
            probs.resize(tpos.len(), unreported);
            let mut idx = 0usize;
            for (i, f) in fields.enumerate() {
                let rel: usize = std::str::from_utf8(f).ok()?.trim().parse().ok()?;
                idx += rel + 1;
                let v = *ml.get(ml_off + i)? as u16;
                if idx > probs.len() {
                    // R would fail on this record (more calls than T's); skip it
                    cnt.mm_overflow += 1;
                    return None;
                }
                probs[idx - 1] = v;
            }
        }
        ml_off += n_rel;
    }
    if !found_b {
        return None;
    }

    // --- process.signal: map to reference, keep positions in [start, end] ---
    // Iterate T's in ascending BAM query order so reference positions are ascending.
    let signal = &mut scratch.signal;
    signal.clear();
    let shift = read_length - rec.l_seq as i64; // mod_pos = read_length - (p - 1) with p = l_seq - j
    let mut bi = 0usize;
    let n_t = tpos.len();
    for k in 0..n_t {
        let kk = if minus { n_t - 1 - k } else { k };
        let pr = probs[kk];
        if pr == NONE {
            continue;
        }
        let j = tpos[kk] as i64 + 1; // 1-based BAM query index
        let q = if minus { j + shift } else { j };
        while bi < blocks.len() && blocks[bi].q + blocks[bi].len <= q {
            bi += 1;
        }
        if bi == blocks.len() {
            break;
        }
        let b = &blocks[bi];
        if q < b.q {
            continue;
        }
        let r = b.r + (q - b.q);
        if r >= start && r <= end {
            signal.push((r, pr as u8));
        }
    }
    if signal.is_empty() {
        return None;
    }

    // --- binning (signalbin) and medians ---
    let mut hist = [0u32; 256];
    for &(_, c) in signal.iter() {
        hist[c as usize] += 1;
    }
    let med_signal = r_median_ml(&hist, signal.len());

    let bs = p.bin_size;
    let mut bins: Vec<(f64, f64)> = Vec::new();
    let vals = &mut scratch.vals;
    vals.clear();
    let mut cur_bin = i64::MIN;
    for &(r, c) in signal.iter() {
        let b = (r - 1).div_euclid(bs) * bs + 1;
        if b != cur_bin {
            if !vals.is_empty() {
                bins.push((cur_bin as f64, r_mean(vals)));
                vals.clear();
            }
            cur_bin = b;
        }
        vals.push(c as f64 / 255.0);
    }
    bins.push((cur_bin as f64, r_mean(vals)));
    let mut bvals: Vec<f64> = bins.iter().map(|b| b.1).collect();
    let med_signalbin = r_median(&mut bvals);

    let m_query = sa.map(|_| {
        blocks
            .iter()
            .map(|b| (b.q as u32, (b.q + b.len) as u32))
            .collect()
    });

    cnt.mappings += 1;
    Some(Mapping {
        read_id: String::from_utf8_lossy(rec.name).into(),
        flag,
        chrom: rec.ref_id as u32,
        minus,
        start,
        end,
        bins,
        med_signal,
        med_signalbin,
        m_query,
    })
}

#[derive(Default)]
pub struct Scratch {
    blocks: Vec<Block>,
    tpos: Vec<u32>,
    probs: Vec<u16>,
    signal: Vec<(i64, u8)>,
    vals: Vec<f64>,
}

/// Size of the intersection of two sorted, disjoint interval lists.
fn overlap(a: &[(u32, u32)], b: &[(u32, u32)]) -> u64 {
    let (mut i, mut j, mut n) = (0, 0, 0u64);
    while i < a.len() && j < b.len() {
        let lo = a[i].0.max(b[j].0);
        let hi = a[i].1.min(b[j].1);
        if hi > lo {
            n += (hi - lo) as u64;
        }
        if a[i].1 < b[j].1 { i += 1 } else { j += 1 }
    }
    n
}

fn total(a: &[(u32, u32)]) -> u64 {
    a.iter().map(|&(s, e)| (e - s) as u64).sum()
}

/// `supp_filter()` applied per read over the whole file, then rows ordered
/// like R's `arrange(chrom, read_id, flag, start)`.
pub fn supp_filter(mut maps: Vec<Mapping>, max_dist: i64) -> (Vec<Mapping>, u64) {
    maps.sort_unstable_by(|a, b| {
        (&a.read_id, a.chrom, a.flag, a.start).cmp(&(&b.read_id, b.chrom, b.flag, b.start))
    });
    let mut keep = vec![true; maps.len()];
    let mut missing_mq = 0u64;
    let mut i = 0;
    while i < maps.len() {
        let mut j = i + 1;
        while j < maps.len() && maps[j].read_id == maps[i].read_id {
            j += 1;
        }
        if j - i > 1 {
            let first = &maps[i];
            // tokeep: read-position overlap with the first mapping is total or null
            for k in i + 1..j {
                keep[k] = match (&maps[k].m_query, &first.m_query) {
                    (Some(x), Some(f)) => {
                        let ov = overlap(x, f);
                        ov == total(f) || ov == 0
                    }
                    _ => {
                        missing_mq += 1;
                        true
                    }
                };
            }
            // distance to the first mapping
            let (s1, e1) = (first.start, first.end);
            for k in i + 1..j {
                if keep[k] {
                    let (s, e) = (maps[k].start, maps[k].end);
                    let d = (s - s1).abs().min((e - s1).abs()).min((s - e1).abs()).min((e - e1).abs());
                    keep[k] = d < max_dist;
                }
            }
        }
        i = j;
    }
    let mut it = keep.iter();
    maps.retain(|_| *it.next().unwrap());
    maps.sort_unstable_by(|a, b| {
        (a.chrom, &a.read_id, a.flag, a.start).cmp(&(b.chrom, &b.read_id, b.flag, b.start))
    });
    (maps, missing_mq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_overlap() {
        let a = [(1, 10), (12, 20)];
        let b = [(5, 15)];
        assert_eq!(overlap(&a, &b), 5 + 3);
        assert_eq!(overlap(&a, &a), total(&a));
        assert_eq!(overlap(&a, &[(30, 40)]), 0);
    }
}
