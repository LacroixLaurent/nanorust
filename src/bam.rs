//! Minimal zero-copy BAM decoding: header, raw records, CIGAR and aux tags.

use anyhow::{Context, Result, bail};
use std::io::{self, Read};

pub struct Reference {
    pub name: String,
    pub len: u32,
}

pub struct Header {
    pub references: Vec<Reference>,
}

fn read_i32<R: Read>(r: &mut R) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

pub fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic).context("reading BAM magic")?;
    if &magic != b"BAM\x01" {
        bail!("not a BAM file (bad magic)");
    }
    let l_text = read_i32(r)? as usize;
    io::copy(&mut r.by_ref().take(l_text as u64), &mut io::sink())?;
    let n_ref = read_i32(r)? as usize;
    let mut references = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let l_name = read_i32(r)? as usize;
        let mut name = vec![0u8; l_name];
        r.read_exact(&mut name)?;
        if name.last() == Some(&0) {
            name.pop();
        }
        let len = read_i32(r)? as u32;
        references.push(Reference { name: String::from_utf8(name)?, len });
    }
    Ok(Header { references })
}

/// Reads the next raw record (without its block_size prefix) into `buf`.
/// Returns false on clean EOF.
pub fn read_raw_record<R: Read>(r: &mut R, buf: &mut Vec<u8>) -> Result<bool> {
    let mut b = [0u8; 4];
    let mut got = 0;
    while got < 4 {
        let n = r.read(&mut b[got..])?;
        if n == 0 {
            if got == 0 {
                return Ok(false);
            }
            bail!("truncated BAM record length");
        }
        got += n;
    }
    let size = u32::from_le_bytes(b) as usize;
    let start = buf.len();
    buf.resize(start + size, 0);
    r.read_exact(&mut buf[start..]).context("truncated BAM record")?;
    Ok(true)
}

pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_SECONDARY: u16 = 0x100;

// CIGAR op codes (BAM encoding "MIDNSHP=X")
pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;

pub struct Record<'a> {
    pub ref_id: i32,
    /// 0-based leftmost position
    pub pos: i32,
    pub flag: u16,
    pub name: &'a [u8],
    pub cigar: Vec<(u8, u32)>,
    pub l_seq: usize,
    seq: &'a [u8],
    pub aux: &'a [u8],
}

impl<'a> Record<'a> {
    pub fn parse(d: &'a [u8]) -> Result<Self> {
        if d.len() < 32 {
            bail!("BAM record too short");
        }
        let le32 = |o: usize| i32::from_le_bytes(d[o..o + 4].try_into().unwrap());
        let le16 = |o: usize| u16::from_le_bytes(d[o..o + 2].try_into().unwrap());
        let ref_id = le32(0);
        let pos = le32(4);
        let l_read_name = d[8] as usize;
        let n_cigar = le16(12) as usize;
        let flag = le16(14);
        let l_seq = le32(16) as usize;
        let mut o = 32;
        let name_end = o + l_read_name;
        let name = &d[o..name_end - 1]; // drop NUL
        o = name_end;
        let cigar_bytes = &d[o..o + 4 * n_cigar];
        o += 4 * n_cigar;
        let seq_len_bytes = l_seq.div_ceil(2);
        let seq = &d[o..o + seq_len_bytes];
        o += seq_len_bytes + l_seq; // skip qual
        let aux = &d[o..];
        let mut cigar: Vec<(u8, u32)> = cigar_bytes
            .chunks_exact(4)
            .map(|c| {
                let v = u32::from_le_bytes(c.try_into().unwrap());
                ((v & 0xf) as u8, v >> 4)
            })
            .collect();
        // Long CIGAR (>65535 ops) is stored in the CG:B:I tag with a kSmN placeholder.
        if cigar.len() == 2 && cigar[0] == (OP_S, l_seq as u32) && cigar[1].0 == OP_N
            && let Some(Tag::Array(b'I', cg)) = find_tag(aux, b"CG") {
                cigar = cg
                    .chunks_exact(4)
                    .map(|c| {
                        let v = u32::from_le_bytes(c.try_into().unwrap());
                        ((v & 0xf) as u8, v >> 4)
                    })
                    .collect();
            }
        Ok(Record { ref_id, pos, flag, name, cigar, l_seq, seq, aux })
    }

    /// 4-bit base code at 0-based query position `i` ("=ACMGRSVTWYHKDBN").
    #[inline]
    pub fn base_code(&self, i: usize) -> u8 {
        let b = self.seq[i >> 1];
        if i & 1 == 0 { b >> 4 } else { b & 0xf }
    }
}

pub const BASE_A: u8 = 1;
pub const BASE_T: u8 = 8;

pub enum Tag<'a> {
    /// Z or H string (without NUL)
    Str(&'a [u8]),
    /// B array: subtype and raw little-endian payload
    Array(u8, &'a [u8]),
    Other,
}

fn scalar_size(t: u8) -> Option<usize> {
    match t {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        _ => None,
    }
}

/// Iterates aux fields and returns the first tag named `key`.
pub fn find_tag<'a>(aux: &'a [u8], key: &[u8; 2]) -> Option<Tag<'a>> {
    let mut o = 0;
    while o + 3 <= aux.len() {
        let k = &aux[o..o + 2];
        let t = aux[o + 2];
        o += 3;
        let (val, next) = match t {
            b'Z' | b'H' => {
                let end = o + aux[o..].iter().position(|&c| c == 0)?;
                (Tag::Str(&aux[o..end]), end + 1)
            }
            b'B' => {
                let sub = aux[o];
                let n = u32::from_le_bytes(aux[o + 1..o + 5].try_into().ok()?) as usize;
                let sz = scalar_size(sub)?;
                let start = o + 5;
                (Tag::Array(sub, &aux[start..start + n * sz]), start + n * sz)
            }
            _ => (Tag::Other, o + scalar_size(t)?),
        };
        if k == key {
            return Some(val);
        }
        o = next;
    }
    None
}
