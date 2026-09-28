//! Parallel BGZF decompression with libdeflate, exposed as an ordered `Read`.
//!
//! A reader thread groups compressed BGZF blocks, worker threads inflate the
//! groups (with CRC32 check), and `ParallelBgzfReader` hands the decompressed
//! stream back in file order. `InlineBgzfReader` is the single-thread variant.

use crate::budget::Budget;
use anyhow::{Context, Result, bail};
use crossbeam_channel::{Receiver, bounded};
use std::collections::BTreeMap;
use std::io::{self, Read};
use std::thread::JoinHandle;

const GROUP_BYTES: usize = 1 << 20; // ~1 MB compressed per work unit

type Group = (u64, Vec<u8>, Vec<(usize, usize)>);
type Inflated = (u64, Result<Vec<u8>>);

pub struct ParallelBgzfReader {
    rx: Receiver<Inflated>,
    pending: BTreeMap<u64, Result<Vec<u8>>>,
    next: u64,
    cur: Vec<u8>,
    pos: usize,
    handles: Vec<JoinHandle<()>>,
    reader: Option<JoinHandle<Result<()>>>,
}

/// Reads one BGZF block; returns false at EOF.
fn read_block<R: Read>(r: &mut R, buf: &mut Vec<u8>) -> Result<bool> {
    let start = buf.len();
    let mut head = [0u8; 12];
    match r.read_exact(&mut head) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e.into()),
    }
    if head[0] != 31 || head[1] != 139 || head[3] & 4 == 0 {
        bail!("not a BGZF file (bad block header)");
    }
    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    r.read_exact(&mut extra)?;
    let mut bsize = None;
    let mut o = 0;
    while o + 4 <= xlen {
        let slen = u16::from_le_bytes([extra[o + 2], extra[o + 3]]) as usize;
        if extra[o] == b'B' && extra[o + 1] == b'C' && slen == 2 {
            bsize = Some(u16::from_le_bytes([extra[o + 4], extra[o + 5]]) as usize);
        }
        o += 4 + slen;
    }
    let total = bsize.context("BGZF block without BSIZE")? + 1;
    buf.extend_from_slice(&head);
    buf.extend_from_slice(&extra);
    let rest = total - 12 - xlen;
    let at = buf.len();
    buf.resize(at + rest, 0);
    r.read_exact(&mut buf[at..]).context("truncated BGZF block")?;
    debug_assert_eq!(buf.len() - start, total);
    Ok(true)
}

fn inflate_group(d: &mut libdeflater::Decompressor, data: &[u8], blocks: &[(usize, usize)]) -> Result<Vec<u8>> {
    let total: usize = blocks
        .iter()
        .map(|&(_, e)| u32::from_le_bytes(data[e - 4..e].try_into().unwrap()) as usize)
        .sum();
    let mut out = Vec::with_capacity(total);
    for &(s, e) in blocks {
        let b = &data[s..e];
        let xlen = u16::from_le_bytes([b[10], b[11]]) as usize;
        let isize = u32::from_le_bytes(b[b.len() - 4..].try_into().unwrap()) as usize;
        let crc = u32::from_le_bytes(b[b.len() - 8..b.len() - 4].try_into().unwrap());
        let cdata = &b[12 + xlen..b.len() - 8];
        let at = out.len();
        out.resize(at + isize, 0);
        if isize > 0 {
            let n = d
                .deflate_decompress(cdata, &mut out[at..])
                .map_err(|e| anyhow::anyhow!("BGZF inflate error: {e:?}"))?;
            if n != isize {
                bail!("BGZF block size mismatch");
            }
        }
        if libdeflater::crc32(&out[at..]) != crc {
            bail!("BGZF CRC mismatch");
        }
    }
    Ok(out)
}

impl ParallelBgzfReader {
    /// `workers` inflate threads, each holding a `budget` permit while inflating.
    pub fn new<R: Read + Send + 'static>(mut inner: R, workers: usize, budget: Budget) -> Self {
        let workers = workers.max(1);
        let (gtx, grx) = bounded::<Group>(workers * 4);
        let (otx, orx) = bounded::<Inflated>(workers * 4);
        let reader = std::thread::spawn(move || -> Result<()> {
            let mut seq = 0u64;
            loop {
                let mut data = Vec::with_capacity(GROUP_BYTES + (64 << 10));
                let mut blocks = Vec::new();
                let mut eof = false;
                while data.len() < GROUP_BYTES {
                    let s = data.len();
                    if !read_block(&mut inner, &mut data)? {
                        eof = true;
                        break;
                    }
                    blocks.push((s, data.len()));
                }
                if !blocks.is_empty() && gtx.send((seq, data, blocks)).is_err() {
                    return Ok(());
                }
                seq += 1;
                if eof {
                    return Ok(());
                }
            }
        });
        let handles = (0..workers)
            .map(|_| {
                let (grx, otx, budget) = (grx.clone(), otx.clone(), budget.clone());
                std::thread::spawn(move || {
                    let mut d = libdeflater::Decompressor::new();
                    for (seq, data, blocks) in grx {
                        let out = {
                            let _permit = budget.acquire();
                            inflate_group(&mut d, &data, &blocks)
                        };
                        if otx.send((seq, out)).is_err() {
                            return;
                        }
                    }
                })
            })
            .collect();
        ParallelBgzfReader {
            rx: orx,
            pending: BTreeMap::new(),
            next: 0,
            cur: Vec::new(),
            pos: 0,
            handles,
            reader: Some(reader),
        }
    }

    /// Loads the next group in order; false at end of stream.
    fn advance(&mut self) -> io::Result<bool> {
        loop {
            if let Some(r) = self.pending.remove(&self.next) {
                self.next += 1;
                self.cur = r.map_err(|e| io::Error::other(e.to_string()))?;
                self.pos = 0;
                return Ok(true);
            }
            match self.rx.recv() {
                Ok((seq, r)) => {
                    self.pending.insert(seq, r);
                }
                Err(_) => {
                    // all workers done: surface a reader error if any
                    if let Some(h) = self.reader.take()
                        && let Ok(Err(e)) = h.join() {
                            return Err(io::Error::other(e.to_string()));
                        }
                    return Ok(false);
                }
            }
        }
    }
}

impl Read for ParallelBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.cur.len() {
            if !self.advance()? {
                return Ok(0);
            }
        }
        let n = buf.len().min(self.cur.len() - self.pos);
        buf[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for ParallelBgzfReader {
    fn drop(&mut self) {
        // unblock workers by dropping the receiver first
        let (_, dummy) = bounded(0);
        drop(std::mem::replace(&mut self.rx, dummy));
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
        if let Some(h) = self.reader.take() {
            let _ = h.join();
        }
    }
}

/// Single-threaded BGZF reader (used with `-t 1`): inflates block by block
/// in the calling thread, no helper threads.
pub struct InlineBgzfReader<R> {
    inner: R,
    d: libdeflater::Decompressor,
    block: Vec<u8>,
    cur: Vec<u8>,
    pos: usize,
}

impl<R: Read> InlineBgzfReader<R> {
    pub fn new(inner: R) -> Self {
        InlineBgzfReader { inner, d: libdeflater::Decompressor::new(), block: Vec::new(), cur: Vec::new(), pos: 0 }
    }
}

impl<R: Read> Read for InlineBgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.cur.len() {
            self.block.clear();
            if !read_block(&mut self.inner, &mut self.block).map_err(|e| io::Error::other(e.to_string()))? {
                return Ok(0);
            }
            let n = self.block.len();
            self.cur = inflate_group(&mut self.d, &self.block, &[(0, n)]).map_err(|e| io::Error::other(e.to_string()))?;
            self.pos = 0;
        }
        let n = buf.len().min(self.cur.len() - self.pos);
        buf[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}
