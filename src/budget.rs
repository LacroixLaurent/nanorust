//! Global CPU budget shared by all compute stages (`-t N`).
//!
//! Every CPU-heavy unit of work (inflating a BGZF group, processing a record
//! batch) holds a permit while it runs, so at most N units execute at once no
//! matter how many threads exist. Permits are never held while blocking on a
//! channel, which rules out deadlocks.

use crossbeam_channel::{Receiver, Sender, bounded};

#[derive(Clone)]
pub struct Budget {
    tx: Sender<()>,
    rx: Receiver<()>,
}

pub struct Permit<'a>(&'a Budget);

impl Budget {
    pub fn new(n: usize) -> Self {
        let n = n.max(1);
        let (tx, rx) = bounded(n);
        for _ in 0..n {
            tx.send(()).unwrap();
        }
        Budget { tx, rx }
    }

    pub fn acquire(&self) -> Permit<'_> {
        // the Budget owns both ends, so the channel can never disconnect
        self.rx.recv().unwrap();
        Permit(self)
    }
}

impl Drop for Permit<'_> {
    fn drop(&mut self) {
        let _ = self.0.tx.send(());
    }
}
