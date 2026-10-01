//! Minimal writer for R's serialization format (version 3, XDR), as produced
//! by `saveRDS()` (gzip-compressed).

use crate::signal::{Mapping, ModSpec};
use std::collections::HashMap;
use std::io::{self, Write};

const NILVALUE_SXP: i32 = 254;
const REFSXP: i32 = 255;
const SYMSXP: i32 = 1;
const LISTSXP: i32 = 2;
const CHARSXP: i32 = 9;
const INTSXP: i32 = 13;
const REALSXP: i32 = 14;
const STRSXP: i32 = 16;
const VECSXP: i32 = 19;
const IS_OBJECT: i32 = 1 << 8;
const HAS_ATTR: i32 = 1 << 9;
const HAS_TAG: i32 = 1 << 10;
const ASCII_MASK: i32 = 1 << 6;
const UTF8_MASK: i32 = 1 << 3;
const NA_INTEGER: i32 = i32::MIN;

pub struct RdsWriter<W: Write> {
    w: W,
    symbols: HashMap<&'static str, i32>,
}

impl<W: Write> RdsWriter<W> {
    pub fn new(mut w: W) -> io::Result<Self> {
        w.write_all(b"X\n")?;
        let mut s = RdsWriter { w, symbols: HashMap::new() };
        s.int(3)?;
        s.int(0x0004_0403)?;
        s.int(0x0003_0500)?;
        s.int(5)?;
        s.w.write_all(b"UTF-8")?;
        Ok(s)
    }

    pub fn finish(self) -> W {
        self.w
    }

    #[inline]
    fn int(&mut self, v: i32) -> io::Result<()> {
        self.w.write_all(&v.to_be_bytes())
    }

    fn len(&mut self, n: usize) -> io::Result<()> {
        if n > i32::MAX as usize {
            self.int(-1)?;
            self.int((n >> 32) as i32)?;
            self.int(n as u32 as i32)
        } else {
            self.int(n as i32)
        }
    }

    fn charsxp(&mut self, s: &str) -> io::Result<()> {
        let gp = if s.is_ascii() { ASCII_MASK } else { UTF8_MASK };
        self.int(CHARSXP | (gp << 12))?;
        self.int(s.len() as i32)?;
        self.w.write_all(s.as_bytes())
    }

    fn symbol(&mut self, name: &'static str) -> io::Result<()> {
        if let Some(&idx) = self.symbols.get(name) {
            return self.int((idx << 8) | REFSXP);
        }
        let idx = self.symbols.len() as i32 + 1;
        self.symbols.insert(name, idx);
        self.int(SYMSXP)?;
        self.charsxp(name)
    }

    fn strsxp<S: AsRef<str>>(&mut self, v: &[S]) -> io::Result<()> {
        self.int(STRSXP)?;
        self.len(v.len())?;
        for s in v {
            self.charsxp(s.as_ref())?;
        }
        Ok(())
    }

    fn intsxp(&mut self, flags: i32, v: impl ExactSizeIterator<Item = i32>) -> io::Result<()> {
        self.int(INTSXP | flags)?;
        self.len(v.len())?;
        for x in v {
            self.int(x)?;
        }
        Ok(())
    }

    fn realsxp(&mut self, v: impl ExactSizeIterator<Item = f64>) -> io::Result<()> {
        self.int(REALSXP)?;
        self.len(v.len())?;
        for x in v {
            self.w.write_all(&x.to_bits().to_be_bytes())?;
        }
        Ok(())
    }

    fn attr_tag(&mut self, name: &'static str) -> io::Result<()> {
        self.int(LISTSXP | HAS_TAG)?;
        self.symbol(name)
    }

    fn end_attrs(&mut self) -> io::Result<()> {
        self.int(NILVALUE_SXP)
    }

    fn tibble_attrs(&mut self, names: &[&str], nrow: usize, compact_rownames: bool) -> io::Result<()> {
        self.attr_tag("row.names")?;
        if compact_rownames {
            self.intsxp(0, [NA_INTEGER, -(nrow as i32)].into_iter())?;
        } else {
            self.intsxp(0, 1..nrow as i32 + 1)?;
        }
        self.attr_tag("names")?;
        self.strsxp(names)?;
        self.attr_tag("class")?;
        self.strsxp(&["tbl_df", "tbl", "data.frame"])?;
        self.end_attrs()
    }

    fn factor(&mut self, codes: impl ExactSizeIterator<Item = i32>, levels: &[String]) -> io::Result<()> {
        self.intsxp(IS_OBJECT | HAS_ATTR, codes)?;
        self.attr_tag("levels")?;
        self.strsxp(levels)?;
        self.attr_tag("class")?;
        self.strsxp(&["factor"])?;
        self.end_attrs()
    }

    /// Écrit la tibble `alldata` regroupant toutes les modifications d'un même read dans une sous-tibble unique `signalbin`.
    pub fn write_alldata(
        &mut self,
        maps: &[Mapping],
        chrom_levels: &[String],
        modspecs: &[ModSpec],
    ) -> io::Result<()> {
        let n = maps.len();
        let num_mods = modspecs.len();

        // 1. Colonnes de la sous-tibble: "positions", "signalB", "signalE", ...
        let signalbin_col_names: Vec<String> = std::iter::once("positions".to_string())
            .chain(modspecs.iter().map(|ms| format!("signal{}", String::from_utf8_lossy(&ms.code))))
            .collect();
        let signalbin_col_refs: Vec<&str> = signalbin_col_names.iter().map(|s| s.as_str()).collect();

        // 2. Colonnes de la tibble principale
        let mut top_col_names: Vec<String> = vec![
            "read_id".into(),
            "flag".into(),
            "chrom".into(),
            "strand".into(),
            "start".into(),
            "end".into(),
            "signalbin".into(),
        ];

        if num_mods == 1 {
            top_col_names.push("med_signal".into());
            top_col_names.push("med_signalbin".into());
        } else {
            for ms in modspecs {
                let code_str = String::from_utf8_lossy(&ms.code);
                top_col_names.push(format!("med_signal{}", code_str));
                top_col_names.push(format!("med_signalbin{}", code_str));
            }
        }
        let top_col_refs: Vec<&str> = top_col_names.iter().map(|s| s.as_str()).collect();

        self.int(VECSXP | IS_OBJECT | HAS_ATTR)?;
        self.len(top_col_names.len())?;

        // 1. read_id
        self.int(STRSXP)?;
        self.len(n)?;
        for m in maps {
            self.charsxp(&m.read_id)?;
        }
        // 2. flag
        self.intsxp(0, maps.iter().map(|m| m.flag as i32))?;
        // 3. chrom
        self.factor(maps.iter().map(|m| m.chrom as i32 + 1), chrom_levels)?;
        // 4. strand
        let strand_levels = ["+".to_string(), "-".to_string(), "*".to_string()];
        self.factor(maps.iter().map(|m| if m.minus { 2 } else { 1 }), &strand_levels)?;
        // 5. start
        self.realsxp(maps.iter().map(|m| m.start as f64))?;
        // 6. end
        self.realsxp(maps.iter().map(|m| m.end as f64))?;
        // 7. signalbin (liste de tibbles multi-colonnes)
        self.int(VECSXP)?;
        self.len(n)?;
        for m in maps {
            self.int(VECSXP | IS_OBJECT | HAS_ATTR)?;
            self.len(1 + num_mods)?;
            self.realsxp(m.bin_positions.iter().copied())?;
            for mod_idx in 0..num_mods {
                let sigs = &m.mod_bin_signals[mod_idx];
                self.realsxp(sigs.iter().copied())?;
            }
            self.tibble_attrs(&signalbin_col_refs, m.bin_positions.len(), false)?;
        }
        // 8.. med_signal / med_signalbin par modification
        for mod_idx in 0..num_mods {
            self.realsxp(maps.iter().map(|m| m.med_signals.get(mod_idx).copied().unwrap_or(f64::NAN)))?;
            self.realsxp(maps.iter().map(|m| m.med_signalbins.get(mod_idx).copied().unwrap_or(f64::NAN)))?;
        }

        self.tibble_attrs(&top_col_refs, n, true)
    }
}
