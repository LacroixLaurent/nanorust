//! Software emulation of the x87 80-bit extended format (`long double` on
//! x86-64 Linux): 64-bit significand, round-to-nearest-even after every
//! operation. Used to reproduce R's `mean()` bit for bit as computed by R on
//! Linux x86-64. Only finite values are supported (all we need).

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F80 {
    neg: bool,
    /// value = mant * 2^exp; mant has bit 63 set unless the value is zero
    mant: u64,
    exp: i32,
}

pub const ZERO: F80 = F80 { neg: false, mant: 0, exp: 0 };

impl F80 {
    #[inline]
    pub fn from_f64(v: f64) -> F80 {
        let bits = v.to_bits();
        let neg = bits >> 63 != 0;
        let be = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & ((1u64 << 52) - 1);
        let (m, e) = if be == 0 {
            if frac == 0 {
                return ZERO;
            }
            (frac, -1074) // subnormal
        } else {
            (frac | (1u64 << 52), be - 1075)
        };
        let lz = m.leading_zeros() as i32;
        F80 { neg, mant: m << lz, exp: e - lz }
    }

    #[cfg(test)]
    pub fn from_u64(n: u64) -> F80 {
        if n == 0 {
            return ZERO;
        }
        let lz = n.leading_zeros() as i32;
        F80 { neg: false, mant: n << lz, exp: -lz }
    }

    #[inline]
    fn is_zero(self) -> bool {
        self.mant == 0
    }

    #[inline]
    pub fn neg(self) -> F80 {
        F80 { neg: !self.neg, ..self }
    }

    /// Rounds `m * 2^e` (m may carry a sticky bit in its LSB) to 64 bits, RNE.
    #[inline]
    fn round(neg: bool, m: u128, e: i32) -> F80 {
        if m == 0 {
            return ZERO;
        }
        let bits = 128 - m.leading_zeros() as i32;
        if bits <= 64 {
            let sh = 64 - bits;
            return F80 { neg, mant: (m as u64) << sh, exp: e - sh };
        }
        let sh = (bits - 64) as u32;
        let mut mant = (m >> sh) as u64;
        let rem = m & ((1u128 << sh) - 1);
        let half = 1u128 << (sh - 1);
        let mut exp = e + sh as i32;
        if rem > half || (rem == half && mant & 1 == 1) {
            mant = mant.wrapping_add(1);
            if mant == 0 {
                mant = 1 << 63;
                exp += 1;
            }
        }
        F80 { neg, mant, exp }
    }

    #[inline]
    pub fn add(self, o: F80) -> F80 {
        if self.is_zero() {
            return o;
        }
        if o.is_zero() {
            return self;
        }
        // larger magnitude first
        let (a, b) = if (self.exp, self.mant) >= (o.exp, o.mant) { (self, o) } else { (o, self) };
        // 62 guard bits below the 64-bit significand
        let ma = (a.mant as u128) << 62;
        let mut mb = (b.mant as u128) << 62;
        let d = (a.exp - b.exp) as u32;
        if d > 0 {
            if d >= 127 {
                mb = 1; // only a sticky bit remains
            } else {
                let lost = mb & ((1u128 << d) - 1) != 0;
                mb = (mb >> d) | lost as u128;
            }
        }
        let e = a.exp - 62;
        if a.neg == b.neg {
            F80::round(a.neg, ma + mb, e)
        } else {
            F80::round(a.neg, ma - mb, e)
        }
    }

    #[inline]
    pub fn sub(self, o: F80) -> F80 {
        self.add(o.neg())
    }

    /// Division by a positive integer (as R's `s /= n` with n converted to long double).
    #[inline]
    pub fn div_u64(self, n: u64) -> F80 {
        if self.is_zero() {
            return ZERO;
        }
        let num = (self.mant as u128) << 64;
        let q = num / n as u128;
        let r = num % n as u128;
        F80::round(self.neg, q | (r != 0) as u128, self.exp - 64)
    }

    /// Rounds to the nearest f64 (ties to even), like `(double) s`.
    pub fn to_f64(self) -> f64 {
        if self.is_zero() {
            return if self.neg { -0.0 } else { 0.0 };
        }
        // value = mant * 2^exp with mant in [2^63, 2^64)
        let unbiased = self.exp + 63; // exponent of the leading bit
        let sign = (self.neg as u64) << 63;
        if unbiased >= -1022 {
            let mut m = self.mant >> 11;
            let rem = self.mant & 0x7ff;
            let mut ue = unbiased;
            if rem > 0x400 || (rem == 0x400 && m & 1 == 1) {
                m += 1;
                if m == 1u64 << 53 {
                    m >>= 1;
                    ue += 1;
                }
            }
            if ue > 1023 {
                return f64::from_bits(sign | 0x7ff0_0000_0000_0000);
            }
            f64::from_bits(sign | (((ue + 1023) as u64) << 52) | (m & ((1u64 << 52) - 1)))
        } else {
            // subnormal result: significand bits available = 52 - (-1022 - unbiased)
            let shift = 11 + (-1022 - unbiased) as u32;
            if shift > 64 {
                return f64::from_bits(sign);
            }
            let m = if shift >= 64 { 0 } else { self.mant >> shift };
            let rem = if shift >= 64 { self.mant as u128 } else { (self.mant & ((1u64 << shift) - 1)) as u128 };
            let half = 1u128 << (shift - 1);
            let m = if rem > half || (rem == half && m & 1 == 1) { m + 1 } else { m };
            f64::from_bits(sign | m)
        }
    }
}

/// R's `mean()` for doubles as computed on x86-64 Linux (LDOUBLE = x87 80-bit):
/// s = sum(x) / n; t = sum(x - s); s += t / n; return (double) s.
pub fn r_mean_x87(x: &[f64]) -> f64 {
    let n = x.len() as u64;
    if n == 0 {
        return f64::NAN;
    }
    let mut s = ZERO;
    for &v in x {
        s = s.add(F80::from_f64(v));
    }
    s = s.div_u64(n);
    let mut t = ZERO;
    for &v in x {
        t = t.add(F80::from_f64(v).sub(s));
    }
    s = s.add(t.div_u64(n));
    s.to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_and_exact_ops() {
        for v in [0.0, 1.0, 0.5, 1.0 / 255.0, 254.0 / 255.0, 1e-300, 123456.789, -3.25] {
            assert_eq!(F80::from_f64(v).to_f64(), v);
        }
        let a = F80::from_f64(1.0).add(F80::from_f64(2.0));
        assert_eq!(a.to_f64(), 3.0);
        assert_eq!(F80::from_f64(1.0).sub(F80::from_f64(1.0)).to_f64(), 0.0);
        assert_eq!(F80::from_f64(1.0).div_u64(4).to_f64(), 0.25);
        assert_eq!(F80::from_u64(10).to_f64(), 10.0);
    }

    #[test]
    fn extended_precision_is_kept() {
        // 1 + 2^-60 is representable in 80-bit but not in f64
        let a = F80::from_f64(1.0).add(F80::from_f64(2f64.powi(-60)));
        assert_eq!(a.sub(F80::from_f64(1.0)).to_f64(), 2f64.powi(-60));
        // 1/3 in 80-bit, rounded to f64, equals f64 1/3
        assert_eq!(F80::from_f64(1.0).div_u64(3).to_f64(), 1.0 / 3.0);
    }

    #[test]
    fn mean_simple() {
        assert_eq!(r_mean_x87(&[1.0, 2.0, 3.0, 4.0]), 2.5);
    }
}

/// Per-run lookup of the 256 possible values (`code/255`, or 0/1 when binarised)
/// and a per-call cache of `value - mean`, to speed up `r_mean_x87` on ML codes.
pub struct CodeMean {
    x: [F80; 256],
    diff: [F80; 256],
    stamp: [u32; 256],
    epoch: u32,
}

impl CodeMean {
    pub fn new(values: &[f64; 256]) -> Self {
        let mut x = [ZERO; 256];
        for (xi, &v) in x.iter_mut().zip(values) {
            *xi = F80::from_f64(v);
        }
        CodeMean { x, diff: [ZERO; 256], stamp: [0; 256], epoch: 0 }
    }

    /// Same arithmetic, same order as `r_mean_x87` over `values[codes[i]]`.
    pub fn mean(&mut self, codes: &[u8]) -> f64 {
        let n = codes.len() as u64;
        if n == 0 {
            return f64::NAN;
        }
        let mut s = ZERO;
        for &c in codes {
            s = s.add(self.x[c as usize]);
        }
        s = s.div_u64(n);
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.stamp = [0; 256];
            self.epoch = 1;
        }
        let mut t = ZERO;
        for &c in codes {
            let c = c as usize;
            if self.stamp[c] != self.epoch {
                self.stamp[c] = self.epoch;
                self.diff[c] = self.x[c].sub(s);
            }
            t = t.add(self.diff[c]);
        }
        s.add(t.div_u64(n)).to_f64()
    }
}
