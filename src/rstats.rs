//! Replicas of R's `mean()` and `median()` numerics.
//!
//! R computes `mean` in two passes with a LDOUBLE accumulator (80-bit on x86
//! Linux, 64-bit on arm64 macOS). We use f64, which is bit-identical to R on
//! arm64 and within ~1e-15 relative of R on x86.

pub fn r_mean(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mut s = 0.0f64;
    for &v in x {
        s += v;
    }
    s /= n;
    if s.is_finite() {
        let mut t = 0.0f64;
        for &v in x {
            t += v - s;
        }
        s += t / n;
    }
    s
}

/// `median(x)` for non-empty, NaN-free input. Sorts `x` in place.
pub fn r_median(x: &mut [f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    x.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let half = n / 2;
    if n % 2 == 1 { x[half] } else { r_mean(&[x[half - 1], x[half]]) }
}

/// Median of values `code/255` given a histogram of the 256 ML codes.
pub fn r_median_ml(hist: &[u32; 256], n: usize) -> f64 {
    let nth = |k: usize| -> f64 {
        // k-th smallest (0-based)
        let mut acc = 0usize;
        for (code, &c) in hist.iter().enumerate() {
            acc += c as usize;
            if acc > k {
                return code as f64 / 255.0;
            }
        }
        unreachable!()
    };
    let half = n / 2;
    if n % 2 == 1 { nth(half) } else { r_mean(&[nth(half - 1), nth(half)]) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_matches_r() {
        assert_eq!(r_median(&mut [3.0, 1.0, 2.0]), 2.0);
        assert_eq!(r_median(&mut [4.0, 1.0, 2.0, 3.0]), 2.5);
        let mut hist = [0u32; 256];
        hist[0] = 2;
        hist[255] = 2;
        assert_eq!(r_median_ml(&hist, 4), 0.5);
        hist[255] = 3;
        assert_eq!(r_median_ml(&hist, 5), 1.0);
    }

    #[test]
    fn mean_two_pass() {
        assert_eq!(r_mean(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert!(r_mean(&[]).is_nan());
    }
}
