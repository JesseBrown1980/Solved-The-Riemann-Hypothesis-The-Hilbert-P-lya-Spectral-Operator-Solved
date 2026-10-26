//! Relic landings: the nested operator (a network inside a network).
//!
//! Fixed by `PREREGISTRATION-RELIC.md`. Reads no zeta data.
//!
//! Landing sizes follow THE RELIC REDISCOVERY: A_0 = 27, A_(n+1) = A_n · 3^(n+2),
//! so 27 → 243 → 6,561. The operator at each landing nests the previous one:
//!
//! ```text
//! T^(0)   = the frozen K = 27 circulant              (crate::family::build(27))
//! T^(n+1) = T^(n) ⊗ I_q + I_(A_n) ⊗ S_q,  q = 3^(n+2)
//! ```
//!
//! S_q is the family circulant of size q (same P_q, G_q, W_q rule) whose
//! couplings take fresh primes: c[d] = (W^d − W^(−d)) · prime[offset + d], where
//! offset counts the couplings of every earlier factor (13 for T^(0), then 4 for
//! S_9, …). No prime is reused, so the inner wave is never a copy of the outer
//! network. Each outer node carries a copy of the inner wave network S_q. Both terms are
//! real skew-symmetric integer matrices, so T^(n+1) is too. Its spectrum is the
//! Kronecker sum {i(a + b) : ia ∈ spec T^(n), ib ∈ spec S_q}, which keeps exactly
//! one zero mode, the shared free centre C = 1 (0 + 0), when no other sum vanishes.

use crate::family::{
    build, first_primes, garner_signs, leading_minors_mod, moduli, mul_mod_u64, ordinate_enclosure,
    pow_mod_u64, FrozenOperator,
};

/// Family circulant of size `k` with couplings on primes offset by `offset`.
pub fn build_with_prime_offset(k: usize, offset: usize) -> FrozenOperator {
    let base = build(k);
    let (p, w) = (base.p, base.w);
    let primes = first_primes(offset + k);
    let m = (k - 1) / 2;
    let mut c = vec![0i64; k];
    for d in 1..=m {
        let wd = pow_mod_u64(w, d as u64, p);
        let wnd = pow_mod_u64(w, (k - d) as u64, p);
        let r = mul_mod_u64((wd + p - wnd) % p, primes[offset + d] % p, p);
        let lifted = if r > p / 2 {
            r as i64 - p as i64
        } else {
            r as i64
        };
        c[d] = lifted;
        c[k - d] = -lifted;
    }
    FrozenOperator { c, ..base }
}

/// The factors of landing `level`: (size, prime offset), outer first.
pub fn factors(level: usize) -> Vec<(usize, usize)> {
    let mut out = vec![(27usize, 0usize)];
    let mut offset = 13usize;
    for n in 0..level {
        let q = 3usize.pow(n as u32 + 2);
        out.push((q, offset));
        offset += (q - 1) / 2;
    }
    out
}

/// All L signed ordinates of a landing as enclosures at scale 10^18; index 0 is
/// the exact zero mode.
#[derive(Clone, Debug)]
pub struct Landing {
    pub level: usize,
    pub size: usize,
    /// (lo, hi, exact_zero)
    pub signed: Vec<(i128, i128, bool)>,
}

/// Signed spectrum of a frozen circulant: γ_0 = 0, γ_j, γ_(K−j) = −γ_j.
pub fn signed_circulant(op: &FrozenOperator) -> Vec<(i128, i128, bool)> {
    let k = op.k;
    let mut out = vec![(0i128, 0i128, true); k];
    for j in 1..=op.m() {
        let (lo, hi) = ordinate_enclosure(op, j);
        out[j] = (lo, hi, false);
        out[k - j] = (-hi, -lo, false);
    }
    out
}

/// Size of landing n: 27, 243, 6561, …
pub fn landing_size(level: usize) -> usize {
    let mut a = 27usize;
    for n in 0..level {
        a *= 3usize.pow(n as u32 + 2);
    }
    a
}

/// The landing at `level`, built by nesting.
pub fn landing(level: usize) -> Landing {
    let mut cur = Landing {
        level: 0,
        size: 27,
        signed: signed_circulant(&build(27)),
    };
    for (n, &(q, offset)) in factors(level).iter().skip(1).enumerate() {
        let inner = signed_circulant(&build_with_prime_offset(q, offset));
        let mut next = Vec::with_capacity(cur.size * q);
        for &(al, ah, az) in &cur.signed {
            for &(bl, bh, bz) in &inner {
                next.push((al + bl, ah + bh, az && bz));
            }
        }
        cur = Landing {
            level: n + 1,
            size: cur.size * q,
            signed: next,
        };
    }
    cur
}

/// Positive ordinates sorted by lower end, as (lo, hi, index). `Err` if an
/// enclosure other than the free centre contains 0.
pub fn positive_ordinates(l: &Landing) -> Result<Vec<(i128, i128, usize)>, String> {
    let mut pos = Vec::new();
    let mut zeros = 0usize;
    for (idx, &(lo, hi, z)) in l.signed.iter().enumerate() {
        if z {
            zeros += 1;
        } else if lo > 0 {
            pos.push((lo, hi, idx));
        } else if hi >= 0 {
            return Err(format!("ordinate {idx} enclosure [{lo},{hi}] contains 0"));
        }
    }
    if zeros != 1 {
        return Err(format!("{zeros} exact zero modes"));
    }
    pos.sort();
    Ok(pos)
}

/// Dense integer generator of a landing (row-major), for the inertia check.
pub fn dense_generator(level: usize) -> (usize, Vec<i64>) {
    let base = build(27);
    let mut n = 27usize;
    let mut t: Vec<i64> = (0..n * n)
        .map(|e| base.c[(e % n + n - e / n) % n])
        .collect();
    for &(q, offset) in factors(level).iter().skip(1) {
        let s = build_with_prime_offset(q, offset);
        let m = n * q;
        let mut u = vec![0i64; m * m];
        for x in 0..n {
            for y in 0..q {
                let row = x * q + y;
                for x2 in 0..n {
                    let v = t[x * n + x2];
                    if v != 0 {
                        u[row * m + x2 * q + y] += v;
                    }
                }
                for y2 in 0..q {
                    u[row * m + x * q + y2] += s.c[(y2 + q - y) % q];
                }
            }
        }
        t = u;
        n = m;
    }
    (n, t)
}

/// Exact number of negative eigenvalues of TᵀT − τI for a dense integer T.
pub fn inertia_dense_gram(n: usize, t: &[i64], tau: i64) -> Result<usize, String> {
    let mut a = vec![0i64; n * n];
    for i in 0..n {
        for j in i..n {
            let mut s: i128 = 0;
            for l in 0..n {
                s += t[l * n + i] as i128 * t[l * n + j] as i128;
            }
            let v = i64::try_from(s).map_err(|_| "gram entry overflow".to_string())?;
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
        a[i * n + i] -= tau;
    }
    let entry = |i: usize, j: usize| -> i64 { a[i * n + j] };
    let mut row_bits = 0usize;
    for i in 0..n {
        let mut norm2: u128 = 0;
        for j in 0..n {
            let v = entry(i, j).unsigned_abs() as u128;
            norm2 += v * v;
        }
        row_bits = row_bits.max((128 - norm2.leading_zeros()).div_ceil(2) as usize);
    }
    let need = (row_bits * n + 2).div_ceil(60) + 1;
    let mut residues: Vec<Vec<u64>> = Vec::with_capacity(need);
    let mut used: Vec<u64> = Vec::with_capacity(need);
    let mut rejected: Vec<u64> = Vec::new();
    while used.len() < need {
        for p in moduli(
            need - used.len(),
            &[used.clone(), rejected.clone()].concat(),
        ) {
            match leading_minors_mod(n, &entry, p) {
                Some(r) => {
                    residues.push(r);
                    used.push(p);
                }
                None => {
                    rejected.push(p);
                    if rejected.len() > 16 {
                        return Err("too many zero pivots".into());
                    }
                }
            }
        }
    }
    let mut prev = 1i8;
    let mut neg = 0usize;
    for (idx, s) in garner_signs(&used, &residues, n).into_iter().enumerate() {
        if s == 0 {
            return Err(format!("leading minor D_{} is 0", idx + 1));
        }
        neg += (s != prev) as usize;
        prev = s;
    }
    Ok(neg)
}
