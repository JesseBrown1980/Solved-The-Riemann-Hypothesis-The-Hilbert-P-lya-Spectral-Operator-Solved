//! The frozen K-family of skew-circulant generators (the experiment builder).
//!
//! Everything here is fixed by `PREREGISTRATION.md` and uses no zeta data.
//!
//! For K = 3^k:
//! - P_K is the smallest prime P ≥ 1,000,000 with P ≡ 1 (mod K);
//! - G_K is the smallest primitive root mod P_K;
//! - W_K = G_K^((P_K − 1)/K), an element of order K;
//! - c[d] = (W^d − W^(−d)) · prime[d] (mod P) for d = 1..(K−1)/2, where prime[d]
//!   is the (d+1)-th prime (prime[0] = 2), and c[K − d] = −c[d];
//! - each c[d] is lifted to its balanced integer representative in (−P/2, P/2];
//! - T[i][j] = c[(j − i) mod K] is a real skew-symmetric integer circulant.
//!
//! At K = 27 this rule gives P = 1,000,081, G = 7, W = 951,846, the constants of
//! the original operator.
//!
//! The ordinates of T are γ_j = 2 Σ_{d=1}^{m} c[d] sin(2π j d / K), j = 1..m,
//! m = (K − 1)/2, and spec(T) = {0, ±iγ_1, …, ±iγ_m}. They are enclosed with
//! integer interval trigonometry ([`crate::interval`]) and cross-checked by exact
//! integer inertia of A = −T² = TᵀT ([`inertia_of_shifted_gram`]).

use crate::interval::{sin_turn, SCALE};

/// Smallest P ≥ `floor` with P prime.
pub const P_FLOOR: u64 = 1_000_000;

pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

#[inline]
pub fn mul_mod_u64(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

pub fn pow_mod_u64(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut acc = 1 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod_u64(acc, b, m);
        }
        b = mul_mod_u64(b, b, m);
        e >>= 1;
    }
    acc
}

fn distinct_prime_factors(mut n: u64) -> Vec<u64> {
    let mut f = Vec::new();
    let mut q = 2;
    while q * q <= n {
        if n % q == 0 {
            f.push(q);
            while n % q == 0 {
                n /= q;
            }
        }
        q += 1;
    }
    if n > 1 {
        f.push(n);
    }
    f
}

/// The first `count` primes.
pub fn first_primes(count: usize) -> Vec<u64> {
    let mut out = Vec::with_capacity(count);
    let mut n = 2u64;
    while out.len() < count {
        if is_prime_u64(n) {
            out.push(n);
        }
        n += 1;
    }
    out
}

/// The frozen operator at size K.
#[derive(Clone, Debug)]
pub struct FrozenOperator {
    pub k: usize,
    pub p: u64,
    pub g: u64,
    pub w: u64,
    /// Balanced integer couplings c[0..K], c[0] = 0, c[K − d] = −c[d].
    pub c: Vec<i64>,
}

impl FrozenOperator {
    pub fn m(&self) -> usize {
        (self.k - 1) / 2
    }
}

/// Build the frozen operator for odd K ≥ 3.
pub fn build(k: usize) -> FrozenOperator {
    assert!(k >= 3 && k % 2 == 1, "K must be odd and at least 3");
    let ku = k as u64;
    let mut p = P_FLOOR;
    while !(p % ku == 1 && is_prime_u64(p)) {
        p += 1;
    }
    let factors = distinct_prime_factors(p - 1);
    let g = (2..p)
        .find(|&g| factors.iter().all(|&q| pow_mod_u64(g, (p - 1) / q, p) != 1))
        .expect("primitive root exists");
    let w = pow_mod_u64(g, (p - 1) / ku, p);
    let primes = first_primes(k);
    let m = (k - 1) / 2;
    let mut c = vec![0i64; k];
    for d in 1..=m {
        let wd = pow_mod_u64(w, d as u64, p);
        let wnd = pow_mod_u64(w, (k - d) as u64, p);
        let diff = (wd + p - wnd) % p;
        let r = mul_mod_u64(diff, primes[d] % p, p);
        let lifted = if r > p / 2 {
            r as i64 - p as i64
        } else {
            r as i64
        };
        c[d] = lifted;
        c[k - d] = -lifted;
    }
    FrozenOperator { k, p, g, w, c }
}

/// Signed enclosure [lo, hi] of γ_j at scale [`SCALE`], for j = 1..=m.
pub fn ordinate_enclosure(op: &FrozenOperator, j: usize) -> (i128, i128) {
    let k = op.k;
    let (mut lo, mut hi) = (0i128, 0i128);
    for d in 1..=op.m() {
        let (sl, sh) = sin_turn(((j * d) % k) as u64, k as u64);
        let c = op.c[d] as i128;
        if c >= 0 {
            lo += c * sl;
            hi += c * sh;
        } else {
            lo += c * sh;
            hi += c * sl;
        }
    }
    (2 * lo, 2 * hi)
}

/// One positive ordinate |γ_j| with its enclosure.
#[derive(Clone, Copy, Debug)]
pub struct Ordinate {
    pub j: usize,
    pub lo: i128,
    pub hi: i128,
}

/// All |γ_j| enclosures sorted by lower end. `Err` if an enclosure contains 0.
pub fn sorted_positive_ordinates(op: &FrozenOperator) -> Result<Vec<Ordinate>, String> {
    let mut out = Vec::with_capacity(op.m());
    for j in 1..=op.m() {
        let (lo, hi) = ordinate_enclosure(op, j);
        let (alo, ahi) = if lo > 0 {
            (lo, hi)
        } else if hi < 0 {
            (-hi, -lo)
        } else {
            return Err(format!("gamma_{j} enclosure [{lo},{hi}] contains 0"));
        };
        out.push(Ordinate {
            j,
            lo: alo,
            hi: ahi,
        });
    }
    out.sort_by_key(|o| (o.lo, o.j));
    Ok(out)
}

/// First row of A = TᵀT: a[s] = Σ_e c[e]·c[(e + s) mod K].
pub fn gram_first_row(op: &FrozenOperator) -> Vec<i64> {
    let k = op.k;
    (0..k)
        .map(|s| {
            (0..k)
                .map(|e| op.c[e] as i128 * op.c[(e + s) % k] as i128)
                .sum::<i128>() as i64
        })
        .collect()
}

/// Moduli for the multi-modular determinant signs: the largest primes below 2^61.
fn moduli(count: usize, skip: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(count);
    let mut n: u64 = (1u64 << 61) - 1;
    while out.len() < count {
        if is_prime_u64(n) && !skip.contains(&n) {
            out.push(n);
        }
        n -= 2;
    }
    out
}

/// Exact number of negative eigenvalues of A − τI, A = TᵀT.
///
/// Counts sign changes in the leading principal minors D_0 = 1, D_1, …, D_K of
/// the integer matrix A − τI (Jacobi/Sylvester). Each D_k is computed modulo
/// enough 61-bit primes to exceed twice its Hadamard bound, and its sign is read
/// from the balanced mixed-radix (Garner) digits: the sign of the last nonzero
/// digit. Returns `Err` if some D_k is exactly 0.
pub fn inertia_of_shifted_gram(op: &FrozenOperator, tau: i64) -> Result<usize, String> {
    let k = op.k;
    let a = gram_first_row(op);
    let entry = |i: usize, j: usize| -> i64 {
        let v = a[(j + k - i) % k];
        if i == j {
            v - tau
        } else {
            v
        }
    };
    // Hadamard: |D_n| ≤ Π_{i<n} ||row_i||; all rows of a circulant shift have the
    // same norm, so one bound covers every leading minor.
    let mut norm2: u128 = 0;
    for j in 0..k {
        let v = entry(0, j).unsigned_abs() as u128;
        norm2 += v * v;
    }
    let row_bits = (128 - norm2.leading_zeros()).div_ceil(2) as usize; // ||row|| < 2^row_bits
    let total_bits = row_bits * k + 2;
    let need = total_bits.div_ceil(60) + 1;

    let mut residues: Vec<Vec<u64>> = Vec::with_capacity(need); // [prime][minor index 1..=k]
    let mut used: Vec<u64> = Vec::with_capacity(need);
    let mut rejected: Vec<u64> = Vec::new();
    while used.len() < need {
        let batch = moduli(
            need - used.len(),
            &[used.clone(), rejected.clone()].concat(),
        );
        for p in batch {
            match leading_minors_mod(k, &entry, p) {
                Some(r) => {
                    residues.push(r);
                    used.push(p);
                }
                None => {
                    rejected.push(p);
                    if rejected.len() > 16 {
                        return Err("too many zero pivots: a leading minor is 0".into());
                    }
                }
            }
        }
    }
    let signs = garner_signs(&used, &residues, k);
    let mut sign_prev = 1i8;
    let mut neg = 0usize;
    for (idx, &s) in signs.iter().enumerate() {
        if s == 0 {
            return Err(format!("leading minor D_{} is 0", idx + 1));
        }
        if s != sign_prev {
            neg += 1;
        }
        sign_prev = s;
    }
    Ok(neg)
}

/// D_1..D_K mod p by elimination without pivoting; `None` on a zero pivot.
fn leading_minors_mod(k: usize, entry: &dyn Fn(usize, usize) -> i64, p: u64) -> Option<Vec<u64>> {
    let red = |v: i64| -> u64 { v.rem_euclid(p as i64) as u64 };
    // Symmetric: keep the upper triangle only, row-major.
    let mut m: Vec<Vec<u64>> = (0..k)
        .map(|i| (i..k).map(|j| red(entry(i, j))).collect())
        .collect();
    let mut out = Vec::with_capacity(k);
    let mut det = 1u64;
    for piv in 0..k {
        let pv = m[piv][0];
        if pv == 0 {
            return None;
        }
        det = mul_mod_u64(det, pv, p);
        out.push(det);
        let inv = pow_mod_u64(pv, p - 2, p);
        let (head, tail) = m.split_at_mut(piv + 1);
        let prow = &head[piv];
        for (off, row) in tail.iter_mut().enumerate() {
            let i = piv + 1 + off;
            // row i (upper part starts at column i); multiplier from column i of pivot row
            let f = mul_mod_u64(prow[i - piv], inv, p);
            if f == 0 {
                continue;
            }
            for (jj, cell) in row.iter_mut().enumerate() {
                let j = i + jj;
                let sub = mul_mod_u64(f, prow[j - piv], p);
                *cell = if *cell >= sub {
                    *cell - sub
                } else {
                    *cell + p - sub
                };
            }
        }
        m[piv].clear();
        // shift: next pivot row's first stored entry is its diagonal already
    }
    Some(out)
}

/// Signs of the integers with residues `res[i][idx]` mod `mods[i]`, for every
/// idx < `n`, via balanced mixed-radix (Garner) digits: with every digit in
/// (−m_i/2, m_i/2], the sign of the value is the sign of its last nonzero digit.
fn garner_signs(mods: &[u64], res: &[Vec<u64>], n: usize) -> Vec<i8> {
    let r = mods.len();
    // radix[i][t] = Π_{s<t} m_s (mod m_i), t ≤ i; inv[i] = radix[i][i]^(-1) mod m_i
    let mut radix: Vec<Vec<u64>> = Vec::with_capacity(r);
    let mut inv: Vec<u64> = Vec::with_capacity(r);
    for i in 0..r {
        let mi = mods[i];
        let mut row = Vec::with_capacity(i + 1);
        let mut acc = 1u64;
        for &mt in mods.iter().take(i) {
            row.push(acc);
            acc = mul_mod_u64(acc, mt % mi, mi);
        }
        row.push(acc);
        inv.push(pow_mod_u64(acc, mi - 2, mi));
        radix.push(row);
    }
    let mut out = Vec::with_capacity(n);
    let mut digits: Vec<i128> = vec![0; r];
    for idx in 0..n {
        for i in 0..r {
            let mi = mods[i];
            let mut acc: u64 = 0;
            for (dig, &rad) in digits.iter().zip(radix[i].iter()).take(i) {
                let dt = dig.rem_euclid(mi as i128) as u64;
                acc = (acc + mul_mod_u64(dt, rad, mi)) % mi;
            }
            let x = res[i][idx] % mi;
            let v = mul_mod_u64((x + mi - acc) % mi, inv[i], mi);
            digits[i] = if v > mi / 2 {
                v as i128 - mi as i128
            } else {
                v as i128
            };
        }
        out.push(match digits.iter().rev().find(|&&d| d != 0) {
            None => 0,
            Some(&d) if d > 0 => 1,
            _ => -1,
        });
    }
    out
}

/// Exact-integer bound: ceil(x / SCALE) for x ≥ 0.
pub fn ceil_units(x: i128) -> i128 {
    (x + SCALE - 1).div_euclid(SCALE)
}

/// Exact-integer bound: floor(x / SCALE).
pub fn floor_units(x: i128) -> i128 {
    x.div_euclid(SCALE)
}
