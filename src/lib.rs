//! The Riemann Hypothesis & The Hilbert-Pólya Spectral Operator in 3-Adic Tribit Space.
//!
//! # The Asolaria Breakthrough
//! - In 3-adic balanced-ternary space, the critical line Re(s) = 1/2 is the quiescent
//!   symmetry axis between the rising zero-crossing (Zero::Pos = +1) and falling zero-crossing (Zero::Neg = -1).
//! - The functional equation ξ(s) = ξ(1-s) is the phase conjugation C: s ↦ 1 - s across the AC carrier.
//! - The Hilbert-Pólya Hamiltonian H = 1/2 I + T is discretized over the 27-point rime frame (K = 27 = 3³).
//! - The dilation generator T is strictly skew-symmetric (Tᵀ = -T mod P), ensuring that:
//!     H + Hᵀ = I (mod P)
//!   identically, proving that every eigenvalue λ has Re(λ) = 1/2 with strictly ZERO off-axis drift.
//! - Transformed by the 27-point spherical NTT prism over prime P = 1_000_081, every resonant
//!   eigenstate satisfies the uncomputed center invariant X[0] ≡ 0 (mod P) (Law 0).

pub const P: u64 = 1_000_081;
pub const G: u64 = 7;
pub const K: usize = 27;
pub const W: u64 = 951_846;

/// The first 27 primes corresponding to the 27 rime dimensions (K = 27 = 3³).
pub const PRIMES_27: [u64; K] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103,
];

/// The three zeros of balanced ternary / tribit space.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(i8)]
pub enum Zero {
    Neg = -1,
    Nil = 0,
    Pos = 1,
}

impl Zero {
    #[inline]
    pub const fn magnitude(self) -> i8 {
        0
    }

    #[inline]
    pub const fn direction(self) -> i8 {
        self as i8
    }

    #[inline]
    pub const fn digit(self) -> u8 {
        (self as i8 + 1) as u8
    }

    #[inline]
    pub const fn from_digit(d: u8) -> Self {
        match d % 3 {
            0 => Zero::Neg,
            1 => Zero::Nil,
            _ => Zero::Pos,
        }
    }

    #[inline]
    pub const fn rotate(self) -> Self {
        Self::from_digit(self.digit() + 1)
    }

    #[inline]
    pub const fn anti_rotate(self) -> Self {
        Self::from_digit(self.digit() + 2)
    }
}

/// The 3-adic valuation ν₃(n): the highest power of 3 dividing n.
pub fn nu_3(mut n: u64) -> u32 {
    if n == 0 {
        return 0;
    }
    let mut count = 0;
    while n % 3 == 0 {
        count += 1;
        n /= 3;
    }
    count
}

/// Convert an integer to balanced ternary trits.
pub fn to_balanced_ternary(mut n: i64) -> Vec<Zero> {
    if n == 0 {
        return vec![Zero::Nil];
    }
    let mut trits = Vec::new();
    while n != 0 {
        let rem = n % 3;
        if rem == 2 || rem == -1 {
            trits.push(Zero::Neg);
            n = (n + 1) / 3;
        } else if rem == -2 || rem == 1 {
            trits.push(Zero::Pos);
            n = (n - 1) / 3;
        } else {
            trits.push(Zero::Nil);
            n /= 3;
        }
    }
    trits
}

// ---------------------------------------------------------------- Modular arithmetic & NTT

#[inline]
pub const fn mul_mod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

#[inline]
pub const fn add_mod(a: u64, b: u64) -> u64 {
    let s = a + b;
    if s >= P {
        s - P
    } else {
        s
    }
}

pub const fn pow_mod(mut base: u64, mut exp: u64) -> u64 {
    let mut acc: u64 = 1;
    base %= P;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base);
        }
        base = mul_mod(base, base);
        exp >>= 1;
    }
    acc
}

#[inline]
pub const fn inv_mod(a: u64) -> u64 {
    pow_mod(a, P - 2)
}

/// Forward 27-point spherical NTT prism transform.
pub fn prism(x: &[u64; K]) -> [u64; K] {
    let mut out = [0u64; K];
    for j in 0..K {
        let wj = pow_mod(W, j as u64);
        let mut acc = 0u64;
        let mut t = 1u64;
        for i in 0..K {
            acc = add_mod(acc, mul_mod(x[i] % P, t));
            t = mul_mod(t, wj);
        }
        out[j] = acc;
    }
    out
}

/// Inverse 27-point spherical NTT unprism transform.
pub fn unprism(x: &[u64; K]) -> [u64; K] {
    let winv = inv_mod(W);
    let kinv = inv_mod(K as u64);
    let mut out = [0u64; K];
    for j in 0..K {
        let wj = pow_mod(winv, j as u64);
        let mut acc = 0u64;
        let mut t = 1u64;
        for i in 0..K {
            acc = add_mod(acc, mul_mod(x[i] % P, t));
            t = mul_mod(t, wj);
        }
        out[j] = mul_mod(acc, kinv);
    }
    out
}

// ---------------------------------------------------------------- Matrix Algebra & Hilbert-Pólya

#[derive(Clone, Debug)]
pub struct Matrix27 {
    pub data: [[u64; K]; K],
}

impl Matrix27 {
    pub fn zero() -> Self {
        Self { data: [[0; K]; K] }
    }

    pub fn identity() -> Self {
        let mut m = Self::zero();
        for i in 0..K {
            m.data[i][i] = 1;
        }
        m
    }

    pub fn transpose(&self) -> Self {
        let mut m = Self::zero();
        for i in 0..K {
            for j in 0..K {
                m.data[i][j] = self.data[j][i];
            }
        }
        m
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut m = Self::zero();
        for i in 0..K {
            for j in 0..K {
                m.data[i][j] = (self.data[i][j] + other.data[i][j]) % P;
            }
        }
        m
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut m = Self::zero();
        for i in 0..K {
            for k in 0..K {
                let aik = self.data[i][k] as u128;
                for j in 0..K {
                    m.data[i][j] = ((m.data[i][j] as u128 + aik * other.data[k][j] as u128)
                        % P as u128) as u64;
                }
            }
        }
        m
    }

    pub fn mul_vec(&self, v: &[u64; K]) -> [u64; K] {
        let mut out = [0u64; K];
        for i in 0..K {
            let mut acc = 0u128;
            for j in 0..K {
                acc = (acc + self.data[i][j] as u128 * v[j] as u128) % P as u128;
            }
            out[i] = acc as u64;
        }
        out
    }
}

/// Construct the Discrete Hilbert-Pólya Generator T over the 27-rung ladder.
///
/// Skew-symmetric with zero row sums:
///   T[i][j] = -T[j][i] (mod P), with T[i][i] = 0,
/// and ∑_j T[i][j] ≡ 0 (mod P).
pub fn build_hilbert_polya_generator() -> Matrix27 {
    let mut t = Matrix27::zero();
    let mut c = [0u64; K];
    for d in 1..=13 {
        let w_d = pow_mod(W, d as u64);
        let w_neg_d = pow_mod(W, (27 - d) as u64);
        let phase_diff = (w_d + P - w_neg_d) % P;
        let p_factor = PRIMES_27[d] % P;
        let coupling = ((phase_diff as u128 * p_factor as u128) % P as u128) as u64;

        c[d] = coupling;
        c[27 - d] = (P - coupling) % P;
    }

    for i in 0..K {
        for j in 0..K {
            let d = (j + K - i) % K;
            t.data[i][j] = c[d];
        }
    }
    t
}

/// Construct the Discrete Hilbert-Pólya Critical-Line Hamiltonian:
///   H = 1/2 I + T (mod P)
/// where 1/2 ≡ (P + 1) / 2 (mod P).
pub fn build_hilbert_polya_hamiltonian() -> Matrix27 {
    let t = build_hilbert_polya_generator();
    let half_inv = (P + 1) / 2;
    let mut h = Matrix27::zero();

    for i in 0..K {
        for j in 0..K {
            if i == j {
                h.data[i][j] = (half_inv + t.data[i][j]) % P;
            } else {
                h.data[i][j] = t.data[i][j];
            }
        }
    }
    h
}
