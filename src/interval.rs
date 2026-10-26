//! Integer interval trigonometry at fixed scale 10^18.
//!
//! Every value is an exact integer enclosure [lo, hi] of x · 10^18. Rounding is
//! outward at every step (floor for lower ends, ceil for upper ends), and the
//! Taylor tail is added as an explicit bound.

/// Fixed-point scale.
pub const SCALE: i128 = 1_000_000_000_000_000_000;

/// π · 10^18 enclosure: π = 3.14159265358979323846…
pub const PI_LO: i128 = 3_141_592_653_589_793_238;
pub const PI_HI: i128 = 3_141_592_653_589_793_239;

const TERMS: u32 = 16;

#[inline]
fn div_floor(a: i128, b: i128) -> i128 {
    a.div_euclid(b)
}

#[inline]
fn div_ceil(a: i128, b: i128) -> i128 {
    -((-a).div_euclid(b))
}

/// Enclosure of sin x or cos x for x ∈ [xl, xh] ⊂ [0, 0.8], scaled.
fn taylor(xl: i128, xh: i128, cosine: bool) -> (i128, i128) {
    // first term: cos -> 1, sin -> x
    let (mut tl, mut th) = if cosine { (SCALE, SCALE) } else { (xl, xh) };
    let mut n: i128 = if cosine { 0 } else { 1 };
    let (mut lo, mut hi) = (tl, th);
    for step in 1..=TERMS {
        // next term magnitude: t · x² / ((n+1)(n+2))
        let den = (n + 1) * (n + 2);
        tl = div_floor(div_floor(div_floor(tl * xl, SCALE) * xl, SCALE), den);
        th = div_ceil(div_ceil(div_ceil(th * xh, SCALE) * xh, SCALE), den);
        n += 2;
        if step == TERMS {
            // alternating series with decreasing terms: tail bounded by this term
            lo -= th;
            hi += th;
        } else if step % 2 == 1 {
            lo -= th;
            hi -= tl;
        } else {
            lo += tl;
            hi += th;
        }
    }
    (lo.max(-SCALE), hi.min(SCALE))
}

/// Enclosure of sin(2π n / k) · 10^18, for n < k.
pub fn sin_turn(n: u64, k: u64) -> (i128, i128) {
    let k = k as i128;
    let mut a = 8 * (n as i128 % k); // angle = 2π a / (8k)
    let mut neg = false;
    if a >= 4 * k {
        neg = true;
        a -= 4 * k;
    }
    if a > 2 * k {
        a = 4 * k - a;
    }
    let (b, cosine) = if a > k { (2 * k - a, true) } else { (a, false) };
    // x = π b / (4k) ∈ [0, π/4]
    let xl = div_floor(PI_LO * b, 4 * k);
    let xh = div_ceil(PI_HI * b, 4 * k);
    let (lo, hi) = if b == 0 {
        if cosine {
            (SCALE, SCALE)
        } else {
            (0, 0)
        }
    } else {
        taylor(xl, xh, cosine)
    };
    if neg {
        (-hi, -lo)
    } else {
        (lo, hi)
    }
}
