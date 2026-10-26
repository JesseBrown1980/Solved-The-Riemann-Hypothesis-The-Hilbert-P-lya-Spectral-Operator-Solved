//! Experiment comparator (PREREGISTRATION.md §4–§5).
//!
//! Reads a builder receipt and the zeta-zero table and writes the comparison:
//!
//! ```text
//! cargo run --release --bin compare_zeta -- receipts/experiment/SPECTRUM-K27.hbp data/zeta/ZEROS-1000.hbp
//! ```
//!
//! Self-contained: it does not link the builder modules. All arithmetic is
//! exact integer arithmetic with 256-bit products.

use std::fs;

const SCALE: u128 = 1_000_000_000_000_000_000;

fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split_whitespace()
        .find_map(|t| t.strip_prefix(key).and_then(|r| r.strip_prefix('=')))
}

fn num(line: &str, key: &str) -> u128 {
    field(line, key)
        .unwrap_or_else(|| panic!("missing {key} in: {line}"))
        .parse()
        .unwrap_or_else(|_| panic!("bad {key} in: {line}"))
}

/// Full 256-bit product of two u128 values as (hi, lo).
fn mul256(a: u128, b: u128) -> (u128, u128) {
    let mask = u64::MAX as u128;
    let (a1, a0) = (a >> 64, a & mask);
    let (b1, b0) = (b >> 64, b & mask);
    let p00 = a0 * b0;
    let p01 = a0 * b1;
    let p10 = a1 * b0;
    let p11 = a1 * b1;
    let mid = (p00 >> 64) + (p01 & mask) + (p10 & mask);
    let lo = (p00 & mask) | ((mid & mask) << 64);
    let hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    (hi, lo)
}

/// a·b ≤ c·d, exactly.
fn le_prod(a: u128, b: u128, c: u128, d: u128) -> bool {
    mul256(a, b) <= mul256(c, d)
}

/// floor(a·b / c) for the descriptive table (values small enough for u128 after division).
fn muldiv_floor(a: u128, b: u128, c: u128) -> u128 {
    // long division of the 256-bit product by c, bit by bit
    let (hi, lo) = mul256(a, b);
    let mut rem: u128 = 0;
    let mut q: u128 = 0;
    for i in (0..256).rev() {
        let bit = if i >= 128 {
            (hi >> (i - 128)) & 1
        } else {
            (lo >> i) & 1
        };
        let carry = rem >> 127;
        rem = (rem << 1) | bit;
        q <<= 1;
        if carry == 1 || rem >= c {
            rem = rem.wrapping_sub(c);
            q |= 1;
        }
    }
    q
}

fn fmt_scaled(x: u128) -> String {
    format!("{}.{:06}", x / SCALE, (x % SCALE) / 1_000_000_000_000)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let spec = fs::read_to_string(&args[1]).expect("spectrum receipt");
    let zeros = fs::read_to_string(&args[2]).expect("zero table");

    let header = spec.lines().find(|l| l.starts_with("K=")).expect("K line");
    let k = num(header, "K");
    let m = num(header, "M") as usize;
    let mut ord: Vec<(u128, u128)> = vec![(0, 0); m];
    for l in spec.lines().filter(|l| l.starts_with("ORD ")) {
        let r = num(l, "rank") as usize;
        ord[r - 1] = (num(l, "lo"), num(l, "hi"));
    }
    assert!(
        ord.iter().all(|&(lo, hi)| 0 < lo && lo <= hi),
        "ordinate enclosures"
    );
    let mut z: Vec<(u128, u128)> = Vec::new();
    for l in zeros.lines().filter(|l| l.starts_with("Z ")) {
        assert_eq!(num(l, "n") as usize, z.len() + 1, "zeros in order");
        z.push((num(l, "lo"), num(l, "hi")));
    }
    assert!(z.len() > m, "need zeros 1..m+1");

    let (g1l, g1h) = ord[0];
    let (z1l, z1h) = z[0];
    println!("HBP ZETA-COMPARISON v1");
    println!(
        "K={k} M={m} anchor=rank1_to_zeta1 zeros_available={}",
        z.len()
    );

    let (mut agree, mut disagree, mut ambiguous) = (0usize, 0usize, 0usize);
    let mut first_disagree: Option<usize> = None;
    let mut max_dev: u128 = 0;
    for kk in 1..=m {
        let (_, uk) = z[kk - 1];
        let (lk1, _) = z[kk];
        let t = (uk + lk1) / 2;
        assert!(uk < t && t < lk1, "threshold {kk} separates zeros");
        let (mut n_in, mut n_amb) = (0usize, 0usize);
        for (n, &(gl, gh)) in ord.iter().enumerate() {
            if n == 0 || le_prod(gh, z1h, t, g1l) {
                n_in += 1;
            } else if !le_prod(gl, z1l, t, g1h) {
                // definitely out: gl·ℓ1 > t·gh1
            } else {
                n_amb += 1;
            }
        }
        let status = if n_amb == 0 && n_in == kk {
            agree += 1;
            "AGREE"
        } else if kk < n_in || kk > n_in + n_amb {
            disagree += 1;
            first_disagree.get_or_insert(kk);
            "DISAGREE"
        } else {
            ambiguous += 1;
            "AMBIGUOUS"
        };
        if n_amb == 0 {
            max_dev = max_dev.max((n_in as i128 - kk as i128).unsigned_abs());
        }
        println!(
            "T k={kk} t={t} N_zeta={kk} N_op_in={n_in} N_op_ambiguous={n_amb} status={status}"
        );
    }
    // descriptive: anchored ordinates γ̂_n = |γ|_(n) · γ1 / |γ|_(1), midpoints
    let z1 = (z1l + z1h) / 2;
    let g1 = (g1l + g1h) / 2;
    let show: Vec<usize> = (1..=m.min(10)).chain(std::iter::once(m)).collect();
    for n in show {
        let gn = (ord[n - 1].0 + ord[n - 1].1) / 2;
        let hat = muldiv_floor(gn, z1, g1);
        let zn = (z[n - 1].0 + z[n - 1].1) / 2;
        println!(
            "ANCHORED n={n} gamma_hat={} zeta={}",
            fmt_scaled(hat),
            fmt_scaled(zn)
        );
    }
    let verdict = if agree == m {
        "MATCH"
    } else if disagree > 0 {
        "MISMATCH"
    } else {
        "INCONCLUSIVE"
    };
    println!(
        "VERDICT K={k} verdict={verdict} agree={agree}/{m} disagree={disagree} ambiguous={ambiguous} first_disagree={} max_abs_dev={max_dev}",
        first_disagree.map_or("-".to_string(), |x| x.to_string())
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul256_matches_small_products() {
        assert_eq!(mul256(3, 5), (0, 15));
        let (hi, lo) = mul256(u128::MAX, u128::MAX);
        assert_eq!((hi, lo), (u128::MAX - 1, 1));
        assert_eq!(
            muldiv_floor(10u128.pow(30), 10u128.pow(30), 10u128.pow(24)),
            10u128.pow(36)
        );
    }
}
