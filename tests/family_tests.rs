use riemann_hilbert_polya_spectral::family::{
    build, gram_first_row, inertia_of_shifted_gram, ordinate_enclosure, sorted_positive_ordinates,
};
use riemann_hilbert_polya_spectral::interval::{sin_turn, SCALE};
use riemann_hilbert_polya_spectral::{generator_couplings_integer, G, P, W};

#[test]
fn k27_rule_reproduces_the_original_operator() {
    let op = build(27);
    assert_eq!((op.p, op.g, op.w), (P, G, W));
    assert_eq!(op.c.as_slice(), generator_couplings_integer().as_slice());
}

#[test]
fn k81_rule_reproduces_the_clasy_x3_field() {
    let op = build(81);
    assert_eq!((op.p, op.w), (1_000_999, 30_722));
}

#[test]
fn sin_enclosures_are_tight_and_consistent() {
    for k in [27u64, 81, 243, 729] {
        for n in 0..k {
            let (lo, hi) = sin_turn(n, k);
            assert!(lo <= hi && hi - lo < 200, "width at {n}/{k}: {}", hi - lo);
            // sin(2π(k−n)/k) = −sin(2πn/k)
            let (l2, h2) = sin_turn((k - n) % k, k);
            assert!(l2 <= -lo && -hi <= h2);
            // sin² + cos² = 1, cos(θ) = sin(θ + π/2) not on the grid for odd k; check |sin| ≤ 1
            assert!(lo >= -SCALE && hi <= SCALE);
        }
    }
    // sin(2π/4) is not on an odd grid; check sin(2π·1/12) = 1/2 on k = 12
    let (lo, hi) = sin_turn(1, 12);
    assert!(lo <= SCALE / 2 && SCALE / 2 <= hi);
    let (lo, hi) = sin_turn(1, 8); // √2/2 = 0.70710678118654752440…
    assert!(lo <= 707_106_781_186_547_524 && 707_106_781_186_547_525 <= hi);
}

#[test]
fn trace_of_gram_matches_sum_of_squared_ordinates() {
    // Tr(TᵀT) = K·a[0] = 2 Σ γ_j²  (exact identity; checked with enclosures)
    for k in [27usize, 81] {
        let op = build(k);
        let a0 = gram_first_row(&op)[0] as i128;
        let trace = k as i128 * a0;
        let (mut lo, mut hi) = (0i128, 0i128);
        for j in 1..=op.m() {
            let (l, h) = ordinate_enclosure(&op, j);
            let (al, ah) = if l > 0 { (l, h) } else { (-h, -l) };
            lo += 2 * (al / SCALE) * (al / SCALE);
            hi += 2 * (ah / SCALE + 1) * (ah / SCALE + 1);
        }
        assert!(lo <= trace && trace <= hi, "K={k}: {lo} <= {trace} <= {hi}");
    }
}

#[test]
fn inertia_matches_the_circulant_formula_at_k27() {
    let op = build(27);
    let ords = sorted_positive_ordinates(&op).unwrap();
    // below the first ordinate: only the zero mode is negative in A − τI
    assert_eq!(inertia_of_shifted_gram(&op, 1).unwrap(), 1);
    for r in 1..=op.m() {
        let a = ords[r - 1].hi / SCALE + 1;
        let b = if r < op.m() {
            ords[r].lo / SCALE
        } else {
            a + 2
        };
        assert!(a < b);
        let tau = (a * b) as i64;
        assert_eq!(
            inertia_of_shifted_gram(&op, tau).unwrap(),
            1 + 2 * r,
            "rank {r}"
        );
    }
}
