use riemann_hilbert_polya_spectral::family::build;
use riemann_hilbert_polya_spectral::interval::SCALE;
use riemann_hilbert_polya_spectral::nested::{
    build_with_prime_offset, dense_generator, factors, inertia_dense_gram, landing, landing_size,
    positive_ordinates,
};

#[test]
fn landing_sizes_follow_the_relic_closed_form() {
    // A_n = 3^(3 + n(n+3)/2)
    for n in 0..3usize {
        assert_eq!(landing_size(n), 3usize.pow(3 + (n * (n + 3) / 2) as u32));
    }
    assert_eq!(
        [landing_size(0), landing_size(1), landing_size(2)],
        [27, 243, 6561]
    );
}

#[test]
fn each_wave_takes_fresh_primes() {
    assert_eq!(factors(2), vec![(27, 0), (9, 13), (27, 17)]);
    assert_eq!(build_with_prime_offset(27, 0).c, build(27).c);
    assert_ne!(build_with_prime_offset(27, 17).c, build(27).c);
}

#[test]
fn every_landing_has_one_free_centre_and_separated_pairs() {
    for level in 0..3 {
        let l = landing(level);
        let pos = positive_ordinates(&l).unwrap();
        assert_eq!(2 * pos.len() + 1, l.size);
        for w in pos.windows(2) {
            assert!(w[0].1 < w[1].0, "separated");
        }
    }
}

#[test]
fn dense_243_generator_inertia_matches_the_kronecker_spectrum_at_spot_ranks() {
    let (n, t) = dense_generator(1);
    assert_eq!(n, 243);
    let pos = positive_ordinates(&landing(1)).unwrap();
    for r in [0usize, 1, 60, 121] {
        let tau = if r == 0 {
            1
        } else {
            let a = pos[r - 1].1 / SCALE + 1;
            let b = if r < pos.len() {
                pos[r].0 / SCALE
            } else {
                a + 2
            };
            (a * b) as i64
        };
        assert_eq!(
            inertia_dense_gram(n, &t, tau).unwrap(),
            1 + 2 * r,
            "rank {r}"
        );
    }
}
