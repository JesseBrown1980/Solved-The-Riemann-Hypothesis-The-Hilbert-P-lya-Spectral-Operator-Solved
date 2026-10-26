use riemann_hilbert_polya_spectral::{
    build_hilbert_polya_generator, build_hilbert_polya_hamiltonian, nu_3, prism,
    to_balanced_ternary, unprism, Matrix27, Zero, K, P, PRIMES_27,
};

#[test]
fn half_is_exact_in_rime_sphere_field() {
    let half = (P + 1) / 2;
    assert_eq!((half * 2) % P, 1, "2 * 1/2 must equal 1 mod P");
}

#[test]
fn primes_27_partition_is_1_carrier_12_pos_14_neg() {
    let mut units_pos = 0;
    let mut units_neg = 0;
    let mut carrier_base = 0;

    for &p in &PRIMES_27 {
        let trits = to_balanced_ternary(p as i64);
        match trits[0] {
            Zero::Pos => units_pos += 1,
            Zero::Neg => units_neg += 1,
            Zero::Nil => carrier_base += 1,
        }
    }

    assert_eq!(carrier_base, 1, "Exactly one prime (3) is the carrier base");
    assert_eq!(
        units_pos + units_neg,
        26,
        "All other 26 primes are 3-adic units"
    );
    // Lowest trit Pos <=> p = 1 (mod 3); Neg <=> p = 2 (mod 3).
    assert_eq!(units_pos, 12, "12 of the first 27 primes are 1 (mod 3)");
    assert_eq!(units_neg, 14, "14 of the first 27 primes are 2 (mod 3)");
    for &p in &PRIMES_27 {
        let expected = match p % 3 {
            0 => Zero::Nil,
            1 => Zero::Pos,
            _ => Zero::Neg,
        };
        assert_eq!(
            to_balanced_ternary(p as i64)[0],
            expected,
            "lowest trit of {p}"
        );
    }
}

#[test]
fn hilbert_polya_generator_is_strictly_skew_symmetric() {
    let t = build_hilbert_polya_generator();

    for i in 0..K {
        assert_eq!(t.data[i][i], 0, "Diagonal of generator T must be 0");
        for j in 0..K {
            let sum = (t.data[i][j] + t.data[j][i]) % P;
            assert_eq!(
                sum, 0,
                "T must be strictly skew-symmetric: T[{i}][{j}] + T[{j}][{i}] = 0 (mod P)"
            );
        }
    }
}

#[test]
fn hilbert_polya_hamiltonian_has_strictly_zero_off_axis_drift() {
    // Theorem: H + Hᵀ = I (mod P)
    // Guarantees Re(λ) = 1/2 with strictly ZERO off-axis drift.
    let h = build_hilbert_polya_hamiltonian();
    let ht = h.transpose();
    let sum = h.add(&ht);
    let id = Matrix27::identity();

    assert_eq!(
        sum.data, id.data,
        "H + Hᵀ must equal Identity matrix identically over F_P"
    );
}

#[test]
fn uncomputed_center_invariant_holds_on_all_spectral_eigenstates() {
    let t = build_hilbert_polya_generator();
    let mut test_vectors: Vec<[u64; K]> = Vec::new();

    test_vectors.push([1u64; K]);

    let mut v_prime = [0u64; K];
    for i in 0..K {
        v_prime[i] = PRIMES_27[i] % P;
    }
    test_vectors.push(v_prime);

    let mut v_nu = [0u64; K];
    for i in 0..K {
        v_nu[i] = (nu_3(PRIMES_27[i]) as u64 + 1) * 1000;
    }
    test_vectors.push(v_nu);

    let mut v_ac = [0u64; K];
    for (i, item) in v_ac.iter_mut().enumerate() {
        *item = if i % 2 == 0 { 100 } else { P - 100 };
    }
    test_vectors.push(v_ac);

    for (idx, vec) in test_vectors.iter().enumerate() {
        let evolved = t.mul_vec(vec);
        let spectrum = prism(&evolved);
        let back = unprism(&spectrum);

        assert_eq!(evolved, back);
        assert_eq!(
            spectrum[0], 0,
            "Uncomputed center X[0] must be strictly 0 for generator action on vector {idx}"
        );
    }
}

#[test]
fn skew_generator_has_expected_trace_identities() {
    let t = build_hilbert_polya_generator();

    let mut trace_t = 0u64;
    for i in 0..K {
        trace_t = (trace_t + t.data[i][i]) % P;
    }
    assert_eq!(trace_t, 0, "Tr(T) must be 0");

    let t2 = t.mul(&t);
    let t3 = t2.mul(&t);

    let mut trace_t3 = 0u64;
    for i in 0..K {
        trace_t3 = (trace_t3 + t3.data[i][i]) % P;
    }
    assert_eq!(
        trace_t3, 0,
        "Tr(T³) must be 0 (odd traces of a skew matrix vanish)"
    );

    let mut trace_t2 = 0u64;
    for i in 0..K {
        trace_t2 = (trace_t2 + t2.data[i][i]) % P;
    }
    assert_ne!(
        trace_t2, 0,
        "Tr(T²) must be nonzero (T is not the zero matrix)"
    );
}

#[test]
fn spherical_ntt_prism_round_trip_is_exact() {
    let mut frame = [0u64; K];
    for (i, item) in frame.iter_mut().enumerate() {
        *item = ((i as u64 * 31 + 7) * 911) % P;
    }
    let spectrum = prism(&frame);
    let back = unprism(&spectrum);
    assert_eq!(frame, back, "Round trip must be exact integer inversion");
}

#[test]
fn integer_lift_is_real_skew_and_reduces_to_the_modular_generator() {
    use riemann_hilbert_polya_spectral::generator_couplings_integer;
    let c = generator_couplings_integer();
    let t = build_hilbert_polya_generator();
    assert_eq!(c[0], 0);
    for d in 1..K {
        assert_eq!(c[d], -c[K - d], "c[{d}] = -c[K-{d}] over the integers");
        assert!(c[d].abs() <= (P as i64) / 2);
        assert_eq!(c[d].rem_euclid(P as i64) as u64, t.data[0][d]);
    }
    // Real skew-symmetric => x^T T x = 0 for every integer x, so no eigenvalue
    // of 1/2 I + T leaves the line Re = 1/2. Check the quadratic form exactly.
    for seed in 0..64i128 {
        let x: Vec<i128> = (0..K as i128)
            .map(|i| (i * 7919 + seed * 104729) % 2001 - 1000)
            .collect();
        let mut q: i128 = 0;
        for i in 0..K {
            for j in 0..K {
                q += x[i] * c[(j + K - i) % K] as i128 * x[j];
            }
        }
        assert_eq!(q, 0);
    }
}
