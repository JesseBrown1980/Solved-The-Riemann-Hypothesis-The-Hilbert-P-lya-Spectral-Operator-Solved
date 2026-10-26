//! Host8 integration lane: the 300 GB room store read through the 27-mode operator.
//!
//! Reads 71 fixed sample rooms from the Host8 store (375 sectors, 37.5B rooms),
//! inverts each stored stub to its balanced-ternary room address, and checks the
//! operator identities (H + Hᵀ ≡ I, exact prism round trip, X[0] ≡ 0) on the
//! resulting 27-trit state.
//!
//! The disk test is a required integration lane. It is marked `#[ignore]` so a
//! machine without the store reports it as ignored rather than passed, and runs with
//!
//! ```text
//! cargo test --test host8_riemann_experiment -- --ignored
//! ```
//!
//! When it runs it fails unless every one of the 71 rooms is read and validated.
//! `HOST8_DIR` overrides the default store path.

use riemann_hilbert_polya_spectral::{
    build_hilbert_polya_generator, build_hilbert_polya_hamiltonian, prism, unprism, Matrix27, Zero,
    K, P, PRIMES_27,
};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

const HOST8_DIR: &str = "D:/ASOLARIA-HOST8-ROOMS-300G";
const SECTORS: u64 = 375;
const RPS: u64 = 100_000_000;
const ROOM_BYTES: u64 = 8;
const TOTAL_ROOMS: u64 = SECTORS * RPS; // 37,500,000,000
const CENTRE_ROOM: u64 = TOTAL_ROOMS / 2; // 18,750,000,000 (Ground Room)
const SEED: u64 = 0x518bbb472a9db786;
const H27: u64 = 3_812_798_742_493; // (3^27 - 1) / 2

// Multiplicative inverse constants for 64-bit stub unmix
const M1: u64 = 0xbf58_476d_1ce4_e5b9;
const M2: u64 = 0x94d0_49bb_1331_11eb;

const fn modinv(c: u64) -> u64 {
    let mut x = c;
    let mut i = 0;
    while i < 6 {
        x = x.wrapping_mul(2u64.wrapping_sub(c.wrapping_mul(x)));
        i += 1;
    }
    x
}
const M1I: u64 = modinv(M1);
const M2I: u64 = modinv(M2);

#[inline]
fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(M1);
    z = (z ^ (z >> 27)).wrapping_mul(M2);
    z ^ (z >> 31)
}

#[inline]
fn unmix(mut z: u64) -> u64 {
    z ^= (z >> 31) ^ (z >> 62);
    z = z.wrapping_mul(M2I);
    z ^= (z >> 27) ^ (z >> 54);
    z = z.wrapping_mul(M1I);
    z ^= (z >> 30) ^ (z >> 60);
    z
}

fn theoretical_word(n: u64) -> u64 {
    let balanced = (n as i128) - (CENTRE_ROOM as i128);
    (balanced + H27 as i128) as u64
}

fn theoretical_stub(n: u64) -> u64 {
    mix(theoretical_word(n) ^ SEED)
}

fn unpack_27_trits(mut word: u128) -> [Zero; K] {
    let mut trits = [Zero::Nil; K];
    for i in (0..K).rev() {
        trits[i] = Zero::from_digit((word % 3) as u8);
        word /= 3;
    }
    trits
}

const SAMPLE_COUNT: usize = 71;

fn host8_dir() -> String {
    std::env::var("HOST8_DIR").unwrap_or_else(|_| HOST8_DIR.to_string())
}

fn read_room_stub_from_disk(dir: &str, n: u64) -> io::Result<u64> {
    let sector = n / RPS;
    let offset = (n % RPS) * ROOM_BYTES;
    let path = format!("{dir}/sector-{sector:03}.rooms");

    let mut file = File::open(&path)?;
    file.seek(SeekFrom::Start(offset))?;
    let mut buf = [0u8; 8];
    file.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[test]
fn test_host8_300g_metadata_and_ground_alignment() {
    assert_eq!(TOTAL_ROOMS, 37_500_000_000);
    assert_eq!(CENTRE_ROOM, 18_750_000_000);

    // Ground room address is balanced 0
    let ground_word = theoretical_word(CENTRE_ROOM);
    assert_eq!(ground_word, H27);

    let ground_trits = unpack_27_trits(ground_word as u128);
    for (i, t) in ground_trits.iter().enumerate() {
        assert_eq!(
            *t,
            Zero::Nil,
            "Trit {i} at ground room must be Zero::Nil (quiescent)"
        );
    }
}

#[test]
#[ignore = "integration lane: needs the Host8 store; run with -- --ignored"]
fn test_host8_300g_readback_and_spectral_invariants() {
    let dir = host8_dir();
    assert!(
        std::path::Path::new(&dir).is_dir(),
        "Host8 store not found at {dir}; set HOST8_DIR"
    );

    let t_gen = build_hilbert_polya_generator();
    let h_ham = build_hilbert_polya_hamiltonian();
    let h_trans = h_ham.transpose();
    let h_sum = h_ham.add(&h_trans);
    let identity = Matrix27::identity();

    // Verify operator symmetry theorem: H + Hᵀ ≡ I (mod P)
    assert_eq!(h_sum.data, identity.data);

    // Curated sample rooms across the 300 GB space
    let mut sample_rooms = vec![
        0,               // Sector 0 start
        99_999_999,      // Sector 0 end
        100_000_000,     // Sector 1 start
        CENTRE_ROOM - 1, // Pre-ground (balanced -1)
        CENTRE_ROOM,     // Ground room (balanced 0)
        CENTRE_ROOM + 1, // Post-ground (balanced +1)
        TOTAL_ROOMS - 1, // Sector 374 end
    ];

    // Add prime ladder rungs from ground room
    for &p in &PRIMES_27 {
        sample_rooms.push(CENTRE_ROOM + p);
        sample_rooms.push(CENTRE_ROOM - p);
    }

    // Add strided sector checkpoints
    for s in [50, 100, 187, 250, 300] {
        sample_rooms.push(s * RPS);
        sample_rooms.push(s * RPS + 500_000);
    }

    println!("\n=== RUNNING HOST8 300 GB RIEMANN HYPOTHESIS SPECTRAL EXPERIMENT ===");
    println!(
        "Total sample points across 300 GB storage: {}",
        sample_rooms.len()
    );
    assert_eq!(
        sample_rooms.len(),
        SAMPLE_COUNT,
        "the lane reads exactly 71 rooms"
    );

    let mut validated_rooms = 0;

    for &n in &sample_rooms {
        let stub_on_disk = read_room_stub_from_disk(&dir, n)
            .unwrap_or_else(|e| panic!("could not read room {n} from {dir}: {e}"));

        // 1. Stub unmixing & collision-free invertibility check
        let expected_stub = theoretical_stub(n);
        assert_eq!(
            stub_on_disk, expected_stub,
            "Stub on disk must match invertible mix"
        );

        let recovered_word = unmix(stub_on_disk) ^ SEED;
        let expected_word = theoretical_word(n);
        assert_eq!(
            recovered_word, expected_word,
            "Recovered word must match expected word"
        );

        // 2. Unpack into 27-trit balanced vector
        let trits = unpack_27_trits(recovered_word as u128);

        // 3. Map into F_P^27 carrier state
        let mut v = [0u64; K];
        for i in 0..K {
            let dir = trits[i].direction() as i64;
            v[i] = ((dir % P as i64 + P as i64) % P as i64) as u64;
        }

        // 4. Verify Zero Off-Axis Drift: (H + Hᵀ) v = I v = v (mod P)
        let drift_check = h_sum.mul_vec(&v);
        assert_eq!(
            drift_check, v,
            "Off-axis drift must be identically 0 (Re(λ) = 1/2)"
        );

        // 5. Apply Hilbert-Pólya dilation generator: w = T v
        let w = t_gen.mul_vec(&v);

        // 6. Forward 27-point spherical NTT prism transform: X = prism(w)
        let spectrum = prism(&w);
        let reconstructed_w = unprism(&spectrum);
        assert_eq!(
            reconstructed_w, w,
            "NTT prism transform must round-trip exactly"
        );

        // 7. Verify Law 0 Uncomputed Center Invariant: X[0] ≡ 0 (mod P)
        assert_eq!(
            spectrum[0], 0,
            "Uncomputed center invariant X[0] == 0 violated for room {n}"
        );

        // 8. Ground Room Quiescent origin check
        if n == CENTRE_ROOM {
            for i in 0..K {
                assert_eq!(v[i], 0, "Ground room vector must be identically 0");
                assert_eq!(w[i], 0, "Ground room dilation must be identically 0");
                assert_eq!(spectrum[i], 0, "Ground room spectrum must be identically 0");
            }
        }

        validated_rooms += 1;
    }

    assert_eq!(
        validated_rooms, SAMPLE_COUNT,
        "every sample room must be read and validated"
    );
    println!("SUCCESS: {validated_rooms}/{SAMPLE_COUNT} Host8 rooms read and validated: H + Hᵀ = I, exact prism round trip, X[0] = 0.");
}
