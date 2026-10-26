# A finite Hilbert–Pólya operator on the 27-point tribit frame

[![Rust](https://img.shields.io/badge/rust-1.81.0-orange.svg)](rust-toolchain.toml)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)]()

Jesse Daniel Brown's finite approach to the Hilbert–Pólya programme: build an
exact, integer-only operator whose spectrum sits on the critical line by
construction, then test whether its ordinates follow the Riemann zeta zeros.

The repository does both halves.
- **The operator.** A skew-circulant generator T on the 27-point (3³) tribit
  frame, built from the order-27 root of unity in F_P and the first 27 primes.
  Its Hamiltonian H = ½I + T lies on the line Re = ½, exactly.
- **The test.** A pre-registered comparison of the operator's ordinates with the
  zeta zeros, at K = 27, 81, 243 and 729. Everything is frozen before any zero
  data is read.

## 1. The operator

Work in F_P with P = 1,000,081. W = 951,846 has multiplicative order 27, and
PRIMES₂₇ = 2, 3, 5, …, 103. For d = 1, …, 13 the couplings are

```text
c[d]      = (W^d − W^(−d)) · PRIMES₂₇[d]   (mod P)
c[27 − d] = −c[d]
T[i][j]   = c[(j − i) mod 27]
H         = ½I + T,   ½ = (P + 1)/2
```

**Theorem 1 (critical-line symmetry).** Tᵀ = −T and T has zero diagonal, so

```text
H + Hᵀ = (½I + T) + (½I − T) = I.
```

Over F_P this is an exact modular identity. Lift each c[d] to its balanced
integer representative in (−P/2, P/2]. The lift is still antisymmetric,
c[27 − d] = −c[d], so the integer T is a real skew-symmetric matrix. Its
spectrum over ℂ is purely imaginary,

```text
spec(T) = {0, ±iγ₁, …, ±iγ₁₃},
```

and every eigenvalue of H = ½I + T has real part exactly ½. The quadratic form
xᵀTx vanishes identically, which the tests check over the integers.

**Theorem 2 (centre invariant).** Each row of T sums to zero, so for every
state v the 27-point prism of Tv has DC component prism(Tv)[0] = 0.

**Theorem 3 (paired ordinates).** Odd powers of T are skew, so
Tr(T^(2k+1)) = 0. Every ordinate +γ comes with −γ, and the eigenvalues of H
pair as ½ ± iγ.

**The prime frame.** In balanced ternary, 3 is the single carrier prime. The
other 26 primes are 3-adic units: 12 with lowest trit Pos (p ≡ 1 mod 3) and 14
with lowest trit Neg (p ≡ 2 mod 3).

## 2. The experiment: do the ordinates follow the zeta zeros?

Theorem 1 puts the spectrum on Re = ½ for this construction. The step from
there to ζ(s) is a spectral correspondence: the ordinates γ would have to be
the zeta ordinates. That is the open research step, and the experiment tests it
directly.

The test is built to be decided by arithmetic alone.
- **Freeze first.** The coupling rule, the sizes K = 27, 81, 243, 729, the
  choice of P_K and W_K, the normalization and the match criterion are fixed
  in [`PREREGISTRATION.md`](PREREGISTRATION.md) before any zero table enters
  the repository.
- **Separate roles.** The builder computes the operator's ordinates and never
  sees zero data. A separate comparator reads the builder's receipts and the
  zero table.
- **Exact counts.** Ordinates are certified integer enclosures, from the
  circulant formula with interval trigonometry. They are cross-checked against
  exact integer inertia of A = −T² = TᵀT.
- **Fixed capacity.** A real skew circulant of odd size K has (K − 1)/2
  positive ordinates: 13, 40, 121 and 364 at K = 27, 81, 243 and 729.

Results are in [§4](#4-results).

## 3. Tests

```bash
cargo test
cargo clippy --all-targets -- -D warnings
```

| Test | What it checks |
| :--- | :--- |
| `half_is_exact_in_rime_sphere_field` | 2 · (P + 1)/2 = 1 in F_P |
| `primes_27_partition_is_1_carrier_12_pos_14_neg` | 3 is the carrier; 12 Pos and 14 Neg units, each prime's lowest trit |
| `hilbert_polya_generator_is_strictly_skew_symmetric` | Tᵀ = −T (mod P), zero diagonal |
| `hilbert_polya_hamiltonian_has_strictly_zero_off_axis_drift` | H + Hᵀ = I over F_P |
| `integer_lift_is_real_skew_and_reduces_to_the_modular_generator` | Integer lift antisymmetric, reduces to T mod P, xᵀTx = 0 |
| `uncomputed_center_invariant_holds_on_all_spectral_eigenstates` | prism(Tv)[0] = 0 |
| `skew_generator_has_expected_trace_identities` | Tr T = Tr T³ = 0, Tr T² ≠ 0 |
| `spherical_ntt_prism_round_trip_is_exact` | Exact prism inversion |
| `test_host8_300g_metadata_and_ground_alignment` | Host8 geometry; ground room is all Nil |

### Host8 integration lane

`test_host8_300g_readback_and_spectral_invariants` reads 71 fixed rooms from
the 300 GB Host8 store: sector ends, the ground room and its neighbours, the
prime ladder around the ground room, and strided checkpoints. Each stored stub
is inverted to its 27-trit room address and pushed through the operator, which
checks H + Hᵀ = I, the exact prism round trip and the centre invariant.

The lane needs the store, so it runs on request:

```bash
HOST8_DIR=/mnt/d/ASOLARIA-HOST8-ROOMS-300G \
  cargo test --test host8_riemann_experiment -- --ignored
```

It passes only when all 71 rooms are read and validated. A missing store or any
failed read fails the lane.

## 4. Results

The pre-registered comparison is recorded in [`PREREGISTRATION.md`](PREREGISTRATION.md)
and its receipts.

## Provenance

Theory: **Jesse Daniel Brown** (`JesseBrown1980`) and **Rayssa Chiqueto**.
Integer-only Rust 1.81.0, no external crates.
