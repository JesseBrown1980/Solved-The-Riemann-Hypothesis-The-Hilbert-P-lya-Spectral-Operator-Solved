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
| `k27_rule_reproduces_the_original_operator` | The frozen K-rule gives P, G, W and the couplings at K = 27 |
| `k81_rule_reproduces_the_clasy_x3_field` | The K-rule gives P = 1,000,999, W = 30,722 at K = 81 |
| `sin_enclosures_are_tight_and_consistent` | Interval sine widths, symmetry and known values |
| `trace_of_gram_matches_sum_of_squared_ordinates` | Tr(TᵀT) = 2 Σ γ_j² inside the enclosures |
| `inertia_matches_the_circulant_formula_at_k27` | Exact inertia counts equal 1 + 2r at every rank |

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

The experiment ran in the order fixed in [`PREREGISTRATION.md`](PREREGISTRATION.md).
Git history shows the freeze first, then the builder receipts, then the zero
table and comparator, then the verdicts.

**Builder.** Every ordinate enclosure is separated from its neighbours. Exact
integer inertia of A − τI agrees with the circulant formula at every checked
rank.

| K | P_K | ordinates | separated | inertia checks agreeing |
|---:|---:|---:|---:|---|
| 27 | 1,000,081 | 13 | 13 | 14 / 14 (all ranks) |
| 81 | 1,000,999 | 40 | 40 | 41 / 41 (all ranks) |
| 243 | 1,002,619 | 121 | 121 | 122 / 122 (all ranks) |
| 729 | 1,006,021 | 364 | 364 | 9 / 9 (ranks 0, 46, 91, 137, 182, 228, 273, 319, 364) |

**Zeros.** The table holds γ₁ … γ₁₀₀₀ at 50 digits, with enclosures of width
2·10⁻¹⁵. Each zero is isolated by a sign change of Hardy's Z function, and
`nzeros` confirms the count is complete at 1000. All 1000 values lie within
Odlyzko's published table (zeros6, stated accuracy 4·10⁻⁹).

**Comparison.** The smallest operator ordinate is anchored to γ₁ = 14.134725…,
and counts are compared between consecutive zeros.

| K | verdict | thresholds agreeing | max \|N_op − N_ζ\| | γ̂₂ vs γ₂ | γ̂_m vs γ_m |
|---:|---|---:|---:|---|---|
| 27 | MISMATCH | 1 / 13 | 11 | 31.04 vs 21.02 | 502.92 vs 59.35 |
| 81 | MISMATCH | 1 / 40 | 24 | 14.76 vs 21.02 | 531.68 vs 122.95 |
| 243 | MISMATCH | 1 / 121 | 83 | 14.80 vs 21.02 | 2457.93 vs 271.49 |
| 729 | MISMATCH | 1 / 364 | 354 | 29.39 vs 21.02 | 81620.23 vs 630.81 |

**Reading, as pre-registered.** The ordinates of this operator family, under
this anchoring, do not follow the zeta ordinates at any of the four sizes.
- After anchoring, the largest operator ordinate exceeds γ_m by a factor of
  8.5, 4.3, 9.1 and 129 at K = 27, 81, 243 and 729.
- The zeta ordinates grow like 2πn / log n.
- The count deviation grows with K.

This measures the coupling rule (W^d − W^(−d)) · prime[d]. Theorem 1, the
critical-line symmetry of the operator, is unaffected. Any next candidate rule
is evaluated the same way: frozen first, then compared on these four sizes.

Receipts: [`receipts/experiment/`](receipts/experiment/), with `SPECTRUM-K*.hbp`
from the builder and `COMPARE-K*.hbp` from the comparator. The zero table is
[`data/zeta/ZEROS-1000.hbp`](data/zeta/ZEROS-1000.hbp).

### The Relic-landing operator

The second experiment follows the lineage of
[THE RELIC REDISCOVERY](https://github.com/JesseBrown1980/THE-RELIC-REDISCOVERY).
At each landing, 27 → 243 → 6,561, the operator nests the previous one and adds
a new inner wave on fresh primes:

```text
T^(n+1) = T^(n) ⊗ I_q + I ⊗ S_q,   q = 3^(n+2)
```

Every landing keeps exactly one zero eigenvalue, the free centre. Every choice
was frozen in [`PREREGISTRATION-RELIC.md`](PREREGISTRATION-RELIC.md) before the
comparison.

**Builder.**
- **Factors.** Every factor passes exact inertia at every rank.
- **Landing 243.** Dense inertia of the full 243 × 243 matrix agrees at all 122
  ranks.
- **Landing 6,561.** One free centre plus 3,280 separated positive ordinates.

**Zeros.** γ₁ … γ₃₃₀₀, built the same way as the first table and complete at
3,300. All lie within Odlyzko's zeros6.

| Landing | verdict | thresholds agreeing | max \|N_op − N_ζ\| | γ̂₂ vs γ₂ | γ̂_m vs γ_m |
|---:|---|---:|---:|---|---|
| 27 | MISMATCH | 1 / 13 | 11 | 31.04 vs 21.02 | 502.92 vs 59.35 |
| 243 | MISMATCH | 1 / 121 | 100 | 33.48 vs 21.02 | 3118.36 vs 271.49 |
| 6,561 | MISMATCH | 1 / 3,280 | 2,930 | 27.70 vs 21.02 | 99889.38 vs 3809.87 |

**Reading, as pre-registered.** Under the same one-anchor normalization, the
nested operator's ordinates do not follow the zeta ordinates at any landing.
- At the top of the spectrum the nested operator overshoots by 11.5× at 243
  and 26× at 6,561. The flat family overshot by 129× at K = 729.
- The low-lying anchored ordinates sit closer to the zeros: γ̂₂ to γ̂₅ are 27.7,
  36.2, 39.9 and 51.4 at 6,561, against 21.0, 25.0, 30.4 and 32.9.
- The counts still separate from the first threshold on.

This result measures the nesting rule frozen in `PREREGISTRATION-RELIC.md`. The
operator's critical-line symmetry, Theorem 1, holds at every landing.

Receipts: [`receipts/relic/`](receipts/relic/), with `SPECTRUM-L*.hbp` and
`COMPARE-L*.hbp`. The zero table is
[`data/zeta/ZEROS-3300.hbp`](data/zeta/ZEROS-3300.hbp).

```bash
cargo run --release --bin build_spectrum -- 243 full > SPECTRUM-K243.hbp
cargo run --release --bin compare_zeta -- SPECTRUM-K243.hbp data/zeta/ZEROS-1000.hbp
```

## Provenance

Theory: **Jesse Daniel Brown** (`JesseBrown1980`) and **Rayssa Chiqueto**.
Integer-only Rust 1.81.0, no external crates.
