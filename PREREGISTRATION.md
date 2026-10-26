# Pre-registration: frozen-operator ordinates against the zeta zeros

This file fixes every choice in the experiment. It is committed together
with the builder, before any zeta-zero table or comparator code enters the
repository. Git history records that order.

## 1. The operator family (builder, `src/family.rs`)

For each K ∈ {27, 81, 243, 729}:

| Item | Rule |
|---|---|
| P_K | smallest prime P ≥ 1,000,000 with P ≡ 1 (mod K) |
| G_K | smallest primitive root mod P_K |
| W_K | G_K^((P_K − 1)/K), of multiplicative order K |
| primes | prime[0] = 2, prime[1] = 3, … (the first K primes) |
| couplings | c[d] = (W^d − W^(−d)) · prime[d] mod P, for d = 1..m, m = (K − 1)/2 |
| lift | c[d] → balanced representative in (−P/2, P/2]; c[K − d] = −c[d]; c[0] = 0 |
| generator | T[i][j] = c[(j − i) mod K], a real skew-symmetric integer circulant |

The rule reproduces the original constants at K = 27 (P = 1,000,081, G = 7,
W = 951,846) and the Clasy X3 field at K = 81 (P = 1,000,999, W = 30,722).
The tests check both.

| K | P_K | G_K | W_K | positive ordinates m |
|---:|---:|---:|---:|---:|
| 27 | 1,000,081 | 7 | 951,846 | 13 |
| 81 | 1,000,999 | 6 | 30,722 | 40 |
| 243 | 1,002,619 | 2 | 990,283 | 121 |
| 729 | 1,006,021 | 2 | 635,719 | 364 |

The construction has no free parameter, and the builder reads no zeta data.

## 2. The operator's ordinates

spec(T) = {0, ±iγ_1, …, ±iγ_m} with

```text
γ_j = 2 Σ_{d=1..m} c[d] · sin(2π j d / K),   j = 1..m.
```

- **Enclosures.** Each γ_j is enclosed at scale 10^18 by integer interval
  trigonometry (`src/interval.rs`). π is enclosed as
  [3.141592653589793238, 3.141592653589793239], Taylor series are summed with
  outward rounding, and the tail is added as a bound.
- **Ordering.** The positive ordinates are the |γ_j|, sorted. Adjacent
  enclosures must be disjoint, and the receipt records `separated` for each
  rank.
- **Inertia cross-check.** A = −T² = TᵀT is an integer positive-semidefinite
  circulant. For each checked rank r, an integer τ is chosen strictly between
  |γ|²₍r₎ and |γ|²₍r+1₎. The exact count of negative eigenvalues of A − τI must
  equal 1 + 2r: one zero mode, plus each ordinate counted twice.
  - The count is the number of sign changes in the leading principal minors.
  - Each minor is computed modulo 61-bit primes, enough that their product
    exceeds twice the Hadamard bound.
  - The sign of each minor is read from balanced mixed-radix digits.
  - This is an independent diagonalization-free route and uses no
    trigonometry.
- **Inertia coverage.**
  - K = 27, 81, 243: every rank r = 0..m.
  - K = 729: rank 0 and ranks ⌈364·i/8⌉ for i = 1..8, which are 46, 91, 137,
    182, 228, 273, 319 and 364.
- **Builder receipt.** `receipts/experiment/SPECTRUM-K{K}.hbp`, written by
  `build_spectrum K full` (or `spot8` for K = 729).

## 3. Zeta-zero data (comparator input)

- **Source.** The ordinates γ^ζ_n, n = 1..1000, are computed with mpmath
  `zetazero` at 50 significant digits. The values for n = 1..10 are checked
  against the published values in Odlyzko's tables.
- **Enclosure.** For each n the stored rational enclosure is
  [ℓ_n, u_n] = [⌊γ̃·10^18⌋ − 1000, ⌈γ̃·10^18⌉ + 1000] / 10^18.
- **Isolation evidence.** Hardy's Z function changes sign between ℓ_n and u_n,
  evaluated at 50 digits.
- **Completeness evidence.** mpmath `nzeros` at the midpoint between zeros 1000
  and 1001 returns 1000.
- **Data file.** `data/zeta/ZEROS-1000.hbp`, which also records its provenance.

## 4. Normalization (one anchor)

The operator's ordinates are in arbitrary units. A single scale is fixed by
anchoring the smallest positive operator ordinate to the first zeta ordinate:

```text
σ_K = |γ|₍1₎ / γ^ζ_1,     γ̂_n = |γ|₍n₎ / σ_K.
```

No other parameter is fitted, and the anchor is the same rule at every K.

## 5. The comparison

For k = 1..m, the threshold is t_k = ⌊(u_k + ℓ_{k+1}) / 2⌋ / 10^18. It lies
strictly between zeta zeros k and k + 1, so N_ζ(t_k) = k exactly. The
comparator checks u_k < t_k < ℓ_{k+1}.

N_op(t_k) is the number of ranks n with γ̂_n ≤ t_k. With interval enclosures:

| Status of rank n | Condition |
|---|---|
| in | hi(|γ|₍n₎) · u_1 ≤ t_k · lo(|γ|₍1₎) |
| out | lo(|γ|₍n₎) · ℓ_1 > t_k · hi(|γ|₍1₎) |
| ambiguous | otherwise |

Rank 1 is in for every k, by the anchor. All arithmetic is exact integer
arithmetic, with 256-bit products.

Threshold k is called:

| Result | Condition |
|---|---|
| AGREE | N_op(t_k) = k with no ambiguous rank |
| DISAGREE | k lies outside [#in, #in + #ambiguous] |
| AMBIGUOUS | otherwise |

**Verdict per K, fixed now:**

| Verdict | Condition |
|---|---|
| MATCH | all m thresholds AGREE |
| MISMATCH | at least one threshold DISAGREES |
| INCONCLUSIVE | otherwise |

The verdict line also reports these descriptive statistics, which are not
criteria:
- the number of agreeing thresholds;
- the first disagreeing k;
- max |N_op(t_k) − k|;
- the anchored γ̂_n beside γ^ζ_n for n ≤ 10, and for n = m.

## 6. Interpretation, fixed now

- **MATCH.** Finite evidence for this particular construction at that K.
- **MISMATCH.** Evidence against this construction's ordinates being the
  zeta ordinates under this normalization. It bears on this operator family
  only.
- **Rule changes.** Any change to the coupling rule made after the comparison
  is labelled exploratory. It is then evaluated on sizes or zeros it has not
  been compared against.

## 7. Roles

- **Builder:** `src/family.rs`, `src/interval.rs` and
  `src/bin/build_spectrum.rs`. It reads no files.
- **Comparator:** `src/bin/compare_zeta.rs`, added after the builder receipts.
  It reads only `receipts/experiment/SPECTRUM-K*.hbp` and
  `data/zeta/ZEROS-1000.hbp`, and does not link the builder modules.
