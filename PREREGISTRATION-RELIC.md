# Pre-registration: the Relic-landing operator against the zeta zeros

This file fixes every choice of the second experiment before any comparison of
the nested operator with zeta data. It is committed together with the frozen
builder. Git history records that order.

The first experiment ([`PREREGISTRATION.md`](PREREGISTRATION.md)) tested the
flat circulant family at K = 27, 81, 243 and 729. Its verdicts stay attached to
that flat operator. This experiment tests a different operator: the nested
construction at the Relic landings.

## 1. Landings

The landing sizes follow THE RELIC REDISCOVERY: A_0 = 27 and
A_(n+1) = A_n · 3^(n+2).

| Level n | Size A_n | Positive ordinates m = (A_n − 1)/2 |
|---:|---:|---:|
| 0 | 27 | 13 |
| 1 | 243 | 121 |
| 2 | 6,561 | 3,280 |

## 2. The operator (builder: `src/nested.rs`, `src/bin/build_relic.rs`)

```text
T^(0)   = the frozen K = 27 circulant of PREREGISTRATION.md
T^(n+1) = T^(n) ⊗ I_q + I_(A_n) ⊗ S_q,      q = 3^(n+2)
index   (x, y) ↦ x·q + y,   x the outer index, y the inner index
```

**Inner waves.** S_q is the family circulant of size q. It uses the same P_q,
G_q and W_q rule as `PREREGISTRATION.md` §1, but its couplings take fresh
primes:

```text
c[d] = (W^d − W^(−d)) · prime[offset + d]  (mod P_q),   d = 1..(q − 1)/2,
```

lifted to the balanced representative in (−P/2, P/2], with c[q − d] = −c[d].
The offset counts every coupling used by the earlier factors, so no prime is
used twice.

| Factor | Size | P | G | W | Prime offset | Primes used |
|---|---:|---:|---:|---:|---:|---|
| outer T^(0) | 27 | 1,000,081 | 7 | 951,846 | 0 | prime[1..13] = 3 … 43 |
| inner S_9 | 9 | 1,000,081 | 7 | 53,001 | 13 | prime[14..17] = 47 … 59 |
| inner S_27 | 27 | 1,000,081 | 7 | 951,846 | 17 | prime[18..30] = 67 … 127 |

**Why fresh primes.** Before this freeze, reusing build(27) as the inner wave at
level 2 made S_27 a copy of T^(0). The copy cancelled exactly and left 27 zero
modes. Fresh primes keep the inner network distinct.

## 3. Ordinates and the free centre

Every factor is real skew-symmetric, so T^(n) is too. Its spectrum is the
Kronecker sum:

```text
spec T^(n+1) = { i(a + b) : ia ∈ spec T^(n), ib ∈ spec S_q }.
```

- **Enclosures.** Each eigenvalue's enclosure is the sum of its factors'
  enclosures, which come from the interval trigonometry at scale 10^18.
- **Free centre.** The zero mode 0 + 0 + … is the single free centre. Every
  other enclosure must exclude 0, and the receipt must report exactly one zero
  mode and (A_n − 1)/2 positive ordinates.
- **Ordering.** The positive ordinates are sorted, and adjacent enclosures must
  be disjoint.

## 4. Exact cross-checks

| Check | Coverage |
|---|---|
| Factor inertia: neg(SᵀS − τI) = 1 + 2r | every rank of each factor (27, 9, 27) |
| Landing 1 inertia: exact inertia of the dense 243 × 243 matrix TᵀT − τI | every rank 0..121 |
| Landing 2 | the Kronecker-sum theorem applied to the inertia-validated factors |

The landing 1 check is computed directly from the dense matrix and does not use
the Kronecker formula. The multi-modular minors are the same as in
`PREREGISTRATION.md` §2.

**Receipts.** `receipts/relic/SPECTRUM-L{0,1,2}.hbp`, written by
`build_relic n full` for n = 0, 1 and `build_relic 2 none`.

## 5. Zero data

`data/zeta/ZEROS-3300.hbp` is built the same way as `ZEROS-1000.hbp`
(`data/zeta/make_zeros.py`, run with N = 3300):

| Item | Rule |
|---|---|
| Source | mpmath `zetazero` at 50 digits |
| Enclosure | [⌊γ·10^18⌋ − 1000, ⌈γ·10^18⌉ + 1000] |
| Isolation | Hardy's Z changes sign across each enclosure |
| Completeness | `nzeros` at the midpoint between zeros 3,300 and 3,301 returns 3,300 |
| Cross-check | all 3,300 values against Odlyzko's zeros6 table, stated accuracy 4·10⁻⁹ |

Landing 2 needs zeros 1 to 3,281.

## 6. Comparison and verdicts

The comparator is `src/bin/compare_zeta.rs`, unchanged, with the rules of
`PREREGISTRATION.md` §4–§5:
- **Anchor.** The smallest positive ordinate is anchored to γ₁. No other
  parameter is fitted.
- **Thresholds.** t_k sits between zeta zeros k and k + 1, for k = 1..m.
- **Counting.** Each rank is in, out or ambiguous, using exact 256-bit
  integer products.

| Status | Meaning |
|---|---|
| AGREE | N_op(t_k) = k with no ambiguous rank |
| DISAGREE | k lies outside [#in, #in + #ambiguous] |
| AMBIGUOUS | otherwise; reported, never rounded away |

| Verdict per landing | Condition |
|---|---|
| MATCH | all m thresholds AGREE |
| MISMATCH | at least one DISAGREE |
| INCONCLUSIVE | no DISAGREE and at least one AMBIGUOUS |

**Receipts.** `receipts/relic/COMPARE-L{0,1,2}.hbp`.

Landing 0 is the original K = 27 operator, so its comparison repeats the first
experiment's K = 27 result. The new information is at landings 1 and 2.

## 7. Interpretation, fixed now

- **MATCH.** Finite evidence for this nested construction at that landing.
- **MISMATCH.** Evidence against this nested construction's ordinates being
  the zeta ordinates under this anchoring. It says nothing about the flat
  family, or about other nestings.
- **Rule changes.** Any change to the nesting rule, the factor sizes, the prime
  assignment or the anchor after the comparison is labelled exploratory. It is
  then evaluated on landings or zeros it has not been compared against.
