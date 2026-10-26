"""Zeta-zero table for the pre-registered comparison (comparator input).

gamma_n for n = 1..1000 via mpmath.zetazero at 50 digits. Enclosure at scale 1e18:
[floor(g*1e18) - 1000, ceil(g*1e18) + 1000]. Isolation: Hardy Z changes sign
across each enclosure. Completeness: mpmath.nzeros at the midpoint between zeros
1000 and 1001 equals 1000.
"""
import sys, mpmath
mpmath.mp.dps = 50
S = 10**18
N = 1000
out = []
gs = []
for n in range(1, N + 2):
    g = mpmath.zetazero(n).imag
    gs.append(g)
lines = ["HBP ZETA-ZEROS v1",
         f"SOURCE mpmath-{mpmath.__version__} zetazero dps=50",
         f"ENCLOSURE scale={S} lo=floor(g*scale)-1000 hi=ceil(g*scale)+1000",
         "ISOLATION hardy_Z_sign_change_across_[lo,hi] dps=50"]
for n in range(1, N + 1):
    g = gs[n - 1]
    lo = int(mpmath.floor(g * S)) - 1000
    hi = int(mpmath.ceil(g * S)) + 1000
    zl = mpmath.siegelz(mpmath.mpf(lo) / S)
    zh = mpmath.siegelz(mpmath.mpf(hi) / S)
    assert zl * zh < 0, (n, zl, zh)
    if n > 1:
        assert lo > prev_hi
    prev_hi = hi
    lines.append(f"Z n={n} lo={lo} hi={hi}")
mid = (gs[N - 1] + gs[N]) / 2
cnt = mpmath.nzeros(mid)
assert cnt == N, cnt
lines.append(f"COMPLETENESS nzeros(t)={cnt} t_between_zero_{N}_and_{N+1}")
sys.stdout.write("\n".join(lines) + "\n")
