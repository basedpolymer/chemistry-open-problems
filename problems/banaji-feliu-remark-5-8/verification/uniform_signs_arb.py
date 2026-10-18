#!/usr/bin/env python3
"""Interval-arithmetic (python-flint / Arb) check of the sign pattern used in the manuscript's Proposition, for R = 64
and n = 1..300 (the Proposition itself is proved by hand for all n; this only guards against slips in the algebra):

  f(u) = log u + 2 sum_{m=1}^n (-1)^m log(R^m + u),   c = M_1 - L/2,   M_1 = L (1 + 2 sum_m (-1)^m m),  L = log R,

  claim: (-1)^(i+1) (f(R^i) - c) >= log 2 - 4/(R-1) > 0 for 1 <= i <= n.
Also checks the closed form e^c = 2^(-6n-3) (n odd), 2^(6n+3) (n even) used in direct_check.py.
Runtime: ~10 s.  Prints PASS/FAIL, exit code 1 on failure.
"""
import sys
from flint import arb, ctx

ctx.prec = 200
R = 64
L = arb(R).log()
ok = True
worst = None
for n in range(1, 301):
    M1 = L * (1 + 2 * sum((-1) ** m * m for m in range(1, n + 1)))
    c = M1 - L / 2
    k = 6 * n + 3
    expected = -k if n % 2 else k                      # log2(e^c)
    ok &= abs(c - expected * arb(2).log()) < arb(10) ** -50
    for i in range(1, n + 1):
        u = arb(R) ** i
        f = u.log() + 2 * sum((-1) ** m * (arb(R) ** m + u).log() for m in range(1, n + 1))
        s = (f - c) * (-1) ** (i + 1)
        margin = arb(2).log() - arb(4) / (R - 1)
        if not (s > margin):                           # rigorous comparison of balls
            ok = False
            print("FAIL at", n, i, s)
        lo = s.lower()
        worst = lo if worst is None or lo < worst else worst
print(("PASS" if ok else "FAIL") + f" R = 64, n = 1..300: (-1)^(i+1)(f(R^i) - c) > log 2 - 4/63 for all i, "
      f"and e^c = 2^(+-(6n+3)); smallest certified value {worst.mid().str(8)} (log 2 - 4/63 = {(arb(2).log() - arb(4) / 63).mid().str(8)})")
sys.exit(0 if ok else 1)
