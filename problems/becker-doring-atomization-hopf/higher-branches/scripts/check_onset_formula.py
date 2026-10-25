#!/usr/bin/env python3
"""check_onset_formula.py -- exact identities behind the onset chart (result.md, Sections 2-3).

Environment: Python 3.13, sympy, python-flint 0.9.  Runtime: about 1 min.
Claims (PASS/FAIL):
  O1  (sympy; symbols s, phi, M, P = phi^M, Q = s^M, A = 1/s) PV's F (71)-(75), with (A phi)^-M = Q/P, satisfies
      F/A^3 = f1 (H1 - s^(M-1) phi^-M H2) - (G1 - phi^M G2) f2 with the s-polynomials of onset_sweep.py.
  O2  (sympy) lambda/A = (1 - s/phi)(phi - 1) for lambda = (A - 1/phi)(phi - 1).
  O3  at s = 0 the s-form equals -phi^3 G(phi) (kinf.py's G), i.e. the limit of Lemma 2.4.
  O4  (Arb) onset_sweep.system agrees with result 8's bdlib.hopf_system at integer M (M = 50, 142, 300; random
      phi, K): F/A^3 (onset chart) = K^3/A^3 * K^-3 F (bdlib), lambda/A = (K/M) Lambda / A; the balls overlap.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from flint import arb, acb, ctx

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

s, phi, M, P, Q = sp.symbols("s phi M P Q")
A = 1 / s; K = A - 1
ApM = Q / P                                   # (A phi)^-M = s^M phi^-M
f1 = (phi - 1) + P * A
f2 = (1 - A * phi) + ApM * A ** 2 * phi
G1 = phi * (1 - A * phi) * (K * phi + (K * M + 1) * (phi - 1))
G2 = (1 - A * phi) * (A * M * (phi - 1) + K * phi + A * phi * (phi - 1))
H1 = (phi - 1) * (K + (K * M + 1) * (1 - A * phi))
H2 = (phi - 1) * ((A + M * A ** 2 * phi) * (1 - A * phi) + K * A * phi)
F = f1 * (H1 - ApM * H2) - (G1 - P * G2) * f2
SM1 = Q / s                                   # s^(M-1)
g1 = P + s * (phi - 1)
g2 = s - phi + SM1 * phi / P                  # s^(M-1) phi^(1-M)
a = M * (1 - s) + s
GG1 = phi * (s - phi) * ((1 - s) * phi + a * (phi - 1))
GG2 = (s - phi) * (M * (phi - 1) + (1 - s) * phi + phi * (phi - 1))
HH1 = (phi - 1) * (s * (1 - s) + a * (s - phi))
HH2 = (phi - 1) * ((s + M * phi) * (s - phi) + s * (1 - s) * phi)
Fs = g1 * (HH1 - SM1 / P * HH2) - (GG1 - P * GG2) * g2
report("O1", sp.simplify(sp.expand(F * s ** 3 - Fs)) == 0, "F/A^3 = s-form (symbolic in s, phi, M, phi^M, s^M)")
lam = (A - 1 / phi) * (phi - 1)
report("O2", sp.simplify(lam / A - (1 - s / phi) * (phi - 1)) == 0, "lambda/A = (1 - s/phi)(phi - 1)")
Gs = -P / phi ** 2 * (phi ** 3 + M * phi ** 2 - 2 * M * phi + M) + (M + 1) * phi - M
F0 = Fs.subs(Q, 0).subs(s, 0)
report("O3", sp.simplify(F0 + phi ** 3 * Gs) == 0, "s = 0 (and s^M = 0): s-form = -phi^3 G")

import onset_sweep as O
import bdlib as bd
ctx.prec = 192
random.seed(1)
ok4 = True
for Mi in (50, 142, 300):
    for _ in range(5):
        Kv = arb(random.uniform(0.5, 50)); ph = acb(1 + random.uniform(-0.05, 0.05), random.uniform(0.02, 0.3))
        sv = 1 / (1 + Kv)
        Fo, Lo = O.system(sv, arb(Mi), ph)
        kap = Kv * arb(Mi).sqrt(); z = (ph - 1) * Mi
        Qe, Lam = bd.hopf_system(kap, z, Mi)
        Av = 1 + Kv
        ok4 &= (Fo.v - Qe.v * Kv ** 3 / Av ** 3).contains(0)
        ok4 &= (Lo.v - Lam.v * Kv / Mi / Av).contains(0)
report("O4", ok4, "onset_sweep.system = bdlib.hopf_system (rescaled) at 15 random points, M = 50, 142, 300")
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
