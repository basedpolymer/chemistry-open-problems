#!/usr/bin/env python3
"""
table1.py -- rigorous enclosures of the two coefficients of Table 1 of Aminian-Camassa-McLaughlin
(Stud. Appl. Math. 141 (2018), p. 407; arXiv:1708.05097 v3 Table 1), for the regular n-gon with
circumradius 1 (the paper's normalization):
   LT coefficient  kappa_n = 3<u g2>/(2<u g1>)^(3/2) = 6<u g1^2>/(2<u g1>)^(3/2)   (scales like length)
   geometric skewness S^G_n = <u^3>/<u^2>^(3/2)                                   (scale invariant)
n = 3 is exact (triangle_exact.py); n >= 4 from lt_certify.certify / sg_certify.certify_sg with larger K.
Reports AGREES/DIFFERS: whether each printed value of the paper lies in our enclosure widened by the paper's
rounding (half a unit of the last printed digit). A DIFFERS is a finding about the paper's numerics, not a failure.
Runtime: several minutes (one core).
"""
import sys, os, math
from flint import arb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lt_certify import certify
from sg_certify import certify_sg

# (n, printed LT coeff, printed SG, K for LT, K for SG); arXiv v3 values (published values are rounded versions)
ROWS = [(4, '0.139', '0.0813', 24, 20), (5, '0.194', '0.042', 24, 20), (6, '0.215', '0.024', 20, 20),
        (7, '0.225', '0.0149', 16, 16), (8, '0.231', '0.0097', 16, 16), (20, '0.243', '3.97e-4', 10, 10),
        (100, '0.245', '8.6e-7', 4, 10)]
ok_all = True


def consistent(enc, printed):
    s = printed
    if 'e' in s:
        mant, ex = s.split('e')
        dec = len(mant.split('.')[1]) if '.' in mant else 0
        half = 0.5 * 10 ** (int(ex) - dec)
    else:
        dec = len(s.split('.')[1])
        half = 0.5 * 10 ** (-dec)
    v = float(s)
    lo, hi = float(enc.lower()), float(enc.upper())
    return lo <= v + half and hi >= v - half


print("n=3 (exact): kappa(R=1) = -789 sqrt(770)/618800 = -0.03538 (paper -0.035), S^G = 4 sqrt(21)/99 = 0.18515 (paper 0.185)")
for n, lt, sg, Kl, Ks in ROWS:
    r = certify(n, K=Kl, verbose=False)
    kap = arb(r['kappa_circumradius1'][0]).union(arb(r['kappa_circumradius1'][1]))
    c1 = consistent(kap, lt)
    line = (f"n={n:3d}: kappa(R=1) in [{float(kap.lower()):.5f}, {float(kap.upper()):.5f}] (paper {lt}) "
            f"{'AGREES' if c1 else 'DIFFERS'}")
    ok_all &= c1
    if Ks is not None:
        rs = certify_sg(n, K=Ks, verbose=False)
        sgv = arb(rs['SG'][0]).union(arb(rs['SG'][1]))
        c2 = consistent(sgv, sg)
        ok_all &= c2
        line += (f";  S^G in [{float(sgv.lower()):.5g}, {float(sgv.upper()):.5g}] (paper {sg}) "
                 f"{'AGREES' if c2 else 'DIFFERS'}")
    print(line, flush=True)
print("All printed values agree with the enclosures up to rounding." if ok_all else
      "Some printed values of the paper lie outside the rigorous enclosures (see DIFFERS above).")
print("PASS: all enclosures computed")   # the comparison itself is information, not a claim of ours
sys.exit(0)
