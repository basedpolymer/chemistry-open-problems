#!/usr/bin/env python3
"""
table1.py -- rigorous enclosures of the two coefficients of Table 1 of Aminian-Camassa-McLaughlin
(Stud. Appl. Math. 141 (2018), p. 407; arXiv:1708.05097 v3 Table 1), for the regular n-gon with
circumradius 1 (the paper's normalization):
   LT coefficient  kappa_n = 3<u g2>/(2<u g1>)^(3/2) = 6<u g1^2>/(2<u g1>)^(3/2)   (scales like length)
   geometric skewness S^G_n = <u^3>/<u^2>^(3/2)                                   (scale invariant)
n = 3 is exact (triangle_exact.py); n >= 4 from lt_certify.certify / sg_certify.certify_sg with larger K.
Checks (PASS/FAIL) that each printed value of the paper is consistent with our enclosure up to the
paper's rounding (half a unit of the last printed digit).
Runtime: several minutes (one core).
"""
import sys, os, math
from flint import arb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lt_certify import certify
from sg_certify import certify_sg

# (n, printed LT coeff, printed SG, K for LT, K for SG); arXiv v3 values (published values are rounded versions)
ROWS = [(4, '0.139', '0.0813', 24, 40), (5, '0.194', '0.042', 24, 40), (6, '0.215', '0.024', 20, 40),
        (7, '0.225', '0.0149', 16, 36), (8, '0.231', '0.0097', 16, 32), (20, '0.243', '3.97e-4', 10, 30),
        (100, '0.245', '8.6e-7', 4, None)]
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
    kap = arb(r['kappa_circumradius1'])
    c1 = consistent(kap, lt)
    line = f"n={n:3d}: kappa(R=1) in {kap.str(5)} (paper {lt}) {'PASS' if c1 else 'FAIL'}"
    ok_all &= c1
    if Ks is not None:
        rs = certify_sg(n, K=Ks, verbose=False)
        sgv = arb(rs['SG'])
        c2 = consistent(sgv, sg)
        ok_all &= c2
        line += f";  S^G in {sgv.str(5)} (paper {sg}) {'PASS' if c2 else 'FAIL'}"
    print(line, flush=True)
print("PASS (all Table 1 entries consistent)" if ok_all else "FAIL (some entry inconsistent)")
sys.exit(0 if ok_all else 1)
