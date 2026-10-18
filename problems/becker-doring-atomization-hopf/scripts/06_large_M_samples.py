#!/usr/bin/env python3
"""06_large_M_samples.py -- certified l1 at sample values of M up to 1e8 and convergence to the limit L_j.

Environment: WSL Fedora, Python 3.12, python-flint 0.9, mpmath.  Runtime: about 1-3 minutes.

For j = 1, 2, 3 and M in {1e4, 3e4, 1e5, ..., 1e8}: Krawczyk enclosure of the Hopf point, ball enclosure of l1.
Claims (PASS/FAIL):
  S1  l1 < 0 (certified) at every sample point;
  S2  (numerical observation, not a proof) the normalized coefficient  l1(q_1 = 1)/M  approaches the limit value
      L_j = -R(R-1)(7R^2-34R+23) / [t(R+1)(5R-3)((R-1)^3+R+1)]  (t = t_j, R = sqrt(1+t^2)) of the formal limit
      problem; the deviation at M = 1e8 is below 1% and decreases like M^(-1/2) (fitted exponent printed).
"""
import sys, math
import numpy as np
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, __import__("os").path.dirname(__file__))
import bdlib as bd

mp.mp.dps = 30
def limit_data(j):
    c = (2*j - mp.mpf(1)/2)*mp.pi
    t = mp.findroot(lambda u: mp.tan(u) - u, c - 1/c)
    R = mp.sqrt(1 + t**2)
    L = -R*(R - 1)*(7*R**2 - 34*R + 23)/(t*(R + 1)*(5*R - 3)*((R - 1)**3 + R + 1))
    return float(t), float(mp.sqrt(R - 1)), float(L)

Ms = [10**4, 3*10**4, 10**5, 3*10**5, 10**6, 3*10**6, 10**7, 3*10**7, 10**8]
fails = 0
for j in (1, 2, 3):
    t, k0, L = limit_data(j)
    # starting guess from the limit, refined by continuation from a moderate M
    X = [arb(k0*(1 + 0.6/math.sqrt(1e4)*j)), arb(0.1), arb(t)]
    devs = []
    print("j=%d  limit: t_j=%.8f kappa_j^0=%.8f L_j=%.10f" % (j, t, k0, L))
    for M in Ms:
        res = bd.certify_point(X, M, prec=192)
        if res is None:
            print("  M=%.0e  KRAWCZYK FAIL" % M); fails += 1; continue
        X = [res["kappa"].mid(), res["x"].mid(), res["y"].mid()]
        l1 = res["l1"]; v = float(res["l1q1"].mid())/M
        ok = bool(l1 < 0) and res["simple"] and bool(res["transv"] > 0)
        fails += (not ok)
        devs.append((M, abs(v/L - 1)))
        print("  M=%.0e kappa=%.10f z=%.7f%+.7fi  l1 in [%.6e, %.6e]  l1(q1=1)/M=%.10f  rel.dev=%.3e  %s"
              % (M, float(res["kappa"].mid()), float(res["x"].mid()), float(res["y"].mid()),
                 float(l1.lower()), float(l1.upper()), v, abs(v/L - 1), "PASS" if ok else "FAIL"), flush=True)
    Mv = np.array([d[0] for d in devs[-5:]], float); dv = np.array([d[1] for d in devs[-5:]])
    slope = np.polyfit(np.log(Mv), np.log(dv), 1)[0]
    s2 = devs[-1][1] < 0.01
    print("  S2 (observation): deviation at M=1e8 = %.3e (<1%%: %s), fitted exponent of the deviation in M: %.3f"
          % (devs[-1][1], s2, slope))
    fails += (not s2)
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
