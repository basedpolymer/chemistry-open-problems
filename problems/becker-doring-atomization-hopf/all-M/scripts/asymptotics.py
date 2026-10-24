#!/usr/bin/env python3
"""asymptotics.py -- the limit M -> infinity of l1(q1=1)/M on branches j = 1, 2, 3 (rigorous enclosures).

Environment: Python 3.13, python-flint 0.9, mpmath.  Runtime: a few seconds.
Claims (PASS/FAIL):
  A1  at eps = 0 the Hopf system is Q(z;kappa) = 0, Re z = 0; its zero near (kappa_j^0, i t_j) is enclosed by the
      Krawczyk test (tan t_j = t_j, cos t_j < 0, kappa_j^0 = (sqrt(1+t_j^2) - 1)^(1/2), PV (89)-(91));
  A2  the eps-regular formula (5.2) of result.md, evaluated at eps = 0, gives Lhat(0) = lim l1(q1=1)/M (continuity of
      the formula on [0, delta], result.md Sec. 6), enclosed to ~1e-50;
  A3  (comparison, not a proof of an identity) result 8's closed-form limit coefficient
      L_j = -R(R-1)(7R^2-34R+23)/[t(R+1)(5R-3)((R-1)^3+R+1)], R = sqrt(1+t^2), t = t_j,
      lies in the enclosure of A2 (agreement to the printed number of digits);
  A4  uniform enclosure of Lhat on the box eps in [0, 1e-5] (all M >= 1e10) by the mean value form, and Lhat < 0.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from flint import arb, acb, ctx
import allm_core as C
import l1reg as R
import h2_lib as H

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

ctx.prec = 256
mp.mp.dps = 80
for j in (1, 2, 3):
    c = (2 * j - mp.mpf(1) / 2) * mp.pi
    t = mp.findroot(lambda u: mp.tan(u) - u, c - 1 / c)
    k0 = mp.sqrt(mp.sqrt(1 + t * t) - 1)
    X0 = [arb(mp.nstr(k0, 70)), arb(0), arb(mp.nstr(t, 70))]
    enc, box = C.certify_hopf(X0, arb(0))
    okA1 = enc is not None and bool(mp.cos(t) < 0)
    report("A1", okA1, "j=%d: kappa_j^0 in %s, t_j in %s" % (j, enc[0].str(30), enc[2].str(30)) if okA1 else "")
    out = R.blocks(arb(0), enc[0], acb(enc[1], enc[2]))
    L0 = out["Lhat"]
    okA2 = L0 is not None and L0.real.rad() < arb("1e-40")
    report("A2", okA2, "j=%d: Lhat(0) = %s" % (j, acb(L0).real.str(40)))
    # result 8's closed form, in ball arithmetic from the enclosure of t_j
    tb = enc[2]
    Rb = (1 + tb * tb).sqrt()
    Lj = -Rb * (Rb - 1) * (7 * Rb ** 2 - 34 * Rb + 23) / (tb * (Rb + 1) * (5 * Rb - 3) * ((Rb - 1) ** 3 + Rb + 1))
    diff = acb(L0).real - Lj
    okA3 = bool(diff.contains(0))
    report("A3", okA3, "j=%d: L_j (result 8) = %s, difference enclosure %s" % (j, Lj.str(40), diff.str(5)))
    # A4: uniform enclosure for eps in [0, 1e-5]
    eb, Xm, rad, enc4, box4 = H.hopf_box(0, "1e-5", [arb(enc[0].mid()), arb(0), arb(enc[2].mid())])
    L = R.Lhat_box(eb, box4[0], box4[1], box4[2]) if enc4 is not None else None
    okA4 = L is not None and bool(L < 0)
    report("A4", okA4, "j=%d: Lhat in [%s, %s] for all M >= 1e10" % (j, arb(L.lower()).str(8), arb(L.upper()).str(8)) if L is not None else "")
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
