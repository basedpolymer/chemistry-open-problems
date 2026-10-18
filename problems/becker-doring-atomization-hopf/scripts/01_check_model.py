#!/usr/bin/env python3
"""01_check_model.py -- model, linearization and reproduction of the paper's numbers.

Environment: WSL Fedora, Python 3.12, numpy/scipy, python-flint 0.9.  Runtime: about 1 minute.
Claims (PASS/FAIL):
  M1  the Jacobian of (2)-(4) at n = A = 1+K equals the paper's matrix B (27) (finite differences);
  M2  f(A + v) = B v + v_1 S v exactly (quadratic, no cubic term), for random v;
  M3  l = (1..M) is a left null vector of B and of S (mass conservation of the linear and quadratic parts);
  M4  eigenvalues for M = 100, K = 3 quoted in the paper: 0.05836 +- 0.2014i, 0.02585 +- 0.3618i, -410.94;
  M5  Table 1 of the paper (kappa_1 = K sqrt(M), Im lambda) for M = 1e2 .. 1e6, from certified Hopf enclosures;
  M6  the closed-form (ball) value of l1 agrees with the dense floating-point Kuznetsov formula (M = 20..200).
"""
import sys, math
import numpy as np
from flint import arb, acb, ctx
sys.path.insert(0, __import__("os").path.dirname(__file__))
import bdlib as bd

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-3s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

rng = np.random.default_rng(12345)
# M1
errs = []
for (M, K) in [(7, 0.3), (12, 0.7), (25, 2.5)]:
    A = 1 + K; n0 = A*np.ones(M); h = 1e-6
    J = np.array([(bd.field(n0 + h*e, K) - bd.field(n0 - h*e, K))/(2*h) for e in np.eye(M)]).T
    errs.append(np.abs(J - bd.Bmat(K, M)).max()/np.abs(bd.Bmat(K, M)).max())
report("M1", max(errs) < 1e-8, "max rel. error %.1e" % max(errs))
# M2
errs = []
for (M, K) in [(7, 0.3), (12, 0.7), (25, 2.5)]:
    A = 1 + K; v = rng.normal(size=M)
    r = bd.field(A + v, K) - bd.Bmat(K, M)@v - v[0]*(bd.Smat(M)@v)
    errs.append(np.abs(r).max())
report("M2", max(errs) < 1e-11, "max residual %.1e" % max(errs))
# M3
errs = []
for (M, K) in [(7, 0.3), (12, 0.7), (25, 2.5)]:
    l = np.arange(1, M + 1)
    errs.append(max(np.abs(l@bd.Bmat(K, M)).max(), np.abs(l@bd.Smat(M)).max()))
report("M3", max(errs) < 1e-10, "max |l B|, |l S| = %.1e" % max(errs))
# M4
w = np.linalg.eigvals(bd.Bmat(3.0, 100))
top = sorted([x for x in w if x.real > 0 and x.imag > 0], key=lambda x: -x.real)
ok4 = (len(top) == 2 and abs(top[0] - (0.05836 + 0.2014j)) < 1e-4 and abs(top[1] - (0.02585 + 0.3618j)) < 1e-4
       and abs(w.real.min() + 410.94) < 0.01)
report("M4", ok4, "unstable: %s ; min eigenvalue %.2f" % (["%.5f%+.4fi" % (x.real, x.imag) for x in top], w.real.min()))
# M5
table = {10**2: (3.9349, 0.021740), 10**3: (2.3722, 3.6176e-4), 10**4: (2.0376, 9.3596e-6),
         10**5: (1.9414, 2.7777e-7), 10**6: (1.9118, 8.6091e-9)}
ok5 = True
X = None
for M in sorted(table):
    if X is None:
        cr = [c for c in bd.dense_crossings(M) if c[3] == 1][0]
        K, lam, phi = bd.dense_hopf(M, cr[0], cr[1], 0)
        z = M*(phi - 1); X = [arb(K*math.sqrt(M)), arb(z.real), arb(z.imag)]
    res = bd.certify_point(X, M)
    X = [res["kappa"].mid(), res["x"].mid(), res["y"].mid()]
    kap = float(res["kappa"].mid()); om = float(res["omega"].mid())
    k_t, o_t = table[M]
    good = abs(kap - k_t) < 1.5e-4*k_t and abs(om - o_t) < 1.5e-4*o_t
    ok5 &= good
    print("     M=%.0e kappa_1 = %.6f (paper %.4f)  Im lambda = %.6e (paper %.4e)  %s" % (M, kap, k_t, om, o_t, good))
report("M5", ok5, "Table 1 reproduced")
# M6
ok6 = True
for M in [20, 50, 100, 200]:
    for (a, b, p, c) in bd.dense_crossings(M):
        K, lam, phi = bd.dense_hopf(M, a, b, p)
        _, l1d, _ = bd.lyap_dense(K, M, lam)
        ctx.prec = 192
        l1c = bd.lyapunov_cf(M, K, acb(phi.real, phi.imag))
        rel = abs(float(l1c.mid()) - l1d)/abs(l1d)
        ok6 &= rel < 1e-7
        print("     M=%d j=%d  l1 dense = %.10e   closed form = %.10e   rel.diff %.1e" % (M, c, l1d, float(l1c.mid()), rel))
report("M6", ok6, "closed form = dense formula")
print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
