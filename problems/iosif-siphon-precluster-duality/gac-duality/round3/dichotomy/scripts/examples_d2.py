"""Agent D2 (round 3): exact symbolic checks of the hand-proved examples of notes.md (sympy).

Environment: Python 3.12 (WSL Fedora 44), sympy.  Usage (from round3/dichotomy/scripts/): python examples_d2.py
Claims:
  WRS  Proposition WR-Sigma on all 18 861 networks of ~/gacD2-data/s2_all.pkl.gz: Yp = Ye*Sigma, and the dual sources
       are pairwise distinct;
  U1   example U: the dual ODE, z' = a[(2k1 - k2) + K z (1 - 2 z)] with z = a/b, a = x2, b = x3, K = k3 x1^3;
  U2   example U: u' = a[(k2 - 2k1) u^2 - K u + 2K] and a' = a^2 (K - k2 u) with u = b/a;
  U3   example U: u* > K/k2  <=>  K < k2^2/k1 (positive root of (2k1 - k2)u^2 + K u - 2K, k2 < 2k1);
  CM   CM example: reduced coefficient on the center manifold at (0,T,0) equals k3 T (k2 T - k1)/(k1 + k2 T);
  L    L example: V = x2 + 2 x3 has V' = -2k1 x3^2 - 2k2 x2^2 + 2k3 x1^2 x2 x3 at the vertex (T,0,0) chart;
  ALL  smallest ALL example: x3'/x3 = k3 x1^2 - 2 k1 x3 (so V = x3 repels near (T,0,0)).
"""
import os
import sys
import gzip
import time
import pickle
import platform
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d2lib  # noqa: E402

fails = 0


def claim(name, ok, msg=""):
    global fails
    print(f"{'PASS' if ok else 'FAIL'}: {name} {msg}")
    fails += 0 if ok else 1


def dual_ode(n, rx, ks, xs):
    N = d2lib.Net(n, rx)
    f = [0] * N.r
    for i in range(N.n):
        rate = ks[i]
        for j in range(N.r):
            rate *= xs[j] ** N.Y[i][j]
        for j in range(N.r):
            f[j] += rate * N.G[i][j]
    return [sp.expand(e) for e in f]


def main():
    t0 = time.time()
    print(f"# examples_d2.py  Python {platform.python_version()} on {platform.system()}; sympy {sp.__version__}; exact")
    # WRS
    recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
    bad = 0
    for q in recs:
        n, rx = q["n"], q["rx"]
        srcs = [a for a, b in rx]
        if not all(b in srcs for a, b in rx):
            bad += 1
            continue
        N = d2lib.Net(n, rx)
        rows = [tuple(N.Y[i]) for i in N.changed]
        if len(set(rows)) != len(rows):
            bad += 1
    claim("WRS", bad == 0, f"(products are sources and dual sources distinct in {len(recs) - bad} of {len(recs)})")
    k1, k2, k3, T = sp.symbols("kappa1 kappa2 kappa3 T", positive=True)
    x1, x2, x3 = sp.symbols("x1 x2 x3", positive=True)
    # example U: 3X3 -> 2X1+X2, X2+2X3 -> 3X3, 2X1+X2 -> X2+2X3
    U = (((0, 0, 3), (2, 1, 0)), ((0, 1, 2), (0, 0, 3)), ((2, 1, 0), (0, 1, 2)))
    f = dual_ode(3, U, [k1, k2, k3], [x1, x2, x3])
    ok = (sp.simplify(f[1] - (-k2 * x2 * x3 + k3 * x1 ** 3 * x2 ** 2)) == 0 and
          sp.simplify(f[2] - (-2 * k1 * x3 ** 2 + 2 * k3 * x1 ** 3 * x2 ** 2)) == 0)
    K = k3 * x1 ** 3
    z = x2 / x3
    zp = f[1] / x3 - z * f[2] / x3
    ok = ok and sp.simplify(zp - x2 * ((2 * k1 - k2) + K * z * (1 - 2 * z))) == 0
    claim("U1", ok, "(z' identity)")
    u = x3 / x2
    up = f[2] / x2 - u * f[1] / x2
    ok = sp.simplify(up - x2 * ((k2 - 2 * k1) * u ** 2 - K * u + 2 * K)) == 0
    ok = ok and sp.simplify(f[1] - x2 ** 2 * (K - k2 * u)) == 0
    claim("U2", ok, "(u' and a' identities)")
    Ks, us = sp.symbols("K u", positive=True)
    ustar = (-Ks + sp.sqrt(Ks ** 2 + 8 * Ks * (2 * k1 - k2))) / (2 * (2 * k1 - k2))
    g = (2 * k1 - k2) * us ** 2 + Ks * us - 2 * Ks
    val = sp.simplify(g.subs(us, Ks / k2))
    ok = sp.simplify(val - Ks * (Ks * 2 * k1 / k2 ** 2 - 2)) == 0
    ok = ok and sp.simplify(g.subs(us, ustar)) == 0
    claim("U3", ok, "(g(K/k2) = K (2 k1 K / k2^2 - 2): u* > K/k2 iff K < k2^2/k1, g increasing on u > 0)")
    # CM example: 2X3 -> X2+X3, X2+X3 -> X1+X2, X1+X2 -> 2X3 ; vertex (0, T, 0), a = x1, b = x3
    C = (((0, 0, 2), (0, 1, 1)), ((0, 1, 1), (1, 1, 0)), ((1, 1, 0), (0, 0, 2)))
    f = dual_ode(3, C, [k1, k2, k3], [x1, x2, x3])
    a, b = sp.symbols("a b")
    sub = {x1: a, x3: b, x2: T - a - b}
    fa, fb = sp.expand(f[0].subs(sub)), sp.expand(f[2].subs(sub))
    beta = fa.coeff(b, 1).subs(a, 0)
    delta = fb.coeff(b, 1).subs(a, 0)
    f2 = fa.subs(b, 0).coeff(a, 2)
    g2 = fb.subs(b, 0).coeff(a, 2)
    ok = fa.subs({a: 0, b: 0}) == 0 and fa.coeff(a, 1).subs(b, 0) == 0 and fb.coeff(a, 1).subs(b, 0) == 0
    c = sp.simplify(f2 + beta * g2 / (-delta))
    ok = ok and sp.simplify(c - k3 * T * (k2 * T - k1) / (k1 + k2 * T)) == 0
    claim("CM", ok, f"(reduced coefficient {c})")
    # L example: 2X3 -> 2X2+X3, 2X2+X3 -> 2X1+X3, 2X1+X3 -> 2X3 ; V = x2 + 2 x3
    L = (((0, 0, 2), (0, 2, 1)), ((0, 2, 1), (2, 0, 1)), ((2, 0, 1), (0, 0, 2)))
    f = dual_ode(3, L, [k1, k2, k3], [x1, x2, x3])
    Vp = sp.expand(f[1] + 2 * f[2])
    ok = sp.simplify(Vp - (-2 * k1 * x3 ** 2 - 2 * k2 * x2 ** 2 + 2 * k3 * x1 ** 2 * x2 * x3)) == 0
    claim("L", ok, f"(V' = {Vp})")
    # smallest ALL example: 2X3 -> X2, X2 -> 2X1+X3, 2X1+X3 -> 2X3
    A = (((0, 0, 2), (0, 1, 0)), ((0, 1, 0), (2, 0, 1)), ((2, 0, 1), (0, 0, 2)))
    f = dual_ode(3, A, [k1, k2, k3], [x1, x2, x3])
    ok = sp.simplify(f[2] / x3 - (k3 * x1 ** 2 - 2 * k1 * x3)) == 0
    claim("ALL", ok, f"(x3' = {f[2]})")
    print(f"# runtime {time.time() - t0:.1f} s")
    print("ALL CLAIMS PASS" if fails == 0 else f"{fails} FAIL")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
