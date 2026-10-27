"""Exact checks of the explicit examples of ../notes.md (agent Q1, round 2).

Environment: Python 3.12 (WSL Fedora 44) with sympy 1.14 (symbolic identities); the rest is standard library and
exact rationals.  Every claim prints PASS/FAIL; exit status 1 on any FAIL.
"""
import os
import sys
import time
import platform

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "computation", "scripts")))
import ecore  # noqa: E402
from s1lib import classify, dual, show  # noqa: E402
from drainsets import drain_sets, drain_species_E  # noqa: E402

FAILS = []


def claim(name, ok):
    print(f"{name}: {'PASS' if ok else 'FAIL'}")
    if not ok:
        FAILS.append(name)


def ode(p, rx, kap, xs):
    f = [0] * p
    for j, (a, b) in enumerate(rx):
        mono = kap[j]
        for i in range(p):
            mono *= xs[i] ** a[i]
        for i in range(p):
            f[i] += (b[i] - a[i]) * mono
    return [sp.expand(e) for e in f]


def H(n, rx):
    N = ecore.CNet(n, rx)
    return N.is_conservative() and N.is_consistent() and N.rows_distinct()


def info(n, rx):
    N = ecore.CNet(n, rx)
    return dict(WR=N.weakly_reversible(), lc=N.linkage_classes(), delta=N.deficiency(), s=N.rank(), H=H(n, rx))


def types(p, drx):
    v, data = classify(p, drx)
    return v, sorted((f["sigma"], f["Z"], f["type"]) for f in data[2])


def main():
    t0 = time.time()
    print(f"# examples.py, Python {platform.python_version()}, sympy {sp.__version__}, {platform.system()}")
    k1, k2, k3, k4 = sp.symbols("k1 k2 k3 k4", positive=True)
    K1, K2, K3 = sp.symbols("kappa1 kappa2 kappa3", positive=True)

    # ---- Example X: N = 2X2 <-> X1+X2 (smallest WR network whose dual is non-persistent for every rate) -------------
    rx = [((0, 2), (1, 1)), ((1, 1), (0, 2))]
    I = info(2, rx)
    claim("X1 N = 2X2 <-> X1+X2 is WR, one linkage class, deficiency 0, s = 1, (H1)-(H3)",
          I == dict(WR=True, lc=1, delta=0, s=1, H=True))
    p, drx = dual(2, rx)
    claim("X2 N* = k2 -> k1, 2k1+k2 -> k1+2k2", drx == [((0, 1), (1, 0)), ((2, 1), (1, 2))])
    f = ode(p, drx, [K1, K2], [k1, k2])
    claim("X3 k2' = k2 (kappa2 k1^2 - kappa1) (identity)", sp.expand(f[1] - k2 * (K2 * k1 ** 2 - K1)) == 0)
    claim("X4 k1' + k2' = 0 (classes k1 + k2 = T)", sp.expand(f[0] + f[1]) == 0)
    # on k1 + k2 = T with kappa2 T^2 < kappa1: k2' <= -(kappa1 - kappa2 T^2) k2  (since 0 < k1 < T)
    T = sp.symbols("T", positive=True)
    claim("X5 for 0<k1<T: kappa2 k1^2 - kappa1 <= kappa2 T^2 - kappa1 (monotone in k1)",
          sp.simplify(sp.diff(K2 * k1 ** 2 - K1, k1) - 2 * K2 * k1) == 0)
    v, ty = types(p, drx)
    claim("X6 face classification: verdict NONE, lower face {k2} of type ESC", v == "NONE" and (1, (1,), "ESC") in ty)
    claim("X7 no drain species and no drain set in N* (the extinction is not of drain type)",
          not drain_species_E(p, drx) and not drain_sets(p, drx))

    # ---- Example M (agent G): 3A <-> 2A+B, dual has an OUT face (drain species) -------------------------------------
    rx = [((3, 0), (2, 1)), ((2, 1), (3, 0))]
    p, drx = dual(2, rx)
    v, ty = types(p, drx)
    claim("M1 dual of 3A <-> 2A+B: verdict NONE with an OUT face", v == "NONE" and any(t[2] == "OUT" for t in ty))
    claim("M2 the OUT face is the drain species of E (k2) and a drain set", drain_species_E(p, drx) == [1] and bool(drain_sets(p, drx)))

    # ---- Example P: a face with |Z| = 2 decides (WR, s = 1) ---------------------------------------------------------
    rx = [((2, 0, 1), (3, 2, 0)), ((3, 2, 0), (2, 0, 1)), ((0, 2, 1), (1, 4, 0)), ((1, 4, 0), (0, 2, 1))]
    I = info(3, rx)
    claim("P1 N = 2X1+X3 <-> 3X1+2X2, 2X2+X3 <-> X1+4X2: WR, s = 1, (H1)-(H3)", I["WR"] and I["s"] == 1 and I["H"])
    p, drx = dual(3, rx)
    v, ty = types(p, drx)
    sing = all(t[2] == "IN" for t in ty if len(t[1]) == 1)
    pair = [t for t in ty if len(t[1]) == 2 and t[2] != "IN"]
    claim("P2 all singleton faces of N* are IN, a two-element face is not, verdict NONE", sing and pair and v == "NONE")
    print("   non-IN faces:", pair)

    # ---- Example S: non-WR, dual persistent exactly for kappa1 > kappa2 (open set) ---------------------------------
    rx = [((1, 1, 0), (0, 2, 2)), ((1, 1, 2), (2, 0, 0))]
    I = info(3, rx)
    claim("S1 N = X1+X2 -> 2X2+2X3, X1+X2+2X3 -> 2X1: s = 1, (H1)-(H3), not WR", I["s"] == 1 and I["H"] and not I["WR"])
    p, drx = dual(3, rx)
    claim("S2 N* = k1+k2 -> 2k2, k1+k2 -> 2k1, 2k2 -> 2k1", drx == [((1, 1), (0, 2)), ((1, 1), (2, 0)), ((0, 2), (2, 0))])
    f = ode(p, drx, [K1, K2, K3], [k1, k2])
    claim("S3 k1' = k2 ((kappa2 - kappa1) k1 + 2 kappa3 k2) (identity)",
          sp.expand(f[0] - k2 * ((K2 - K1) * k1 + 2 * K3 * k2)) == 0)
    v, ty = types(p, drx)
    claim("S4 verdict SOME: face {k2} DOM, face {k1} IN", v == "SOME" and {(t[1], t[2]) for t in ty} == {((1,), "DOM"), ((0,), "IN")})

    # ---- Example B: non-WR, dual persistent only on the hypersurface kappa1 = kappa2 -------------------------------
    rx = [((1, 1), (0, 2)), ((1, 1), (2, 0))]
    p, drx = dual(2, rx)
    f = ode(p, drx, [K1, K2], [k1, k2])
    claim("B1 N = X1+X2 -> 2X2, X1+X2 -> 2X1 (H1)-(H3); N* has k1' = (kappa2 - kappa1) k1 k2",
          H(2, rx) and sp.expand(f[0] - (K2 - K1) * k1 * k2) == 0)
    v, ty = types(p, drx)
    claim("B2 verdict TWO-SIDED (DOM faces on both sides)", v == "TWO-SIDED")

    # ---- Example C: s = 2, molecularity 2, one linkage class, drain SET {k2, k3} but no drain species ------------
    rx = [((0, 0, 2), (0, 1, 1)), ((0, 1, 1), (1, 0, 1)), ((1, 0, 1), (0, 0, 2))]
    I = info(3, rx)
    claim("C1 N = 2X3 -> X2+X3 -> X1+X3 -> 2X3: WR, one linkage class, deficiency 0, s = 2, (H1)-(H3)",
          I == dict(WR=True, lc=1, delta=0, s=2, H=True))
    p, drx = dual(3, rx)
    claim("C2 N* = k3 -> k2, k2 -> k1, 2k1+k2+k3 -> k1+k2+2k3",
          drx == [((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((2, 1, 1), (1, 1, 2))])
    f = ode(p, drx, [K1, K2, K3], [k1, k2, k3])
    V = sp.expand(f[1] + 2 * f[2])
    claim("C3 V = k2 + 2 k3: V' = -kappa1 k3 - kappa2 k2 + 2 kappa3 k1^2 k2 k3 (identity)",
          sp.expand(V - (-K1 * k3 - K2 * k2 + 2 * K3 * k1 ** 2 * k2 * k3)) == 0)
    claim("C4 N* conservative (k1+k2+k3)", sp.expand(f[0] + f[1] + f[2]) == 0)
    ds = drain_sets(p, drx)
    claim("C5 drain set {k2, k3} found, no drain species", (1, 2) in [Z for Z, d, c in ds] and not drain_species_E(p, drx))

    # ---- Example Q: s = 2, N pure (not WR): N* non-persistent unless kappa2 = kappa3 --------------------------------
    rx = [((1, 0, 0), (0, 1, 0)), ((1, 0, 0), (0, 0, 1)), ((0, 1, 1), (2, 0, 0))]
    I = info(3, rx)
    pure = all(a[i] == 0 or b[i] == 0 for a, b in rx for i in range(3))
    claim("Q1 N = X1 -> X2, X1 -> X3, X2+X3 -> 2X1: pure, s = 2, (H1)-(H3), not WR", pure and I["s"] == 2 and I["H"] and not I["WR"])
    p, drx = dual(3, rx)
    claim("Q2 N* = k1+k2 -> 2k3, k3 -> k1, k3 -> k2 (pure too)",
          drx == [((1, 1, 0), (0, 0, 2)), ((0, 0, 1), (1, 0, 0)), ((0, 0, 1), (0, 1, 0))])
    f = ode(p, drx, [K1, K2, K3], [k1, k2, k3])
    claim("Q3 (k2 - k1)' = (kappa3 - kappa2) k3 and k1+k2+k3 conserved (identities)",
          sp.expand(f[1] - f[0] - (K3 - K2) * k3) == 0 and sp.expand(f[0] + f[1] + f[2]) == 0)

    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
