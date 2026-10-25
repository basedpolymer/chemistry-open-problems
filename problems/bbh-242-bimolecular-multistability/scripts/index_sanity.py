#!/usr/bin/env python3
"""Floating-point sanity check of Theorem thm:index (NOT part of the proof).

Random quadratic vector fields F = L z + q(z) with F1(0,y) >= 0, F2(x,0) >= 0 (all planar mass-action systems with
reactant molecularity <= 2 have this form).  (a) If L is Metzler and Hurwitz (origin a hyperbolic sink), the sum of
sign det DF over the positive equilibria must be <= 0 (Theorem thm:index).  (b) Control: with a saddle at the origin,
positive sums do occur, so test (a) is not vacuous.

Runtime ~1 min.  Environment: Python 3.12, numpy 2.x (WSL Fedora 44).  Prints PASS/FAIL; exit 1 on FAIL.
"""
import sys
import numpy as np
rng = np.random.default_rng(1)
P = np.polynomial.polynomial

def solve(cf, cg):
    f2 = np.array([cf[5]]); f1 = np.array([cf[2], cf[4]]); f0 = np.array([cf[0], cf[1], cf[3]])
    g2 = np.array([cg[5]]); g1 = np.array([cg[2], cg[4]]); g0 = np.array([cg[0], cg[1], cg[3]])
    a = P.polysub(P.polymul(f2, g0), P.polymul(g2, f0))
    b = P.polysub(P.polymul(f2, g1), P.polymul(g2, f1))
    c = P.polysub(P.polymul(f1, g0), P.polymul(g1, f0))
    R = P.polysub(P.polymul(a, a), P.polymul(b, c))
    xs = P.polyroots(R)
    sols = []
    for xr in xs:
        if abs(xr.imag) > 1e-9 or xr.real <= 1e-9: continue
        x = xr.real
        for (q2, q1, q0) in [(cf[5], cf[2] + cf[4] * x, cf[0] + cf[1] * x + cf[3] * x * x)]:
            ys = np.roots([q2, q1, q0]) if abs(q2) > 1e-14 else ([-q0 / q1] if abs(q1) > 1e-14 else [])
            for y in ys:
                if abs(np.imag(y)) > 1e-9 or np.real(y) <= 1e-9: continue
                y = np.real(y)
                gv = cg[0] + cg[1]*x + cg[2]*y + cg[3]*x*x + cg[4]*x*y + cg[5]*y*y
                T = abs(cg[0]) + abs(cg[1]*x) + abs(cg[2]*y) + abs(cg[3]*x*x) + abs(cg[4]*x*y) + abs(cg[5]*y*y)
                if abs(gv) < 1e-8 * T and not any(abs(x - s[0]) + abs(y - s[1]) < 1e-7 * (1 + x + y) for s in sols):
                    sols.append((x, y))
    return sols



def run(stable_origin, n):
    rng = np.random.default_rng(1 if stable_origin else 2)
    viol = 0; tested = 0; hist = {}
    for it in range(n):
        L = rng.normal(size=(2, 2)); L[0, 1] = abs(L[0, 1]); L[1, 0] = abs(L[1, 0])
        if stable_origin and not (np.trace(L) < 0 and np.linalg.det(L) > 0):
            continue
        if not stable_origin and not (np.linalg.det(L) < 0):
            continue
        q = rng.normal(size=(2, 3)) * 10.0 ** rng.uniform(-1, 1, size=(2, 3))
        q[0, 2] = abs(q[0, 2]); q[1, 0] = abs(q[1, 0])
        cf = np.array([0, L[0, 0], L[0, 1], q[0, 0], q[0, 1], q[0, 2]])
        cg = np.array([0, L[1, 0], L[1, 1], q[1, 0], q[1, 1], q[1, 2]])
        sols = solve(cf, cg)
        tested += 1
        S = 0
        for (x, y) in sols:
            J = np.array([[cf[1] + 2*cf[3]*x + cf[4]*y, cf[2] + cf[4]*x + 2*cf[5]*y],
                          [cg[1] + 2*cg[3]*x + cg[4]*y, cg[2] + cg[4]*x + 2*cg[5]*y]])
            S += int(np.sign(np.linalg.det(J)))
        hist[(len(sols), S)] = hist.get((len(sols), S), 0) + 1
        if S > 0:
            viol += 1
    return tested, viol, hist


if __name__ == "__main__":
    t, v, h = run(True, 400000)
    print("stable origin: fields tested", t, "| (number of positive equilibria, index sum) histogram", sorted(h.items()))
    ok1 = (v == 0 and t > 10000)
    print(("PASS" if ok1 else "FAIL"), "index sum <= 0 whenever the origin is a hyperbolic sink (%d violations)" % v)
    t2, v2, h2 = run(False, 400000)
    print("saddle origin (control): fields tested", t2, "| histogram", sorted(h2.items()))
    ok2 = v2 > 0
    print(("PASS" if ok2 else "FAIL"), "control: positive index sums occur when the origin is a saddle (%d cases)" % v2)
    sys.exit(0 if (ok1 and ok2) else 1)
