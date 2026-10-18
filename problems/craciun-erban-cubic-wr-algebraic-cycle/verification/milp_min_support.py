#!/usr/bin/env python3
"""Verifier's side check (not needed for the theorem): smallest number of reactions of a network in W_3 whose
mass-action system has the circle h = (x-2)^2 + (y-2)^2 - 1 = 0 as an invariant curve free of equilibria.

MILP (HiGHS, 1 thread). Unknowns: rate k_r for each of the 90 reactions u -> v between the 10 complexes of
molecularity <= 3, binaries z_r (reaction present), circulation phi_r, cofactor s (6 coefficients, degree <= 2).
  * f h_x + g h_y = s h (15 linear equations);
  * z_r <= k_r <= KMAX z_r   (time rescaling allows min rate >= 1; KMAX bounds the rate ratio);
  * z_r <= phi_r <= PHIMAX z_r, flow conservation at every complex: a strictly positive circulation on the support,
    i.e. weak reversibility;
  * sign * P(theta_j) >= 1 at NS sample points of the circle (P = tangential speed; sign = +1 / -1 for
    counter-clockwise / clockwise orientation). This is a relaxation of "no equilibrium on the circle".
Minimise sum z_r. Hyperbolicity is NOT imposed, so the optimum is a lower bound for the number of reactions of a
network with this circle as a limit cycle (under the KMAX/PHIMAX bounds).
Runtime ~ 20 s (WSL Fedora). Environment: Python 3.12, highspy (HiGHS), sympy, numpy.
"""
import sys
import itertools
import numpy as np
import sympy as sp
import highspy

KMAX, PHIMAX, NS = 1.0e5, 1.0e3, 96
x, y = sp.symbols("x y")
cx = [(a, b) for a in range(4) for b in range(4) if a + b <= 3]
R = [(u, v) for u in cx for v in cx if u != v]
nR = len(R)
h = (x - 2)**2 + (y - 2)**2 - 1
hx, hy = sp.diff(h, x), sp.diff(h, y)
mons = [(i, j) for i in range(5) for j in range(5) if i + j <= 4]
smons = [(i, j) for i in range(3) for j in range(3) if i + j <= 2]


def coeff_vec(expr):
    p = sp.Poly(sp.expand(expr), x, y)
    d = dict(p.terms())
    return [float(d.get(m, 0)) for m in mons]


# columns: k (nR), z (nR), phi (nR), c (6)
cols_eq = []
for (u, v) in R:
    fr = (v[0] - u[0]) * x**u[0] * y**u[1]
    gr = (v[1] - u[1]) * x**u[0] * y**u[1]
    cols_eq.append(coeff_vec(fr * hx + gr * hy))
for (i, j) in smons:
    cols_eq.append(coeff_vec(-(x**i * y**j) * h))
Aeq = np.array(cols_eq).T          # rows: monomials; columns: k..., c...

thetas = np.linspace(0, 2 * np.pi, NS, endpoint=False)
Prow = []
for (u, v) in R:
    fr = (v[0] - u[0]) * x**u[0] * y**u[1]
    gr = (v[1] - u[1]) * x**u[0] * y**u[1]
    Pr = sp.lambdify((x, y), -(y - 2) * fr + (x - 2) * gr, "numpy")
    Prow.append(np.broadcast_to(Pr(2 + np.cos(thetas), 2 + np.sin(thetas)), thetas.shape).astype(float))
Prow = np.array(Prow).T            # NS x nR


def solve(sign):
    hs = highspy.Highs()
    hs.setOptionValue("output_flag", False)
    hs.setOptionValue("threads", 1)
    hs.setOptionValue("time_limit", 900.0)
    inf = highspy.kHighsInf
    nv = 3 * nR + len(smons)
    lo = [0.0] * nR + [0.0] * nR + [0.0] * nR + [-inf] * len(smons)
    up = [KMAX] * nR + [1.0] * nR + [PHIMAX] * nR + [inf] * len(smons)
    cost = [0.0] * nR + [1.0] * nR + [0.0] * nR + [0.0] * len(smons)
    for i in range(nv):
        hs.addVar(lo[i], up[i])
        hs.changeColCost(i, cost[i])
    for i in range(nR, 2 * nR):
        hs.changeColIntegrality(i, highspy.HighsVarType.kInteger)

    def row(idx, vals, lb, ub):
        hs.addRow(lb, ub, len(idx), np.array(idx, dtype=np.int32), np.array(vals, dtype=float))

    kc = list(range(nR)) + list(range(3 * nR, nv))
    for m in range(Aeq.shape[0]):
        row(kc, Aeq[m], 0.0, 0.0)
    for r in range(nR):
        row([r, nR + r], [1.0, -1.0], 0.0, inf)            # k_r >= z_r
        row([r, nR + r], [1.0, -KMAX], -inf, 0.0)          # k_r <= KMAX z_r
        row([2 * nR + r, nR + r], [1.0, -1.0], 0.0, inf)   # phi_r >= z_r
        row([2 * nR + r, nR + r], [1.0, -PHIMAX], -inf, 0.0)
    for c in cx:
        idx, vals = [], []
        for r, (u, v) in enumerate(R):
            if v == c:
                idx.append(2 * nR + r); vals.append(1.0)
            if u == c:
                idx.append(2 * nR + r); vals.append(-1.0)
        row(idx, vals, 0.0, 0.0)
    for j in range(NS):
        row(list(range(nR)), sign * Prow[j], 1.0, inf)
    hs.run()
    st = hs.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        return str(st), None, None
    sol = np.array(hs.getSolution().col_value)
    sup = [R[r] for r in range(nR) if sol[nR + r] > 0.5]
    return "optimal", round(hs.getInfo().objective_function_value), [(u, v, sol[r]) for r in range(nR)
                                                                    for (u2, v2) in [R[r]] if (u2, v2) in sup
                                                                    for u, v in [(u2, v2)]]


for sign in (+1, -1):
    st, obj, sup = solve(sign)
    print(f"orientation {'counter-clockwise' if sign > 0 else 'clockwise'}: status {st}, min #reactions = {obj}")
    if sup:
        for u, v, k in sup:
            print(f"   {u} -> {v}   k = {k:.6g}")
