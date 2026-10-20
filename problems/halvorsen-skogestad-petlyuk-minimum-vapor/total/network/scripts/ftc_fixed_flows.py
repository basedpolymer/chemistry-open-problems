#!/usr/bin/env python3
"""Finite-stage FTC (generalized Petlyuk) at FIXED flows F_delta, increasing stage numbers (Agent N, attainment).

Flows: result 12's attainment point (every column [i,m] with m - i >= 2 at its preferred split, binary columns with
excess e_r = max H - H_r; theory/scripts/check_ftc_achievability.py), plus a vapor/liquid circulation of total size
delta from the reboiler to the condenser spread over all monotone lattice paths of the FTC grid (adds the same
amount to V and L of a section, so every net flow D_s = V_s - L_s and every junction balance is unchanged, every
section gets a positive margin, and R = max_j H_j + delta exactly).  With these flows FIXED, the finite-stage steady
state is computed (own sparse Newton solver `stagenet.py`) for n = n1 < n2 < ... stages per section; the script
prints the product impurities and R.  Claim tested (network.md §10, step (C)): the impurities tend to 0 as n grows
while R stays max H + delta.  Floating point; evidence only.

Usage: python ftc_fixed_flows.py N delta n1,n2,... [seed]
Environment: Python 3, numpy, scipy, mpmath.  Runtime: seconds to minutes.
"""
import sys
import itertools
import numpy as np
import mpmath as mp
from stagenet import Net, solve_lse

mp.mp.dps = 40


def attainment_point(al, z, q):
    """Columns (i,j) 0-based -> dict with d, b, VT, VB (result 12, Theorem achieve)."""
    N = len(al)
    al = [mp.mpf(a) for a in al]
    z = [mp.mpf(x) for x in z]
    q = mp.mpf(q)
    th = []
    for r in range(N - 1):
        f = lambda t: sum(a * x / (a - t) for a, x in zip(al, z)) - (1 - q)
        a_, b_ = al[r + 1], al[r]
        for _ in range(300):
            m = (a_ + b_) / 2
            if f(m) > 0:
                b_ = m
            else:
                a_ = m
        th.append((a_ + b_) / 2)
    H = [sum(al[p] * z[p] / (al[p] - th[r]) for p in range(r + 1)) - (1 - q) for r in range(N - 1)]
    C = max(H)
    d, b, VT, VB = {}, {}, {}, {}
    cols = sorted([(i, j) for i in range(N) for j in range(i + 1, N)], key=lambda c: -(c[1] - c[0]))
    for (i, j) in cols:
        if (i, j) == (0, N - 1):
            f = {p: z[p] for p in range(N)}
            u = 1 - q
        else:
            f = {p: mp.mpf(0) for p in range(i, j + 1)}
            u = mp.mpf(0)
            if j + 1 <= N - 1:
                for p in range(i, j + 1):
                    f[p] += d[(i, j + 1)].get(p, 0)
                u += VT[(i, j + 1)]
            if i - 1 >= 0:
                for p in range(i, j + 1):
                    f[p] += b[(i - 1, j)].get(p, 0)
                u -= VB[(i - 1, j)]
        if j - i >= 2:
            def coef(p):
                num = mp.fprod([al[p] - th[r] for r in range(i, j)])
                den = mp.fprod([al[p] - al[pp] for pp in range(i, j) if pp != p])
                return num / den / al[p]
            vt = f[i] / coef(i)
            dd = {p: vt * coef(p) for p in range(i, j)}
            bb = {p: f[p] - dd.get(p, 0) for p in range(i + 1, j + 1)}
            vb = vt - u
        else:
            r = i
            dd = {i: f[i]}
            bb = {j: f[j]}
            yT = al[i] * f[i] / (al[i] - th[r])
            yB = al[j] * f[j] / (th[r] - al[j])
            e = C - H[r]
            vt, vb = yT + e, yB + e
        d[(i, j)], b[(i, j)], VT[(i, j)], VB[(i, j)] = dd, bb, vt, vb
    return d, b, VT, VB, [float(h) for h in H], float(C)


def build_ftc(al, z, q, n, delta, nstages_fn=None, forced=None):
    """forced = (eta_A, eta_C) (N = 3 only): distillate z_A - eta_A, bottoms z_C - eta_C, side product
    z_B + eta_A + eta_C (the middle product collects prescribed impurities; liquids of the main column adjusted)."""
    N = len(al)
    d, b, VT, VB, H, C = attainment_point(al, z, q)
    net = Net(al, z, q)
    J = {}
    for i in range(N):
        for j in range(i + 1, N):
            J[(i, j)] = net.add_stage(f"J{i}{j}")
    P = {p: net.add_stage(f"P{p}") for p in range(N)}
    # sections: list of (name, bottom unit, top unit, V, L)
    secs = []
    for (i, j) in J:
        top_unit = J[(i, j - 1)] if j - 1 > i else P[i]
        bot_unit = J[(i + 1, j)] if i + 1 < j else P[j]
        vt = float(VT[(i, j)]); vb = float(VB[(i, j)])
        Dt = float(sum(d[(i, j)].values())); Bb = float(sum(b[(i, j)].values()))
        secs.append([("T", i, j), J[(i, j)], top_unit, vt, vt - Dt])
        secs.append([("B", i, j), bot_unit, J[(i, j)], vb, vb + Bb])
    # circulation delta from reboiler P[N-1] to condenser P[0] over all monotone paths (bottom -> top)
    out_edges = {}
    for k, s in enumerate(secs):
        out_edges.setdefault(s[1], []).append(k)
    # count paths: number of paths from unit u to condenser
    memo = {}

    def npaths(u):
        if u == P[0]:
            return 1
        if u in memo:
            return memo[u]
        memo[u] = sum(npaths(secs[k][2]) for k in out_edges.get(u, []))
        return memo[u]
    tot = npaths(P[N - 1])
    # paths from reboiler to u
    memo2 = {P[N - 1]: 1}
    order = sorted(set(s[1] for s in secs) | set(s[2] for s in secs), key=lambda u: -npaths(u) if u != P[0] else 0)
    # topological order by repeated relaxation (DAG)
    into = {}
    for k, s in enumerate(secs):
        into.setdefault(s[2], []).append(k)

    def npaths_from(u):
        if u == P[N - 1]:
            return 1
        if u in memo2:
            return memo2[u]
        memo2[u] = sum(npaths_from(secs[k][1]) for k in into.get(u, []))
        return memo2[u]
    for k, s in enumerate(secs):
        w = npaths_from(s[1]) * npaths(s[2]) / tot
        s[3] += delta * w
        s[4] += delta * w
    prod = [float(v) for v in z]
    if forced is not None:
        assert N == 3
        eA, eC = forced
        for s in secs:
            if s[0] in (("T", 0, 1), ("B", 0, 1)):
                s[4] += eA
            if s[0] in (("T", 1, 2), ("B", 1, 2)):
                s[4] -= eC
        prod = [prod[0] - eA, prod[1] + eA + eC, prod[2] - eC]
    # build stages
    for k, s in enumerate(secs):
        name, bu, tu, V, L = s
        ns = n if nstages_fn is None else nstages_fn(name)
        st = [net.add_stage(f"{name}.{t}") for t in range(ns)]
        if ns == 0:
            net.add_vap(bu, tu, V); net.add_liq(tu, bu, L)
            continue
        net.add_vap(bu, st[0], V)
        for t in range(ns - 1):
            net.add_vap(st[t], st[t + 1], V)
            net.add_liq(st[t + 1], st[t], L)
        net.add_vap(st[-1], tu, V)
        net.add_liq(tu, st[-1], L)
        net.add_liq(st[0], bu, L)
    for p in range(N):
        net.add_liq(P[p], ('P', p), prod[p])
    net.add_feed(J[(0, N - 1)], 1.0)
    return net, P, H, C


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    delta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
    ns = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "5,10,20,40").split(",")]
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    if seed == 0 and N == 3:
        al, z, q = [4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0
    else:
        rng = np.random.default_rng(seed)
        al = list(np.cumprod(np.r_[1.0, rng.uniform(1.4, 3.0, N - 1)])[::-1])
        z = list(rng.dirichlet(np.ones(N) * 3))
        q = float(rng.uniform(0.0, 1.0))
    print(f"N = {N}, alpha = {np.round(al, 4)}, z = {np.round(z, 4)}, q = {q:.3f}, delta = {delta}")
    prev = None  # (names, log-compositions, n) of the previous solve, for continuation in n
    out = []
    for n in ns:
        net, P, H, C = build_ftc(al, z, q, n, delta)
        Lout, Vout, Lin, Vin, Fl, Fv, dlt, massres = net.flows()
        assert np.max(np.abs(massres)) < 1e-10 and np.min(Lout) > -1e-12
        R = float(np.sum(np.maximum(dlt, 0)))
        cands = []
        # (a) fresh solve (successive substitution + Newton)
        Xa, ra = net.solve(maxit=100)
        cands.append((ra, np.log(np.maximum(Xa, 1e-300))))
        # (b) continuation from the previous n (interpolated profiles), plain Newton
        if prev is not None:
            pn, pxi, pnst = prev
            idx = {nm: k for k, nm in enumerate(pn)}
            xi0 = np.zeros((net.nst, net.N))
            for k, nm in enumerate(net.names):
                if nm in idx:
                    xi0[k] = pxi[idx[nm]]
                else:
                    sec, t = nm.rsplit(".", 1)
                    t_old = min(pnst - 1, int(round(int(t) * (pnst - 1) / max(n - 1, 1))))
                    xi0[k] = pxi[idx[f"{sec}.{t_old}"]]
            X0 = np.exp(xi0)
            X0 = X0 / X0.sum(axis=1, keepdims=True)
            Xb, rb = net.solve(X0=X0, maxit=100)
            cands.append((rb, np.log(np.maximum(Xb, 1e-300))))
        # (c) log-sum-exp Newton polish of the best candidate (kept only if it improves and stays normalized)
        cands.sort(key=lambda t: t[0])
        best_r, best_xi = cands[0]
        xi_c, rel_c = solve_lse(net, best_xi)
        Xc = np.exp(xi_c)
        res_c = float(np.max(np.abs(net.residual(Xc.ravel()))))
        if np.all(np.isfinite(Xc)) and np.max(np.abs(Xc.sum(axis=1) - 1)) < 1e-9 and res_c <= best_r:
            xi, res = xi_c, res_c
        else:
            xi, res = best_xi, best_r
        prev = (list(net.names), xi, n)
        # impurity of product p = sum of the other mole fractions (computed in log space)
        imp = []
        for p in range(N):
            v = np.delete(xi[P[p]], p)
            m = v.max()
            imp.append(float(np.exp(m) * np.sum(np.exp(v - m))))
        out.append((n, res, imp))
        print(f"  n = {n:4d}: residual {res:.1e}; R = {R:.6f} (max H = {C:.6f}); "
              f"impurities of products 1..N: " + " ".join(f"{v:.2e}" for v in imp), flush=True)
    return out


if __name__ == "__main__":
    main()
