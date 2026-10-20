#!/usr/bin/env python3
"""Structure of finite-stage ternary Petlyuk (FTC, N = 3) steady states at fixed flows F_delta, as n grows
(Agent N, attainment, network.md §11): junction compositions, net component flows of every section, product
impurities.  Floating point (log-sum-exp Newton), exploration only.
Usage: python ftc3_structure.py delta n1,n2,...   Environment: Python 3, numpy, scipy, mpmath.
"""
import sys
import numpy as np
from ftc_fixed_flows import build_ftc
from stagenet import solve_lse

delta = float(sys.argv[1]); ns = [int(v) for v in sys.argv[2].split(",")]
al, z, q = [4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0
prev = None
for n in ns:
    net, P, H, C = build_ftc(al, z, q, n, delta)
    Xa, ra = net.solve(maxit=100)
    xi = np.log(np.maximum(Xa, 1e-300))
    if prev is not None:
        pn, pxi, pnst = prev
        idx = {nm: k for k, nm in enumerate(pn)}
        xi0 = np.zeros((net.nst, 3))
        for k, nm in enumerate(net.names):
            if nm in idx:
                xi0[k] = pxi[idx[nm]]
            else:
                sec, t = nm.rsplit(".", 1)
                t_old = min(pnst - 1, int(round(int(t) * (pnst - 1) / max(n - 1, 1))))
                xi0[k] = pxi[idx[f"{sec}.{t_old}"]]
        X0 = np.exp(xi0); X0 /= X0.sum(axis=1, keepdims=True)
        Xb, rb = net.solve(X0=X0, maxit=100)
        if rb < ra:
            xi, ra = np.log(np.maximum(Xb, 1e-300)), rb
    xi2, rel = solve_lse(net, xi)
    X = np.exp(xi2)
    if np.max(np.abs(X.sum(axis=1) - 1)) < 1e-9:
        xi = xi2
    prev = (list(net.names), xi, n)
    X = np.exp(xi)
    a = np.array(al)
    def nm(s): return net.names.index(s)
    J1, J3, F, J2 = X[nm('J01')], X[nm('J12')], X[nm('J02')], X[nm('P1')]
    xD, xR = X[nm('P0')], X[nm('P2')]
    # net upward flows of sections: w = V y(x_bottom-end) - L x_(stage 1) ... use stage 0 of each section
    def netflow(sec, bottom_unit):
        V = [f for d, f in net.vap[bottom_unit] if d == nm(f"{sec}.0")][0]
        L = [f for d, f in net.liq[nm(f"{sec}.0")] if d == bottom_unit][0]
        xb = X[bottom_unit]; y = a * xb / (a @ xb)
        return V * y - L * X[nm(f"{sec}.0")]
    wPT = netflow("('T', 0, 2)", nm('J02'))
    wPB = netflow("('B', 0, 2)", nm('J12'))
    wM1 = netflow("('T', 0, 1)", nm('J01'))
    wM2 = netflow("('B', 0, 1)", nm('P1'))
    wM3 = netflow("('T', 1, 2)", nm('J12'))
    wM4 = netflow("('B', 1, 2)", nm('P2'))
    print(f"n={n} res={ra:.1e}")
    print(f"  J1 = {J1}  xi=A/(A+B)={J1[0]/(J1[0]+J1[1]):.6f}")
    print(f"  J3 = {J3}  zeta=B/(B+C)={J3[1]/(J3[1]+J3[2]):.6f}")
    print(f"  F  = {F}")
    print(f"  J2 = {J2};  xD = {xD};  xReb = {xR}")
    print(f"  wPT = {wPT}  (eps_C = {wPT[2]:.3e});  wPB = {wPB} (eps_A = {-wPB[0]:.3e})")
    print(f"  wM1 = {wM1}; wM2 = {wM2}")
    print(f"  wM3 = {wM3}; wM4 = {wM4}", flush=True)
