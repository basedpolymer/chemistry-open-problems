#!/usr/bin/env python3
"""Ternary Petlyuk (FTC, N = 3) at fixed flows with PRESCRIBED forced impurities of the side product (Agent N,
attainment, network.md §11): distillate z_A - eta_A, bottoms z_C - eta_C, side product z_B + eta_A + eta_C.
Reports, as n grows: product impurities, the junction ratios xi (A/(A+B) at J1) and zeta (B/(B+C) at J3), the net
flows of the main-column middle sections and of the prefractionator, and the smallest singular value of the
steady-state Jacobian (does the neutral mode disappear?).  Floating point; exploration only.
Usage: python ftc3_forced.py delta eta_A eta_C n1,n2,... [svd_max_n]
"""
import sys
import numpy as np
from ftc_fixed_flows import build_ftc
from stagenet import solve_lse

delta, eA, eC = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
ns = [int(v) for v in sys.argv[4].split(",")]
svd_max = int(sys.argv[5]) if len(sys.argv) > 5 else 40
al, z, q = [4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0
a = np.array(al)
prev = None
for n in ns:
    net, P, H, C = build_ftc(al, z, q, n, delta, forced=(eA, eC))
    Lout, Vout, Lin, Vin, Fl, Fv, dlt, massres = net.flows()
    assert np.max(np.abs(massres)) < 1e-10 and np.min(Lout) > 0
    R = float(np.sum(np.maximum(dlt, 0)))
    Xa, ra = net.solve(maxit=100)
    best = (ra, np.log(np.maximum(Xa, 1e-300)))
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
        if rb < best[0]:
            best = (rb, np.log(np.maximum(Xb, 1e-300)))
    xi, res = best[1], best[0]
    xi2, rel = solve_lse(net, xi)
    X2 = np.exp(xi2)
    r2 = float(np.max(np.abs(net.residual(X2.ravel()))))
    if np.max(np.abs(X2.sum(axis=1) - 1)) < 1e-9 and r2 <= res:
        xi, res = xi2, r2
    prev = (list(net.names), xi, n)
    X = np.exp(xi)
    def nm(s): return net.names.index(s)
    def netflow(sec, bu):
        V = [f for d, f in net.vap[bu] if d == nm(f"{sec}.0")][0]
        L = [f for d, f in net.liq[nm(f"{sec}.0")] if d == bu][0]
        y = a * X[bu] / (a @ X[bu])
        return V * y - L * X[nm(f"{sec}.0")]
    J1, J3 = X[nm('J01')], X[nm('J12')]
    imp = [1 - X[P[0]][0], 1 - X[P[1]][1], 1 - X[P[2]][2]]
    wPT = netflow("('T', 0, 2)", nm('J02')); wPB = netflow("('B', 0, 2)", nm('J12'))
    wM2 = netflow("('B', 0, 1)", nm('P1')); wM3 = netflow("('T', 1, 2)", nm('J12'))
    line = (f"n={n:4d} res={res:.1e} R={R:.5f} (maxH={C:.5f}) imp(D,S,B)=({imp[0]:.2e},{imp[1]:.3e},{imp[2]:.2e}) "
            f"xi={J1[0]/(J1[0]+J1[1]):.6f} zeta={J3[1]/(J3[1]+J3[2]):.6f} C@J1={J1[2]:.2e} A@J3={J3[0]:.2e} "
            f"epsC={wPT[2]:.2e} epsA={-wPB[0]:.2e} wM2={np.array2string(wM2, precision=2)} "
            f"wM3={np.array2string(wM3, precision=2)}")
    if n <= svd_max:
        Jm = net.jacobian(X.ravel()).toarray()
        s = np.linalg.svd(Jm, compute_uv=False)
        line += f" smin={s[-1]:.2e}"
    print(line, flush=True)
