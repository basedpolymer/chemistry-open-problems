#!/usr/bin/env python3
"""crosscheck_L.py -- Agent C.  Cross-validation of two independent implementations of the exact finite-stage system.

My finite-stage steady states (petnet.py: unknowns = all unit compositions and all section net flows, product form
with polynomial left eigenvectors; profiles re-checked stage by stage) are mapped into Agent L's variables (net
splits d_c, junction compositions X_c; limit/scripts/limitsys.py, ratio form with L's own eigen-solver and trace
normalization), and L's residual G_n is evaluated there with L's code (read-only import).  L's Design object is
filled from MY operating specification at full precision.  Also compares L's recomputed limit configuration v_inf
with mine.  Expected: L's |G_n(my v_n)| at the level of the working precision; v_inf identical.
Cases: N = 3 (BT, network.md 11.3 data), N = 4 (BBB), N = 5 (BBBB): n = 24, 48; N = 6 (BTBTB, mixed data,
rho = 0.90, pre-asymptotic at n = 24): n = 48, 72.
Prints PASS/FAIL.  Numerical evidence only.
Usage: python crosscheck_L.py
Environment: Python 3.12, mpmath 1.3 (WSL Fedora sci env).  Runtime: a few minutes.
"""
import os
import sys
import time
import mpmath as mp
from petnet import Problem, Design, limit_config, solve_audited

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'limit', 'scripts'))
import limitsys as ls  # noqa: E402

t0 = time.time()
fails = []
CASES = [
    ('N3 BT', [4, 2, 1], [1, 1, 1], 1, '0.05', '0.025', 'BT', '0.003'),
    ('N4 BBB', [8, 4, 2, 1], [1, 1, 1, 1], 1, '0.05', '0.005', 'BBB', '0.003'),
    ('N5 BBBB', [16, 8, 4, 2, 1], [1] * 5, '0.5', '0.05', '0.002', 'BBBB', '0.002'),
    ('N6 BTBTB', [7, 5, 3.4, 2.3, 1.5, 1], [0.2, 0.1, 0.25, 0.15, 0.1, 0.2], '0.8', '0.03', '0.002', 'BTBTB', '0.001'),
]


def to_L(des):
    """L's Design object filled from my specification."""
    pb = des.pb
    N = pb.N
    D = ls.Design()
    D.N, D.al, D.z, D.q, D.H, D.R = N, list(pb.al), list(pb.z), pb.q, list(pb.H), des.R
    D.pattern = {p: des.col[(p, p + 1)]['kind'] for p in range(N - 1)}
    D.cols = sorted([(i, m) for i in range(N) for m in range(i + 1, N)], key=lambda c: (-(c[1] - c[0]), c[0]))
    D.d, D.b, D.f, D.VT, D.VB, D.LT, D.LB, D.chi = {}, {}, {}, {}, {}, {}, {}, {}
    for c in D.cols:
        r = des.col[c]
        D.d[c], D.b[c], D.f[c] = list(r['d']), list(r['b']), list(r['f'])
        D.VT[c], D.VB[c], D.LT[c], D.LB[c] = r['VT'], r['VB'], r['LT'], r['LB']
        if c[1] - c[0] >= 2:
            D.chi[c] = [r['chi'][j] for j in range(c[0], c[1])]
    D.P = list(des.P)
    return D


for name, al, z, q, delta, enb, pat, eb in CASES:
    N = len(al)

    def build(al=al, z=z, q=q, delta=delta, enb=enb, pat=pat, eb=eb):
        pb = Problem(al, z, mp.mpf(q))
        des = Design(pb, mp.mpf(delta), eta_nb_default=mp.mpf(enb),
                     eta_bin={(p, p + 1): (pat[p], mp.mpf(eb)) for p in range(len(al) - 1)})
        return pb, des, limit_config(des)
    for n in ([48, 72] if N == 6 else [24, 48]):
        res = solve_audited(build, n, 80)
        if not res['phys']:
            print(f"[FAIL] {name} n = {n}: my solution not physical")
            fails.append(f"{name} n={n}")
            continue
        des, sn, v = res['des'], res['sn'], res['v']
        X, W = sn.unpack(v)
        D = to_L(des)
        SL = ls.System(D, n)
        # map my solution into L's variables
        vals = []
        for (c, typ, k, kind, comp) in SL.vars:
            if typ == 'd':
                vals.append(W[('T',) + c][k])
            elif typ == 'b':
                vals.append(-W[('B',) + c][k])
            else:
                vals.append(X[c][k])
        G = SL.G(vals)
        gmax = max(abs(g) for g in G)
        # limit configurations
        SLi = ls.System(D, None)
        lim = res['lim']
        dX = max(max(abs(a - b) for a, b in zip(SLi.X_inf[c], lim['X'][c])) for c in D.cols)
        Gi = max(abs(g) for g in SLi.G(SLi.v_inf))
        ok = gmax < mp.mpf(10) ** (-(mp.mp.dps // 2)) and dX < mp.mpf(10) ** -40
        print(f"[{'PASS' if ok else 'FAIL'}] {name} n = {n} (dps {mp.mp.dps}): L's |G_n| at my steady state = "
              f"{mp.nstr(gmax, 3)} ({len(G)} equations, {len(SL.vars)} unknowns); L's v_inf vs mine: max |dX| = "
              f"{mp.nstr(dX, 3)}; L's |G_inf(L's v_inf)| = {mp.nstr(Gi, 3)}", flush=True)
        if not ok:
            fails.append(f"{name} n={n}")
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
