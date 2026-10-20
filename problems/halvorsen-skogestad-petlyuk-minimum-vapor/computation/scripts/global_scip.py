"""Global minimisation (SCIP spatial branch-and-bound) of the total reboiler vapor of ONE configuration at fixed
parameters, to cross-check the multistart ipopt results (nlp.py).

Reformulation: u_{c,r,p} = 1/(alpha_p - theta_{c,r}) (bilinear constraint u*(alpha_p - theta) = 1), so the Underwood
feed equations and minimum-vapor inequalities become bilinear.  Roots are kept at relative distance >= eps from
the volatilities (eps = 1e-6 of the interval width), as in the literature (TG use 1e-4); vapor flows are bounded by
vmax.  SCIP's gap tolerance gives a certified lower bound only up to SCIP's floating-point feasibility tolerances
(not an exact certificate).

Environment: WSL Fedora, Python 3, pyscipopt 6.3 (SCIP), 1 thread.
Runtime: 0.05 s to a few minutes per configuration (N <= 4).
"""
from __future__ import annotations

import sys
import time

from pyscipopt import Model, quicksum

from nlp import structure


def scip_min(cfg, alpha, F, q, vmax=None, eps=1e-6, gap=1e-7, timelimit=600, verbose=False, cutoff=None,
             permissive=False):
    N = cfg['N']
    al = {p: float(alpha[p - 1]) for p in range(1, N + 1)}
    Fd = {p: float(F[p - 1]) for p in range(1, N + 1)}
    Ft = sum(Fd.values())
    cols, info, reb, side = structure(cfg)
    if vmax is None:
        vmax = 50.0 * Ft
    m = Model()
    if not verbose:
        m.hideOutput()
    m.setParam('limits/gap', gap)
    m.setParam('limits/time', timelimit)
    m.setParam('parallel/maxnthreads', 1)
    m.setParam('numerics/feastol', 1e-9)
    d, b, VT, VB, th, u = {}, {}, {}, {}, {}, {}
    for s in cols:
        I = info[s]
        for p in range(I['i'], I['k'] + 1):
            d[s, p] = m.addVar(lb=0, ub=Fd[p])
        for p in range(I['l'], I['j'] + 1):
            b[s, p] = m.addVar(lb=0, ub=Fd[p])
        VT[s] = m.addVar(lb=0, ub=vmax)
        VB[s] = m.addVar(lb=0, ub=vmax)
        for r in range(I['l'] - 1, I['k'] + 1):
            w = al[r] - al[r + 1]
            lo, hi = al[r + 1] + eps * w, al[r] - eps * w
            th[s, r] = m.addVar(lb=lo, ub=hi)
            for p in range(I['i'], I['j'] + 1):
                a1, a2 = al[p] - hi, al[p] - lo  # range of alpha_p - theta
                ulo, uhi = sorted([1.0 / a1, 1.0 / a2])
                u[s, r, p] = m.addVar(lb=ulo, ub=uhi)
                m.addCons(u[s, r, p] * (al[p] - th[s, r]) == 1)

    def feed(s, p):
        if s == (1, N):
            return Fd[p]
        I = info[s]
        e = 0
        if I['top'] is not None:
            e = e + d[I['top'], p]
        if I['bot'] is not None:
            e = e + b[I['bot'], p]
        return e

    def vfeed(s):
        if s == (1, N):
            return (1 - q) * Ft
        I = info[s]
        e = 0
        if I['top'] is not None:
            t = I['top']
            e = e + (quicksum(d[t, p] for p in range(info[t]['i'], info[t]['k'] + 1)) if I['kind'] == 'HX'
                     else VT[t])
        if I['bot'] is not None and I['kind'] != 'HX':
            e = e - VB[I['bot']]
        return e

    for s in cols:
        I = info[s]
        for p in range(I['i'], I['j'] + 1):
            rhs = 0
            if p <= I['k']:
                rhs = rhs + d[s, p]
            if p >= I['l']:
                rhs = rhs + b[s, p]
            m.addCons(feed(s, p) - rhs == 0)
        m.addCons(VT[s] - VB[s] - vfeed(s) == 0)
        tgt = VT[s]
        if permissive and I['l'] <= I['k'] - 1:
            tgt = m.addVar(lb=0, ub=vmax)
            m.addCons(VT[s] - tgt >= 0)
        for r in range(I['l'] - 1, I['k'] + 1):
            m.addCons(quicksum(al[p] * feed(s, p) * u[s, r, p] for p in range(I['i'], I['j'] + 1))
                      - vfeed(s) == 0)
            g = quicksum(al[p] * d[s, p] * u[s, r, p] for p in range(I['i'], I['k'] + 1))
            if I['l'] <= r <= I['k'] - 1:
                m.addCons(tgt - g == 0)
            else:
                m.addCons(tgt - g >= 0)
    for (t, uu) in side:
        m.addCons(VT[t] - VB[uu] == 0)
    obj = m.addVar(lb=0, ub=vmax * len(reb))
    m.addCons(obj - quicksum(VB[x] for x in reb) == 0)
    if cutoff is not None:
        m.addCons(obj <= cutoff)
    m.setObjective(obj, 'minimize')
    t0 = time.time()
    m.optimize()
    st = m.getStatus()
    try:
        primal = m.getObjVal() if m.getNSols() > 0 else None
    except Exception:  # noqa: BLE001
        primal = None
    return {'status': st, 'primal': primal, 'dual': m.getDualbound(), 'time': time.time() - t0,
            'nodes': m.getNNodes()}


if __name__ == '__main__':
    from hsmodel import all_configurations, config_name, petlyuk_value
    N = int(sys.argv[1])
    al = [float(x) for x in sys.argv[2].split(',')]
    z = [float(x) for x in sys.argv[3].split(',')]
    q = float(sys.argv[4])
    ids = [int(x) for x in sys.argv[5:]]
    vt, vp = petlyuk_value(al, z, q)
    print('Petlyuk reboiler vapor', float(vp))
    C = all_configurations(N)
    for i in ids:
        r = scip_min(C[i], al, z, q, vmax=20 * float(vt))
        print(i, config_name(C[i]), r, flush=True)
