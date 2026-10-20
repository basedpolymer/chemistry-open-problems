"""Adversarial search: for a configuration c and a peak index j, minimise jointly over the parameters
(alpha, z, q) and the operating variables the ratio  V_c / (P_j - (1-q)),  where V_c is the total reboiler
vapor of c and P_j the j-th peak of the feed's Vmin diagram (P_j - (1-q) <= Petlyuk reboiler vapor).
A value < 1 (beyond noise) would be a candidate counterexample.

Environment: WSL Fedora, Python 3, pyomo 6.10, ipopt 3.14.  Floating point exploration only.
Runtime: about 1-3 min per configuration (N = 4, 8 starts x 3 peaks), single thread.
Usage: python adversarial.py N qlo qhi starts seed outfile [config ids...]
"""
from __future__ import annotations

import json
import os
import math
import random
import sys
import time

import pyomo.environ as pyo
from pyomo.opt import TerminationCondition

from hsmodel import all_configurations, config_name, petlyuk_value
from nlp import build_model, random_start, solver, structure, max_violation


def build_joint(cfg, jpk, qlo, qhi, rmin=1.02, rmax=30.0, zmin=1e-3, dth=1e-7, permissive=False):
    N = cfg['N']
    cols, info, reb, side = structure(cfg)
    m = pyo.ConcreteModel()
    P = range(1, N + 1)
    m.al = pyo.Var(P, bounds=(1.0, rmax ** (N - 1)))
    m.z = pyo.Var(P, bounds=(zmin, 1.0))
    m.q = pyo.Var(bounds=(qlo, qhi))
    m.c = pyo.ConstraintList()
    m.c.add(m.al[N] == 1.0)
    for r in range(1, N):
        m.c.add(m.al[r] >= rmin * m.al[r + 1])
        m.c.add(m.al[r] <= rmax * m.al[r + 1])
    m.c.add(sum(m.z[p] for p in P) == 1.0)
    al, Fd, q = m.al, m.z, m.q
    m.d = pyo.Var([(s, p) for s in cols for p in range(info[s]['i'], info[s]['k'] + 1)], bounds=(0, 1))
    m.b = pyo.Var([(s, p) for s in cols for p in range(info[s]['l'], info[s]['j'] + 1)], bounds=(0, 1))
    m.VT = pyo.Var(cols, bounds=(0, 1e4))
    m.VB = pyo.Var(cols, bounds=(0, 1e4))
    th_idx = [(s, r) for s in cols for r in range(info[s]['l'] - 1, info[s]['k'] + 1)]
    m.th = pyo.Var(th_idx)
    m.Th = pyo.Var(range(1, N))  # feed roots

    def thc(var, r):
        m.c.add(var >= al[r + 1] + dth * (al[r] - al[r + 1]))
        m.c.add(var <= al[r] - dth * (al[r] - al[r + 1]))
    for (s, r) in th_idx:
        thc(m.th[s, r], r)
    for r in range(1, N):
        thc(m.Th[r], r)
        m.c.add(sum(al[p] * Fd[p] / (al[p] - m.Th[r]) for p in P) == 1 - q)

    def feed(s, p):
        if s == (1, N):
            return Fd[p]
        d = info[s]
        e = 0
        if d['top'] is not None:
            e = e + m.d[d['top'], p]
        if d['bot'] is not None:
            e = e + m.b[d['bot'], p]
        return e

    def vfeed(s):
        if s == (1, N):
            return 1 - q
        d = info[s]
        e = 0
        if d['top'] is not None:
            t = d['top']
            e = e + (sum(m.d[t, p] for p in range(info[t]['i'], info[t]['k'] + 1)) if d['kind'] == 'HX'
                     else m.VT[t])
        if d['bot'] is not None and d['kind'] != 'HX':
            e = e - m.VB[d['bot']]
        return e

    for s in cols:
        d = info[s]
        for p in range(d['i'], d['j'] + 1):
            rhs = 0
            if p <= d['k']:
                rhs = rhs + m.d[s, p]
            if p >= d['l']:
                rhs = rhs + m.b[s, p]
            m.c.add(feed(s, p) == rhs)
        m.c.add(m.VT[s] - m.VB[s] == vfeed(s))
        has_int = d['l'] <= d['k'] - 1
        tgt = m.VT[s]
        if permissive and has_int:
            ups = pyo.Var(bounds=(0, 1e4))
            m.add_component(f"ups_{s[0]}_{s[1]}", ups)
            m.c.add(m.VT[s] >= ups)
            tgt = ups
        for r in range(d['l'] - 1, d['k'] + 1):
            m.c.add(sum(al[p] * feed(s, p) / (al[p] - m.th[s, r]) for p in range(d['i'], d['j'] + 1))
                    == vfeed(s))
            g = sum(al[p] * m.d[s, p] / (al[p] - m.th[s, r]) for p in range(d['i'], d['k'] + 1))
            if d['l'] <= r <= d['k'] - 1:
                m.c.add(tgt == g)
            else:
                m.c.add(tgt >= g)
    for (t, u) in side:
        m.c.add(m.VT[t] == m.VB[u])
    m.V = pyo.Expression(expr=sum(m.VB[u] for u in reb))
    m.Pj = pyo.Expression(expr=sum(al[p] * Fd[p] / (al[p] - m.Th[jpk]) for p in range(1, jpk + 1)) - (1 - q))
    m.c.add(m.Pj >= float(os.environ.get('PMIN', '0.05')))
    m.obj = pyo.Objective(expr=m.V / m.Pj, sense=pyo.minimize)
    return m


def sample_params(N, rng, qlo, qhi, rmin=1.02, rmax=30.0, zmin=1e-3):
    rat = [math.exp(rng.uniform(math.log(rmin), math.log(rmax))) for _ in range(N - 1)]
    al = [1.0]
    for r in reversed(rat):
        al.insert(0, al[0] * r)
    w = [rng.expovariate(1.0) ** rng.choice([1, 2, 3]) for _ in range(N)]
    s = sum(w)
    z = [max(x / s, 2 * zmin) for x in w]
    s = sum(z)
    z = [x / s for x in z]
    q = rng.uniform(qlo, qhi)
    return al, z, q


import os
PERMISSIVE = os.environ.get('PERMISSIVE') == '1'


def run(N, qlo, qhi, starts, seed, outfile, ids=None):
    cfgs = all_configurations(N)
    if ids:
        cfgs = [c for c in cfgs if c['id'] in ids]
    rng = random.Random(seed)
    out = []
    t0 = time.time()
    for c in cfgs:
        bestc = None
        for jpk in range(1, N):
            mj = build_joint(c, jpk, qlo, qhi, permissive=PERMISSIVE)
            for st in range(starts):
                al, z, q = sample_params(N, rng, qlo, qhi)
                # feasible start from fixed-parameter model
                mf = build_model(c, al, z, q, permissive=PERMISSIVE)
                random_start(mf, rng)
                try:
                    r0 = solver().solve(mf, load_solutions=False)
                    if r0.solver.termination_condition == TerminationCondition.optimal:
                        mf.solutions.load_from(r0)
                except Exception:  # noqa: BLE001
                    pass
                for p in range(1, N + 1):
                    mj.al[p].value = al[p - 1]
                    mj.z[p].value = z[p - 1]
                mj.q.value = q
                from hsmodel import feed_roots
                th = feed_roots(al, z, q, dps=20)
                for r in range(1, N):
                    mj.Th[r].value = float(th[r - 1])
                for k, v in mf.d.items():
                    mj.d[k].value = v.value
                for k, v in mf.b.items():
                    mj.b[k].value = v.value
                for k, v in mf.VT.items():
                    mj.VT[k].value = v.value
                for k, v in mf.VB.items():
                    mj.VB[k].value = v.value
                for k, v in mf.th.items():
                    mj.th[k].value = v.value
                for comp in mf.component_objects(pyo.Var):
                    if comp.name.startswith('ups_'):
                        getattr(mj, comp.name).value = comp.value
                try:
                    res = solver().solve(mj, load_solutions=False)
                except Exception:  # noqa: BLE001
                    continue
                if res.solver.termination_condition != TerminationCondition.optimal:
                    continue
                mj.solutions.load_from(res)
                viol = max_violation(mj)
                if viol > 1e-7:
                    continue
                val = pyo.value(mj.obj)
                pars = ([pyo.value(mj.al[p]) for p in range(1, N + 1)],
                        [pyo.value(mj.z[p]) for p in range(1, N + 1)], pyo.value(mj.q))
                if bestc is None or val < bestc['ratio']:
                    bestc = {'id': c['id'], 'name': config_name(c), 'j': jpk, 'ratio': val,
                             'alpha': pars[0], 'z': pars[1], 'q': pars[2], 'viol': viol}
        if bestc is not None:
            vt, vp = petlyuk_value(bestc['alpha'], bestc['z'], bestc['q'])
            bestc['petlyuk_reboiler'] = float(vp)
        out.append(bestc)
        print(c['id'], None if bestc is None else (round(bestc['ratio'], 9), bestc['j']),
              round(time.time() - t0, 1), flush=True)
        with open(outfile, 'w') as fh:
            json.dump(out, fh, indent=1)


if __name__ == '__main__':
    N = int(sys.argv[1])
    qlo, qhi = float(sys.argv[2]), float(sys.argv[3])
    starts, seed = int(sys.argv[4]), int(sys.argv[5])
    outfile = sys.argv[6]
    ids = [int(x) for x in sys.argv[7:]] or None
    run(N, qlo, qhi, starts, seed, outfile, ids)
