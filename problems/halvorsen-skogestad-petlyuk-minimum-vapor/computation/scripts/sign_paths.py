"""Sign-graph test of a sufficient condition for R >= P_j - (1-q) (Agent C).

At theta = theta_j, every section e (lower node x, upper node y) gets the sign of Psi_e(theta_j).  Add an arc
x -> y if Psi_e > 0 and y -> x if Psi_e < 0.  If no directed path leads from a light-only node (stream within
[1,j]) to a heavy-only node (stream within [j+1,N]), then U = {nodes reachable from the light nodes} is a
partition whose crossing sections all enter the cut identity with a nonnegative sign, which proves
R - (P_j - (1-q)) >= 0 at that operating point.  This script tests the condition at
  (a) optimal operating points, (b) feasible points obtained by minimising random positive objectives.

Environment: WSL Fedora, Python 3, pyomo + ipopt.  Floating point exploration.
Runtime: about 10-20 min for 60 trials (N = 4), single thread.
Usage: python sign_paths.py N ntrials seed
"""
import random
import sys

import pyomo.environ as pyo
from pyomo.opt import TerminationCondition

from hsmodel import all_configurations, config_name, feed_roots, petlyuk_peaks
from nlp import build_model, random_start, solver, structure, max_violation, snapshot
from cut_identity import psi_terms


def sign_graph_bad(cfg, sol, alpha, F, q, j, th, tol=1e-9):
    N = cfg['N']
    cols, info, psiT, psiB, hx = psi_terms(cfg, sol, alpha, F, q, th)
    # nodes: streams (mixtures and pure products); section edges
    edges = []  # (lower, upper, psi)
    for c in cols:
        I = info[c]
        dist = (I['i'], I['k'])
        res = (I['l'], I['j'])
        edges.append((c, dist, psiT[c]))   # top section: from feed node c up to distillate node
        edges.append((res, c, psiB[c]))    # bottom section: from residue node up to feed node c
    adj = {}
    for (x, y, v) in edges:
        if v > tol:
            adj.setdefault(x, []).append(y)
        elif v < -tol:
            adj.setdefault(y, []).append(x)
    light = lambda s: s[1] <= j
    heavy = lambda s: s[0] >= j + 1
    nodes = set([e[0] for e in edges] + [e[1] for e in edges])
    stack = [s for s in nodes if light(s)]
    seen = set(stack)
    while stack:
        x = stack.pop()
        for y in adj.get(x, []):
            if y not in seen:
                seen.add(y)
                stack.append(y)
    bad = [s for s in seen if heavy(s)]
    return bad, {str(c): (round(psiT[c], 6), round(psiB[c], 6)) for c in cols}


def random_feasible(cfg, alpha, F, q, rng):
    m = build_model(cfg, alpha, F, q)
    cols = list(m.VT.keys())
    w = {c: rng.random() for c in cols}
    w2 = {c: rng.random() for c in cols}
    m.obj.deactivate()
    m.obj2 = pyo.Objective(expr=sum(w[c] * m.VB[c] + w2[c] * m.VT[c] for c in cols))
    random_start(m, rng, 2.0)
    res = solver().solve(m, load_solutions=False)
    if res.solver.termination_condition != TerminationCondition.optimal:
        return None
    m.solutions.load_from(res)
    if max_violation(m) > 1e-7:
        return None
    return pyo.value(m.obj), snapshot(m)


if __name__ == '__main__':
    N = int(sys.argv[1])
    ntr = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]))
    C = all_configurations(N)
    nbad = 0
    ntot = 0
    for tr in range(ntr):
        c = rng.choice(C)
        rat = [rng.uniform(1.05, 5) for _ in range(N - 1)]
        al = [1.0]
        for r in reversed(rat):
            al.insert(0, al[0] * r)
        wz = [rng.random() + 0.02 for _ in range(N)]
        z = [x / sum(wz) for x in wz]
        q = rng.choice([1.0, rng.random()])
        ths = [float(x) for x in feed_roots(al, z, q)]
        P = [float(x) for x in petlyuk_peaks(al, z, q)]
        for mode in ('opt', 'rand'):
            if mode == 'opt':
                from nlp import solve_config
                best, _ = solve_config(c, al, z, q, starts=3, seed=tr)
            else:
                best = random_feasible(c, al, z, q, rng)
            if best is None:
                continue
            R, sol = best
            if mode == 'rand':
                R = sum(sol['VB'][str(u)] for u in structure(c)[2])
            for j in range(1, N):
                bad, signs = sign_graph_bad(c, sol, al, z, q, j, ths[j - 1])
                ntot += 1
                if bad:
                    nbad += 1
                    print('BAD-PATH', mode, c['id'], config_name(c), 'j', j, 'alpha', [round(a, 4) for a in al],
                          'z', [round(x, 4) for x in z], 'q', round(q, 4), 'gap', R - (P[j - 1] - (1 - q)),
                          'reach heavy', bad, signs, flush=True)
    print('total', ntot, 'bad', nbad)
