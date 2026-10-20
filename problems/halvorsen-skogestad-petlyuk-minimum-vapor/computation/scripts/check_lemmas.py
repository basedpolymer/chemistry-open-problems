"""Numerical checks requested by Agent T, at optimal and random feasible points (strict and permissive models):
  (L1) column sign rule: for every mixed column, not (A_c < 0 and B_c < 0) at theta_j  [= Agent T's Lemma 1];
  (RO) recovery ordering: within each column, top recoveries d_p/f_p are non-increasing in p over the
       distributing components (p = l..k);
  (S5) strict model only: j interior => not (A_c > 0 and B_c > 0);
  (BP) sign graph has no light->heavy path (sign_paths.sign_graph_bad).
Environment: WSL Fedora, Python 3, pyomo + ipopt.  Usage: python check_lemmas.py N ntrials seed permissive(0/1)
Quick check (about 15 s): python check_lemmas.py 4 20 7 0   -> prints PASS/FAIL, exit code 1 on FAIL.
Runtime: about 10-20 min for 60 trials (N = 4), single thread.
"""
import random
import sys

import pyomo.environ as pyo
from pyomo.opt import TerminationCondition

from hsmodel import all_configurations, config_name, feed_roots, petlyuk_peaks
from nlp import build_model, random_start, solver, structure, max_violation, snapshot, solve_config
from cut_identity import psi_terms
from sign_paths import sign_graph_bad


def random_feasible(cfg, alpha, F, q, rng, permissive):
    m = build_model(cfg, alpha, F, q, permissive=permissive)
    cols = list(m.VT.keys())
    m.obj.deactivate()
    w = {c: rng.uniform(-0.2, 1.0) for c in cols}
    w2 = {c: rng.random() for c in cols}
    wd = {k: rng.uniform(-1, 1) for k in m.d.keys()}
    m.obj2 = pyo.Objective(expr=sum(w[c] * m.VB[c] + w2[c] * m.VT[c] for c in cols)
                           + sum(wd[k] * m.d[k] for k in m.d.keys()))
    random_start(m, rng, 2.0)
    res = solver().solve(m, load_solutions=False)
    if res.solver.termination_condition != TerminationCondition.optimal:
        return None
    m.solutions.load_from(res)
    if max_violation(m) > 1e-7:
        return None
    return snapshot(m)


if __name__ == '__main__':
    N, ntr, seed, perm = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), bool(int(sys.argv[4]))
    rng = random.Random(seed)
    C = all_configurations(N)
    stats = {'L1': 0, 'RO': 0, 'S5': 0, 'BP': 0, 'tests': 0, 'points': 0}
    worst_ro = 1.0
    for tr in range(ntr):
        c = rng.choice(C)
        rat = [rng.uniform(1.05, 6) for _ in range(N - 1)]
        al = [1.0]
        for r in reversed(rat):
            al.insert(0, al[0] * r)
        wz = [rng.random() + 0.01 for _ in range(N)]
        z = [x / sum(wz) for x in wz]
        q = rng.choice([1.0, rng.random(), rng.uniform(-0.5, 1.5)])
        ths = [float(x) for x in feed_roots(al, z, q)]
        sols = []
        best, _ = solve_config(c, al, z, q, starts=2, seed=tr, permissive=perm)
        if best:
            sols.append(best[1])
        for _ in range(3):
            s = random_feasible(c, al, z, q, rng, perm)
            if s:
                sols.append(s)
        cols, info, reb, side = structure(c)
        for sol in sols:
            stats['points'] += 1
            d = {eval(k): v for k, v in sol['d'].items()}
            b = {eval(k): v for k, v in sol['b'].items()}
            for s in cols:
                I = info[s]
                rec = []
                for p in range(I['l'], I['k'] + 1):
                    f = d.get((s[0], s[1], p), 0) + b.get((s[0], s[1], p), 0)
                    if f > 1e-9:
                        rec.append(d.get((s[0], s[1], p), 0) / f)
                for a1, a2 in zip(rec, rec[1:]):
                    worst_ro = min(worst_ro, a1 - a2)
                    if a1 < a2 - 1e-7:
                        stats['RO'] += 1
                        print('RO-VIOLATION', c['id'], config_name(c), s, rec, al, z, q, flush=True)
            for j in range(1, N):
                _, _, psiT, psiB, _ = psi_terms(c, sol, al, z, q, ths[j - 1])
                stats['tests'] += 1
                for s in cols:
                    I = info[s]
                    if not (I['i'] <= j < I['j']):
                        continue
                    if psiT[s] < -1e-8 and psiB[s] < -1e-8:
                        stats['L1'] += 1
                        print('L1-VIOLATION', c['id'], s, j, psiT[s], psiB[s], flush=True)
                    if (not perm) and I['l'] <= j <= I['k'] - 1 and psiT[s] > 1e-8 and psiB[s] > 1e-8:
                        stats['S5'] += 1
                        print('S5-VIOLATION', c['id'], s, j, psiT[s], psiB[s], flush=True)
                bad, signs = sign_graph_bad(c, sol, al, z, q, j, ths[j - 1], tol=1e-8)
                if bad:
                    stats['BP'] += 1
                    print('BAD-PATH', 'perm' if perm else 'strict', c['id'], config_name(c), j, al, z, q, signs,
                          flush=True)
    print(stats, 'min consecutive recovery difference', worst_ro)
    # PASS criterion: no violation of the sign rule (L1), of the recovery ordering (RO) and, in the strict model,
    # of (S5) and of the sign-graph property (BP).  In the permissive model bad paths are allowed (informational).
    bad = stats['L1'] + stats['RO'] + (0 if perm else stats['S5'] + stats['BP'])
    print(('PASS' if bad == 0 else 'FAIL') + f": lemma checks, N={N}, {'permissive' if perm else 'strict'} model, "
          f"{stats['points']} points, {stats['tests']} (point, j) tests")
    sys.exit(0 if bad == 0 else 1)
