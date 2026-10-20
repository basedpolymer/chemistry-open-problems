"""Parallel driver for the PERMISSIVE-model bad-path study (see badpath_search.py), N = 5.

Targets: every TC/HX variant of the basic configurations that admit a sign-consistent light->heavy path in the
permissive model AND whose mixed sub-network for that j contains a merge node (a mixed submixture produced by
two splits = mixed side draw) — the case left open by Agent G's Theorem G3.  For each (variant, j, path):
maximise the margin t of the path's sign requirements over parameters (alpha, z; q fixed) and operating points;
if t* > 1e-3, minimise R / (P_j - (1-q)) inside the bad-path region (margin >= 1e-5).  A ratio < 1 would be a
counterexample to the conjecture in the permissive model.

Environment: WSL Fedora, Python 3, pyomo + ipopt, multiprocessing (W workers).  Floating point exploration.
Runtime: about 2-3 h with 4 workers (1248 tasks, 4 starts each).
Usage: python badpath_parallel.py outfile.jsonl starts seed q W [maxvariants_per_basic]
"""
import json
import multiprocessing as mpc
import random
import sys
import time

import badpath_search as bp
from hsmodel import all_configurations

bp.STRICT = False
_C = None


def _init():
    global _C
    bp.STRICT = False
    _C = all_configurations(5)


def has_merge(c, j):
    for s, prods in c['producers'].items():
        if s[0] == s[1] or s == (1, 5):
            continue
        if s[0] <= j < s[1] and len(prods) == 2:
            return True
    return False


def task(args):
    cid, j, reqs, starts, seed, q = args
    import pyomo.environ as pyo
    from pyomo.opt import TerminationCondition
    from adversarial import sample_params
    from hsmodel import feed_roots
    from nlp import build_model, random_start, solver, max_violation
    c = _C[cid]
    rng = random.Random(seed)
    m = bp.margin_model(c, j, reqs, q, q)
    m.t.setub(0.02)  # small margins keep vapors moderate, so that phase 2 starts from a reasonable point
    best_t, best_gap, minratio, best_pt = -1e9, None, float('inf'), None
    for st in range(starts):
        al, z, qq = sample_params(5, rng, q, q)
        mf = build_model(c, al, z, qq)
        random_start(mf, rng, 2.0)
        for p in range(1, 6):
            m.al[p].value = al[p - 1]
            m.z[p].value = z[p - 1]
        m.q.value = qq
        ths = feed_roots(al, z, qq, dps=20)
        for r in range(1, 5):
            m.Th[r].value = float(ths[r - 1])
        for comp in ('d', 'b', 'VT', 'VB', 'th'):
            for k, v in getattr(mf, comp).items():
                getattr(m, comp)[k].value = v.value
        m.t.value = -1.0
        try:
            res = solver().solve(m, load_solutions=False)
        except Exception:  # noqa: BLE001
            continue
        if res.solver.termination_condition != TerminationCondition.optimal:
            continue
        m.solutions.load_from(res)
        if max_violation(m) > 1e-7:
            continue
        tv = pyo.value(m.t)
        if tv > best_t:
            best_t, best_gap = tv, pyo.value(m.Rgap)
        if tv > 1e-3:
            m.t.setlb(1e-5)
            m.mobj.deactivate()
            m.obj.activate()
            try:
                r2 = solver().solve(m, load_solutions=False)
                if r2.solver.termination_condition == TerminationCondition.optimal:
                    m.solutions.load_from(r2)
                    if max_violation(m) < 1e-7 and pyo.value(m.obj) < minratio:
                        minratio = pyo.value(m.obj)
                        best_pt = {'alpha': [pyo.value(m.al[p]) for p in range(1, 6)],
                                   'z': [pyo.value(m.z[p]) for p in range(1, 6)], 'q': pyo.value(m.q)}
            except Exception:  # noqa: BLE001
                pass
            m.obj.deactivate()
            m.mobj.activate()
            m.t.setlb(-10)
    return {'id': cid, 'j': j, 'reqs': [[list(x[0]), x[1], x[2]] for x in reqs], 't': best_t, 'gap': best_gap,
            'min_ratio_in_badpath_region': None if minratio == float('inf') else minratio, 'point': best_pt}


if __name__ == '__main__':
    out, starts, seed, q, W = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    maxv = int(sys.argv[6]) if len(sys.argv) > 6 else 10 ** 9
    _init()
    rng = random.Random(seed)
    tasks = []
    per_basic = {}
    for c in _C:
        for j in (2, 3):
            if not has_merge(c, j):
                continue
            P = bp.potential_paths(c, j)
            if not P:
                continue
            key = (c['basic_id'], j)
            per_basic.setdefault(key, 0)
            if per_basic[key] >= maxv:
                continue
            per_basic[key] += 1
            for reqs in P:
                tasks.append((c['id'], j, reqs, starts, rng.randrange(10 ** 9), q))
    print('tasks', len(tasks), 'basic (b,j) pairs', len(per_basic), flush=True)
    t0 = time.time()
    worst = float('inf')
    with mpc.Pool(W, initializer=_init) as pool, open(out, 'a') as fh:
        for n, r in enumerate(pool.imap_unordered(task, tasks)):
            fh.write(json.dumps(r) + '\n')
            fh.flush()
            mr = r['min_ratio_in_badpath_region']
            if mr is not None:
                worst = min(worst, mr)
            if (mr is not None and mr < 1 - 1e-7) or n % 25 == 0:
                print(n, r['id'], r['j'], f"t*={r['t']:.2e}", mr, 'worst so far', worst,
                      'CANDIDATE' if (mr is not None and mr < 1 - 1e-7) else '', round(time.time() - t0),
                      flush=True)
    print('done; worst min ratio in bad-path regions', worst)
