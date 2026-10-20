"""Adversarial search for 'bad sign paths' (see sign_paths.py) over parameters AND operating points.

For a configuration c and split index j, enumerate the simple paths, in the potential sign graph, from a light
node to a heavy node.  Each path imposes sign requirements on A_x = Psi_top(x)(theta_j), B_x = Psi_bot(x)(theta_j).
For each path we maximise the margin t subject to all model constraints (parameters free in a box) and
sign * Psi >= t.  t* > 0 means a feasible operating point exists where the sign-graph argument fails at j;
we then also report R - (P_j - (1-q)) at that point (negative = counterexample to the conjecture).

Environment: WSL Fedora, Python 3, pyomo + ipopt.  Floating point exploration.
Runtime: about 30 s per (configuration, path) with 10 starts, single thread.
Usage: python badpath_search.py N starts seed outfile [qlo qhi] [config ids...]
"""
import json
import random
import sys
import time

import pyomo.environ as pyo
from pyomo.opt import TerminationCondition

from adversarial import build_joint, sample_params
from hsmodel import all_configurations, config_name, feed_roots
from nlp import build_model, random_start, solver, structure, max_violation


STRICT = True  # False: TG literal (permissive) model, where (S5) and the root-interior rule are not available


def consistent(info, N, j, reqs):
    """Necessary sign conditions at theta_j valid at every feasible point (lemmas S3-S5 in PROGRESS.md):
    (S3) k < j => A_c >= 0 ;  l > j+1 => B_c >= 0 ;
    (S4) never A_c < 0 and B_c < 0 ;
    (S5) j interior for c (l <= j <= k-1) => not (A_c > 0 and B_c > 0) ;
    (root) A = B at the main column (its feed function vanishes at every feed root), hence A = B >= 0, and
           A = B = 0 if j is interior for the main column."""
    sg = {}
    for (c, w, s) in reqs:
        if sg.get((c, w), s) != s:
            return False
        sg[(c, w)] = s
    for (c, w), s in sg.items():
        I = info[c]
        if w == 'A' and s < 0 and I['k'] < j:
            return False
        if w == 'B' and s < 0 and I['l'] > j + 1:
            return False
        if sg.get((c, 'A')) == -1 and sg.get((c, 'B')) == -1:
            return False
        if STRICT and I['l'] <= j <= I['k'] - 1 and sg.get((c, 'A')) == 1 and sg.get((c, 'B')) == 1:
            return False
        if c == (1, N):
            if s < 0:
                return False
            if STRICT and I['l'] <= j <= I['k'] - 1:
                return False
    return True


def potential_paths(cfg, j, maxlen=12):
    """Simple paths light -> heavy with sign requirements [(col, 'A'|'B', +1|-1), ...]."""
    N = cfg['N']
    cols, info, reb, side = structure(cfg)
    light = lambda s: s[1] <= j
    heavy = lambda s: s[0] >= j + 1
    moves = {}  # node -> list of (next node, requirement)

    def add(x, y, req):
        moves.setdefault(x, []).append((y, req))
    for c in cols:
        I = info[c]
        dist, res = (I['i'], I['k']), (I['l'], I['j'])
        # along the top section of c (between node c and node dist)
        add(c, dist, (c, 'A', +1))   # upward, A_c > 0
        add(dist, c, (c, 'A', -1))   # downward, A_c < 0
        # along the bottom section of c (between node res and node c)
        add(res, c, (c, 'B', +1))    # upward, B_c > 0
        add(c, res, (c, 'B', -1))    # downward, B_c < 0
    paths = []
    starts = set(x for x in moves if light(x))

    def dfs(node, visited, reqs):
        if len(reqs) > maxlen:
            return
        for (y, req) in moves.get(node, []):
            if y in visited:
                continue
            # consistency: same quantity with opposite signs, or A<0 and B<0 on same column
            bad = False
            for (cc, w, sg) in reqs:
                if cc == req[0] and w == req[1] and sg != req[2]:
                    bad = True
                if cc == req[0] and w != req[1] and sg == -1 and req[2] == -1:
                    bad = True
            if bad or not consistent(info, N, j, reqs + [req]):
                continue
            nr = reqs + [req]
            if heavy(y):
                paths.append(nr)
                continue
            if light(y):
                continue  # restarting from a light node is covered by other starts
            dfs(y, visited | {y}, nr)
    for s in starts:
        dfs(s, {s}, [])
    # deduplicate requirement sets
    uniq = {}
    for p in paths:
        key = tuple(sorted(set(p)))
        uniq[key] = p
    return list(uniq.keys())


def margin_model(cfg, j, reqs, qlo, qhi):
    m = build_joint(cfg, j, qlo, qhi, permissive=not STRICT)
    m.obj.deactivate()
    cols, info, reb, side = structure(cfg)
    th = m.Th[j]
    m.t = pyo.Var(bounds=(-10, 1.0))
    al = m.al

    def A(c):
        I = info[c]
        return m.VT[c] - sum(al[p] * m.d[c, p] / (al[p] - th) for p in range(I['i'], I['k'] + 1))

    def B(c):
        I = info[c]
        return m.VB[c] + sum(al[p] * m.b[c, p] / (al[p] - th) for p in range(I['l'], I['j'] + 1))
    m.req = pyo.ConstraintList()
    for (c, w, sg) in reqs:
        e = A(c) if w == 'A' else B(c)
        m.req.add(sg * e >= m.t)
    m.Rgap = pyo.Expression(expr=m.V - m.Pj)
    m.mobj = pyo.Objective(expr=-m.t)
    return m


def run(N, starts, seed, outfile, qlo, qhi, ids):
    C = all_configurations(N)
    if ids:
        C = [c for c in C if c['id'] in ids]
    rng = random.Random(seed)
    out = []
    t0 = time.time()
    for c in C:
        for j in range(1, N):
            P = potential_paths(c, j)
            for reqs in P:
                best = None
                minratio = float('inf')
                m = margin_model(c, j, reqs, qlo, qhi)
                for st in range(starts):
                    al, z, q = sample_params(N, rng, qlo, qhi)
                    mf = build_model(c, al, z, q)
                    random_start(mf, rng, 2.0)
                    for p in range(1, N + 1):
                        m.al[p].value = al[p - 1]
                        m.z[p].value = z[p - 1]
                    m.q.value = q
                    ths = feed_roots(al, z, q, dps=20)
                    for r in range(1, N):
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
                    ratio = None
                    if tv > 1e-3:
                        # phase 2: minimise R/(P_j-(1-q)) inside the bad-path region (margin >= 1e-4)
                        m.t.setlb(1e-4)
                        m.mobj.deactivate()
                        m.obj.activate()
                        try:
                            r2 = solver().solve(m, load_solutions=False)
                            if r2.solver.termination_condition == TerminationCondition.optimal:
                                m.solutions.load_from(r2)
                                if max_violation(m) < 1e-7:
                                    ratio = pyo.value(m.obj)
                        except Exception:  # noqa: BLE001
                            pass
                        m.obj.deactivate()
                        m.mobj.activate()
                        m.t.setlb(-10)
                    if ratio is not None:
                        minratio = min(minratio, ratio)
                    if best is None or tv > best['t']:
                        best = {'id': c['id'], 'j': j, 'reqs': [[str(x[0]), x[1], x[2]] for x in reqs], 't': tv,
                                'gap': pyo.value(m.Rgap),
                                'alpha': [pyo.value(m.al[p]) for p in range(1, N + 1)],
                                'z': [pyo.value(m.z[p]) for p in range(1, N + 1)], 'q': pyo.value(m.q)}
                if best is not None:
                    best['min_ratio_in_badpath_region'] = minratio if minratio < float('inf') else None
                    out.append(best)
                    flag = 'BADPATH' if best['t'] > 1e-6 else ''
                    if best['t'] > 1e-6:
                        flag += ' COUNTEREXAMPLE?' if best['gap'] < -1e-7 else ''
                    print(c['id'], j, len(reqs), f"t*={best['t']:.3e}", f"gap={best['gap']:.4e}",
                          f"minratio={minratio:.9f}", flag,
                          round(time.time() - t0), flush=True)
        with open(outfile, 'w') as fh:
            json.dump(out, fh, indent=1)


if __name__ == '__main__':
    N = int(sys.argv[1])
    starts, seed = int(sys.argv[2]), int(sys.argv[3])
    outfile = sys.argv[4]
    qlo = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    qhi = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    ids = [int(x) for x in sys.argv[7:]]
    import os
    if os.environ.get('PERMISSIVE') == '1':
        STRICT = False
        globals()['STRICT'] = False
    run(N, starts, seed, outfile, qlo, qhi, ids)
