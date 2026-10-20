"""Sign-relaxation check (Agent T).  For each configuration and each split j, decide whether the column sign
rule (Lemma 1) alone implies R >= H_j, i.e. whether  sum_{light exits} t + sum kappa >= 0  holds for every
choice of signs/values allowed by:  a_c = t_c + b_c,  (a_c >= 0 and t_c >= 0) or (a_c <= 0 and b_c <= 0),
a_root = 0, child inflow = (t_parent - kappa) + (b_parent' + rho), kappa, rho >= 0 (HX links only).
Also imposed: the unconditional part of the rule (k < j => t >= 0, l > j+1 => b <= 0).
The feasible set is a finite union of polyhedral cones, so the bound fails iff some cone contains a point with
objective <= -1 (LP feasibility per sign pattern).  Exact enough: LP with HiGHS; patterns are few.
Environment: WSL Fedora, Python 3, scipy (HiGHS), read-only import of computation/scripts/hsmodel.py.
Runtime: 'certify4' about 20 s; full listing for N = 5 about 6 minutes.
Usage: python signcheck.py certify4 | python signcheck.py N"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '../../computation/scripts')
from hsmodel import all_configurations, basic_configurations, config_name, candidates

def mixed_structure(cfg, j):
    N = cfg['N']
    M = [s for s in cfg['streams'] if s[0] <= j < s[1]]
    return M

def sign_hard(cfg, j, hx):
    """hx: dict stream->'TC'/'HX' for single-produced submixtures. Returns (hard, witness or None)."""
    N = cfg['N']
    M = mixed_structure(cfg, j)
    idx = {}
    nv = 0
    for s in M:
        idx[('t', s)] = nv; nv += 1
        idx[('b', s)] = nv; nv += 1
    # HX variables on links INTO mixed streams
    hxvars = []
    for s in M:
        if s == (1, N):
            continue
        prods = cfg['producers'][s]
        if len(prods) == 1 and hx.get(s) == 'HX':
            kind, par = prods[0]
            idx[('h', s)] = nv; nv += 1
            hxvars.append((s, kind))
    def inflow(s):
        """linear expression (dict var->coef) of a_s"""
        e = {}
        if s == (1, N):
            return e
        for kind, par in cfg['producers'][s]:
            if kind == 'top':
                e[idx[('t', par)]] = e.get(idx[('t', par)], 0) + 1
            else:
                e[idx[('b', par)]] = e.get(idx[('b', par)], 0) + 1
        if ('h', s) in idx:
            kind = cfg['producers'][s][0][0]
            e[idx[('h', s)]] = -1 if kind == 'top' else +1
        return e
    light_exit = [s for s in M if cfg['splits'][s][0] <= j]
    obj = np.zeros(nv)
    for s in light_exit:
        obj[idx[('t', s)]] += 1
    for (s, kind) in hxvars:
        if kind == 'top':
            obj[idx[('h', s)]] += 1
    nonroot = [s for s in M if s != (1, N)]
    for pattern in itertools.product([1, -1], repeat=len(nonroot)):
        A_ub, b_ub, A_eq, b_eq = [], [], [], []
        def row(e):
            r = np.zeros(nv)
            for k, v in e.items():
                r[k] += v
            return r
        # objective <= -1
        A_ub.append(obj.copy()); b_ub.append(-1.0)
        bounds = [(None, None)] * nv
        for (s, kind) in hxvars:
            bounds[idx[('h', s)]] = (0, None)
        # root: t>=0, b<=0, t+b=0
        r0 = (1, N)
        bounds[idx[('t', r0)]] = (0, None)
        bounds[idx[('b', r0)]] = (None, 0)
        A_eq.append(row({idx[('t', r0)]: 1, idx[('b', r0)]: 1})); b_eq.append(0.0)
        for s in M:   # unconditional part of the sign rule: k < j => A >= 0 (t >= 0); l > j+1 => B >= 0 (b <= 0)
            k_s, l_s = cfg['splits'][s]
            if k_s < j:
                lo, hi = bounds[idx[('t', s)]]
                bounds[idx[('t', s)]] = (0, hi)
            if l_s > j + 1:
                lo, hi = bounds[idx[('b', s)]]
                bounds[idx[('b', s)]] = (lo, 0)
        for s, sg in zip(nonroot, pattern):
            a = inflow(s)
            # a = t + b
            e = dict(a)
            e[idx[('t', s)]] = e.get(idx[('t', s)], 0) - 1
            e[idx[('b', s)]] = e.get(idx[('b', s)], 0) - 1
            A_eq.append(row(e)); b_eq.append(0.0)
            if sg == 1:   # a >= 0, t >= 0
                A_ub.append(-row(a)); b_ub.append(0.0)
                bounds[idx[('t', s)]] = (0, bounds[idx[('t', s)]][1])
            else:         # a <= 0, b <= 0
                A_ub.append(row(a)); b_ub.append(0.0)
                bounds[idx[('b', s)]] = (bounds[idx[('b', s)]][0], 0)
        res = linprog(np.zeros(nv), A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      A_eq=np.array(A_eq) if A_eq else None, b_eq=np.array(b_eq) if b_eq else None,
                      bounds=bounds, method='highs')
        if res.status == 0:
            return True, (dict(zip(nonroot, pattern)), res.x)
    return False, None

def certify4():
    """Cross-check of the case analysis of Theorem 2 (N = 4): the (configuration, j) pairs NOT settled by the
    column sign rule alone are exactly the 16 TC/HX variants of family K at j = 2."""
    A = all_configurations(4)
    K = {(1, 4), (1, 3), (2, 4), (1, 2), (3, 4)}
    hard = [(c['id'], j) for c in A for j in (1, 2, 3) if sign_hard(c, j, c['hx'])[0]]
    expect = [(c['id'], 2) for c in A if set(c['streams']) == K]
    ok = sorted(hard) == sorted(expect) and len(expect) == 16
    print('N = 4: sign-rule-hard pairs', len(hard), '; family K variants at j = 2:', len(expect),
          '->', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    if sys.argv[1] == 'certify4':
        sys.exit(0 if certify4() else 1)
    N = int(sys.argv[1])
    A = all_configurations(N)
    hard = {}
    for c in A:
        for j in range(1, N):
            h, w = sign_hard(c, j, c['hx'])
            if h:
                hard.setdefault(j, []).append(c['id'])
    print('N =', N, 'configurations:', len(A))
    for j in range(1, N):
        L = hard.get(j, [])
        print(f' j={j}: sign-hard configurations: {len(L)}')
        basics = sorted(set(A[i]['basic_id'] for i in L))
        for b in basics:
            ex = [i for i in L if A[i]['basic_id'] == b]
            print(f'   basic {b}: {len(ex)} variants, e.g. {config_name(A[ex[0]])}')

def dump(N, path):
    """Write {j: [ids]} of sign-hard (configuration, j) pairs and the union, as JSON."""
    import json
    A = all_configurations(N)
    hard = {j: [] for j in range(1, N)}
    for c in A:
        for j in range(1, N):
            if sign_hard(c, j, c['hx'])[0]:
                hard[j].append(c['id'])
    union = sorted(set().union(*[set(v) for v in hard.values()]))
    json.dump({'N': N, 'n_configs': len(A), 'hard': hard, 'union': union,
               'basic_of_union': sorted(set(A[i]['basic_id'] for i in union))}, open(path, 'w'), indent=0)
    print(N, 'configs', len(A), 'hard per j', {j: len(v) for j, v in hard.items()}, 'union', len(union))
