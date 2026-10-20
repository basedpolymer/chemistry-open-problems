"""Numerical check of the 'cut identity' (Agent C, 2026-10-20) on optimal operating points.

For a feasible operating point of a configuration, a feed root theta_j, and Psi_s(theta) = V_s - sum_p alpha_p
w_{s,p}/(alpha_p - theta) for every section s (w = net upward component flows: d >= 0 in top sections, -b in
bottom sections), the canonical top cut U_j = {nodes whose stream is contained in [1, j]} gives

   R - (P_j - (1-q)) = FR(U_j) + FC(complement of U_j) + sum_{s in E_j} Psi_s(theta_j)          (*)

where R = total reboiler vapor, E_j = top sections of the columns [i,m] (i <= j < m) whose distillate [i,k]
has k <= j.  This script solves the minimum-vapor NLP for random configurations/parameters and checks (*).

Environment: WSL Fedora, Python 3, pyomo, ipopt.  Floating point check only.
Runtime: about 1 min for 25 trials, single thread.
"""
import random
import sys

from hsmodel import all_configurations, config_name, feed_roots, petlyuk_peaks
from nlp import solve_config, structure


def psi_terms(cfg, sol, alpha, F, q, th):
    N = cfg['N']
    cols, info, reb, side = structure(cfg)
    al = {p: alpha[p - 1] for p in range(1, N + 1)}
    d = {eval(k): v for k, v in sol['d'].items()}
    b = {eval(k): v for k, v in sol['b'].items()}
    VT = {eval(k): v for k, v in sol['VT'].items()}
    VB = {eval(k): v for k, v in sol['VB'].items()}
    dd = lambda s, p: d.get((s[0], s[1], p), 0.0)
    bb = lambda s, p: b.get((s[0], s[1], p), 0.0)
    psiT = {s: VT[s] - sum(al[p] * dd(s, p) / (al[p] - th) for p in range(info[s]['i'], info[s]['k'] + 1))
            for s in cols}
    psiB = {s: VB[s] + sum(al[p] * bb(s, p) / (al[p] - th) for p in range(info[s]['l'], info[s]['j'] + 1))
            for s in cols}
    # heat exchangers: list of (node stream, 'FC'|'FR', duty)
    hx = []
    for s in cols:
        if s == (1, N):
            continue
        dct = info[s]
        if dct['kind'] == 'HX':
            if dct['top'] is not None:
                t = dct['top']
                D = sum(dd(t, p) for p in range(info[t]['i'], info[t]['k'] + 1))
                hx.append((s, 'FC', VT[t] - D))
            else:
                hx.append((s, 'FR', VB[dct['bot']]))
    for p in range(1, N + 1):
        prods = cfg['producers'][(p, p)]
        top = [x for (t, x) in prods if t == 'top']
        bot = [x for (t, x) in prods if t == 'bot']
        if top and not bot:
            hx.append(((p, p), 'FC', VT[top[0]]))
        if bot and not top:
            hx.append(((p, p), 'FR', VB[bot[0]]))
    return cols, info, psiT, psiB, hx


def check(cfg, alpha, F, q, starts=4, seed=0):
    N = cfg['N']
    best, _ = solve_config(cfg, alpha, F, q, starts=starts, seed=seed)
    if best is None:
        return None
    R, sol = best
    ths = [float(x) for x in feed_roots(alpha, F, q)]
    P = [float(x) for x in petlyuk_peaks(alpha, F, q)]
    out = []
    for j in range(1, N):
        cols, info, psiT, psiB, hx = psi_terms(cfg, sol, alpha, F, q, ths[j - 1])
        inU = lambda s: s[1] <= j
        FRU = sum(v for (s, t, v) in hx if t == 'FR' and inU(s))
        FCc = sum(v for (s, t, v) in hx if t == 'FC' and not inU(s))
        Ej = [s for s in cols if s[0] <= j < s[1] and info[s]['k'] <= j]
        rhs = FRU + FCc + sum(psiT[s] for s in Ej)
        lhs = R - (P[j - 1] - (1 - q))
        out.append((j, lhs, rhs, FRU, FCc, {s: psiT[s] for s in Ej}))
    return R, out


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    rng = random.Random(5)
    C = all_configurations(N)
    worst = 0.0
    for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 30):
        c = rng.choice(C)
        rat = [rng.uniform(1.1, 4) for _ in range(N - 1)]
        al = [1.0]
        for r in reversed(rat):
            al.insert(0, al[0] * r)
        w = [rng.random() + 0.05 for _ in range(N)]
        z = [x / sum(w) for x in w]
        q = rng.choice([1.0, rng.random()])
        res = check(c, al, z, q)
        if res is None:
            continue
        R, out = res
        for (j, lhs, rhs, FRU, FCc, ps) in out:
            worst = max(worst, abs(lhs - rhs))
            neg = {f'{s}': round(v, 4) for s, v in ps.items() if v < -1e-9}
            print(c['id'], j, f'{lhs:.6f}', f'{rhs:.6f}', 'FRU', round(FRU, 4), 'FCc', round(FCc, 4), 'negPsi', neg)
    print('max |lhs-rhs| =', worst)
