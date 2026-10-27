#!/usr/bin/env python3
"""
rand_S7.py -- numerical exploration for E2: random reaction networks with reactant complexes {E, D+E, 2C, A+B}
and R reactions (Rmin..Rmax), flux sampled in ker Gamma (LP), x* random; class solved with the resultant solver
of enum_S7.py.  Prints hits (two exponentially stable states in one class) with counts.
Usage: python rand_S7.py seed seconds Rmin Rmax
"""
import numpy as np, itertools, sys, time
from scipy.optimize import linprog
sys.argv_saved = sys.argv
seed, TMAX, RMIN, RMAX = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
sys.argv = ["x", "4", "1", "0", "0", "1"]
src = open("enum_S7.py").read().split("prods = [[p")[0]
exec(src)
rng = np.random.default_rng(seed)
prods = [[p for p in COMP if not np.array_equal(p, r)] for r in REAC]
t0 = time.time(); hits = 0; tried = 0
while time.time() - t0 < TMAX:
    R_ = rng.integers(RMIN, RMAX + 1)
    alloc = [1, 1, 1, 1]
    for _ in range(R_ - 4):
        alloc[rng.integers(4)] += 1
    rx = []
    for a in range(4):
        for p in rng.choice(20, size=alloc[a], replace=False):
            rx.append((a, prods[a][p]))
    G = np.array([p - REAC[a] for a, p in rx], float).T
    if np.linalg.matrix_rank(G) != 3:
        continue
    r = G.shape[1]
    res = linprog(np.zeros(r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, None)] * r, method="highs")
    if res.status != 0:
        continue
    tried += 1
    W = null(G.T)
    Yr = np.array([REAC[a] for a, _ in rx], float).T
    Bs = np.linalg.svd(G)[0][:, :3]
    for s in range(15):
        res = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, 50)] * r, method="highs")
        lam = res.x
        xs = np.exp(rng.uniform(-3, 3, n))
        sols = class_solutions(W, xs)
        if not sols or len(sols) < 3:
            continue
        nst = 0; info = []
        for x in sols:
            t = np.linalg.lstsq(Vt, np.log(x / xs), rcond=None)[0]
            beta = Yr[:, 0] @ Vt
            J = G @ np.diag(lam * np.exp(beta @ t)) @ Yr.T @ np.diag(1 / x)
            ev = np.linalg.eigvals(Bs.T @ J @ Bs)
            nst += np.max(ev.real) < -1e-9 * np.max(np.abs(ev))
            info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist()))
        if nst >= 2:
            hits += 1
            comps = set(tuple(REAC[a]) for a, _ in rx) | set(tuple(p) for _, p in rx)
            print("HIT R=%d complexes=%d" % (r, len(comps)), [(REAC[a].tolist(), p.tolist()) for a, p in rx],
                  "lam", np.round(lam, 4).tolist(), "xs", np.round(xs, 4).tolist(), info, flush=True)
            break
print("tried", tried, "hits", hits)
