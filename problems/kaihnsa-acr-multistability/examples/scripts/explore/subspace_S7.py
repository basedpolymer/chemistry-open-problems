#!/usr/bin/env python3
"""
subspace_S7.py -- numerical exploration for E2, structure S7 = {E, D+E, 2C, A+B}.
Enumerates stoichiometric subspaces S spanned by 3 bimolecular reaction vectors from these reactant complexes;
for each S, takes the MAXIMAL network (all generators lying in S), screens x* for classes with 3 positive steady
states (this depends only on S), then samples fluxes for two stable states.
Usage: python subspace_S7.py part nparts xsamples lamsamples seed
"""
import numpy as np, itertools, sys, time
from scipy.optimize import linprog
part, nparts, XS, LS, seed = map(int, sys.argv[1:6])
sys.argv = ["x", "4", "1", "0", "0", "1"]
exec(open("enum_S7.py").read().split("prods = [[p")[0])
rng = np.random.default_rng(seed)
gens = []
for a in range(4):
    for p in COMP:
        if not np.array_equal(p, REAC[a]):
            gens.append((a, p, (p - REAC[a]).astype(float)))
seen = set()
t0 = time.time(); nS = 0; three = 0; hits = 0
for idx, (i, j, k) in enumerate(itertools.combinations(range(len(gens)), 3)):
    if idx % nparts != part:
        continue
    M = np.array([gens[i][2], gens[j][2], gens[k][2]])
    if np.linalg.matrix_rank(M) < 3:
        continue
    W = null(M)  # 5 x 2, S^perp
    inS = [g for g in range(len(gens)) if np.max(np.abs(gens[g][2] @ W)) < 1e-9]
    key = tuple(inS)
    if key in seen:
        continue
    seen.add(key)
    if set(gens[g][0] for g in inS) != {0, 1, 2, 3}:
        continue
    G = np.array([gens[g][2] for g in inS]).T
    r = G.shape[1]
    res = linprog(np.zeros(r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, None)] * r, method="highs")
    if res.status != 0:
        continue
    nS += 1
    Yr = np.array([REAC[gens[g][0]] for g in inS], float).T
    Bs = np.linalg.svd(G)[0][:, :3]
    found3 = 0
    for s in range(XS):
        xs = np.exp(rng.uniform(-4, 4, n))
        sols = class_solutions(W, xs)
        if not sols or len(sols) < 3:
            continue
        found3 += 1
        for q in range(LS):
            res = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(rng.uniform(0, 1) < 0.5) * 0 + 1e-3, 100] and [(1e-3, 100)] * r, method="highs")
            if res.status != 0:
                continue
            lam = res.x
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
                print("HIT r=%d" % r, [(REAC[gens[g][0]].tolist(), gens[g][1].tolist()) for g in inS],
                      "W^T", np.round(W.T, 4).tolist(), "lam", np.round(lam, 4).tolist(), "xs",
                      np.round(xs, 4).tolist(), info, flush=True)
                break
        if found3 >= 3:
            break
    if found3:
        three += 1
        print("THREE r=%d" % r, [(REAC[gens[g][0]].tolist(), gens[g][1].tolist()) for g in inS], flush=True)
print("part", part, "subspaces", nS, "with 3-state classes", three, "hits", hits, "time %.0f" % (time.time() - t0),
      flush=True)
