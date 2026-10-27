#!/usr/bin/env python3
"""
subspace_gen.py -- numerical exploration for E2 (4 reactant complexes, dim S = 3), any structure.
Valid networks need S = Im N, i.e. every reaction vector lies in the 3-dim kinetic subspace.  Enumerates
subspaces S spanned by 3 bimolecular reaction vectors from the 4 reactant complexes, takes the maximal network
(all generators in S), screens x* for classes with >= 3 positive steady states (depends only on S), and samples
fluxes for two exponentially stable states.  The class is solved by multistart Newton in t (k = n - 3 dims).
Usage: python subspace_gen.py struct part nparts xsamples lamsamples seed
"""
import numpy as np, itertools, sys, time
from scipy.optimize import linprog, fsolve
STRUCT = {
    'S3': [(0, 0, 0, 0, 1), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    'S7': [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)],
    'S9': [(0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 0, 0)],
}
name = sys.argv[1]
part, nparts, XS, LS, seed = map(int, sys.argv[2:7])
REAC = [np.array(r, float) for r in STRUCT[name]]
n = len(REAC[0])
COMP = [np.array(c, float) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
rng = np.random.default_rng(seed)


def null(M, tol=1e-10):
    u, s, vt = np.linalg.svd(M)
    r = int(np.sum(s > tol * max(1, s[0])))
    return vt[r:].T


A = np.array([REAC[a] - REAC[3] for a in range(3)])
Vt = null(A)
k = Vt.shape[1]
gens = [(a, p, p - REAC[a]) for a in range(4) for p in COMP if not np.array_equal(p, REAC[a])]
seen = set()
t0 = time.time(); nS = 0; three = 0; hits = 0
for idx, (i, j, l) in enumerate(itertools.combinations(range(len(gens)), 3)):
    if idx % nparts != part:
        continue
    M = np.array([gens[i][2], gens[j][2], gens[l][2]])
    if np.linalg.matrix_rank(M) < 3:
        continue
    W = null(M)
    inS = tuple(g for g in range(len(gens)) if np.max(np.abs(gens[g][2] @ W)) < 1e-9)
    if inS in seen:
        continue
    seen.add(inS)
    if set(gens[g][0] for g in inS) != {0, 1, 2, 3}:
        continue
    G = np.array([gens[g][2] for g in inS]).T
    r = G.shape[1]
    res = linprog(np.zeros(r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, None)] * r, method="highs")
    if res.status != 0:
        continue
    nS += 1
    Yr = np.array([REAC[gens[g][0]] for g in inS]).T
    Bs = np.linalg.svd(G)[0][:, :3]
    f3 = 0
    for s in range(XS):
        xs = np.exp(rng.uniform(-4, 4, n))
        F = lambda t: W.T @ (xs * np.exp(Vt @ t) - xs)
        sols = [np.zeros(k)]
        for st in rng.uniform(-6, 6, size=(25, k)):
            t, info, ier, msg = fsolve(F, st, full_output=True)
            if ier == 1 and np.all(np.isfinite(t)) and np.max(np.abs(F(t))) < 1e-10 * (1 + np.max(xs)):
                if all(np.max(np.abs(t - q)) > 1e-6 for q in sols):
                    sols.append(t)
        if len(sols) < 3 or len(sols) > 6:
            continue
        # isolated solutions only: Jacobian of the class equations must be nonsingular at every solution
        okiso = True
        for t in sols:
            x = xs * np.exp(Vt @ t)
            DF = W.T @ np.diag(x) @ Vt
            if abs(np.linalg.det(DF)) < 1e-8 * np.prod(np.linalg.norm(DF, axis=0) + 1e-300):
                okiso = False
        if not okiso:
            continue
        f3 += 1
        for q in range(LS):
            res = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(1e-3, 100)] * r, method="highs")
            if res.status != 0:
                continue
            lam = res.x
            beta = Yr[:, 0] @ Vt
            nst = 0; info = []
            for t in sols:
                x = xs * np.exp(Vt @ t)
                J = G @ np.diag(lam * np.exp(beta @ t)) @ Yr.T @ np.diag(1 / x)
                ev = np.linalg.eigvals(Bs.T @ J @ Bs)
                if np.min(np.abs(ev)) < 1e-7 * np.max(np.abs(ev)):
                    nst -= 10
                nst += np.max(ev.real) < -1e-9 * np.max(np.abs(ev))
                info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist()))
            if nst >= 2:
                hits += 1
                print("HIT r=%d" % r, [(REAC[gens[g][0]].astype(int).tolist(), gens[g][1].astype(int).tolist())
                                      for g in inS], "lam", np.round(lam, 4).tolist(), "xs",
                      np.round(xs, 4).tolist(), info, flush=True)
                break
        if f3 >= 3:
            break
    if f3:
        three += 1
        print("THREE r=%d" % r, [(REAC[gens[g][0]].astype(int).tolist(), gens[g][1].astype(int).tolist())
                                 for g in inS], flush=True)
print("struct", name, "part", part, "subspaces", nS, "with 3-state classes", three, "hits", hits,
      "time %.0f" % (time.time() - t0), flush=True)
