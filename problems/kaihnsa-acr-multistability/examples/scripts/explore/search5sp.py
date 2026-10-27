#!/usr/bin/env python3
"""
search5sp.py -- numerical exploration (not a certificate) for E2 with n = 5 species, 4 reactant complexes, dim 3,
using the 8 reactant-set structures with an ACR species found by stage1_n.py.
Random reactions (1-3 per reactant complex), random flux lambda in ker Gamma (positive), random x*.
Steady states of the class of x*: t in R^2 with W^T (x* o exp(V t) - x*) = 0 (W = conservation laws).
Usage: python search5sp.py seed seconds
"""
import numpy as np, itertools, sys, time
from scipy.optimize import linprog, fsolve

rng = np.random.default_rng(int(sys.argv[1]))
TMAX = float(sys.argv[2])
n = 5
COMP = [np.array(c) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
STRUCT = [
    [(0, 0, 0, 0, 0), (0, 0, 0, 0, 1), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 0), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 1, 0, 1), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 1, 0, 1), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 1, 1, 0, 0), (1, 0, 1, 0, 0)],
]


def null(M, tol=1e-9):
    u, s, vt = np.linalg.svd(M)
    r = int(np.sum(s > tol * max(1, s[0])))
    return vt[r:].T


t0 = time.time()
tried = 0
hits = 0
multi = 0
while time.time() - t0 < TMAX:
    R = [np.array(r) for r in STRUCT[rng.integers(len(STRUCT))]]
    A = np.array([R[a] - R[3] for a in range(3)], float)
    V = null(A)
    rx = []
    for a in range(4):
        k = rng.integers(1, 4)
        for p in rng.choice(len(COMP), size=k, replace=False):
            if not np.array_equal(COMP[p], R[a]):
                rx.append((a, COMP[p]))
    if len(set(a for a, _ in rx)) < 4:
        continue
    G = np.array([p - R[a] for a, p in rx], float).T
    if np.linalg.matrix_rank(G) != 3:
        continue
    r = G.shape[1]
    res = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, 20)] * r, method="highs")
    if res.status != 0:
        continue
    tried += 1
    W = null(G.T)
    Yr = np.array([R[a] for a, _ in rx], float).T
    for rep in range(10):
        lam = res.x * np.exp(rng.uniform(-0.5, 0.5, r))
        # re-project lam onto ker G keeping positivity
        res2 = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, 20)] * r, method="highs")
        if res2.status == 0:
            lam = rng.uniform(0.1, 1) * res.x + rng.uniform(0.1, 1) * res2.x
        else:
            lam = res.x
        xs = np.exp(rng.uniform(-2, 2, n))
        Tv = W.T @ xs

        def F(t):
            return W.T @ (xs * np.exp(V @ t)) - Tv
        sols = [np.zeros(2)]
        for start in rng.uniform(-6, 6, size=(60, 2)):
            t, info, ier, msg = fsolve(F, start, full_output=True)
            if ier == 1 and np.max(np.abs(F(t))) < 1e-9 * (1 + np.max(np.abs(Tv))):
                if all(np.max(np.abs(t - s)) > 1e-5 for s in sols):
                    sols.append(t)
        if len(sols) < 3:
            continue
        multi += 1
        nst = 0
        info = []
        beta = Yr[:, 0] @ V
        Bs = np.linalg.svd(G)[0][:, :3]  # orthonormal basis of S
        for t in sols:
            x = xs * np.exp(V @ t)
            J = G @ np.diag(lam * np.exp(beta @ t)) @ Yr.T @ np.diag(1 / x)
            JS = Bs.T @ J @ Bs
            ev = np.linalg.eigvals(JS)
            DF = W.T @ np.diag(x) @ V
            if abs(np.linalg.det(DF)) < 1e-8 * (1 + np.max(x)) ** 2 or np.min(np.abs(ev)) < 1e-7 * np.max(np.abs(ev)):
                st = False
                nst -= 10
            else:
                st = np.max(ev.real) < -1e-9 * np.max(np.abs(ev))
            nst += st
            info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist(), bool(st)))
        if nst >= 2:
            hits += 1
            print("HIT", [(R[a].tolist(), p.tolist()) for a, p in rx], "lam", np.round(lam, 4).tolist(), info,
                  flush=True)
            break
print("tried", tried, "classes with >=3 steady states", multi, "hits", hits)
