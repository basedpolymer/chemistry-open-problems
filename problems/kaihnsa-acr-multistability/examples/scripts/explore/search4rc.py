#!/usr/bin/env python3
"""
search4rc.py -- numerical exploration (not a certificate): random search for bimolecular networks with n species,
exactly 4 reactant complexes, dim S = 3, structural ACR (e_i in the row space of A, rows y_a - y_4) and two
exponentially stable positive steady states in one compatibility class.

Steady states for kappa: the toric coset x = x* o exp(V w), V basis of ker A (dim n-3).  Fluxes along the coset:
lambda_r(x) = lambda_r(x*) * exp(beta . w) with beta = y_r . V (the same for all reactions since y_a - y_4 in row(A)).
For n = 4 the coset is a curve and the class condition is a 1-variable exponential sum.

Usage: python search4rc.py seed seconds n
"""
import numpy as np, itertools, sys, time

rng = np.random.default_rng(int(sys.argv[1]))
TMAX = float(sys.argv[2])
n = int(sys.argv[3]) if len(sys.argv) > 3 else 4
COMP = [np.array(c) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]


def null(M, tol=1e-9):
    u, s, vt = np.linalg.svd(M)
    r = int(np.sum(s > tol))
    return vt[r:].T


def flux_sample(G):
    """random strictly positive vector in ker G (or None) via LP-free trick: project random positive vectors."""
    from scipy.optimize import linprog
    r = G.shape[1]
    c = rng.uniform(-1, 1, r)
    res = linprog(c, A_eq=G, b_eq=np.zeros(G.shape[0]), bounds=[(1, 10) for _ in range(r)], method="highs")
    if res.status != 0:
        return None
    lam = res.x
    # randomize inside the cone
    for _ in range(3):
        c = rng.uniform(-1, 1, r)
        res2 = linprog(c, A_eq=G, b_eq=np.zeros(G.shape[0]), bounds=[(1, 10) for _ in range(r)], method="highs")
        if res2.status == 0:
            lam = lam + rng.uniform(0, 2) * res2.x
    return lam


t0 = time.time()
tried = 0
hits = 0
while time.time() - t0 < TMAX:
    idx = rng.choice(len(COMP), size=4, replace=False)
    R = [COMP[i] for i in idx]
    A = np.array([R[a] - R[3] for a in range(3)], float)
    if np.linalg.matrix_rank(A) < 3:
        continue
    V = null(A)  # n x (n-3)
    # ACR species: rows of V equal to zero
    acr = [i for i in range(n) if np.allclose(V[i], 0)]
    if not acr:
        continue
    # reactions: each reactant complex gets 1..2 products
    rx = []
    for a in range(4):
        k = rng.integers(1, 3)
        prods = rng.choice(len(COMP), size=k, replace=False)
        for p in prods:
            if np.array_equal(COMP[p], R[a]):
                continue
            rx.append((a, COMP[p]))
    if len(set(a for a, _ in rx)) < 4:
        continue
    G = np.array([p - R[a] for a, p in rx], float).T  # n x r
    if np.linalg.matrix_rank(G) != 3:
        continue
    W = null(G.T)  # conservation laws n x (n-3)
    tried += 1
    lam = flux_sample(G)
    if lam is None:
        continue
    Yr = np.array([R[a] for a, _ in rx], float).T  # n x r
    for rep in range(20):
        xs = np.exp(rng.uniform(-3, 3, n))
        lam2 = lam * np.exp(rng.uniform(-0.0, 0.0, len(lam)))
        if n == 4:
            v = V[:, 0]
            w = W[:, 0]
            beta = Yr[:, 0] @ v

            def phi(t):
                return w @ (xs * np.exp(t * v))
            ts = np.linspace(-12, 12, 4001)
            ph = np.array([phi(t) for t in ts])
            dph = np.diff(ph)
            sc = np.sum(np.sign(dph[1:]) != np.sign(dph[:-1]))
            if sc < 2:
                continue
            # pick T between the two extreme values
            ext = np.nonzero(np.sign(dph[1:]) != np.sign(dph[:-1]))[0]
            lo, hi = sorted([ph[ext[0] + 1], ph[ext[1] + 1]])
            for T in np.linspace(lo, hi, 7)[1:-1]:
                g = ph - T
                roots = []
                for k in range(len(ts) - 1):
                    if g[k] == 0 or g[k] * g[k + 1] < 0:
                        a_, b_ = ts[k], ts[k + 1]
                        for _ in range(60):
                            m = 0.5 * (a_ + b_)
                            if (phi(a_) - T) * (phi(m) - T) <= 0:
                                b_ = m
                            else:
                                a_ = m
                        roots.append(0.5 * (a_ + b_))
                nst = 0
                info = []
                for t in roots:
                    x = xs * np.exp(t * v)
                    J = G @ np.diag(lam2 * np.exp(t * beta)) @ Yr.T @ np.diag(1 / x)
                    ev = np.linalg.eigvals(J)
                    ev = ev[np.argsort(np.abs(ev))][1:]  # drop the zero eigenvalue
                    st = np.max(ev.real) < -1e-9
                    nst += st
                    info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist()))
                if nst >= 2:
                    hits += 1
                    print("HIT", [(R[a].tolist(), p.tolist()) for a, p in rx], "lam", np.round(lam2, 4).tolist(),
                          "acr", acr, "T", T, info, flush=True)
                    break
print("tried", tried, "hits", hits)
