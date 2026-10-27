#!/usr/bin/env python3
"""
general_direct.py -- numerical exploration for E2: n species, 4 reactant complexes (structures S1..S9 from
stage1_n.py).  Columns N_1..N_3 random in the bimolecular cones, x* random, c_a = m_a(x*)/m_4(x*),
N_4 = -sum c_a N_a must lie in cone(y_4).  Class of x*: t in R^k, W^T (x* o exp(V t) - x*) = 0, solved by
multistart Newton.  Records (number of steady states, number of stable ones).
Usage: python general_direct.py struct_index seed seconds
"""
import numpy as np, sys, time
from scipy.optimize import fsolve

STRUCT = [
    [(0, 0, 0, 0, 0), (0, 0, 0, 0, 1), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 0), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 1, 0, 1), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 1, 0, 1), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)],
    [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 1, 1, 0, 0), (1, 0, 1, 0, 0)],
    [(0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 0, 0)],
]
si = int(sys.argv[1])
rng = np.random.default_rng(int(sys.argv[2]))
TMAX = float(sys.argv[3])
R = [np.array(r, float) for r in STRUCT[si]]
n = len(R[0])
Y = np.array(R).T  # n x 4


def null(M, tol=1e-10):
    u, s, vt = np.linalg.svd(M)
    r = int(np.sum(s > tol * max(1, s[0])))
    return vt[r:].T


A = np.array([R[a] - R[3] for a in range(3)])
V = null(A)
k = V.shape[1]


def in_cone(w, y, tol=1e-12):
    off = [i for i in range(n) if y[i] == 0]
    if any(w[i] < -tol for i in off):
        return False
    if y.sum() == 2 and w.sum() > tol:
        return False
    return True


def rand_col(y):
    w = rng.exponential(1, n) * (rng.random(n) < 0.5)
    sup = [i for i in range(n) if y[i] > 0]
    for i in sup:
        w[i] = rng.uniform(-3, 1)
    if y.sum() == 2 and w.sum() > 0:
        i = sup[rng.integers(len(sup))]
        w[i] -= w.sum() * rng.uniform(1, 2)
    if y.sum() == 0:
        w = np.abs(w)
    return w


t0 = time.time()
hist = {}
hits = 0
tried = 0
while time.time() - t0 < TMAX:
    xs = np.exp(rng.uniform(-2, 2, n))
    m = np.exp(Y.T @ np.log(xs))
    cols = [rand_col(R[a]) for a in range(3)]
    c = m[:3] / m[3]
    N4 = -sum(c[a] * cols[a] for a in range(3))
    if not in_cone(N4, R[3]):
        continue
    N = np.column_stack(cols + [N4])
    if np.linalg.matrix_rank(N) != 3:
        continue
    tried += 1
    W = null(N.T)
    Bs = np.linalg.svd(N)[0][:, :3]

    def F(t):
        return W.T @ (xs * np.exp(V @ t) - xs)
    sols = [np.zeros(k)]
    for st in rng.uniform(-5, 5, size=(40, k)):
        try:
            t, info, ier, msg = fsolve(F, st, full_output=True)
        except Exception:
            continue
        if ier == 1 and np.all(np.isfinite(t)) and np.max(np.abs(F(t))) < 1e-10 * (1 + np.max(xs)):
            if all(np.max(np.abs(t - s)) > 1e-6 for s in sols):
                sols.append(t)
    nst = 0
    deg = False
    info = []
    for t in sols:
        x = xs * np.exp(V @ t)
        mm = np.exp(Y.T @ np.log(x))
        J = N @ np.diag(mm) @ Y.T @ np.diag(1 / x)
        JS = Bs.T @ J @ Bs
        ev = np.linalg.eigvals(JS)
        DF = W.T @ np.diag(x) @ V
        if abs(np.linalg.det(DF)) < 1e-9 * (1 + np.max(x)) ** k or np.min(np.abs(ev)) < 1e-8 * np.max(np.abs(ev)):
            deg = True
        s_ = np.max(ev.real) < -1e-9 * np.max(np.abs(ev))
        nst += s_
        info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist()))
    key = (len(sols), nst, deg)
    hist[key] = hist.get(key, 0) + 1
    if nst >= 2 and not deg:
        hits += 1
        print("HIT", np.round(N, 4).tolist(), np.round(xs, 4).tolist(), info, flush=True)
print("struct", si, "tried", tried, "hist (n_ss, n_stable, degenerate):", hist, "hits", hits)
