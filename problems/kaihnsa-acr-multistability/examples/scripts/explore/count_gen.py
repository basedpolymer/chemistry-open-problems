#!/usr/bin/env python3
"""count_gen.py struct seed nsamples -- for each valid realizable subspace (valid_<struct>.json), sample x* and
count isolated positive steady states in the class of x* (multistart Newton in t); degenerate subspaces (class
Jacobian singular at x*) are reported separately."""
import numpy as np, json, sys
from fractions import Fraction
from scipy.optimize import fsolve, linprog
STRUCT = {
    'S3': [(0, 0, 0, 0, 1), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    'S7': [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)],
    'S9': [(0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 0, 0)],
}
name = sys.argv[1]; rng = np.random.default_rng(int(sys.argv[2])); NS = int(sys.argv[3])
REAC = [np.array(r, float) for r in STRUCT[name]]
n = len(REAC[0])
u, s, vt = np.linalg.svd(np.array([REAC[a] - REAC[3] for a in range(3)]))
Vt = vt[3:].T
k = Vt.shape[1]
valid = json.load(open("valid_%s.json" % name))
hist = {}; degen = 0
kdeg = 0; stab_hits = 0
for v in valid:
    W = np.array([[float(Fraction(x)) for x in row] for row in v['W']]).T
    xs0 = np.exp(rng.uniform(-1, 1, n))
    DF = W.T @ np.diag(xs0) @ Vt
    if abs(np.linalg.det(DF)) < 1e-10:
        degen += 1
        continue
    rx = v['rx']
    G = np.array([np.array(p, float) - np.array(r, float) for r, p in rx]).T
    ridx = [[tuple(R_) for R_ in STRUCT[name]].index(tuple(r)) for r, p in rx]
    Yr = np.array([np.array(r, float) for r, p in rx]).T
    r_ = G.shape[1]
    lams = []
    for q in range(20):
        res = linprog(rng.uniform(-1, 1, r_), A_eq=G, b_eq=np.zeros(n), bounds=[(1e-3, 100)] * r_, method="highs")
        if res.status != 0:
            continue
        Nagg = np.zeros((n, 4))
        for j, a in enumerate(ridx):
            Nagg[:, a] += res.x[j] * G[:, j]
        if np.linalg.matrix_rank(Nagg, tol=1e-8 * np.max(np.abs(Nagg))) == 3:
            lams.append(res.x)
    if not lams:
        kdeg += 1
        continue
    Bs = np.linalg.svd(G)[0][:, :3]
    mx = 0
    for s_ in range(NS):
        xs = np.exp(rng.uniform(-5, 5, n))
        F = lambda t: W.T @ (xs * np.exp(Vt @ t) - xs)
        sols = [np.zeros(k)]
        for st in rng.uniform(-8, 8, size=(30, k)):
            with np.errstate(all='ignore'):
                t, info, ier, msg = fsolve(F, st, full_output=True)
                if ier == 1 and np.all(np.isfinite(t)) and np.max(np.abs(F(t))) < 1e-10 * (1 + np.max(xs)):
                    if all(np.max(np.abs(t - q)) > 1e-6 for q in sols):
                        sols.append(t)
        mx = max(mx, len(sols))
        if len(sols) >= 3:
            beta = Yr[:, 0] @ Vt
            for lam in lams:
                nst = 0
                for t in sols:
                    x = xs * np.exp(Vt @ t)
                    J = G @ np.diag(lam * np.exp(beta @ t)) @ Yr.T @ np.diag(1 / x)
                    ev = np.linalg.eigvals(Bs.T @ J @ Bs)
                    if np.min(np.abs(ev)) > 1e-7 * np.max(np.abs(ev)) and np.max(ev.real) < -1e-9 * np.max(np.abs(ev)):
                        nst += 1
                if nst >= 2:
                    stab_hits += 1
                    print("STABLE-HIT", v['W'], xs.tolist(), lam.tolist(), flush=True)
    hist[mx] = hist.get(mx, 0) + 1
    if mx >= 3:
        print(v['W'], mx, v['rx'], flush=True)
print(name, "valid subspaces", len(valid), "kinetically degenerate (rank N < 3 for all sampled fluxes)", kdeg,
      "two-stable hits", stab_hits, "degenerate (class Jacobian singular)", degen,
      "max #steady states per class (sampled):", hist, flush=True)
