#!/usr/bin/env python3
"""
agg5_sample.py -- numerical exploration: 3 species, exactly 5 reactant complexes, columns of N = random nonnegative
combinations of the bimolecular reaction vectors of each reactant complex (i.e. any number of reactions).
Counts how often g' changes sign >= 2 times on the cone ker N cap R^5_{>0} (needed for 3 positive steady states).
Usage: python agg5_sample.py seed samples_per_set kmax
"""
import numpy as np, itertools, sys

rng = np.random.default_rng(int(sys.argv[1]))
NS = int(sys.argv[2])
KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
COMP = [c for c in itertools.product(range(3), repeat=3) if sum(c) <= 2]
PERMS = list(itertools.permutations(range(3)))


def canon():
    out = []
    for sub in itertools.combinations(COMP, 5):
        key = tuple(sorted(sub))
        if key == min(tuple(sorted(tuple(c[p[i]] for i in range(3)) for c in sub)) for p in PERMS):
            out.append(key)
    return out


def analyse(N, Y):
    u, s, vt = np.linalg.svd(N)
    if s[2] < 1e-9:
        return None
    K = vt[3:].T  # 5x2
    # cone: angles
    ang = np.arctan2(K[:, 1], K[:, 0])
    # feasible theta: K theta > 0 ; test directions
    cand = []
    for a in ang:
        for d in (a + np.pi / 2, a - np.pi / 2):
            th = np.array([np.cos(d), np.sin(d)])
            if np.all(K @ th >= -1e-12):
                cand.append(th)
    if len(cand) < 2:
        return None
    best = None
    for i in range(len(cand)):
        for j in range(i + 1, len(cand)):
            th = cand[i] + cand[j]
            if np.all(K @ th > 1e-9):
                best = (cand[i], cand[j])
    if best is None:
        return None
    r1, r2 = best
    A = np.vstack([Y, np.ones(5)])
    uu, ss, vv = np.linalg.svd(A)
    c = vv[-1]
    s_ = np.linspace(1e-7, 1 - 1e-7, 4001)
    mu = np.outer(1 - s_, K @ r1) + np.outer(s_, K @ r2)
    v = K @ r2 - K @ r1
    gp = (mu ** -1 * v * c).sum(axis=1)
    sg = np.sign(gp)
    sg = sg[np.abs(gp) > 1e-10]
    return int(np.sum(sg[1:] != sg[:-1]))


subs = canon()
tot = 0
hist = {}
for R in subs:
    Y = np.array(R, float).T
    gens = [[np.array(p) - np.array(r) for p in COMP if p != r] for r in R]
    for _ in range(NS):
        cols = []
        ks = []
        dbl = rng.integers(5)
        for gi, g in enumerate(gens):
            k = (rng.integers(1, KMAX + 1) if KMAX > 0 else (2 if gi == dbl else 1))
            ks.append(int(k))
            idx = rng.choice(len(g), size=k, replace=False)
            w = np.exp(rng.uniform(-3, 3, k))
            cols.append(sum(w[i] * g[idx[i]] for i in range(k)))
        N = np.array(cols, float).T
        sc = analyse(N, Y)
        if sc is None:
            continue
        tot += 1
        hist[sc] = hist.get(sc, 0) + 1
        if sc >= 2:
            print("SC>=2", R, "ks", ks, np.round(N, 4).tolist(), flush=True)
print("consistent samples", tot, "hist", hist)
