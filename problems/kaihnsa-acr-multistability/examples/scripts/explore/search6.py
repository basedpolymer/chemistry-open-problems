#!/usr/bin/env python3
"""
search6.py -- numerical exploration (not a certificate): 3 species, 6 reactions, exactly 5 reactant complexes,
structural ACR of "pair" type: two reactant complexes a, b with y_a - y_b = e_j, and the reactions of the other
three reactant complexes have linearly dependent reaction vectors (rank <= 2), so that lambda_a ~ lambda_b on
ker N for every kappa and x_j is constant on all positive steady states.
Reports networks with two exponentially stable positive steady states, with the number of complexes.
Usage: python search6.py seed seconds
"""
import numpy as np, itertools, sys, time

rng = np.random.default_rng(int(sys.argv[1]))
TMAX = float(sys.argv[2])
COMP = [np.array(c) for c in itertools.product(range(3), repeat=3) if sum(c) <= 2]
PAIRS = [(a, b, j) for a in range(10) for b in range(10) for j in range(3)
         if np.array_equal(COMP[a] - COMP[b], np.eye(3, dtype=int)[j])]


def analyse(N, Y, want_acr):
    """N 3x5 real; returns list of (x, eig) of steady states for a level with 3 roots, best over levels."""
    u, s, vt = np.linalg.svd(N)
    if s[2] < 1e-9 * s[0]:
        return None
    K = vt[3:].T
    ang = np.arctan2(K[:, 1], K[:, 0])
    cand = []
    for a in ang:
        for d in (a + np.pi / 2, a - np.pi / 2):
            th = np.array([np.cos(d), np.sin(d)])
            if np.all(K @ th >= -1e-12):
                cand.append(th)
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
    c = np.linalg.svd(A)[2][-1]
    u0 = K @ r1
    v0 = K @ r2 - u0
    ss = np.linspace(1e-6, 1 - 1e-6, 6001)
    mu = np.outer(1 - ss, u0) + np.outer(ss, K @ r2)
    if np.any(mu <= 0):
        mu = np.maximum(mu, 1e-300)
    g = np.log(mu) @ c
    gp = (v0 / mu) @ c
    scale = (np.abs(v0) / mu) @ np.abs(c)
    good = np.abs(gp) > 1e-7 * scale
    if not np.any(good):
        return None
    sg = np.sign(gp[good])
    idxg = np.nonzero(good)[0]
    chg = np.nonzero(sg[1:] != sg[:-1])[0]
    if len(chg) < 2 or len(chg) > 3:
        return None
    ch = idxg[chg]
    out = []
    for k in range(len(ch) - 1):
        lo, hi = sorted([g[ch[k] + 1], g[ch[k + 1] + 1]])
        for L in np.linspace(lo, hi, 9)[1:-1]:
            h = g - L
            roots = []
            for q in np.nonzero(np.sign(h[1:]) != np.sign(h[:-1]))[0]:
                a_, b_ = ss[q], ss[q + 1]
                for _ in range(50):
                    m = 0.5 * (a_ + b_)
                    hm = np.log((1 - m) * u0 + m * (K @ r2)) @ c - L
                    ha = np.log((1 - a_) * u0 + a_ * (K @ r2)) @ c - L
                    if ha * hm <= 0:
                        b_ = m
                    else:
                        a_ = m
                roots.append(0.5 * (a_ + b_))
            if len(roots) < 3:
                continue
            beta = np.exp(L * c / (c @ c))
            Np = N * beta  # column scaling: level shift
            sts = []
            for sr in roots:
                muv = (1 - sr) * u0 + sr * (K @ r2)
                if np.any(muv <= 0):
                    continue
                rhs = np.log(muv) - np.log(beta)
                if not np.all(np.isfinite(rhs)):
                    continue
                M = np.hstack([Y.T, -np.ones((5, 1))])
                sol = np.linalg.lstsq(M, rhs, rcond=None)[0]
                if not np.all(np.isfinite(sol)) or np.max(np.abs(sol[:3])) > 300:
                    continue
                x = np.exp(sol[:3])
                m = np.exp(Y.T @ np.log(x))
                f = Np @ m
                if np.max(np.abs(f)) > 1e-6 * np.max(np.abs(Np) @ m):
                    continue
                J = Np @ np.diag(m) @ Y.T @ np.diag(1 / x)
                ev = np.linalg.eigvals(J)
                sts.append((x, ev))
            nst = sum(np.max(e.real) < -1e-9 for _, e in sts)
            if nst >= 2:
                xs = np.array([x for x, _ in sts])
                acr = [j for j in range(3) if np.ptp(np.log(xs[:, j])) < 1e-6]
                if (not want_acr) or acr:
                    out.append((L, sts, acr, Np))
                    return out
    return out or None


t0 = time.time()
hits = 0
tried = 0
seen = set()
while time.time() - t0 < TMAX:
    a, b, j = PAIRS[rng.integers(len(PAIRS))]
    others = rng.choice([i for i in range(10) if i not in (a, b)], size=3, replace=False)
    reac = [a, b] + list(others)
    dbl = rng.integers(5)
    rx = []
    for idx, rc in enumerate(reac):
        k = 2 if idx == dbl else 1
        prods = rng.choice([p for p in range(10) if p != rc], size=k, replace=False)
        for p in prods:
            rx.append((idx, rc, p))
    # dependence condition on the other three complexes' reaction vectors
    V_oth = np.array([COMP[p] - COMP[rc] for idx, rc, p in rx if idx >= 2])
    if np.linalg.matrix_rank(V_oth) > 2:
        continue
    G = np.array([COMP[p] - COMP[rc] for idx, rc, p in rx]).T
    if np.linalg.matrix_rank(G) < 3:
        continue
    tried += 1
    Y = np.array([COMP[rc] for rc in reac], float).T
    for rep in range(30):
        kap = np.exp(rng.uniform(-3, 3, len(rx)))
        N = np.zeros((3, 5))
        for (idx, rc, p), kk in zip(rx, kap):
            N[:, idx] += kk * (COMP[p] - COMP[rc])
        res = analyse(N, Y, True)
        if res:
            L, sts, acr, Np = res[0]
            comps = set()
            for idx, rc, p in rx:
                comps.add(rc)
                comps.add(p)
            key = tuple(sorted((rc, p) for idx, rc, p in rx))
            hits += 1
            if key not in seen:
                seen.add(key)
                print("HIT complexes=%d" % len(comps), [(COMP[rc].tolist(), COMP[p].tolist()) for idx, rc, p in rx],
                      "kappa", np.round(kap, 4).tolist(), "acr", acr,
                      "x", [np.round(x, 4).tolist() for x, _ in sts], flush=True)
            break
print("tried", tried, "hits", hits)
