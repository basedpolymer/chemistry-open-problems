#!/usr/bin/env python3
"""
enum_S7.py -- numerical exploration for E2: structure S7, reactant complexes {E, D+E, 2C, A+B} (species A..E),
exactly R reactions (R = 4: one per reactant complex; R = 5: one doubled).  For each network with rank 3 and a
positive flux vector, sample steady states x* (any x* > 0 is a steady state for suitable kappa) and solve the
class of x* exactly enough: with u = exp(t1), w = exp(t2/2) the toric coset is
x = (a w^2/u, b u, c w, d, e w^2) and the two conservation laws give two quadratics in u; resultant -> w.
Reports networks with a class containing >= 2 exponentially stable steady states.
Usage: python enum_S7.py R samples seed [part nparts]
"""
import numpy as np, itertools, sys, time

RN = int(sys.argv[1])
SAMP = int(sys.argv[2])
rng = np.random.default_rng(int(sys.argv[3]))
part = int(sys.argv[4]) if len(sys.argv) > 4 else 0
nparts = int(sys.argv[5]) if len(sys.argv) > 5 else 1
n = 5
COMP = [np.array(c) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
REAC = [np.array(r) for r in [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)]]
Vt = np.array([[-1, 1, 0, 0, 0], [1, 0, 0.5, 0, 1]], float).T  # x = x* o exp(V t)


def null(M, tol=1e-10):
    u, s, vt = np.linalg.svd(M)
    r = int(np.sum(s > tol * max(1, s[0])))
    return vt[r:].T


def class_solutions(W, xs):
    """all positive steady states in the class of xs; returns list of x."""
    a, b, c, d, e = xs
    K = W.T @ xs  # 2-vector
    # row l: W_A a w^2/u + W_B b u + W_C c w + W_D d + W_E e w^2 = K_l ; times u:
    # (W_B b) u^2 + (W_C c w + W_E e w^2 + W_D d - K_l) u + W_A a w^2 = 0
    P = []
    for l in range(2):
        A2 = np.array([W[1, l] * b])
        A1 = np.array([W[3, l] * d - K[l], W[2, l] * c, W[4, l] * e])
        A0 = np.array([0, 0, W[0, l] * a])
        P.append((A2, A1, A0))
    from numpy.polynomial import polynomial as Pl
    (a2, a1, a0), (b2, b1, b0) = P
    t1 = Pl.polysub(Pl.polymul(a2, b0), Pl.polymul(b2, a0))
    t2 = Pl.polysub(Pl.polymul(a2, b1), Pl.polymul(b2, a1))
    t3 = Pl.polysub(Pl.polymul(a1, b0), Pl.polymul(b1, a0))
    res = Pl.polysub(Pl.polymul(t1, t1), Pl.polymul(t2, t3))
    res = np.trim_zeros(res, 'b')
    if len(res) < 2:
        return None
    out = []
    for w in Pl.polyroots(res):
        if abs(w.imag) > 1e-7 * (1 + abs(w)) or w.real <= 1e-9:
            continue
        w = w.real
        den = Pl.polyval(w, t2)
        if abs(den) < 1e-14:
            continue
        u = -Pl.polyval(w, t1) / den
        if u <= 1e-12:
            continue
        x = np.array([a * w * w / u, b * u, c * w, d, e * w * w])
        if np.max(np.abs(W.T @ x - K)) > 1e-7 * (1 + np.max(np.abs(K)) + np.max(x)):
            continue
        if all(np.max(np.abs(x - y) / (1 + np.abs(y))) > 1e-7 for y in out):
            out.append(x)
    return out


prods = [[p for p in COMP if not np.array_equal(p, r)] for r in REAC]
t0 = time.time()
count = 0
hits = 0
hist = {}
choices = []
if RN == 4:
    it = itertools.product(range(20), repeat=4)
    mk = lambda tup: [(a, prods[a][tup[a]]) for a in range(4)]
else:
    it = ((dbl, pp, tup) for dbl in range(4) for pp in itertools.combinations(range(20), 2)
          for tup in itertools.product(range(20), repeat=3))

    def mk(item):
        dbl, pp, tup = item
        rx = []
        j = 0
        for a in range(4):
            if a == dbl:
                rx += [(a, prods[a][pp[0]]), (a, prods[a][pp[1]])]
            else:
                rx.append((a, prods[a][tup[j]]))
                j += 1
        return rx
for idx, item in enumerate(it):
    if idx % nparts != part:
        continue
    rx = mk(item)
    G = np.array([p - REAC[a] for a, p in rx], float).T
    if np.linalg.matrix_rank(G) != 3:
        continue
    L = null(G)
    r = G.shape[1]
    # positive flux: need lambda > 0 in ker G
    if L.shape[1] == 1:
        lam = L[:, 0] * np.sign(L[0, 0])
        if not np.all(lam > 1e-12):
            continue
        lams = [lam]
    elif L.shape[1] == 2:
        ang = np.arctan2(L[:, 1], L[:, 0])
        cand = []
        for a_ in ang:
            for d_ in (a_ + np.pi / 2, a_ - np.pi / 2):
                th = np.array([np.cos(d_), np.sin(d_)])
                if np.all(L @ th >= -1e-12):
                    cand.append(th)
        rays = None
        for i_ in range(len(cand)):
            for j_ in range(i_ + 1, len(cand)):
                if np.all(L @ (cand[i_] + cand[j_]) > 1e-9):
                    rays = (L @ cand[i_], L @ cand[j_])
        if rays is None:
            continue
        lams = rays
    else:
        from scipy.optimize import linprog
        res = linprog(np.zeros(r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, None)] * r, method="highs")
        if res.status != 0:
            continue
        lams = None
    count += 1
    W = null(G.T)
    Yr = np.array([REAC[a] for a, _ in rx], float).T
    Bs = np.linalg.svd(G)[0][:, :3]
    best = 0
    for s in range(SAMP):
        if lams is None:
            res = linprog(rng.uniform(-1, 1, r), A_eq=G, b_eq=np.zeros(n), bounds=[(1, 30)] * r, method="highs")
            lam = res.x
        elif len(lams) == 2:
            al = rng.uniform(0.02, 0.98)
            lam = al * lams[0] / np.max(lams[0]) + (1 - al) * lams[1] / np.max(lams[1])
        else:
            lam = lams[0]
        xs = np.exp(rng.uniform(-3, 3, n))
        sols = class_solutions(W, xs)
        if not sols or len(sols) < 2:
            continue
        nst = 0
        info = []
        for x in sols:
            t = np.linalg.lstsq(Vt, np.log(x / xs), rcond=None)[0]
            beta = Yr[:, 0] @ Vt
            J = G @ np.diag(lam * np.exp(beta @ t)) @ Yr.T @ np.diag(1 / x)
            ev = np.linalg.eigvals(Bs.T @ J @ Bs)
            st = np.max(ev.real) < -1e-9 * np.max(np.abs(ev))
            nst += st
            info.append((np.round(x, 4).tolist(), np.round(ev, 4).tolist()))
        key = (len(sols), nst)
        hist[key] = hist.get(key, 0) + 1
        if nst >= 2:
            hits += 1
            print("HIT", [(REAC[a].tolist(), p.tolist()) for a, p in rx], "lam", np.round(lam, 4).tolist(),
                  "xs", np.round(xs, 4).tolist(), info, flush=True)
            break
print("networks", count, "hist", hist, "hits", hits, "time %.0f" % (time.time() - t0))
