#!/usr/bin/env python3
"""
design_fast.py -- numerical exploration (not a certificate), vectorized version of search_design6.py.
Integral-feedback design on X, Y, Z: X+Y -> 2Y (k1), Y -> W (k2) + (NR-2) extra reactions keeping
dY/dt = y*(alpha x - beta).  Enumerates all designs with <= MAXRC reactant complexes and <= MAXC complexes and
samples kappa (log-uniform) in batch; reports designs with two exponentially stable positive steady states.
Usage: python design_fast.py NR MAXRC MAXC samples seed part nparts
"""
import numpy as np, itertools, sys, time

NR, MAXRC, MAXC, SAMP = map(int, sys.argv[1:5])
rng = np.random.default_rng(int(sys.argv[5]))
part, nparts = int(sys.argv[6]), int(sys.argv[7])
free = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (2, 0, 0), (0, 0, 2), (1, 0, 1)]
cat = [(0, 1, 0), (1, 1, 0), (0, 1, 1)]
lib = [(a, b) for a in free for b in free if a != b] + [(a, b) for a in cat for b in cat if a != b]
lib += [((0, 1, 0), w) for w in free]
lib += [((1, 1, 0), (1, 0, 0)), ((1, 1, 0), (0, 0, 1)), ((1, 1, 0), (1, 0, 1))]
lib = sorted(set(lib))


def batch(rx, K):
    S = K.shape[0]
    alpha = np.zeros(S); beta = np.zeros(S)
    for j, (a, b) in enumerate(rx):
        dy = b[1] - a[1]
        if a == (1, 1, 0):
            alpha += K[:, j] * dy
        elif a == (0, 1, 0):
            beta -= K[:, j] * dy
    ok = (alpha > 0) & (beta > 0)
    if not np.any(ok):
        return 0, None
    K = K[ok]; alpha = alpha[ok]; beta = beta[ok]; S = K.shape[0]
    c = beta / alpha
    PX0 = np.zeros((S, 3)); PX1 = np.zeros((S, 2)); PZ0 = np.zeros((S, 3)); PZ1 = np.zeros((S, 2))
    for j, (a, b) in enumerate(rx):
        v = np.array(b) - np.array(a)
        coef = K[:, j] * c ** a[0]
        if a[1] == 0:
            PX0[:, a[2]] += coef * v[0]; PZ0[:, a[2]] += coef * v[2]
        else:
            PX1[:, a[2]] += coef * v[0]; PZ1[:, a[2]] += coef * v[2]
    R = np.zeros((S, 4))
    for i in range(3):
        for j in range(2):
            R[:, i + j] += PX0[:, i] * PZ1[:, j] - PZ0[:, i] * PX1[:, j]
    good = np.abs(R[:, 3]) > 1e-12 * np.max(np.abs(R), axis=1)
    if not np.any(good):
        return 0, None
    K, c, R, PX0, PX1 = K[good], c[good], R[good], PX0[good], PX1[good]
    S = K.shape[0]
    comp = np.zeros((S, 3, 3))
    comp[:, 1, 0] = 1; comp[:, 2, 1] = 1
    comp[:, :, 2] = -(R[:, :3] / R[:, 3:4])
    fin = np.all(np.isfinite(comp.reshape(S, -1)), axis=1)
    if not np.all(fin):
        K, c, R, PX0, PX1, comp = K[fin], c[fin], R[fin], PX0[fin], PX1[fin], comp[fin]
        S = K.shape[0]
        if S == 0:
            return 0, None
    zr = np.linalg.eigvals(comp)  # S x 3
    real = np.abs(zr.imag) < 1e-9 * (1 + np.abs(zr))
    zz = np.where(real, zr.real, -1.0)
    den = PX1[:, 0:1] + PX1[:, 1:2] * zz
    yy = -(PX0[:, 0:1] + PX0[:, 1:2] * zz + PX0[:, 2:3] * zz ** 2) / np.where(np.abs(den) > 1e-300, den, 1e-300)
    pos = (zz > 1e-12) & (yy > 1e-12) & real
    cand = np.nonzero(pos.sum(axis=1) >= 3)[0]  # two stable states need >= 3 (det J alternates)
    if len(cand) == 0:
        return 0, None
    best = None
    for s in cand:
        X = np.array([[c[s], yy[s, q], zz[s, q]] for q in range(3) if pos[s, q]])
        if np.max(X) / np.min(X) > 1e6:
            continue
        Js = []
        for x in X:
            fval = np.zeros(3); fsc = 0.0
            for j, (a, b) in enumerate(rx):
                mon = K[s, j] * np.prod(x ** np.array(a))
                fval += mon * (np.array(b) - np.array(a)); fsc += abs(mon)
            if np.max(np.abs(fval)) > 1e-8 * fsc:
                break
            J = np.zeros((3, 3))
            for j, (a, b) in enumerate(rx):
                a_ = np.array(a); v = np.array(b) - a_
                for i in range(3):
                    if a_[i] > 0:
                        J[:, i] += v * K[s, j] * a_[i] * np.prod(x ** (a_ - np.eye(3, dtype=int)[i]))
            Js.append(J)
        Js = np.array(Js)
        if not np.all(np.isfinite(Js)):
            continue
        ev = np.linalg.eigvals(Js)
        st = np.max(ev.real, axis=1) < -1e-10 * np.max(np.abs(ev), axis=1)
        if st.sum() >= 2:
            return 1, (K[s], X, ev)
    return 0, None


t0 = time.time()
count = 0
hits = 0
for wi, W in enumerate(free):
    base = [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), W)]
    for extra in itertools.combinations(range(len(lib)), NR - 2):
        if hash((wi, extra)) % nparts != part:
            continue
        rx = base + [lib[i] for i in extra]
        if len(set(rx)) < len(rx):
            continue
        comps = set(); reac = set()
        for a, b in rx:
            comps.add(a); comps.add(b); reac.add(a)
        if len(reac) > MAXRC or len(comps) > MAXC:
            continue
        V = np.array([np.array(b) - np.array(a) for a, b in rx])
        if np.linalg.matrix_rank(V) < 3:
            continue
        count += 1
        K = np.exp(rng.uniform(-4, 4, size=(SAMP, len(rx))))
        h, info = batch(rx, K)
        if h:
            hits += 1
            Ks, X, ev = info
            print("HIT complexes=%d reactants=%d" % (len(comps), len(reac)), rx, "k=", np.round(Ks, 4).tolist(),
                  "ss=", np.round(X, 4).tolist(), flush=True)
print("part", part, "networks tested", count, "hits", hits, "time %.0f" % (time.time() - t0), flush=True)
