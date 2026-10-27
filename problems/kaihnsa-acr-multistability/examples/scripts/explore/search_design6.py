#!/usr/bin/env python3
"""
search_design6.py -- numerical exploration (not a certificate).
Integral-feedback design on species X, Y, Z: X+Y -> 2Y (k1), Y -> W (k2), plus extra reactions that keep
dY/dt = y * (affine function of x) (Y-neutral reactions, extra Y -> W' reactions, X+Y -> X).
ACR in X is then structural.  Enumerates ALL designs with a given total number of reactions whose
(reactant complexes, complexes) counts are below given thresholds, and samples kappa to find two exponentially
stable positive steady states.
Usage: python search_design6.py nreac max_reactant_complexes max_complexes samples seed [W-index]
"""
import numpy as np, itertools, sys, time

NR = int(sys.argv[1])
MAXRC = int(sys.argv[2])
MAXC = int(sys.argv[3])
SAMP = int(sys.argv[4])
rng = np.random.default_rng(int(sys.argv[5]))
free = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (2, 0, 0), (0, 0, 2), (1, 0, 1)]
cat = [(0, 1, 0), (1, 1, 0), (0, 1, 1)]
lib = [(a, b) for a in free for b in free if a != b] + [(a, b) for a in cat for b in cat if a != b]
lib += [((0, 1, 0), w) for w in free]          # extra Y -> W'
lib += [((1, 1, 0), (1, 0, 0))]                # X+Y -> X
lib += [((1, 1, 0), w) for w in [(0, 0, 1), (1, 0, 1)] ]  # X+Y -> Z, X+Y -> X+Z: change Y by -1 with xy monomial
base_W = free


def steady(rx, K):
    """rx list of (a,b); K (S,r) rates. Reaction 0 must be X+Y->2Y. Returns per-sample list of (x, stable)."""
    S = K.shape[0]
    # dY/dt = y*(alpha x - beta): alpha = sum over reactions with reactant X+Y of k*(dY), beta = -sum over Y of k*dY
    alpha = np.zeros(S)
    beta = np.zeros(S)
    for j, (a, b) in enumerate(rx):
        dy = b[1] - a[1]
        if a == (1, 1, 0):
            alpha += K[:, j] * dy
        elif a == (0, 1, 0):
            beta -= K[:, j] * dy
        else:
            assert dy == 0, (a, b)
    ok = (alpha > 0) & (beta > 0)
    c = np.where(ok, beta / np.where(alpha > 0, alpha, 1), 1.0)
    PX0 = np.zeros((S, 3)); PX1 = np.zeros((S, 2)); PZ0 = np.zeros((S, 3)); PZ1 = np.zeros((S, 2))
    for j, (a, b) in enumerate(rx):
        v = np.array(b) - np.array(a)
        coef = K[:, j] * c ** a[0]
        dz = a[2]
        if a[1] == 0:
            PX0[:, dz] += coef * v[0]; PZ0[:, dz] += coef * v[2]
        else:
            PX1[:, dz] += coef * v[0]; PZ1[:, dz] += coef * v[2]
    R = np.zeros((S, 4))
    for i in range(3):
        for j in range(2):
            R[:, i + j] += PX0[:, i] * PZ1[:, j] - PZ0[:, i] * PX1[:, j]
    res = []
    for s in np.nonzero(ok)[0]:
        co = R[s, ::-1]
        nz = np.nonzero(np.abs(co) > 1e-14 * np.max(np.abs(co)))[0]
        if len(nz) == 0:
            continue
        z = np.roots(co[nz[0]:])
        z = z[(np.abs(z.imag) < 1e-9 * (1 + np.abs(z))) & (z.real > 1e-12)].real
        if len(z) < 2:
            continue
        sols = []
        for zz in z:
            den = PX1[s, 0] + PX1[s, 1] * zz
            if abs(den) < 1e-14:
                continue
            y = -(PX0[s, 0] + PX0[s, 1] * zz + PX0[s, 2] * zz * zz) / den
            if y <= 1e-12:
                continue
            x = np.array([c[s], y, zz])
            J = np.zeros((3, 3)); f = np.zeros(3)
            for j, (a, b) in enumerate(rx):
                a_ = np.array(a); v = np.array(b) - a_
                f += K[s, j] * np.prod(x ** a_) * v
                for i in range(3):
                    if a_[i] > 0:
                        J[:, i] += v * K[s, j] * a_[i] * np.prod(x ** (a_ - np.eye(3, dtype=int)[i]))
            sc = np.abs(K[s]).max() * (1 + np.max(x)) ** 2
            if np.max(np.abs(f)) > 1e-7 * sc:
                continue
            ev = np.linalg.eigvals(J)
            sols.append((x, bool(np.max(ev.real) < -1e-10 * np.max(np.abs(ev)))))
        if sum(st for _, st in sols) >= 2:
            res.append((s, sols))
    return res


t0 = time.time()
count = 0
hits = 0
Wsel = [int(sys.argv[6])] if len(sys.argv) > 6 else range(len(base_W))
for wi in Wsel:
    W = base_W[wi]
    base = [((1, 1, 0), (0, 2, 0)), ((0, 1, 0), W)]
    for extra in itertools.combinations(range(len(lib)), NR - 2):
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
        res = steady(rx, K)
        if res:
            hits += 1
            s, sols = res[0]
            print("HIT complexes=%d reactants=%d" % (len(comps), len(reac)), rx, "k=",
                  np.round(K[s], 4).tolist(), "ss=", [(np.round(x, 4).tolist(), st) for x, st in sols], flush=True)
print("networks tested", count, "hits", hits, "time %.0f" % (time.time() - t0))
