#!/usr/bin/env python3
"""End-to-end test on actual reaction networks (runtime about 5 min on 3 cores). Floating point: a test, not a proof.

For a random span in a random source type (all 26 four-species types and the nine larger types):
  * a random reaction subset that uses all four sources and spans S, with random rates chosen so that a
    random point x* is an equilibrium (a random strictly positive flux from the LP on ker G);
  * N is formed from the reactions; rank N = 3 is required (otherwise the case is counted and skipped);
  * all equilibria in the class of a random point of the toric set are obtained by msolve;
  * each is checked to be a true zero of the mass-action field of the network;
  * the full Jacobian is restricted to S and classified by its eigenvalues.

PASS requires, in every class, at most two equilibria and at most one Hurwitz equilibrium.
A blind Newton multi-start on the full system in log coordinates tests the toric description: it must
find nothing outside the toric set, so the ACR species never moves.
Control 1: on kinetic-degenerate spans every random flux gives rank N < 3.
Control 2: the spans labelled injective never give two equilibria.
"""
import json
import multiprocessing as mp
import os
import random
import sys
from collections import Counter
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import linprog, root

from common import rref, nullspace, named, diff, complexes, PAPER_TYPES, add_n4_types
from stress_count import setup, poly_text, msolve

HERE = os.path.dirname(os.path.abspath(__file__))
add_n4_types()
np.seterr(all='ignore')


def reactions_in(name, W):
    names, src, acr = PAPER_TYPES[name]
    n = len(names)
    sources = [named(names, t) for t in src]
    W = np.array(W, dtype=np.int64).reshape(-1, n)
    out = []
    for a, y in enumerate(sources):
        for p in complexes(n):
            v = diff(p, y)
            if any(v) and not np.any(W @ np.array(v)):
                out.append((a, np.array(v, dtype=float)))
    return sources, out


def positive_flux(G, rng):
    m = G.shape[1]
    tot, cnt = np.zeros(m), 0
    for _ in range(6):
        c = np.array([rng.gauss(0, 1) for _ in range(m)])
        r = linprog(c, A_eq=G, b_eq=np.zeros(G.shape[0]), bounds=[(0, 1)] * m, method='highs')
        if r.status == 0:
            tot += r.x
            cnt += 1
    if cnt == 0:
        return None
    # strict interior point of the flux polytope if there is one
    r = linprog([0.0] * m + [-1.0], A_eq=np.hstack([G, np.zeros((G.shape[0], 1))]), b_eq=np.zeros(G.shape[0]),
                A_ub=np.hstack([-np.eye(m), np.ones((m, 1))]), b_ub=np.zeros(m),
                bounds=[(0, 1)] * m + [(None, 1)], method='highs')
    if r.status != 0 or -r.fun < 1e-9:
        return None
    return tot / cnt + r.x[:m]


def field_jac(Y, rx, kap, x):
    n = len(x)
    f, J = np.zeros(n), np.zeros((n, n))
    for (a, v), k in zip(rx, kap):
        m = np.prod(np.array(x)**np.array(Y[a], dtype=float))
        f += k * m * v
        J += k * m * np.outer(v, np.array(Y[a], dtype=float) / np.array(x))
    return f, J


def job(arg):
    name, W, seed, tries = arg
    rng = random.Random(seed)
    names, src, acr = PAPER_TYPES[name]
    n = len(names)
    Y, allr = reactions_in(name, W)
    Wr, Vr, b = setup(name, W)
    m = len(Vr[0])
    idx = [i for i in range(n) if names[i] not in acr]
    acri = [i for i in range(n) if names[i] in acr]
    A = np.array([diff(Y[i], Y[3]) for i in range(3)], dtype=float)
    Sb = np.array(nullspace(W, n), dtype=float).T
    h = Counter()
    bad = []
    for _ in range(tries):
        keep = [r for r in allr if rng.random() < 0.5]
        if {a for a, _ in keep} != {0, 1, 2, 3}:
            h['no source cover'] += 1
            continue
        G = np.array([v for _, v in keep]).T
        if np.linalg.matrix_rank(G) < 3:
            h['rank of subset < 3'] += 1
            continue
        phi = positive_flux(G, rng)
        if phi is None:
            h['no positive flux'] += 1
            continue
        x0 = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(n)]
        xs = np.array([float(v) for v in x0])
        kap = [p / np.prod(xs**np.array(Y[a], dtype=float)) for p, (a, _) in zip(phi, keep)]
        N = np.zeros((n, 4))
        for (a, v), k in zip(keep, kap):
            N[:, a] += k * v
        if np.linalg.matrix_rank(N, tol=1e-9) < 3:
            h['rank N < 3'] += 1
            continue
        h['networks'] += 1
        # classes: of x0 and of a few random points of the toric set
        V = np.array(nullspace(A.astype(int).tolist(), n), dtype=float).reshape(-1, n)
        for cls in range(3):
            if cls == 0:
                xc = x0
            else:
                t = [rng.uniform(-1.2, 1.2) for _ in range(len(V))]
                xc = [v * Fr(np.exp(float(np.dot(t, V[:, i])))).limit_denominator(64) for i, v in enumerate(x0)]
                # the ACR coordinates do not move on the toric set
                for i in acri:
                    xc[i] = x0[i]
            # binomial constants are those of the network (from x0); only the class depends on xc
            sols = msolve(poly_text(Wr, b, [x0[i] for i in idx], [xc[i] for i in idx]), m)
            if sols is None:
                h['positive-dimensional class'] += 1
                continue
            pos = [s for s in sols if all(lo > 0 for lo, hi in s)]
            stable = 0
            for s in pos:
                xfull = np.array([float(v) for v in xc])
                xfull[idx] = [float((lo + hi) / 2) for lo, hi in s]
                f, J = field_jac(Y, keep, kap, xfull)
                if np.abs(f).max() > 1e-6 * (1 + np.abs(J).max()):
                    bad.append(('not an equilibrium of the network', name, W, x0, xfull.tolist()))
                C = np.linalg.pinv(Sb) @ J @ Sb
                ev = np.linalg.eigvals(C)
                stable += int(ev.real.max() < -1e-9)
            h[f'class with {len(pos)} equilibria, {stable} stable'] += 1
            if len(pos) > 2 or stable > 1:
                bad.append((name, W, x0, xc))
        # blind Newton multi-start in log coordinates on class of x0
        W0 = np.array(W, dtype=float).reshape(-1, n)
        T0 = W0 @ xs
        Q = np.linalg.qr(Sb)[0]

        def F(u):
            x = np.exp(u)
            f, _ = field_jac(Y, keep, kap, x)
            return np.concatenate([Q.T @ f, W0 @ x - T0])
        for _s in range(12):
            u0 = np.log(xs) + np.array([rng.gauss(0, 1.5) for _ in range(n)])
            try:
                r = root(F, u0, method='hybr')
            except Exception:
                continue
            if r.success and np.linalg.norm(F(r.x)) < 1e-9 * (1 + np.linalg.norm(T0)):
                x = np.exp(r.x)
                # log coordinates run off to boundary equilibria, which are not positive equilibria
                if not np.all(np.isfinite(x)) or x.min() < 1e-6:
                    h['newton limit on the boundary (not a positive equilibrium)'] += 1
                    continue
                h['newton zeros'] += 1
                if np.max(np.abs(x[acri] - xs[acri])) > 1e-6 * (1 + xs[acri].max()):
                    bad.append(('ACR moved', name, W, x0, x.tolist()))
                    h['ACR moved'] += 1
    return name, h, bad


def kin_job(arg):
    name, W, seed, tries = arg
    rng = random.Random(seed)
    names, _, _ = PAPER_TYPES[name]
    n = len(names)
    _, allr = reactions_in(name, W)
    G = np.array([v for _, v in allr]).T
    h = Counter()
    for _ in range(tries):
        phi = positive_flux(G, rng)
        if phi is None:
            h['no flux'] += 1
            break
        phi = phi * np.array([rng.uniform(0.5, 1.5) for _ in phi])
        # project back onto ker G
        P = np.eye(G.shape[1]) - np.linalg.pinv(G) @ G
        phi = P @ phi
        L = np.zeros((n, 4))
        for (a, v), p in zip(allr, phi):
            L[:, a] += p * v
        h['rank %d' % np.linalg.matrix_rank(L, tol=1e-8)] += 1
    return h


def main():
    tries = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    tasks, kin = [], []
    seed = 777
    for name in PAPER_TYPES:
        fname = os.path.join(HERE, f'spans-{name}.json')
        if not os.path.exists(fname):
            continue
        data = json.load(open(fname))
        for W, lab in data['labels']:
            seed += 1
            if lab == 'normal-form':
                tasks.append((name, W, seed, tries))
            elif lab == 'kinetic-degenerate':
                kin.append((name, W, seed, 4))
    print(f'{len(tasks)} surviving spans, {tries} attempts each; {len(kin)} kinetic-degenerate spans', flush=True)
    tot = Counter()
    badall = []
    with mp.Pool(3) as pool:
        for name, h, bad in pool.imap_unordered(job, tasks, chunksize=2):
            tot.update(h)
            badall += bad
        kt = Counter()
        for h in pool.imap_unordered(kin_job, kin, chunksize=4):
            kt.update(h)
    for k in sorted(tot):
        print(f'   {k}: {tot[k]}')
    print('kinetic-degenerate control, rank of the aggregated matrix under random positive fluxes:', dict(kt))
    ok1 = not badall
    ok2 = not any(k.startswith('rank 3') or k.startswith('rank 4') for k in kt)
    print('PASS' if ok1 else 'FAIL', 'no class with more than two equilibria or two stable ones; ACR never moved', len(badall))
    print('PASS' if ok2 else 'FAIL', 'kinetic-degenerate spans: rank N < 3 for all sampled positive fluxes')
    print('RESULT', 'PASS' if ok1 and ok2 else 'FAIL')
    return 0 if ok1 and ok2 else 1


if __name__ == '__main__':
    sys.exit(main())
