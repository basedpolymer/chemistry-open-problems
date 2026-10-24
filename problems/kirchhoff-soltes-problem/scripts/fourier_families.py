#!/usr/bin/env python3
"""fourier_families.py -- floating-point scan (EXPLORATION, not a proof) of large vertex-transitive /
two-orbit families via characters and the local deletion formula (Proposition 9 of result.md).

Families
  tori4 : Cay(Z^2/Lambda, {+-e1, +-e2}) for every sublattice Lambda of index n (Hermite normal form rows
          (a,b),(0,d), ad = n, 0 <= b < d), with e1, e2, 2e1, 2e2, e1+e2, e1-e2 not in Lambda (simple
          4-regular graph).  Includes the cyclic quotients (circulants C_n(1,k)).
  tri6  : Cay(Z^2/Lambda, {+-e1, +-e2, +-(e1+e2)}) with the six elements distinct and non-zero.
  gp    : generalized Petersen graphs GP(n,k), 1 <= k < n/2, Delta at an outer and at an inner vertex.
For each graph prints nothing unless |Delta|/Kf < THRESH; at the end prints counts, the smallest |Delta|/Kf
and the near-miss list.  Usage: fourier_families.py FAMILY NMIN NMAX [THRESH]
Validation mode: fourier_families.py validate  (compares with dense numpy computations; PASS/FAIL).
"""
import math
import sys

import numpy as np


def local_delta(n, nbrs, gam, gam2, trG):
    """Proposition 9 in floating point. nbrs: list of neighbour labels; gam(x,y), gam2(x,y) callables with
    vertex labels; v is label 'v'."""
    gvv = gam('v', 'v')
    g2vv = gam2('v', 'v')
    d = len(nbrs)
    Z = np.zeros((d + 1, d + 1))
    W = np.zeros((d + 1, d + 1))
    gav = [gam(a, 'v') for a in nbrs]
    g2av = [gam2(a, 'v') for a in nbrs]
    for i, a in enumerate(nbrs):
        for j, b in enumerate(nbrs):
            Z[i, j] = gam(a, b) - gav[i] - gav[j] + gvv - (1.0 if i == j else 0.0)
            W[i, j] = gam2(a, b) - g2av[i] - g2av[j] + g2vv + n * (gav[i] - gvv) * (gav[j] - gvv)
        Z[i, d] = Z[d, i] = n * (gvv - gav[i])
        W[i, d] = W[d, i] = -n * (g2av[i] - g2vv) - n * n * gvv * (gav[i] - gvv)
    Z[d, d] = n - 1 + n * n * gvv
    W[d, d] = n * n * g2vv + n ** 3 * gvv ** 2
    tr = np.trace(np.linalg.solve(Z, W))
    Rt = trG + n * gvv
    return n * trG, (n - 1) * (Rt - 1 - tr)


def lattice_quotient(a, b, d, gens):
    n = a * d
    i = np.arange(a)[:, None]
    j = np.arange(d)[None, :]
    th1 = ((i * d - b * j) / n).ravel()
    th2 = np.broadcast_to(j / d, (a, d)).ravel()
    lam = np.zeros(n)
    for s in gens:
        lam += 1 - np.cos(2 * np.pi * (th1 * s[0] + th2 * s[1]))
    return th1, th2, lam


def in_lattice(x, a, b, d):
    if x[0] % a:
        return False
    p = x[0] // a
    return (x[1] - p * b) % d == 0


def cayley_lattice_delta(a, b, d, gens):
    n = a * d
    th1, th2, lam = lattice_quotient(a, b, d, gens)
    mask = np.ones(n, bool)
    mask[0] = False                     # theta = 0 is index 0 (i=j=0)
    if np.any(lam[mask] < 1e-12):
        return None
    inv1 = np.zeros(n)
    inv1[mask] = 1 / lam[mask]
    inv2 = inv1 ** 2
    cache = {}

    def G(x, which):
        key = (x[0], x[1], which)
        if key not in cache:
            c = np.cos(2 * np.pi * (th1 * x[0] + th2 * x[1]))
            cache[key] = float(np.dot(c, inv1 if which == 1 else inv2)) / n
        return cache[key]
    lab = {'v': (0, 0)}
    for idx, s in enumerate(gens):
        lab[idx] = s

    def gam(p, q):
        x, y = lab[p], lab[q]
        return G((x[0] - y[0], x[1] - y[1]), 1)

    def gam2(p, q):
        x, y = lab[p], lab[q]
        return G((x[0] - y[0], x[1] - y[1]), 2)
    trG = n * G((0, 0), 1)
    return local_delta(n, list(range(len(gens))), gam, gam2, trG)


def gp_delta(n, k):
    """Delta at outer vertex (0,0) and inner vertex (1,0) of GP(n,k)."""
    th = 2 * np.pi * np.arange(n) / n
    A11 = 3 - 2 * np.cos(th)
    A22 = 3 - 2 * np.cos(k * th)
    det = A11 * A22 - 1.0
    # pseudo-inverse blocks; j = 0 block [[1,-1],[-1,1]] has pinv [[1,-1],[-1,1]]/4
    P11 = np.where(det > 1e-14, A22 / np.where(det > 1e-14, det, 1), 0.25)
    P22 = np.where(det > 1e-14, A11 / np.where(det > 1e-14, det, 1), 0.25)
    P12 = np.where(det > 1e-14, 1.0 / np.where(det > 1e-14, det, 1), -0.25)
    if np.any(det[1:] < 1e-12):
        return None
    Q11 = P11 * P11 + P12 * P12
    Q22 = P22 * P22 + P12 * P12
    Q12 = P11 * P12 + P12 * P22
    blocks1 = {(0, 0): P11, (1, 1): P22, (0, 1): P12, (1, 0): P12}
    blocks2 = {(0, 0): Q11, (1, 1): Q22, (0, 1): Q12, (1, 0): Q12}
    cache = {}

    def G(al, be, dx, which):
        key = (al, be, dx % n, which)
        if key not in cache:
            c = np.cos(th * dx)
            cache[key] = float(np.dot(c, (blocks1 if which == 1 else blocks2)[(al, be)])) / n
        return cache[key]
    trG = n * (G(0, 0, 0, 1) + G(1, 1, 0, 1))
    res = []
    for (vt, nb) in (((0, 0), [(0, 1), (0, -1), (1, 0)]), ((1, 0), [(1, k), (1, -k), (0, 0)])):
        lab = {'v': vt}
        for idx, w in enumerate(nb):
            lab[idx] = w
        gam = lambda p, q: G(lab[p][0], lab[q][0], lab[p][1] - lab[q][1], 1)
        gam2 = lambda p, q: G(lab[p][0], lab[q][0], lab[p][1] - lab[q][1], 2)
        res.append(local_delta(2 * n, list(range(3)), gam, gam2, trG))
    return res


def dense_delta(adj_matrix, v):
    A = adj_matrix
    L = np.diag(A.sum(1)) - A
    n = len(L)
    ev = np.linalg.eigvalsh(L)
    kf = n * np.sum(1 / ev[1:])
    idx = [i for i in range(n) if i != v]
    B = A[np.ix_(idx, idx)]
    LB = np.diag(B.sum(1)) - B
    evb = np.linalg.eigvalsh(LB)
    return kf, (n - 1) * np.sum(1 / evb[1:])


def validate():
    ok = True
    for (a, b, d) in [(3, 0, 4), (4, 1, 5), (2, 1, 6), (5, 2, 5), (1, 3, 11), (6, 0, 6)]:
        for gens in ([(1, 0), (-1, 0), (0, 1), (0, -1)], [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]):
            elems = []
            n = a * d
            # enumerate group elements by reduction: representatives via characters is enough: build Cayley graph directly
            reps = {}
            for x in range(-n, n + 1):
                for y in range(-n, n + 1):
                    key = None
                    # canonical representative of (x,y) mod Lambda: x mod a then y mod d after subtracting
                    p = x // a
                    xr = x - p * a
                    yr = (y - p * b) % d
                    reps.setdefault((xr, yr), (xr, yr))
            elems = sorted(reps)
            if len(elems) != n:
                print('bad enumeration')
                ok = False
                continue
            idx = {e: i for i, e in enumerate(elems)}

            def red(x, y):
                p = x // a
                return (x - p * a, (y - p * b) % d)
            A = np.zeros((n, n))
            simple = True
            for e in elems:
                imgs = [red(e[0] + s[0], e[1] + s[1]) for s in gens]
                if len(set(imgs)) < len(gens) or e in imgs:
                    simple = False
                for t in imgs:
                    A[idx[e], idx[t]] = 1
            if not simple:
                continue
            kf, kfh = dense_delta(A, idx[(0, 0)])
            res = cayley_lattice_delta(a, b, d, gens)
            good = abs(res[0] - kf) < 1e-8 * kf and abs(res[1] - kfh) < 1e-8 * kf
            ok &= good
            print(('PASS' if good else 'FAIL'), 'lattice', (a, b, d), len(gens), kf, res[0], kfh, res[1])
    for (n, k) in [(5, 2), (7, 2), (8, 3), (10, 2), (12, 5), (13, 5)]:
        A = np.zeros((2 * n, 2 * n))
        for i in range(n):
            for (u, w) in [(i, (i + 1) % n), (i, n + i), (n + i, n + (i + k) % n)]:
                A[u, w] = A[w, u] = 1
        r = gp_delta(n, k)
        for vtx, rr in zip((0, n), r):
            kf, kfh = dense_delta(A, vtx)
            good = abs(rr[0] - kf) < 1e-8 * kf and abs(rr[1] - kfh) < 1e-8 * kf
            ok &= good
            print(('PASS' if good else 'FAIL'), 'GP', (n, k), vtx, kf, rr[0], kfh, rr[1])
    print('PASS validation' if ok else 'FAIL validation')
    return ok


def main():
    if sys.argv[1] == 'validate':
        sys.exit(0 if validate() else 1)
    fam, nmin, nmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    thresh = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-9
    near = []
    count = 0
    best = (1e9, None)
    signs = [0, 0]
    for n in range(nmin, nmax + 1):
        if fam in ('tori4', 'tri6'):
            gens = [(1, 0), (-1, 0), (0, 1), (0, -1)] + ([(1, 1), (-1, -1)] if fam == 'tri6' else [])
            bad = [(1, 0), (0, 1), (2, 0), (0, 2), (1, 1), (1, -1)] if fam == 'tori4' else \
                  [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 2), (1, -1), (2, 1), (1, 2)]
            for a in range(1, n + 1):
                if n % a:
                    continue
                d = n // a
                for b in range(d):
                    if any(in_lattice(x, a, b, d) for x in bad):
                        continue
                    r = cayley_lattice_delta(a, b, d, gens)
                    if r is None:
                        continue
                    count += 1
                    rel = (r[0] - r[1]) / r[0]
                    signs[int(rel > 0)] += 1
                    if abs(rel) < best[0]:
                        best = (abs(rel), (n, a, b, d))
                    if abs(rel) < thresh:
                        near.append((n, a, b, d, rel))
                        print('NEAR', n, a, b, d, rel, flush=True)
        elif fam == 'gp':
            for k in range(1, (n - 1) // 2 + 1):
                r = gp_delta(n, k)
                if r is None:
                    continue
                count += 1
                rels = [(x[0] - x[1]) / x[0] for x in r]
                m = max(abs(rels[0]), abs(rels[1]))      # both must vanish for KS
                signs[int(rels[0] > 0)] += 1
                if m < best[0]:
                    best = (m, (n, k))
                if m < thresh:
                    near.append((n, k, rels))
                    print('NEAR', n, k, rels, flush=True)
    print(f'family {fam} n={nmin}..{nmax}: {count} graphs; Delta<0 (at first test vertex): {signs[0]}, Delta>0: {signs[1]}; '
          f'min max|Delta|/Kf = {best[0]:.3e} at {best[1]}; near-misses below {thresh}: {len(near)}')


if __name__ == '__main__':
    main()
