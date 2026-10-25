#!/usr/bin/env python3
"""Floating-point exploration (sanity check only, NOT part of the proof).

Equilibrium-first sampling: for every dynamically nontrivial class with >= 3 distinct sources, draw a point
z* in the open quadrant and a vector m in ker(Gamma) cap R^4_{>0}; the rate constants kappa = m / z*^A make
z* a positive equilibrium (every kappa with a positive equilibrium arises this way).  We then test, in
floating point, whether z* is linearly stable and whether a boundary equilibrium is asymptotically stable at
the same time (origin: hyperbolic test or, for one zero eigenvalue, the sign of the second-order coefficient on
the centre manifold; axis points: both eigenvalues negative).  Expected outcome if the theorem is right:
hits only in the 10 classes of BBH (2024), and in each of them.

Usage: python explore_numeric.py [samples_per_class=20000]      runtime ~16 min on 12 cores
Environment: Python 3.12, numpy 2.x (WSL Fedora 44).  Prints PASS/FAIL; exit 1 on FAIL.
"""
import sys
import numpy as np
from multiprocessing import Pool
from networks import *

MON = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
IDX = {mm: i for i, mm in enumerate(MON)}


def kernel_rays(net):
    import sympy as sp
    from math import gcd
    G = sp.Matrix([[vec(r)[0] for r in net], [vec(r)[1] for r in net]])
    b1, b2 = G.nullspace()
    cands = set()
    for i in range(4):
        a, b = b1[i], b2[i]
        if a == 0 and b == 0:
            continue
        w = b * b1 - a * b2
        w = w * sp.lcm([sp.fraction(x)[1] for x in w])
        wi = [int(x) for x in w]
        g = 0
        for x in wi:
            g = gcd(g, abs(x))
        wi = [x // g for x in wi]
        for s in (wi, [-x for x in wi]):
            if all(x >= 0 for x in s):
                cands.add(tuple(s))
    cands = sorted(cands)
    assert len(cands) == 2
    return np.array(cands[0], float), np.array(cands[1], float)


def coeffs(net, k):
    cf = np.zeros(6); cg = np.zeros(6)
    for kj, r in zip(k, net):
        v = vec(r)
        cf[IDX[r[0]]] += kj * v[0]
        cg[IDX[r[0]]] += kj * v[1]
    return cf, cg


def jac(cf, cg, x, y):
    return np.array([[cf[1] + 2 * cf[3] * x + cf[4] * y, cf[2] + cf[4] * x + 2 * cf[5] * y],
                     [cg[1] + 2 * cg[3] * x + cg[4] * y, cg[2] + cg[4] * x + 2 * cg[5] * y]])


def origin_status(cf, cg):
    if cf[0] != 0 or cg[0] != 0:
        return 'noeq'
    J = jac(cf, cg, 0.0, 0.0)
    sc = max(1e-300, np.max(np.abs(np.concatenate([cf, cg]))))
    ev = np.linalg.eigvals(J).real
    rho, mu = max(ev), min(ev)
    tol = 1e-10 * sc
    if rho < -tol:
        return 'AS'
    if rho > tol:
        return 'unstable'
    if mu > -tol:
        return 'nilpotent'
    w, V = np.linalg.eig(J)
    i = int(np.argmax(w.real))
    e0 = np.real(V[:, i]); e0 = e0 if e0.sum() >= 0 else -e0
    wl, Wl = np.linalg.eig(J.T)
    j = int(np.argmax(wl.real))
    l0 = np.real(Wl[:, j]); l0 = l0 / (l0 @ e0)
    x, y = e0
    q = np.array([cf[3] * x * x + cf[4] * x * y + cf[5] * y * y, cg[3] * x * x + cg[4] * x * y + cg[5] * y * y])
    c = l0 @ q
    if c < -1e-9 * sc:
        return 'AS_cm'
    if c > 1e-9 * sc:
        return 'unstable'
    return 'cm_c0'


def axis_AS(cf, cg):
    out = []
    for axis in (0, 1):
        A, B = (cf, cg) if axis == 0 else (cg, cf)
        iu, ivl = ([0, 1, 3], [2, 4]) if axis == 0 else ([0, 2, 5], [1, 4])
        if any(B[i] != 0 for i in iu):
            continue
        p = [A[iu[0]], A[iu[1]], A[iu[2]]]
        if not any(p[1:]):
            continue
        for r in np.polynomial.polynomial.polyroots(p):
            if abs(r.imag) > 1e-12 or r.real <= 0:
                continue
            u = r.real
            dp = p[1] + 2 * p[2] * u
            lt = B[ivl[0]] + B[ivl[1]] * u
            if dp < 0 and lt < 0:
                out.append(axis)
    return out


def work(args):
    net, nsamp, seed = args
    rng = np.random.default_rng(seed)
    u, v = kernel_rays(net)
    A = np.array([r[0] for r in net], float)
    hits = {}
    for _ in range(nsamp):
        t = rng.uniform(0, 1)
        if rng.random() < 0.3:
            e = 10.0 ** rng.uniform(-6, 0)
            t = e if rng.random() < 0.5 else 1 - e
        m = (1 - t) * u + t * v
        if np.any(m <= 0):
            continue
        z = 10.0 ** rng.uniform(-3, 3, size=2)
        k = m / np.prod(z ** A, axis=1)
        cf, cg = coeffs(net, k)
        J = jac(cf, cg, *z)
        tr, det = np.trace(J), np.linalg.det(J)
        sc = np.max(np.abs(J))
        if not (det > 1e-9 * sc * sc and tr < -1e-9 * sc):
            continue
        os_ = origin_status(cf, cg)
        if os_ in ('AS', 'AS_cm'):
            hits['origin_' + os_] = hits.get('origin_' + os_, 0) + 1
        if axis_AS(cf, cg):
            hits['axis'] = hits.get('axis', 0) + 1
    return net_str(net), hits


def main():
    nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    cl = classes()
    jobs = []
    for i, (key, mem) in enumerate(sorted(cl.items(), key=lambda kv: net_str(kv[1][0]))):
        net = mem[0]
        if cone_canon([vec(r) for r in net]) != ('plane',) or len(key) < 3:
            continue
        jobs.append((net, nsamp, 7000 + i))
    with Pool(12) as pool:
        res = pool.map(work, jobs, chunksize=4)
    orange = set(canon_class(parse_net(s)) for s in ORANGE10)
    hit_classes = set()
    print("classes sampled:", len(res))
    for name, hits in res:
        if hits:
            c = canon_class(parse_net(name))
            hit_classes.add(c)
            print("  hit:", name, hits, "(orange)" if c in orange else "(NOT in the list of 10)")
    ok1 = hit_classes <= orange
    ok2 = hit_classes == orange
    print(("PASS" if ok1 else "FAIL"), "every class with a numerical hit is one of the 10")
    print(("PASS" if ok2 else "FAIL"), "every one of the 10 classes has a numerical hit")
    sys.exit(0 if (ok1 and ok2) else 1)


if __name__ == "__main__":
    main()
