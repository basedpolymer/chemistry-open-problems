"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: seeded random tests (corroboration of proved
statements, not proofs), independent of the authors' code.

  RW1  Theorem S1-WR (Q1): WR rank-one N  =>  no DOM face of N* (random, molecularity up to 8, n = 2..4).
  RW2  Theorem S1-WR under the weaker hypothesis (the two extreme complexes of each linkage class are sources),
       random non-WR rank-one networks.
  RW3  Theorem W (D2): WR N of rank 2 and 3  =>  no facet set Z of N* (siphon or not) is DOM.
  RW4  Prop. 2.1 (G): q = beta* + delta* on random networks (any rank, self-loops excluded).
  RW5  Prop. 2.2 (G): N WR with (H3)  =>  [N* WR  <=>  source-row set = product-row set], and then
       m* = n, l* = c(pi), delta* = q - c(pi), each pi-orbit is a conservation law.
  RW6  Prop. WR-Sigma (D2): k_j' = G_{y'_j} - G_{y_j} for random WR N (sympy identity).
Usage: python v3_random.py [seed]   (exit 1 on FAIL).  Environment: Python 3.12, sympy; exact arithmetic.
"""
import sys
import time
import platform
import random
from fractions import Fraction
from itertools import permutations
from math import gcd
import sympy as sp

import v3lib as L

FAILS = 0


def claim(name, ok, msg=""):
    global FAILS
    print(("PASS" if ok else "FAIL") + ": " + name + ((" -- " + msg) if msg else ""))
    FAILS += 0 if ok else 1


def rand_sc_digraph(rnd, k, extra):
    """Random strongly connected digraph on 0..k-1: random Hamiltonian cycle plus random extra arcs."""
    order = list(range(k))
    rnd.shuffle(order)
    E = {(order[i], order[(i + 1) % k]) for i in range(k)}
    for _ in range(extra):
        u, v = rnd.sample(range(k), 2)
        E.add((u, v))
    return sorted(E)


def rand_rank1_wr(rnd, n, maxmol, weak=False):
    while True:
        a = [rnd.randint(-3, 3) for _ in range(n)]
        if not (any(x > 0 for x in a) and any(x < 0 for x in a)):
            continue
        g = 0
        for x in a:
            g = gcd(g, abs(x))
        a = [x // g for x in a]
        break
    R = []
    used = set()
    for _ in range(rnd.randint(1, 3)):
        for _try in range(50):
            y0 = [rnd.randint(0, maxmol) for _ in range(n)]
            ts = [t for t in range(-maxmol, maxmol + 1)
                  if all(y0[i] + t * a[i] >= 0 for i in range(n)) and 0 < sum(y0[i] + t * a[i] for i in range(n)) <= maxmol]
            if len(ts) >= 2:
                break
        else:
            continue
        k = rnd.randint(2, min(5, len(ts)))
        pos = sorted(rnd.sample(ts, k))
        pts = [tuple(y0[i] + t * a[i] for i in range(n)) for t in pos]
        if any(p in used for p in pts):
            continue
        used.update(pts)
        if not weak:
            E = rand_sc_digraph(rnd, k, rnd.randint(0, 3))
        else:
            # extreme positions 0 and k-1 must be sources; others arbitrary; random arcs, not necessarily WR
            E = set()
            E.add((0, rnd.randint(1, k - 1)))
            E.add((k - 1, rnd.randint(0, k - 2)))
            for _ in range(rnd.randint(0, 4)):
                u, v = rnd.sample(range(k), 2)
                E.add((u, v))
            E = sorted(E)
        R.extend((pts[u], pts[v]) for u, v in E)
    return R


def rand_wr(rnd, n, maxmol, ncl, kmax):
    R = []
    used = set()
    for _ in range(ncl):
        k = rnd.randint(2, kmax)
        pts = []
        tries = 0
        while len(pts) < k and tries < 500:  # bounded: the pool of complexes may be exhausted
            tries += 1
            y = tuple(rnd.randint(0, maxmol) for _ in range(n))
            if 0 < sum(y) <= maxmol and y not in used and y not in pts:
                pts.append(y)
        if len(pts) < 2:
            break
        k = len(pts)
        used.update(pts)
        E = rand_sc_digraph(rnd, k, rnd.randint(0, 2))
        R.extend((pts[u], pts[v]) for u, v in E)
    return R


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 271027
    rnd = random.Random(seed)
    t0 = time.time()
    print("# v3_random.py seed=%d  Python %s on %s; sympy %s; exact" % (seed, platform.python_version(),
                                                                       platform.system(), sp.__version__))
    # RW1
    tested = dom = 0
    verdicts = {'ALL': 0, 'NONE': 0, 'DOM-only': 0}
    for _ in range(4000):
        n = rnd.randint(2, 4)
        R = rand_rank1_wr(rnd, n, rnd.randint(2, 8))
        if not R or len(R) > 10 or L.rank(L.Gamma(n, R)) != 1 or not L.is_wr(R):
            continue
        nd, Rd = L.dual(n, R)
        v, ft = L.verdict_rank_one(nd, Rd)
        verdicts[v] += 1
        tested += 1
        if any(t == 'DOM' for _, _, t in ft):
            dom += 1
            print("   DOM witness:", L.show(n, R))
    claim("RW1 S1-WR: no DOM face in the dual of %d random WR rank-one networks (verdicts %s)" % (tested, verdicts),
          dom == 0 and tested > 1000)
    # RW2
    tested = dom = 0
    for _ in range(4000):
        n = rnd.randint(2, 4)
        R = rand_rank1_wr(rnd, n, rnd.randint(2, 8), weak=True)
        if not R or len(R) > 10 or L.rank(L.Gamma(n, R)) != 1:
            continue
        nd, Rd = L.dual(n, R)
        ft = L.face_types(nd, Rd)
        tested += 1
        if any(t == 'DOM' for _, _, t in ft):
            dom += 1
            print("   DOM witness (weak hypothesis):", L.show(n, R))
    claim("RW2 S1-WR weak hypothesis (extreme complexes are sources): no DOM face in %d random duals" % tested,
          dom == 0 and tested > 1000)
    # RW3
    stats = {2: [0, 0], 3: [0, 0]}
    dom = 0
    for _ in range(3000):
        n = rnd.randint(3, 5)
        R = rand_wr(rnd, n, rnd.randint(2, 4), rnd.randint(1, 3), 4)
        s = L.rank(L.Gamma(n, R))
        if s not in (2, 3) or len(R) > 12:
            continue
        ft = L.facet_types_any_rank(n, R)
        stats[s][0] += 1
        stats[s][1] += len(ft)
        if any(t == 'DOM' for _, t in ft):
            dom += 1
            print("   DOM facet witness:", L.show(n, R))
    claim("RW3 Theorem W: no DOM facet set (siphon or not) in the duals of random WR networks: s=2: %d networks / "
          "%d facet sets, s=3: %d / %d" % (stats[2][0], stats[2][1], stats[3][0], stats[3][1]),
          dom == 0 and stats[2][0] > 300 and stats[3][0] > 300)
    # RW4
    ok = True
    cnt = 0
    for _ in range(2000):
        n = rnd.randint(1, 5)
        r = rnd.randint(1, 6)
        R = []
        for _j in range(r):
            s_ = tuple(rnd.randint(0, 2) for _ in range(n))
            p_ = tuple(rnd.randint(0, 2) for _ in range(n))
            if s_ != p_:
                R.append((s_, p_))
        if not R:
            continue
        nd, Rd = L.dual(n, R)
        md, ld, sd, dd = L.counts(nd, Rd)
        beta = n - (md - ld)  # dim ker of the incidence matrix of the dual graph (n dual reactions, loops allowed)
        q = n - L.rank(L.Gamma(n, R))
        cnt += 1
        if q != beta + dd or sd != L.rank(L.Gamma(n, R)):
            ok = False
            print("   RW4 failure", L.show(n, R))
    claim("RW4 q = beta* + delta* on %d random networks" % cnt, ok)
    # RW5
    ok = True
    nwr = nwrd = 0
    for it in range(3000):
        n = rnd.randint(2, 5)
        if it % 2 == 0:
            R = rand_wr(rnd, n, rnd.randint(1, 3), rnd.randint(1, 3), 4)
        else:
            # coordinate-permutation construction: complexes y, Py, P^2 y, ... (dual WR by Prop. 2.2)
            P = list(range(n))
            rnd.shuffle(P)
            R = []
            used = set()
            for _ in range(rnd.randint(1, 3)):
                y = tuple(rnd.randint(0, 2) for _ in range(n))
                orbit = [y]
                while True:
                    z = tuple(orbit[-1][P[i]] for i in range(n))
                    if z == y:
                        break
                    orbit.append(z)
                if len(orbit) < 2 or any(o in used for o in orbit) or sum(y) == 0:
                    continue
                used.update(orbit)
                R.extend((orbit[i], orbit[(i + 1) % len(orbit)]) for i in range(len(orbit)))
        if not R or not L.is_wr(R) or not L.H3(n, R):
            continue
        nwr += 1
        nd, Rd = L.dual(n, R)
        Yer = {tuple(row) for row in L.Ye(n, R)}
        Ypr = {tuple(row) for row in L.Yp(n, R)}
        dwr = L.is_wr(Rd)
        if dwr != (Yer == Ypr):
            ok = False
            print("   RW5 equivalence failure", L.show(n, R))
            continue
        if dwr:
            nwrd += 1
            ye = [tuple(row) for row in L.Ye(n, R)]
            yp = [tuple(row) for row in L.Yp(n, R)]
            pi = [ye.index(yp[i]) for i in range(n)]
            seen, cyc = set(), 0
            orbits = []
            for i in range(n):
                if i not in seen:
                    cyc += 1
                    o = []
                    j = i
                    while j not in seen:
                        seen.add(j); o.append(j); j = pi[j]
                    orbits.append(o)
            md, ld, sd, dd = L.counts(nd, Rd)
            q = n - L.rank(L.Gamma(n, R))
            G = L.Gamma(n, R)
            laws = all(all(sum(G[i][j] for i in o) == 0 for j in range(len(R))) for o in orbits)
            if not (md == n and ld == cyc and dd == q - cyc and laws):
                ok = False
                print("   RW5 count failure", L.show(n, R))
    claim("RW5 Prop 2.2 on %d WR networks with (H3) (%d with WR dual)" % (nwr, nwrd), ok and nwrd > 50)
    # RW6
    ok = True
    for _ in range(40):
        n = rnd.randint(2, 4)
        R = rand_wr(rnd, n, 3, rnd.randint(1, 2), 3)
        r = len(R)
        ks = sp.symbols('k0:%d' % r, positive=True)
        kap = sp.symbols('c0:%d' % n, positive=True)
        nd, Rd = L.dual(n, R)
        f = [0] * r
        for i, (s_, p_) in enumerate(Rd):
            rate = kap[i]
            for j in range(r):
                rate *= ks[j] ** s_[j]
            for j in range(r):
                f[j] += rate * (p_[j] - s_[j])
        Cx = L.complexes(R)
        Pc = {c: sp.Mul(*[ks[j] for j in range(r) if R[j][0] == c]) for c in Cx}
        Psym = {c: sp.Symbol('P%d' % idx, positive=True) for idx, c in enumerate(Cx)}
        Psi = sum(kap[i] * sp.Mul(*[Psym[c] ** c[i] for c in Cx if R and any(R[j][0] == c for j in range(r))])
                  for i in range(n))

        def Gc(c):
            return (Psym[c] * sp.diff(Psi, Psym[c])).subs({Psym[d]: Pc[d] for d in Cx})
        for j in range(r):
            if sp.expand(f[j] - (Gc(R[j][1]) - Gc(R[j][0]))) != 0:
                ok = False
    claim("RW6 WR-Sigma: k_j' = G_{y'_j}(P) - G_{y_j}(P) on 40 random WR networks", ok)
    print("# runtime %.1f s" % (time.time() - t0))
    print("ALL CLAIMS PASS" if FAILS == 0 else "%d FAIL" % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
