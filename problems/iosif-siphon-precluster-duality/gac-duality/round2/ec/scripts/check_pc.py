"""Theorem PC (round2/ec/notes.md §8): proportional-conflict networks (two reactions sharing a reactant species have
positively proportional source complexes; equivalently Ye is block diagonal with rank-one positive blocks, one per
Petri cluster - a self-dual class containing EC and its dual class BH).  Agent Q3, Claude Opus 5.5 (effort xhigh).

Environment: Python 3.12, standard library only, exact rationals (crn.py).  Usage (from round2/ec/):
python scripts/check_pc.py NSPECIES NNETS SEED.  Exit 1 on FAIL.  Runtime printed.
Random PC networks (seeded): species partitioned into blocks; block B gets a primitive shape vector b_B (entries
1..2 on B) and its reactions have sources a_t * b_B (a_t in 1..2); products random nonzero complexes of
molecularity <= 4.  Kept if connected, (H1), (H2), s = #clusters - 1 and the network is PC but neither EC nor block-homogeneous
('genuine'; at most 40*NNETS draws); all kept networks are tested:
  PC1  network and its dual are PC with the same hypotheses (self-duality);
  PC2  at two random rational points x and random rational k: J|_S Hurwitz (exact charpoly, Routh array) and the
       matrix G = Q P is Metzler with a positive right null vector u(a o w) (exact kernel of A_{a o w});
  PC3  the same for the dual N* (at random points and kappa);
  PC4  no critical siphon in N and N* (Theorem C^gen).
Also the fixed example NP = {A+2B -> C, 2A+4B -> C+D, C -> A+2B, D -> A+2B}.
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import product

import eclib
from eclib import Network


def clusters(N):
    return [c for c in eclib.petri.clusters(N) if c[1]]


def is_pc(N):
    for S, T in clusters(N):
        ys = [N.reactions[t][0] for t in T]
        y0 = ys[0]
        for y in ys[1:]:
            # positively proportional: y * y0[i] == y0 * y[i] componentwise for the pivot i
            i = next(j for j in range(N.n) if y0[j])
            if any(y[j] * y0[i] != y0[j] * y[i] for j in range(N.n)) or y[i] == 0:
                return False
    return True


def wf_pc(N):
    return (is_pc(N) and eclib.petri.is_connected(N) and N.rank() == len(clusters(N)) - 1
            and N.is_consistent() and N.is_conservative())


def jac_certificate(N, rng):
    cl = clusters(N)
    q = len(cl)
    x = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.n)]
    k = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.r)]
    w = [k[t] * eclib.monomial(x, N.reactions[t][0]) for t in range(N.r)]
    # shape vectors b_c (primitive) and multiples a_t
    a = [None] * N.r
    b = {}
    for ci, (S, T) in enumerate(cl):
        T = sorted(T)
        y0 = N.reactions[T[0]][0]
        from math import gcd
        g = 0
        for v in y0:
            g = gcd(g, v)
        bc = [v // g for v in y0]
        b[ci] = bc
        i = next(j for j in range(N.n) if bc[j])
        for t in T:
            a[t] = Fraction(N.reactions[t][0][i], bc[i])
    G = N.Gamma
    rcl = {}
    for ci, (S, T) in enumerate(cl):
        for t in T:
            rcl[t] = ci
    P = [[sum((G[p][t] * a[t] * w[t] for t in range(N.r) if rcl[t] == c), Fraction(0)) for c in range(q)]
         for p in range(N.n)]
    Q = [[Fraction(b[c][p]) / x[p] for p in range(N.n)] for c in range(q)]
    Gm = eclib.matmul(Q, P)
    if any(Gm[c][d] < 0 for c in range(q) for d in range(q) if c != d):
        return False
    ker = eclib.nullspace(P, q)
    if len(ker) != 1 or not (all(v > 0 for v in ker[0]) or all(v < 0 for v in ker[0])):
        return False
    J = eclib.jacobian(N, k, x)
    cp = eclib.charpoly(J)
    low = N.n - N.rank()
    if any(cp[i] != 0 for i in range(low)):
        return False
    return cp[low] != 0 and eclib.is_hurwitz(cp[low:])


def random_pc(n, rng):
    cpx = [c for c in product(range(5), repeat=n) if 0 < sum(c) <= 4]
    for _ in range(100000):
        # random block partition into 2..n blocks
        nb = rng.randint(2, n)
        lab = [rng.randrange(nb) for _ in range(n)]
        blocks = [[i for i in range(n) if lab[i] == c] for c in range(nb)]
        blocks = [B for B in blocks if B]
        rx = set()
        for B in blocks:
            bvec = [0] * n
            for i in B:
                bvec[i] = rng.randint(1, 2)
            if sum(bvec) > 4:
                bvec = [1 if v else 0 for v in bvec]
            for _ in range(rng.randint(1, 2)):
                amult = rng.randint(1, 2)
                if amult * sum(bvec) > 4:
                    amult = 1
                y = tuple(amult * v for v in bvec)
                yp = cpx[rng.randrange(len(cpx))]
                if yp != y:
                    rx.add((y, yp))
        N = Network(n, sorted(rx))
        if N.r >= 2 and wf_pc(N):
            return N
    return None


def main():
    n, nnets, seed = (int(a) for a in sys.argv[1:4])
    t0 = time.time()
    print(f"# check_pc.py n={n} nets={nnets} seed={seed}; Python {platform.python_version()} on {platform.system()}")
    rng = random.Random(seed)
    nets = [Network(4, [((1, 2, 0, 0), (0, 0, 1, 0)), ((2, 4, 0, 0), (0, 0, 1, 1)), ((0, 0, 1, 0), (1, 2, 0, 0)),
                        ((0, 0, 0, 1), (1, 2, 0, 0))], names=["A", "B", "C", "D"])]
    tries = 0
    while len(nets) < nnets + 1 and tries < 40 * nnets:
        tries += 1
        N = random_pc(n, rng)
        if N is not None and not eclib.ECStructure(N).is_ec() and not eclib.is_balanced(N):
            nets.append(N)            # keep only 'genuine' PC networks: neither EC nor block-homogeneous
    tally = {}
    bad = 0
    ngen = 0
    for N in nets:
        D = N.dual()
        genuine = not eclib.ECStructure(N).is_ec() and not eclib.is_balanced(N)
        ngen += genuine
        res = {
            "PC1 self-duality": wf_pc(D),
            "PC2 J|_S Hurwitz, G Metzler with positive null vector (N)": all(jac_certificate(N, rng) for _ in range(2)),
            "PC3 same for N*": all(jac_certificate(D, rng) for _ in range(2)),
            "PC4 no critical siphon in N and N*": not eclib.petri.critical_siphons(N) and not eclib.petri.critical_siphons(D),
        }
        for kname, v in res.items():
            h, f = tally.get(kname, (0, 0))
            tally[kname] = (h + v, f + (not v))
            if not v:
                print(f"FAIL  {kname}: {N.describe()}")
    print(f"# networks tested: {len(nets)} (neither EC nor balanced: {ngen}); first: {nets[0].describe()}")
    for kname, (h, f) in sorted(tally.items()):
        bad += f > 0
        print(f"{'PASS' if f == 0 else 'FAIL'}  {kname}: {h}/{h + f}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
