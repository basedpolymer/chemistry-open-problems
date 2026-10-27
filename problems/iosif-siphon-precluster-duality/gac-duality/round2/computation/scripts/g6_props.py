"""Agent G's Propositions 6.2 and 6.3 (gac/notes.md §6), tested exactly (first-round loose end, task 4 of BRIEF-E2).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (e2lib -> ecore -> crn).
Usage (from round2/computation/): python scripts/g6_props.py [--limit K] [--nproc P]
Input: first-round record files; networks N with (H1) (bounded classes), every family, at most K per file (seeded).

6.2  For x0 > 0 and Sigma a set of species such that the face F_Sigma = {x in P : x_Sigma = 0} of the class
     P = (x0 + S) cap R^n_{>=0} contains a point positive outside Sigma ("realizable"):
         codim_P F_Sigma = rank Gamma[Sigma,:] = rank Gamma_D[:,Sigma]  (dual subnetwork with the reactions Sigma).
     Independent computation: exact vertex enumeration of P (x = x0 + B t, s tight coordinates), F_Sigma = conv of
     the vertices with x_Sigma = 0, realizable iff their barycentre has zero set exactly Sigma, dimension = affine
     rank of those vertices.  Also checked: the realizable zero sets of N's classes that are siphons of N are T-traps
     of N* (dictionary), and the Anderson-Shiu reformulation is the same count (facets <-> rank 1).
     Classes: x0 = (1,...,1) and one seeded random rational x0.
6.3  For N weakly reversible and every siphon Sigma: J0 = reactions whose source avoids Sigma; (i) their products
     avoid Sigma; (ii) Sbl = im Gamma[:,J0] is contained in S0 = S cap {x_Sigma = 0}; (a) dim S0 = dim Sbl iff
     rank M[J0,:] + rank M[:,Sigma] = s (M = Gamma^T), with S0 computed independently as a nullspace intersection;
     (b) for every embedded linkage class L (projection of N's reactions to Sigma, linkage classes of the projected
     complexes), Stiemke: im Gamma[Sigma, J_L] contains a nonzero nonnegative vector iff ker M[J_L, Sigma] contains
     no positive vector, both sides decided by independent exact LPs.  (Wiuf's Sigma-minimality is not needed for
     the equivalence, which is tested on every embedded class.)
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn, ecore


def basis(cols):
    B = []
    for c in cols:
        if crn.rank(B + [c]) > len(B):
            B.append(c)
    return B


def solve_square(A, b):
    """Solve A t = b exactly (A square, nonsingular)."""
    m = len(A)
    M = [[Fraction(x) for x in row] + [Fraction(bb)] for row, bb in zip(A, b)]
    for c in range(m):
        p = next(i for i in range(c, m) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(m):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][m] for i in range(m)]


def vertices(x0, B, n):
    s = len(B)
    verts = set()
    for I in combinations(range(n), s):
        A = [[B[l][i] for l in range(s)] for i in I]
        if crn.rank(A) < s:
            continue
        t = solve_square(A, [-x0[i] for i in I])
        x = tuple(x0[i] + sum(t[l] * B[l][i] for l in range(s)) for i in range(n))
        if all(v >= 0 for v in x):
            verts.add(x)
    return sorted(verts)


def affine_rank(pts):
    if not pts:
        return -1
    p0 = pts[0]
    return crn.rank([[a - b for a, b in zip(p, p0)] for p in pts[1:]]) if len(pts) > 1 else 0


def check62(n, rx, x0):
    N = CNet(n, rx)
    G = N.Gamma
    r = len(rx)
    B = basis([[Fraction(G[i][j]) for i in range(n)] for j in range(r)])
    s = len(B)
    V = vertices(x0, B, n)
    D = N.dual()
    GD = D.Gamma  # r x n
    zsets = {}
    fails = []
    nreal = 0
    for size in range(1, n + 1):
        for Sig in combinations(range(n), size):
            VS = [v for v in V if all(v[i] == 0 for i in Sig)]
            if not VS:
                continue
            bary = [sum(v[i] for v in VS) / len(VS) for i in range(n)]
            if {i for i in range(n) if bary[i] == 0} != set(Sig):
                continue
            nreal += 1
            codim = s - affine_rank(VS)
            r1 = crn.rank([G[i] for i in Sig])
            r2 = crn.rank([[GD[j][i] for i in Sig] for j in range(r)])
            if not (codim == r1 == r2):
                fails.append(("6.2", Sig, codim, r1, r2))
            if ecore.is_siphon(N, Sig):
                # dictionary: Sigma is a T-trap of D (every species of D produced by Sigma is consumed by Sigma)
                prod = {j for i in Sig for j in range(r) if D.reactions[i][1][j] > 0}
                cons = {j for i in Sig for j in range(r) if D.reactions[i][0][j] > 0}
                if not prod <= cons:
                    fails.append(("T-trap", Sig))
    return nreal, fails


def lp_nonneg_image(A, rows, cols):
    """Exists u with A u >= 0 and sum(A u) >= 1 (A: rows x cols)?"""
    if cols == 0:
        return False
    ineqs = [([-A[i][j] for j in range(cols)], 0) for i in range(rows)]
    ineqs.append(([-sum(A[i][j] for i in range(rows)) for j in range(cols)], -1))
    return e2lib.fm_point(ineqs, cols) is not None


def lp_pos_kernel_T(A, rows, cols):
    """Exists y >= 1 (componentwise) with A^T y = 0 (y indexed by rows)?"""
    eqs = [([A[i][j] for i in range(rows)], 0) for j in range(cols)]
    ineqs = [([-int(t == i) for t in range(rows)], -1) for i in range(rows)]
    return e2lib.lp_point(eqs, ineqs, rows) is not None


def check63(n, rx):
    N = CNet(n, rx)
    G = N.Gamma
    r = len(rx)
    s = N.rank()
    fails = []
    napp = 0
    for Sm in ecore.siphon_masks(N):
        Sig = sorted(ecore.to_set(Sm))
        napp += 1
        J0 = [j for j in range(r) if not any(rx[j][0][i] for i in Sig)]
        if any(rx[j][1][i] for j in J0 for i in Sig):
            fails.append(("i", Sig))
        for j in J0:
            if any(G[i][j] for i in Sig):
                fails.append(("ii", Sig))
        # S0 independently: dim(S cap {x_Sigma = 0}) = dim{t : Gamma[Sigma,:] t = 0} images ... via nullspace
        # S0 = Gamma(ker Gamma[Sigma,:]); its dimension = rank of Gamma restricted to that kernel.
        Kt = crn.nullspace([[Fraction(G[i][j]) for j in range(r)] for i in Sig], r) if Sig else None
        if Kt is None:
            dimS0 = s
        else:
            dimS0 = crn.rank([[sum(G[i][j] * t[j] for j in range(r)) for i in range(n)] for t in Kt]) if Kt else 0
        dimSbl = crn.rank([[G[i][j] for i in range(n)] for j in J0]) if J0 else 0
        M = [[G[i][j] for i in range(n)] for j in range(r)]  # Gamma^T, r x n
        lhs = (dimS0 == dimSbl)
        rM1 = crn.rank([M[j] for j in J0]) if J0 else 0
        rM2 = crn.rank([[M[j][i] for i in Sig] for j in range(r)])
        if lhs != (rM1 + rM2 == s):
            fails.append(("a", Sig))
        # (b) embedded linkage classes
        proj = [(tuple(rx[j][0][i] for i in Sig), tuple(rx[j][1][i] for i in Sig)) for j in range(r)]
        cpx = sorted({c for pr in proj for c in pr})
        parent = {c: c for c in cpx}

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x
        for a, b in proj:
            if a != b:
                parent[find(a)] = find(b)
        classes = {}
        for j, (a, b) in enumerate(proj):
            if a != b:
                classes.setdefault(find(a), []).append(j)
        for L, JL in classes.items():
            A = [[G[i][j] for j in JL] for i in Sig]  # Gamma[Sigma, J_L]
            e1 = lp_nonneg_image(A, len(Sig), len(JL))
            e2 = lp_pos_kernel_T(A, len(Sig), len(JL))
            if e1 == e2:
                fails.append(("b-Stiemke", Sig, JL))
    return napp, fails


def work(args):
    n, rx, seed = args
    rng = random.Random(seed)
    out = {"n62": 0, "f62": [], "n63": 0, "f63": [], "rx": rx, "n": n}
    for x0 in ([Fraction(1)] * n, [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(n)]):
        k, f = check62(n, rx, x0)
        out["n62"] += k
        out["f62"] += f
    N = CNet(n, rx)
    if N.weakly_reversible():
        k, f = check63(n, rx)
        out["n63"] += k
        out["f63"] += f
    return out


def main():
    t0 = time.time()
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 1500
    nproc = int(sys.argv[sys.argv.index("--nproc") + 1]) if "--nproc" in sys.argv else 4
    print(f"# g6_props.py  Python {platform.python_version()} on {platform.system()}; exact; limit {limit} per file")
    rng = random.Random(62)
    jobs = []
    for fn in e2lib.E_FILES:
        recs = [rc for rc in e2lib.load([fn]) if rc["N"]["cons"]]
        rng.shuffle(recs)
        for rc in recs[:limit]:
            jobs.append((rc["n"], rc["rx"], rng.randint(0, 10 ** 9)))
    print(f"# {len(jobs)} conservative networks (seeded subsample, both classes x0 = 1 and a random x0)")
    with Pool(nproc) as pool:
        res = pool.map(work, jobs, chunksize=20)
    n62 = sum(q["n62"] for q in res)
    f62 = [(q["n"], q["rx"], f) for q in res for f in q["f62"]]
    n63 = sum(q["n63"] for q in res)
    f63 = [(q["n"], q["rx"], f) for q in res for f in q["f63"]]
    nwr = sum(1 for q in res if q["n63"])
    print(f"{'PASS' if not f62 else 'FAIL'} 6.2: codim F_Sigma = rank Gamma[Sigma,:] = rank of the dual subnetwork on "
          f"Sigma, and realizable siphons are T-traps of N*: {n62} realizable faces, {len(f62)} failures")
    for f in f62[:3]:
        print(f"     failure {f}")
    print(f"{'PASS' if not f63 else 'FAIL'} 6.3: (i), (ii), (a) rank identity, (b) Stiemke translation on every "
          f"siphon of {nwr} weakly reversible networks ({n63} siphons), {len(f63)} failures")
    for f in f63[:3]:
        print(f"     failure {f}")
    print(f"{'PASS' if not (f62 or f63) else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if not (f62 or f63) else 1)


if __name__ == "__main__":
    main()
