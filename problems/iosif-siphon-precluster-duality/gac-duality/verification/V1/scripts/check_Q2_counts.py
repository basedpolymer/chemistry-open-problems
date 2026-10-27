"""V1 independent recount of rows of Q2's table (round2/preclusters/notes.md section 7) for the two-species families.

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Standard library, exact. For n = 2 every network with
(H1) has s <= 1; with s = 1, Gamma = d gamma^T and A_V = d f(V)^T with f_c(V) = sum_{t in c} gamma_t V_t, so
CF is decided EXACTLY: c ~ c' in CF iff f_{c''} vanishes on every extreme flux ray for c'' not in {c, c'} and f_c is
sign-definite (and not identically zero) on the extreme rays. The dual N* has the same rank s = 1.

Families (as E's enum_all.py / result 17): sets of at most R distinct reactions between distinct nonzero complexes of
molecularity <= m, every species occurring, up to swapping X1, X2, with (H1) and (H2).
"""
import sys, time, platform, random
from itertools import combinations, product
from math import gcd
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
RNG = random.Random(7)
FAILS = []


def cf_rank_one(N):
    """Exact CF for a network with (H2) and rank Gamma = 1 (or 0)."""
    s = srank(N)
    keys, cls = source_classes(N)
    m = len(cls)
    if s == 0:
        return [sorted(c) for c in cls] if m else []
    # gamma: Gamma = d gamma^T
    col = next(j for j in range(N.r) if any(N.G[i][j] for i in range(N.n)))
    i0 = next(i for i in range(N.n) if N.G[i][col])
    gamma = [Fr(N.G[i0][t]) for t in range(N.r)]  # row i0 of Gamma is proportional to gamma
    rays = [v for _, v in flux_rays(N)]
    F = [[sum(gamma[t] * v[t] for t in cls[c]) for v in rays] for c in range(m)]  # f_c on each ray
    zero = [all(x == 0 for x in F[c]) for c in range(m)]
    sd = [not zero[c] and (all(x >= 0 for x in F[c]) or all(x <= 0 for x in F[c])) for c in range(m)]
    edges = []
    for a, b in combinations(range(m), 2):
        if sd[a] and sd[b] and all(zero[c] for c in range(m) if c not in (a, b)):
            edges.append((a, b))
    comps = components(m, edges)
    return canon([[t for c in comp for t in cls[c]] for comp in comps])


def enumerate_n2(mol, R):
    C = [y for y in product(range(mol + 1), repeat=2) if 0 < sum(y) <= mol]
    reacs = [(a, b) for a in C for b in C if a != b]
    bydir = {}
    for a, b in reacs:
        d = (b[0] - a[0], b[1] - a[1])
        if d[0] * d[1] >= 0:
            continue  # (H1) with n = 2 forces every reaction vector to be +-(p,-q), p, q > 0
        g = gcd(abs(d[0]), abs(d[1]))
        prim = (abs(d[0]) // g, abs(d[1]) // g)
        bydir.setdefault(prim, []).append((a, b))
    seen = set()
    out = []
    for prim, rs in bydir.items():
        for k in range(2, R + 1):
            for S in combinations(rs, k):
                signs = {(b[0] - a[0]) > 0 for a, b in S}
                if len(signs) < 2:
                    continue
                key = tuple(sorted(S))
                sw = tuple(sorted(((a[1], a[0]), (b[1], b[0])) for a, b in S))
                can = min(key, sw)
                if can in seen:
                    continue
                seen.add(can)
                Ye = [[S[j][0][i] for j in range(k)] for i in range(2)]
                Yp = [[S[j][1][i] for j in range(k)] for i in range(2)]
                N = Net(Ye, Yp)
                if is_H1(N) and is_H2(N):
                    out.append(N)
    return out


def row(label, nets):
    c = dict(n=0, h3=0, PneqPcs=0, PcsneqPi=0, pineqCF_N=0, pineqCF_D=0, kgenneqCF=0, iso_N=0, iso_D=0, q2=0,
             q2h3=0, sound=0)
    for N in nets:
        D = N.dual()
        c["n"] += 1
        h3 = is_H3(N)
        c["h3"] += h3
        P, Pcs, pi = preclusters_def19(N), class_span_closure(N), flat_closure_partition(N)
        CF, CFD = cf_rank_one(N), cf_rank_one(D)
        piD = flat_closure_partition(D)
        if not (refines(pi, CF) and refines(piD, CFD)):
            c["sound"] += 1
        c["PneqPcs"] += P != Pcs
        c["PcsneqPi"] += Pcs != pi
        c["pineqCF_N"] += pi != CF
        c["pineqCF_D"] += piD != CFD
        V = interior_flux(N, [RNG.randint(1, 10 ** 6) for _ in flux_rays(N)])
        c["kgenneqCF"] += kappa_V(N, V) != CF
        c["iso_N"] += all(any(set(S) <= set(b) for b in CF) for S, _ in flux_rays(N))
        c["iso_D"] += all(any(set(S) <= set(b) for b in CFD) for S, _ in flux_rays(D))
        q2 = all(any(set(s) <= set(b) for b in CFD) for s in min_siphons(N))
        c["q2"] += q2
        c["q2h3"] += q2 and h3
    print("   %s: networks %d ((H3) %d) | P!=Pcs %d | Pcs!=pi %d | pi!=CF N/N* %d/%d | kappa_gen!=CF %d | "
          "isolation N/N* %d/%d | species-side Q2 %d; %d of %d with (H3) | pi not<= CF: %d"
          % (label, c["n"], c["h3"], c["PneqPcs"], c["PcsneqPi"], c["pineqCF_N"], c["pineqCF_D"], c["kgenneqCF"],
             c["iso_N"], c["iso_D"], c["q2"], c["q2h3"], c["h3"], c["sound"]))
    return c


expected = {  # Q2's table, round2/preclusters/notes.md section 7
    (2, 5): dict(n=100, PneqPcs=0, PcsneqPi=10, pineqCF_N=2, pineqCF_D=0, kgenneqCF=0, iso_N=20, iso_D=100, q2=100, q2h3=100),
    (2, 6): dict(n=116, PneqPcs=0, PcsneqPi=10, pineqCF_N=2, pineqCF_D=0, kgenneqCF=0, iso_N=20, iso_D=116, q2=116, q2h3=116),
    (3, 4): dict(n=2784, PneqPcs=27, PcsneqPi=63, pineqCF_N=38, pineqCF_D=0, kgenneqCF=5, iso_N=224, iso_D=2784, q2=2784,
                 q2h3=None),
}
if __name__ == "__main__":
    for (mol, R), exp in expected.items():
        t = time.time()
        nets = enumerate_n2(mol, R)
        c = row("all (2,%d,%d)" % (mol, R), nets)
        diffs = {k: (c[k], v) for k, v in exp.items() if v is not None and c[k] != v}
        ok = not diffs and c["sound"] == 0
        print(("PASS" if ok else "FAIL") + " Q2 table row all (2,%d,%d) reproduced exactly (CF decided exactly, s = 1)%s; %.1f s"
              % (mol, R, "" if not diffs else " -- differences (mine, Q2): %s" % diffs, time.time() - t))
        if not ok:
            FAILS.append((mol, R))
    print("runtime %.1f s" % (time.time() - T0))
    sys.exit(1 if FAILS else 0)
