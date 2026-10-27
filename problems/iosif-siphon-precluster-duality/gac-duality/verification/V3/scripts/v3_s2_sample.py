"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: independent exact necessary-condition test of
agent D2's ALL verdicts at s = 2 (round3/dichotomy), with attention to classes of very small / very large mass.

For each network certified ALL by D2 (verdicts recomputed with D2's d2cert.analyse, imported unchanged, only to
obtain the list to test) and for sampled rate vectors kappa and classes k0 * 10^e (e = -6..6 and anisotropic
scalings), every vertex of the class polygon of N* is computed exactly (own code, v3lib), and if it is an equilibrium
its Jacobian on S* is tested for Hurwitz stability (tr < 0, det > 0, exact).  A Hurwitz vertex equilibrium attracts
an open set of the class, i.e. N* would be non-persistent for that kappa: a contradiction with ALL.
The same test is run on D2's NONE networks for information (how often a plain scaling exhibits a stable vertex).

Usage: python v3_s2_sample.py [--nproc P] [--nk K]   Environment: Python 3.12, standard library, exact arithmetic.
"""
import os
import sys
import gzip
import time
import pickle
import random
import platform
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v3lib as L  # noqa: E402

D2DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round3", "dichotomy", "scripts"))


def verdict(q):
    sys.path.insert(0, D2DIR)
    import d2cert  # noqa: E402
    return d2cert.analyse(q["n"], q["rx"])["verdict"]


def vertex_test(n, rx, kap, k0):
    """Return list of Hurwitz vertex equilibria (zero sets) of the class of k0 for N* with rates kap."""
    r, Rd = L.dual(n, rx)
    G = L.Gamma(n, rx)  # rows = dual reaction vectors pi_i (in R^r)
    rows = [[Fraction(x) for x in G[i]] for i in range(n)]
    i1 = next(i for i in range(n) if any(rows[i]))
    i2 = next(i for i in range(n) if L.rank([rows[i1], rows[i]]) == 2)
    B = [[rows[i1][j], rows[i2][j]] for j in range(r)]  # r x 2 basis of S*
    out = []
    seen = set()
    for j1, j2 in combinations(range(r), 2):
        det = B[j1][0] * B[j2][1] - B[j1][1] * B[j2][0]
        if det == 0:
            continue
        # solve k0_j + (B z)_j = 0 for j1, j2
        b1, b2 = -k0[j1], -k0[j2]
        z0 = (b1 * B[j2][1] - B[j1][1] * b2) / det
        z1 = (B[j1][0] * b2 - b1 * B[j2][0]) / det
        k = [k0[j] + B[j][0] * z0 + B[j][1] * z1 for j in range(r)]
        if any(x < 0 for x in k):
            continue
        key = tuple(k)
        if key in seen:
            continue
        seen.add(key)
        # field
        f = [Fraction(0)] * r
        for i, (s, p) in enumerate(Rd):
            rate = kap[i]
            for j in range(r):
                rate *= k[j] ** s[j]
            for j in range(r):
                f[j] += rate * (p[j] - s[j])
        if any(f):
            continue
        # Jacobian Df(k) B, then coordinates in basis B (use rows a, b of B with nonzero minor)
        Df = [[Fraction(0)] * r for _ in range(r)]
        for i, (s, p) in enumerate(Rd):
            for l in range(r):
                if s[l] == 0:
                    continue
                d = kap[i] * s[l]
                for j in range(r):
                    e = s[j] - (1 if j == l else 0)
                    d *= k[j] ** e
                if d == 0:
                    continue
                for j in range(r):
                    Df[j][l] += d * (p[j] - s[j])
        M = [[sum(Df[j][l] * B[l][c] for l in range(r)) for c in range(2)] for j in range(r)]
        a, b = j1, j2
        Bab = [[B[a][0], B[a][1]], [B[b][0], B[b][1]]]
        dd = Bab[0][0] * Bab[1][1] - Bab[0][1] * Bab[1][0]
        inv = [[Bab[1][1] / dd, -Bab[0][1] / dd], [-Bab[1][0] / dd, Bab[0][0] / dd]]
        JS = [[sum(inv[x][y] * M[[a, b][y]][c] for y in range(2)) for c in range(2)] for x in range(2)]
        # consistency: B JS == M
        assert all(sum(B[j][x] * JS[x][c] for x in range(2)) == M[j][c] for j in range(r) for c in range(2))
        tr = JS[0][0] + JS[1][1]
        de = JS[0][0] * JS[1][1] - JS[0][1] * JS[1][0]
        if tr < 0 and de > 0:
            out.append(tuple(j for j in range(r) if k[j] == 0))
    return out


def work(args):
    q, nk, seed = args
    rnd = random.Random(seed)
    n, rx = q["n"], q["rx"]
    r = len(rx)
    hits = []
    for _ in range(nk):
        kap = [Fraction(rnd.randint(1, 30), rnd.randint(1, 30)) for _ in range(n)]
        base = [Fraction(rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(r)]
        scalings = [[Fraction(10) ** e] * r for e in range(-6, 7, 2)]
        for _a in range(4):
            scalings.append([Fraction(10) ** rnd.randint(-5, 5) for _ in range(r)])
        for sc in scalings:
            k0 = [b * s for b, s in zip(base, sc)]
            h = vertex_test(n, rx, kap, k0)
            if h:
                hits.append((kap, k0, h))
                break
    return hits


def main():
    argv = sys.argv[1:]
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 6
    nk = int(argv[argv.index("--nk") + 1]) if "--nk" in argv else 3
    t0 = time.time()
    print("# v3_s2_sample.py  Python %s on %s; exact; nk=%d" % (platform.python_version(), platform.system(), nk))
    recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
    with Pool(nproc) as P:
        ver = P.map(verdict, recs, chunksize=20)
    from collections import Counter
    print("# D2 verdicts recomputed:", dict(Counter(ver)))
    allr = [q for q, v in zip(recs, ver) if v == "ALL"]
    noner = [q for q, v in zip(recs, ver) if v == "NONE"]
    with Pool(nproc) as P:
        resA = P.map(work, [(q, nk, 1000 + i) for i, q in enumerate(allr)], chunksize=20)
    bad = [(q, h) for q, h in zip(allr, resA) if h]
    for q, h in bad[:5]:
        print("   ALL network with a Hurwitz vertex equilibrium:", L.show(q["n"], q["rx"]), h[0])
    ok = not bad
    print(("PASS" if ok else "FAIL") + ": no Hurwitz boundary vertex equilibrium found in %d ALL networks x %d rate "
          "vectors x 11 class scalings (10^-6..10^6, anisotropic)" % (len(allr), nk))
    rnd = random.Random(7)
    sample = rnd.sample(noner, min(3000, len(noner)))
    with Pool(nproc) as P:
        resN = P.map(work, [(q, nk, 5000 + i) for i, q in enumerate(sample)], chunksize=20)
    print("# information: among %d sampled NONE networks, %d show a Hurwitz vertex equilibrium for some sampled "
          "(kappa, scaling) (other NONE mechanisms: non-hyperbolic vertices, edges, theta-scalings)" %
          (len(sample), sum(1 for h in resN if h)))
    print("# runtime %.1f s" % (time.time() - t0))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
