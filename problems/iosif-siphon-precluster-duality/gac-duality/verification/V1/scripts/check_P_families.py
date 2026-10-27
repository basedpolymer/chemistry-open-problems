"""V1 independent re-enumeration of two of agent P's EC families (petri/logs/enum_ec_*.log) and re-check of
Theorems A, B, C, C', D and Prop F on them.

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Standard library, exact; independent code (v1lib).
Family definition (from P's enumerate_ec.py docstring, re-implemented): networks up to species permutation whose
distinct source complexes have pairwise disjoint supports covering all species (coefficients >= 1, molecularity
<= m; ordinary: coefficients 1), each source with a nonempty set of distinct products (nonzero complexes of
molecularity <= m different from the source; ordinary: 0/1 complexes), at most R reactions, connected, (H1), (H2).
Usage: python check_P_families.py n m R ordinary(0/1)
"""
import sys, time, platform
from itertools import combinations, permutations, product
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
from check_P import thmA, thmB_species, thmC_Cp, thmD, propF

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
n, m, R, ordinary = map(int, sys.argv[1:5])


def set_partitions(xs):
    if not xs:
        yield []
        return
    first, rest = xs[0], xs[1:]
    for p in set_partitions(rest):
        for k in range(len(p)):
            yield p[:k] + [[first] + p[k]] + p[k + 1:]
        yield [[first]] + p


coef = [1] if ordinary else list(range(1, m + 1))
cpx = [c for c in product(range(2) if ordinary else range(m + 1), repeat=n) if 0 < sum(c) <= m]
perms = list(permutations(range(n)))


def canonical(rx):
    best = None
    for p in perms:
        t = tuple(sorted((tuple(a[p[i]] for i in range(n)), tuple(b[p[i]] for i in range(n))) for a, b in rx))
        if best is None or t < best:
            best = t
    return best


seen = set()
nets = []
for part in set_partitions(list(range(n))):
    src_choices = []
    for block in part:
        opts = []
        for cs in product(coef, repeat=len(block)):
            y = [0] * n
            for i, c in zip(block, cs):
                y[i] = c
            if sum(y) <= m:
                opts.append(tuple(y))
        src_choices.append(opts)
    for sources in product(*src_choices):
        prods = [[c for c in cpx if c != s] for s in sources]

        def rec(k, budget, acc):
            if k == len(sources):
                yield list(acc)
                return
            rest = len(sources) - k - 1
            for size in range(1, budget - rest + 1):
                for P in combinations(prods[k], size):
                    yield from rec(k + 1, budget - size, acc + [(sources[k], p) for p in P])
        for rx in rec(0, R, []):
            can = canonical(rx)
            if can in seen:
                continue
            seen.add(can)
            r = len(rx)
            N = Net([[rx[j][0][i] for j in range(r)] for i in range(n)], [[rx[j][1][i] for j in range(r)] for i in range(n)])
            if is_connected(N) and is_H1(N) and is_H2(N):
                nets.append(N)

cnt = dict(EC=0, H3=0, WF=0, WFFC=0, WR=0, pureWF=0)
bad = {}
for N in nets:
    assert is_EC(N)
    cnt["EC"] += 1
    h3, wf, wr = is_H3(N), srank(N) == n_sources(N) - 1, is_WR(N)
    cnt["H3"] += h3; cnt["WF"] += wf; cnt["WR"] += wr
    D = N.dual()
    for side in (N, D):
        if not thmA(side):
            bad.setdefault("Theorem A", N)
        if not thmB_species(side):
            bad.setdefault("Theorem B", N)
    if wf:
        noC, noCp, ms, mt = thmC_Cp(N)
        if is_pure(N):
            cnt["pureWF"] += 1
            if not (noC and ms):
                bad.setdefault("Theorem C (pure)", N)
        if not (noCp and mt):
            bad.setdefault("Theorem C'", N)
        if is_ordinary(N):
            cnt["WFFC"] += 1
            if not thmD(N):
                bad.setdefault("Theorem D", N)
    if wr:
        if not propF(N):
            bad.setdefault("Prop F", N)
print("family n=%d m<=%d r<=%d ordinary=%d: %s" % (n, m, R, ordinary, cnt))
for k in ("Theorem A", "Theorem B", "Theorem C (pure)", "Theorem C'", "Theorem D", "Prop F"):
    print(("PASS " if k not in bad else "FAIL ") + k + ("" if k not in bad else " -- %r" % bad[k]))
print("runtime %.1f s" % (time.time() - T0))
sys.exit(1 if bad else 0)
