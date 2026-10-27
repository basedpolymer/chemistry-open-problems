"""K8 strong form, species side: random search (agent Q3, Claude Opus 5.5, effort xhigh; candidates Q3-K8s).

Looks for two networks N1, N2 on the same species set with (H1)-(H3), rows and columns of Gamma pairwise linearly
independent (both), the same twin classes and the same minimal siphons (up to a common relabelling of the species)
and different P(N*).  Networks are drawn at random (seeded) with a forced twin pair {X1, X2}: every source complex
has equal coefficients on X1 and X2.  Exact arithmetic (crn.py).  Environment: Python 3.12, standard library only.
Usage: python k8s_search.py N MAXMOL RMIN RMAX NSAMPLES SEED [NPROC].  Exit 0 (search).  Runtime printed.
"""
import sys
import time
import random
import platform
from itertools import permutations, product
from multiprocessing import Pool

import eclib
from eclib import Network


def twin_partition(N):
    Ye = N.Ye
    cl = {}
    for i in range(N.n):
        cl.setdefault(tuple(Ye[i]) if any(Ye[i]) else ("zero", i), []).append(i)
    return frozenset(frozenset(b) for b in cl.values())


def draw(n, maxmol, r, rng):
    cpx = [c for c in product(range(maxmol + 1), repeat=n) if 0 < sum(c) <= maxmol]
    srcs = [c for c in cpx if c[0] == c[1]]
    rx = set()
    while len(rx) < r:
        a = rng.choice(srcs)
        b = rng.choice(cpx)
        if a != b:
            rx.add((a, b))
    return tuple(sorted(rx))


def worker(args):
    n, maxmol, rmin, rmax, nsamp, seed = args
    rng = random.Random(seed)
    out = []
    for _ in range(nsamp):
        rx = draw(n, maxmol, rng.randint(rmin, rmax), rng)
        N = Network(n, rx)
        if not (N.rows_pairwise_independent() and N.columns_pairwise_independent()):
            continue
        if not (N.is_conservative() and N.is_consistent()):
            continue
        tw = twin_partition(N)
        if all(len(b) == 1 for b in tw):
            continue
        MS = frozenset(frozenset(s) for s in N.minimal_siphons())
        P = frozenset(frozenset(b) for b in N.dual().preclusters())
        out.append((rx, tw, MS, P))
    return out


def permset(fam, pi):
    return frozenset(frozenset(pi[x] for x in B) for B in fam)


def main():
    n, maxmol, rmin, rmax, nsamp, seed = (int(a) for a in sys.argv[1:7])
    nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 8
    t0 = time.time()
    print(f"# k8s_search.py n={n} mol<={maxmol} r in [{rmin},{rmax}] samples={nsamp}x{nproc} seed={seed}; "
          f"Python {platform.python_version()} on {platform.system()}")
    with Pool(nproc) as pool:
        res = pool.map(worker, [(n, maxmol, rmin, rmax, nsamp, seed + 1000 * i) for i in range(nproc)])
    nets = {}
    for chunk in res:
        for rx, tw, MS, P in chunk:
            nets.setdefault(rx, (tw, MS, P))
    print(f"# strong (H1)-(H3) networks with a twin pair: {len(nets)}")
    perms = list(permutations(range(n)))
    table = {}
    wit = None
    for rx, (tw, MS, P) in nets.items():
        for pi in perms:
            key = (permset(tw, pi), permset(MS, pi))
            d = table.setdefault(key, {})
            pp = permset(P, pi)
            if pp not in d:
                d[pp] = rx
                if len(d) >= 2 and wit is None:
                    wit = (key, list(d.values())[:2])
    print(f"# labelled keys: {len(table)}; witness: {'YES' if wit else 'none'}")
    if wit:
        key, (a, b) = wit
        print(f"#   N1 = {Network(n, a).describe()}")
        print(f"#   N2 = {Network(n, b).describe()}")
        print(f"#   twins {sorted(map(sorted, key[0]))}, minimal siphons {sorted(map(sorted, key[1]))}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
