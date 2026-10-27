"""Exhaustive enumeration of ordinary equal-conflict (= free-choice in agent P's sense) networks with (H1), (H2).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python enum_fc.py N_SPECIES MAX_REACTIONS [--nproc K] [--out PATH]

Networks: every complex is a 0/1 vector (ordinary); the source complexes are the indicator vectors of pairwise
disjoint nonempty blocks of species (this is exactly equal conflict for ordinary networks: two reactions sharing a
reactant have the same source); each block is the source of at least one reaction; products are nonzero 0/1 vectors
different from the source; at most R reactions; every species occurs; one representative per species-permutation
orbit.  A record is kept when (H1) and (H2) hold.  Claim: C1 every kept network is ordinary and EC.
"""
import sys
import time
import gzip
import pickle
import platform
from itertools import combinations, product
from multiprocessing import Pool

from ecore import CNet, canonical, record

ARGS = None


def partial_partitions(n):
    """All sets of pairwise disjoint nonempty blocks of {0..n-1} (as tuples of sorted tuples)."""
    out = []

    def rec(i, blocks):
        if i == n:
            out.append(tuple(tuple(b) for b in blocks))
            return
        rec(i + 1, blocks)  # species i in no source
        for b in blocks:
            b.append(i)
            rec(i + 1, blocks)
            b.pop()
        blocks.append([i])
        rec(i + 1, blocks)
        blocks.pop()
    rec(0, [])
    return [p for p in out if p]


def worker(pp):
    n, R = ARGS
    vecs = [v for v in product((0, 1), repeat=n) if any(v)]
    srcs = [tuple(1 if i in b else 0 for i in range(n)) for b in pp]
    choices = [[v for v in vecs if v != s] for s in srcs]
    out = []
    k = len(srcs)
    if k > R:
        return out

    def rec(t, rx, budget):
        if t == k:
            rxs = tuple(sorted(rx))
            if not all(any(a[i] or b[i] for a, b in rxs) for i in range(n)):
                return
            for i in range(n):  # necessary for (H1)+(H2): each species gained and lost, or neither
                if any(b[i] < a[i] for a, b in rxs) != any(b[i] > a[i] for a, b in rxs):
                    return
            if canonical(n, rxs) != rxs:
                return
            N = CNet(n, rxs)
            if N.is_consistent() and N.is_conservative():
                out.append(record(rxs, n))
            return
        rest = k - t - 1  # each later block needs >= 1 reaction
        for j in range(1, budget - rest + 1):
            for prods in combinations(choices[t], j):
                rec(t + 1, rx + [(srcs[t], p) for p in prods], budget - j)
    rec(0, [], R)
    return out


def main():
    global ARGS
    argv = sys.argv[1:]
    n, R = int(argv[0]), int(argv[1])
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 8
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    ARGS = (n, R)
    t0 = time.time()
    print(f"# enum_fc.py n={n} R={R} nproc={nproc}; Python {platform.python_version()} on {platform.system()}")
    pps = partial_partitions(n)
    recs = []
    with Pool(nproc) as pool:
        for part in pool.imap_unordered(worker, pps):
            recs.extend(part)
    recs.sort(key=lambda r: r["rx"])
    for r in recs:
        r["family"] = ("fc", n, R)
    ok = len({r["rx"] for r in recs}) == len(recs)
    for r in recs:
        rx = r["rx"]
        ok &= all(x in (0, 1) for a, b in rx for x in a + b)
        ok &= all(a == c for (a, b), (c, d) in combinations(rx, 2) if any(x and y for x, y in zip(a, c)))
    nH = sum(r["N"]["rows_distinct"] for r in recs)
    print(f"# {len(recs)} ordinary EC networks with (H1),(H2); {nH} with (H3); "
          f"{sum(r['N']['WR'] for r in recs)} weakly reversible")
    if out:
        with gzip.open(out, "wb") as fh:
            pickle.dump(recs, fh, protocol=5)
        print(f"# records written to {out}")
    print(f"{'PASS' if ok else 'FAIL'} C1: kept networks pairwise distinct, ordinary and equal conflict")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
