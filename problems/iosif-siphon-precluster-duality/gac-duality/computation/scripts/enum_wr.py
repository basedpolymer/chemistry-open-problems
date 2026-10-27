"""Exhaustive (or seeded-sample) enumeration of weakly reversible networks and of all objects of N and N*.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python enum_wr.py N_SPECIES MAX_MOLECULARITY MAX_REACTIONS [--zero | --homog] [--nproc K]
                         [--sample K --seed S] [--out PATH]
--homog: only complexes of molecularity exactly m (all such networks are conservative).

Networks: reaction sets whose reaction graph is weakly reversible (each linkage class strongly connected), complexes
of molecularity <= m (the zero complex only with --zero), at most R reactions, every species occurring, one
representative per orbit of the species permutations (canonical form of ecore.canonical).  Generation: disjoint
unions of linkage classes, each a labelled strongly connected spanning digraph on a set of complexes.

Phase 1 enumerates all orbit representatives (cheap).  Phase 2 computes ecore.record for every representative, or
for a seeded random sample of K of them (--sample), and pickles the list of records (gzip) to --out.  The summary
printed below is by stoichiometric dimension s.  Claims checked here (PASS/FAIL):
  C1  every generated network is weakly reversible (crn.Network.weakly_reversible) and has <= R reactions;
  C2  representatives are pairwise non-isomorphic (distinct canonical forms) - true by construction, re-checked;
  C3  Lemma 1 of result 17 on every network: Gamma(N*) = Gamma^T, rank equal, N** = N, (H1)(N) <=> (H2)(N*),
      (H2)(N) <=> (H1)(N*), (H3)(N*) <=> columns of Gamma distinct.
"""
import sys
import time
import gzip
import pickle
import random
import platform
from multiprocessing import Pool

import ecore
from ecore import CNet, canonical, complexes, lc_instances, wr_networks_from, record

INST = None
ARGS = None


def phase1(firsts):
    n, R = ARGS["n"], ARGS["R"]
    out = []
    for f in firsts:
        for rx in wr_networks_from(f, INST, R, n):
            if canonical(n, rx) == rx:
                out.append(rx)
    return out


def phase2(batch):
    n = ARGS["n"]
    out = []
    bad = []
    for rx in batch:
        rec = record(rx, n)
        N = CNet(n, rx)
        D = N.dual()
        ok = N.weakly_reversible() and len(rx) <= ARGS["R"]
        G = N.Gamma
        GT = [list(c) for c in zip(*G)]
        ok3 = (D.Gamma == GT and D.rank() == N.rank() and
               [tuple(map(tuple, x)) for x in ecore.Network.dual(D).reactions] == [tuple(map(tuple, x)) for x in N.reactions] and
               rec["N"]["cons"] == rec["D"]["consist"] and rec["N"]["consist"] == rec["D"]["cons"] and
               rec["D"]["rows_distinct"] == N.columns_distinct())
        if not ok:
            bad.append(("C1", rx))
        if not ok3:
            bad.append(("C3", rx))
        out.append(rec)
    return out, bad


def main():
    global INST, ARGS
    argv = sys.argv[1:]
    n, m, R = int(argv[0]), int(argv[1]), int(argv[2])
    zero = "--zero" in argv
    homog = "--homog" in argv
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 8
    sample = int(argv[argv.index("--sample") + 1]) if "--sample" in argv else None
    seed = int(argv[argv.index("--seed") + 1]) if "--seed" in argv else 0
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    ARGS = dict(n=n, m=m, R=R)
    t0 = time.time()
    print(f"# enum_wr.py n={n} m={m} R={R} zero={zero} homog={homog} nproc={nproc} sample={sample} seed={seed}")
    print(f"# Python {platform.python_version()} on {platform.system()} {platform.release()}")
    C = complexes(n, m, zero)
    if homog:  # only complexes of molecularity exactly m (then (1,...,1) is a conservation law)
        C = [c for c in C if sum(c) == m]
    INST = lc_instances(C, R)
    print(f"# {len(C)} complexes, {len(INST)} linkage-class instances")
    chunks = [list(range(i, len(INST), nproc * 8)) for i in range(nproc * 8)]
    reps = []
    with Pool(nproc) as pool:
        for part in pool.imap_unordered(phase1, chunks):
            reps.extend(part)
    reps.sort()
    t1 = time.time()
    print(f"# phase 1: {len(reps)} orbit representatives in {t1 - t0:.1f} s")
    c2 = len(set(reps)) == len(reps)
    if sample is not None and sample < len(reps):
        rng = random.Random(seed)
        reps = sorted(rng.sample(reps, sample))
        print(f"# seeded random sample of {len(reps)} representatives (seed {seed})")
    batches = [reps[i:i + 50] for i in range(0, len(reps), 50)]
    recs, bad = [], []
    with Pool(nproc) as pool:
        for rr, bb in pool.imap(phase2, batches):
            recs.extend(rr)
            bad.extend(bb)
    t2 = time.time()
    print(f"# phase 2: {len(recs)} records in {t2 - t1:.1f} s")
    for r in recs:
        r["family"] = (n, m, R, "zero" if zero else ("homog" if homog else ""))
    if out:
        with gzip.open(out, "wb") as fh:
            pickle.dump(recs, fh, protocol=5)
        print(f"# records written to {out}")
    # summary by s
    smax = max((r["s"] for r in recs), default=0)
    hdr = ["s", "nets", "H1", "H2", "H3", "H", "H&H*", "def0", "H&def0", "WR*", "H&WR*", "H&def0&WR*",
           "H&crit", "H&EC", "H&EFC", "H&FC"]
    print("# counts by stoichiometric dimension s (H = (H1)-(H3) for N; H* = (H1)-(H3) for N*; crit = some critical "
          "minimal siphon; EC/EFC/FC = equal conflict / extended free choice / free choice)")
    print("  ".join(f"{h:>10s}" for h in hdr))
    tot = [0] * (len(hdr) - 1)
    for s in range(0, smax + 1):
        rs = [r for r in recs if r["s"] == s]
        if not rs:
            continue
        def h(r):
            return r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]
        def hd(r):
            return r["D"]["cons"] and r["D"]["consist"] and r["D"]["rows_distinct"]
        row = [len(rs),
               sum(r["N"]["cons"] for r in rs), sum(r["N"]["consist"] for r in rs),
               sum(r["N"]["rows_distinct"] for r in rs), sum(h(r) for r in rs),
               sum(h(r) and hd(r) for r in rs), sum(r["N"]["deficiency"] == 0 for r in rs),
               sum(h(r) and r["N"]["deficiency"] == 0 for r in rs), sum(r["D"]["WR"] for r in rs),
               sum(h(r) and r["D"]["WR"] for r in rs),
               sum(h(r) and r["N"]["deficiency"] == 0 and r["D"]["WR"] for r in rs),
               sum(h(r) and len(r["N"]["CS"]) > 0 for r in rs),
               sum(h(r) and r["N"]["EC"] for r in rs), sum(h(r) and r["N"]["EFC"] for r in rs),
               sum(h(r) and r["N"]["FC"] for r in rs)]
        tot = [a + b for a, b in zip(tot, row)]
        print("  ".join(f"{x:>10}" for x in [s] + row))
    print("  ".join(f"{x:>10}" for x in ["all"] + tot))
    ok1 = not any(b[0] == "C1" for b in bad)
    ok3 = not any(b[0] == "C3" for b in bad)
    print(f"{'PASS' if ok1 else 'FAIL'} C1: all {len(recs)} networks weakly reversible with <= {R} reactions")
    print(f"{'PASS' if c2 else 'FAIL'} C2: orbit representatives pairwise distinct")
    print(f"{'PASS' if ok3 else 'FAIL'} C3: Gamma(N*) = Gamma^T, same rank, N** = N, (H1)/(H2) swap, (H3)(N*) <=> "
          f"distinct columns, on all {len(recs)} networks")
    for b in bad[:10]:
        print("  failure:", b)
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if (ok1 and c2 and ok3) else 1)


if __name__ == "__main__":
    main()
