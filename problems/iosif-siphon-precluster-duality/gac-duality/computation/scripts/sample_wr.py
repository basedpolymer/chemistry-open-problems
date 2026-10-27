"""Seeded random sample of weakly reversible networks in intermediate dimensions (beyond exhaustive reach).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python sample_wr.py N_SPECIES MAX_MOLECULARITY RMIN RMAX COUNT SEED [--homog] [--s 2,3,4] [--cons]
                           [--nproc K] [--out PATH]

Sampling scheme (not uniform; deterministic given the seed): repeatedly draw a reaction budget R in [RMIN, RMAX],
then add linkage classes until the budget is used: a class has k in [2, min(5, budget)] complexes and e in
[k, min(budget, k+2)] reactions, its arc set is a uniformly chosen labelled strongly connected digraph with k
vertices and e arcs (ecore.sc_patterns) and its complexes are k distinct unused complexes chosen uniformly.
Keep the network if every species occurs, s is in the target set, it is conservative (with --cons), and its
canonical form is new.  Stop after COUNT networks or 400*COUNT draws.  Records as in enum_wr.py.
Claim checked: C1 every sampled network is weakly reversible, has the target s and (with --cons) satisfies (H1).
"""
import sys
import time
import gzip
import pickle
import random
import platform
from multiprocessing import Pool

from ecore import CNet, canonical, complexes, sc_patterns, record


def draw(rng, C, pats, Rmin, Rmax):
    R = rng.randint(Rmin, Rmax)
    used = set()
    rx = []
    budget = R
    while budget >= 2:
        k = rng.randint(2, min(5, budget))
        e = rng.randint(k, min(budget, k + 2))
        P = pats.get((k, e))
        if not P:
            continue
        free = [i for i in range(len(C)) if i not in used]
        if len(free) < k:
            break
        sub = rng.sample(free, k)
        pat = rng.choice(P)
        used.update(sub)
        rx.extend((C[sub[i]], C[sub[j]]) for i, j in pat)
        budget -= e
        if rng.random() < 0.25:
            break
    return tuple(sorted(rx))


def rec_worker(args):
    rx, n = args
    return record(rx, n)


def main():
    argv = sys.argv[1:]
    n, m, Rmin, Rmax, count, seed = map(int, argv[:6])
    homog = "--homog" in argv
    cons = "--cons" in argv
    starget = {int(x) for x in argv[argv.index("--s") + 1].split(",")} if "--s" in argv else {2, 3, 4}
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 8
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    t0 = time.time()
    print(f"# sample_wr.py n={n} m={m} R in [{Rmin},{Rmax}] count={count} seed={seed} homog={homog} cons={cons} "
          f"s in {sorted(starget)}; Python {platform.python_version()} on {platform.system()}")
    C = complexes(n, m)
    if homog:
        C = [c for c in C if sum(c) == m]
    pats = {}
    for k in range(2, 6):
        allp = sc_patterns(k, k + 2)
        for p in allp:
            pats.setdefault((k, len(p)), []).append(p)
    rng = random.Random(seed)
    seen = set()
    reps = []
    draws = 0
    while len(reps) < count and draws < 400 * count:
        draws += 1
        rx = draw(rng, C, pats, Rmin, Rmax)
        if not rx or not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
            continue
        cf = canonical(n, rx)
        if cf in seen:
            continue
        N = CNet(n, cf)
        if N.rank() not in starget:
            continue
        if cons and not N.is_conservative():
            continue
        seen.add(cf)
        reps.append(cf)
    print(f"# {len(reps)} distinct networks kept after {draws} draws ({time.time() - t0:.1f} s)")
    with Pool(nproc) as pool:
        recs = pool.map(rec_worker, [(rx, n) for rx in reps], chunksize=20)
    for r in recs:
        r["family"] = ("sample", n, m, Rmin, Rmax, seed, "homog" if homog else "", "cons" if cons else "")
    ok = all(r["N"]["WR"] and r["s"] in starget and (r["N"]["cons"] or not cons) for r in recs)
    bys = {}
    for r in recs:
        h = r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]
        d = bys.setdefault(r["s"], [0, 0, 0, 0])
        d[0] += 1
        d[1] += h
        d[2] += h and r["N"]["deficiency"] == 0
        d[3] += h and r["D"]["WR"]
    for s in sorted(bys):
        print(f"#   s={s}: {bys[s][0]} networks, {bys[s][1]} with (H1)-(H3), {bys[s][2]} of these with deficiency 0, "
              f"{bys[s][3]} with WR dual")
    if out:
        with gzip.open(out, "wb") as fh:
            pickle.dump(recs, fh, protocol=5)
        print(f"# records written to {out}")
    print(f"{'PASS' if ok else 'FAIL'} C1: all {len(recs)} sampled networks WR, s in target"
          + (", conservative" if cons else ""))
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
