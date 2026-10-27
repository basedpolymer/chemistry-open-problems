"""Networks with the cluster rank equation s = #clusters - 1 outside EC (agent Q3, Claude Opus 5.5, effort xhigh;
round2/ec/notes.md §7.3).  Environment: Python 3.12, standard library only, exact (crn.py).
Usage: python cluster_rank.py FILES...  (agent E's record files, WSL ~/gacE-data/all_*.pkl.gz)

Selects networks N with (H1)-(H3), connected, s = (number of Petri clusters with a reaction) - 1 that are not EC.
For each, deletes the pure catalysts (species with a zero row of Gamma) and classifies the reduced network: EC,
block-homogeneous (BH: TEFC with every source a multiple of a 0/1 vector), or neither ('genuine').  Reports counts
and the smallest genuine examples with their minimal siphons, critical siphons of N and of N*, and the rank data.
Exit 0 (report).
"""
import sys
import time
import gzip
import pickle
import platform

import eclib
from eclib import Network


def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                seen.setdefault((r["n"], tuple(r["rx"])), r)
    return list(seen.values())


def classify(N):
    srcs = sorted({a for a, b in N.reactions})
    supp = [frozenset(i for i in range(N.n) if y[i]) for y in srcs]
    if all(not (supp[i] & supp[j]) for i in range(len(srcs)) for j in range(i + 1, len(srcs))):
        return "EC"
    tefc = all(not (supp[i] & supp[j]) or supp[i] == supp[j] for i in range(len(srcs)) for j in range(i + 1, len(srcs)))
    bh = all(len({c for c in y if c}) == 1 for y in srcs)
    if tefc and bh:
        return "BH"
    return "TEFC" if tefc else "other"


def delete_catalysts(N):
    G = N.Gamma
    keep = [i for i in range(N.n) if any(G[i])]
    rx = [(tuple(a[i] for i in keep), tuple(b[i] for i in keep)) for a, b in N.reactions]
    return Network(len(keep), rx), len(keep) < N.n


def main():
    t0 = time.time()
    print(f"# cluster_rank.py; Python {platform.python_version()} on {platform.system()}")
    recs = load(sys.argv[1:])
    counts = {}
    genuine = []
    for r in recs:
        if not (r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]):
            continue
        N = Network(r["n"], r["rx"])
        if not eclib.petri.is_connected(N):
            continue
        ncl = len([c for c in eclib.petri.clusters(N) if c[1]])
        if N.rank() != ncl - 1:
            continue
        k0 = classify(N)
        if k0 == "EC":
            continue
        R, had_cat = delete_catalysts(N)
        k1 = classify(R)
        key = (k0, k1 if had_cat else "no catalyst")
        counts[key] = counts.get(key, 0) + 1
        if k1 in ("TEFC", "other"):
            genuine.append((N.n + N.r, N, r))
    for key in sorted(counts):
        print(f"REPORT  kind {key[0]:6s} -> after deleting pure catalysts: {key[1]:12s} {counts[key]}")
    genuine.sort(key=lambda t: (t[0], t[1].describe()))
    print(f"# genuine (not EC/BH even after deleting pure catalysts): {len(genuine)}")
    for size, N, r in genuine[:8]:
        print(f"#   {N.describe()}  s={N.rank()}  critical siphons N: {len(eclib.petri.critical_siphons(N))}, "
              f"N*: {len(eclib.petri.critical_siphons(N.dual()))}  record siphons_eq_SF N/D: "
              f"{r['N']['siphons_eq_SF']}/{r['D']['siphons_eq_SF']}")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
