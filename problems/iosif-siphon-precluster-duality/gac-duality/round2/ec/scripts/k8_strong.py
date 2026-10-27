"""K8, strong form (agent Q3, Claude Opus 5.5, effort xhigh; round2/ec/candidates.md Q3-K8s).

Search, among the first-round record files of agent E (WSL ~/gacE-data/*.pkl.gz, produced by
computation/scripts/enum_all.py, enum_wr.py, enum_fc.py, sample_wr.py), for two networks with
  * (H1)-(H3), rows of Gamma pairwise linearly independent and columns of Gamma pairwise linearly independent
    (record fields rows_pi, cols_pi), on both networks;
  * side N: the same number of species, the same twin classes (equal nonzero rows of Ye) and the same minimal
    siphons (up to one common relabelling of the species), but different P(N*);
  * side D: the same number of reactions, the same doubling classes and the same minimal T-traps (minimal siphons
    of N*), up to a common relabelling of the reactions, but different P(N).
Method: bucket by a permutation-invariant signature, then decide equality of the labelled keys exactly by trying
all permutations inside a bucket.  Environment: Python 3.12, standard library only.  Usage:
[PURE=1] python k8_strong.py FILES...   Exit status 0 (search, no claim).  Runtime printed.
"""
import os
import sys
import time
import gzip
import pickle
import platform
from itertools import permutations


def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                seen.setdefault((r["n"], tuple(r["rx"])), r)
    return list(seen.values())


def twin_partition(rows):
    cl = {}
    for i, row in enumerate(rows):
        cl.setdefault(tuple(row) if any(row) else ("zero", i), []).append(i)
    return frozenset(frozenset(b) for b in cl.values())


def permset(fam, pi):
    return frozenset(frozenset(pi[x] for x in B) for B in fam)


def signature(m, tw, MS, P):
    # invariant of (tw, MS) under relabelling: for each element, (size of its twin block, sorted sizes of the
    # minimal siphons containing it); plus the global multisets
    per = sorted((len(next(B for B in tw if x in B)), tuple(sorted(len(S) for S in MS if x in S))) for x in range(m))
    return (m, tuple(per), tuple(sorted(len(S) for S in MS)))


def main():
    t0 = time.time()
    print(f"# k8_strong.py  Python {platform.python_version()} on {platform.system()}")
    recs = load(sys.argv[1:])
    print(f"# {len(recs)} distinct networks in the record files")
    strong = [r for r in recs if r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]
              and r.get("rows_pi") and r.get("cols_pi")]
    if os.environ.get("PURE") == "1":
        strong = [r for r in strong if all(not (a[i] and b[i]) for a, b in r["rx"] for i in range(r["n"]))]
        print("# PURE=1: only pure networks (no species on both sides of a reaction)")
    print(f"# {len(strong)} with (H1)-(H3) and rows, columns of Gamma pairwise linearly independent")
    for side in ("N", "D"):
        buckets = {}
        for idx, r in enumerate(strong):
            rx = r["rx"]
            if side == "N":
                m = r["n"]
                tw = twin_partition([tuple(a[i] for a, b in rx) for i in range(m)])
                MS = frozenset(frozenset(S) for S in r["N"]["MS"])
                P = frozenset(frozenset(B) for B in r["D"]["P"])
            else:
                m = r["r"]
                tw = twin_partition([a for a, b in rx])
                MS = frozenset(frozenset(S) for S in r["D"]["MS"])
                P = frozenset(frozenset(B) for B in r["N"]["P"])
            if all(len(B) == 1 for B in tw):
                continue            # no twins / doubling pair: P of the other side is discrete (manuscript Cor. 9)
            buckets.setdefault(signature(m, tw, MS, P), []).append((idx, m, tw, MS, P))
        cand = 0
        witnesses = []
        for sig, items in buckets.items():
            if len(items) < 2:
                continue
            m = sig[0]
            if m > 6:
                continue
            # exact: labelled keys under all permutations of the first item of each P-class
            table = {}
            for idx, m_, tw, MS, P in items:
                key0 = (tw, MS)
                for pi in permutations(range(m)):
                    key = (permset(tw, pi), permset(MS, pi))
                    table.setdefault(key, {}).setdefault(permset(P, pi), idx)
                    if len(table) > 2_000_000:
                        break
            for key, d in table.items():
                ids = sorted(set(d.values()))
                if len(ids) >= 2:
                    witnesses.append((ids[0], ids[1], key))
                    cand += 1
        print(f"# side {side}: {len(buckets)} buckets; colliding labelled keys: {cand}")
        shown = set()
        for a, b, key in sorted(witnesses, key=lambda w: (strong[w[0]]["n"] + strong[w[0]]["r"]
                                                           + strong[w[1]]["n"] + strong[w[1]]["r"], w[0], w[1])):
            if (a, b) in shown:
                continue
            shown.add((a, b))
            for lab, i in (("first", a), ("second", b)):
                r = strong[i]
                print(f"#   {lab}: n={r['n']} r={r['r']} rx={r['rx']} family={r.get('family')} WR={r['N']['WR']}")
            print(f"#   key: classes {sorted(map(sorted, key[0]))}, minimal sets {sorted(map(sorted, key[1]))}")
            if len(shown) >= 3:
                break
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
