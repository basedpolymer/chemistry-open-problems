"""Agent P's K8: are preclusters of N* determined by the twin classes and the minimal siphons of N?

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (records of
enum_all.py / enum_wr.py).  Usage: python k8.py RECORD_FILES...

Side N: among networks with (H1)-(H3), labelled on the species set {0..n-1} (every orbit representative under every
permutation of the species), look for two networks with the same n, the same partition into twin classes (equal
nonzero rows of Ye) and the same family of minimal siphons, but different P(N*).
Side D: the same on the dual side: networks N with (H1)-(H3) on the reaction set {0..r-1} (every permutation of
the reactions), same doubling classes (equal nonzero sources), same minimal siphons of N* (minimal T-traps of N),
different P(N).  Prints the smallest witness pair of each side (by total size) and its count of colliding keys.
Claim checked (agent P, proved by K3(e)): no witness pair consists of two WF-FC networks -- reported, the WF-FC
status of witnesses is printed.  The script never fails (K8 is a search); it exits 0.
"""
import sys
import time
import gzip
import pickle
import platform
from itertools import permutations, combinations

import ecore
from ecore import CNet


def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                seen.setdefault((r["n"], r["rx"]), r)
    return list(seen.values())


def twin_partition(rows):
    cl = {}
    for i, row in enumerate(rows):
        if any(row):
            cl.setdefault(row, []).append(i)
        else:
            cl.setdefault(("zero", i), []).append(i)
    return frozenset(frozenset(b) for b in cl.values())


def permset(fam, pi):
    return frozenset(frozenset(pi[x] for x in B) for B in fam)


def main():
    t0 = time.time()
    print(f"# k8.py  Python {platform.python_version()} on {platform.system()} {platform.release()}")
    recs = [r for r in load(sys.argv[1:])
            if r["N"]["cons"] and r["N"]["consist"] and r["N"]["rows_distinct"]]
    print(f"# {len(recs)} networks with (H1)-(H3)")
    for side in ("N", "D"):
        table = {}
        for idx, r in enumerate(recs):
            n, rx = r["n"], r["rx"]
            if side == "N":
                m = n
                rows = [tuple(a[i] for a, b in rx) for i in range(n)]
                tw = twin_partition(rows)
                MS = frozenset(r["N"]["MS"])
                P = frozenset(r["D"]["P"])
            else:
                m = r["r"]
                rows = [a for a, b in rx]
                tw = twin_partition(rows)
                MS = frozenset(r["D"]["MS"])
                P = frozenset(r["N"]["P"])
            if m > 6:
                continue
            for pi in permutations(range(m)):
                key = (m, permset(tw, pi), permset(MS, pi))
                val = permset(P, pi)
                d = table.setdefault(key, {})
                if val not in d:
                    d[val] = idx
        wit = []
        for key, d in table.items():
            if len(d) > 1:
                idxs = sorted(set(d.values()))
                if len(idxs) >= 2:
                    a, b = idxs[0], idxs[1]
                    size = recs[a]["n"] + recs[a]["r"] + recs[b]["n"] + recs[b]["r"]
                    wit.append((size, a, b, key))
                else:
                    # the same network under two labellings gives different P: impossible (equivariance)
                    wit.append((0, idxs[0], idxs[0], key))
        print(f"# side {side}: {len(table)} keys, {len(wit)} keys with two different preclusterings")
        if wit:
            wit.sort(key=lambda t: (t[0], t[1], t[2]))
            size, a, b, key = wit[0]
            for lab, i in (("first", a), ("second", b)):
                r = recs[i]
                N = CNet(r["n"], r["rx"])
                sd = None
                print(f"#   {lab}: N = {N.describe()}   N* = {N.dual().describe()}   WR={r['N']['WR']}")
            print(f"#   common key (ground-set size, twin/doubling classes, minimal siphons): "
                  f"{key[0]}, {sorted(map(sorted, key[1]))}, {sorted(map(sorted, key[2]))}")
            nwffc = 0
            sys.path.insert(0, __file__.rsplit("/", 1)[0])
            import tests
            for size, a, b, key in wit:
                fa = tests.Side(recs[a], side).WFFC()
                fb = tests.Side(recs[b], side).WFFC()
                nwffc += fa and fb
            print(f"{'PASS' if nwffc == 0 else 'FAIL'} side {side}: no witness pair of two WF-FC networks "
                  f"({nwffc} found)")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
