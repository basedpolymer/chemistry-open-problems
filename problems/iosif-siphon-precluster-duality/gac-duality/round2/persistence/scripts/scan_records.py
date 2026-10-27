"""Scan agent E's weakly reversible record files: drain species (E.4) versus drain sets (Theorem D) in the dual, and
the rank-one classification of the dual (Theorems S1, S1-WR) (agent Q1, round 2).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact arithmetic.  Record files: agent E's pickles
in the WSL directory ~/gacE-data (SHA-256 in ../../../computation/logs/data_sha256.txt).
Usage: python scan_records.py FILE [FILE ...]   (paths of wr_*.pkl.gz record files)

Checks (PASS/FAIL, exit 1 on FAIL):
  R1  every drain species of E is a drain set (|Z| = 1) of our implementation;
  R2  Theorem D implies non-persistence for all rates; for s = 1 duals this must agree with the face
      classification: a dual with a drain set has verdict NONE;
  R3  Theorem S1-WR: no DOM face in the dual of a WR network with s = 1.
Reports, for WR networks with (H1)-(H3), by s: duals with a drain species, with a drain set, with a drain set but
no drain species (new witnesses; smallest printed), and the s = 1 verdict counts.
"""
import gzip
import os
import pickle
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from drainsets import drain_sets, drain_species_E  # noqa: E402
from s1lib import classify, dual, show  # noqa: E402


def H(rec):
    N = rec["N"]
    return N["cons"] and N["consist"] and N["rows_distinct"]


def main():
    t0 = time.time()
    print(f"# scan_records.py, Python {platform.python_version()} on {platform.system()}")
    fails = {"R1": 0, "R2": 0, "R3": 0}
    tot = {}
    newest = None
    for path in sys.argv[1:]:
        L = pickle.load(gzip.open(path))
        fam = os.path.basename(path)
        cnt = {}
        for rec in L:
            if not (rec["N"]["WR"] and H(rec)):
                continue
            n, rx, s = rec["n"], [tuple(map(tuple, x)) for x in rec["rx"]], rec["s"]
            p, drx = dual(n, rx)
            ds = drain_species_E(p, drx)
            dsets = drain_sets(p, drx, maxsize=4, first_only=False)
            singles = {Z[0] for Z, d, ch in dsets if len(Z) == 1}
            if not set(ds) <= singles:
                fails["R1"] += 1
                print("FAIL R1", show(n, rx))
            c = cnt.setdefault(s, {"nets": 0, "drainsp": 0, "drainset": 0, "new": 0, "ALL": 0, "NONE": 0, "other": 0})
            c["nets"] += 1
            c["drainsp"] += bool(ds)
            c["drainset"] += bool(dsets)
            if dsets and not ds:
                c["new"] += 1
                size = (len(rx), s, sum(sum(a) + sum(b) for a, b in rx))
                if newest is None or size < newest[0]:
                    newest = (size, n, rx, dsets[0])
            if s == 1:
                v, data = classify(p, drx)
                if v in ("ALL", "NONE"):
                    c[v] += 1
                else:
                    c["other"] += 1
                if data and any(f["type"] == "DOM" for f in data[2]):
                    fails["R3"] += 1
                    print("FAIL R3", show(n, rx))
                if dsets and v != "NONE":
                    fails["R2"] += 1
                    print("FAIL R2", show(n, rx))
        for s, c in sorted(cnt.items()):
            print(f"{fam:28s} s={s} WR+H={c['nets']:6d} drain-species={c['drainsp']:5d} drain-set={c['drainset']:5d} "
                  f"set-not-species={c['new']:5d}" + (f" s1: ALL={c['ALL']} NONE={c['NONE']} other={c['other']}" if s == 1 else ""))
            T = tot.setdefault(s, {})
            for k, x in c.items():
                T[k] = T.get(k, 0) + x
    print("# totals by s:", {s: T for s, T in sorted(tot.items())})
    if newest:
        size, n, rx, ex = newest
        print(f"# smallest WR network whose dual has a drain set but no drain species: N = {show(n, rx)}; "
              f"N* = {show(*dual(n, rx), names=[f'k{j+1}' for j in range(len(rx))])}; drain set Z (0-based) = {ex[0]}, d = {ex[1]}")
    for k, v in fails.items():
        print(f"{k}: {'PASS' if v == 0 else 'FAIL'} ({v} failures)")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if any(fails.values()) else 0)


if __name__ == "__main__":
    main()
