#!/usr/bin/env python3
"""Tamper test for verify_certificate.py: corrupt a certificate in several ways and check that the verifier
rejects every corrupted version (exit code != 0).  Prints PASS/FAIL per corruption.
Corruptions: (1) increase certR of one leaf with positive outR; (2) decrease outR of one leaf;
(3) delete one leaf (tiling hole); (4) duplicate one leaf (overlap); (5) change the theta level of one leaf.
Environment: Python 3.13.  Runtime: a few seconds per corruption on a small certificate (reproduce_small.sh).
Usage: python3 tamper_test.py LEAVES.gz SUMMARY.json HEXCERT TRIG_TABLE
"""
import sys, gzip, os, subprocess, tempfile, json
HERE = os.path.dirname(os.path.abspath(__file__))

def write(rows, path):
    with gzip.open(path, "wt") as f:
        for r in rows: f.write(" ".join(map(str, r)) + "\n")

def main():
    leaves, summary, binp, table = sys.argv[1:5]
    rows = [list(map(int, l.split())) for l in gzip.open(leaves, "rt")]
    i = next(k for k, r in enumerate(rows) if r[7] > 0 and r[8] < r[7])
    variants = []
    r1 = [list(r) for r in rows]; r1[i][8] += 1; variants.append(("certR + 1", r1))
    r2 = [list(r) for r in rows]; r2[i][7] -= 1; variants.append(("outR - 1", r2))
    r3 = [list(r) for r in rows]; del r3[i]; variants.append(("leaf deleted", r3))
    r4 = [list(r) for r in rows]; r4.append(list(rows[i])); variants.append(("leaf duplicated", r4))
    r5 = [list(r) for r in rows]; r5[i][6] += 1; variants.append(("theta level changed", r5))
    ok = True
    with tempfile.TemporaryDirectory() as td:
        for name, rr in variants:
            p = os.path.join(td, "t.leaves.gz"); write(rr, p)
            idx = [str(i)]
            res = subprocess.run([sys.executable, os.path.join(HERE, "verify_certificate.py"), p, summary, binp, table,
                                  "--only", str(i)], capture_output=True, text=True)
            rejected = res.returncode != 0
            print(("PASS" if rejected else "FAIL"), "tampered certificate rejected:", name)
            ok = ok and rejected
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
