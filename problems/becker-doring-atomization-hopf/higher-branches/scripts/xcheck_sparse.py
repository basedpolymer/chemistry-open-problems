#!/usr/bin/env python3
"""xcheck_sparse.py -- separate numerical cross-check of l1 at selected integer points, with the O(M) sparse code
of result 8's verifier (../../verification/v03_sparse_float.py, v04_sparse_arb.py; imported, not modified).

That code works in reduced coordinates on the mass hyperplane (tridiagonal + rank-one Jacobian), re-locates the
Hopf point by a bordered Newton iteration for the eigenpair and a secant in K, and evaluates Kuznetsov's formula,
all in Arb midpoint arithmetic at 320 bits.  It shares no code with bdlib / allm_core / l1reg.  This is a
NUMERICAL cross-check (midpoint iteration, no enclosure), used for the points where l1 > 0 is certified.

Usage:  python xcheck_sparse.py PW.csv M [M ...]
Claim (PASS/FAIL): at each M, kappa agrees with the certificate to 1e-20 (relative), l1 in both normalizations
to 1e-12 (relative), and the sign of l1 is the certified one.
Environment: python-flint 0.9, numpy, scipy.  Runtime: seconds to minutes per point.
"""
import sys, os, csv
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "verification"))
os.environ.setdefault("V04_PREC", "320")
from flint import arb
import v03_sparse_float as v3
import v04_sparse_arb as v4

def main():
    pw = sys.argv[1]; Ms = [int(a) for a in sys.argv[2:]]
    rows = {int(r["M"]): r for r in csv.DictReader(open(pw))}
    bad = 0
    for M in Ms:
        r = rows[M]; j = int(r["j"])
        # v4.run compares with the stored l1 and returns a flag for result 8's sign convention; recompute here
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            v4.run(M, j, {M: r})
        line = buf.getvalue().strip()
        print(line, flush=True)
        # parse "rel" values from the printed line
        import re
        rels = [float(x) for x in re.findall(r"rel ([0-9.e+-]+)", line)]
        l1u = float(re.search(r"l1_unit=([-0-9.e+]+)", line).group(1))
        cert_pos = float(r["l1_unit_lo"]) > 0
        ok = len(rels) == 3 and rels[0] < 1e-20 and rels[1] < 1e-12 and rels[2] < 1e-12 and ((l1u > 0) == cert_pos)
        print("M=%d j=%d: certified sign %s, sparse code sign %s, agreement %s -> %s" % (
            M, j, "+" if cert_pos else "-", "+" if l1u > 0 else "-", rels, "PASS" if ok else "FAIL"), flush=True)
        bad += not ok
    print("SUMMARY: %d points -> %s" % (len(Ms), "PASS" if bad == 0 else "FAIL"))
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
