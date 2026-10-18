#!/usr/bin/env python3
"""04_verify_certificates.py -- independent re-check of the sweep certificates.

Usage: python 04_verify_certificates.py STRIDE FILE [FILE ...]      (files may be .csv or .csv.gz)
       STRIDE = 1 re-checks every line; STRIDE = k re-checks every k-th line plus all lines with M < 200.
Environment: WSL Fedora, Python 3.12, python-flint 0.9.  Runtime: ~20 ms per re-checked line.

For each file (one branch j), the script checks
  coverage : every integer M from the first to the last line is present; the only line without a Hopf point is the
             first one (no crossing for that M, i.e. M = 14, 43, 86 for j = 1, 2, 3; float scan, K <= 1e5);
  claims   : from the stored midpoints ONLY (no continuation), Newton + Krawczyk re-prove a unique Hopf point in a
             small box, l1 is re-enclosed, and
               - for (M, j) = (15, 1):  l1 > 0  (subcritical),
               - for every other line:  l1 < 0  (supercritical),
             together with simplicity (0 not in dQ/dz) and transversality (Re dlambda/dkappa > 0).
Prints PASS/FAIL per file and exits non-zero on any failure.
"""
import sys, csv, gzip, time
from flint import arb, ctx
sys.path.insert(0, __import__("os").path.dirname(__file__))
import bdlib as bd

def rows(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        yield from csv.DictReader(fh)

def main():
    stride = int(sys.argv[1]); files = sys.argv[2:]
    total_fail = 0
    for fn in files:
        t0 = time.time(); nfail = 0; nchk = 0; Ms = []; first = True; jj = None
        for r in rows(fn):
            M, j = int(r["M"]), int(r["j"]); Ms.append(M); jj = j
            if r["status"] == "NO_CROSSING":
                if not first:
                    nfail += 1; print("  unexpected NO_CROSSING at M=%d" % M)
                first = False
                continue
            first = False
            if r["status"] not in ("PASS", "FAIL"):
                nfail += 1; print("  M=%d status %s" % (M, r["status"])); continue
            if not (M < 200 or M % stride == 0):
                continue
            nchk += 1
            X0 = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
            res = bd.certify_point(X0, M, prec=192)
            if res is None:
                nfail += 1; print("  M=%d: Krawczyk failed on re-check" % M); continue
            # the re-proved point must be the stored one
            same = abs(float(res["kappa"].mid()) - float(r["kappa_mid"])) < 1e-12*abs(float(r["kappa_mid"]))
            sign_ok = (res["l1"] > 0) if (M, j) == (15, 1) else (res["l1"] < 0)
            ok = same and sign_ok and res["simple"] and bool(res["transv"] > 0)
            if not ok:
                nfail += 1
                print("  M=%d j=%d: FAIL (same=%s, l1=%s, simple=%s, transv=%s)" % (M, j, same, res["l1"].str(5),
                                                                                  res["simple"], res["transv"].str(5)))
        cover = Ms == list(range(Ms[0], Ms[-1] + 1))
        if not cover:
            nfail += 1; print("  coverage gap in", fn)
        total_fail += nfail
        print("%s: j=%d, M in [%d, %d], %d lines re-checked, coverage %s: %s (%.0fs)"
              % (fn, jj, Ms[0], Ms[-1], nchk, cover, "PASS" if nfail == 0 else "FAIL", time.time() - t0), flush=True)
    print("SUMMARY:", "PASS" if total_fail == 0 else "FAIL")
    sys.exit(0 if total_fail == 0 else 1)

if __name__ == "__main__":
    main()
