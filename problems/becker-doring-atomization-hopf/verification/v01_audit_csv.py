#!/usr/bin/env python3
"""v01_audit_csv.py -- verifier's audit of the six sweep certificate files (no re-proof, data consistency only).

Usage: python v01_audit_csv.py ../certificates/sweep_j*.csv.gz
Checks, per file:
  - coverage of every integer M, single NO_CROSSING line at the start;
  - stored enclosures: omega_lo > 0 (needed: Kuznetsov's formula assumes omega > 0), y = Im z > 0 (phi non-real),
    l1 (unit) and l1 (q1=1) of the same strict sign, the sign claimed (positive only at (15,1)), simple == 1,
    transv_lo > 0, status consistent with these;
  - smoothness of the branch in M (a jump to another branch would show as an outlier in the second differences
    of kappa and Im z);
  - Im z inside the window ((2j-1) pi, 2 j pi) used by the solver as a branch label aid.
Environment: WSL Fedora, Python 3.12. Runtime: ~1 min.
"""
import sys, csv, gzip, math
import numpy as np

def main():
    bad = 0
    for fn in sys.argv[1:]:
        rows = list(csv.DictReader(gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)))
        Ms = [int(r["M"]) for r in rows]
        j = int(rows[0]["j"])
        cover = Ms == list(range(Ms[0], Ms[-1] + 1))
        nc = [r for r in rows if r["status"] == "NO_CROSSING"]
        data = [r for r in rows if r["status"] in ("PASS", "FAIL")]
        other = [r for r in rows if r["status"] not in ("PASS", "FAIL", "NO_CROSSING")]
        errs = []
        kap, y, MM = [], [], []
        for r in data:
            M = int(r["M"])
            om_lo = float(r["omega_lo"]); l_lo = float(r["l1_unit_lo"]); l_hi = float(r["l1_unit_hi"])
            q_lo = float(r["l1_q1_lo"]); q_hi = float(r["l1_q1_hi"])
            yy = float(r["y_mid"]); ry = float(r["box_rad_y"])
            sub = (M, j) == (15, 1)
            ok = om_lo > 0 and yy - ry > 0 and r["simple"] == "1" and float(r["transv_lo"]) > 0
            ok &= (l_lo > 0 and q_lo > 0) if sub else (l_hi < 0 and q_hi < 0)
            ok &= (r["status"] == ("FAIL" if sub else "PASS"))
            if not ok:
                errs.append(M)
            kap.append(float(r["kappa_mid"])); y.append(yy); MM.append(M)
        kap = np.array(kap); y = np.array(y); MM = np.array(MM)
        win = np.all(((2*j - 1)*math.pi < y) & (y < 2*j*math.pi))
        mono = np.all(np.diff(kap) < 0)
        # relative second differences (branch jumps would be outliers)
        d2k = np.abs(np.diff(kap, 2))/kap[1:-1]
        d2y = np.abs(np.diff(y, 2))
        # normalise by the local scale of the first differences
        d1k = np.abs(np.diff(kap))/kap[1:]
        ratio = np.max(d2k[5:]/np.maximum(d1k[5:-1], 1e-300)) if len(d2k) > 10 else float("nan")
        # small files start with the single NO_CROSSING line; the "_b" files continue at M = 10001 without one
        start_ok = (len(nc) == 1 and rows[0]["status"] == "NO_CROSSING") or (len(nc) == 0 and Ms[0] == 10001)
        lines_ok = (not errs) and cover and start_ok and not other
        bad += not lines_ok
        print("%s: j=%d M=[%d,%d] lines=%d no_cross=%s errors=%d cover=%s | Im z window %s, kappa decreasing %s, "
              "max|d2 kappa|/|d1 kappa| (M>=first+6) = %.3g, max|d2 y| = %.3g | %s"
              % (fn, j, Ms[0], Ms[-1], len(data), [int(r["M"]) for r in nc], len(errs), cover, win, mono,
                 ratio, d2y.max(), "PASS" if lines_ok else "FAIL"))
        if errs:
            print("   lines with problems:", errs[:20])
        i0 = np.argmax(np.diff(kap) >= 0) if not mono else None
        if i0 is not None:
            print("   kappa not decreasing at M =", MM[i0], kap[i0:i0+3])
    print("SUMMARY:", "PASS" if bad == 0 else "FAIL")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
