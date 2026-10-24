#!/usr/bin/env python3
"""h2j_sweep.py -- spectral certificate (H2 and unstable count) at the j-th Hopf point on the eps-boxes of an
l1_sweep certificate file, re-proving the Krawczyk box from the stored data.

Usage:  python h2j_sweep.py J IN.csv OUT.csv R0 Y cL cR [EPS_STOP]
  IN.csv : l1_sweep.py output for branch J (boxes eps_lo, eps_hi, midpoints and radii of the Hopf box)
  R0, Y, cL, cR : contour parameters of h2j_lib (localization radius, half-height, left/right level values)
  EPS_STOP : stop after the first box with eps_hi > EPS_STOP (default: all boxes)
Environment: Python 3.13, python-flint 0.9.  1 thread.  Runtime about 2 s per box.
Claim per box (PASS): the Krawczyk box is re-proved for all eps in the box, and the checks of h2j_lib hold
(K < 1, localization, monotonicity, left/right chains, N(P) = 2J + 2, N(P-) = 4).  The output keeps the l1 columns
and adds the H2 columns.  SUMMARY line; exit 0 iff all processed boxes pass.
"""
import sys, os, csv, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import h2j_lib as HJ

def main():
    j = int(sys.argv[1]); fin, fout = sys.argv[2], sys.argv[3]
    R0 = float(sys.argv[4]); Y = float(sys.argv[5]); cL = float(sys.argv[6]); cR = float(sys.argv[7])
    estop = float(sys.argv[8]) if len(sys.argv) > 8 else float("inf")
    ctx.prec = 192
    rows = list(csv.DictReader(open(fin)))
    nf = n = 0; t0 = time.time()
    with open(fout, "w", newline="") as fh:
        fields = list(rows[0].keys()) + ["R0", "Y", "cL", "cR", "loc_mode", "mono_ub", "NP", "NPm", "h2_status"]
        wr = csv.DictWriter(fh, fieldnames=fields); wr.writeheader()
        for r in rows:
            if float(r["eps_lo"]) > estop:
                break
            n += 1
            eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
            Xm = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
            rad = [arb(r["rad_kappa"]), arb(r["rad_x"]), arb(r["rad_y"])]
            enc, box = C.krawczyk(Xm, rad, eb)
            if enc is None:
                res = dict(ok=False, fail="krawczyk")
            else:
                res = HJ.h2j_checks(j, eb, box, R0, Y, cL, cR)
            ok = res.get("ok", False)
            nf += (not ok)
            out = dict(r); out.update(R0=R0, Y=Y, cL=cL, cR=cR, loc_mode=res.get("loc_mode", ""),
                                      mono_ub="%.4g" % res["mono_ub"] if "mono_ub" in res else "",
                                      NP=res.get("NP", ""), NPm=res.get("NPm", ""),
                                      h2_status="PASS" if ok else "FAIL:" + res.get("fail", ""))
            wr.writerow(out)
            if not ok or n % 50 == 0:
                fh.flush()
                print("box %d [%s, %s] kappa=%.6f NP=%s NPm=%s %s (%.0fs)" % (n, r["eps_lo"], r["eps_hi"],
                      float(Xm[0].mid()), res.get("NP"), res.get("NPm"), "PASS" if ok else "FAIL " + res.get("fail", ""),
                      time.time() - t0), flush=True)
    last_hi = rows[min(n, len(rows)) - 1]["eps_hi"]
    print("SUMMARY j=%d: %d boxes processed (eps up to %s), %d failures -> %s (%.0fs)" % (
        j, n, last_hi, nf, "PASS" if nf == 0 else "FAIL", time.time() - t0))
    sys.exit(0 if nf == 0 else 1)

if __name__ == "__main__":
    main()
