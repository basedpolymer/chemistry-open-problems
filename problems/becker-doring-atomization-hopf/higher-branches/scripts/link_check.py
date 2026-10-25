#!/usr/bin/env python3
"""link_check.py -- the integer points of pointwise.py lie on the certified j-th Hopf curve of onset_sweep.py.

Usage:  python link_check.py PW.csv ONSET.csv
Environment: Python 3.13, python-flint 0.9.  One thread.  Runtime: about 10-60 s.

For every line of PW.csv (integer M, Krawczyk box X_M of result 8's chart: kappa, x, y; unique zero in X_M):
the image of X_M in the onset chart (s = 1/(1 + kappa M^-1/2), u + i v = (x + i y)/M) is computed in ball
arithmetic.  A line of ONSET.csv is selected -- an M-chart line whose M-range contains M, or else an s-chart line
whose s-range contains the image -- and re-proved (onset_sweep.krawczyk on the stored box: for every parameter in
its range exactly one zero in the box); the image must lie in that box.  Then the zero of X_M is the point of the
onset curve at that M.  Claim (PASS/FAIL): every line of PW.csv is linked and the M of PW.csv are contiguous.
"""
import sys, os, csv, time, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import onset_sweep as O

def main():
    pw, on = sys.argv[1], sys.argv[2]
    ctx.prec = 192
    t0 = time.time()
    P = list(csv.DictReader(open(pw)))
    S = list(csv.DictReader(open(on)))
    Mrows = [r for r in S if r["chart"] == "M"]
    srows = [r for r in S if r["chart"] == "s"]
    Mlo = [float(r["p_lo"]) for r in Mrows]
    slo = [float(r["p_lo"]) for r in srows]
    proved = {}
    def prove(row):
        key = (row["chart"], row["p_lo"])
        if key not in proved:
            Pb = C.eps_interval(row["p_lo"], row["p_hi"])
            Xm = [arb(row["X0_mid"]), arb(row["u_mid"]), arb(row["v_mid"])]
            rad = [arb(row["rad_X0"]), arb(row["rad_u"]), arb(row["rad_v"])]
            img, box = O.krawczyk(row["chart"], Xm, rad, Pb)
            proved[key] = (Pb, box if img is not None else None)
        return proved[key]
    nf = 0
    Ms = [int(r["M"]) for r in P]
    contiguous = Ms == list(range(Ms[0], Ms[-1] + 1))
    for r in P:
        M = int(r["M"])
        X = [arb(r["kappa_mid"], r["rad_kappa"]), arb(r["x_mid"], r["rad_x"]), arb(r["y_mid"], r["rad_y"])]
        Mb = arb(M)
        s_img = 1 / (1 + X[0] / Mb.sqrt())
        u_img = X[1] / Mb; v_img = X[2] / Mb
        ok = False; where = ""
        k = bisect.bisect_right(Mlo, M) - 1
        for kk in (k, k - 1, k + 1):
            if 0 <= kk < len(Mrows):
                Pb, box = prove(Mrows[kk])
                if box is not None and arb(Pb.lower()) <= Mb and Mb <= arb(Pb.upper()) and \
                        box[0].contains(s_img) and box[1].contains(u_img) and box[2].contains(v_img):
                    ok = True; where = "M-chart line %d" % kk; break
        if not ok:
            k = bisect.bisect_right(slo, float(s_img.mid())) - 1
            for kk in (k, k - 1, k + 1):
                if 0 <= kk < len(srows):
                    Pb, box = prove(srows[kk])
                    if box is not None and arb(Pb.lower()) <= arb(s_img.lower()) and arb(s_img.upper()) <= arb(Pb.upper()) \
                            and box[0].contains(Mb) and box[1].contains(u_img) and box[2].contains(v_img):
                        ok = True; where = "s-chart line %d" % kk; break
        nf += not ok
        if not ok or M % 25 == 0 or M == Ms[0]:
            print("M=%d: s = %s, linked: %s (%s)" % (M, s_img.str(8), ok, where), flush=True)
    okall = nf == 0 and contiguous
    print("SUMMARY %s: %d integer points (M = %d..%d, contiguous: %s) linked to the curve of %s, %d failures -> %s (%.0fs)"
          % (os.path.basename(pw), len(P), Ms[0], Ms[-1], contiguous, os.path.basename(on), nf,
             "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

if __name__ == "__main__":
    main()
