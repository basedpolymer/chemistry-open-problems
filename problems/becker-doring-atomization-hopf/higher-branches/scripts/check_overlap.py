#!/usr/bin/env python3
"""check_overlap.py -- consistency of the two first-Lyapunov-coefficient implementations on the overlap of the
eps-chart with the integer range: result 8's closed form (bdlib.lyapunov_cf, integer M) and the all-M eps-regular
formula (l1reg, eps-boxes).

Usage:  python check_overlap.py L1.csv MFROM MTO
For every integer M in [MFROM, MTO] (with M^-1/2 inside the eps-range of L1.csv): bdlib's Krawczyk test at M
started from the l1 box encloses the Hopf point; its enclosure must lie inside the l1 Krawczyk box of the eps-box
containing M^-1/2 (same zero, by uniqueness there), and bdlib's l1(q_1 = 1)/M must lie inside that line's
enclosure [Lhat_lo, Lhat_hi].  Claim (PASS/FAIL): all M pass.
Environment: Python 3.13, python-flint 0.9.  Runtime: about 0.1 s per M.
"""
import sys, os, csv, bisect, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import bdlib as bd

def main():
    fn = sys.argv[1]; M0 = int(sys.argv[2]); M1 = int(sys.argv[3])
    ctx.prec = 192
    rows = list(csv.DictReader(open(fn)))
    his = [float(r["eps_hi"]) for r in rows]
    t0 = time.time(); nf = 0; n = 0
    for M in range(M0, M1 + 1):
        e = float(arb(1) / arb(M).sqrt())
        k = bisect.bisect_left(his, e)
        r = rows[k]
        eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
        if not eb.contains(arb(1) / arb(M).sqrt()):
            k += 1; r = rows[k]; eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
        Xm = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
        rad = [arb(r["rad_kappa"]), arb(r["rad_x"]), arb(r["rad_y"])]
        enc_e, box = C.krawczyk(Xm, rad, eb)
        res = bd.certify_point(Xm, M, prec=192)
        ok = enc_e is not None and res is not None
        if ok:
            ok = all(box[i].contains(res["box"][i]) for i in range(3))
            L = res["l1q1"] / M
            ok = ok and bool(arb(r["Lhat_lo"]) <= L) and bool(L <= arb(r["Lhat_hi"]))
        nf += not ok; n += 1
        if not ok or M % 50 == 0:
            print("M=%d: bdlib l1(q1)/M = %s, l1 line [%s, %s]: %s" % (M, (res["l1q1"] / M).str(10) if res else None,
                  r["Lhat_lo"], r["Lhat_hi"], "PASS" if ok else "FAIL"), flush=True)
    print("SUMMARY %s M in [%d, %d]: %d points, %d failures -> %s (%.0fs)" % (os.path.basename(fn), M0, M1, n, nf,
          "PASS" if nf == 0 else "FAIL", time.time() - t0))
    sys.exit(0 if nf == 0 else 1)

if __name__ == "__main__":
    main()
