#!/usr/bin/env python3
"""l1_sweep.py -- certified sign of the first Lyapunov coefficient for ALL M in a range, by eps-boxes (eps = M^-1/2).

Usage:  python l1_sweep.py J EPS_LO EPS_HI OUT.csv [--h2] [W0] [WMAX]
  J        branch (1, 2, 3)
  EPS_LO   lower end (0 allowed: then the box [0, w] contains the limit M = infinity)
  EPS_HI   upper end (e.g. 0.0031623 > 1e5^-1/2, or 0.0788111 > 161^-1/2)
  --h2     also run the H2 checks of h2_lib (j = 1 only) on the same Krawczyk box
Environment: Python 3.13, python-flint 0.9.  1 thread.  Runtime ~0.1 s per box (0.4 s with --h2).

For each box E: parametric Krawczyk box X for the Hopf point (kappa_j, z) valid for all eps in E (allm_core),
simplicity (dQ/dz != 0), transversality, Im z > 0, f(phi) != 0, and the enclosure of
     Lhat = eps^2 l1(q_1 = 1) = l1(q_1 = 1)/M
over E x X by the mean value form (l1reg.Lhat_box).  Claim per box: Lhat < 0 (hence l1 < 0, any normalization).
PASS/FAIL per box in OUT.csv; exit 0 iff all boxes pass and [EPS_LO, EPS_HI] is covered.
"""
import sys, os, csv, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import mpmath as mp
import allm_core as C
import h2_lib as H
import l1reg as R

def limit_point(j):
    mp.mp.dps = 30
    c = (2 * j - mp.mpf(1) / 2) * mp.pi
    t = mp.findroot(lambda u: mp.tan(u) - u, c - 1 / c)
    k0 = mp.sqrt(mp.sqrt(1 + t * t) - 1)
    return [arb(str(k0)), arb(0), arb(str(t))]

FIELDS = ["j", "eps_lo", "eps_hi", "kappa_mid", "x_mid", "y_mid", "rad_kappa", "rad_x", "rad_y", "simple",
          "transv_lo", "fphi_lo", "imz_pos", "Lhat_lo", "Lhat_hi", "h2", "status"]

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_h2 = "--h2" in sys.argv
    j = int(args[0]); elo = args[1]; ehi = float(args[2]); out = args[3]
    w = float(args[4]) if len(args) > 4 else 1e-5
    wmax = float(args[5]) if len(args) > 5 else 2e-4
    ctx.prec = 192
    X0 = limit_point(j)
    if float(elo) > 0:
        X0 = C.newton(X0, arb(0))
        # continuation from eps = 0 to eps_lo in small steps
        e = 0.0
        while e < float(elo):
            e = min(e + 2e-4, float(elo)); X0 = C.newton(X0, arb(e))
    lo = elo
    nbox = nfail = 0; t0 = time.time(); prev = None
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(FIELDS)
        while float(lo) < ehi:
            hi = "%.17g" % min(float(lo) + w, ehi)
            ok = False; info = ""
            if do_h2:
                res = H.certify_box(lo, hi, X0)
                eb = C.eps_interval(lo, hi); box = res["box"]; Xm = res["Xm"]; rad = res["rad"]
                h2ok = res["ok"]
            else:
                eb, Xm, rad, enc, box = H.hopf_box(lo, hi, X0)
                h2ok = None
            if box is not None and rad is not None:
                simple, tr, fphi, laim = C.hopf_checks(box, eb)
                imz = bool(arb(box[2].lower()) > 0)
                L = R.Lhat_box(eb, box[0], box[1], box[2])
                ok = (L is not None and bool(L < 0) and simple and bool(tr > 0) and bool(fphi.abs_lower() > 0)
                      and imz and (h2ok is None or h2ok))
                tight = L is not None and float(L.rad()) < 0.3 * abs(float(L.mid()))
            else:
                L = None; tight = False
            if not ok or not tight:
                if w > 1e-10:
                    w /= 2; continue
                nfail += 1
                print("box [%s, %s] FAIL" % (lo, hi), flush=True)
                if box is None or rad is None:
                    break
            wr.writerow([j, lo, hi, Xm[0].mid().str(30, radius=False), Xm[1].mid().str(30, radius=False),
                         Xm[2].mid().str(30, radius=False), "%.6e" % float(rad[0]), "%.6e" % float(rad[1]),
                         "%.6e" % float(rad[2]), int(simple), "%.6g" % float(arb(tr.lower())),
                         "%.6g" % float(fphi.abs_lower()), int(imz), "%.12g" % float(L.lower()), "%.12g" % float(L.upper()),
                         "" if h2ok is None else int(h2ok), "PASS" if ok else "FAIL"])
            nbox += 1
            if nbox % 50 == 0 or not ok:
                fh.flush()
                print("box %d [%s, %s] kappa=%.8f Lhat in [%.6g, %.6g] %s (%.0fs)" % (
                    nbox, lo, hi, float(Xm[0].mid()), float(L.lower()), float(L.upper()), "PASS" if ok else "FAIL",
                    time.time() - t0), flush=True)
            if prev is not None:
                X0 = [arb(2 * float(Xm[i].mid()) - float(prev[i].mid())) for i in range(3)]
            else:
                X0 = [arb(Xm[i].mid()) for i in range(3)]
            prev = Xm
            lo = hi
            w = min(w * 1.3, wmax)
    cover = float(lo) >= ehi
    okall = nfail == 0 and cover
    print("SUMMARY j=%d eps in [%s, %.10g]: %d boxes, %d failures, covered: %s -> %s (%.0fs)" % (
        j, elo, ehi, nbox, nfail, cover, "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

if __name__ == "__main__":
    main()
