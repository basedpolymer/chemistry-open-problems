#!/usr/bin/env python3
"""h2_sweep.py -- certified spectral hypothesis H2 at kappa_1(M) for ALL M >= MMIN (branch j = 1).

Usage:  python h2_sweep.py OUT.csv [MMIN=161] [W0=1e-4] [WMAX=4e-4]
Environment: Python 3.13, python-flint 0.9 (Arb).  1 thread.  Runtime: about 0.3-0.5 s per eps-box
(a few hundred boxes cover eps in [0, MMIN^-1/2]).

Covers the closed interval eps in [0, MMIN^-1/2] (eps = M^-1/2) by consecutive boxes; for each box all the checks
of h2_lib.certify_box are made (Krawczyk Hopf box, simplicity, transversality, localization, region, winding
number = 4).  Each line of OUT.csv is one box; status PASS iff all checks hold.  Prints PASS/FAIL per box
(every 25th box and failures) and a final SUMMARY; exit code 0 iff every box passes and the boxes cover
[0, MMIN^-1/2] without gaps.
"""
import sys, os, csv, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import h2_lib as H

def main():
    out = sys.argv[1]
    Mmin = int(sys.argv[2]) if len(sys.argv) > 2 else 161
    w = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
    wmax = float(sys.argv[4]) if len(sys.argv) > 4 else 4e-4
    ctx.prec = 192
    eps_end = 1 / arb(Mmin).sqrt()
    e_end = float(eps_end.upper()) * (1 + 1e-12)
    X0 = [arb("1.8982"), arb("0.0"), arb("4.4934")]
    lo = "0"
    nbox = 0; nfail = 0; t0 = time.time()
    prev = None
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(H.CSV_FIELDS)
        while True:
            flo = float(lo)
            if flo >= e_end:
                break
            fhi = min(flo + w, e_end)
            hi = "%.17g" % fhi
            res = H.certify_box(lo, hi, X0)
            if not res["ok"]:
                if w > 1e-9:
                    w /= 2
                    continue
                nfail += 1
                print("box [%s, %s] FAIL (%s)" % (lo, hi, res.get("fail")), flush=True)
                if res["rad"] is None:
                    break
            wr.writerow(H.row(res)); nbox += 1
            if nbox % 25 == 0 or not res["ok"]:
                fh.flush()
                print("box %d [%s, %s] M in [%.0f, %s] kappa=%.8f count=%s loc=%s rhoI=%.3f %s (%.0fs)" % (
                    nbox, lo, hi, 1 / fhi ** 2, ("%.0f" % (1 / flo ** 2)) if flo > 0 else "inf",
                    float(res["Xm"][0].mid()), res["count"], res["loc_mode"], res["rhoI"],
                    "PASS" if res["ok"] else "FAIL", time.time() - t0), flush=True)
            # continuation
            Xm = res["Xm"]
            if prev is not None:
                X0 = [arb(2 * float(Xm[i].mid()) - float(prev[i].mid())) for i in range(3)]
            else:
                X0 = [arb(Xm[i].mid()) for i in range(3)]
            prev = Xm
            lo = hi
            w = min(w * 1.25, wmax)
    # coverage: boxes are [lo_k, hi_k] with hi_k = lo_{k+1} (same decimal string, rounded outwards), first lo = 0,
    # last hi >= MMIN^-1/2
    cover = float(lo) >= float(eps_end.upper())
    ok = nfail == 0 and cover
    print("SUMMARY: %d boxes, %d failures, coverage of [0, %d^-1/2]: %s -> %s (%.0fs)" % (
        nbox, nfail, Mmin, cover, "PASS" if ok else "FAIL", time.time() - t0))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
