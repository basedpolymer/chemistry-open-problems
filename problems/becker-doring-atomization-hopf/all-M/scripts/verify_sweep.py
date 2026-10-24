#!/usr/bin/env python3
"""verify_sweep.py -- re-proves, from the stored data alone, every line of an eps-box certificate file.

Usage:  python verify_sweep.py FILE.csv EPS_LO EPS_HI [--h2] [STRIDE]
  FILE.csv  written by l1_sweep.py (columns eps_lo, eps_hi, kappa_mid, x_mid, y_mid, rad_*, Lhat_*, h2, status)
            or by h2_sweep.py (no Lhat columns; then only the Hopf/H2 claims are re-proved)
  --h2      re-prove the H2 claims (localization, region, winding = 4) with R0 = 8, sigma = 0.3, a1 = 0.1
Environment: Python 3.13, python-flint 0.9.  1 thread.  Runtime about the same as the sweep (0.1-0.6 s per line).
For each line (every STRIDE-th line; STRIDE = 1 re-proves all): the box X = mid +- rad is re-checked by the
Krawczyk test uniformly in eps in [eps_lo, eps_hi] (so X contains a unique zero of the Hopf system for every such
eps), simplicity/transversality/Im z > 0/f(phi) != 0, and (if present) Lhat < 0 on E x X by the mean value form.
Also checks coverage: first eps_lo = EPS_LO, consecutive lines share endpoints (identical strings), last
eps_hi >= EPS_HI.  Prints a PASS/FAIL line per checked line (failures, every 100th) and a SUMMARY; exit code 0 iff
all PASS.
"""
import sys, os, csv, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import h2_lib as H
import l1reg as R

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    fn, elo, ehi = args[0], args[1], float(args[2])
    stride = int(args[3]) if len(args) > 3 else 1
    do_h2 = "--h2" in sys.argv
    ctx.prec = 192
    rows = list(csv.DictReader(open(fn)))
    cover = rows[0]["eps_lo"] == elo and float(rows[-1]["eps_hi"]) >= ehi
    for a, b in zip(rows[:-1], rows[1:]):
        cover &= (a["eps_hi"] == b["eps_lo"])
    has_L = "Lhat_lo" in rows[0]
    nf = 0; n = 0; t0 = time.time()
    for k, r in enumerate(rows):
        if k % stride and k != len(rows) - 1:
            continue
        n += 1
        eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
        Xm = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
        rad = [arb(r["rad_kappa"]), arb(r["rad_x"]), arb(r["rad_y"])]
        enc, box = C.krawczyk(Xm, rad, eb)
        ok = enc is not None and r["status"] == "PASS"
        Ls = ""
        if enc is not None:
            simple, tr, fphi, laim = C.hopf_checks(box, eb)
            ok &= simple and bool(tr > 0) and bool(fphi.abs_lower() > 0) and bool(arb(box[2].lower()) > 0)
            if has_L:
                L = R.Lhat_box(eb, box[0], box[1], box[2])
                ok &= L is not None and bool(L < 0)
                Ls = " Lhat in [%.6g, %.6g]" % (float(L.lower()), float(L.upper())) if L is not None else " Lhat None"
            if do_h2:
                h = H.h2_checks(eb, box)
                ok &= h["ok"]
                Ls += " H2 count=%s loc=%s" % (h.get("count"), h.get("loc_mode"))
        nf += (not ok)
        if not ok or n % 100 == 0:
            print("line %d [%s, %s]%s %s (%.0fs)" % (k + 1, r["eps_lo"], r["eps_hi"], Ls, "PASS" if ok else "FAIL",
                                                      time.time() - t0), flush=True)
    okall = nf == 0 and cover
    print("SUMMARY %s: %d lines re-proved (stride %d), %d failures, coverage of [%s, %g]: %s -> %s (%.0fs)" % (
        os.path.basename(fn), n, stride, nf, elo, ehi, cover, "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

if __name__ == "__main__":
    main()
