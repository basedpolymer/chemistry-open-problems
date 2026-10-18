#!/usr/bin/env python3
"""03_certify_range.py -- certified sign of the first Lyapunov coefficient on a range of M, for one branch j.

Usage:  python 03_certify_range.py J MMIN MMAX OUT.csv
Environment: WSL Fedora, Python 3.12, python-flint 0.9 (Arb ball arithmetic), numpy/scipy for starting guesses.
Runtime: about 15-30 ms per (M, j) for M <= 1e4 (plus a dense float scan for M <= 60).

For every integer M in [MMIN, MMAX]:
  1. a starting guess for the j-th Hopf point (kappa, z) is obtained (dense float scan for M <= 60,
     continuation in M afterwards);
  2. Newton + Krawczyk in ball arithmetic proves that the scaled Hopf system has a unique zero in a small box
     and encloses it;
  3. at the enclosure, the first Lyapunov coefficient l1 is enclosed with the closed-form formula
     (q normalized by q^H q = 1) -- claim: l1 < 0;
  4. nondegeneracy: 0 not in dQ/dz (simple eigenvalue), Re dlambda/dkappa > 0 (transversality).
  5. (non-rigorous identification aid) for M <= 300, and every 50th M up to 2000, a float eigenvalue count
     at kappa_j (1 -/+ 1e-7) must give j-1 / j unstable pairs (the point is the j-th crossing).
Each line of OUT.csv records the certified data; PASS iff all claims hold for all M.
"""
import sys, time, csv, math
import numpy as np
from flint import arb, acb, ctx
sys.path.insert(0, __import__("os").path.dirname(__file__))
import bdlib as bd

def main():
    j, Mmin, Mmax, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    ctx.prec = 192
    hist = []
    if len(sys.argv) > 5:          # optional seed: previous sweep file, its last two points start the continuation
        import gzip
        op = gzip.open if sys.argv[5].endswith(".gz") else open
        with op(sys.argv[5], "rt") as fh:
            prev = [r for r in csv.DictReader(fh) if r["status"] in ("PASS", "FAIL")][-2:]
        hist = [(int(r["M"]), [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]) for r in prev]
        assert hist[-1][0] == Mmin - 1, "seed file must end at MMIN-1"
    nfail = 0
    t0 = time.time()
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["M", "j", "kappa_mid", "x_mid", "y_mid", "box_rad_kappa", "box_rad_x", "box_rad_y",
                     "omega_lo", "omega_hi", "l1_unit_lo", "l1_unit_hi", "l1_q1_lo", "l1_q1_hi",
                     "simple", "transv_lo", "float_id_check", "status"])
        for M in range(Mmin, Mmax + 1):
            res = None
            use_dense = M <= 60 or len(hist) < 2
            if not use_dense:
                (M1, X1), (M2, X2) = hist[-2], hist[-1]
                if abs(float(X2[0]) - float(X1[0])) > 0.02*abs(float(X2[0])):
                    use_dense = True
                else:
                    X0 = [arb(2*float(X2[i]) - float(X1[i])) for i in range(3)]
                    res = bd.certify_point(X0, M, prec=192)
                    if res is not None and abs(float(res["kappa"].mid()) - float(X0[0])) > 0.01*abs(float(X0[0])):
                        res = None          # converged far from the predictor: redo from a dense scan
                    if res is None:
                        use_dense = True
            if use_dense:
                cr = [c for c in bd.dense_crossings(M) if c[3] == j and c[2] == j - 1]
                if not cr:
                    wr.writerow([M, j] + [""]*15 + ["NO_CROSSING"]); fh.flush()
                    print("M=%d j=%d: no crossing found (float scan K<=1e5)" % (M, j), flush=True)
                    hist = []
                    continue
                a, b, _, _ = cr[0]
                K, lam, phi = bd.dense_hopf(M, a, b, j - 1)
                z = M*(phi - 1)
                X0 = [arb(K*math.sqrt(M)), arb(z.real), arb(z.imag)]
                res = bd.certify_point(X0, M, prec=192)
            if res is None:
                nfail += 1
                wr.writerow([M, j] + [""]*15 + ["KRAWCZYK_FAIL"]); fh.flush()
                print("M=%d j=%d: KRAWCZYK FAIL" % (M, j), flush=True)
                hist = []
                continue
            kap, x, y = res["kappa"], res["x"], res["y"]
            # window sanity check on Im z (branch j lies in ((2j-1)pi, 2j pi) in all our data)
            ywin = (2*j - 1)*math.pi < float(y.mid()) < 2*j*math.pi
            fid = ""
            if M <= 300 or (M <= 2000 and M % 50 == 0):
                Kc = float(kap.mid())/math.sqrt(M)
                c_lo = bd.unstable_pairs(Kc*(1 - 1e-7), M); c_hi = bd.unstable_pairs(Kc*(1 + 1e-7), M)
                fid = "ok" if (c_lo == j - 1 and c_hi == j) else "MISMATCH(%d,%d)" % (c_lo, c_hi)
            l1 = res["l1"]; l1q = res["l1q1"]
            ok = (l1 < 0) and (l1q < 0) and res["simple"] and (res["transv"] > 0) and ywin and (fid in ("", "ok"))
            if not ok:
                nfail += 1
            lo = lambda a: float(a.lower()); hi = lambda a: float(a.upper())
            wr.writerow([M, j, kap.mid().str(40, radius=False), x.mid().str(40, radius=False),
                         y.mid().str(40, radius=False),
                         "%.3e" % float(res["box"][0].rad()), "%.3e" % float(res["box"][1].rad()),
                         "%.3e" % float(res["box"][2].rad()),
                         "%.17g" % lo(res["omega"]), "%.17g" % hi(res["omega"]),
                         "%.17g" % lo(l1), "%.17g" % hi(l1), "%.17g" % lo(l1q), "%.17g" % hi(l1q),
                         int(res["simple"]), "%.6g" % lo(res["transv"]), fid, "PASS" if ok else "FAIL"])
            if M % 500 == 0 or M <= Mmin + 2 or not ok:
                fh.flush()
                print("M=%d j=%d kappa=%.10f z=%.6f%+.6fi l1=%.6e l1(q1=1)/M=%.6e simple=%s transv>%.3g id=%s %s (%.0fs)"
                      % (M, j, float(kap.mid()), float(x.mid()), float(y.mid()), float(l1.mid()),
                         float(l1q.mid())/M, res["simple"], lo(res["transv"]), fid, "PASS" if ok else "FAIL",
                         time.time() - t0), flush=True)
            hist.append((M, [kap.mid(), x.mid(), y.mid()]))
            hist = hist[-2:]
    print("SUMMARY j=%d M in [%d,%d]: %s (%d failures) in %.0fs" % (j, Mmin, Mmax, "PASS" if nfail == 0 else "FAIL",
                                                                    nfail, time.time() - t0))
    sys.exit(0 if nfail == 0 else 1)

if __name__ == "__main__":
    main()
