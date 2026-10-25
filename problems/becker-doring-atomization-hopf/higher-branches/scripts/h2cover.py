#!/usr/bin/env python3
"""h2cover.py -- spectral certificate (hypothesis H2 and the unstable count) at the j-th Hopf point for every eps
covered by an l1_sweep certificate file, on merged parameter boxes.

Usage:  python h2cover.py J L1.csv OUT.csv R0 Y cL cR [KREL [cL2 cR2]]
  L1.csv : output of l1_sweep.py for branch J (consecutive eps-boxes E_i with Krawczyk boxes X_i)
  R0, Y, cL, cR : contour parameters of h2j_lib (all-M Theorem 4.5)
  KREL   : maximal relative width of the merged kappa range (default 2e-3)
  cL2 cR2: second pair of level values, tried when the first fails (narrow chains, needed near eps = 0 for
           large J, where the zeros next to the Hopf zero are at distance ~ log(t_{J+1}/t_J) from the axis)
Environment: Python 3.13, python-flint 0.9.  One thread.  Runtime about 5-10 s per merged box.

Method.  Every Krawczyk box X_i of L1.csv is first re-proved (allm_core.krawczyk, all eps in E_i): it contains the
unique zero (kappa_j(eps), z_H(eps)) for eps in E_i.  Consecutive boxes are merged into groups G = {i0..i1};
for a group, E_G = union of the E_i (an interval) and the parameter box  kappa in K_G = hull of the kappa-ranges of
the X_i, z_H in Z_G = hull of their z-ranges.  The checks of h2j_lib.h2j_checks (localization, monotonicity,
region, left/right chains, winding counts N(P) = 2J + 2 and N(P-) = 4) are valid for ALL (eps, kappa) in
E_G x K_G, in particular at (eps, kappa_j(eps)) for every eps in E_G; the containment |Im z_H| < Y is checked on Z_G.
By all-M Theorem 4.5 this gives, for every M with M^-1/2 in E_G: at kappa_j(M) the only eigenvalues of B on iR
are 0, +-i omega, and exactly 2J - 2 zeros of F (with multiplicity) give eigenvalues with Re lambda > 0.
If a group fails it is split in two (down to single boxes).  Claim: every group PASSes and the groups cover the
eps-range of L1.csv.  Exit 0 iff so.
"""
import sys, os, csv, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import allm_core as C
import h2j_lib as HJ

def main():
    j = int(sys.argv[1]); fin, fout = sys.argv[2], sys.argv[3]
    R0, Y, cL, cR = (float(a) for a in sys.argv[4:8])
    krel = float(sys.argv[8]) if len(sys.argv) > 8 else 2e-3
    levels = [(cL, cR)]
    if len(sys.argv) > 10:
        levels.append((float(sys.argv[9]), float(sys.argv[10])))
    ctx.prec = 192
    rows = list(csv.DictReader(open(fin)))
    t0 = time.time()
    cover = rows[0]["eps_lo"] == "0"
    for a, b in zip(rows[:-1], rows[1:]):
        cover &= a["eps_hi"] == b["eps_lo"]
    # 1) re-prove the Krawczyk boxes
    boxes = []
    nbad = 0
    for r in rows:
        eb = C.eps_interval(r["eps_lo"], r["eps_hi"])
        Xm = [arb(r["kappa_mid"]), arb(r["x_mid"]), arb(r["y_mid"])]
        rad = [arb(r["rad_kappa"]), arb(r["rad_x"]), arb(r["rad_y"])]
        enc, box = C.krawczyk(Xm, rad, eb)
        if enc is None:
            nbad += 1
        boxes.append((r["eps_lo"], r["eps_hi"], box, enc is not None))
    print("re-proved %d Krawczyk boxes, %d failures (%.0fs)" % (len(boxes), nbad, time.time() - t0), flush=True)
    # 2) groups
    groups = []
    i = 0
    while i < len(boxes):
        k_lo = float(boxes[i][2][0].lower()); k_hi = float(boxes[i][2][0].upper())
        i1 = i
        while i1 + 1 < len(boxes):
            lo2 = min(k_lo, float(boxes[i1 + 1][2][0].lower())); hi2 = max(k_hi, float(boxes[i1 + 1][2][0].upper()))
            if hi2 - lo2 > krel * hi2:
                break
            k_lo, k_hi = lo2, hi2; i1 += 1
        groups.append((i, i1))
        i = i1 + 1
    print("%d groups" % len(groups), flush=True)
    nf = 0; ndone = 0
    with open(fout, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["j", "eps_lo", "eps_hi", "first_line", "last_line", "kappa_lo", "kappa_hi", "R0", "Y", "cL", "cR",
                     "loc_mode", "mono_ub", "region_ub", "NP", "NPm", "status"])
        stack = list(reversed(groups))
        while stack:
            a, b = stack.pop()
            if not all(boxes[k][3] for k in range(a, b + 1)):
                res = dict(ok=False, fail="krawczyk")
            else:
                eb = C.eps_interval(boxes[a][0], boxes[b][1])
                kap = boxes[a][2][0]; xr = boxes[a][2][1]; yr = boxes[a][2][2]
                for k in range(a + 1, b + 1):
                    kap = kap.union(boxes[k][2][0]); xr = xr.union(boxes[k][2][1]); yr = yr.union(boxes[k][2][2])
                for (cl, cr) in levels:
                    res = HJ.h2j_checks(j, eb, [kap, xr, yr], R0, Y, cl, cr)
                    if res.get("ok", False):
                        break
            if not res.get("ok", False) and b > a:
                m = (a + b) // 2
                stack.append((m + 1, b)); stack.append((a, m))
                continue
            ok = res.get("ok", False)
            nf += not ok; ndone += 1
            wr.writerow([j, boxes[a][0], boxes[b][1], a + 1, b + 1,
                         "%.12g" % float(kap.lower()) if ok else "", "%.12g" % float(kap.upper()) if ok else "",
                         R0, Y, res.get("cL", cL), res.get("cR", cR), res.get("loc_mode", ""), "%.4g" % res["mono_ub"] if "mono_ub" in res else "",
                         "%.4g" % res["region_ub"] if "region_ub" in res else "", res.get("NP", ""), res.get("NPm", ""),
                         "PASS" if ok else "FAIL:" + res.get("fail", "")])
            fh.flush()
            if not ok or ndone % 20 == 0:
                print("group %d [%s, %s] lines %d-%d NP=%s NPm=%s %s (%.0fs)" % (ndone, boxes[a][0], boxes[b][1], a + 1,
                      b + 1, res.get("NP"), res.get("NPm"), "PASS" if ok else "FAIL " + res.get("fail", ""),
                      time.time() - t0), flush=True)
    okall = nf == 0 and nbad == 0 and cover
    print("SUMMARY j=%d: %d groups covering eps in [%s, %s] (contiguous: %s), %d failures -> %s (%.0fs)" % (
        j, ndone, rows[0]["eps_lo"], rows[-1]["eps_hi"], cover, nf, "PASS" if okall else "FAIL", time.time() - t0))
    sys.exit(0 if okall else 1)

if __name__ == "__main__":
    main()
