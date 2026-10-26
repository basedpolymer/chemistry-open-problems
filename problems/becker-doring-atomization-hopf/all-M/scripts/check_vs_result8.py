#!/usr/bin/env python3
"""check_vs_result8.py -- the Hopf points certified by result 8 (pointwise in M) are the Hopf points of the eps-box
certificates of this folder.

Usage:  python check_vs_result8.py J OURS.csv R8_FILE1.csv.gz [R8_FILE2.csv.gz ...]
Environment: Python 3.13, python-flint 0.9.  Runtime: about 1 minute per 1e5 lines.
For every line of result 8's sweep files with M >= 161 (status PASS/FAIL, i.e. a certified Hopf box): locate the
eps-box [eps_lo, eps_hi] of OURS.csv containing eps = M^-1/2 (exact comparison eps_lo^2 <= 1/M <= eps_hi^2), and check
that result 8's box (kappa, x, y) = mid +- rad lies inside our Krawczyk box X = mid +- rad (which contains a unique
zero of the same Hopf system for this eps).  Then both certificates concern the same Hopf point, and result 8's
sign of l1 and ours (Lhat < 0 on the whole box) agree.  Also reports the agreement of our Lhat enclosure with result
8's value l1(q1=1)/M.  PASS iff every M is matched.
"""
import sys, os, csv, gzip, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx

def main():
    j = int(sys.argv[1]); ours = list(csv.DictReader(open(sys.argv[2])))
    ctx.prec = 192
    los = [float(r["eps_lo"]) for r in ours]
    n = nf = nL = 0
    seen = set()
    for fn in sys.argv[3:]:
        for r in csv.DictReader(gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)):
            if r["status"] not in ("PASS", "FAIL") or int(r["M"]) < 161 or int(r["j"]) != j:
                continue
            M = int(r["M"]); n += 1
            if M in seen or r["status"] != "PASS":
                nf += 1
            seen.add(M)
            e = 1 / M ** 0.5
            k = bisect.bisect_right(los, e) - 1
            cand = [i for i in (k - 1, k, k + 1) if 0 <= i < len(ours)]
            ok = False
            for i in cand:
                b = ours[i]
                if not (arb(b["eps_lo"]) ** 2 * M <= 1 and arb(b["eps_hi"]) ** 2 * M >= 1):
                    continue
                inside = True
                for key, rk in (("kappa", "rad_kappa"), ("x", "rad_x"), ("y", "rad_y")):
                    ours_box = arb(b[key + "_mid"], b[rk])
                    r8 = arb(r[key + "_mid"], r["box_rad_" + key])
                    inside &= bool(ours_box.contains(r8))
                if inside:
                    ok = True
                    if "Lhat_lo" in b and r["l1_q1_lo"]:
                        v = (float(r["l1_q1_lo"]) + float(r["l1_q1_hi"])) / 2 / M
                        nL += (float(b["Lhat_lo"]) <= v <= float(b["Lhat_hi"]))
                    break
            nf += (not ok)
            if not ok:
                print("M=%d: not matched" % M, flush=True)
    nf += seen != set(range(161, 100001))
    print("SUMMARY j=%d: %d result-8 Hopf points with M >= 161, %d matched inside our Krawczyk boxes, "
          "%d with result 8's l1(q1=1)/M inside our Lhat enclosure -> %s" % (j, n, n - nf, nL, "PASS" if nf == 0 else "FAIL"))
    sys.exit(0 if nf == 0 else 1)

if __name__ == "__main__":
    main()
