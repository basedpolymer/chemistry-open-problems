#!/usr/bin/env python3
"""Independent check of the branch-and-bound certificate.

Checks (PASS/FAIL):
  V1  the leaves are dyadic boxes of the chart [0,2]x[0,1]^2 with L1 = L2 + 1 = L3 + 1, pairwise
      non-nested (hence with disjoint interiors), and their total volume is exactly 2 (so they tile the chart);
  V2  recomputing every leaf (or a random sample with --sample K) with the certifier reproduces the stored
      counts outR, certR, outL, certL exactly, with zero 'unknot' anomalies;
  V3  the exact dyadic sums U = sum w(D)(outR+outL), L = sum w(D)(certR+certL) (formula (3) of result.md)
      equal the values in the summary file; prints U and L to 12 digits and as fractions num/2^EXP;
  V4  L < Hake's bound: U < (14 - 3 pi)/192 (checked with exact rationals and a rational upper bound of pi... the
      comparison uses pi > 3.14159265358979, i.e. (14-3pi)/192 < (14 - 3*3.14159265358979)/192).
Environment: Python 3.13; certifier binary built by scripts/build.sh.  Runtime: full recomputation ~ the
original run time (hours, one thread); --sample 2000 takes a few minutes.
Usage: python3 verify_certificate.py LEAVES.gz SUMMARY.json HEXCERT TRIG_TABLE [--sample K] [--seed S]
"""
import sys, gzip, json, random, argparse, subprocess
from fractions import Fraction
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from bnb_driver import weight_num, EXP

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("leaves"); ap.add_argument("summary"); ap.add_argument("bin"); ap.add_argument("table")
    ap.add_argument("--sample", type=int, default=0); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rows = [tuple(map(int, line.split())) for line in gzip.open(a.leaves, "rt")]
    summ = json.load(open(a.summary))
    ok = True
    def chk(name, cond, info=""):
        nonlocal ok
        print(("PASS" if cond else "FAIL"), name, info, flush=True); ok = ok and cond
    # V1 tiling
    keys = set()
    good = True
    vol = 0
    for r in rows:
        L1, a1, L2, a2, L3, a3 = r[:6]
        if not (L1 == L2 + 1 == L3 + 1 and 0 <= a1 < 2 ** L1 and 0 <= a2 < 2 ** L2 and 0 <= a3 < 2 ** L3):
            good = False
        keys.add((L2, a1, a2, a3))
        vol += Fraction(2, 2 ** (L1 + L2 + L3))
    for (L, a1, a2, a3) in keys:
        b1, b2, b3, l = a1, a2, a3, L
        while l > 0:
            b1, b2, b3, l = b1 // 2, b2 // 2, b3 // 2, l - 1
            if (l, b1, b2, b3) in keys:
                good = False; break
    chk("V1 leaves are non-nested dyadic boxes", good and len(keys) == len(rows), "(%d leaves)" % len(rows))
    chk("V1 total chart volume = 2", vol == 2, "(%s)" % vol)
    # V2 recomputation
    idx = list(range(len(rows)))
    if a.sample:
        random.Random(a.seed).shuffle(idx); idx = sorted(idx[:a.sample])
    p = subprocess.Popen([a.bin, a.table, "1"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
    mism = 0; anom = 0
    for k in range(0, len(idx), 200):
        part = idx[k:k + 200]
        for i in part:
            p.stdin.write("%d %d %d %d %d %d %d %d\n" % ((i,) + rows[i][:7]))
        p.stdin.flush()
        for i in part:
            f = list(map(int, p.stdout.readline().split()))
            assert f[0] == i
            Lth, outR, inR, certR, outL, inL, certL = f[1:8]
            unk = f[12] + f[13]
            if (outR, certR, outL, certL) != rows[i][7:11] or Lth != rows[i][6]:
                mism += 1
            anom += unk
    p.stdin.close(); p.wait()
    chk("V2 recomputed counts match", mism == 0, "(%d leaves recomputed, %d mismatches)" % (len(idx), mism))
    chk("V2 no component certified unknotted", anom == 0)
    # V3 sums
    up = sum(weight_num(*r[:7]) * (r[7] + r[9]) for r in rows)
    lo = sum(weight_num(*r[:7]) * (r[8] + r[10]) for r in rows)
    chk("V3 exact sums equal the summary", str(up) == summ["upper_num"] and str(lo) == summ["lower_num"])
    U = Fraction(up, 1 << EXP); L = Fraction(lo, 1 << EXP)
    print("   upper U = %.12e   lower L = %.12e   U/L = %.6f" % (float(U), float(L), float(U / L) if L else float("inf")))
    hake_upper_rational = (Fraction(14) - 3 * Fraction("3.14159265358979")) / 192   # > (14-3pi)/192
    chk("V4 0 < L <= U < (14-3pi)/192 (Hake)", 0 < L <= U and U < (Fraction(14) - 3 * Fraction("3.14159265358980")) / 192)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
