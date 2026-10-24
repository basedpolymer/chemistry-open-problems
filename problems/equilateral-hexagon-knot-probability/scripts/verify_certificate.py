#!/usr/bin/env python3
"""Independent check of the branch-and-bound certificate.

Checks (PASS/FAIL):
  V1  the leaves are dyadic boxes of the chart [0,2]x[0,1]^2 forming a k-d tree partition (recursive bisection
      check) and their total volume is exactly 2 (so they tile the chart);
  V2  recomputing every leaf (or a random sample with --sample K) with the certifier reproduces the stored
      counts outR, certR, outL, certL exactly, with zero 'unknot' anomalies;
  V3  the exact dyadic sums U = sum w(D)(outR+outL), L = sum w(D)(certR+certL) (formula (3) of result.md)
      equal the values in the summary file; prints U and L to 12 digits and as fractions num/2^EXP;
  V4  0 < L <= U < (14 - 3 pi)/192 (Hake's bound), checked with exact rationals using pi < 3.1415926535898,
      so that (14 - 3*3.1415926535898)/192 < (14 - 3 pi)/192.
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
    ap.add_argument("--only", default="", help="comma-separated leaf indices to recompute (testing)")
    ap.add_argument("--workers", type=int, default=1, help="parallel certifier processes")
    a = ap.parse_args()
    rows = [tuple(map(int, line.split())) for line in gzip.open(a.leaves, "rt")]
    summ = json.load(open(a.summary))
    ok = True
    def chk(name, cond, info=""):
        nonlocal ok
        print(("PASS" if cond else "FAIL"), name, info, flush=True); ok = ok and cond
    # V1 tiling: recursive bisection check (every node of the k-d tree is either a leaf or is bisected along an
    # axis on which all leaves inside it are finer, and both halves are tiled), plus exact volume = 2
    sys.setrecursionlimit(10000)
    good = all(0 <= r[1] < 2 ** r[0] and 0 <= r[3] < 2 ** r[2] and 0 <= r[5] < 2 ** r[4] for r in rows)
    boxes = [r[:6] for r in rows]
    def tiled(node, items):
        if len(items) == 1 and boxes[items[0]] == node:
            return True
        if not items:
            return False
        for k in range(3):
            Lk, ak = node[2 * k], node[2 * k + 1]
            if all(boxes[i][2 * k] > Lk for i in items):
                halves = ([], [])
                for i in items:
                    Li, ai = boxes[i][2 * k], boxes[i][2 * k + 1]
                    halves[(ai >> (Li - Lk - 1)) - 2 * ak].append(i)
                for h in range(2):
                    child = list(node); child[2 * k] = Lk + 1; child[2 * k + 1] = 2 * ak + h
                    if not tiled(tuple(child), halves[h]):
                        return False
                return True
        return False
    good = good and tiled((0, 0, 0, 0, 0, 0), list(range(len(rows))))
    vol = sum(Fraction(2, 2 ** (r[0] + r[2] + r[4])) for r in rows)
    chk("V1 leaves tile the chart [0,2]x[0,1]^2 (k-d tree check)", good, "(%d leaves)" % len(rows))
    chk("V1 total chart volume = 2", vol == 2, "(%s)" % vol)
    # V2 recomputation
    idx = list(range(len(rows)))
    if a.sample:
        random.Random(a.seed).shuffle(idx); idx = sorted(idx[:a.sample])
    if a.only:
        idx = sorted(int(x) for x in a.only.split(","))
    from concurrent.futures import ThreadPoolExecutor
    def work(sub):
        p = subprocess.Popen([a.bin, a.table, "1"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
        mm = 0; an = 0
        for k in range(0, len(sub), 100):
            part = sub[k:k + 100]
            for i in part:
                p.stdin.write("%d %d %d %d %d %d %d %d\n" % ((i,) + tuple(rows[i][:7])))
            p.stdin.flush()
            for i in part:
                f = list(map(int, p.stdout.readline().split()[:14]))
                assert f[0] == i
                Lth, outR, inR, certR, outL, inL, certL = f[1:8]
                if (outR, certR, outL, certL) != tuple(rows[i][7:11]) or Lth != rows[i][6]:
                    mm += 1
                an += f[12] + f[13]
        p.stdin.close(); p.wait()
        return mm, an
    K = max(1, a.workers)
    parts = [idx[j::K] for j in range(K)]
    with ThreadPoolExecutor(K) as ex:
        res = list(ex.map(work, parts))
    mism = sum(r[0] for r in res); anom = sum(r[1] for r in res)
    chk("V2 recomputed counts match", mism == 0, "(%d leaves recomputed, %d mismatches)" % (len(idx), mism))
    chk("V2 no component certified unknotted", anom == 0)
    # V3 sums
    up = sum(weight_num(*r[:7]) * (r[7] + r[9]) for r in rows)
    lo = sum(weight_num(*r[:7]) * (r[8] + r[10]) for r in rows)
    chk("V3 exact sums equal the summary", summ["EXP"] == EXP and str(up) == summ["upper_num"] and str(lo) == summ["lower_num"])
    U = Fraction(up, 1 << EXP); L = Fraction(lo, 1 << EXP)
    print("   upper U = %.12e   lower L = %.12e   U/L = %.6f" % (float(U), float(L), float(U / L) if L else float("inf")))
    chk("V4 0 < L <= U < (14-3pi)/192 (Hake)", 0 < L <= U and U < (Fraction(14) - 3 * Fraction("3.14159265358980")) / 192)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
