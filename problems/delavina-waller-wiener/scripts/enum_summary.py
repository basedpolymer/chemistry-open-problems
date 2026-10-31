#!/usr/bin/env python3
"""Aggregate the geng/enum_stream outputs in ~/dw into a table of max W by (n, diameter) and check the
conjecture for n = 7, 9, 11. PASS/FAIL printed; exit non-zero on FAIL.  Needs out_N.txt from
run_enum_small.sh and out_11_{0..3}.txt from run_enum_11.sh. Expected total for n=11: 1006700565
connected graphs (OEIS A001349)."""
import glob, os, re, sys
from collections import defaultdict

base = os.path.expanduser("~/dw")
hist = defaultdict(lambda: defaultdict(int))
sample = defaultdict(dict)
totals = defaultdict(int)
for f in sorted(glob.glob(os.path.join(base, "out_*.txt"))):
    n = None
    d = None
    for line in open(f):
        m = re.match(r"n=(\d+) total_connected_read=(\d+)", line)
        if m:
            n = int(m.group(1))
            totals[n] += int(m.group(2))
            continue
        m = re.match(r"diam=(\d+) count=(\d+) maxW=(\d+)", line)
        if m:
            d = int(m.group(1))
            continue
        m = re.match(r"\s+W=(\d+) count=(\d+) :\s*(.*)", line)
        if m:
            w, c = int(m.group(1)), int(m.group(2))
            hist[(n, d)][w] += c
            sample[(n, d, w)] = m.group(3).split()[:4]

expected = {5: 21, 6: 112, 7: 853, 8: 11117, 9: 261080, 10: 11716571, 11: 1006700565}
fail = 0
for n in sorted(totals):
    ok = totals[n] == expected.get(n)
    print(("PASS" if ok else "FAIL"), f"n={n}: {totals[n]} connected graphs read, expected {expected.get(n)}")
    fail += not ok
print("max W by (n, diameter), with the three largest values and the number of graphs attaining them")
for key in sorted(hist):
    tops = sorted(hist[key], reverse=True)[:3]
    print(key, [(w, hist[key][w]) for w in tops])
for d in (3, 4, 5):
    n = 2 * d + 1
    Wc = n * d * (d + 1) // 2
    top = max(max(hist[(n, dd)]) for dd in range(1, d + 1) if (n, dd) in hist)
    ok = top == Wc
    print(("PASS" if ok else "FAIL"), f"d={d}: max W over connected graphs of order {n} and diameter <= {d} is {top}, W(C_{n})={Wc}; "
          f"graphs of diameter {d} attaining it: {hist[(n, d)][Wc]} {sample.get((n, d, Wc))}")
    fail += not ok
sys.exit(1 if fail else 0)
