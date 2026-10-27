"""Exploratory: ALL-certified networks having an OUT (or ESC) siphon facet (Conjecture F), smallest first."""
import os
import sys
import gzip
import pickle
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d2cert  # noqa: E402
import d2lib  # noqa: E402

out = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/try1.pkl.gz")))
C = Counter()
shown = 0
for n, rx, v, r in out:
    if v != "ALL":
        continue
    N = d2lib.Net(n, rx)
    F = d2cert.facets(N)
    ts = sorted({f["type"] for f in F if f["siphon"]})
    C[tuple(ts)] += 1
    if "OUT" in ts and shown < 4:
        shown += 1
        print("ALL with OUT facet:", d2lib.describe(n, rx))
        print("   dual:", d2lib.describe(N.r, N.dual_rx(), [f"k{j+1}" for j in range(N.r)]))
        print("   facets:", [(f["Z"], f["type"], f["fl"], f["mu"]) for f in F if f["siphon"]])
        print("   vertices:", r["rv"])
print(C)
