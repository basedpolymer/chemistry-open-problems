"""Exploratory: what remains undecided after d2cert (no claim)."""
import os
import sys
import gzip
import pickle
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d2lib  # noqa: E402

out = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/try1.pkl.gz")))
und = [o for o in out if o[2] == "UNDECIDED"]
C = Counter()
ex = {}
for n, rx, v, r in und:
    key = (tuple(sorted({f"{t}/{fl}/mu{mu}" for Z, t, fl, mu in r["uf"]})), len(r["uv"]) > 0, r["cycle"])
    C[key] += 1
    ex.setdefault(key, (n, rx, r))
for k, c in C.most_common(60):
    n, rx, r = ex[k]
    print(c, k, "|", d2lib.describe(n, rx), "| uv", r["uv"], "uf", r["uf"])
