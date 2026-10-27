"""Exploratory: consistency of d2cert verdicts with E2's certified ALL / NONE networks (no claim)."""
import os
import sys
import gzip
import pickle
import time
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d2cert  # noqa: E402
import d2lib  # noqa: E402


def work(q):
    r = d2cert.analyse(q["n"], q["rx"])
    return (q["cat"], r["verdict"], q["n"], q["rx"])


if __name__ == "__main__":
    t0 = time.time()
    recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
    dec = [q for q in recs if q["cat"] != "UNDECIDED"]
    with Pool(11) as P:
        out = P.map(work, dec, chunksize=10)
    C = Counter((a, b) for a, b, n, rx in out)
    print(sorted(C.items()), round(time.time() - t0, 1))
    for a, b, n, rx in out:
        if (a == "ALL" and b == "NONE") or (a.startswith("NONE") and b == "ALL"):
            print("CONFLICT", a, b, d2lib.describe(n, rx))
