"""Exploratory run of d2cert.analyse on E2's undecided networks (no claim)."""
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
    try:
        r = d2cert.analyse(q["n"], q["rx"])
        r.pop("F")
        r.pop("V")
    except Exception as e:  # noqa: BLE001
        return (q["n"], q["rx"], "ERR " + repr(e), None)
    return (q["n"], q["rx"], r["verdict"], r)


if __name__ == "__main__":
    t0 = time.time()
    recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
    und = [q for q in recs if q["cat"] == "UNDECIDED"]
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else len(und)
    with Pool(11) as P:
        out = P.map(work, und[:lim], chunksize=10)
    C = Counter(o[2] if not o[2].startswith("ERR") else "ERR" for o in out)
    print(C, round(time.time() - t0, 1))
    for o in out:
        if o[2].startswith("ERR"):
            print(d2lib.describe(o[0], o[1]), o[2])
            break
    with gzip.open(os.path.expanduser("~/gacD2-data/try1.pkl.gz"), "wb") as fh:
        pickle.dump(out, fh)
    sm = {}
    for o in out:
        sm.setdefault(o[2], o)
    for k, o in sm.items():
        print(k, d2lib.describe(o[0], o[1]))
        if o[3]:
            print("   ", {a: o[3][a] for a in ("none", "uf", "uv", "cycle")})
