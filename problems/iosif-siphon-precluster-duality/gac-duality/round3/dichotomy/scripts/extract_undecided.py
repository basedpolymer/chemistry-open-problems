"""Rebuild E2's s = 2 classification (round2/computation/scripts/s2_persistence.py, imported unchanged) and save the
18 861 WR networks with (H1)-(H3) and s = 2, with E2's verdict, to ~/gacD2-data/s2_all.pkl.gz.

Environment: Python 3.12 (WSL Fedora 44), standard library + E2's modules.  Exact (no numerics).
Usage (from round3/dichotomy/scripts/): python extract_undecided.py [--nproc P]
"""
import os
import sys
import gzip
import time
import pickle
import random
import platform
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
E2DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round2", "computation", "scripts"))
sys.path.insert(0, E2DIR)
import e2lib  # noqa: E402
import s2_persistence as s2p  # noqa: E402

OUT = os.path.expanduser("~/gacD2-data")


def main():
    t0 = time.time()
    argv = sys.argv[1:]
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 8
    print(f"# extract_undecided.py  Python {platform.python_version()} on {platform.system()}; exact")
    rng = random.Random(22)  # same seed and order as s2_persistence.main
    jobs = []
    for fn in e2lib.E_FILES:
        recs = [rc for rc in e2lib.load([fn]) if rc["s"] == 2 and e2lib.H(rc) and rc["N"]["WR"]]
        rng.shuffle(recs)
        for rc in recs:
            jobs.append((rc["n"], rc["rx"], not rc["D"]["CS"], rng.randint(0, 10 ** 9),
                         rc["N"].get("deficiency"), rc["N"].get("lc", None)))
    seen = set()
    uj = []
    for j in jobs:
        if (j[0], j[1]) not in seen:
            seen.add((j[0], j[1]))
            uj.append(j)
    uj.sort(key=lambda j: (j[0] + len(j[1]), len(j[1]), j[1]))
    print(f"# {len(uj)} WR networks with (H1)-(H3) and s = 2")
    with Pool(nproc) as pool:
        res = pool.map(s2p.classify, [(n, rx, ce, sd, False) for n, rx, ce, sd, _, _ in uj], chunksize=20)
    for q, j in zip(res, uj):
        q["deficiency"] = j[4]
    cnt = {}
    for q in res:
        cnt[q["cat"]] = cnt.get(q["cat"], 0) + 1
    print(f"# E2 categories: {dict(sorted(cnt.items()))}")
    ok = cnt.get("UNDECIDED") == 14825 and len(res) == 18861
    print(f"{'PASS' if ok else 'FAIL'}: reproduces E2's counts (18 861 networks, 14 825 undecided)")
    os.makedirs(OUT, exist_ok=True)
    with gzip.open(os.path.join(OUT, "s2_all.pkl.gz"), "wb") as fh:
        pickle.dump(res, fh)
    print(f"# saved ~/gacD2-data/s2_all.pkl.gz; runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
