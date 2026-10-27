"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: how many of D2's NONE verdicts rest only on the
bicoloured case of Prop. E (whose written proof has a gap, see verification.md), split by the minimal order mu of
the facet.  Uses D2's d2cert.analyse (imported unchanged, read-only) on ~/gacD2-data/s2_all.pkl.gz.
Usage: python v3_propE_split.py.  Environment: Python 3.12, standard library."""
import os, sys, gzip, pickle, time
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
D2DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round3", "dichotomy", "scripts"))

def work(q):
    sys.path.insert(0, D2DIR)
    import d2cert
    r = d2cert.analyse(q["n"], q["rx"])
    if r["verdict"] != "NONE":
        return None
    kinds = {c[0] for c in r["none"]}
    if kinds != {"E"}:
        return ("other",)
    wit = [c for c in r["none"]]
    if any(c[3][0] == "zero-flow" for c in wit):
        return ("E-zero",)
    mus = []
    for c in wit:
        Z = c[1]
        f = next(f for f in r["F"] if f["Z"] == Z)
        mus.append(f["mu"])
    return ("E-bi-only", min(mus), q["cat"])

t0 = time.time()
recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
with Pool(6) as P:
    out = P.map(work, recs, chunksize=20)
c = Counter(o[0] for o in out if o)
print("NONE verdicts by reliance:", dict(c))
mu = Counter(o[1] for o in out if o and o[0] == "E-bi-only")
print("bicoloured-only, by smallest mu among the witnesses:", dict(sorted(mu.items())))
cat = Counter(o[2] for o in out if o and o[0] == "E-bi-only")
print("bicoloured-only, E2 category:", dict(cat))
print("# runtime %.1f s" % (time.time() - t0))
