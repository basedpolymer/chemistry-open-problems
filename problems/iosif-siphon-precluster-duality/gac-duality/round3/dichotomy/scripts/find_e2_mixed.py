"""Identify E2's two 'certificate for some sampled kappa only' networks (exploratory)."""
import os, sys, gzip, pickle, random
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
E2DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "round2", "computation", "scripts"))
sys.path.insert(0, E2DIR); sys.path.insert(0, HERE)
import e2lib, s2_persistence as s2p, d2lib
recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
und = [q for q in recs if q["cat"] == "UNDECIDED"][:200]
for q in und:
    n, rx = q["n"], q["rx"]
    p, d = len(rx), e2lib.dual_rx(n, rx)
    rng = random.Random(7)
    res = []
    for _ in range(4):
        kappa = [Fraction(10) ** rng.randint(-2, 2) * Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in d]
        u = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(p)]
        cert = None
        for T in s2p.SCALES:
            v = s2p.stable_vertex(p, d, kappa, [T * x for x in u])
            if v is not None:
                cert = (T, v); break
        res.append((kappa, u, cert))
        if False:
            break  # E2's rng then draws LSODA-free; sequence identical only for cert prefixes
    nc = sum(r[2] is not None for r in res)
    if 0 < nc < len(res) or (nc and len(res) < 4):
        print(d2lib.describe(n, rx), "|", [(r[0], r[2]) for r in res])
