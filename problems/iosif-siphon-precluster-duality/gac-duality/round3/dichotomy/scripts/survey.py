"""Survey of the generic boundary structure of N* on E2's undecided s = 2 networks (exploratory, no claim)."""
import os, sys, gzip, pickle, time
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d2lib
from itertools import combinations

recs = pickle.load(gzip.open(os.path.expanduser("~/gacD2-data/s2_all.pkl.gz")))
und = [q for q in recs if q["cat"] == "UNDECIDED"]
print(len(und))
C = Counter()
ex = {}
for q in und:
    N = d2lib.Net(q["n"], q["rx"])
    es = [N.edge(j) for j in range(N.r)]
    sip = [e for e in es if e["siphon"]]
    key_e = []
    for e in sip:
        cs = {1 if c > 0 else -1 for i, c in e["flow"]}
        fl = "zero" if not cs else ("uni" if len(cs) == 1 else "bi")
        key_e.append(f"{e['type']}/{fl}")
    veq = [(j, l) for j, l in combinations(range(N.r), 2) if N.independent(j, l) and N.vertex(j, l)["eq"]] if True else []
    key = (tuple(sorted(key_e)), len(veq))
    C[key] += 1
    ex.setdefault(key, (q["n"], q["rx"]))
from itertools import combinations
tot = Counter()
for (ke, nv), c in C.items():
    tot[("hasOUT" if any(k.startswith("OUT") for k in ke) else "") + ("hasMIX" if any(k.startswith("MIX") for k in ke) else "") + f" veq>0={nv>0}"] += c
for k, v in sorted(tot.items()):
    print(k, v)
for (ke, nv), c in C.most_common(40):
    n, rx = ex[(ke, nv)]
    print(c, ke, nv, "|", d2lib.describe(n, rx))
