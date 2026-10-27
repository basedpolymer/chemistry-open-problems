#!/usr/bin/env python3
"""count_S7.py -- for each valid S7 subspace (valid_S7.json), sample x* = (a,b,c,d,e) massively and record the
maximal number of positive steady states found in one class (class solver: resultant in w, see enum_S7.py)."""
import numpy as np, json, sys, time
sys.argv_ = sys.argv
seed = int(sys.argv[1]); NS = int(sys.argv[2])
sys.argv = ["x", "4", "1", "0", "0", "1"]
exec(open("enum_S7.py").read().split("prods = [[p")[0])
rng = np.random.default_rng(seed)
valid = json.load(open("valid_S7.json"))
from fractions import Fraction
hist = {}
worst = []
for v in valid:
    W = np.array([[float(Fraction(x)) for x in row] for row in v['W']]).T  # 5 x 2
    mx = 0
    for s in range(NS):
        xs = np.exp(rng.uniform(-6, 6, n))
        sols = class_solutions(W, xs)
        if sols is None:
            mx = max(mx, 99); break
        mx = max(mx, len(sols))
    hist[mx] = hist.get(mx, 0) + 1
    if mx >= 3:
        worst.append((v['W'], mx))
print("max #positive steady states per class (sampled):", hist)
for w in worst[:20]:
    print(w)
