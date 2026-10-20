#!/usr/bin/env python3
"""Sign graph of the permissive (TG-literal) counterexample of result 12 at theta_3 (Agent N, network.md, Sec. 6).

Reads ../../../verification/permissive_counterexample_N5.json (the verifier's certified point; 40-digit centre of a
1e-30 box), computes the Underwood fluxes Lambda_s(theta_3) of all 18 sections, builds the sign graph (arc x -> y
when the flux flows from x to y) and prints every light -> heavy sign path.  By Proposition 5 of network.md, in a
physical steady state the node sign sigma = sign u(x; theta_3) can only go - -> + along an arc, and light products
are +, heavy products -, so such a path cannot exist; the script shows where the permissive point needs a + -> -
transition.  Floating point (mpmath, 40 digits) on data certified elsewhere; illustration only.
Environment: Python 3, mpmath.  Runtime: < 1 s.
"""
import json, os, sys
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 40
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(here, '..', '..', '..', 'verification', 'permissive_counterexample_N5.json')))
alpha = [mp.mpf(Fraction(a).numerator) / Fraction(a).denominator for a in d['alpha']]
th = [mp.mpf(Fraction(t).numerator) / Fraction(t).denominator for t in d['feed_roots']]
c = {k: mp.mpf(v) for k, v in d['center'].items()}
j = 3
theta = th[j - 1]
cols = {(1, 5): ((1, 4), (2, 5)), (1, 4): ((1, 3), (2, 4)), (2, 5): ((2, 4), (3, 5)), (1, 3): ((1, 2), (2, 3)),
        (2, 4): ((2, 3), (4, 4)), (3, 5): ((3, 3), (4, 5)), (1, 2): ((1, 1), (2, 2)), (2, 3): ((2, 2), (3, 3)),
        (4, 5): ((4, 4), (5, 5))}
name = lambda lab: 'ABCDE'[lab[0] - 1:lab[1]]
def flows(col, kind):
    out = {}
    for k, v in c.items():
        if k.startswith(f"{kind}({col[0]}, {col[1]})_"):
            out[int(k.split('_')[-1])] = v
    return out
arcs = []
print(f"theta_{j} = {mp.nstr(theta, 12)}")
for col, (dist, res) in cols.items():
    dT = flows(col, 'd'); bB = flows(col, 'b')
    VT = c[f"VT({col[0]}, {col[1]})"]; VB = c[f"VB({col[0]}, {col[1]})"]
    LT = VT - sum(alpha[p - 1] * w / (alpha[p - 1] - theta) for p, w in dT.items())
    LB = VB + sum(alpha[p - 1] * w / (alpha[p - 1] - theta) for p, w in bB.items())
    print(f"  column {name(col):6s} top -> {name(dist):5s} Lambda_T = {mp.nstr(LT, 8):>14s}   "
          f"bottom from {name(res):5s} Lambda_B = {mp.nstr(LB, 8):>14s}")
    # top section: bottom end at col, top end at dist
    arcs.append((name(col), name(dist), LT) if LT > 0 else (name(dist), name(col), -LT))
    arcs.append((name(res), name(col), LB) if LB > 0 else (name(col), name(res), -LB))
def region(n):
    i = 'ABCDE'.index(n[0]) + 1; m = 'ABCDE'.index(n[-1]) + 1
    return 'light' if m <= j else ('heavy' if i >= j + 1 else 'mixed')
adj = {}
for a, b, w in arcs:
    if w > mp.mpf('1e-25'):
        adj.setdefault(a, []).append(b)
paths = []
def dfs(v, path):
    if region(v) == 'heavy' and len(path) > 1:
        paths.append(list(path)); return
    for w in adj.get(v, []):
        if w not in path:
            path.append(w); dfs(w, path); path.pop()
for s in sorted({a for a, _, _ in arcs} | {b for _, b, _ in arcs}):
    if region(s) == 'light':
        dfs(s, [s])
print("light -> heavy sign paths:")
for p in paths:
    print("   " + " -> ".join(f"{n}({region(n)[0]})" for n in p))
print(("PASS" if paths else "FAIL") + "  the permissive point has a light -> heavy sign path (impossible physically)")
sys.exit(0 if paths else 1)
