#!/usr/bin/env python3
"""list_subspaces.py struct -- list all valid maximal networks (S spanned by reaction vectors, all 4 reactant
complexes used, positive flux) for a 4-reactant-complex structure, with exact conservation laws (rref)."""
import numpy as np, itertools, sys, json
import sympy as sp
from scipy.optimize import linprog
STRUCT = {
    'S3': [(0, 0, 0, 0, 1), (0, 0, 0, 0, 2), (0, 0, 1, 1, 0), (1, 1, 0, 0, 0)],
    'S7': [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)],
    'S9': [(0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1), (0, 0, 1, 1, 0, 0), (1, 1, 0, 0, 0, 0)],
}
name = sys.argv[1]
REAC = [np.array(r) for r in STRUCT[name]]
n = len(REAC[0])
COMP = [np.array(c) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
gens = [(a, p, p - REAC[a]) for a in range(4) for p in COMP if not np.array_equal(p, REAC[a])]
seen = {}
for (i, j, l) in itertools.combinations(range(len(gens)), 3):
    M = np.array([gens[i][2], gens[j][2], gens[l][2]], float)
    if np.linalg.matrix_rank(M) < 3:
        continue
    Ms = sp.Matrix(M.astype(int))
    Wb = Ms.nullspace()
    Wm = sp.Matrix.hstack(*Wb).T.rref()[0]
    key = tuple(Wm)
    if key in seen:
        continue
    W = np.array(Wm.tolist(), float).T
    inS = [g for g in range(len(gens)) if np.max(np.abs(gens[g][2] @ W)) < 1e-9]
    if set(gens[g][0] for g in inS) != {0, 1, 2, 3}:
        seen[key] = None
        continue
    G = np.array([gens[g][2] for g in inS], float).T
    if np.linalg.matrix_rank(G) != 3:
        seen[key] = None
        continue
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=np.zeros(n), bounds=[(1, None)] * G.shape[1], method="highs")
    if res.status != 0:
        seen[key] = None
        continue
    seen[key] = dict(W=[[str(x) for x in row] for row in Wm.tolist()],
                     rx=[(REAC[gens[g][0]].tolist(), gens[g][1].tolist()) for g in inS])
valid = [v for v in seen.values() if v]
print(name, "distinct 3-dim S spanned by reaction vectors:", len(seen), " valid (all reactant complexes, positive flux):",
      len(valid))
json.dump(valid, open("valid_%s.json" % name, "w"))
for v in valid[:400]:
    print(v['W'], len(v['rx']))
