#!/usr/bin/env python3
"""stage1_n.py n -- 4 bimolecular reactant complexes covering all n species, rank A = 3, ACR species present."""
import itertools, sys
import numpy as np
import sympy as sp

n = int(sys.argv[1])
COMP = [c for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
seen = set()
out = []
for sub in itertools.combinations(COMP, 4):
    cover = np.array(sub).sum(axis=0)
    if np.any(cover == 0):
        continue
    A = np.array([np.array(sub[a]) - np.array(sub[3]) for a in range(3)], float)
    if np.linalg.matrix_rank(A) < 3:
        continue
    As = sp.Matrix(A.astype(int))
    V = As.nullspace()
    Vm = sp.Matrix.hstack(*V)
    zero_rows = [i for i in range(n) if all(Vm[i, k] == 0 for k in range(Vm.shape[1]))]
    if not zero_rows:
        continue
    # canonical form under species permutations (only for printing distinct structures)
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(sorted(tuple(c[p[i]] for i in range(n)) for c in sub))
        if best is None or key < best:
            best = key
    if best in seen:
        continue
    seen.add(best)
    out.append((best, zero_rows))
print("n =", n, "distinct reactant sets (up to S_n) with rank A = 3 and an ACR species:", len(out))
for b, z in out:
    A = sp.Matrix([[b[a][k] - b[3][k] for k in range(n)] for a in range(3)])
    V = sp.Matrix.hstack(*A.nullspace())
    zr = [i for i in range(n) if all(V[i, k] == 0 for k in range(V.shape[1]))]
    print(b, "ACR species", zr, "V^T =", V.T.tolist())
