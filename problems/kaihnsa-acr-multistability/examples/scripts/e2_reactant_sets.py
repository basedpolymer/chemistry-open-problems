#!/usr/bin/env python3
"""
e2_reactant_sets.py -- exact enumeration behind Section 4 of report.md (task E2, exactly 4 reactant complexes).

Runtime: about 10-15 minutes (n = 8 dominates; WSL, Python 3.12, sympy).   Usage: python e2_reactant_sets.py [nmax]

For n = 4..8 species, enumerates all sets of 4 at-most-bimolecular reactant complexes that involve every species
(a species absent from all reactant complexes is never consumed, so no positive steady state exists), with
rank A = 3 (A = rows y_a - y_4; rank A < 3 makes every positive steady state degenerate in a class of dimension 3),
and with an ACR species (a zero row of a basis of ker A; with 4 reactant complexes and rank N = 3 the positive
steady states form the coset x* o exp(ker A), so ACR in X_i <=> e_i is orthogonal to ker A).
Claims checked (PASS/FAIL):
  n = 4: for every such set, the entries of v (ker A = span v) off the ACR species are NOT three distinct nonzero
         numbers  (=> at most 3 distinct exponents in the exponential sum of a class => <= 2 steady states);
  n = 5: exactly 8 such sets up to permutation of species (S1..S8 of report.md);
  n = 6: exactly 1 (S9);  n = 7, 8: none.
"""
import itertools, sys, time
import numpy as np
import sympy as sp

T0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
FAIL = []


def check(c, msg):
    print(("PASS " if c else "FAIL ") + msg, flush=True)
    if not c:
        FAIL.append(msg)


def structures(n):
    COMP = [c for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
    seen = {}
    for sub in itertools.combinations(COMP, 4):
        cover = np.array(sub).sum(axis=0)
        if np.any(cover == 0):
            continue
        A = np.array([np.array(sub[a]) - np.array(sub[3]) for a in range(3)])
        if np.linalg.matrix_rank(A) < 3:
            continue
        As = sp.Matrix(A.tolist())
        if As.rank() < 3:
            continue
        V = sp.Matrix.hstack(*As.nullspace())
        zr = [i for i in range(n) if all(V[i, j] == 0 for j in range(V.shape[1]))]
        if not zr:
            continue
        best = min(tuple(sorted(tuple(c[p[i]] for i in range(n)) for c in sub)) for p in itertools.permutations(range(n)))
        Ab = sp.Matrix([[best[a][i] - best[3][i] for i in range(n)] for a in range(3)])
        Vb = sp.Matrix.hstack(*Ab.nullspace())
        zb = [i for i in range(n) if all(Vb[i, j] == 0 for j in range(Vb.shape[1]))]
        seen[best] = (Vb, zb)
    return seen


expected = {5: 8, 6: 1, 7: 0, 8: 0}
for n in range(4, NMAX + 1):
    st = structures(n)
    print("n = %d: %d reactant sets (up to species permutation) with rank A = 3 and an ACR species  [%.0f s]"
          % (n, len(st), time.time() - T0), flush=True)
    if n == 4:
        bad = []
        for key, (V, zr) in st.items():
            v = list(V[:, 0])
            off = [v[i] for i in range(n) if i not in zr]
            if len(zr) == 1 and len(set(off)) == 3 and all(x != 0 for x in off):
                bad.append(key)
        check(not bad, "n = 4: no reactant set has an ACR species with three distinct nonzero exponents off it "
                       "(%d sets checked)" % len(st))
    else:
        check(len(st) == expected[n], "n = %d: number of structures = %d" % (n, expected[n]))
        for key, (V, zr) in sorted(st.items()):
            print("    ", key, "ACR species (0-based):", zr, " ker A basis^T:", V.T.tolist())
print("runtime %.0f s" % (time.time() - T0))
sys.exit(1 if FAIL else 0)
