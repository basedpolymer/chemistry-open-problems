#!/usr/bin/env python3
"""
enum4rc.py -- exploration for E2 (n = 4 species, exactly 4 reactant complexes, dim S = 3).
Stage 1 (exact): reactant sets with rank A = 3, v = ker A with exactly one zero entry (the ACR species) and the
other three entries distinct (necessary for 3 positive steady states in a class).
Stage 2 (exact): one reaction per reactant complex: rank Gamma = 3, ker Gamma positive, conservation law omega
with omega_j != 0 off the ACR species and a sign pattern allowing 3 sign changes of the exponential sum.
Prints counts and writes the survivors.
"""
import itertools, sys, json
import sympy as sp

n = 4
COMP = [c for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
PERMS = list(itertools.permutations(range(n)))


def pc(c, p):
    return tuple(c[p[i]] for i in range(n))


def canon_sets():
    out = []
    for sub in itertools.combinations(COMP, 4):
        key = tuple(sorted(sub))
        if key == min(tuple(sorted(pc(c, p) for c in sub)) for p in PERMS):
            out.append(key)
    return out


def intvec(v):
    den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
    w = [int(x * den) for x in v]
    g = 0
    for x in w:
        g = sp.igcd(g, x)
    return [x // g for x in w]


def alternation_possible(v, om, i):
    """3 sign changes possible for sum_{j!=i} om_j x_j e^{t v_j} + c0 (c0 free sign, exponent 0)?"""
    terms = [(v[j], 1 if om[j] > 0 else -1) for j in range(n) if j != i]
    for c0 in (1, -1):
        seq = sorted(terms + [(0, c0)])
        signs = [s for _, s in seq]
        if all(signs[k] != signs[k + 1] for k in range(3)):
            return True
    return False


sets = canon_sets()
st1 = []
for R in sets:
    A = sp.Matrix([[R[a][k] - R[3][k] for k in range(n)] for a in range(3)])
    if A.rank() < 3:
        continue
    v = intvec(list(A.nullspace()[0]))
    zeros = [k for k in range(n) if v[k] == 0]
    if len(zeros) != 1:
        continue
    i = zeros[0]
    others = [v[k] for k in range(n) if k != i]
    if len(set(others)) < 3:
        continue
    st1.append((R, v, i))
print("canonical reactant sets:", len(sets), " stage-1 survivors:", len(st1))
for R, v, i in st1:
    print("  ", R, "v =", v, "ACR species", i)

surv = []
cnt = dict(total=0, rank3=0, consistent=0, omega_ok=0)
for R, v, i in st1:
    prods = [[c for c in COMP if c != R[a]] for a in range(4)]
    for P in itertools.product(*prods):
        cnt['total'] += 1
        G = sp.Matrix([[P[a][k] - R[a][k] for a in range(4)] for k in range(n)])
        # quick: every species row must have + and - entries? not needed for conserved combos; use rank+kernel
        if G.rank() != 3:
            continue
        cnt['rank3'] += 1
        lam = G.nullspace()[0]
        if not (all(x > 0 for x in lam) or all(x < 0 for x in lam)):
            continue
        cnt['consistent'] += 1
        om = intvec(list(G.T.nullspace()[0]))
        if any(om[j] == 0 for j in range(n) if j != i):
            continue
        if alternation_possible(v, om, i):
            cnt['omega_ok'] += 1
            surv.append(dict(R=[list(r) for r in R], P=[list(p) for p in P], v=v, acr=i, omega=om))
print(cnt)
json.dump(surv, open("enum4rc_survivors.json", "w"))
for s in surv[:50]:
    print(s)
