"""Verifier V2: exact re-check of Q3's K8 strong pair M3/M4 (round2/ec/notes.md §6) and of the small pair M1/M2.

Environment: WSL Fedora 44, Python 3.12.  Independent precluster implementation (Definition 19, one-step span rule)
compared with the repository library scripts/crn.py (all three variants).  PASS/FAIL, exit 1 on FAIL.
"""
import sys, os, time
from itertools import permutations, combinations
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
from v2lib import *  # noqa

import crn

T0 = time.time()
FAILS = []


def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + info) if info else ""))
    if not cond:
        FAILS.append(name)


def weakly_reversible(N):
    cs = []
    for a, b in N:
        for c in (a, b):
            if c not in cs:
                cs.append(c)
    idx = {c: i for i, c in enumerate(cs)}
    m = len(cs)
    R = [[i == j for j in range(m)] for i in range(m)]
    for a, b in N:
        R[idx[a]][idx[b]] = True
    for k in range(m):
        for i in range(m):
            if R[i][k]:
                for j in range(m):
                    if R[k][j]:
                        R[i][j] = True
    return all(R[idx[b]][idx[a]] for a, b in N)


def pure(N):
    return all(not (x and y) for a, b in N for x, y in zip(a, b))


def pairwise_indep(rows):
    return all(rank([list(rows[i]), list(rows[j])]) == 2 for i, j in combinations(range(len(rows)), 2))


def doubling_classes(N):
    cl = {}
    for t, (a, _) in enumerate(N):
        cl.setdefault(a, []).append(t)
    return sorted(tuple(v) for v in cl.values())


def min_T_traps(N):
    """Minimal T-traps of N = minimal siphons of N* (as reaction index sets)."""
    return sorted(tuple(S) for S in minimal_siphons(dual(N)))


def preclusters_v2(N):
    """Definition 19 (one step): graph on reactions with edges {i,j} for equal sources; for each such edge and each
    s with n_s in span{n_i, n_j} add {i,s},{j,s}; n_t = row t of the matrix of extreme flux rays."""
    rays = flux_rays(N)
    r = len(N)
    nv = [[v[t] for v in rays] for t in range(r)]
    G = [(i, j) for i, j in combinations(range(r), 2) if N[i][0] == N[j][0]]
    E = set(G)
    for i, j in G:
        base = rank([nv[i], nv[j]])
        for s in range(r):
            if s not in (i, j) and rank([nv[i], nv[j], nv[s]]) == base:
                E.add(tuple(sorted((i, s))))
                E.add(tuple(sorted((j, s))))
    parent = list(range(r))

    def f(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for i, j in E:
        parent[f(i)] = f(j)
    blocks = {}
    for t in range(r):
        blocks.setdefault(f(t), []).append(t)
    return sorted(tuple(b) for b in blocks.values())


def crn_net(N):
    return crn.Network(nsp(N), N)


def relabel(blocks, pi):
    return sorted(tuple(sorted(pi[t] for t in b)) for b in blocks)


def pair_check(tag, A, B, expect_matches=None):
    for nm, N in ((tag + "a", A), (tag + "b", B)):
        D = dual(N)
        check(f"{nm}: (H1)-(H3) on N and N*, WR, connected", H123(N) and H123(D) and weakly_reversible(N)
              and connected(N))
        cr = crn_net(N)
        pv = preclusters_v2(N)
        pc = [sorted(tuple(sorted(b)) for b in cr.preclusters(variant=v)) for v in ("onestep", "thesis", "closure")]
        check(f"{nm}: independent preclusters agree with crn.py (3 variants)", all(sorted(p) == pv for p in pc),
              f"{pv}")
        print(f"     {nm}: pure={pure(N)} rows-indep={pairwise_indep(Gamma(N))} "
              f"cols-indep={pairwise_indep(transpose(Gamma(N)))} doubling={doubling_classes(N)} "
              f"minTtraps={min_T_traps(N)} P={pv} s={rank(Gamma(N))}")
    dA, dB = doubling_classes(A), doubling_classes(B)
    tA, tB = min_T_traps(A), min_T_traps(B)
    pA, pB = preclusters_v2(A), preclusters_v2(B)
    r = len(A)
    matches = 0
    differ = 0
    for pi in permutations(range(r)):
        if relabel(dB, pi) == dA and relabel(tB, pi) == tA:
            matches += 1
            if relabel(pB, pi) != pA:
                differ += 1
    check(f"{tag}: relabellings matching doubling classes and minimal T-traps: {matches}; preclusters differ under all",
          matches > 0 and differ == matches and (expect_matches is None or matches == expect_matches))


sp_ = ["X1", "X2", "X3", "X4"]
M3 = parse("X4->2X3, X4->X2+X3, 2X3->2X2, X2+X3->2X1, 2X2->X4, 2X1->2X3", sp_)
M4 = parse("X4->X2+X3, 2X3->X4, 2X3->2X2, X2+X3->2X1, 2X2->X4, 2X1->2X3", sp_)
print("== M3 / M4")
pair_check("M3/M4", M3, M4, expect_matches=8)
for nm, N in (("M3", M3), ("M4", M4)):
    check(f"{nm}: pure, rows and columns of Gamma pairwise linearly independent",
          pure(N) and pairwise_indep(Gamma(N)) and pairwise_indep(transpose(Gamma(N))))
    D = dual(N)
    check(f"{nm}*: rows and columns pairwise independent, (H1)-(H3)", H123(D) and pairwise_indep(Gamma(D))
          and pairwise_indep(transpose(Gamma(D))))
print("== M1 / M2 (Q3's smallest; E2's smallest strong pair)")
M1 = parse("2X3->2X2, 2X3->X1+X2, 2X2->X1+X3, X1+X3->2X3, X1+X2->2X2", ["X1", "X2", "X3"])
M2 = parse("2X4->X3, X3->X2+X4, X3->X1+X4, X2+X4->2X4, X1+X4->2X4", sp_)
pair_check("M1/M2", M1, M2)
for nm, N in (("M1", M1), ("M2", M2)):
    print(f"     {nm}: pairwise-indep rows/cols N: {pairwise_indep(Gamma(N))}/{pairwise_indep(transpose(Gamma(N)))}"
          f"; N*: {pairwise_indep(Gamma(dual(N)))}/{pairwise_indep(transpose(Gamma(dual(N))))}")
print(f"runtime {time.time() - T0:.1f} s; {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
