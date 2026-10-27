"""Exact checks of the Petri-net dictionary and of the well-formed free-choice obstruction (agent P).

Environment: Python 3.12, standard library only (crn.py + petri.py, exact rationals).  Runtime: a few seconds.
Prints PASS/FAIL per claim and exits with status 1 on any failure.

Claims checked (numbering of petri/notes.md):
  D1  N* is the weighted reverse-dual of N (Pre and Post transposed), N** = N.
  D2  siphons of N = T-traps of N*, traps of N = T-siphons of N*, and the same with N, N* exchanged.
  D3  clusters of N = clusters of N* (same partition of species + reactions).
  D4  equal-conflict sets of N* = twin classes of N; N is EC & balanced-joins iff N* is;
      for EC & balanced-joins networks the equal-conflict sets are the reaction parts of the clusters with
      a nonzero source and the twin classes are the species parts of the clusters.
  O1  R, T, W, W*, F are connected ordinary free-choice networks satisfying (H1), (H2), (H3) and the rank equation,
      hence well-formed (Rank Theorem); their duals are too (Duality Theorem); token-game liveness confirmed.
  O2  every formulation of Table 1 of the manuscript fails on one of R, T, W, W*, F.
  O3  U, D, N3 are not free choice; D is EC but not well-formed (rank equation fails), N3 is not EC.
  P1  pairwise span rule: in the 2-branch and 3-branch well-formed free-choice networks below, the join reaction
      lies in the precluster of the choice for 2 branches and not for 3 branches, while the class-span closure
      (Conradi-Flockerzi ratio closure) contains it in both.
"""
import sys
import time
import platform

from petri import (Network, is_reverse_dual, clusters, ec_sets, twin_classes, is_EC, is_free_choice,
                   is_balanced_joins, is_ordinary, is_connected, rank_equation_clusters, rank_equation_ec,
                   well_formed_by_rank_theorem, is_t_trap, is_t_siphon, is_trap, t_traps, minimal_t_traps,
                   is_live_and_bounded, ts_marking, refines, span_closure_preclusters)
from itertools import combinations

T0 = time.time()
print(f"# check_dictionary_examples.py  Python {platform.python_version()} on {platform.system()}")
FAILS = []


def claim(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILS.append(name)


def cx(n, *terms):
    v = [0] * n
    for c, i in terms:
        v[i] += c
    return tuple(v)


def net(n, rx):
    return Network(n, [(cx(n, *a), cx(n, *b)) for a, b in rx])


one = lambda i: (1, i)  # noqa: E731
R = net(2, [([one(0)], [one(1)]), ([one(1)], [one(0)])])
T = net(3, [([one(0)], [one(1)]), ([one(1)], [one(2)]), ([one(2)], [one(0)])])
W = net(3, [([one(0)], [one(1)]), ([one(0)], [one(2)]), ([one(1)], [one(0)]), ([one(2)], [one(0)])])
Wd = W.dual()
U = net(4, [([one(3)], [one(1), one(2)]), ([one(2), one(3)], [one(0)]), ([one(0), one(1)], [(2, 3)])])
D = net(3, [([one(0)], [one(1)]), ([one(0)], [one(2)]), ([one(1), one(2)], [(2, 0)])])
# F: species S0,S1,E,F,ES0,FS1 = 0..5
F = net(6, [([one(0), one(2)], [one(4)]), ([one(4)], [one(0), one(2)]), ([one(4)], [one(1), one(2)]),
            ([one(1), one(3)], [one(5)]), ([one(5)], [one(1), one(3)]), ([one(5)], [one(0), one(3)])])
N3 = net(2, [([(2, 0)], [one(0), one(1)]), ([one(0), one(1)], [(2, 0)]), ([one(0), one(1)], [(2, 1)])])
N412 = net(2, [([(2, 0)], [one(0), one(1)]), ([one(0), one(1)], [(2, 0)])])
# k-branch choice: X -> Y_a (a=1..k), Y_a -> Z, Z -> X ; species X=0, Y_a=a, Z=k+1
def branch(k):
    n = k + 2
    rx = [([one(0)], [one(a)]) for a in range(1, k + 1)]
    rx += [([one(a)], [one(k + 1)]) for a in range(1, k + 1)]
    rx += [([one(k + 1)], [one(0)])]
    return net(n, rx)


B2, B3 = branch(2), branch(3)
NETS = {"R": R, "T": T, "W": W, "W*": Wd, "U": U, "D": D, "F": F, "N3": N3, "N4.12": N412, "B2": B2, "B3": B3}

# ---------------------------------------------------------------- D1-D4
for name, N in NETS.items():
    Ns = N.dual()
    claim(f"D1 {name}: N* is the weighted reverse-dual of N, N** = N",
          is_reverse_dual(N, Ns) and Ns.dual().reactions == N.reactions)
    sub_sp = [set(S) for k in range(1, N.n + 1) for S in combinations(range(N.n), k)]
    ok2 = all(N.is_siphon(S) == is_t_trap(Ns, S) and is_trap(N, S) == is_t_siphon(Ns, S) for S in sub_sp)
    sub_rx = [set(U_) for k in range(1, N.r + 1) for U_ in combinations(range(N.r), k)]
    ok2 = ok2 and all(Ns.is_siphon(U_) == is_t_trap(N, U_) and is_trap(Ns, U_) == is_t_siphon(N, U_)
                      for U_ in sub_rx)
    claim(f"D2 {name}: siphons/traps of N = T-traps/T-siphons of N*, and conversely", ok2)
    cN = clusters(N)
    cD = {(c[1], c[0]) for c in clusters(Ns)}  # species of N* = reactions of N
    claim(f"D3 {name}: clusters of N = clusters of N*", set(cN) == cD)
    claim(f"D4 {name}: equal-conflict sets of N* = twin classes of N", ec_sets(Ns) == twin_classes(N))
    claim(f"D4 {name}: N EC+balanced joins iff N* EC+balanced joins",
          (is_EC(N) and is_balanced_joins(N)) == (is_EC(Ns) and is_balanced_joins(Ns)))
    if is_EC(N) and is_balanced_joins(N):
        ecs = sorted(sorted(c[1]) for c in cN if c[1])
        tws = sorted(sorted(c[0]) for c in cN if c[0] and any(N.Ye[i][j] for i in c[0] for j in range(N.r)))
        claim(f"D4 {name}: (EC, balanced) EC sets / twin classes = reaction / species parts of clusters",
              ecs == sorted(sorted(b) for b in ec_sets(N)) and
              tws == sorted(sorted(b) for b in twin_classes(N) if any(N.Ye[i][j] for i in b for j in range(N.r))))

# ---------------------------------------------------------------- O1, O3
def H(N):
    return N.is_conservative() and N.is_consistent() and N.rows_distinct()


for name in ["R", "T", "W", "W*", "F", "B2", "B3"]:
    N = NETS[name]
    Ns = N.dual()
    wf = (is_free_choice(N) and is_connected(N) and H(N) and rank_equation_clusters(N))
    wfd = (is_free_choice(Ns) and is_connected(Ns) and Ns.is_conservative() and Ns.is_consistent()
           and rank_equation_clusters(Ns))
    ok_b, live = is_live_and_bounded(N, ts_marking(N))
    ok_bd, lived = is_live_and_bounded(Ns, ts_marking(Ns))
    claim(f"O1 {name}: ordinary FC, connected, (H1)-(H3), rank = #clusters-1 = {len(clusters(N)) - 1}", wf)
    claim(f"O1 {name}: dual is well-formed FC too (Duality Theorem) and both are live+bounded at M0", wfd and ok_b
          and live and ok_bd and lived)

claim("O3 U is not EC (hence not free choice)", not is_EC(U))
claim("O3 D is EC, balanced joins, not ordinary, (H1)-(H3), rank 2 != #EC sets - 1 = 1: not well-formed",
      is_EC(D) and is_balanced_joins(D) and not is_ordinary(D) and H(D) and not rank_equation_ec(D)
      and well_formed_by_rank_theorem(D) is False)
okD, liveD = is_live_and_bounded(D, ts_marking(D))
claim("O3 D: token game from M0 = (1,1,1) is bounded and not live", okD and liveD is False)
claim("O3 N3 is not EC", not is_EC(N3))

# ---------------------------------------------------------------- O2
def formulations(N):
    Ns = N.dual()
    out = {}
    for side, (P, M, S, C) in {
        "A": (Ns.preclusters(), N.mips_ray(), N.minimal_siphons(), None),
        "B": (N.preclusters(), Ns.mips_ray(), Ns.minimal_siphons(), N.internal_cycles()),
    }.items():
        out[side + "1 P=M"] = P == M
        out[side + "2 P<=M"] = refines(P, M)
        out[side + "3 M<=P"] = refines(M, P)
        out[side + "4 |P|=|M|"] = len(P) == len(M)
        out["Q1" + side + " P=S"] = sorted(map(sorted, P)) == sorted(map(sorted, S))
        out["Q2" + side + " S<=P"] = refines(S, P)
        out["Q3" + side + " P<=S"] = refines(P, S)
        net_ = N if side == "A" else Ns
        out["Q4" + side + " every P block is a siphon"] = all(net_.is_siphon(b) for b in P)
        out["Q5" + side + " every P block is a union of min siphons"] = all(
            set().union(*([s for s in S if s <= b] or [set()])) == set(b) for b in P)
        out["Q6" + side + " |P|=|S|"] = len(P) == len(S)
        if C is not None:
            out["R23 cycles<=P"] = refines(C, P)
    out["C |P(N)|=|M(N)|"] = len(N.preclusters()) == len(N.mips_ray())
    return out


fails = {}
for name in ["R", "T", "W", "W*", "F"]:
    for st, val in formulations(NETS[name]).items():
        if not val:
            fails.setdefault(st, []).append(name)
all_st = sorted(formulations(T).keys())
for st in all_st:
    claim(f"O2 '{st}' fails on a well-formed FC network: {fails.get(st, [])}", st in fails)

# ---------------------------------------------------------------- P1
for name, N, expect_join in [("B2", B2, True), ("B3", B3, False)]:
    P = N.preclusters()
    join = N.r - 1
    choice = 0
    same = any(join in b and choice in b for b in P)
    claim(f"P1 {name}: Iosif precluster of the choice contains the join reaction: {same} (expected {expect_join});"
          f" P = {[sorted(b) for b in P]}", same == expect_join)
    Q = span_closure_preclusters(N)
    claim(f"P1 {name}: class-span (ratio) closure is one block: {[sorted(b) for b in Q]}", len(Q) == 1)
    claim(f"P1 {name}: all three Iosif variants agree",
          N.preclusters("onestep") == N.preclusters("thesis") == N.preclusters("closure"))

print(f"# runtime {time.time() - T0:.2f} s; {len(FAILS)} failure(s)")
sys.exit(1 if FAILS else 0)
