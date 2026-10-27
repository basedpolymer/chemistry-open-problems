"""Exact check of the K8 strong witness, reaction side (agent Q3, Claude Opus 5.5, effort xhigh; notes.md §6).

Environment: Python 3.12+, standard library only, exact rationals (crn.py).  Run from round2/ec/:
python scripts/check_k8s.py.  Runtime < 5 s.  PASS/FAIL per claim, exit 1 on FAIL.

M1 (n = 3, r = 5) and M2 (n = 4, r = 5), found by scripts/k8_strong.py in agent E's records:
  M1: 2X3 -> 2X2, 2X3 -> X1+X2, 2X2 -> X1+X3, X1+X3 -> 2X3, X1+X2 -> 2X2
  M2: 2X4 -> X3, X3 -> X2+X4, X3 -> X1+X4, X2+X4 -> 2X4, X1+X4 -> 2X4
Second (pure) pair on the same species set, found by PURE=1 scripts/k8_strong.py:
  M3: X4 -> 2X3, X4 -> X2+X3, 2X3 -> 2X2, X2+X3 -> 2X1, 2X2 -> X4, 2X1 -> 2X3
  M4: X4 -> X2+X3, 2X3 -> X4, 2X3 -> 2X2, X2+X3 -> 2X1, 2X2 -> X4, 2X1 -> 2X3
Claims (each pair): both weakly reversible, (H1)-(H3) on M and on M*, rows and columns of Gamma pairwise linearly independent;
there is a relabelling of the reactions carrying the doubling classes and the minimal T-traps (minimal siphons of
M*) of M1 onto those of M2; for EVERY such relabelling the preclusters P(M1), P(M2) (all three variants of crn.py)
do not correspond.  So P(N) is not a function of the doubling classes and the minimal T-traps, even for networks
with pairwise linearly independent rows and columns (dually: P(N*) is not a function of the twin classes and the
minimal siphons of N, with N = M*).
"""
import sys
import time
import platform
from itertools import permutations

import eclib
from eclib import Network

t0 = time.time()
print(f"# check_k8s.py; Python {platform.python_version()} on {platform.system()}")
bad = 0


def claim(name, ok, extra=""):
    global bad
    print(("PASS" if ok else "FAIL") + f"  {name}" + (f"  {extra}" if extra else ""))
    bad += not ok


M1 = Network(3, [((0, 0, 2), (0, 2, 0)), ((0, 0, 2), (1, 1, 0)), ((0, 2, 0), (1, 0, 1)),
                 ((1, 0, 1), (0, 0, 2)), ((1, 1, 0), (0, 2, 0))])
M2 = Network(4, [((0, 0, 0, 2), (0, 0, 1, 0)), ((0, 0, 1, 0), (0, 1, 0, 1)), ((0, 0, 1, 0), (1, 0, 0, 1)),
                 ((0, 1, 0, 1), (0, 0, 0, 2)), ((1, 0, 0, 1), (0, 0, 0, 2))])


def doubling(N):
    cl = {}
    for t, (a, b) in enumerate(N.reactions):
        cl.setdefault(a, []).append(t)
    return frozenset(frozenset(v) for v in cl.values())


def permset(fam, pi):
    return frozenset(frozenset(pi[x] for x in B) for B in fam)


def check_pair(lab1, A, lab2, B):
    info = {}
    for name, M in ((lab1, A), (lab2, B)):
        D = M.dual()
        claim(f"{name} = {M.describe()}: weakly reversible", M.weakly_reversible())
        claim(f"{name}: (H1)-(H3) on M and on M*", M.is_conservative() and M.is_consistent() and M.rows_distinct()
              and D.is_conservative() and D.is_consistent() and D.rows_distinct())
        claim(f"{name}: rows and columns of Gamma pairwise linearly independent",
              M.rows_pairwise_independent() and M.columns_pairwise_independent())
        P = {v: frozenset(frozenset(b) for b in M.preclusters(v)) for v in ("onestep", "thesis", "closure")}
        claim(f"{name}: the three precluster variants agree", len(set(P.values())) == 1)
        MT = frozenset(frozenset(x) for x in D.minimal_siphons())
        info[name] = (doubling(M), MT, P["onestep"])
        print(f"#   {name}: pure={eclib.is_pure(M)}, doubling classes {sorted(map(sorted, doubling(M)))}, minimal "
              f"T-traps {sorted(map(sorted, MT))}, P = {sorted(map(sorted, P['onestep']))}, s = {M.rank()}, "
              f"deficiency = {M.deficiency()}")
    d1, t1, p1 = info[lab1]
    d2, t2, p2 = info[lab2]
    r = A.r
    matching = [pi for pi in permutations(range(r)) if permset(d1, pi) == d2 and permset(t1, pi) == t2]
    claim(f"[{lab1},{lab2}] some relabelling of the reactions carries (doubling classes, minimal T-traps) onto "
          f"each other", len(matching) > 0, f"{len(matching)} relabellings")
    claim(f"[{lab1},{lab2}] for every such relabelling the preclusters differ",
          all(permset(p1, pi) != p2 for pi in matching))


check_pair("M1", M1, "M2", M2)
# pure witness on one species set (found by PURE=1 k8_strong.py)
M3 = Network(4, [((0, 0, 0, 1), (0, 0, 2, 0)), ((0, 0, 0, 1), (0, 1, 1, 0)), ((0, 0, 2, 0), (0, 2, 0, 0)),
                 ((0, 1, 1, 0), (2, 0, 0, 0)), ((0, 2, 0, 0), (0, 0, 0, 1)), ((2, 0, 0, 0), (0, 0, 2, 0))])
M4 = Network(4, [((0, 0, 0, 1), (0, 1, 1, 0)), ((0, 0, 2, 0), (0, 0, 0, 1)), ((0, 0, 2, 0), (0, 2, 0, 0)),
                 ((0, 1, 1, 0), (2, 0, 0, 0)), ((0, 2, 0, 0), (0, 0, 0, 1)), ((2, 0, 0, 0), (0, 0, 2, 0))])
claim("M3, M4 pure", eclib.is_pure(M3) and eclib.is_pure(M4))
check_pair("M3", M3, "M4", M4)

print(f"# runtime {time.time() - t0:.2f} s")
sys.exit(1 if bad else 0)
