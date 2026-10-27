"""Exact re-check of the smallest strong K8 witness pair found by k8_strong.py (reaction side).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact (crn.py via ecore).  Usage (from
round2/computation/): python scripts/k8_pair_check.py.  Claim (PASS/FAIL): both networks satisfy (H1)-(H3), are
weakly reversible, have rows and columns of Gamma pairwise linearly independent; some relabelling of the reactions of
the second matches the doubling classes and the minimal T-traps (minimal siphons of the dual) of the first, and for
some such relabelling the preclusters P (Definition 19, crn.py) differ (a labelled witness pair for K8, reaction
side; the two networks have 3 and 4 species and the same 5 reaction labels).
"""
import sys
import time
import platform
from itertools import permutations

import e2lib
from e2lib import CNet, ecore

N1 = (((0, 0, 2), (0, 2, 0)), ((0, 0, 2), (1, 1, 0)), ((0, 2, 0), (1, 0, 1)), ((1, 0, 1), (0, 0, 2)),
      ((1, 1, 0), (0, 2, 0)))
# N2 = 2X4 -> X3, X3 -> X2 + X4, X3 -> X1 + X4, X2 + X4 -> 2X4, X1 + X4 -> 2X4 (smallest pair of k8_strong.py)
N2 = (((0, 0, 0, 2), (0, 0, 1, 0)), ((0, 0, 1, 0), (0, 1, 0, 1)), ((0, 0, 1, 0), (1, 0, 0, 1)),
      ((0, 1, 0, 1), (0, 0, 0, 2)), ((1, 0, 0, 1), (0, 0, 0, 2)))


def props(rx):
    N = CNet(len(rx[0][0]), rx)
    D = N.dual()
    dbl = {}
    for j, (a, b) in enumerate(rx):
        dbl.setdefault(a, set()).add(j)
    return dict(H=N.is_conservative() and N.is_consistent() and N.rows_distinct(), WR=N.weakly_reversible(),
                strong=N.rows_pairwise_independent() and N.columns_pairwise_independent(),
                dbl=frozenset(frozenset(v) for v in dbl.values()),
                mt=frozenset(ecore.to_set(S) for S in ecore.minimal_masks(ecore.siphon_masks(D))),
                P=frozenset(frozenset(b) for b in N.preclusters("onestep")), desc=N.describe(), s=N.rank(),
                d=N.deficiency())


t0 = time.time()
print(f"# k8_pair_check.py  Python {platform.python_version()} on {platform.system()}")
p1, p2 = props(N1), props(N2)
for lab, p in (("N1", p1), ("N2", p2)):
    print(f"# {lab} = {p['desc']}: H={p['H']} WR={p['WR']} strong={p['strong']} s={p['s']} deficiency={p['d']}")
    print(f"#   doubling classes {sorted(map(sorted, p['dbl']))}, minimal T-traps {sorted(map(sorted, p['mt']))}, "
          f"P = {sorted(map(sorted, p['P']))}")
match = differ = 0
for pi in permutations(range(5)):
    rel = lambda fam: frozenset(frozenset(pi[j] for j in B) for B in fam)  # noqa: E731
    if rel(p2["dbl"]) == p1["dbl"] and rel(p2["mt"]) == p1["mt"]:
        match += 1
        differ += rel(p2["P"]) != p1["P"]
ok = all(p[k] for p in (p1, p2) for k in ("H", "WR", "strong")) and match > 0 and differ > 0
print(f"{'PASS' if ok else 'FAIL'}: relabellings matching doubling classes and minimal T-traps: {match}; with "
      f"different P: {differ}")
print(f"# runtime {time.time() - t0:.2f} s")
sys.exit(0 if ok else 1)
