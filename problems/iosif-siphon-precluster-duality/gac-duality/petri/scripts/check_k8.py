"""Candidate K8 (agent P): the preclusters of N* are NOT determined by the twin classes, the species parts of the
clusters and the minimal siphons of N outside the well-formed free-choice class.
Environment: Python 3.12+, standard library only (crn.py, petri.py), exact rationals.  Runtime < 1 s.
The pair was found by `k8_search.py 4 2 3 6 360 1` (random search, seed 1).  PASS/FAIL per claim, exit 1 on failure.
"""
import sys
import platform
from petri import Network, twin_classes, clusters, is_EC

print(f"# check_k8.py  Python {platform.python_version()} on {platform.system()}")
FAILS = []


def claim(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILS.append(name)


N1 = Network(4, [((0, 0, 1, 0), (0, 0, 0, 1)), ((0, 1, 1, 0), (1, 1, 0, 0)), ((1, 0, 0, 1), (0, 0, 2, 0))])
N2 = Network(4, [((0, 0, 2, 0), (0, 0, 0, 2)), ((0, 1, 1, 0), (0, 2, 0, 0)), ((0, 1, 1, 0), (2, 0, 0, 0)),
                 ((1, 0, 0, 1), (0, 0, 2, 0))])
for name, N in (("N1", N1), ("N2", N2)):
    print(f"#   {name}: {N.describe()}; extreme conservation laws {N.conservation_rays()}; P(N*) = "
          f"{[sorted(b) for b in N.dual().preclusters()]}; EC: {is_EC(N)}")
    claim(f"{name} satisfies (H1), (H2), (H3)", N.is_conservative() and N.is_consistent() and N.rows_distinct())
key = lambda part: sorted(sorted(b) for b in part)  # noqa: E731
claim("same twin classes", key(twin_classes(N1)) == key(twin_classes(N2)))
claim("same species parts of the clusters",
      key([c[0] for c in clusters(N1)]) == key([c[0] for c in clusters(N2)]))
claim("same minimal siphons", key(N1.minimal_siphons()) == key(N2.minimal_siphons()))
claim("different preclusters of the dual", key(N1.dual().preclusters()) != key(N2.dual().preclusters()))
claim("the three precluster variants agree on both duals",
      all(len({str(key(N.dual().preclusters(v))) for v in ("onestep", "thesis", "closure")}) == 1 for N in (N1, N2)))
# reaction side (found by check_general.py on the family n=3, molecularity <= 2, r <= 4): P(N) is not determined by
# the doubling classes and the minimal T-traps (= minimal siphons of N*); no species has a zero row of Gamma.
M1 = Network(3, [((0, 0, 1), (1, 1, 0)), ((0, 1, 1), (1, 0, 1)), ((1, 1, 0), (0, 0, 1)), ((1, 1, 0), (0, 2, 0))])
M2 = Network(3, [((0, 0, 2), (0, 1, 1)), ((1, 0, 1), (1, 1, 0)), ((1, 1, 0), (0, 0, 2)), ((1, 1, 0), (2, 0, 0))])
for name, N in (("M1", M1), ("M2", M2)):
    print(f"#   {name}: {N.describe()}; extreme flux rays {N.flux_rays()}; P(N) = "
          f"{[sorted(b) for b in N.preclusters()]}; EC: {is_EC(N)}")
    claim(f"{name} satisfies (H1), (H2), (H3); no zero row of Gamma (columns distinct: {N.columns_distinct()})",
          N.is_conservative() and N.is_consistent() and N.rows_distinct() and all(any(row) for row in N.Gamma))
claim("M1, M2: same doubling classes (= twin classes of the duals)",
      key(twin_classes(M1.dual())) == key(twin_classes(M2.dual())))
claim("M1, M2: same minimal T-traps (= minimal siphons of the duals)",
      key(M1.dual().minimal_siphons()) == key(M2.dual().minimal_siphons()))
claim("M1, M2: different preclusters", key(M1.preclusters()) != key(M2.preclusters()))
sys.exit(1 if FAILS else 0)
