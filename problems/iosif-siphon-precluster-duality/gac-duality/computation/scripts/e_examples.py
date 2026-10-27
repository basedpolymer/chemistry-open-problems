"""Explicit examples of agent E (counterexamples and witnesses quoted in notes.md), checked exactly.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Runtime: a few seconds.

Claims (PASS/FAIL):
  E1  Y5: A->B, B->C, B->D, B->E, C->A, D->A, E->A (monomolecular, one linkage class, WR, deficiency 0,
      (H1)-(H3)): the join of P(N) and M(N*) is not the linkage-class partition (A->B is a singleton of both); a
      monomolecular illustration of the failure of E-c1' (tests.py finds smaller ones).
  E2  D = 3U+2V -> 2U+3V, V -> U (dual of G's M): for EVERY rate constants k1, k2 > 0 and every class u+v = T > 0,
      v' = v (k1 u^3 v - k2) <= -(k2/2) v whenever 0 < v <= v* = min(T, k2/(2 k1 T^3)), so v -> 0: D is
      non-persistent for every choice of rate constants (checked: the bound k1 T^3 v* <= k2/2 on a grid, and the
      'drain species' certificate of endo.py); M is reversible with one linkage class (Anderson 2011: globally
      attracting for all rate constants, here even explicitly).
  E3  Smallest drain-species witnesses found by tests.py (X-drainD) re-checked: N WR, one linkage class,
      deficiency 0, (H1)-(H3) (so complex balanced and globally attracting for all rate constants: Anderson 2011),
      N* with (H1)-(H3) has a drain species (non-persistent for all rate constants, notes.md Prop. E.4).
"""
import sys
import time
import platform
from fractions import Fraction as F

import ecore
import endo
from ecore import CNet, refines

t0 = time.time()
print(f"# e_examples.py  Python {platform.python_version()} on {platform.system()} {platform.release()}")
OK = True


def claim(tag, cond, doc):
    global OK
    OK &= bool(cond)
    print(f"{'PASS' if cond else 'FAIL'} {tag}: {doc}")


def H(N):
    return N.is_conservative() and N.is_consistent() and N.rows_distinct()


def e(i, n=5):
    return tuple(int(j == i) for j in range(n))


A, B, C, D_, E = range(5)
Y5 = CNet(5, [(e(A), e(B)), (e(B), e(C)), (e(B), e(D_)), (e(B), e(E)), (e(C), e(A)), (e(D_), e(A)), (e(E), e(A))],
          names=list("ABCDE"))
P = Y5.preclusters("onestep")
M = Y5.dual().mips_ray()
blocks = set(frozenset(b) for b in P) | set(frozenset(b) for b in M)
# join
parent = list(range(Y5.r))


def find(x):
    while parent[x] != x:
        x = parent[x]
    return x


for b in blocks:
    b = sorted(b)
    for x in b[1:]:
        parent[find(x)] = find(b[0])
join = {}
for x in range(Y5.r):
    join.setdefault(find(x), set()).add(x)
claim("E1.1", Y5.weakly_reversible() and Y5.deficiency() == 0 and H(Y5), "Y5 is WR, deficiency 0, (H1)-(H3)")
claim("E1.2", frozenset({0}) in set(P) and frozenset({0}) in set(M) and len(join) == 2,
      f"P(Y5) = {[sorted(b) for b in P]}, M(Y5*) = {[sorted(b) for b in M]}: the join has 2 blocks, Y5 one "
      f"linkage class")

# E2
okb = True
for k1 in (F(1), F(1, 3), F(7, 2)):
    for k2 in (F(1), F(5, 4), F(1, 9)):
        for T in (F(1, 2), F(2), F(9)):
            vs = min(T, k2 / (2 * k1 * T ** 3))
            okb &= k1 * T ** 3 * vs <= k2 / 2 and vs > 0
claim("E2.1", okb, "k1 T^3 v* <= k2/2 with v* = min(T, k2/(2 k1 T^3)) > 0 on a grid of (k1, k2, T)")
Dn = CNet(2, [((3, 2), (2, 3)), ((0, 1), (1, 0))], names=["U", "V"])
claim("E2.2", endo.drain_species(2, Dn.reactions) == [1] and H(Dn) and not Dn.weakly_reversible(),
      "V is a drain species of D (source V of V->U, producing reaction has V-coefficient 2 > 1); D has (H1)-(H3)")
Mn = Dn.dual()
claim("E2.3", Mn.weakly_reversible() and Mn.linkage_classes() == 1 and Mn.deficiency() == 0 and H(Mn)
      and endo.drain_species(2, Mn.reactions) == [], f"M = D* = {Mn.describe()}: reversible, one linkage class, "
      "deficiency 0, (H1)-(H3), no drain species")

# E3: witnesses (filled from the test logs)
WITNESSES = [
    (2, [((3, 0), (2, 1)), ((2, 1), (3, 0))]),  # G's M (in the family (2,3,6))
    (3, [((0, 0, 3), (0, 1, 0)), ((0, 1, 0), (1, 0, 2)), ((1, 0, 2), (0, 0, 3))]),  # smallest n = 3, s = 2 (3,3,5)
    (4, [((0, 0, 0, 1), (0, 0, 3, 0)), ((0, 0, 3, 0), (1, 2, 0, 0)), ((0, 3, 0, 0), (0, 0, 0, 1)),
         ((1, 2, 0, 0), (0, 3, 0, 0))]),  # s = 3: X4 -> 3X3 -> X1+2X2 -> 3X2 -> X4 (family (4,3,4))
]
for n, rx in WITNESSES:
    N = CNet(n, rx)
    Dd = N.dual()
    claim("E3", N.weakly_reversible() and H(N) and N.deficiency() == 0 and N.linkage_classes() == 1 and H(Dd)
          and endo.drain_species(Dd.n, Dd.reactions) != [],
          f"N = {N.describe()} WR with (H1)-(H3); N* = {Dd.describe()} has drain species "
          f"{endo.drain_species(Dd.n, Dd.reactions)}; s = {N.rank()}")
print(f"# runtime {time.time() - t0:.1f} s")
sys.exit(0 if OK else 1)
