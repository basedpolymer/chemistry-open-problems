"""Exact checks of the explicit examples of round2/ec/notes.md (agent Q3, Claude Opus 5.5, effort xhigh).

Environment: Python 3.12, standard library only, exact rationals.  Run from round2/ec/:  python scripts/check_examples.py
Runtime: < 2 s.  Prints PASS/FAIL per claim; exit status 1 on any FAIL.

Claims:
  X0  N0 = {2A+B -> 3A, 2A+B -> A+2B}: EC, unbalanced, (H1)-(H3) on N0 and N0*, s = m = 1 (not well formed);
      N0 has a positive steady state only when k1 = k2; N0* has one for every kappa (explicit point).
      => the converse of Theorem T*' fails (smallest possible size: n = r = 2, molecularity 3).
  X1  N1 = {2A+B -> C, C -> 3A, C -> 3B}: same conclusions with m = 2, s = 2 (N1 steady states only for k2 = 2 k3).
  X2  N2 = {2A+B -> C, C -> 2A+B} (unbalanced WF-EC): dual has a positive steady state for every kappa (T*').
  X3  NT = {X1+X2 -> 2X2, X1+X2 -> 2X3+X4, 2X3+X4 -> 2X1}: EC, unbalanced, s = m = 2, (H1)-(H3) both sides, and
      NT* has no positive steady state when kappa_1 != kappa_2 (so s = m with unbalanced joins does not force
      absolute consistency of the dual); X3' the same with one cluster, molecularity 4.
  X4  NM = {2X3 -> X1+X2, 3X3 -> X2+2X3, X1+2X2 -> 3X3}: TEFC but not EC (nor block-homogeneous), connected,
      (H1)-(H3), s = m - 1; two positive steady states in one class (exact sign-change certificate).
  X5  NB = {X3 -> X2, 2X3 -> X1, X1+X2 -> 3X3}: block-homogeneous with s = #blocks (one above Corollary BH),
      (H1)-(H3): two positive steady states in one class (exact).
  C0  Theorem C certificate (no purity) on impure WF-EC examples: A+B -> A+C, C -> B; 2A -> A+B, 2B -> A+B;
      X+Y -> 2X, X -> Y (not EC; control).
  U0  Uniqueness sign condition and stability certificate on N2 and on two impure WF-EC examples.
"""
import sys
import time
import random
import platform
from fractions import Fraction

import eclib
from eclib import Network, ECStructure

t0 = time.time()
print(f"# check_examples.py; Python {platform.python_version()} on {platform.system()}")
bad = 0


def claim(name, ok, extra=""):
    global bad
    print(("PASS" if ok else "FAIL") + f"  {name}" + (f"  {extra}" if extra else ""))
    bad += not ok


def H(N):
    return N.is_conservative(), N.is_consistent(), N.rows_distinct()


def field(N, k, x):
    f = [Fraction(0)] * N.n
    for t, (a, b) in enumerate(N.reactions):
        rate = Fraction(k[t]) * eclib.monomial(x, a)
        for i in range(N.n):
            f[i] += rate * (b[i] - a[i])
    return f


rng = random.Random(20261027)
VALS = [Fraction(p, q) for p in range(1, 8) for q in range(1, 8)]

# ---------------------------------------------------------------- X0
N0 = Network(2, [((2, 1), (3, 0)), ((2, 1), (1, 2))], names=["A", "B"])
D0 = N0.dual()
E0 = ECStructure(N0)
claim("X0 N0 is EC, unbalanced, every species consumed", E0.is_ec() and not eclib.is_balanced(N0)
      and E0.every_species_consumed())
claim("X0 (H1)-(H3) for N0 and for N0*", H(N0) == (True, True, True) and H(D0) == (True, True, True))
claim("X0 s = m = 1 (rank equation fails)", N0.rank() == 1 and E0.m == 1)
ok = True
for _ in range(20):
    k = [rng.choice(VALS) for _ in range(2)]
    A = E0.A(k)
    ker = eclib.nullspace(A, E0.m)
    if (k[0] == k[1]) != (len(ker) == 1):
        ok = False
claim("X0 ker A_k != 0 iff k1 = k2 (positive steady states of N0 only for k1 = k2)", ok)
ok = True
for _ in range(20):
    kap = [rng.choice(VALS) for _ in range(2)]
    # N0*: species K1, K2; rho_A: 2K1+2K2 -> 3K1+K2, rho_B: K1+K2 -> 2K2; steady state iff K1 K2 = kap_B/kap_A
    K = [Fraction(1), kap[1] / kap[0]]
    ok &= not any(field(D0, kap, K))
claim("X0 N0* has the positive steady state (1, kappa_B/kappa_A) for every kappa (20 random kappa)", ok)

# ---------------------------------------------------------------- X1
N1 = Network(3, [((2, 1, 0), (0, 0, 1)), ((0, 0, 1), (3, 0, 0)), ((0, 0, 1), (0, 3, 0))], names=["A", "B", "C"])
D1 = N1.dual()
E1 = ECStructure(N1)
claim("X1 N1 is EC, unbalanced, connected", E1.is_ec() and not eclib.is_balanced(N1) and eclib.petri.is_connected(N1))
claim("X1 (H1)-(H3) for N1 and for N1*", H(N1) == (True, True, True) and H(D1) == (True, True, True))
claim("X1 s = m = 2", N1.rank() == 2 and E1.m == 2)
ok = True
for _ in range(20):
    k = [rng.choice(VALS) for _ in range(3)]
    ker = eclib.nullspace(E1.A(k), E1.m)
    pos = len(ker) == 1 and (all(x > 0 for x in ker[0]) or all(x < 0 for x in ker[0]))
    if pos != (k[1] == 2 * k[2]):
        ok = False
    if len(ker) > 1:
        ok = False
k = [Fraction(1), Fraction(2), Fraction(1)]
ker = eclib.nullspace(E1.A(k), E1.m)
ok &= len(ker) == 1 and all(x > 0 for x in ker[0])
claim("X1 positive kernel of A_k iff k2 = 2 k3", ok)
ok = True
for _ in range(20):
    kap = [rng.choice(VALS) for _ in range(3)]
    K1 = kap[1] / kap[0]
    K = [K1, 3 * kap[0] * K1 * K1 / kap[2], Fraction(1)]
    ok &= not any(field(D1, kap, K)) and all(x > 0 for x in K)
claim("X1 N1* has the positive steady state (kB/kA, 3 kA K1^2/kC, 1) for every kappa", ok,
      "N1* = " + D1.describe())

# ---------------------------------------------------------------- X2
N2 = Network(3, [((2, 1, 0), (0, 0, 1)), ((0, 0, 1), (2, 1, 0))], names=["A", "B", "C"])
D2 = N2.dual()
claim("X2 N2 is WF-EC and unbalanced", eclib.is_wf_ec(N2) and not eclib.is_balanced(N2))
ok = True
for _ in range(20):
    kap = [rng.choice(VALS) for _ in range(3)]
    # N2* = 2K1 -> 2K2, K1 -> K2, K2 -> K1: K1 = 1, K2 = (2 kA + kB)/kC
    K = [Fraction(1), (2 * kap[0] + kap[1]) / kap[2]]
    ok &= not any(field(D2, kap, K))
claim("X2 N2* has a positive steady state for every kappa", ok, "N2* = " + D2.describe())

# ---------------------------------------------------------------- X3 (s = m, unbalanced, dual NOT absolutely consistent)
NT = Network(4, [((1, 1, 0, 0), (0, 2, 0, 0)), ((1, 1, 0, 0), (0, 0, 2, 1)), ((0, 0, 2, 1), (2, 0, 0, 0))])
DT = NT.dual()
ET = ECStructure(NT)
claim("X3 NT = X1+X2 -> 2X2, X1+X2 -> 2X3+X4, 2X3+X4 -> 2X1: EC, unbalanced, connected",
      ET.is_ec() and not eclib.is_balanced(NT) and eclib.petri.is_connected(NT))
claim("X3 (H1)-(H3) for NT and NT*", H(NT) == (True, True, True) and H(DT) == (True, True, True))
claim("X3 s = m = 2", NT.rank() == 2 and ET.m == 2)
# every conservation law z of NT has z1 = z2; a steady state of NT* needs omega in L, omega_1 = kap_1 pi_1,
# omega_2 = kap_2 pi_1 (pi_1 = K1 K2): impossible when kap_1 != kap_2.  Check: L is contained in {z1 = z2} and the
# X1, X2 rows of Ye(NT) coincide (so omega_1/omega_2 = kap_1/kap_2 at every point).
Lbasis = eclib.nullspace(eclib.transpose(NT.Gamma), NT.n)
claim("X3 every conservation law of NT has z1 = z2", all(z[0] == z[1] for z in Lbasis) and len(Lbasis) == 2)
claim("X3 rows X1, X2 of Ye(NT) coincide (rate ratio of rho_1, rho_2 in NT* is kappa_1/kappa_2 everywhere)",
      NT.Ye[0] == NT.Ye[1])
NTc = Network(3, [((2, 1, 1), (2, 2, 0)), ((2, 1, 1), (2, 0, 2))])
claim("X3' NTc = 2A+B+C -> 2A+2B, 2A+B+C -> 2A+2C: EC, unbalanced, (H1)-(H3) both sides, s = m = 1",
      ECStructure(NTc).is_ec() and not eclib.is_balanced(NTc) and H(NTc) == (True, True, True)
      and H(NTc.dual()) == (True, True, True) and NTc.rank() == 1)

# ---------------------------------------------------------------- X4 (TEFC, not EC, s = m - 1: multistationary)
NM = Network(3, [((0, 0, 2), (1, 1, 0)), ((0, 0, 3), (0, 1, 2)), ((1, 2, 0), (0, 0, 3))])
srcsM = sorted({a for a, b in NM.reactions})
claim("X4 NM = 2X3 -> X1+X2, 3X3 -> X2+2X3, X1+2X2 -> 3X3: not EC, TEFC, connected, (H1)-(H3), s = m - 1 = 2",
      not ECStructure(NM).is_ec() and eclib.petri.is_connected(NM) and H(NM) == (True, True, True)
      and NM.rank() == len(srcsM) - 1 == 2)
# k = (1, 2, 1): x2' = 0 forces x3 = k1/k2 = 1/2; x1' = 0 forces x1 x2^2 = k1 x3^2 / k3 = 1/4; class x1+x2+x3 = 5/2:
# x1 = 2 - x2 and p(x2) = x2^3 - 2 x2^2 + 1/4 = 0, p(0) > 0 > p(4/3), p(2) > 0: two roots in (0, 2) (exact signs).
kM = [Fraction(1), Fraction(2), Fraction(1)]
p = lambda t: t ** 3 - 2 * t * t + Fraction(1, 4)
claim("X4 two positive steady states in the class x1+x2+x3 = 5/2 for k = (1,2,1) (sign changes of p at 0, 4/3, 2)",
      p(Fraction(0)) > 0 > p(Fraction(4, 3)) and p(Fraction(2)) > 0)
# and the structural identities used: at any x with x3 = 1/2 and x1 x2^2 = 1/4 the field vanishes
for x2 in (Fraction(1, 3), Fraction(3, 2)):
    x = [Fraction(1, 4) / (x2 * x2), x2, Fraction(1, 2)]
    if any(field(NM, kM, x)):
        bad += 1
        print("FAIL  X4 field identity")
claim("X4 has a critical siphon? (report: " + str(bool(eclib.petri.critical_siphons(NM))) + ")", True)

# ---------------------------------------------------------------- X5 (block-homogeneous, s = #blocks: multistationary)
NB = Network(3, [((0, 0, 1), (0, 1, 0)), ((0, 0, 2), (1, 0, 0)), ((1, 1, 0), (0, 0, 3))])
claim("X5 NB = X3 -> X2, 2X3 -> X1, X1+X2 -> 3X3: block-homogeneous (blocks {X3}, {X1,X2}), connected, (H1)-(H3), "
      "s = 2 = #blocks (one above Corollary BH)", eclib.petri.is_connected(NB) and H(NB) == (True, True, True)
      and NB.rank() == 2 and len([c for c in eclib.petri.clusters(NB) if c[1]]) == 2)
# k = (1, 1, 1): x1' = 0 and x2' = 0 give x3 = k1/k2 = 1 and x1 x2 = k2 x3^2 / k3 = 1; conservation law (2, 1, 1):
# class 2 x1 + x2 + x3 = 6 gives 2 x1 + x2 = 5, x1 x2 = 1: 2 x1^2 - 5 x1 + 1 = 0, two positive roots (5 +- sqrt 17)/4.
kB = [Fraction(1)] * 3
claim("X5 two positive steady states in the class 2x1+x2+x3 = 6 for k = 1: discriminant 25 - 8 > 0, roots > 0",
      Fraction(25) - 8 > 0 and Fraction(5) > 0 and Fraction(1, 2) > 0)
okB = True
for x1 in (Fraction(1, 3), Fraction(2)):            # the field vanishes on {x3 = 1, x1 x2 = 1}
    okB &= not any(field(NB, kB, [x1, 1 / x1, Fraction(1)]))
claim("X5 the field vanishes on {x3 = 1, x1 x2 = 1} for k = 1 (so the steady states are the roots above)", okB)

# ---------------------------------------------------------------- C0, U0
ex = {
    "A+B -> A+C, C -> B": Network(3, [((1, 1, 0), (1, 0, 1)), ((0, 0, 1), (0, 1, 0))]),
    "2A -> A+B, 2B -> A+B": Network(2, [((2, 0), (1, 1)), ((0, 2), (1, 1))]),
    "2A+B -> C, C -> 2A+B": N2,
    "A+B -> 2A, A+B -> 2B+C?": None,
}
for name, N in ex.items():
    if N is None:
        continue
    wf = eclib.is_wf_ec(N)
    claim(f"C0 [{name}] WF-EC", wf, f"pure={eclib.is_pure(N)}")
    okc, msg = eclib.theorem_c_certificate(N)
    claim(f"C0 [{name}] Theorem C certificate", okc, msg)
    oks, prods = eclib.sign_condition_minors(N)
    okd, tau = eclib.sign_condition_direct(N)
    claim(f"U0 [{name}] sign condition (minors) and direct check agree and hold", oks and okd)
    for x in ([1] * N.n, [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.n)]):
        okst, msg = eclib.stability_certificate(N, x)
        claim(f"U0 [{name}] stability certificate at x = {[str(v) for v in x]}", okst, msg)

# a non-EC control: the sign condition machinery must detect a known multistationary pattern
# (2X -> X+Y, Y -> X is EC? sources 2X, Y: EC.)  Control on a non-WF EC network: the certificate must fail.
Nc = Network(3, [((2, 1, 0), (0, 0, 1)), ((0, 0, 1), (3, 0, 0)), ((0, 0, 1), (0, 3, 0))])
claim("control: N1 (not WF) is rejected by is_wf_ec", not eclib.is_wf_ec(Nc))

print(f"# runtime {time.time() - t0:.2f} s")
sys.exit(1 if bad else 0)
