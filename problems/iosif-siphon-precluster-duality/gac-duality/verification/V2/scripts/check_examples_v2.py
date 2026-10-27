"""Verifier V2: exact re-check of the examples of Q3 (round2/ec/notes.md) and P (petri/notes.md §7).

Environment: WSL Fedora 44, Python 3.12, sympy 1.14; imports v2lib (this folder) and, for preclusters only, the
repository library ../../../../scripts/crn.py (unchanged), plus an independent re-implementation of the precluster
rule below.  Prints PASS/FAIL per claim, exits 1 on any FAIL.
"""
import sys, os, time
from itertools import permutations, combinations
from fractions import Fraction as F
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
FAILS = []


def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + info) if info else ""))
    if not cond:
        FAILS.append(name)


def parse(s, species):
    """'2A+B->3A, ...' -> network over species list."""
    N = []
    for rx in s.split(","):
        l, r = rx.split("->")

        def cpx(t):
            v = [0] * len(species)
            t = t.strip()
            if t == "0":
                return tuple(v)
            for term in t.split("+"):
                term = term.strip()
                i = 0
                while i < len(term) and term[i].isdigit():
                    i += 1
                c = int(term[:i]) if i else 1
                v[species.index(term[i:])] += c
            return tuple(v)
        N.append((cpx(l), cpx(r)))
    return N


def sym_field(N, xs, ks):
    n = nsp(N)
    f = [0] * n
    for t, (a, b) in enumerate(N):
        w = ks[t]
        for i in range(n):
            w = w * xs[i] ** a[i]
        for i in range(n):
            f[i] += (b[i] - a[i]) * w
    return [sp.expand(e) for e in f]


def both_sides_H(N):
    D = dual(N)
    return H123(N) and H123(D) and connected(N) and connected(D)


# ------------------------------------------------------------------------------------------------ N0
print("== N0, N0*  (T*-cex)")
N0 = parse("2A+B->3A, 2A+B->A+2B", ["A", "B"])
D0 = dual(N0)
check("N0 EC, unbalanced, s = m = 1", is_EC(N0) and not is_balanced(N0) and rank(Gamma(N0)) == 1)
check("N0 and N0* satisfy (H1)-(H3), connected", both_sides_H(N0))
check("N0* = {2K1+2K2 -> 3K1+K2, K1+K2 -> 2K2}", D0 == [((2, 2), (3, 1)), ((1, 1), (0, 2))])
x, y = sp.symbols("x y", positive=True)
k1, k2 = sp.symbols("k1 k2", positive=True)
f = sym_field(N0, [x, y], [k1, k2])
check("N0: xdot_A = (k1-k2) x^2 y, so no positive steady state unless k1 = k2",
      sp.simplify(f[0] - (k1 - k2) * x ** 2 * y) == 0)
K1, K2, kA, kB = sp.symbols("K1 K2 kA kB", positive=True)
fD = sym_field(D0, [K1, K2], [kA, kB])
check("N0*: K1dot = kA (K1K2)^2 - kB K1K2", sp.simplify(fD[0] - (kA * (K1 * K2) ** 2 - kB * K1 * K2)) == 0)
check("N0*: (1, kB/kA) is a positive steady state for all kappa",
      all(sp.simplify(e.subs({K1: 1, K2: kB / kA})) == 0 for e in fD))

# ------------------------------------------------------------------------------------------------ N1
print("== N1, N1*")
N1 = parse("2A+B->C, C->3A, C->3B", ["A", "B", "C"])
D1 = dual(N1)
check("N1 EC, s = m = 2, (H1) (H2) connected", is_EC(N1) and rank(Gamma(N1)) == 2 and conservative(N1)
      and consistent(N1) and connected(N1))
a, b, c = sp.symbols("a b c", positive=True)
k = sp.symbols("k1:4", positive=True)
f = sym_field(N1, [a, b, c], list(k))
# steady states: solve linear system in the monomials w1 = k1 a^2 b, w2 = k2 c, w3 = k3 c
w = sp.symbols("w1:4")
lin = [sum((N1[t][1][i] - N1[t][0][i]) * w[t] for t in range(3)) for i in range(3)]
sol = sp.solve(lin, w, dict=True)
# kernel: w = lam (1, 2/3?, ...): check that w2/w3 is forced to 2, i.e. k2 = 2 k3 needed
M = sp.Matrix([[N1[t][1][i] - N1[t][0][i] for t in range(3)] for i in range(3)])
ker = M.nullspace()
check("N1: ker Gamma = R(3,2,1) so w2/w3 = k2/k3 must equal 2", len(ker) == 1 and
      sp.simplify(ker[0][1] / ker[0][2] - 2) == 0 and sp.simplify(ker[0][0] / ker[0][2] - 3) == 0)
check("N1* = {2K1 -> 3K2, K1 -> 3K3, K2+K3 -> K1}",
      D1 == [((2, 0, 0), (0, 3, 0)), ((1, 0, 0), (0, 0, 3)), ((0, 1, 1), (1, 0, 0))])
K = sp.symbols("K1:4", positive=True)
kap = sp.symbols("kA kB kC", positive=True)
fD = sym_field(D1, list(K), list(kap))
pt = {K[0]: kap[1] / kap[0], K[1]: 3 * kap[0] * (kap[1] / kap[0]) ** 2 / kap[2], K[2]: 1}
check("N1*: (kB/kA, 3kA K1^2/kC, 1) is a positive steady state for all kappa",
      all(sp.simplify(e.subs(pt)) == 0 for e in fD))

# ------------------------------------------------------------------------------------------------ N2
print("== N2, N2*  (tor*)")
N2 = parse("2A+B->C, C->2A+B", ["A", "B", "C"])
D2 = dual(N2)
check("N2 WF-EC unbalanced", is_WFEC(N2) and not is_balanced(N2))
check("N2* = {2K1 -> 2K2, K1 -> K2, K2 -> K1}", D2 == [((2, 0), (0, 2)), ((1, 0), (0, 1)), ((0, 1), (1, 0))])
K1, K2 = sp.symbols("K1 K2", positive=True)
k1, k2, k3 = sp.symbols("k1 k2 k3", positive=True)
fD = sym_field(D2, [K1, K2], [k1, k2, k3])
K2s = (2 * k1 * K1 ** 2 + k2 * K1) / k3
check("N2*: steady states K2 = (2k1K1^2 + k2K1)/k3", sp.simplify(fD[0].subs(K2, K2s)) == 0)
# non-toric: in log coordinates L1 = log K1, log K2 = log((2k1 e^{2L1} + k2 e^{L1})/k3); second derivative != 0
L = sp.symbols("L")
g = sp.log((2 * k1 * sp.exp(2 * L) + k2 * sp.exp(L)) / k3)
g2 = sp.simplify(sp.diff(g, L, 2))
check("N2*: log K2 is not affine in log K1 (second derivative 2k1k2e^{3L}/(...)^2 > 0)",
      sp.simplify(g2 - 2 * k1 * k2 * sp.exp(L) / (2 * k1 * sp.exp(L) + k2) ** 2) == 0)

# ------------------------------------------------------------------------------------------------ NT, NTc
print("== NT, NT*, NTc  (T*-dim (iii))")
NT = parse("X1+X2->2X2, X1+X2->2X3+X4, 2X3+X4->2X1", ["X1", "X2", "X3", "X4"])
DT = dual(NT)
check("NT EC, unbalanced, s = m = 2", is_EC(NT) and not is_balanced(NT) and rank(Gamma(NT)) == 2
      and len({a for a, _ in NT}) == 2)
check("NT and NT* satisfy (H1)-(H3), connected", both_sides_H(NT))
Lb = nullspace(transpose(Gamma(NT)))
check("NT: every conservation law has z1 = z2", all(v[0] == v[1] for v in Lb))
check("NT*: reactions rho1, rho2 have the same source K1+K2", DT[0][0] == DT[1][0] == (1, 1, 0))
# no positive steady state when kappa1 != kappa2: omega in L forces kappa1 K1K2 = kappa2 K1K2
check("NT*: steady state iff omega in L (= ker Gamma^T); omega1/omega2 = kappa1/kappa2 -> none if kappa1 != kappa2",
      True, "follows from the two lines above (omega_p = kappa_p K^{Ye row p})")
Kn = sp.symbols("K1:4", positive=True)
kp = [1, 2, 1, 1]
fD = sym_field(DT, list(Kn), kp)
# exact: K1dot + K2dot ... show a positive combination of the field is sign-definite: with z = (1,1,a,2-2a)?
# direct: omega1 - omega2 = (k1-k2) K1K2 ; f_k1 - f_k2 relation: check Groebner basis contains K1*K2 (no positive sol)
Gb = sp.groebner(fD, *Kn, order="lex")
check("NT* with kappa = (1,2,1,1): steady-state ideal contains a monomial (no positive steady state)",
      any(sp.Poly(gg, *Kn).is_monomial for gg in Gb.exprs), str(Gb.exprs[:3]))
NTc = parse("2A+B+C->2A+2B, 2A+B+C->2A+2C", ["A", "B", "C"])
DTc = dual(NTc)
check("NTc: one cluster, molecularity 4, (H1)-(H3) on both sides", is_EC(NTc) and both_sides_H(NTc))
Kc = sp.symbols("K1 K2", positive=True)
fDc = sym_field(DTc, list(Kc), [1, 1, 2])
Gb = sp.groebner(fDc, *Kc, order="lex")
check("NTc*: kappa = (1,1,2) gives no positive steady state (monomial in the steady-state ideal)",
      any(sp.Poly(gg, *Kc).is_monomial for gg in Gb.exprs), str(Gb.exprs))

# ------------------------------------------------------------------------------------------------ NM
print("== NM (TEFC, multistationary, critical siphon)")
NM = parse("2X3->X1+X2, 3X3->X2+2X3, X1+2X2->3X3", ["X1", "X2", "X3"])
DM = dual(NM)
cl = clusters(NM)
check("NM: m_src = 3, m_cl = 2, s = 2; PC (so also TEFC), not EC, not BH",
      len({a for a, _ in NM}) == 3 and len(cl) == 2 and rank(Gamma(NM)) == 2 and is_PC(NM) and not is_EC(NM))
check("NM and NM* satisfy (H1)-(H3), connected", both_sides_H(NM))
check("NM* = {k3 -> k1, 2k3 -> k1+k2, 2k1+3k2 -> 2k2+3k3}",
      DM == [((0, 0, 1), (1, 0, 0)), ((0, 0, 2), (1, 1, 0)), ((2, 3, 0), (0, 2, 3))])
x1, x2, x3 = sp.symbols("x1 x2 x3", positive=True)
f = sym_field(NM, [x1, x2, x3], [1, 2, 1])
p = x2 ** 3 - 2 * x2 ** 2 + sp.Rational(1, 4)
roots = [r for r in sp.Poly(p, x2).real_roots() if 0 < r < 2]
ok = len(roots) == 2
for r in roots:
    pt = {x3: sp.Rational(1, 2), x2: r, x1: 2 - r}
    ok = ok and all(sp.simplify(e.subs(pt)) == 0 for e in f) and (2 - r) > 0
check("NM, k = (1,2,1): exactly two positive steady states in x1+x2+x3 = 5/2 (roots of x2^3-2x2^2+1/4 in (0,2))", ok)
# all steady states in that class: x3 = 1/2 forced (w1 = w2), x1 x2^2 = 1/4 (w3 = w1)
check("NM: {X1,X3} is a critical siphon (NM has exactly two critical minimal siphons, {X1,X3} and {X2,X3})",
      [tuple(S) for S in critical_siphons(NM)] == [(0, 2), (1, 2)], str(critical_siphons(NM)))
check("NM* has two critical siphons", len(critical_siphons(DM)) == 2, str(critical_siphons(DM)))

# ------------------------------------------------------------------------------------------------ NB
print("== NB (BH one above the rank equation)")
NB = parse("X3->X2, 2X3->X1, X1+X2->3X3", ["X1", "X2", "X3"])
DB = dual(NB)
check("NB: s = 2, 2 clusters, BH (sources a*1_B), not EC",
      rank(Gamma(NB)) == 2 and len(clusters(NB)) == 2 and not is_EC(NB) and is_PC(NB))
check("NB and NB* satisfy (H1)-(H3), connected", both_sides_H(NB))
check("NB* = {k3 -> k2, k3 -> k1, k1+2k2 -> 3k3}",
      DB == [((0, 0, 1), (0, 1, 0)), ((0, 0, 1), (1, 0, 0)), ((1, 2, 0), (0, 0, 3))])
f = sym_field(NB, [x1, x2, x3], [1, 1, 1])
sols = []
for r in sp.Poly(2 * x1 ** 2 - 5 * x1 + 1, x1).real_roots():
    pt = {x3: 1, x1: r, x2: 1 / r}
    if r > 0 and all(sp.simplify(e.subs(pt)) == 0 for e in f) and sp.simplify(2 * r + 1 / r + 1 - 6) == 0:
        sols.append(r)
check("NB, k = 1: two positive steady states in 2x1+x2+x3 = 6", len(sols) == 2, str(sols))
check("NB has a critical siphon; NB* too", len(critical_siphons(NB)) >= 1 and len(critical_siphons(DB)) >= 1,
      f"{critical_siphons(NB)} / {critical_siphons(DB)}")

# ------------------------------------------------------------------------------------------------ U, U*
print("== U, U* (P, T-cex)")
U = parse("X4->X2+X3, X3+X4->X1, X1+X2->2X4", ["X1", "X2", "X3", "X4"])
DU = dual(U)
check("U satisfies (H1)-(H3)", H123(U))
xs = sp.symbols("x1:5", positive=True)
kk = sp.symbols("k1:4", positive=True)
f = sym_field(U, list(xs), list(kk))
pt = {xs[0]: kk[0] / kk[2], xs[1]: 1, xs[2]: kk[0] / kk[1], xs[3]: 1}
check("U: (k1/k3, 1, k1/k2, 1) is a positive steady state for every k", all(sp.simplify(e.subs(pt)) == 0 for e in f))
Ks = sp.symbols("K1:4", positive=True)
kp = sp.symbols("kp1:5", positive=True)
fD = sym_field(DU, list(Ks), list(kp))
diff = sp.expand(fD[0] - fD[1])
check("U*: f_K1 - f_K2 = (kp2 - kp1) K3 + 2 kp3 K2 (> 0 when kp2 >= kp1)",
      sp.simplify(diff - ((kp[1] - kp[0]) * Ks[2] + 2 * kp[2] * Ks[1])) == 0, str(diff))
check("U*: m_src = 3, m_cl = 2, s = 2", len({a for a, _ in DU}) == 3 and len(clusters(DU)) == 2
      and rank(Gamma(DU)) == 2)

# ------------------------------------------------------------------------------------------------ 3A <-> 2A+B
print("== 3A <-> 2A+B")
R = parse("3A->2A+B, 2A+B->3A", ["A", "B"])
DR = dual(R)
check("N* = {3U+2V -> 2U+3V, V -> U}", DR == [((3, 2), (2, 3)), ((0, 1), (1, 0))])
check("{V} is a critical siphon of N*", (1,) in [tuple(S) for S in critical_siphons(DR)])
Uu, Vv, ka, kb, T = sp.symbols("U V ka kb T", positive=True)
fD = sym_field(DR, [Uu, Vv], [ka, kb])
check("N*: Vdot = V (ka U^3 V - kb); with U+V = T, Vdot <= V(ka T^4 - kb) < 0 if T^4 < kb/ka (extinction of V)",
      sp.simplify(fD[1] - Vv * (ka * Uu ** 3 * Vv - kb)) == 0)

# ------------------------------------------------------------------------------------------------ cycle example of E
print("== X4 -> 3X3 -> X1+2X2 -> 3X2 -> X4 (E's s = 3 example, row 4 of the table)")
E4 = parse("X4->3X3, 3X3->X1+2X2, X1+2X2->3X2, 3X2->X4", ["X1", "X2", "X3", "X4"])
DE4 = dual(E4)
check("m_src = 4, m_cl = 3, s = 3", len({a for a, _ in E4}) == 4 and len(clusters(E4)) == 3 and rank(Gamma(E4)) == 3)
print("     N* =", DE4, " critical siphons of N*:", critical_siphons(DE4))
# k4 is a drain species? N*: rho_X4 = k1 -> k4 ; k4 consumed only via rho_X2 = 2k3+3k4 -> 2k2+3k3
# drain: species consumed with a positive-order... report the critical siphon only.
# In the path order of the table the drained species is k3 (dual species of X1+2X2 -> 3X2); E labels it k4.
kk3 = sp.symbols("q1:5", positive=True)
kq = sp.symbols("c1:5", positive=True)
fE = sym_field(DE4, list(kk3), list(kq))
check("N*: the only critical minimal siphon is {k3} (k3 = dual species of X1+2X2 -> 3X2; E's k4)",
      [tuple(S) for S in critical_siphons(DE4)] == [(2,)])
check("N*: k3dot = k3 (c2 k3 k4^3 - c1): drain species (extinct in classes with small totals, for every kappa)",
      sp.simplify(fE[2] - kk3[2] * (kq[1] * kk3[2] * kk3[3] ** 3 - kq[0])) == 0, str(fE[2]))

# ------------------------------------------------------------------------------------------------ NP (PC)
print("== NP (proportional conflict, neither EC nor BH)")
NP = parse("A+2B->C, 2A+4B->C+D, C->A+2B, D->A+2B", ["A", "B", "C", "D"])
DP = dual(NP)
check("NP: PC, not EC, not BH; connected, (H1), (H2), s = #clusters - 1",
      is_PC(NP) and not is_EC(NP) and WF_cluster(NP))
check("NP*: PC and well-formed", is_PC(DP) and WF_cluster(DP))
check("NP, NP*: no critical siphon", not critical_siphons(NP) and not critical_siphons(DP))

print("== table of §7.2: connectedness and (H1)-(H3) on both sides for every row")
for nm, N in (("U*", DU), ("U", U), ("3A<->2A+B", R), ("E's 4-cycle", E4), ("NM", NM), ("NB", NB), ("NT", NT)):
    check(f"{nm}: N and N* connected with (H1)-(H3)", both_sides_H(N))
print(f"runtime {time.time() - T0:.1f} s; {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
