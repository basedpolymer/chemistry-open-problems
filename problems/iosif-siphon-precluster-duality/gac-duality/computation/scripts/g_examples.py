"""Red team of agent G's explicit examples and finite-family counts (gac/notes.md, sections 3, 5, 10).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (fractions.Fraction;
crn.py via ecore.py).  Runtime: a few seconds.

Claims checked (PASS/FAIL):
  G3a  M: 3A -> 2A+B, 2A+B -> 3A and D = M* = 3U+2V -> 2U+3V, V -> U: both Gamma = [[-1,1],[1,-1]], both satisfy
       (H1)-(H3); M reversible, deficiency 0, one linkage class; D has 4 complexes, 2 linkage classes, deficiency 1,
       is not WR; direction w = (0,1) witnesses that D is not endotactic.
  G3b  with unit rate constants, D: v' = u^3 v^2 - v; on u + v = 2, h(v) = (2-v)^3 v satisfies h' = (2-v)^2 (2-4v)
       (polynomial identity), h(1/8) = 3375/4096 < 1, h(1) = 1 (positive equilibrium (1,1) in that class), so
       v' <= -(721/4096) v on (0, 1/8]: extinction of V from (15/8, 1/8).  M: a' = a^2 (b - a).
  G3c  the class u + v = 2 of D has exactly two positive equilibria (h = 1 has one root in (1/8, 1/2) and the root 1)
       (an observation of E: the trajectory from (15/8,1/8) lies below the smaller one).
  G5   H: 2A <-> 2B, 2A <-> A+B: WR, (H1)-(H3), only minimal siphon {A,B} = support of a+b (ADLS holds); its dual
       has (H1)-(H3) and the siphon {U4} supports no conservation law (ADLS fails).
  G4.1 Prop. 4.1 examples: the 5-species two-cycle example (s=4, q=1, deficiency 4, dual one 5-cycle, critical
       siphon {X1,X2,X3}, strongly endotactic) and the 4-species example showing primality is needed.
  G8.1 Example 8.1: Gamma, J_S = [[-15,27],[24,-25]], det -273, alternative realisation with the same Gamma.
  G10  the 705 labelled networks of section 10 (3 species, the six complexes of molecularity two; unions of 1-3
       reversible pairs: 575; directed 3- and 4-cycles with rotation fixed: 130): 702 satisfy (H3); among them the
       ADLS criterion holds on both sides for 18, only for N for 177, only for N* for 0, for neither 507; 13 have a
       WR dual; and q = beta* + delta* on all 702.
"""
import sys
import time
import platform
from fractions import Fraction as F
from itertools import combinations, permutations

import ecore
from ecore import CNet, siphon_masks, minimal_masks, to_set, semiflow_supports

t0 = time.time()
print(f"# g_examples.py  Python {platform.python_version()} on {platform.system()} {platform.release()}")
OK = True


def claim(tag, cond, doc):
    global OK
    OK &= bool(cond)
    print(f"{'PASS' if cond else 'FAIL'} {tag}: {doc}")


def H(N):
    return N.is_conservative() and N.is_consistent() and N.rows_distinct()


def adls(N):
    sf = [sum(1 << i for i in S) for S in semiflow_supports(N)]
    ms = minimal_masks(siphon_masks(N))
    return all(any((f & S) == f for f in sf) for S in ms)


# ---------------- G3
M = CNet(2, [((3, 0), (2, 1)), ((2, 1), (3, 0))], names=["A", "B"])
D = M.dual()
DD = CNet(2, [((3, 2), (2, 3)), ((0, 1), (1, 0))], names=["U", "V"])
claim("G3a.1", D.reactions == DD.reactions, f"dual of M is {DD.describe()}")
claim("G3a.2", M.Gamma == [[-1, 1], [1, -1]] and D.Gamma == [[-1, 1], [1, -1]], "both Gamma = [[-1,1],[1,-1]]")
claim("G3a.3", H(M) and H(D), "M and D satisfy (H1)-(H3)")
claim("G3a.4", M.weakly_reversible() and M.deficiency() == 0 and M.linkage_classes() == 1,
      "M weakly reversible, deficiency 0, one linkage class")
claim("G3a.5", len(D.complexes()) == 4 and D.linkage_classes() == 2 and D.deficiency() == 1
      and not D.weakly_reversible(), "D: 4 complexes, 2 linkage classes, deficiency 1, not WR")
w = (0, 1)
ess = [(a, b) for a, b in D.reactions if sum(wi * (bi - ai) for wi, ai, bi in zip(w, a, b)) != 0]
hmax = max(sum(wi * ai for wi, ai in zip(w, a)) for a, b in ess)
bad = [(a, b) for a, b in ess if sum(wi * ai for wi, ai in zip(w, a)) == hmax
       and sum(wi * (bi - ai) for wi, ai, bi in zip(w, a, b)) > 0]
claim("G3a.6", len(bad) > 0, "w=(0,1): a w-maximal source among w-essential reactions points outward (D not endotactic)")


# polynomials as dicts {degree: Fraction}
def pmul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i + j] = out.get(i + j, 0) + a * b
    return {k: v for k, v in out.items() if v != 0}


def padd(p, q, c=1):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def pder(p):
    return {k - 1: k * v for k, v in p.items() if k > 0}


def pev(p, x):
    return sum(v * x ** k for k, v in p.items())


two_minus_v = {0: F(2), 1: F(-1)}
h = pmul(pmul(pmul(two_minus_v, two_minus_v), two_minus_v), {1: F(1)})
hp = pmul(pmul(two_minus_v, two_minus_v), {0: F(2), 1: F(-4)})
claim("G3b.1", pder(h) == hp, "h'(v) = (2-v)^2 (2-4v) for h(v) = (2-v)^3 v")
claim("G3b.2", pev(h, F(1, 8)) == F(3375, 4096) and pev(h, F(1)) == 1, "h(1/8) = 3375/4096, h(1) = 1")
claim("G3b.3", all(pev(hp, F(k, 800)) > 0 for k in range(0, 101)) and F(3375, 4096) - 1 == F(-721, 4096),
      "h' > 0 on [0,1/8] (sign of (2-v)^2(2-4v): 2-4v >= 3/2 there) and h(1/8) - 1 = -721/4096")
# vector fields with unit rate constants, symbolic check at random rational points
okf = True
for (u, v) in [(F(3, 2), F(1, 3)), (F(15, 8), F(1, 8)), (F(7, 5), F(2, 9))]:
    rates = [u ** 3 * v ** 2, v]
    vdot = sum(D.Gamma[1][j] * rates[j] for j in range(2))
    okf &= vdot == u ** 3 * v ** 2 - v
    a, b = u, v
    rM = [a ** 3, a ** 2 * b]
    adot = sum(M.Gamma[0][j] * rM[j] for j in range(2))
    okf &= adot == a ** 2 * (b - a)
claim("G3b.4", okf, "mass action with unit rates: v' = u^3 v^2 - v for D, a' = a^2 (b - a) for M")
# G3c: roots of h(v) = 1 on (0,2): h - 1 changes sign on (1/8, 1/2) and has root 1; h increasing on (0,1/2),
# decreasing on (1/2, 2), so exactly two roots
claim("G3c", pev(h, F(1, 8)) < 1 < pev(h, F(1, 2)) and pev(h, F(2)) == 0 and pev(h, F(1)) == 1,
      "u+v=2: h increases on (0,1/2) from 0 to 27/16 and decreases to 0 on (1/2,2): exactly two positive "
      "equilibria, v=1 and one v0 in (1/8,1/2)")

# ---------------- G5
Hn = CNet(2, [((2, 0), (0, 2)), ((0, 2), (2, 0)), ((2, 0), (1, 1)), ((1, 1), (2, 0))], names=["A", "B"])
Hd = Hn.dual()
ms = [to_set(S) for S in minimal_masks(siphon_masks(Hn))]
claim("G5.1", Hn.weakly_reversible() and H(Hn) and ms == [frozenset({0, 1})] and adls(Hn),
      "H is WR with (H1)-(H3); only minimal siphon {A,B}; ADLS holds")
claim("G5.2", Hd.describe() == "2k1 + 2k3 + k4 -> 2k2 + k3 + 2k4, 2k2 + k4 -> 2k1 + k3",
      f"dual of H: {Hd.describe()}")
claim("G5.3", H(Hd) and Hd.is_siphon([3]) and not adls(Hd) and Hd.Gamma[3] == [1, -1],
      "dual has (H1)-(H3), {U4} is a siphon with row (1,-1): no conservation law on it; ADLS fails")

# ---------------- G10
cx = sorted([c for c in [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]])
nets = []
pairs = list(combinations(cx, 2))
for k in (1, 2, 3):
    for ps in combinations(pairs, k):
        rx = []
        for a, b in ps:
            rx += [(a, b), (b, a)]
        nets.append(rx)
n_rev = len(nets)
for k in (3, 4):
    for sub in combinations(cx, k):
        for p in permutations(sub[1:]):
            cyc = (sub[0],) + p
            nets.append([(cyc[t], cyc[(t + 1) % k]) for t in range(k)])
claim("G10.1", n_rev == 575 and len(nets) == 705, f"575 reversible unions + 130 cycles = {len(nets)} networks")
h3 = [rx for rx in nets if CNet(3, rx).rows_distinct()]
both = only_n = only_d = neither = wrd = 0
qok = True
for rx in h3:
    N = CNet(3, rx)
    Dn = N.dual()
    a, b = adls(N), adls(Dn)
    both += a and b
    only_n += a and not b
    only_d += b and not a
    neither += (not a) and (not b)
    wrd += Dn.weakly_reversible()
    q = 3 - N.rank()
    qok &= q == (3 - len(Dn.complexes()) + Dn.linkage_classes()) + Dn.deficiency()
claim("G10.2", len(h3) == 702, f"{len(h3)} satisfy (H3)")
claim("G10.3", (both, only_n, only_d, neither) == (18, 177, 0, 507),
      f"ADLS both / only N / only N* / neither = {both} / {only_n} / {only_d} / {neither}")
claim("G10.4", wrd == 13, f"{wrd} have a WR dual")
claim("G10.5", qok, "q = beta* + delta* on all 702")
# ---------------- G4.1 examples and G8.1
import endo  # noqa: E402


def mod5(i):
    return (i - 1) % 5


rx5 = []
for i in range(1, 6):
    for c in (1, 2):
        a = [0] * 5
        b = [0] * 5
        a[mod5(i)] += c
        a[mod5(i + 2)] += c
        b[mod5(i + 1)] += c
        b[mod5(i + 3)] += c
        rx5.append((tuple(a), tuple(b)))
N5 = CNet(5, rx5)
D5 = N5.dual()
ms5 = [to_set(S) for S in minimal_masks(siphon_masks(N5))]
crit5 = [S for S in ms5 if not any(F <= S for F in semiflow_supports(N5))]
claim("G4.1a", N5.r == 10 and len(N5.complexes()) == 10 and N5.linkage_classes() == 2 and N5.rank() == 4
      and N5.deficiency() == 4 and N5.weakly_reversible() and H(N5),
      "5-species example: 10 reactions, 10 complexes, 2 linkage classes, s=4, deficiency 4, WR, (H1)-(H3)")
claim("G4.1b", D5.weakly_reversible() and D5.deficiency() == 0 and D5.linkage_classes() == 1
      and len(D5.complexes()) == 5 and D5.r == 5, f"its dual is one 5-cycle with deficiency 0: {D5.describe()}")
claim("G4.1c", frozenset({0, 1, 2}) in crit5, f"{{X1,X2,X3}} is a critical minimal siphon (critical: {crit5})")
claim("G4.1d", endo.strongly_endotactic(5, N5.reactions), "the 5-species example is strongly endotactic (exact)")
A4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
rx4 = [(A4[0], A4[1]), (A4[1], A4[2]), (A4[2], A4[3]), (A4[3], A4[0]), ((1, 0, 1, 0), (0, 1, 0, 1)),
       ((0, 1, 0, 1), (1, 0, 1, 0))]
N4 = CNet(4, rx4)
claim("G4.1e", N4.weakly_reversible() and N4.rows_distinct() and 4 - N4.rank() == 1 and N4.dual().weakly_reversible()
      and not endo.strongly_endotactic(4, rx4),
      "A->B->C->D->A, A+C<->B+D: WR, (H3), q=1, dual WR, NOT strongly endotactic (primality needed)")
w = (3, 2, 1, 2)
hts = [sum(x * y for x, y in zip(w, a)) for a, b in rx4]
mx = max(hts)
claim("G4.1f", all(sum(x * (bb - aa) for x, aa, bb in zip(w, a, b)) == 0 for (a, b), h in zip(rx4, hts) if h == mx)
      and mx == 4, "w=(3,2,1,2): the maximal sources (height 4) carry only w-orthogonal reactions")
# Example 8.1
R81 = CNet(3, [((3, 0, 2), (2, 3, 0)), ((2, 3, 0), (0, 2, 3)), ((0, 2, 3), (3, 0, 2))])
G81 = R81.Gamma
claim("G8.1a", G81 == [[-1, -2, 3], [3, -1, -2], [-2, 3, -1]] and R81.rank() == 2 and R81.deficiency() == 0
      and R81.dual().weakly_reversible() and R81.dual().deficiency() == 0, "Gamma, rank 2, deficiency 0, dual WR "
      "with deficiency 0")
Gd = [[10, 0, 2], [2, 1, 0], [0, 10, 1]]
J = [[sum(G81[i][k] * Gd[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
# restriction to S = {a+b+c=0}, basis f1 = e1-e3, f2 = e2-e3; coordinates of x in S: (x1, x2)
Jf = []
for f in ((1, 0, -1), (0, 1, -1)):
    img = [sum(J[i][j] * f[j] for j in range(3)) for i in range(3)]
    Jf.append((img[0], img[1]))
JS = [[Jf[0][0], Jf[1][0]], [Jf[0][1], Jf[1][1]]]
claim("G8.1b", JS == [[-15, 27], [24, -25]] and JS[0][0] * JS[1][1] - JS[0][1] * JS[1][0] == -273,
      f"J_S = {JS}, det = -273 (saddle)")
alt = CNet(3, [((10, 0, 2), (9, 3, 0)), ((2, 1, 0), (0, 0, 3)), ((0, 10, 1), (3, 8, 0))])
claim("G8.1c", alt.Gamma == G81, "10A+2C -> 9A+3B, 2A+B -> 3C, 10B+C -> 3A+8B has the same Gamma")
print(f"# runtime {time.time() - t0:.1f} s")
sys.exit(0 if OK else 1)
