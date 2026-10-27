"""V1 miscellaneous checks: precluster variants (Def. 19 vs thesis Alg. 4.4) on families and random networks,
loops and CF (Q2 section 3 remark), the agent-G example quoted in petri/notes.md section 9, and E's reaction-side
sentence on "MIPS refine preclusters" (computation/notes.md section 3(a)).

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Exact, standard library only.
"""
import sys, time, random, platform
from itertools import combinations
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
from check_P import family
from check_Q2 import random_wr

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
FAILS = []


def claim(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def gbar(N):
    """closure of G under adjunction through ANY edge, without the component restriction (manuscript G-bar)"""
    nv = nvec(N)
    E = set(doubling_pairs(N))
    changed = True
    while changed:
        changed = False
        for (a, b) in list(E):
            for s in range(N.r):
                if s in (a, b):
                    continue
                if in_span(nv[s], [nv[a], nv[b]]):
                    for e in ((min(a, s), max(a, s)), (min(b, s), max(b, s))):
                        if e not in E:
                            E.add(e); changed = True
    return components(N.r, list(E))


# 1. variants of the precluster definition
cnt, bad = 0, []
for gen in (family(2, 2, 4), family(3, 2, 3, zero=False), family(2, 3, 3, zero=False)):
    for N in gen:
        if not is_H2(N):
            continue
        cnt += 1
        a, b, c = preclusters_def19(N), preclusters_thesis(N), gbar(N)
        if not (a == b == c):
            bad.append(N)
claim("Def 19 (one step) = thesis Alg 4.4 = G-bar on all %d (H2) networks of the three small families" % cnt, not bad,
      repr(bad[:1]))
rng = random.Random(5)
cnt, bad = 0, []
for _ in range(3000):
    n = rng.randint(2, 4)
    cpx = [tuple(rng.randint(0, 2) for _ in range(n)) for _ in range(rng.randint(3, 5))]
    cpx = [c for c in set(cpx) if any(c)]
    if len(cpx) < 2:
        continue
    reacs = set()
    for _ in range(rng.randint(3, 7)):
        a, b = rng.sample(cpx, 2)
        reacs.add((a, b))
    reacs = sorted(reacs)
    Ye = [[a[i] for a, b in reacs] for i in range(n)]
    Yp = [[b[i] for a, b in reacs] for i in range(n)]
    N = Net(Ye, Yp)
    if not is_H2(N) or len(doubling_pairs(N)) == 0:
        continue
    cnt += 1
    if not (preclusters_def19(N) == preclusters_thesis(N) == gbar(N)):
        bad.append(N)
claim("Def 19 = thesis Alg 4.4 = G-bar on %d random (H2) networks with a doubling pair (r <= 7)" % cnt, not bad,
      repr(bad[:1]))

# 2. loops change P but not pi or kappa_V (Q2 section 3 remark)
R = net_from_str("X1 -> X2, X2 -> X1, X1 -> X1", species_order=["X1", "X2"])
R0 = net_from_str("X1 -> X2, X2 -> X1", species_order=["X1", "X2"])
V = interior_flux(R, [3, 5])
claim("loop X1 -> X1 added to R: P changes ({R1},{R2} -> one block); pi and kappa_V one block in both",
      preclusters_def19(R0) == [[0], [1]] and preclusters_def19(R) == [[0, 1, 2]] and
      len(flat_closure_partition(R)) == 1 and len(kappa_V(R, V)) == 1 and len(flat_closure_partition(R0)) == 1)

# 3. agent G's example quoted by P (section 9): 2A <-> 2B, 2A <-> A+B satisfies ADLS; dual has a critical siphon
G = net_from_str("2A <-> 2B, 2A <-> A+B")
def adls(N):
    sups = [frozenset(S) for S, _ in cons_rays(N)]
    return all(any(Z <= Sx for Z in sups) for Sx in siphons(N))
claim("G example 2A<->2B, 2A<->A+B: (H1); every siphon contains a semiflow support; dual is conservative and has a "
      "critical siphon; neither EC", is_H1(G) and adls(G) and is_H1(G.dual()) and not adls(G.dual()) and not is_EC(G),
      "dual: %r, critical minimal siphons %s" % (G.dual(), [sorted(x) for x in critical_siphons(G.dual())]))

# 4. E section 3(a): 'MIPS refine preclusters holds exactly when the MIPS are trivial' -- reaction side reading
W = net_from_str("X1 -> X2, X1 -> X3, X2 -> X1, X3 -> X1")
Ms, P = full_coupling(W), preclusters_def19(W)
print("   W (WR, deficiency 0, (H1)-(H3)): M(W*) = %s, P(W) = %s" % (Ms, P))
claim("E 3(a) sentence holds on the species side only: for W (WR, H) M(W*) is not discrete but M(W*) <= P(W)",
      refines(Ms, P) and Ms != [[j] for j in range(W.r)] and is_WR(W) and is_H1(W) and is_H3(W))
rng = random.Random(11)
cnt = cnt2 = 0
for _ in range(2000):
    N = random_wr(rng, rng.randint(2, 4), 2)
    if N is None or not is_WR(N) or not (is_H1(N) and is_H3(N)):
        continue
    cnt += 1
    Ms = full_coupling(N)
    if refines(Ms, preclusters_def19(N)) and Ms != [[j] for j in range(N.r)]:
        cnt2 += 1
print("   random WR networks with (H1),(H3): %d; with M(N*) non-discrete and M(N*) <= P(N): %d" % (cnt, cnt2))


# 5. Theorem B and the definition of twins: two species occurring in no reaction
N = Net([[1, 0], [0, 1], [0, 0], [0, 0]], [[0, 1], [1, 0], [0, 0], [0, 0]])  # X1 <-> X2; X3, X4 absent
D = N.dual()
PD = preclusters_def19(D)
ms = min_siphons(N)
lhs = all(any(set(b) <= S for S in ms) for b in PD)
nonzero_twins = [(i, j) for i, j in twin_pairs(N) if any(N.Ye[i])]
rhs_P = (not nonzero_twins) and all(any(i in S for S in ms) for i in range(N.n))
rhs_M = (not twin_pairs(N)) and all(any(i in S for S in ms) for i in range(N.n))
print("   X1<->X2 with absent X3, X4: P(N*) = %s, minimal siphons %s" % (PD, [sorted(x) for x in ms]))
claim("Theorem B needs twins = equal rows of Ye INCLUDING zero rows (manuscript definition): with P's 'nonzero rows' "
      "the criterion fails when two species occur in no reaction", lhs == rhs_M and lhs != rhs_P)
print("runtime %.1f s" % (time.time() - T0))
sys.exit(1 if FAILS else 0)
