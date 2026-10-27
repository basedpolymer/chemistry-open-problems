"""Exact certificates for the explicit networks of the manuscript.

Environment: Python 3.12 (WSL Fedora 44), standard library only; exact rational arithmetic (crn.py).
Runtime: < 2 s.

Every claim is a comparison between an object computed from scratch (extreme rays by exact linear algebra,
siphons by exhaustive search, preclusters by the three variants of the definition, MIPS by the three readings)
and the value stated in the manuscript.  Prints PASS/FAIL per claim; exits with status 1 if any claim fails.
"""
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from crn import Network, refines, fmt, mips_set_q2  # noqa: E402

FAILS = 0
NCLAIMS = 0


def claim(text, ok):
    global FAILS, NCLAIMS
    NCLAIMS += 1
    if not ok:
        FAILS += 1
    print(f"{'PASS' if ok else 'FAIL'}: {text}")


def fs(*blocks):
    """Partition / family from blocks given with 1-based indices."""
    return {frozenset(i - 1 for i in b) for b in blocks}


def H(N):
    return N.is_conservative() and N.is_consistent() and N.rows_distinct()


def check_all_variants(name, N, expected, what="P"):
    for v in ("onestep", "thesis", "closure"):
        claim(f"{name}: preclusters ({v}) = {sorted(sorted(i + 1 for i in b) for b in expected)}",
              set(N.preclusters(v)) == expected)


def e(*coeffs):
    return tuple(coeffs)


# ------------------------------------------------------------------------------------------------
print("== Consistency of the implementation with the examples of arXiv:2603.08767v2 and of the thesis")
# Example 6 / 15: A + E <-> AE -> Ap + E, species (A, E, AE, Ap) = (x1, x2, x3, x4)
N6 = Network(4, [(e(1, 1, 0, 0), e(0, 0, 1, 0)), (e(0, 0, 1, 0), e(1, 1, 0, 0)), (e(0, 0, 1, 0), e(0, 1, 0, 1))])
claim("Example 15: siphons are {1,2,3,4},{1,2,3},{1,3,4},{2,3,4},{1,3},{2,3}",
      set(N6.siphons()) == fs({1, 2, 3, 4}, {1, 2, 3}, {1, 3, 4}, {2, 3, 4}, {1, 3}, {2, 3}))
claim("Example 15: minimal siphons are {1,3},{2,3}", set(N6.minimal_siphons()) == fs({1, 3}, {2, 3}))
# Example 25: x3x4 -> x1x2, x1 -> x3
N25 = Network(4, [(e(0, 0, 1, 1), e(1, 1, 0, 0)), (e(1, 0, 0, 0), e(0, 0, 1, 0))])
D25 = N25.dual()
claim("Example 25: conservation laws x1+x3, x2+x4 are the extreme rays of the conservation cone",
      set(N25.conservation_rays()) == {(1, 0, 1, 0), (0, 1, 0, 1)})
claim("Example 25: minimal siphons {x4},{x1,x3}", set(N25.minimal_siphons()) == fs({4}, {1, 3}))
claim("Example 25: MIPS (ray reading) {x1,x3},{x2,x4}", set(N25.mips_ray()) == fs({1, 3}, {2, 4}))
claim("Example 25: MIPS (set reading) {x1,x3},{x2,x4}", set(mips_set_q2(N25)) == fs({1, 3}, {2, 4}))
claim("Example 25: dual network is 1 -> k1 (x2), k1 -> 1 (x4), k1 -> k2 (x3), k2 -> k1 (x1)",
      D25.reactions == [((0, 1), (1, 0)), ((0, 0), (1, 0)), ((1, 0), (0, 1)), ((1, 0), (0, 0))])
check_all_variants("Example 25: N*", D25, fs({1, 2, 3, 4}))
claim("Example 25: internal cycles of N* are {x1,x3},{x2,x4}", set(D25.internal_cycles()) == fs({1, 3}, {2, 4}))
claim("Example 25: N has no positive steady state (not consistent), so it is outside Conjecture 22",
      not N25.is_consistent() and N25.is_conservative())
# Thesis Example 1.22: 2X1 <-> X1 + X2 -> 2X2 (this is network N3 of thesis Example 4.1 with x = X1, y = X2)
N3 = Network(2, [(e(2, 0), e(1, 1)), (e(1, 1), e(2, 0)), (e(1, 1), e(0, 2))])
claim("thesis Example 1.22: internal cycles {1,2},{2,3}", set(N3.internal_cycles()) == fs({1, 2}, {2, 3}))
claim("thesis Example 4.1: the only doubling set of N3 is {2,3}", N3.doubling_pairs() == [(1, 2)])
check_all_variants("thesis Example 4.5: N3", N3, fs({1, 2, 3}))
# Thesis Example 4.12: x^2 <-> xy
N412 = Network(2, [(e(2, 0), e(1, 1)), (e(1, 1), e(2, 0))])
check_all_variants("thesis Example 4.12", N412, fs({1}, {2}))

# ------------------------------------------------------------------------------------------------
print("\n== Counterexample T: the 3-cycle X1 -> X2 -> X3 -> X1")
T = Network(3, [(e(1, 0, 0), e(0, 1, 0)), (e(0, 1, 0), e(0, 0, 1)), (e(0, 0, 1), e(1, 0, 0))])
TD = T.dual()
claim("T satisfies H (conservative, consistent, rows of Gamma distinct)", H(T))
claim("T: rows and columns of Gamma pairwise linearly independent",
      T.rows_pairwise_independent() and T.columns_pairwise_independent())
claim("T: weakly reversible, deficiency 0", T.weakly_reversible() and T.deficiency() == 0)
claim("T*: k1 -> k3, k2 -> k1, k3 -> k2 (a 3-cycle again) and T* satisfies H",
      TD.reactions == [((1, 0, 0), (0, 0, 1)), ((0, 1, 0), (1, 0, 0)), ((0, 0, 1), (0, 1, 0))] and H(TD))
for side, (net, base) in {"species side (P(T*) vs MIPS(T))": (TD, T), "reaction side (P(T) vs MIPS(T*))": (T, TD)}.items():
    check_all_variants(f"T, {side}", net, fs({1}, {2}, {3}))
    claim(f"T, {side}: MIPS ray reading = one block", set(base.mips_ray()) == fs({1, 2, 3}))
    claim(f"T, {side}: MIPS semiflow reading = one block", set(base.mips_semiflow()) == fs({1, 2, 3}))
    claim(f"T, {side}: MIPS set reading = one block", set(mips_set_q2(base)) == fs({1, 2, 3}))
    claim(f"T, {side}: minimal siphons = one set {{1,2,3}}", set(base.minimal_siphons()) == fs({1, 2, 3}))
    claim(f"T, {side}: internal cycles of the precluster network = {{1,2,3}}",
          set(net.internal_cycles()) == fs({1, 2, 3}))
    P = net.preclusters("thesis")
    for reading, M in (("ray", base.mips_ray()), ("semiflow", base.mips_semiflow()), ("set", mips_set_q2(base))):
        claim(f"T, {side}, {reading}: P != M, M does not refine P, |P| = 3 != 1 = |M|",
              set(P) != set(M) and not refines(M, P) and len(P) != len(M))
    S = base.minimal_siphons()
    claim(f"T, {side}: P != minimal siphons, a minimal siphon is not inside a precluster, no precluster is a "
          f"siphon, |P| != #minimal siphons",
          set(P) != set(S) and not refines(S, P) and not any(base.is_siphon(b) for b in P) and len(P) != len(S))
    claim(f"T, {side}: the internal cycle {{1,2,3}} is not contained in a precluster (Remark 23 fails)",
          not refines(net.internal_cycles(), P))

# ------------------------------------------------------------------------------------------------
print("\n== Counterexample R: X1 <-> X2 (smallest network satisfying H)")
R = Network(2, [(e(1, 0), e(0, 1)), (e(0, 1), e(1, 0))])
RD = R.dual()
claim("R satisfies H; R* = k1 <-> k2 satisfies H", H(R) and H(RD) and RD.reactions == [((1, 0), (0, 1)), ((0, 1), (1, 0))])
check_all_variants("R (reaction side)", R, fs({1}, {2}))
check_all_variants("R* (species side)", RD, fs({1}, {2}))
claim("R: MIPS(R) = {X1,X2} and MIPS(R*) = {k1,k2} (ray, semiflow and set readings)",
      set(R.mips_ray()) == set(RD.mips_ray()) == set(R.mips_semiflow()) == set(RD.mips_semiflow())
      == set(mips_set_q2(R)) == set(mips_set_q2(RD)) == fs({1, 2}))

# ------------------------------------------------------------------------------------------------
print("\n== Counterexample W: X1 -> X2, X1 -> X3, X2 -> X1, X3 -> X1 (preclusters do not refine MIPS)")
W = Network(3, [(e(1, 0, 0), e(0, 1, 0)), (e(1, 0, 0), e(0, 0, 1)), (e(0, 1, 0), e(1, 0, 0)), (e(0, 0, 1), e(1, 0, 0))])
WD = W.dual()
claim("W satisfies H, weakly reversible, deficiency 0, rows of Gamma pairwise linearly independent, columns distinct",
      H(W) and W.weakly_reversible() and W.deficiency() == 0 and W.rows_pairwise_independent() and W.columns_distinct())
claim("W: flux rays (1,0,1,0) and (0,1,0,1)", set(W.flux_rays()) == {(1, 0, 1, 0), (0, 1, 0, 1)})
check_all_variants("W (reaction side)", W, fs({1, 2, 3, 4}))
claim("W: MIPS(W*) (ray, semiflow, set) = {1,3},{2,4}",
      set(WD.mips_ray()) == set(WD.mips_semiflow()) == set(mips_set_q2(WD)) == fs({1, 3}, {2, 4}))
claim("W: P(W) does not refine MIPS(W*) (any reading)", not refines(W.preclusters("thesis"), WD.mips_ray())
      and not refines(W.preclusters("thesis"), WD.mips_semiflow()))
claim("W* = K1 + K2 -> K3 + K4, K3 -> K1, K4 -> K2 satisfies H",
      WD.reactions == [((1, 1, 0, 0), (0, 0, 1, 1)), ((0, 0, 1, 0), (1, 0, 0, 0)), ((0, 0, 0, 1), (0, 1, 0, 0))] and H(WD))
check_all_variants("W* (species side: preclusters of W** = W)", W, fs({1, 2, 3, 4}))
claim("W*: MIPS(W*) = {K1,K3},{K2,K4}, so P((W*)*) does not refine MIPS(W*)",
      set(WD.mips_ray()) == fs({1, 3}, {2, 4}) and not refines(W.preclusters("thesis"), WD.mips_ray()))

# ------------------------------------------------------------------------------------------------
print("\n== The thesis network N3 (2X1 <-> X1 + X2 -> 2X2) and thesis Example 4.12 (2X1 <-> X1 + X2)")
N3D = N3.dual()
claim("N3 satisfies H", H(N3))
claim("N3: MIPS(N3*) = {1},{2},{3} (no two reactions fully coupled)", set(N3D.mips_ray()) == fs({1}, {2}, {3}))
claim("N3: P(N3) = {1,2,3} does not refine MIPS(N3*)", not refines(N3.preclusters("thesis"), N3D.mips_ray()))
claim("N3: P(N3*) = {X1},{X2} while MIPS(N3) = {X1,X2}",
      set(N3D.preclusters("thesis")) == fs({1}, {2}) and set(N3.mips_ray()) == fs({1, 2}))
N412D = N412.dual()
claim("Example 4.12 satisfies H", H(N412))
claim("Example 4.12: internal cycle {1,2} is not contained in a precluster (Remark 23 fails)",
      set(N412.internal_cycles()) == fs({1, 2}) and not refines(N412.internal_cycles(), N412.preclusters("thesis")))
claim("Example 4.12: MIPS(N*) = {1,2} while P(N) = {1},{2}", set(N412D.mips_ray()) == fs({1, 2}))

# ------------------------------------------------------------------------------------------------
print("\n== Networks where an equality holds (Theorem 3(c))")
Q = Network(2, [(e(2, 0), e(1, 1)), (e(1, 1), e(2, 0)), (e(1, 0), e(0, 1)), (e(0, 1), e(1, 0))])
QD = Q.dual()
claim("Q = {2X1 <-> X1 + X2, X1 <-> X2} satisfies H and is weakly reversible (deficiency 1)",
      H(Q) and Q.weakly_reversible() and Q.deficiency() == 1)
claim("Q: P(Q) = MIPS(Q*) = discrete partition of the 4 reactions",
      set(Q.preclusters("thesis")) == set(QD.mips_ray()) == fs({1}, {2}, {3}, {4}))
claim("Q: P(Q*) = {X1},{X2} but MIPS(Q) = {X1,X2}: the species-side equality fails for Q",
      set(QD.preclusters("thesis")) == fs({1}, {2}) and set(Q.mips_ray()) == fs({1, 2}))

Dn = Network(3, [(e(1, 0, 0), e(0, 1, 0)), (e(1, 0, 0), e(0, 0, 1)), (e(0, 1, 1), e(2, 0, 0))])
DnD = Dn.dual()
claim("D = {X1 -> X2, X1 -> X3, X2 + X3 -> 2X1} satisfies H; D* satisfies H", H(Dn) and H(DnD))
claim("D: flux cone is the ray (1,1,1); doubling pair {1,2}",
      Dn.flux_rays() == [(1, 1, 1)] and Dn.doubling_pairs() == [(0, 1)])
claim("D: P(D) = MIPS(D*) = {1,2,3} (all variants, all readings)",
      all(set(Dn.preclusters(v)) == fs({1, 2, 3}) for v in ("onestep", "thesis", "closure"))
      and set(DnD.mips_ray()) == set(DnD.mips_semiflow()) == set(mips_set_q2(DnD)) == fs({1, 2, 3}))
claim("D: X2, X3 have identical reactant rows; P(D*) = MIPS(D) = {X1,X2,X3}",
      DnD.doubling_pairs() == [(1, 2)]
      and all(set(DnD.preclusters(v)) == fs({1, 2, 3}) for v in ("onestep", "thesis", "closure"))
      and set(Dn.mips_ray()) == set(Dn.mips_semiflow()) == set(mips_set_q2(Dn)) == fs({1, 2, 3}))
claim("D: minimal siphons of D are {X1,X2},{X1,X3}, so even here P(D*) differs from the minimal siphons",
      set(Dn.minimal_siphons()) == fs({1, 2}, {1, 3}))
claim("D: rows and columns of Gamma pairwise linearly independent; D* = K1 + K2 -> 2K3, K3 -> K1, K3 -> K2 "
      "(isomorphic to D)",
      Dn.rows_pairwise_independent() and Dn.columns_pairwise_independent()
      and DnD.reactions == [((1, 1, 0), (0, 0, 2)), ((0, 0, 1), (1, 0, 0)), ((0, 0, 1), (0, 1, 0))])
claim("D: Question 24, formulation Q3 fails on both sides: the single precluster lies in no minimal siphon",
      not refines(DnD.preclusters("thesis"), Dn.minimal_siphons())
      and not refines(Dn.preclusters("thesis"), DnD.minimal_siphons()))

# ------------------------------------------------------------------------------------------------
print("\n== Further entries of Table 1")
claim("T: |P(T)| = 3 != 1 = |MIPS(T)| and |P(R)| = 2 != 1 = |MIPS(R)| (formulation C, no dualization)",
      len(T.preclusters("thesis")) == 3 and len(T.mips_ray()) == 1
      and len(R.preclusters("thesis")) == 2 and len(R.mips_ray()) == 1)
claim("T: no precluster (of T or T*) is a union of minimal siphons (Q5)",
      not any(b in set(T.minimal_siphons()) or not b for b in T.dual().preclusters("thesis"))
      and all(not T.is_siphon(b) for b in T.dual().preclusters("thesis"))
      and all(not T.dual().is_siphon(b) for b in T.preclusters("thesis")))
claim("R: the internal cycle {R1,R2} lies in no precluster (Remark 23); P(R*) = singletons, MIPS(R) = one block",
      not refines(R.internal_cycles(), R.preclusters("thesis"))
      and set(RD.preclusters("thesis")) == fs({1}, {2}) and set(R.mips_ray()) == fs({1, 2}))
claim("N3: P(N3) = {1,2,3} refines no reading of MIPS(N3*) (ray, semiflow, set)",
      not refines(N3.preclusters("thesis"), N3D.mips_ray())
      and not refines(N3.preclusters("thesis"), N3D.mips_semiflow())
      and not refines(N3.preclusters("thesis"), mips_set_q2(N3D)))
claim("Ex. 4.12 is weakly reversible with deficiency zero; W, R, T too",
      all(net.weakly_reversible() and net.deficiency() == 0 for net in (N412, W, R, T)))

# ------------------------------------------------------------------------------------------------
print("\n== The one-site futile cycle F (species S0, S1, E, F, ES0, FS1)")
Fn = Network(6, [(e(1, 0, 1, 0, 0, 0), e(0, 0, 0, 0, 1, 0)), (e(0, 0, 0, 0, 1, 0), e(1, 0, 1, 0, 0, 0)),
                 (e(0, 0, 0, 0, 1, 0), e(0, 1, 1, 0, 0, 0)), (e(0, 1, 0, 1, 0, 0), e(0, 0, 0, 0, 0, 1)),
                 (e(0, 0, 0, 0, 0, 1), e(0, 1, 0, 1, 0, 0)), (e(0, 0, 0, 0, 0, 1), e(1, 0, 0, 1, 0, 0))],
             names=["S0", "S1", "E", "F", "ES0", "FS1"])
FD = Fn.dual()
claim("F satisfies H; F* satisfies H", H(Fn) and H(FD))
claim("F: extreme flux rays r1+r2, r4+r5, r1+r3+r4+r6",
      set(Fn.flux_rays()) == {(1, 1, 0, 0, 0, 0), (0, 0, 0, 1, 1, 0), (1, 0, 1, 1, 0, 1)})
claim("F: extreme conservation laws E+ES0, F+FS1, S0+S1+ES0+FS1",
      set(Fn.conservation_rays()) == {(0, 0, 1, 0, 1, 0), (0, 0, 0, 1, 0, 1), (1, 1, 0, 0, 1, 1)})
check_all_variants("F (reaction side)", Fn, fs({1, 2, 3, 4, 5, 6}))
check_all_variants("F* (species side)", FD, fs({1, 2, 3, 4, 5, 6}))
claim("F: MIPS(F*) = {r3,r6} and singletons; MIPS(F) = {S0,S1} and singletons",
      set(FD.mips_ray()) == fs({1}, {2}, {3, 6}, {4}, {5}) and set(Fn.mips_ray()) == fs({1, 2}, {3}, {4}, {5}, {6}))
claim("F: P(F) does not refine MIPS(F*) and P(F*) does not refine MIPS(F) (ray and semiflow readings)",
      not refines(Fn.preclusters("thesis"), FD.mips_ray()) and not refines(FD.preclusters("thesis"), Fn.mips_ray())
      and not refines(Fn.preclusters("thesis"), FD.mips_semiflow())
      and not refines(FD.preclusters("thesis"), Fn.mips_semiflow()))
claim("F: minimal siphons of F are {E,ES0},{F,FS1},{S0,S1,ES0,FS1}; of F* are {r1,r2},{r4,r5},{r1,r3,r4,r6}",
      set(Fn.minimal_siphons()) == fs({3, 5}, {4, 6}, {1, 2, 5, 6})
      and set(FD.minimal_siphons()) == fs({1, 2}, {4, 5}, {1, 3, 4, 6}))
claim("F: the single precluster is contained in no minimal siphon (both sides)",
      not refines(FD.preclusters("thesis"), Fn.minimal_siphons())
      and not refines(Fn.preclusters("thesis"), FD.minimal_siphons()))
claim("F: P(F*) = one block != MIPS(F) and P(F) = one block != MIPS(F*) (A1, B1)",
      set(FD.preclusters("thesis")) != set(Fn.mips_ray()) and set(Fn.preclusters("thesis")) != set(FD.mips_ray()))
claim("F: internal cycles of F* (semiflow supports) differ from MIPS(F) (ray reading)",
      set(FD.internal_cycles()) != set(Fn.mips_ray()))

# ------------------------------------------------------------------------------------------------
print("\n== X1 + 2X2 <-> X3: the readings of Definition 17 differ")
Y = Network(3, [(e(1, 2, 0), e(0, 0, 1)), (e(0, 0, 1), e(1, 2, 0))])
claim("Y satisfies H", H(Y))
claim("Y: extreme conservation laws X1+X3 and X2+2X3", set(Y.conservation_rays()) == {(1, 0, 1), (0, 1, 2)})
claim("Y: MIPS ray reading {X1},{X2},{X3}; set reading {X1,X2},{X3}; semiflow {X1,X3},{X2,X3}",
      set(Y.mips_ray()) == fs({1}, {2}, {3}) and set(mips_set_q2(Y)) == fs({1, 2}, {3})
      and set(Y.mips_semiflow()) == fs({1, 3}, {2, 3}))
claim("Y: P(Y*) = {X1},{X2},{X3}: equal to MIPS in the ray reading, not in the set or semiflow readings",
      set(Y.dual().preclusters("thesis")) == fs({1}, {2}, {3}))

claim("Y: the minimal siphons {X1,X3},{X2,X3} are not contained in MIPS (ray reading)",
      set(Y.minimal_siphons()) == fs({1, 3}, {2, 3}) and not refines(Y.minimal_siphons(), Y.mips_ray()))

# ------------------------------------------------------------------------------------------------
print("\n== Counterexample U: X4 -> X2 + X3, X3 + X4 -> X1, X1 + X2 -> 2X4 (twins X1, X2 not coupled)")
U = Network(4, [(e(0, 0, 0, 1), e(0, 1, 1, 0)), (e(0, 0, 1, 1), e(1, 0, 0, 0)), (e(1, 1, 0, 0), e(0, 0, 0, 2))])
UD = U.dual()
claim("U and U* satisfy H; rows and columns of Gamma(U) pairwise linearly independent",
      H(U) and H(UD) and U.rows_pairwise_independent() and U.columns_pairwise_independent())
claim("U: flux cone is the ray (1,1,1); extreme conservation laws (1,1,0,1) and (2,0,1,1)",
      U.flux_rays() == [(1, 1, 1)] and set(U.conservation_rays()) == {(1, 1, 0, 1), (2, 0, 1, 1)})
claim("U: X1 and X2 are twins (equal rows of Ye): the only doubling pair of U* is {X1,X2}",
      UD.doubling_pairs() == [(0, 1)])
check_all_variants("U* (species side of U)", UD, fs({1, 2, 3, 4}))
claim("U: MIPS(U) discrete (ray and set readings); semiflow supports {X1,X2,X4},{X1,X3,X4}",
      set(U.mips_ray()) == set(mips_set_q2(U)) == fs({1}, {2}, {3}, {4})
      and set(U.mips_semiflow()) == fs({1, 2, 4}, {1, 3, 4}))
claim("U: P(U*) refines no reading of MIPS(U)",
      not refines(UD.preclusters("thesis"), U.mips_ray()) and not refines(UD.preclusters("thesis"), U.mips_semiflow())
      and not refines(UD.preclusters("thesis"), mips_set_q2(U)))

# ------------------------------------------------------------------------------------------------
print("\n== Duality and the setting of the Global Attractor Conjecture (§5)")
for name, net in (("R", R), ("T", T), ("W", W), ("N3", N3), ("F", Fn), ("Y", Y), ("U", U), ("D", Dn)):
    claim(f"{name}: N** = N, and N* has r species, n reactions and the same rank s as N",
          net.dual().dual().reactions == net.reactions and net.dual().n == net.r and net.dual().r == net.n
          and net.dual().rank() == net.rank())
claim("W: weakly reversible with deficiency zero, but W* is not weakly reversible",
      W.weakly_reversible() and W.deficiency() == 0 and not WD.weakly_reversible())
claim("T: weakly reversible with deficiency zero, and T* = T is weakly reversible",
      T.weakly_reversible() and T.deficiency() == 0 and T.dual().weakly_reversible())

# ------------------------------------------------------------------------------------------------
print("\n== Explicit rate constants and positive steady states for every counterexample")
for name, net in (("R", R), ("T", T), ("W", W), ("W*", WD), ("N3", N3), ("Ex. 4.12", N412), ("F", Fn), ("F*", FD),
                  ("Y", Y), ("U", U), ("U*", UD), ("D", Dn)):
    k = [sum(v[j] for v in net.flux_rays()) for j in range(net.r)]   # positive since net is consistent
    x = [1] * net.n
    G = net.Gamma
    rhs = [sum(G[i][j] * k[j] for j in range(net.r)) for i in range(net.n)]  # x = 1 gives x^{Ye} = 1
    claim(f"{name}: k = {k}, x = {x} is a positive steady state (Gamma diag(k) x^Ye = 0 exactly)",
          all(c > 0 for c in k) and all(c == 0 for c in rhs))

print(f"\n{NCLAIMS - FAILS}/{NCLAIMS} claims PASS")
sys.exit(1 if FAILS else 0)
