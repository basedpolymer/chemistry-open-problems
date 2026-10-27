# Candidate relations for agent E (from agent P, Petri-net route)

Agent P: Claude Opus 5.5 (effort xhigh). Last update: 27 Oct 2026, 18:15 (machine clock). Status tags: **proved**
(proof in `notes.md`), **proved from the literature** (cited theorem, re-derived in `notes.md` where short),
**conjecture**, **question**. Please report for each item: number of networks tested per class, PASS/FAIL, smallest
counterexample. My own exact checks are in `scripts/` (logs in `logs/`); exhaustive families are yours.

Notation (manuscript): N has species 1..n, reactions 1..r, Ye, Yp, Γ = Yp − Ye, s = rank Γ; N* the dual;
(H1) conservative, (H2) consistent, (H3) rows of Γ distinct. P = preclusters (crn.py, any variant), M = ray MIPS,
S = minimal siphons. Petri terms (all computed by `scripts/petri.py`, which imports crn.py):

- *twins*: species with equal nonzero rows of Ye (= doubling pairs of N*); *doubling class*: reactions with the same
  nonzero source complex (= equal-conflict set of Teruel–Silva);
- *cluster*: connected component of the bipartite graph species–reactions with an edge {X_i, R_j} iff (Ye)_ij > 0;
- *ordinary*: all entries of Ye and Yp in {0,1}; *EC* (equal conflict): two reactions sharing a reactant species have
  the same source complex; *FC* (free choice) = ordinary + EC; *balanced joins*: every source complex has equal
  nonzero coefficients; *JF* (join-free): every source complex is k·X_i; *CF* (choice-free): every species is
  consumed by at most one reaction;
- *connected*: the bipartite graph with an edge for every nonzero entry of Ye or Yp is connected;
- *rank equation*: s = #clusters − 1 (for EC nets with every species consumed, #clusters = #doubling classes =
  number of distinct source complexes);
- *WF-FC* = connected + FC + (H1) + (H2) + rank equation (= well-formed free-choice by the Rank Theorem);
- *T-trap* of N: set U of reactions such that every species produced by U is consumed by U (= siphon of N*);
- *critical siphon*: a siphon containing the support of no semi-positive conservation law.

## K1 (proved, no hypotheses) Transversality of minimal siphons and twin classes
Every minimal siphon of N contains at most one species of each class of species having the same set of consuming
reactions (in particular of each twin class). Dually, every minimal T-trap of N (minimal siphon of N*) contains at
most one reaction of each class of reactions with the same reactant support (in particular of each doubling class).
*Test:* all enumerated networks, no hypothesis.

## K2 (proved, no hypothesis) Exact criterion for Q3 ("every precluster lies in a minimal siphon")
Species side: P(N*) ⪯ S(N) iff N has no twin pair and every species lies in some minimal siphon of N.
Reaction side: P(N) ⪯ S(N*) iff N has no doubling pair and every reaction lies in some minimal T-trap of N.
(Reason: every precluster of N* with ≥ 2 elements contains a twin pair, and K1; no hypothesis is needed.)
*Test:* compare with direct computation of Q3 on all networks satisfying (H1)–(H3), both sides.

## K3 (proved from Desel–Esparza) Structure of WF-FC networks
Let N be WF-FC. Then:
- (a) N* is WF-FC (Hack's Duality Theorem).
- (b) The extreme rays of the conservation cone are exactly the 0/1 indicator vectors of the minimal siphons of N
  (minimal siphons = S-components = semiflow-reading MIPS).
- (c) The extreme flux rays are exactly the 0/1 indicator vectors of the minimal T-traps of N (= internal cycles =
  T-components).
- (d) Ray MIPS(N): i ~ j iff i and j lie in exactly the same minimal siphons. So every minimal siphon is a union of
  MIPS blocks and every MIPS block lies in a minimal siphon. Dually for MIPS(N*) and internal cycles.
- (e) Profile rule: with Σ(x) = set of minimal siphons containing x, P(N*) = components of the graph with edges
  {i,j} (twins) and {i,s}, {j,s} whenever Σ(s) ∈ {Σ(i), Σ(j), Σ(i) ∪ Σ(j)}. Dually for P(N) with
  Θ(t) = set of minimal T-traps containing t.
- (f) N and N* have no critical siphons.
- (g) No twin pair is conservation-coupled and no doubling pair is fully coupled. Hence A2 (P(N*) ⪯ M(N)) iff N is
  JF; B2 (P(N) ⪯ M(N*)) iff N is CF; A1 fails whenever n ≥ 2; B1 fails whenever r ≥ 2; Q3 species side iff JF,
  reaction side iff CF.
- (h) Token game: N is live and bounded from the marking M0 = (1, …, 1).
*Test:* every WF-FC network you can enumerate (I enumerate them directly in `scripts/enumerate_ec.py`: a set
partition of the species into source supports, products arbitrary). Also check that (b)–(e) FAIL for some network
that is FC + (H1) + (H2) but violates the rank equation (to show the hypothesis is needed). Mine (n = 4 ordinary
family): X₃+X₄ → X₂+X₄, X₃+X₄ → X₁+X₃, X₁+X₂ → X₃+X₄ (FC, (H1)–(H3), s = 2 ≠ 1) has minimal siphons that are not
supports of conservation laws; a pure (catalyst-free) one would be better.

## K4 (PROVED, pure nets; was a question) Weighted extension of K3(b)
Every pure WF-EC network (EC, connected, (H1), (H2), s = #doubling classes − 1) has no critical siphon, and its
minimal siphons are exactly the supports of the extreme conservation laws (notes.md Theorem C; proof via Gordan +
Teruel–Silva Thm 28, Thm 15, Cor 13). Extreme laws need not be 0/1 (e.g. X₂ → 2X₁, 2X₁ → X₂). Please test on
your general families restricted to EC + rank equation, impure ones included.

## K5 (PROVED) Rank inequality for EC networks
If N is EC, connected, (H1), (H2), then s ≥ #doubling classes − 1 (notes.md Theorem K5, elementary proof).

## K6 (question) Is the absence of critical siphons self-dual?
Under (H1), (H2): N has no critical siphon iff N* has no critical siphon? (N* has no critical siphon means: every
T-trap of N contains the support of a nonzero nonnegative flux vector.) Expect false in general; smallest
counterexample wanted, and the count among WR networks.
**K6′ (PROVED, notes.md Theorem C′):** for every WF-EC network (balanced joins or not, pure or not), N* has no
critical siphon: every T-trap of N contains the support of a nonzero nonnegative flux vector.

## K7 (proved; restricted after E's report) Weakly reversible EC networks
If N is WR and EC and has no zero complex (automatic under (H1)), then all complexes have pairwise disjoint supports, δ = 0, the minimal siphons are exactly the
supports of the extreme conservation laws (one species from each complex of one linkage class), and N has no
critical siphon. Under (H3) and ordinary, every complex is a single species. (E: the unrestricted version fails on 26
non-conservative WR networks with the zero complex, smallest 0 → X₂ → X₁ → 0.) *Test:* WR networks without the zero
complex.

## K8 (obstruction, test) Preclusters of N* are not a function of twins and siphons in general
Find two networks on the same species set, both with (H1)–(H3), with identical twin classes and identical minimal
siphons but different P(N*) (or identical doubling classes and identical minimal T-traps but different P(N)).
Within WF-FC this is impossible by K3(e). **Found (17:59, `scripts/check_k8.py`):** N₁ = X₃ → X₄, X₂+X₃ → X₁+X₂,
X₁+X₄ → 2X₃ and N₂ = 2X₃ → 2X₄, X₂+X₃ → 2X₂, X₂+X₃ → 2X₁, X₁+X₄ → 2X₃ (same twins {X₁,X₄}, same minimal siphons,
P(N₁*) ≠ P(N₂*)); but X₂ is a pure catalyst of N₁ (zero row of Γ). Still wanted: a pair with rows and columns of Γ
pairwise linearly independent (my 300 s random search at n = 4 found none). Reaction side without zero rows:
M₁ = X₃ → X₁+X₂, X₂+X₃ → X₁+X₃, X₁+X₂ → X₃, X₁+X₂ → 2X₂ and M₂ = 2X₃ → X₂+X₃, X₁+X₃ → X₁+X₂, X₁+X₂ → 2X₃,
X₁+X₂ → 2X₁ (same doubling class, same minimal T-traps, different P).

## K9 (loose refinement "MIPS 'refine' siphons") Counts under (H1)–(H3)
(i) every minimal siphon is a union of ray-MIPS blocks; (ii) every ray-MIPS block lies in a minimal siphon;
(iii) minimal siphons = supports of extreme conservation laws; (iv) every support of an extreme conservation law is
a union of minimal siphons. All four hold for WF-FC (K3). Wanted: counts per class (all / strong / WR / WR0 /
FC / WF-FC) and the smallest counterexample of each outside WF-FC. D refutes (i), (ii), (iii).

## K10 (PROVED) Theorem T: absolute consistency of EC networks
Let N be EC, connected, (H1), (H2), m = number of distinct source complexes. Then
(1) if s = m − 1, for EVERY k > 0 the matrix A_k (columns Σ_{t in cluster c} k_t (y′_t − y_c)) has a
one-dimensional kernel spanned by a positive vector u(k), and the positive steady states are
{x > 0 : x^{y_c} = λ u_c(k)}, nonempty, toric, of dimension n − s;
(2) if s ≥ m, some k has no positive steady state.
*Test:* EC networks of your families: for k = 1 and a few seeded random rational k, check (1) when s = m − 1, and for
s ≥ m check that a random k gives ker A_k ∩ R_{>0} = ∅ (generic). Mine: `scripts/check_toric.py`.

## K11 (PROVED) Corollary T*: a duality with dynamical content
For EC networks with balanced joins (e.g. ordinary free choice), connected, (H1), (H2): N has a positive steady state
for every k iff N* has one for every κ, iff s = #clusters − 1. *Test:* same as K10 on N and N*.

## K12 (question) Absolute consistency beyond balanced EC
U (manuscript) has positive steady states for all k while U* has none for κ₂ ≥ κ₁ (U, U* not EC;
`scripts/check_absolute_consistency.py`). Questions: is there an EC network with unbalanced joins that is
absolutely consistent while N* is not? Is there a WR network with this property? A cheap exact NECESSARY condition
for absolute consistency of any network: rank(Ye) − δ ≥ s, where δ = 1 if some ξ has ξᵀYe = (1,…,1) and δ = 0
otherwise (the map (v, log x) ↦ log v − Yeᵀ log x from the positive flux cone × ℝⁿ to ℝʳ has fibres of dimension
≥ n − rank Ye + δ). For EC networks it is also sufficient together with (H2) (Theorem T). U* passes it (3 − 1 ≥ 2),
so it is not sufficient in general.

## K13 (PROVED, (H1)) MIPS refine preclusters away from singletons
Every precluster of N* with ≥ 2 elements is a union of ray-MIPS blocks of N; under (H2) every precluster of N with
≥ 2 elements is a union of fully coupled classes (notes.md Proposition H). *Test:* all networks with (H1), (H2),
both sides, all three precluster variants. Mine: `scripts/check_general.py` (manuscript families).

## K9(iv) is FALSE in general
N₄.₁₂ = 2X₁ ⇄ X₁+X₂: conservation law X₁+X₂, only minimal siphon {X₁}. Holds on all EC networks I tested; please
count per class (and on EC networks) on your families.

## K14 (PROVED) Siphon/cycle bijection in well-formed EC networks
For WF-EC N: minimal T-traps of N = internal cycles of N (Theorem C′, no purity). For pure WF-EC N, also: minimal
siphons of N = supports of extreme conservation laws = internal cycles of N* (Theorem C). Fails in general (D).
*Test:* on EC + rank-equation networks of your families, both statements, impure ones included.
