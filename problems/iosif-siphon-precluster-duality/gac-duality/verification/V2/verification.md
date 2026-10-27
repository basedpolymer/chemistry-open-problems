# Verification V2: well-formed equal-conflict networks (Q3), Theorem T and Corollary T* (P §7), Props E2.3–E2.5

Verifier V2: **Claude Opus 5.5 (effort xhigh)**, fresh context, 27 October 2026. Scope: `../../BRIEF-V2.md`.
Notes verified: `round2/ec/notes.md` (agent Q3, all of it), `petri/notes.md` §7 (agent P), `round2/computation/notes.md`
Props. E2.3–E2.5 and their use (agent E2). Every proof below was re-derived by hand; every example was re-checked
exactly by my own scripts (`scripts/`, own library `v2lib.py`, no import of the authors' code; the repository library
`scripts/crn.py` is imported only to compare preclusters). Logs in `logs/`. Nothing here quotes the author's e-mail.

Notation as in the notes: N has species X₁..Xₙ, reactions y_t → y′_t, Γ = Yp − Ye, s = rank Γ, S = im Γ,
L = S^⊥ (conservation laws); for an EC network the distinct sources y_1..y_m have disjoint supports S_c, T_c is the set
of reactions with source y_c, a_p = (y_{c(p)})_p, A_k = n × m matrix with columns Σ_{t∈T_c} k_t(y′_t − y_c),
φ_c(x) = x^{y_c}. WF-EC = EC + connected + (H1) + (H2) + s = m − 1. q = number of Petri clusters.

## 0. Verdicts

| claim | verdict |
|---|---|
| Q3 Theorem C (no purity) | **CONFIRMED** |
| Q3 Corollary CP (persistence of N and N*) | **CONFIRMED** |
| Q3 Theorem C^gen and Prop. K5gen | **CONFIRMED** (attribution: ordinary-net siphon part essentially in Desel–Esparza Thm 10.16–10.17, see §8) |
| Q3 Theorem U (a) uniqueness per class, (b) Hurwitz on S at every positive point, (c) sign condition | **CONFIRMED** (main statement); the "cluster kinetics" generalisation remark needs a correction (§2.4) |
| Q3 Theorem GS (global stability, max/min Lyapunov functions) | **CONFIRMED** (boundary, ties of the maximiser, analyticity step all checked, §2.5); remark (iii) needs the same correction as U |
| Q3 Theorem T*′ and its Corollary | **CONFIRMED** |
| Q3 Prop. T*-cex (N0, N1) | **CONFIRMED** (exact) |
| Q3 Prop. T*-dim (i)–(iii), NT, NTc, minimality of NT | **CONFIRMED** (exact for the examples; the minimality sketch re-derived) |
| Q3 Theorem U* | **CONFIRMED** |
| Q3 Prop. GS*≤2 | **CONFIRMED** |
| Q3 GS* (s ≥ 3, unbalanced) | open, as stated; my simulations agree (40 networks) |
| Q3 Prop. tor* (N2*) | **CONFIRMED** (exact) |
| Q3 Theorem PC (proportional conflict, self-dual) | **CONFIRMED** |
| Q3 Corollary BH | **CONFIRMED** |
| Q3 K8s (M3/M4, M1/M2) | **CONFIRMED** (exact; 8 matching relabellings for M3/M4, 4 for M1/M2, preclusters differ under all) |
| Q3 examples NM, NB, NT and the table of §7.2 | **CONFIRMED WITH CORRECTIONS** (labelling of the drained species in row 4; NM has two critical siphons) |
| Q3 Conjecture CR | **REFUTED** in its dynamical part (uniqueness and stability): exact counterexamples NC3 (n = 4, r = 4, molecularity 3) and NC, each with three positive steady states in one class, one a saddle (§2.12). Siphon part = Theorem C^gen (proved); absolute-consistency part still open |
| P Lemma CF, Theorem K5, Theorem T (1)(2) | **CONFIRMED**; credits to PMDSC 2012 and Campos–Chiola–Silva 1991 checked against the sources: accurate |
| P Corollary T* | **CONFIRMED** |
| P Theorem C′ (in §7) | **CONFIRMED** (limit argument checked; also subsumed by C^gen) |
| E2 Prop. E2.3 and its use on WF-EC families | **CONFIRMED** (it is PMDSC 2012 Thm 5.5 specialised by Theorem T, as E2 says) |
| E2 Prop. E2.4 (log-cone criterion) | **CONFIRMED** |
| E2 Prop. E2.5 and Cor. E2.5′ | **CONFIRMED** |

Overall: the theory of well-formed equal-conflict networks and their duals holds up. I found no error in any proof
marked "proved"; the substantive new finding is that **Conjecture CR is false** beyond proportional conflict, which
makes Theorem PC the right boundary for the dynamical statements (the siphon statement C^gen goes further).

## 1. Standing lemmas (Q3 §0, P §5, §7)

All re-derived; correct as written.

- **Under (H2) every species occurring is consumed** (if row i of Ye is 0 then row i of Γ is ≥ 0 and annihilates
  v > 0, so it is 0). Hence for EC networks clusters ↔ sources, m = q.
- **Lemma SC** (connected + (H1) + (H2) ⇒ Petri net strongly connected): the terminal-component computation with z > 0,
  v > 0 is right; both sums are ≥ 0 and vanish, so no arc enters the terminal component.
- **(Z) and (Σ)**: for p ∈ S_c, d ≠ c, A_k[p, d] = Σ_{t∈T_d} k_t (y′_t)_p ≥ 0; A_k1 = Γk. For general networks the same
  holds for Ā_k (columns Γ(k ⊙ 1_{T_c})) because a species of S_c is a reactant of no reaction outside T_c.
- **Lemma CF / CF′**: the support of a nonnegative kernel vector is closed under successors in the cluster graph
  (row p ∈ S_c has a positive coefficient Σ_{t∈T_c} k_t(y_t)_p on u_c), strong connectivity gives u > 0, and the
  "largest ε" argument gives a one-dimensional kernel. Correct for impure networks (own-cluster terms only enter the
  left side with the consumption coefficient; production terms of the own cluster go to the right side).
- **Lemma M** (Perron–Frobenius for irreducible Metzler G with a positive right or left null vector): correct
  (G + cI is primitive because its diagonal is positive; a positive eigenvector forces the Perron root). The
  reducible variant (B1 = 0 ⇒ a nonzero ξ ≥ 0 with ξᵀB = 0) is also correct: B + cI has all row sums c, so ρ = c.
- **Lemma IM** (WF-EC: im A_k = S for every k): rank A_k = m − 1 = s and im A_k ⊆ S. Correct.
- **Lemma R** (cluster version for connected + (H1) + (H2) + s = q − 1): closed–open argument with CF′; correct.

## 2. Agent Q3 (`round2/ec/notes.md`)

### 2.1 Theorem C without purity — CONFIRMED

Re-derived. Transversality is P's Theorem A, which is stated for species with the same consumer set (not only twins),
so it applies to the species of one cluster even when their coefficients a_p differ (unbalanced joins): the citation
is correct. Steps 1–3 (Metzler, vanishing outside C_R by the siphon property, zero row sums by (Σ)), the reducible
Perron–Frobenius step and Lemma IM give zᵀΓ = 0 with z ≥ 0, z ≠ 0 supported in R. The consequence (minimal siphons =
supports of extreme laws) is the standard argument. No purity, no balanced joins used. Remark (iii) (splitting a
self-loop reaction through a new species Q_t) is also correct (rank and cluster count both increase by one per split).

### 2.2 Corollary CP — CONFIRMED

ADLS 2007, Theorem 2 (journal version on disk, p. 605–606) states: if the associated Petri net is conservative and
every siphon contains the support of a P-semiflow, the network is persistent, for the kinetics of their §2 (mass
action included). N is conservative with no critical siphon (C); N* is conservative (N consistent) with no critical siphon (P's C′, or
C^gen applied to N*).

### 2.3 Theorem C^gen and Proposition K5gen — CONFIRMED

Re-derived. The only place where the argument differs from C is the use of Gordan's alternative on the non-square
B = Ā_v[R, C_R]: if Bw > 0, at a maximiser c* of w and p ∈ R ∩ S_{c*}, (Bw)_p = Σ_{d≠c*} B[p,d](w_d − w_{c*}) ≤ 0,
contradiction. Then zᵀĀ_v = 0 and im Ā_v = S (rank q − 1 = s by CF′ with u = 1). N* satisfies the same hypotheses
(same clusters with the roles exchanged, same rank, connected, (H1) ↔ (H2)). K5gen: rank Ā_v = q − 1 ≤ s.
Independent checks (exact, `logs/random_wfec_v2_seed{3,4,5,6}.log`): no critical siphon on either side for 4 × 120
WF-EC, 4 × 120 PC and 4 × 120 general (sharing, non-proportional sources) cluster-rank networks generated with
n = 3..6, and for the 437 cluster-rank networks (14 not EC) among 4 × 240 unconstrained random networks with
connected, (H1), (H2) (n ≤ 4, molecularity ≤ 3, r ≤ 5); K5gen held on all 960 unconstrained ones.

### 2.4 Theorem U — CONFIRMED (main statement); correction to the "cluster kinetics" remark

(a) Uniqueness: the argmax/argmin argument on x − x′ = A_v w is correct (row sums of A_v vanish; at an argmax cluster
x ≤ x′ on S_{c*}, so δ = Σ a_p d_p ≤ 0; δ ≥ 0 at an argmin; then x = x′ on argmax clusters and the argmax set is closed
under predecessors). Existence: Brouwer on time-τ maps plus "the zero set of a boundary equilibrium is a siphon",
which contains the support of a conservation law (C): correct.
(b) **Hurwitz at every positive point**: J(x) = A_k Dφ(x); G(x) = Dφ(x)A_k is Metzler and irreducible at every x > 0
(G[c,d] = Σ_{p∈S_c} a_p φ_c(x)/x_p · A_k[p,d] > 0 iff d → c) and G(x)u(k) = 0 with u(k) independent of x. Lemma M
gives a simple 0 and m − 1 eigenvalues with negative real part; det(tIₙ − A_kDφ) = t^{n−m} det(tI_m − DφA_k), so
det(tI − J) = t^{n−m+1}q(t) = t^{n−s}q(t); J maps ℝⁿ into S and induces 0 on ℝⁿ/S (dimension n − s), so the
characteristic polynomial of J|_S is exactly q. Correct, at every positive x, steady or not.
(c) The same argument with (s, d) of equal sign vectors: correct.
Independent exact checks (exact characteristic polynomial by Faddeev–LeVerrier over ℚ, Hurwitz determinants):
U(b) at random rational non-steady points, U(c) by maximal minors, on my exhaustive re-enumeration of P's families
(2,3,≤5) (71 WF-EC), (3,2,≤5) (958), (3,3,≤4) (6 995, of which 52 unbalanced) and (4,2,≤4) (692), 8 716 in all,
and on 480 random
WF-EC networks (n ≤ 6; in seed 3, 57 of 120 unbalanced,
100 of 120 impure): all PASS. Control: on EC networks with s ≥ m, at a positive steady state the test fails as it must
(0 has multiplicity n − m + 1 > n − s): 4 × 30 controls, all detected.

**Correction (remark "Generalisation (cluster kinetics)", §2, and Remark (iii) of §3).** For (b) (and GS) strict
monotonicity of φ_c is not enough: the irreducibility of G needs ∂φ_c/∂x_p > 0 for p ∈ S_c at every positive point.
Example: one-species cluster with φ(x) = (x − 1)³ + 1 (C¹, ≥ 0 on [0,∞), vanishing only at 0, strictly increasing);
at x = 1 its derivative is 0, the corresponding row of G vanishes, G is reducible and J|_S has an extra zero eigenvalue.
Corrected statement: "(a) holds for C¹ cluster rates φ_c that depend on S_c only, are strictly increasing in each
variable and vanish iff a variable vanishes; (b), and Theorem GS when φ_c is moreover analytic, hold if in addition
∂φ_c/∂x_p > 0 on the open orthant for every p ∈ S_c."

### 2.5 Theorem GS — CONFIRMED

The brief asked for particular care here; each step was re-derived.
1. *Identity.* ψ̇ = Dφ(x)ẋ = G(x)ψ with G(x)u = 0 and u = u(k) fixed. Equivalently, with h_c = ψ_c/u_c,
   ḣ_c = Σ_{d≠c} G[c,d](x) (u_d/u_c)(h_d − h_c): a nonlinear consensus system with state-dependent nonnegative
   weights (this is why max and min are monotone; see §8 for credit).
2. *Non-unique maximiser.* M = max_c h_c is a maximum of finitely many C¹ functions, so D⁺M(t) = max over the active
   indices of ḣ_c(t) ≤ 0; a continuous function with D⁺M ≤ 0 everywhere is nonincreasing. Ties are therefore
   harmless; no smoothness of M is needed. Same for μ.
3. *Boundary.* ψ_c(t) ≥ μ(0)u_c > 0 for all t ≥ 0 and all c, and x is bounded above by (H1); since every species
   lies in some S_c, x_p(t)^{a_p} ≥ μ(0)u_c / Π_{q∈S_c, q≠p} B^{a_q} > 0. So no trajectory can approach the boundary:
   the forward orbit lies in a compact K of the open orthant (and exists for all t ≥ 0). This is exactly the point
   where a "two monotone functions" argument could fail, and it does not.
4. *LaSalle step.* W = log M − log μ is continuous on the open orthant, nonincreasing along trajectories, W = 0 iff
   φ(x) ∈ ℝ_{>0}u(k) iff x is a steady state (ker A_k = ℝu). On the ω-limit set Ω ⊆ K, W ≡ W*. Along a full orbit in Ω,
   log M is nonincreasing, log μ nondecreasing and their difference constant, so both are constant.
5. *Analyticity step.* h_c is real-analytic in t (polynomial field), so each set {t : h_c(t) = M̄} is discrete or all
   of ℝ; the sets cover ℝ, so some h_c ≡ M̄; for such c, 0 = ḣ_c = Σ_{d≠c} G[c,d](u_d/u_c)(h_d − M̄) with all terms ≤ 0
   and G[c,d] > 0 for every predecessor d of c at every interior point, so the set I of such c is closed under
   predecessors, hence everything: W = 0, contradiction. So W* = 0, Ω consists of steady states of the class, and Ω is
   the unique steady state of Theorem U(a).
No gap. Numerical red team (`logs/gs_sim_v2.log`, floating point, not a proof): 40 random WF-EC networks with
n = 4..6, starts at random points, in the same class, and at 10⁻⁶ from a face: largest relative increase of M(t) on a
dense grid 1.3·10⁻¹², largest relative decrease of μ(t) 9.4·10⁻¹³ (rounding level), end points of one class coincide to
9·10⁻¹³, residual ≤ 7·10⁻¹⁴; along trajectories started at relative distance 10⁻⁶ from a face, min_t min_p x_p(t)
stays at that scale (1·10⁻⁷ relative) and they converge like the others. The dual N* (GS* is open for s = 3)
behaved the same way on these 40 networks: end points of one class within 5·10⁻¹³.
Remark (ii) (GAC for weakly reversible EC networks without the zero complex): correct — complexes are pairwise
disjoint, each linkage class is a WF-EC network on its own species (z_p = 1/(a_p|S_c|) is a positive law), and the
class of the whole network is the product of the classes of the linkage classes.

### 2.6 Theorem T*′ — CONFIRMED

Re-derived. ω ∈ L iff ω ⊥ im A_v = S (Lemma IM), i.e. H_c(π) = ω(π)·a_c(v) = 0 for all c; Σ_c H_c ≡ 0; H is
cooperative by (Z); H_c ≥ 0 on {π_c = 0}; H_c ≥ −Cπ_c on the simplex; π ↦ π + εH(π) maps the simplex into itself;
a Brouwer fixed point is a zero of H; its zero set Z is closed under successors in the cluster graph, hence empty.
Then k with Π_{t∈T_c} k_t = π_c exists because the T_c are disjoint. The Corollary (N absolutely consistent ⇒ s = m − 1
by P's T(2) ⇒ N* absolutely consistent) is correct.
Independent check (`logs/tstar_v2.log`): 300 random WF-EC networks (99 unbalanced): for m = 2 (112 networks) the
root of H₁(π₁, 1 − π₁) in (0, 1) is an algebraic number and Γᵀω(k) reduces to 0 modulo its minimal polynomial
(exact); for m = 3, 4 (188 networks) the cooperative flow π̇ = H(π) converges to a positive zero with residual
< 10⁻⁸ (numerical). All PASS.

### 2.7 Prop. T*-cex (N0, N1) and Prop. T*-dim — CONFIRMED

Exact (`logs/check_examples_v2.log`): N0 = {2A+B → 3A, 2A+B → A+2B} is EC, unbalanced, s = m = 1, (H1)–(H3) and
connected on both sides, N0* = {2K₁+2K₂ → 3K₁+K₂, K₁+K₂ → 2K₂}; ẋ_A = (k₁ − k₂)x_A²x_B; (1, κ_B/κ_A) is a positive steady
state of N0* for all κ. N1: ker Γ = ℝ(3, 2, 1), so positive steady states need k₂ = 2k₃; N1* = {2K₁ → 3K₂, K₁ → 3K₃,
K₂+K₃ → K₁} has the stated steady state for all κ. T*-dim (i), (ii): the dimension count of
{κ_p = w_p π_{c(p)}^{−a_p}} is right ((n − s) + m < n if s ≥ m + 1; with balanced joins the free ℝ_{>0}-action
(w, π) ↦ (tw, t^{1/a}π) lowers it to n − 1 when s = m; Sard). (iii): NT and NT* are connected with (H1)–(H3), every law
of NT has z₁ = z₂, ρ₁ and ρ₂ of NT* share the source K₁+K₂; for κ = (1, 2, 1, 1) the steady-state ideal of NT*
contains K₁K₂ (Gröbner basis), and for NTc* with κ = (1, 1, 2) it contains K₁K₂. Minimality of NT among molecularity
≤ 3: I re-derived the case analysis (n = 3, m = 2, unbalanced: the 3 × 3 log-linear system has determinant
a₃(a₁ − a₂) ≠ 0, always solvable; n = 3, m = 1: an unbalanced one-cluster source with the needed coincidence of
exponents has molecularity ≥ 4; n = 2 always solvable; r = 2 forces m = 1). Correct.

### 2.8 Theorem U* — CONFIRMED

(a) J* = Γᵀdiag(aω)E_S · E_Tᵀdiag(1/k); G* = Q*P* has G*[c,d] = Σ_{t∈T_c} k_t⁻¹ Σ_{p∈S_d} (y′_t)_p a_p ω_p ≥ 0 for
d ≠ c, irreducible; ξᵀG* = (A_{1/k}ξ)ᵀdiag(aω)E_S, so ξ = u(1/k) > 0 (Theorem T for N with rates 1/k) is a positive left
null vector. Lemma M and the multiplicity count r − m + 1 = r − s give J*|_{S*} Hurwitz at every positive k.
(b) Degree argument: classes compact (N consistent), no boundary steady state (C′), the homotopy λF + (1 − λ)(k̂ − k)
has no zero on the relative boundary (F is in the tangent cone, k̂ − k points strictly inward), every zero has index
(−1)^s, so there is exactly one. Correct.
Independent exact check: J*|_{S*} Hurwitz at random positive points with random κ, on all networks of §2.4 (8 716
enumerated, including the 52 unbalanced ones of (3,3,≤4), + 480 random WF-EC): all PASS.

### 2.9 Prop. GS*≤2 and Prop. tor* — CONFIRMED

s = 1: scalar field on a segment with one hyperbolic interior zero and inward-pointing ends. s = 2: trace(J*|_{S*}) < 0
at every interior point (U*(a)), Bendixson–Dulac excludes periodic orbits and cycles of orbits in the simply connected
interior, ω-limit sets avoid the boundary (CP), Poincaré–Bendixson with a unique equilibrium. Correct. tor*: N2 is
WF-EC and unbalanced, N2* = {2K₁ → 2K₂, K₁ → K₂, K₂ → K₁}, steady states K₂ = (2κ₁K₁² + κ₂K₁)/κ₃, and
d²(log K₂)/d(log K₁)² = 2κ₁κ₂K₁/(2κ₁K₁ + κ₂)² ≠ 0: not toric (exact).

### 2.10 Theorem PC and Corollary BH — CONFIRMED

PC definition and its rank-one block description: correct (sharing is transitive inside a cluster, positive
proportionality is transitive; all sources of a cluster have the same support S_c; integrality of a_t follows from
primitivity of b_c). Self-duality: Ye(N*) = Yeᵀ has the transposed rank-one blocks on the same clusters. (b) The
use of Lemma R on N* with routing z (im Γᵀ spanned by Γᵀ(z ⊙ 1_{S_c})), cooperativity of H (for d ≠ c only
productions (y′_t)_p ≥ 0 appear), Brouwer and "zero set closed under predecessors": correct. (c) J = PQ, G = QP Metzler
and irreducible with positive null vector u(a ∘ w) (Lemma R): correct. (d) Degree argument: correct.
Independent exact checks: on 4 × 120 random PC networks (41 of them neither EC nor BH): N* PC and well-formed, J|_S and
J*|_{S*} Hurwitz at random points, no critical siphon on either side: all PASS. NP and NP* checked exactly.
BH: K* is EC with clusters = blocks (rows q, q′ of Ye(K) coincide when q, q′ lie in one block); then T*′, U*, C′ and
GS*≤2 applied to N = K*. K is EC iff the multiples are constant on each block. Correct.

### 2.11 K8s — CONFIRMED (`logs/check_k8s_v2.log`)

M3 and M4: (H1)–(H3) on N and N*, weakly reversible, pure, rows and columns of Γ pairwise linearly independent on
both networks and on their duals; doubling classes, minimal T-traps and preclusters as stated (my own precluster code
agrees with `crn.py` in its three variants); **8** relabellings match the doubling classes and the minimal T-traps and
the preclusters differ under all 8. M1/M2 (Q3's "smallest", also E2's): impure, WR, (H1)–(H3) both sides, strong
non-degeneracy on both sides; 4 matching relabellings, preclusters differ under all.

### 2.12 The examples of §7.2 and Conjecture CR

**Table of §7.2 — CONFIRMED WITH CORRECTIONS** (`logs/check_examples_v2.log`). Every row is connected with
(H1)–(H3) on N and N*, and the stated (m_src, m_cl, s) are right. NM: two positive steady states in
x₁ + x₂ + x₃ = 5/2 for k = (1, 2, 1) (the two roots of x₂³ − 2x₂² + 1/4 in (0, 2), exact); NM is in fact PC (2X₃ and 3X₃
are proportional), so it is a PC network one above the rank equation. It is "mixed" (one cluster of BH type, one of
EC type) and hence neither EC nor BH at molecularity 3; Q3's remark §4b (iv) ("genuine PC networks need molecularity
≥ 6") is right only if "genuine" means a single cluster with two proportional sources of non-0/1 shape. NB: two positive steady states in
2x₁ + x₂ + x₃ = 6, x₁ = (5 ± √17)/4 (exact). U*: f_{K₁} − f_{K₂} = (κ₂ − κ₁)K₃ + 2κ₃K₂ (exact).
Corrections: (1) NM has **two** critical minimal siphons, {X₁,X₃} and {X₂,X₃} (the text names one; the later sentence
"each with two critical siphons, as N" is right). (2) Row 4: with the reactions numbered in the path order written in
the table (X₄ → 3X₃, 3X₃ → X₁+2X₂, X₁+2X₂ → 3X₂, 3X₂ → X₄), the drained species of N* is **k₃** (the dual species of
X₁+2X₂ → 3X₂, with k̇₃ = k₃(κ₂k₃k₄³ − κ₁) and {k₃} the only critical minimal siphon); "k₄" is agent E's label, which
numbers the reactions in another order. Write "N* loses the species dual to X₁+2X₂ → 3X₂".

**Conjecture CR — REFUTED (dynamical part).** CR states: connected, (H1), (H2), s = m_cl − 1 ⇒ N (and N*) absolutely
consistent with exactly one, exponentially stable, positive steady state per class, for every rate vector. My random
family of cluster-rank networks with sharing, non-proportional sources produced a network whose Jacobian is not
Hurwitz on S at some point (`logs/random_wfec_v2_seed3.log`, `logs/cr_probe.log`); deleting a pure catalyst and one
reaction gives the exact counterexample (`scripts/cr_counterexample.py`, `logs/cr_counterexample.log`):

NC = {C → B, 2A → D, 2B+D → 2A+2B, 2B+2D → 2A+2C+D} (species A, B, C, D).

- connected, (H1) (laws A/2 + D and B + C), (H2) (flux rays (0,1,1,0), (2,1,0,1)), (H3); s = 2; Petri clusters
  {C}|{C → B}, {A}|{2A → D}, {B, D}|{2B+D → …, 2B+2D → …}: q = 3, so s = q − 1; TEFC (2B+D and 2B+2D have the
  same support) but not PC (they are not proportional), not EC; NC and NC* have no critical siphon (as C^gen says).
- k = (2/3, 5/18, 24, 64): x = (3, 1/2, 3, 1/4) is a positive steady state (exact) with det(J(x)|_S) = −88/9 < 0, i.e.
  index −1, whereas (−1)^s = +1.
- Degree argument (the one of U*(b), valid here because the class is a compact polytope without boundary steady
  states by C^gen): the class B + C = 7/2, A + 2D = 7/2 contains further positive steady states of total index +2.
  An exact resultant (real roots isolated exactly, then evaluated to 25 digits) shows exactly three:
  (A, B, C, D) ≈ (3.4728, 3.1485, 0.35145, 0.013589) and (2.7616, 0.34709, 3.1529, 0.36921), both stable on S, and
  the saddle (3, 1/2, 3, 1/4) (eigenvalue ≈ +0.371 on S).

So NC is multistationary (bistable: two stable steady states and a saddle in one class) and has an unstable positive
steady state although it satisfies every hypothesis of CR; it is even TEFC, so "TEFC + cluster rank equation" does not
suffice either and proportional conflict is exactly what Theorem PC needs;
it is also absolutely consistent (k₁C = 2k₄B²D², k₂A² = B²D(k₃ + k₄D) are solved by choosing B, D and then C, A), so
the absolute-consistency part of CR is not refuted by it. The speculation of §7.3 ("such N is always, up to pure
catalysts, EC or block-power") is also false (NC has no pure catalyst and is not PC). The scan of §7.3 (all non-EC
cluster-rank networks of E's families are EC or BH up to catalysts) is not contradicted: those families have
molecularity ≤ 3 with n ≤ 3, or n = 4 with molecularity ≤ 2, and NC has a source of molecularity 4.

**A smaller counterexample, molecularity 3** (found by the targeted search `scripts/cr_search2.py`, certified by
`scripts/cr_counterexample2.py`, `logs/cr_counterexample2.log`):

NC3 = {D → B+C, C → A, A+B → B+C, 2A+B → C+D} (species A, B, C, D; n = 4, r = 4).

- connected, (H1)–(H3) on NC3 and on NC3* (laws B + D and A + C + D; flux rays (0,1,1,0), (1,2,0,1)); s = 2; Petri
  clusters {D}, {C}, {A, B} (the last with sources A+B and 2A+B): q = 3 = s + 1; TEFC, not PC, not EC; no critical
  siphon on either side.
- k = (1/6, 1/3, 6, 6): x = (1/2, 1/3, 6, 3) is a positive steady state (exact) with det(J(x)|_S) = −7/9 < 0.
- In its class B + D = 10/3, A + C + D = 19/2, msolve (exact real root isolation of the steady-state equations with
  the two conservation equations) finds exactly three real solutions, all positive: (0.25380, 1.00437, 6.91724,
  2.32897) and (2.07954, 0.021275, 4.10840, 3.31206), both stable on S, and the saddle (1/2, 1/3, 6, 3) (eigenvalue
  ≈ +0.096 on S).

NC3 lies outside E's general families (n = 4 there only with molecularity ≤ 2), which explains why the scans of
Q3 and E2 did not meet it. The targeted search found several other witnesses with n = 4, r = 4, molecularity 3
(`logs/cr_search2_*.log`): about 4 % of the non-PC cluster-rank networks sampled have such a steady state; none
with n = 3 was met (§6).
Corrected status for the paper: **siphon part: Theorem C^gen (proved); uniqueness/stability: false outside PC
(NC); absolute consistency for every connected (H1)–(H2) network with s = q − 1: open.** Theorem PC is the right
dynamical boundary.

## 3. Agent P, `petri/notes.md` §7

- **Lemma CF, Theorem K5**: confirmed (§1).
- **Theorem T** (1): the closed–open argument (𝒢 closed by normalisation; open because at k ∈ 𝒢 the rank is m − 1,
  stays m − 1 nearby, and the kernel line moves continuously through a strictly positive vector) and the toric
  description (the y_c are linearly independent, so the m − 1 log-linear equations are independent; dimension
  n − m + 1 = n − s): correct. (2): k_t = v_t/u_{c(t)} with v ∈ ker Γ ∩ ℝʳ_{>0}, u ∈ ℝᵐ_{>0} (φ is onto ℝᵐ_{>0}
  because the sources have disjoint supports), image of dimension ≤ (r − s) + m − 1 < r when s ≥ m: correct.
  No purity, no (H3), no balanced joins needed, as stated. **CONFIRMED.**
- **Corollary T***: N* is EC with balanced joins, connected, (H1) ↔ (H2), same rank and cluster count; apply T to
  both sides. **CONFIRMED.**
- **Theorem C′** (limit argument): the bookkeeping of the ε-flux into R_U, the propagation u_{c′} = O(ε) along paths
  to R_U and the limit flux supported in U are correct. **CONFIRMED** (and subsumed by C^gen).
- **U / U***: U has the steady state (k₁/k₃, 1, k₁/k₂, 1) for every k; U* has f_{K₁} − f_{K₂} > 0 when κ₂ ≥ κ₁ (exact).
- **Credits, checked against the sources on disk.** Pérez Millán–Dickenstein–Shiu–Conradi 2012 (arXiv:1102.1590):
  Theorem 3.8 (existence of positive toric steady states under Conditions 3.1, 3.4, 3.6) and Corollary 3.12 (a
  partition valid for all rate constants) are as P describes; P's reduction (Σ has columns A_k[·,c] on sources and 0
  on non-sources; partition {sources} ∪ singletons) is right. Their Remark 3.9 also contains the Brouwer +
  no-boundary-steady-state existence argument used in U(a) and U*(b). Campos–Chiola–Silva 1991 (author's copy):
  Lemma 5.1 (N_R mono-T-semiflow for live, structurally bounded FC nets), Theorem 5.3 (positive visit ratios
  independent of the marking) and Theorem 5.4 with its proof (under the rank condition rank(C) = m − 1 − (a − n) the routing system (5) is solvable
  for every conflict-resolution policy) are as P describes; the rates there are conflict-resolution
  probabilities/rates of a stochastic interpretation (P says "rational rates": acceptable). Add (L2): Recalde–Júlvez–
  Silva IFAC 2002 and Júlvez–Recalde–Silva Automatica 2005 (see §8).

## 4. Agent E2, Props. E2.3–E2.5

- **E2.3** (WF-EC multistationary for some k iff some nonzero d ∈ S, w ∈ D^⊥ have equal signs; E2 writes L for D):
  both directions re-derived (the witness x₂ = d/(e^w − 1), x₁ = x₂e^w, u_c = x₁^{y_c}, k_t = v_t/u_{c(t)} works because
  w·y_c is constant). It is PMDSC 2012 Theorem 5.5 ((5.8)–(5.11), checked in the source) specialised by Theorem T, as
  E2 says. Its use: with Theorem U(c) the criterion never fires; my independent maximal-minor test agrees on all
  8 716 enumerated and 480 random WF-EC networks (0 multistationary). **CONFIRMED.**
- **E2.4** (log-cone criterion): steady state ⇔ v = k ∘ x^{Ye} ∈ V = ker Γ ∩ ℝʳ_{>0} ⇔ log v − log k ∈ im Yeᵀ = K^⊥;
  absolute consistency ⇔ h(v) = (c^l·log v)_l onto ℝ^d. For d = 1: relint of the flux cone = positive combinations of
  all extreme rays, h = Φ(log λ) + O(1) uniformly, Φ piecewise linear on the order cones generated by ±1 and the
  indicators of initial segments; h onto iff Φ takes both signs. Correct. **CONFIRMED.**
- **E2.5 / E2.5′** (Birch): for v⁰ ∈ V and any k, Birch's theorem gives a unique v > 0 with v − v⁰ ∈ K and
  log v − log k ∈ K^⊥; K ⊆ ker Γ puts v in V. For X = N*: cᵀYe = 0 ⇒ cᵀYp = 0 when every product complex is in the span
  of the sources ⇒ cᵀΓ = 0. Correct; N* has (H2) because N has (H1). **CONFIRMED.** (Elementary; I did not search for
  a prior statement of E2.5 in this form.)

## 5. Exact re-checks of the examples (own scripts)

`scripts/check_examples_v2.py` (sympy + fractions, 1.6 s, 0 FAIL): N0, N0*, N1, N1*, N2, N2*, NT, NT*, NTc, NTc*, NM,
NM*, NB, NB*, U, U*, 3A ⇄ 2A+B and its dual (drain species V: V̇ = V(κ_AU³V − κ_B)), E's 4-cycle and its dual, NP, NP*,
and the (H1)–(H3)/connectedness of both sides for every row of §7.2.
`scripts/check_k8s_v2.py` (0.1 s, 0 FAIL): M3/M4, M1/M2.
`scripts/cr_counterexample.py` (30 s, 0 FAIL): NC. `scripts/cr_counterexample2.py` (62 s, 0 FAIL): NC3.

## 6. Independent computations and agreement of counts

- **Exhaustive re-enumeration with P's conventions** (`scripts/enum_wfec_v2.py`, own generator and canonical form):
  (n, molecularity, r) = (2, 3, ≤5): EC 75, WF-EC 71 (62 impure, S ≠ D 34); (3, 2, ≤5): EC 1 143, WF-EC 958 (793
  impure, S ≠ D 415); (3, 3, ≤4): EC 9 399, WF-EC 6 995 (6 230 impure, 52 unbalanced, S ≠ D 5 303); (4, 2, ≤4):
  EC 1 057, WF-EC 692 (436 impure, S ≠ D 486). These agree with P's logs (EC 75 / WF 71), Q3 (71 with 62 impure; 958
  with 793 impure; 6 995 with 6 230 impure and 52 unbalanced; 692 with 436 impure) and E2's table (L ≠ S: 34, 415,
  5 303, 486). On every one of the 8 716 WF-EC networks: T (two random k), C, C′, U(b) (two random non-steady
  points), U(c), U*(a) (two random points) PASS. Runtime 2 s, 21 s, 310 s, 198 s.
- **Random families** (`scripts/random_wfec_v2.py`, seeds 3–6): see §2.3, §2.4, §2.8, §2.10.
- **T*′**: §2.6. **GS**: §2.5.
- **Search for small CR counterexamples.** A first untargeted search (`scripts/cr_search.py`, 3 × 600 s) met only 65
  non-PC cluster-rank networks and no witness (logs superseded). The targeted generator `scripts/cr_search2.py`
  (sources sharing species inside a block, products constrained to keep s ≤ q − 1, molecularity ≤ 3 or 4, 60 exact
  steady-state samples per network), 3 × 420 s: 3 772 networks (n ≤ 4, molecularity ≤ 3) with 157 witnesses, 3 094
  (n ≤ 5, molecularity ≤ 3) with 102, 3 679 (n ≤ 4, molecularity ≤ 4) with 125 (`logs/cr_search2_*.log`). A witness is
  a rational steady state with det(J|_S) of sign (−1)^{s+1}, which by the degree argument forces at least three
  positive steady states in its class when all are nondegenerate. Smallest found: n = 4, r = 4, molecularity 3 (NC3
  is one of them); no witness with n = 3 was met (with n = 3, non-PC forces q = 2 and s = 1).

## 7. Consistency of terms across P's, Q3's and E2's notes

- *Equal conflict*: P (Teruel–Silva Def. 6: •t ∩ •t′ ≠ ∅ ⇒ Pre[P,t] = Pre[P,t′]) and Q3 (two reactions sharing a
  reactant have the same source; equivalently distinct sources have disjoint supports) agree. Teruel–Silva's
  "equal conflict relation" (Pre[P,t] = Pre[P,t′] ≠ 0) gives the doubling classes; under (H2) these are the clusters.
- *Balanced joins*: P's definition (every source has all nonzero coefficients equal) is Teruel–Silva's (a join is
  balanced when all its input arcs carry the same weight, p. 273). Consistent.
- *Pure*: no species on both sides of a reaction (Teruel–Silva's Pre·Post = 0). Note that **Teruel–Silva assume
  purity without loss of generality** (p. 273), so "WF-EC ⇔ Petri well-formed" (P §0, via TS Thm 28) is literally a
  statement about pure nets; Q3 and the theorems above use WF-EC as the *algebraic* definition (EC + connected +
  (H1) + (H2) + s = m − 1), which needs no purity. The paper should define WF-EC algebraically and mention Petri
  well-formedness only as a remark for pure nets.
- *Rank equation*: P defines it as s = #clusters − 1 and notes that for EC with every species consumed it is
  s = m − 1; Q3 uses m for EC and q = m_cl for C^gen and PC; §7.2 separates m_src and m_cl. Consistent. (TS Thm 28
  writes rank(C) = |𝔈| − 1 with 𝔈 the equal-conflict sets.)
- *Well-formed PC*: "connected, (H1), (H2), s = q − 1" in §1.1, §4b and the summary: consistent.
- **Notation clash to fix in the paper**: E2 writes L = span{y_c − y_c′} (Prop. E2.3), while Q3 writes D for that space
  and L = (im Γ)^⊥ for the conservation laws.
- Q3 §4.1 cites "1 989 unbalanced EC networks with s = m of (3,3,4,0)"; E2's notes report "2 003 (n ≤ 3)". Different
  scopes, not a contradiction, but the paper should quote one source.
- P's summary still lists as "open" Theorem C without purity, Birch-type uniqueness and absolute consistency of N*
  for unbalanced joins; these are now Q3's C, U and T*′ (confirmed here). P's Corollary P ("N persistent if pure or
  balanced") is superseded by Q3's CP.

## 8. Attribution

- **Johnston's translated networks** (coordinator's question; arXiv:1305.5845, Bull. Math. Biol. 76 (2014); I read
  §4.1 with Def. 4.1, Def. 4.2 and Thm 4.1). A translation keeps a bijection of reactions and **keeps every reaction
  vector** (Def. 4.1(1)); reactions with the same source go to the same translated source (Def. 4.1(2)); Thm 4.1 needs
  a weakly reversible translation (proper, or strongly resolvable improper) with structural and kinetic deficiency
  δ̃ = δ̃_K = 0. Since Γ = Ỹ·I (I the incidence matrix of the translated graph), δ̃ = 0 means ker Γ = ker I, the cycle
  space of a directed graph, whose circuits are signed cycles with all nonzero entries ±1. **Counterexample:**
  W = {A → C, 2C → 2A} is WF-EC (s = 1, m = 2, (H1)–(H3)), and ker Γ = ℝ(2, 1) is not a cycle space; so W has no
  translation with δ̃ = 0 and Thm 4.1 does not apply, while Theorems T, U, GS do (exact, `logs/translation_check.log`).
  Among 300 random WF-EC networks, 207 have a circuit of ker Γ with unequal absolute entries and hence no δ̃ = 0
  translation. When a δ̃ = 0 translation exists, Johnston's Thm 4.1 gives the toric parametrisation (overlapping the
  toric half of Theorem T, already credited to PMDSC) and uniqueness per class only under σ(S) = σ(S̃) and
  (+,…,+) ∈ σ(S^⊥), a stronger hypothesis than U(c); it says nothing about the Jacobian at non-steady points or global
  stability. Not checked: the later extensions (Tonello–Johnston, Bull. Math. Biol. 2018; Johnston–Müller–Pantea,
  Bull. Math. Biol. 2019) that replace δ̃ = 0 by an effective deficiency; they could cover more WF-EC networks for the
  steady-state parametrisation (Theorem T), not for U(b) or GS.
- **Júlvez–Recalde–Silva, Automatica 41 (2005)** (on disk) and Recalde–Júlvez–Silva IFAC 2002. I confirm the
  coordinator's reading: infinite-server semantics f[t] = λ[t]·min_{p∈•t} m[p]/Pre[p,t] (piecewise linear, not mass
  action); the routing relation (18) f[tᵢ]/(Pre[P,tᵢ]·λ[tᵢ]) = f[tⱼ]/(Pre[P,tⱼ]·λ[tⱼ]) for transitions in continuous
  equal conflict and Definition 7 (MTSR: visit ratios independent of the initial marking) are the analogue of the flux
  half of Theorem T; the paper states that its results apply to systems that, after a transient, reach a steady
  state (p. 607, after Definition 1), i.e. convergence is assumed, not proved; Fig. 2 shows non-monotone
  throughput. I did not locate the exact sentence on non-unique steady-state markings in §3.1, which does not change
  the conclusion: **no predecessor for Theorems U or GS under mass action there.** Credit them (with Campos–Chiola–Silva) for the flux/routing half of T.
- **Desel–Esparza Thm 10.16–10.17** (book on disk): I checked the statements (10.16: weakly connected, positive S-
  and T-invariants, Rank = |C_N| − 1 ⇒ well-formed, credited to Desel; 10.17: then a marking is live and bounded iff it
  marks all proper siphons). That the ordinary-net siphon part of C^gen follows ("every siphon contains an
  S-component") uses the proof of 10.17 (its step (i)); I did not re-derive that step. I agree with L2 that the
  ordinary case is essentially known; the weighted statement, the rate-independent persistence reading and K5gen were
  not found.
- **Theorem GS**: the monotonicity of max/min of h_c = ψ_c/u_c is the classical consensus/compartmental argument (my
  identity ḣ_c = Σ_{d≠c} G[c,d](u_d/u_c)(h_d − h_c) makes this explicit; Maeda–Kodama–Ohta 1978, Moreau 2004/2005 as
  listed by L2 — not opened by me). The application to weighted EC mass-action networks (non-WR, non-complex-balanced)
  was not found by L2 or by me (no search beyond the sources on disk and Johnston). Theorem U(b)'s Metzler/null-vector
  argument is the same Perron–Frobenius device.
- Donnell–Banaji 2013 and node-balanced steady states: Q3's statement that they do not contain WF-EC is plausible;
  not re-checked in depth.

## 9. Corrected statements for the paper

1. Theorems C, C^gen, K5gen, U, GS, T*′, U*, GS*≤2, PC, Corollaries CP and BH, Props. T*-cex, T*-dim, tor*: as in
   `round2/ec/notes.md`, unchanged.
2. Cluster-kinetics generalisation (Q3 §2 and §3 Remark (iii)): add "∂φ_c/∂x_p > 0 on the open orthant for p ∈ S_c"
   for (b) and for GS (strict monotonicity suffices for (a) only).
3. Table of §7.2, row 4: "N* loses the species dual to X₁+2X₂ → 3X₂" (k₃ in the path order of the table, k₄ in E's
   order). NM: "critical siphons {X₁,X₃} and {X₂,X₃}"; NM is PC.
4. Conjecture CR: replace by "Proposition (siphon part) = Theorem C^gen; Proposition (CR fails for uniqueness and
   stability): NC3 = {D → B+C, C → A, A+B → B+C, 2A+B → C+D} is TEFC, connected with (H1)–(H3) on both sides,
   s = q − 1 = 2, and for k = (1/6, 1/3, 6, 6) has exactly three positive steady states in the class B + D = 10/3,
   A + C + D = 19/2, two stable and the saddle (1/2, 1/3, 6, 3). (Also NC = {C → B, 2A → D, 2B+D → 2A+2B,
   2B+2D → 2A+2C+D}, k = (2/3, 5/18, 24, 64), class B + C = A + 2D = 7/2.) Open: absolute consistency under the
   cluster rank equation." Theorem PC is then the boundary of the dynamical results.
5. Notation: one letter for span{y_c − y_c′} (D) across Q3 and E2.
6. Define WF-EC algebraically; "Petri well-formed" only for pure nets (TS assume purity).
7. Credits: add Recalde–Júlvez–Silva 2002 / Júlvez–Recalde–Silva 2005 (min semantics) to Theorem T's flux half;
   Desel–Esparza Thm 10.16–10.17 (Desel 1994) for the ordinary case of C^gen; consensus/compartmental systems for the
   GS technique; state that Johnston 2014 does not cover WF-EC in general (W).

## 10. What I could not check

- GS* for s ≥ 3 with unbalanced joins and GS for general PC networks: open; only simulations (mine: 40 networks).
- The large counts of E2 (464 331 networks for E2.4, 11 185 WF-EC for E2.3) and Q3's runs on E's general families
  were not re-run; I re-implemented the checks on independent families (§6) instead.
- Q3's §7.3 scan count (737 non-EC cluster-rank networks all EC or BH up to catalysts) was not re-run.
- The proof step in Desel–Esparza 10.17 giving "every siphon contains an S-component"; Desel 1994 itself.
- The later translation frameworks (Tonello–Johnston 2018, Johnston–Müller–Pantea 2019); a full literature search
  for Theorems U, GS, PC (beyond L2's table and the sources on disk).

## Files

- `scripts/v2lib.py` — own exact library (rational linear algebra, flux/conservation rays by support enumeration,
  siphons, clusters, EC/PC tests, A_k, u(k), exact Jacobians, Faddeev–LeVerrier, Hurwitz determinants, maximal minors,
  random WF-EC generator).
- `scripts/check_examples_v2.py` → `logs/check_examples_v2.log`; `scripts/check_k8s_v2.py` → `logs/check_k8s_v2.log`.
- `scripts/random_wfec_v2.py` → `logs/random_wfec_v2_seed{3,4,5,6}.log`.
- `scripts/enum_wfec_v2.py` → `logs/enum_wfec_n2_m3_r5.log`, `logs/enum_wfec_n3_m2_r5.log`, `logs/enum_wfec_n3_m3_r4.log`,
  `logs/enum_wfec_n4_m2_r4.log`.
- `scripts/tstar_v2.py` → `logs/tstar_v2.log`; `scripts/gs_sim_v2.py` → `logs/gs_sim_v2.log` (numerical).
- `scripts/cr_probe.py`, `scripts/cr_counterexample.py` (NC), `scripts/cr_counterexample2.py` (NC3, msolve),
  `scripts/cr_search.py` (superseded), `scripts/cr_search2.py` → `logs/cr_probe.log`, `logs/cr_counterexample.log`,
  `logs/cr_counterexample2.log`, `logs/cr_search2_*.log`.
- `scripts/translation_check.py` → `logs/translation_check.log`.
- Johnston's paper was read from a scratch copy outside the repository (arXiv:1305.5845).
