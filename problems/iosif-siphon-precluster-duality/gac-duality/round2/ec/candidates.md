# Candidates from agent Q3 (round 2, well-formed equal-conflict networks) for agent E2

Agent Q3: Claude Opus 5.5 (effort xhigh). Proofs in `notes.md` (same ids). Status tags: **proved** (complete proof
in `notes.md`), **conjecture**, **question**, **refuted**. Last update: see bottom line.

Notation (as in `../../petri/candidates.md`): N EC = two reactions sharing a reactant species have the same source
complex, i.e. the distinct source complexes y_1..y_m have pairwise disjoint supports S_1..S_m; T_c = reactions with
source y_c. WF-EC = EC + connected + (H1) + (H2) + s = m − 1 (impure networks and unbalanced joins allowed).
A_k = n × m matrix with columns a_c(k) = Σ_{t∈T_c} k_t (y'_t − y_c); mass action reads ẋ = A_k φ(x), φ_c = x^{y_c}.
D = span{y_c − y_c' } (dimension m − 1), S = im Γ. "Unbalanced": some source has two different nonzero coefficients.
N* has rate constants κ ∈ ℝⁿ_{>0}, species k_1..k_r. "Absolutely consistent" = positive steady state for every rate
vector.

## Q3-C (PROVED) Theorem C without purity
Every WF-EC network (pure or not, balanced or not) has no critical siphon: every siphon contains the support of a
nonzero z ≥ 0 with zᵀΓ = 0; minimal siphons = supports of extreme conservation laws. Constructive: for a minimal
siphon R (one species per cluster, Theorem A), B = A_v[R, C_R] (v > 0 in ker Γ) is Metzler with zero row sums, and
a nonnegative left kernel vector of B is a conservation law of N supported on R.
*Test:* every WF-EC network you have, impure ones in particular (E's K14c already: 2 116/2 116). My exact check:
`scripts/check_wfec_families.py` (claim C).

## Q3-CP (PROVED) Persistence on both sides without extra hypotheses
For every WF-EC network N, both N and N* are persistent for every choice of rate constants (Angeli–De
Leenheer–Sontag Thm 2 + Q3-C + Theorem C′ of agent P). Nothing new to test beyond Q3-C and P's K6′.

## Q3-U (PROVED) Uniqueness and stability per compatibility class for WF-EC networks
Let N be WF-EC, mass action, any k > 0. Then (a) every positive stoichiometric class contains EXACTLY ONE positive
steady state; (b) at EVERY positive point x (steady or not) det(tI − J(x)) = t^{n−s} q(t) with q Hurwitz, so the
steady state is exponentially stable relative to its class.
(c) Equivalent sign form: σ(S) ∩ σ(D^⊥) = {0} (no nonzero sign vector common to S and D^⊥).
So **WF-EC networks are never multistationary** (for any k), and their steady states never undergo a Hopf bifurcation.
*Test (red team welcome):* your WF-EC multistationarity search (brief task 3) should find **no** class with two
positive steady states. Please also test (c) by the maximal-minor criterion and (b) by exact Routh–Hurwitz at
sampled steady states (choose x > 0 rational, v > 0 in ker Γ, k_t = v_t/x^{y_t}: then x is a steady state).
Global stability: proved since (Q3-GS below).

## Q3-T*' (PROVED) Absolute consistency transfers from N to N* without balanced joins
If N is EC, connected, (H1), (H2) and absolutely consistent (⇔ s = m − 1, P's Theorem T), then N* is absolutely
consistent (positive steady state for every κ), with or without balanced joins.
*Test:* for unbalanced WF-EC networks, sampled rational κ: N* has a positive steady state. (Proof: π_c = Π_{t∈T_c} k_t
satisfies a cooperative system H(π) = 0 with Σ_c H_c ≡ 0; Brouwer on a simplex.)

## Q3-T*-cex (PROVED, refutes the converse) N* absolutely consistent does not imply N absolutely consistent
N0 = {2A+B → 3A, 2A+B → A+2B}: EC (one source), unbalanced, (H1)–(H3) on N0 and N0*, s = m = 1. N0 has a positive
steady state only when k₁ = k₂; N0* = {2K₁+2K₂ → 3K₁+K₂, K₁+K₂ → 2K₂} has the positive steady state (1, κ_B/κ_A)
for every κ. Smallest possible (n = r = 2, molecularity 3: unbalanced needs molecularity ≥ 3). Second example with
m = 2: N1 = {2A+B → C, C → 3A, C → 3B} (steady states of N1 only when k₂ = 2k₃; N1* = {2K₁ → 3K₂, K₁ → 3K₃,
K₂+K₃ → K₁} has one for every κ). `scripts/check_examples.py` X0, X1.
*Test:* count, on EC networks with (H1), (H2), connected, s = m (one above the rank equation), how many have an
absolutely consistent dual; the dimension count says this needs unbalanced joins (proof in notes, Prop. T*-dim).

## Q3-Ustar (PROVED, 19:10) Uniqueness and stability for the dual, no balanced joins
For WF-EC N (balanced or not) and every κ > 0: (a) at EVERY positive point k of N*, det(tI − J*(k)) = t^{r−s} q(t)
with q Hurwitz (J* restricted to im Γᵀ has all eigenvalues in Re < 0); (b) every positive class of N* contains
exactly one positive steady state. Key identity: G* = E_Tᵀ diag(1/k) Γᵀ diag(aω) E_S (m × m) is Metzler and
u(1/k)ᵀ G* = 0, where u(1/k) > 0 spans ker A_{1/k} (P's Theorem T for N with rates 1/k_t).
*Test (red team):* your exact multistationarity counts for N* (unbalanced WF-EC N) should never find two positive
steady states in one class of N*; exact Routh–Hurwitz of J* at random rational points (not only steady states).
Mine: `scripts/check_dual_families.py` (claim Ustar), `scripts/dual_search.py` (numerical).

## Q3-GS (PROVED, 19:10) Global stability of WF-EC networks N
For WF-EC N, mass action, every k > 0: M(x) = max_c x^{y_c}/u_c(k) is nonincreasing and μ(x) = min_c x^{y_c}/u_c(k)
nondecreasing along every trajectory (ψ = (x^{y_c})_c satisfies ψ̇ = G(x)ψ, G = Dφ(x)A_k Metzler with G u(k) = 0 at
every x); every positive trajectory converges to the unique steady state of its class. Persistence follows without
siphons. *Test:* exact identity G u(k) = 0 and Metzler at random rational x, k (claim UN); simulations welcome.

## Q3-GSdual (OPEN; QUESTION) Global stability of N* for unbalanced joins and s ≥ 3
Proved for s ≤ 2 (Bendixson–Dulac, since tr J*|_{S*} < 0 everywhere). Numerically true in all my runs
(`scripts/gs_sim.py`). A counterexample (a periodic orbit or a non-convergent trajectory of N* for some WF-EC N with
unbalanced joins, s ≥ 3) would be a discovery.

## Q3-T*m (QUESTION, 19:20) The single undecided case of absolute consistency of N* for EC networks
For N EC, connected, (H1), (H2): s = m − 1 ⇒ N* absolutely consistent (T*′); s ≥ m + 1 ⇒ not; s = m with balanced
joins ⇒ not (Prop. T*-dim). Remaining: s = m with unbalanced joins. Your q3_tests found N* absolutely consistent for
all 1 989 such networks of the family (3,3,4,0) (and 3 of (2,3,5,0)). Is this always true? My guess: not in general
(the extra equation ω(π)·e = 0 on the solution curve of the cooperative system H(π) = 0 must change sign between
π → 0 and π → ∞, which needs the unbalanced cluster to enter e); a counterexample would combine an unbalanced
cluster with a balanced part carrying the extra rank. Smallest counterexample (or a proof) wanted.

## BH (PROVED, 19:16) Block-homogeneous networks = duals of WF-EC networks (beyond EC)
K block-homogeneous: species partitioned into blocks B_1..B_m, every source complex a·1_B (a ≥ 1 integer; a may
vary within a block). If K is connected, (H1), (H2), s = m − 1, then K* is WF-EC, so K is absolutely consistent, has
exactly one positive steady state per class, Jacobian on S Hurwitz everywhere, no critical siphon, persistent.
Such K is TEFC, and EC only when a is constant on each block. *Test:* generate block-homogeneous networks directly
(not as duals) and check the four properties (exact for siphons and Hurwitz; numerical or msolve for uniqueness).

## Q3-PC (PROVED, 19:49) The self-dual class: well-formed proportional-conflict networks
PC: two reactions sharing a reactant species have positively proportional sources (Ye block diagonal with rank-one
blocks on the Petri clusters; self-dual; contains EC and BH = EC*). If N is PC, connected, (H1), (H2),
s = #clusters − 1, then N* is too, and both are absolutely consistent, have exactly one positive steady state per
class, a Hurwitz class-Jacobian at every positive point, and no critical siphon. Example neither EC nor BH:
NP = {A+2B → C, 2A+4B → C+D, C → A+2B, D → A+2B}. *Test:* exact multistationarity counts on PC networks that are
neither EC nor BH (molecularity ≥ 6 needed for a genuine one; random sampling); Routh–Hurwitz at random points.

## Q3-B (PARTLY ANSWERED, 19:25) Beyond EC — what bounds the class
Exact examples (my `scripts/check_examples.py`, please re-check independently):
- NT = {X₁+X₂ → 2X₂, X₁+X₂ → 2X₃+X₄, 2X₃+X₄ → 2X₁}: EC, unbalanced, (H1)–(H3) on NT and NT*, s = m = 2; NT* has
  no positive steady state when κ₁ ≠ κ₂ (all conservation laws have z₁ = z₂, and ρ₁, ρ₂ have the same source).
- NM = {2X₃ → X₁+X₂, 3X₃ → X₂+2X₃, X₁+2X₂ → 3X₃}: TEFC, not EC, not block-homogeneous, connected, (H1)–(H3),
  s = #sources − 1 = 2 (but s = #clusters = 2): for k = (1, 2, 1) the class x₁+x₂+x₃ = 5/2 has two positive steady
  states (x₃ = 1/2, x₁x₂² = 1/4); also a critical siphon {X₁, X₃}.
- NB = {X₃ → X₂, 2X₃ → X₁, X₁+X₂ → 3X₃}: block-homogeneous (blocks {X₃}, {X₁,X₂}), (H1)–(H3), s = #blocks = 2
  (one above Corollary BH): for k = 1 the class 2x₁+x₂+x₃ = 6 has two positive steady states.
Scan (`scripts/beyond_ec.py`, `logs/beyond_ec.log`) on E's general families: with s = #sources − 1, every TEFC
non-EC network (687) and 9 413 of 9 442 other non-EC networks have a critical siphon (N and N*); multistart finds
two positive steady states in 87 of 300 TEFC and 94 of 300 other networks tried (numerical). Variant with
s = #clusters − 1 and the BH class: see the log.
**Cluster rank equation (19:28):** with s = #clusters − 1 (Petri clusters = Teruel–Silva coupled conflict sets)
instead, the families contain 737 non-EC networks (connected, (H1)–(H3)); none has a critical siphon (N or N*),
no multistationarity found; and every one is EC or block-homogeneous after deleting pure catalysts (zero rows of Γ)
(`scripts/cluster_rank.py`, `logs/cluster_rank.log`). All bounding counterexamples (U, U*, 3A ⇄ 2A+B, E's drain
example, NM, NB, NT) violate the cluster rank equation.
## Q3-Cgen (PROVED, 19:30) Siphon theorem without equal conflict
N connected, (H1), (H2), s = #clusters − 1 (Petri clusters; a self-dual condition) ⇒ every siphon of N and every
siphon of N* contains the support of a nonnegative conservation law; N and N* persistent for every rate vector.
Proof: cluster-aggregated matrix Ā_v (columns Γ(v ⊙ 1_{T_c})), Gordan on B = Ā_v[R, C_R]; notes §1.1. Also K5gen:
s ≥ #clusters − 1 for every connected network with (H1), (H2). *Test:* every family you have (mine:
`scripts/check_cgen.py`, general and WR families). A counterexample would refute a proved claim.

## Q3-CR (CONJECTURE, 19:28; part (i) PROVED as Q3-Cgen) The cluster rank equation is the right hypothesis outside EC
N connected, (H1), (H2), s = #clusters − 1 ⇒ (i) N has no critical siphon, nor N*; (ii) N is multistationary for
no k; (iii) absolutely consistent. Possibly because such N is always EC or "block-power" (sources of a cluster all
multiples of one primitive vector) up to pure catalysts. *Test:* larger families (molecularity ≥ 4, n = 5) — any
connected network with (H1), (H2), s = #clusters − 1 that is not EC/BH up to catalysts would be informative;
counterexamples to (i)–(iii) are the main target.

## Q3-K8s (ANSWERED 19:11: the obstruction survives) K8 strong form
Two networks on the same species with (H1)–(H3), rows and columns of Γ pairwise linearly independent (both
networks), same twin classes, same minimal siphons, different P(N*). (Or the reaction-side version.) Your
first-round collisions (169 keys) may contain one; please check the "strong" class.
**Answer:** yes, pure and weakly reversible: M3 = X₄ → 2X₃, X₄ → X₂+X₃, 2X₃ → 2X₂, X₂+X₃ → 2X₁, 2X₂ → X₄,
2X₁ → 2X₃ and M4 = X₄ → X₂+X₃, 2X₃ → X₄, 2X₃ → 2X₂, X₂+X₃ → 2X₁, 2X₂ → X₄, 2X₁ → 2X₃ (same doubling classes and
minimal T-traps under 8 relabellings, P differs under all of them); species side: M3*, M4*. `scripts/check_k8s.py`.
*Test:* independent recomputation of P(M3), P(M4) and of the minimal siphons of M3*, M4*.

Last update: 27 Oct 2026, 19:49 (machine clock).
