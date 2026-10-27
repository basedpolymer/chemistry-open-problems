# Verification V3: persistence results and obstructions of the gac-duality work

**Verifier V3: Claude Opus 5.5 (effort xhigh), fresh context, 27 October 2026.** Brief: `../../BRIEF-V3.md`
(scope extended by the coordinator, with the owner's approval, to agent A2's notes `round3/dichotomy-astra/`).
I did not take part in the work. Paths below are relative to `gac-duality/` unless stated otherwise.

Method (PROTOCOL §4): every proof marked "proved" in scope was re-derived by hand; every explicit example used as a
counterexample or certificate was re-checked with my own exact code (`verification/V3/scripts/v3lib.py`, written
here, imports no author code; sympy for identities); the authors' scripts were re-run with logs written only in
this folder; the finite counts of Q1 were recomputed with an independent generator and classifier; the cited
sources on disk were compared with the statements that cite them. Environment: WSL Fedora 44, Python 3.12.15,
sympy 1.14, scipy (numerical illustrations only, never used as evidence for a "proved" claim).

## 0. Verdicts at a glance

| Agent | Claim | Verdict |
|---|---|---|
| G | Props 2.1, 2.2 (rank, boundedness, q = β* + δ*, dual WR by a row permutation, δ* = q − c(π)) | CONFIRMED |
| G | §2 endotacticity row tests; complex balance remarks | CONFIRMED |
| G | §3 obstruction M = 3A ⇄ 2A+B, D = M* (unit rates, class T = 2) | CONFIRMED |
| G | §4 small codimension (cycle, rank-preserving operations), Prop. 4.1 (prime n), 5-species example, 4-species counterexample, Gale-dual remarks | CONFIRMED |
| G | §5 ADLS criterion not self-dual (H, H*) | CONFIRMED |
| G | Props 6.1, 6.2, 6.3 (siphons = T-traps, face codimension = dual subnetwork rank, Wiuf's (a), (⋄) in dual matrices) | CONFIRMED (one reading note on Wiuf's text) |
| G | Prop. 7.1 (first-order dual class) | CONFIRMED |
| G | §8 Example 8.1 (J_S = [[−15,27],[24,−25]], det −273) and the Blanchini–Giordano reading | CONFIRMED |
| E | Prop. E.4 (drain species: extinction for every rate vector) | CONFIRMED |
| E | Example X₄ → 3X₃ → X₁+2X₂ → 3X₂ → X₄ (s = 3, dual loses k₄ for all rates); strengthening of G §3 to all rates and all classes | CONFIRMED |
| Q1 | Lemmas 1.1, 1.2; Theorem S1 (1)–(3) | CONFIRMED |
| Q1 | Theorem S1 (4) ("decided by finitely many polynomial inequalities") | CONFIRMED WITH CORRECTIONS (wording) |
| Q1 | S1-N, S1-WR (dichotomy at s = 1), C-2 | CONFIRMED |
| Q1 | C-pure | CONFIRMED WITH CORRECTIONS ((H1) must be assumed; 0 ⇄ A) |
| Q1 | Theorem D (drain sets) and its translation to N | CONFIRMED |
| Q1 | Examples X, M, P, S, B, Q, C | CONFIRMED (exact) |
| Q1 | rank-one table (8 families) | CONFIRMED (independent recount, all numbers equal) |
| E2 | Lemma E2.0, Prop. E2.1 (exact decision at s = 1) | CONFIRMED |
| E2 | Prop. E2.7 (vertex scaling) and its example | CONFIRMED |
| E2 | refutation of Q1-6 (X₄ ⇄ X₂+X₃, 3X₃ ⇄ X₁+2X₄); generic SOME example | CONFIRMED |
| D2 | Prop. WR-Σ | CONFIRMED |
| D2 | **Theorem W** (no DOM facet of N* for WR N, every rank) | CONFIRMED (true for every facet set, siphon or not; only "products of Z-reactions are sources" is used) |
| D2 | Lemma 2.1, Lemma A | CONFIRMED |
| D2 | Lemma BM | CONFIRMED (= stable half of Butler–Waltman 1986 Thm 4.1 with M = {e}) |
| D2 | **Theorem P2** (ALL certificate) | CONFIRMED WITH CORRECTIONS ((A2) needs V_e(e) = 0; isolation step to be written out) |
| D2 | Prop. R, Prop. H, Prop. B | CONFIRMED |
| D2 | Prop. E, bicoloured case | GAP in the written proof (sign change along one edge of one class not justified); repaired here for μ = 1, which covers all 451 NONE verdicts that rest on it |
| D2 | Prop. E zero-flow case (μ = 1) | CONFIRMED |
| D2 | Props CM, L, S | CONFIRMED at sketch level (standard local arguments; their example computations re-derived exactly) |
| D2 | NONE certificates E (example), H, CM, L, S, B examples; Example U | CONFIRMED (exact) |
| D2 | counts (ALL 1 885 / NONE 12 769 / undecided 171 of 14 825) | CONFIRMED by re-run; ALL verdicts pass my independent exact necessary test |
| D2 | Conjecture F (status "checked", open for one-way flows) | status accurate, but the conjecture is **REFUTED** by A2's Theorem A |
| A2 | Theorem A (ALL dual with an OUT facet) | CONFIRMED |
| A2 | Theorem AF (family N_(u,R)) | CONFIRMED |
| A2 | Theorem E / E+ | CONFIRMED; existence part = E2's Cor. E2.5′ (same hypothesis) |
| A2 | Proposition L and boundary-equilibrium remark | CONFIRMED |
| A2 | **Theorem O** (unique equilibrium per class, hyperbolic source, periodic orbit) | CONFIRMED |
| A2 | §7 rejected SOME candidate (NONE by a Hurwitz boundary vertex) | CONFIRMED |

No claim marked "proved" in my scope is refuted. Two statements need a hypothesis or wording fix (C-pure, S1(4)),
one theorem needs a missing condition made explicit (P2), and one sketch has a gap that I could close for every case
in which it is used (Prop. E, bicoloured). D2's Conjecture F is false (A2).

## 1. Agent G (`gac/notes.md`)

**Prop. 2.1 — CONFIRMED.** (1) transposition and rank–nullity. (2) Recession cone S ∩ ℝⁿ≥0 = {0} iff S^⊥ meets
ℝⁿ>0 (Stiemke/Gordan); 0 ⇄ A has unbounded classes, its dual V → U is conservative (checked, `G2.1(2)`).
(3) With n dual reactions, rank I* = m* − ℓ* (loops give zero columns), so β* = n − m* + ℓ* and β* + δ* = n − s = q
identically; in a WR graph each linkage class carries a cycle, so β* ≥ ℓ*; q = 1 forces ℓ* = 1, δ* = 0 and a
strongly connected graph with E = V, i.e. one directed cycle. Counts (2,1,1,0) → (4,2,1,1) for M → D and (2,1,1,0)
for (A ⇄ B)* re-checked exactly. Random test RW4: 1 953 random networks (`logs/v3_random.log`).

**Prop. 2.2 — CONFIRMED.** If N is WR every product is a source, so Yp = Ye·Σ (D2's WR-Σ) and equal rows of Ye give
equal rows of Yp, hence equal rows of Γ, excluded by (H3); the same for Yp via the reversed network. In-/out-degrees
of the dual graph are then ≤ 1, and WR ⇔ every vertex has in = out = 1 ⇔ the two row sets coincide; m* = n,
ℓ* = c(π), δ* = q − c(π); 1_C (C a π-orbit) is a conservation law since Σ_{i∈C}(b_i − a_i) = 0. Random test RW5
(2 237 WR networks with (H3), half of them built from a coordinate permutation; 1 270 have a WR dual).

**Endotacticity tests (§2) — CONFIRMED.** G's row tests are equivalent to Gopalkrishnan–Miller–Shiu, Def. 3.14 (read
in `gac/certificates/sources/gms-2014.txt`, p. 11–12): "all maximisers over all sources have d ≤ 0, one has d < 0" is
equivalent to "endotactic along z and some w-maximal source over all sources points strictly inward".

**§3 — CONFIRMED.** Exact re-check (`v3_examples.py`, claims G3.*): D = M* = {3U+2V → 2U+3V, V → U}; both Γ equal
[[−1,1],[1,−1]]; (H1)–(H3) on both sides; M has (m,ℓ,s,δ) = (2,1,1,0), D has (4,2,1,1); D not WR; w = (0,1)
witnesses non-endotacticity; v′ = u³v² − v; h′(v) = (2−v)²(2−4v); h(1/8) = 3375/4096, so v′ ≤ −(721/4096)v on
(0,1/8]; (1,1) is an equilibrium of the class T = 2, which has exactly two positive equilibria (v = 1 and
v ≈ 0.16071 ∈ (1/8,1/2)). E's strengthening (every rate vector, every class): v′ = v(k₁u³v − k₂) ≤ −(k₂/2)v for
v ≤ min(T, k₂/(2k₁T³)) — CONFIRMED. M: a′ = a²(k₂b − k₁a), globally attracting on every segment.
Note for the paper: Q1's X = 2X₂ ⇄ X₁+X₂ is a smaller instance (see §3 below).

**§4 — CONFIRMED.** The monomolecular n-cycle has q = 1, s = n − 1 and a dual n-cycle. Prop. 4.1: π is one n-cycle,
every linkage class is a full T-orbit (orbit length divides the prime n and is not 1 because no reaction is zero),
Φ_n is irreducible (Eisenstein on Φ_n(t+1)), so d, Td, …, T^{n−2}d are independent and each linkage class spans
S = V; a w not orthogonal to S is non-constant on every class; following the cycle from a global maximiser gives a
strictly inward reaction at the global maximal height and every reaction from that height points weakly inward
(its target is a source). With G's equivalent form of GMS Def. 3.14 this is strong endotacticity, and GMS Thm 1.1/1.2
give permanence and GAC. (H3) is indeed redundant (dual kernel = span(1)). Five-species example: 10 reactions,
(m,ℓ,s,δ) = (10,2,4,4), dual one 5-cycle with δ* = 0, both classes span S, {X₁,X₂,X₃} a critical siphon realizable at
(0,0,0,2,3) with face codimension 3, J₀ = ∅ and dim S₀ = 1 (Wiuf (a) fails) — all exact; strong endotacticity also
passes 3 000 random integer w (corroboration). Four-species example: q = 1, (H3), N and N* WR; w = (3,2,1,2) gives
heights (3,2,1,2,4,4) and displacements (−1,−1,1,1,0,0): no strictly inward reaction at the maximal height, so not
strongly endotactic — exact. Gale-dual remarks: elementary and correctly limited.

**§5 — CONFIRMED.** H* = {2U₁+2U₃+U₄ → 2U₂+U₃+2U₄, 2U₂+U₄ → 2U₁+U₃}; H is WR with (H1)–(H3), its only minimal
siphon {A,B} supports a + b; {U₄} is a siphon of H* with row (1,−1) of Γ*, H* has (H1)–(H3).

**Prop. 6.1 — CONFIRMED** (re-expression of Angeli–De Leenheer–Sontag 2007, Thm 2, read in `adls-2007.txt`; GAC
under complex balance from persistence is classical). **Prop. 6.2 — CONFIRMED**; the Anderson–Shiu reformulation is
right because Anderson–Shiu index faces by the maximal zero set, int(F_W) = {x : x_i = 0 ⇔ i ∈ W} (their (2.6)–(2.7),
Thm 4.6), which is G's "realizable". **Prop. 6.3 — CONFIRMED.** Reading note: in `gac/certificates/sources/wiuf-2026.txt`
(pdftotext) condition (⋄) appears as "S_ℓ ∩ ℝ^{|Σ|}_{≥0} = {0}" because the "≠" is lost; the PDF (p. 5, rendered and
read) has **S_ℓ ∩ ℝ^{|Σ|}_{≥0} ≠ {0}**, and Remark 3.3 there says (⋄) means that no Σ-minimal class has a strictly
positive conservation law of its own. G's (6.4) "ker M[J_L,Σ] ∩ ℝ^Σ_{>0} = ∅" is therefore the correct translation
(Stiemke). The identification J₀ = Wiuf's Σ-blind reactions uses WR (G assumes it).

**Prop. 7.1 — CONFIRMED** (each species lies in at most one source with coefficient 1; no zero complex in a
nontrivial conservative WR class; out-degree 1 + WR gives disjoint cycles on disjoint groups; the siphon/conservation
law and (H3) arguments are correct).

**§8 — CONFIRMED.** Example 8.1: Γ, rank 2, δ = 0, support condition (Ye)_ij > 0 ⇔ Γ_ij < 0, dual WR with δ* = 0,
alternative realisation with the same Γ, G = Dg(1,1,1) = [[10,0,2],[2,1,0],[0,10,1]], and in the basis
(e₁ − e₃, e₂ − e₃) J_S = [[−15,27],[24,−25]] with Γ·G·B = B·J_S and det −273 — exact. G's description of
Blanchini–Giordano's Assumptions 3–4 matches their p. 2 (`certificates/sources/Blanchini-Giordano-2022-LCSS.txt`).

## 2. Agent E (`computation/notes.md` §4–5)

**Prop. E.4 — CONFIRMED.** Production terms of X_i have source exponent ≥ d₀ + 1, the draining reaction contributes
−c₀x_i^{d₀}, so x_i′ ≤ x_i^{d₀}(Cx_i − c₀) on a bounded class; below c₀/(2C) the decay is strict and the limit is 0.
Remark (not in E): the hypothesis can never hold for a WR network N itself, by a direct argument — the cycle through
the source d₀e_i returns to it, so some reaction raises the X_i-coefficient to d₀ from a source with coefficient < d₀,
contradicting (b). E's "X-drainN" check is thus a theorem.

**Example X₄ → 3X₃ → X₁+2X₂ → 3X₂ → X₄ — CONFIRMED.** (m,ℓ,s,δ) = (4,1,3,0), (H1) (w = (1,1,1,3)), (H2), (H3); the dual
is E's up to relabelling of the dual species and satisfies (H1)–(H3); exactly one drain species (d₀ = 1); N has
none. Global attraction of N: Anderson 2011 (single linkage class) and Pantea 2012 (s ≤ 3), both statements checked
in the sources on disk. E's re-run logs: `logs/rerun_E_g_examples.log`, `logs/rerun_E_e_examples.log` (all PASS).

## 3. Agent Q1 (`round2/persistence/notes.md`)

**Lemmas 1.1, 1.2 and Theorem S1 (1)–(3) — CONFIRMED.** I re-derived the endpoint expansion, the converse in
Lemma 1.2 (F is a polynomial in ε along the class), (1) (all IN ⇒ σL < 0; not IN ⇒ large κ_{j₀} at
x⁺ = (0 on Z, 1 off Z)), (2) (OUT; ESC via a separating w and x⁺ = e^{tw}), and (3) (weighted AM–GM; strict linear
inequalities in κ; their complement contains an open set). The sign conventions (upper face Z ⊆ N₋, inward ⇔
σc_j < 0) are consistent throughout and with E2's "good/bad" convention.

**Theorem S1 (4) — CONFIRMED WITH CORRECTIONS.** Only example B is given. The sentence "persistence for some κ is then
decided by finitely many polynomial inequalities" should read: "the set of κ for which M is persistent is
semialgebraic (Tarski–Seidenberg; the condition is a first-order formula in κ and the boundary point); it can be a
hypersurface (example B)".

**S1-N — CONFIRMED.** With Γ = abᵀ, the dual has g = b, c = a; a face is a set of same-orientation reactions of N; the
Z-order is the load ℓ_Z(i); inward ⇔ X_i is a net product of Z.

**S1-WR — CONFIRMED.** I re-derived the proof: on each linkage class (a segment y_L + ta) the affine function
φ_o − φ_v has slope a_o − a_v < 0, the extreme forward position is the only source of reactions outside Z, every
reaction of Z leaves a position where φ_o − φ_v > 0, contradicting equal minimal loads. Neither (H1) nor (H3) is
used; self-loops do not affect the argument. Independent evidence: no DOM face among the duals of the 8 026
networks of Q1's eight families (my generator), nor among 2 669 random WR rank-one networks (molecularity ≤ 8) and
2 745 random non-WR networks satisfying only the weaker hypothesis (RW1, RW2).

**C-pure — CONFIRMED WITH CORRECTIONS.** The proof needs species of both signs in a, i.e. (H1). Without (H1) it is
false: N = 0 ⇄ A is pure, weakly reversible, consistent, s = 1; its dual k₂ → k₁ is conservative and loses k₂ for
every rate (OUT face) (`V3.Cpure`, exact). Corrected statement: *if s = 1, N satisfies (H1) and (H2), and no
reaction of N has one of its net products in its source, then N* is persistent for every rate vector.*

**C-2 — CONFIRMED** (one species per sign class, loads add over Z, so singletons decide).

**Examples — CONFIRMED (exact, `v3_examples.py`).** X: dual k₂ → k₁, 2k₁+k₂ → k₁+2k₂, k₂′ = k₂(κ₂k₁² − κ₁), face {k₂}
ESC, verdict NONE, no drain species and no drain set (the only pure-source dual reaction is k₂ → k₁ and (D2) fails
with ratio 1). X is the smallest WR network with (H1)–(H3) and a NONE dual (n = 1 is impossible under (H1) with
s ≥ 1; X₁ ⇄ X₂ is the only molecularity-1 candidate and is ALL; my (2,2,4) enumeration finds X first). P: every
singleton face IN, the face of the two forward reactions ESC, verdict NONE. S: k₁′ = k₂((κ₂ − κ₁)k₁ + 2κ₃k₂), DOM
only at one end. B: k₁′ = (κ₂ − κ₁)k₁k₂, DOM at both ends. Q: (k₂ − k₁)′ = (κ₃ − κ₂)k₃ (the Barbalat step is
correct). C: V′ = −κ₁k₃ − κ₂k₂ + 2κ₃k₁²k₂k₃, {k₂,k₃} a drain set (w = (1,2), d = (1,1)), no drain species.
Q1's own script re-run: `logs/rerun_Q1_examples.log` (all PASS).

**Theorem D — CONFIRMED.** The estimate x_Z^β ≤ U^θ ≤ U^{θ−1}Σu_i (u_i = x_i^{d_i}, θ > 1), forward invariance of
{V ≤ δ}, V → 0, and "Z contains no semiflow support / ω-limit zero sets are critical siphons" are correct; the
translation to N (private reactant of multiplicity d_j, condition on net products) is correct.

**Rank-one table — CONFIRMED by an independent recount** (`scripts/v3_rank1_enum.py`, own generator, own canonical
form, own classifier; `logs/rank1_*.log`): (2,2,4) 12/4/8/0, (2,3,6) 460/115/345/267, (2,4,6) 3 544/381/3 163/1 833,
(2,5,5) 2 742/662/2 080/1 601, (3,2,6) 12/4/8/0, (3,3,6) 426/101/325/130, (3,4,5) 829/245/584/291, (4,3,4) 1/1/0/0
(networks / ALL / NONE / NONE with an OUT face), 0 SOME — identical to Q1's table. Drain-set counts (D-count) not
re-checked; the s = 2 drain-species count 199 was re-checked independently (`logs/v3_drain_count.log`).

## 4. Agent E2 (`round2/computation/notes.md`)

**Lemma E2.0, Prop. E2.1 — CONFIRMED.** Same content as S1 (IN = (1), OUT/ESC = (2), DOM = (3)); the Newton-vertex
argument in (2) is correct (a vertex of the hull of the tied exponents that is a bad and not a good exponent). The
HYPERPLANE case is immediate.

**Prop. E2.7 — CONFIRMED.** (iii) gives a vertex v with zero set exactly Z in the class of u = v − εd; (ii) makes v
an equilibrium for every κ; only sources e_j + w_k contribute to the linearisation, which is Σκ_kT^{|w_k|}v^{w_k}M_k
in the class of Tu; trace and the mixed-determinant expansion (det M_k = 0) are basis-invariant; the extreme group
dominates. Example 2X₃ → X₂ → X₁+X₃ → 2X₃: linearisation [[−κ₂, κ₁ + κ₃T²],[0, κ₃T² − κ₁]], eigenvalues −κ₂ and
κ₃T² − κ₁ (exact).

**Refutation of Q1-6 — CONFIRMED.** N = X₄ ⇄ X₂+X₃, 3X₃ ⇄ X₁+2X₄ is pure, WR, (H1) (w = (1,1,3,4)), (H2), (H3), s = 2,
δ = 0; its dual (equal to E2's up to relabelling) has a two-species drain set with d = (1,1) (weights (3,1) in E2's
labels; my search finds it). **Generic SOME example — CONFIRMED**: s = 1, (H1), (H2), not WR; the only non-IN face is
the DOM face {k₂}; there σL = −y₁(3y₃² − 2y₃ + 3)/2 < 0 for all y₃ at κ = (1,1,3/2,3/2) and σL = −y₁(y₃² − 3y₃ + 1) > 0
at y₃ = 1 for κ = (3,3,1,1).

## 5. Agent D2 (`round3/dichotomy/notes.md`)

**Prop. WR-Σ — CONFIRMED** (also the identity k_j′ = G_{y′_j} − G_{y_j} on 40 random WR networks, RW6, sympy).

**Theorem W — CONFIRMED.** Proof re-derived: φ(y) = y_o − Σλ_qy_q vanishes on sources outside Z and drops by
b_j(a_o − Σλ_qa_q) < 0 along each reaction of Z; a source with φ ≤ 0 would start an infinite strictly φ-decreasing
chain of distinct complexes; hence φ > 0 on all sources of Z, contradicting Σ_Zφ(y_j) = μ − μ = 0. The proof uses
neither the siphon property nor full WR, only that every product of a reaction of Z is a source. Checks: D2's claim W
re-run (0 DOM siphon facets in 18 861 duals); my RW3 test on random WR networks, all facet sets, siphons or not:
525 networks of rank 2 (2 150 facet sets) and 1 088 of rank 3 (7 903 facet sets), no DOM (`logs/v3_random.log`).

**Lemma 2.1, Lemma A — CONFIRMED.** In Lemma A the siphon edge {k_Z = 0} is invariant in both time directions for
the polynomial field (every monomial of a reaction consuming Z vanishes on the hyperplane; reactions avoiding Z do not
produce Z), so the total trajectory of p stays on it.

**Lemma BM — CONFIRMED.** It is the "stable" half of Butler–Waltman 1986, Theorem 4.1 (with their Lemma 2.1), for
the isolated invariant set M = {e}, and of Hale–Waltman 1989, Lemma 4.3/Remark 4.4 (semiflows) — read in
`certificates/sources/`. D2's proof is BW's; forward continuity plus backward uniqueness of the ODE gives
t_n − s_n → ∞, which is enough for the semiflow on P.

**Theorem P2 — CONFIRMED WITH CORRECTIONS.** The chain argument is a planar version of Butler–Waltman's acyclicity
theorem (BW Thm 3.1, HW Thm 4.1). It cannot be quoted from BW, who assume an invariant boundary; Lemma A is the right
replacement. Corrections: (i) in (A2) add **V_e(e) = 0** (or V_e → 0 at e); without it "V_e increasing ⇒ no interior
point converges to e" fails. D2's functions a^λb^μ and λa + μb satisfy it. (ii) The isolation of e (no compact
invariant set in cl U_e other than {e}) should be written out: for an interior point x of such a set K, V is
monotone along its total orbit, so α(x) and ω(x) contain no interior point; boundary orbits in cl U_e are monotone
along the edges (finitely many boundary equilibria), so α(x) = ω(x) = {e}, impossible with V strictly increasing.
(iii) The theorem as stated allows OUT facets with a one-way flow (the code is stricter: it requires IN); this is
correct, and A2's Theorem A is exactly such a case with a non-monomial V_e. Conclusion is persistence (not uniform
persistence).

**Prop. R, Prop. H, Prop. B — CONFIRMED.** (H: Cauchy–Binet for a sum of rank-one 2×2 matrices; a Hurwitz vertex
attracts an open subset of its class; B: decoupled blocks.)

**Prop. E — GAP in the bicoloured case; zero-flow case CONFIRMED.** The written proof finds a sign change of the edge
flow F "along exp(tθ + u), u in a segment"; these points lie on edges of *different* classes, so the sign change of
F along the edge of one fixed class is not established (the edge of a class through exp(tθ) may be too short in the
needed scale, and its direction d need not be transverse to {F = 0}). Repair (valid for μ = 1): by the intermediate
value theorem on that path there is p* with zero set exactly Z, F(p*) = 0 and normal coefficient L(p*) < 0; p* is an
equilibrium of its own class (Lemma 2.1); the edge is invariant, so the linearisation in (edge, normal) coordinates is
triangular with eigenvalues ∂_τF(p*) and L(p*) < 0; the (strong) stable manifold of the eigenvalue L(p*) is
transverse to the edge and enters the interior, whatever the sign of ∂_τF(p*). So some interior trajectory converges
to p*, for every κ. For μ ≥ 2 a centre-manifold or Ważewski argument with a genuine sign change along the edge would
be needed. Impact (`logs/v3_propE_split.log`): of D2's 13 817 NONE verdicts, 451 rest only on bicoloured Prop. E and
all have μ = 1 (36 rest on the zero-flow case); so every NONE verdict survives with the repaired argument.

**Props CM, L, S — CONFIRMED at sketch level.** Standard centre-manifold, LaSalle and sector arguments; I re-derived
the worked examples exactly: CM (a′ = κ₂Tb − κ₃Ta² + …, b′ = −(κ₁+κ₂T)b + 2κ₃Ta² + …, reduced coefficient
κ₃T(κ₂T − κ₁)/(κ₁ + κ₂T)), L (V′ = −2κ₁k₃² − 2κ₂k₂² + 2κ₃k₁²k₂k₃, negative definite iff κ₃²T⁴ < 4κ₁κ₂), E example
(facet {k₃} ESC, zero edge flow, μ = 1), ALL example (k₃′ = k₃(κ₃k₁² − 2κ₁k₃); edge siphons: only {k₃}; vertex siphon
{k₂,k₃}). I did not re-prove the general tropical conditions implemented in `d2cert.py`.

**Example U — CONFIRMED.** Exact identities z′ = a[(2κ₁ − κ₂) + Kz(1 − 2z)], u′ = a[(κ₂ − 2κ₁)u² − Ku + 2K],
a′ = a²(K − κ₂u), b′ = b²(−2κ₁ + 2Kz²); the three regime arguments re-derived (κ₂ > 2κ₁ in every class; κ₂ < 2κ₁ and
κ₂ = 2κ₁ for κ₃T³ small); log-coordinate numerics agree. The "3-cycles" remark (i) re-derived (det = ΣK_iK_i′(1 − r_i
− r_i′)); (ii) not checked (it is labelled evidence by D2).

**Counts — CONFIRMED.** `certify_s2.py` re-run from a copy that writes only into this folder
(`scripts/certify_s2_v3copy.py`, diff = output path only): output and undecided list identical to D2's (W, K, E2a,
E2b, F, D all PASS; ALL 1 885 / NONE 12 769 / undecided 171). Independent exact necessary test of the ALL verdicts
(`scripts/v3_s2_sample.py`): for all 4 873 ALL networks × 3 random rate vectors × 11 class scalings (10^±6 isotropic
and random anisotropic), no boundary vertex equilibrium is Hurwitz in its class (PASS); the same test detects a
Hurwitz vertex in 796 of 3 000 sampled NONE networks, so it has power.

**Conjecture F, Q1-7 — status statements accurate.** D2 proved F for two-coloured flows and zero flows with μ = 1 and
left one-way flows open; A2's Theorem A (one-way OUT edge, ALL) refutes F in that open case. There is no conflict with
D2's finite data: A's network and every ALL member of A2's family contain the complex X₂+2X₃+RX₄ (R ≥ 3) of
molecularity ≥ 6, outside E's record families (molecularity ≤ 3).

## 6. A2 results (`round3/dichotomy-astra/notes.md`, GPT-6-Astra)

Re-run of A2's certificates with logs here (`logs/rerun_A2_*.log`, all PASS; `reproduce.py` itself was not run
because it writes into A2's folder) and my own checks (`scripts/v3_a2.py`, `logs/v3_a2.log`, all PASS).

**Theorem O — CONFIRMED (priority item).**
- Uniqueness: b′ − 2c′ = (2βc − α)b + 9δc³ (exact); at equilibrium 0 < c < α/(2β), b = 9δc³/(α − 2βc),
  a³ = 3δc(α + βc)/(γ(α − 2βc)), and these solve all three equations (sympy). Monotonicity: b′(c) =
  9δc²(3α − 4βc)/(α − 2βc)² > 0 and (a³)′(c) = 3δ(α² + 2αβc − 2β²c²)/(γ(α − 2βc)²) > 0 on the interval
  (with x = βc/α ∈ (0,1/2): 3 − 4x > 0, 1 + 2x − 2x² > 0). a + b + c runs from 0 to ∞, so each class has exactly
  one positive equilibrium, for every rate vector.
- At (α,β,γ,δ) = (17,400,70000,10⁶), (1,1,1/100) is an equilibrium (fluxes (17,4,7,1)), T = 201/100; Jacobian in
  (b,c) with a = T − b − c: [[−59,3658],[−25,79]], trace 20, det 86 789, eigenvalues 10 ± i√86689; H =
  [[93655,−290457],[−290457,13471234]] is positive definite (det 1 177 283 151 421) and JᵀH + HJ = 3 471 560·I — all
  exact.
- Logic: Theorem A ⇒ every ω-limit set of an interior point is a compact subset of the open triangle; the equilibrium
  is a hyperbolic source, hence an isolated invariant set with trivial stable set, so by Butler–McGehee it lies in no
  other ω-limit set; by uniqueness ω(x) contains no equilibrium; Poincaré–Bendixson in the invariant triangle gives a
  nonconstant periodic orbit. Correct.
- Numerical illustration (`scripts/v3_orbit_numerics.py`, logarithmic coordinates, Radau): from (1, 0.99, 0.02) the
  trajectory settles on a cycle of period ≈ 2.5126·10⁶ along which min c ≈ 10^−12.15 and min b ≈ 10^−19.5 (near the
  vertex (T,0,0)). The orbit is persistent but passes extremely close to the boundary; integration in linear
  coordinates with absolute tolerance 10⁻¹⁴ wrongly suggests convergence to (T,0,0) — a warning for any numerical
  persistence test on these duals (E2's LSODA remarks point the same way).
- Blanchini–Giordano 2022: **no conflict.** Their Theorem 1 transfers *structural* local stability in the
  stoichiometric class (S·G Hurwitz for *every* G in the class of admissible derivative matrices), under Assumption 3
  (a rate depends on x_i iff S_ij < 0). Both N (X₂ is a catalyst in X₁+X₂ → X₂+2X₃+3X₄) and N* (b is a catalyst in
  b+c → a+b) violate Assumption 3, and stability of one mass-action Jacobian is much weaker than structural stability
  (G's Example 8.1). The coordinator's reading is right.

**Theorem A — CONFIRMED.** Network certificates exact ((11,1,4,1) with value 12 on each complex, flux (1,1,1),
(H3), rank 2, δ = 0); dual and equations exact; §2.2: a′ = αb + βbc > 0 on the open edge a = 0, b′ = 2γa³c² + 3δc³ > 0
on b = 0, nonzero field at (0,T,0) and (0,0,T), open points of c = 0 have b(t) = b(0)e^{−αt} backward-unbounded — so
Ω ∩ ∂P ⊆ {e_T}; the only siphon facet is {c} and it is OUT. §2.3: the identity (A.1) is exact, the margin
γT³/8 − Eε ≥ γT³/16 holds (checked also on 3 000 random exact points of U with random rates and T), and the
last-entry argument is correct (entry points x(s_j) with V → 0 accumulate on the open edge c = 0 at b = ε, excluded by
§2.2). This is Theorem P2's scheme with the non-monomial repelling function V = c·e^{−βb/α}.

**Theorem AF — CONFIRMED.** For u = 3..7, R = 1..7: hypotheses ((uR+2, u−2, R+1, u−2) gives u(R+1) on each complex),
dual, identity (AF.1), facet {c} OUT; R ≥ 3 margin ≥ L/2 (sampled exactly); R = 2: c′ ≤ −δc² on T = min(1,
δ/((u−2)γ)); R = 1: c′ ≤ −(δ/2)c on T = min(1, δ/(2(u−2)γ)) (symbolic reduction and exact sampling; T^u ≤ T for
T ≤ 1). For R = 1, k₃ (= c) is also a drain species (E.4). Together with A this refutes Conjecture F; the dependence
on R is structural, as A2 says.

**Theorem E / E+ — CONFIRMED.** The coercivity argument (some π_i·v > 0 for v ≠ 0 because Πᵀw = 0 with w > 0 and Π has
full column rank), strict convexity and the identification Πᵀv = 0 ⇔ Γᵀv = 0 are correct; only (H1) is used.
Attribution: the hypothesis im Γ ⊆ im Ye is the same as E2's (im Yp ⊆ im Ye), and the existence part is exactly E2's
Corollary E2.5′ (absolute consistency of N*), proved there by Birch's theorem. New in A2: the description of the
whole positive equilibrium set as a smooth global graph of dimension r − s in logarithmic coordinates (classical
Birch/moment-map structure, as A2 notes). For the Theorem-A dual (r − s = 1) the equilibrium curve meets each class
once, consistent with Theorem O.

**Proposition L — CONFIRMED** (opposite active vectors on a siphon edge would give a nonnegative conservation law of N
vanishing on a complex, hence on the whole linkage class). Also the boundary-equilibrium remark. Random test: 1 393
single-linkage rank-two networks, no counterexample. Consistent with D2's 18 undecided "IN, two-way flow" networks,
which all have two linkage classes.

**§7 candidate — CONFIRMED (NONE).** Own reconstruction: six species, nine reactions (three 3-cycles), WR, every
complex of total 6, (H2), (H3), rank 2; p = (0,1,0,1,1,1,1,1,1) satisfies p + Γᵀz/4 > 0 for z = (0,−1,0,−2,0,0), is
an equilibrium for every κ, and its Jacobian in (k₁,k₃) is [[−κ₂,0],[κ₂,−κ₄]] (Hurwitz) — exact.

## 7. Attribution and novelty

- **Zhang–Gao–Dochain, IEEE TAC 2023 (DOI 10.1109/TAC.2023.3275880)**, read in `literature2/sources/z.txt` (Def. 2.5,
  2.11, 3.1, Thm 3.10). Their Thm 3.10: every 1-dimensional mass-action system that is w_I-endotactic for the 0/1
  indicator of every siphon I is persistent (for every rate vector; their w-endotacticity is two-sided). Comparison
  with S1: in a conservative rank-one network, non-siphon faces are automatically IN (a minimal-order reaction with
  source avoiding Z can only produce Z), so ZGD's hypothesis implies "every face IN". The converse fails: my example
  A+B → 2A, 2A → A+B, 3A+B → 4A is conservative of rank one with every face IN (persistent for all rates by S1(1);
  a′ = a(κ₁b − κ₂a + κ₃a²b)), but the w_{A}-maximal source 3A+B points outward, so it is not W_I-endotactic
  (`V3.ZGD`, exact). Hence: the "if" half of S1(1) is in print in a weaker form (ZGD 3.10, which also covers
  non-conservative 1-d networks); new in Q1/E2: the exact characterisation (the "only if" half of S1(1)), OUT/ESC ⇒
  non-persistent for every rate vector (S1(2)), the DOM rate dependence (S1(3)–(4)), S1-N/C-2 on N, and the
  dichotomy S1-WR. I agree with L2's table, with this refinement. No mention of duals in ZGD.
- Drain species (E.4), drain sets (Theorem D), Theorem W, the non-self-duality results: I know of no prior statement;
  L2 found none. Agarwal–Craciun–Deshpande–Jin (arXiv:2601.08801) cannot apply to N* as L2 says (N* is consistent,
  and a consistent δ = 0 network is WR).
- Lemma BM and P2: credit Butler–Waltman 1986 (Thm 4.1, Thm 3.1) and Hale–Waltman 1989 (Lemma 4.3, Thm 4.1); D2 does.
- Theorem E existence: credit E2's Cor. E2.5′ (and Birch); A2 credits Müller–Regensburger for the convex argument.
- Q1's "smallest example" X should be cited instead of G's M as the minimal non-self-duality example; G's M remains the
  example with an explicit unit-rate trajectory.

## 8. Corrected statements for the paper

1. **C-pure.** If s = 1, N satisfies (H1) and (H2), and no reaction of N has one of its net products in its source
   (e.g. N pure), then N* is persistent for every rate vector. ((H1) is necessary: 0 ⇄ A.)
2. **S1 (4).** If DOM faces occur at both ends, the set of rate vectors for which M is persistent is semialgebraic; it
   can be a hypersurface (example B).
3. **Theorem P2, (A2).** Every vertex equilibrium e admits a neighbourhood U_e and a C¹ function V_e ≥ 0 with
   V_e(e) = 0, V_e > 0 on int P and V_e′ > 0 on U_e ∩ int P. (Add the isolation argument of §5 above.)
4. **Prop. E (bicoloured).** State it for μ = 1 with the proof: the θ-LP gives, for every κ, a point p* with zero set Z,
   F(p*) = 0 and L(p*) < 0; p* is an equilibrium of its class whose stable manifold for the eigenvalue L(p*) enters the
   interior. (All 451 uses in D2's counts have μ = 1.)
5. **Theorem W.** May be stated for every facet set Z (not only siphon facets) and under the weaker hypothesis "every
   product of a reaction of Z is a source complex".
6. **Conjecture F** is false (A2, Theorem A; family AF): a realizable OUT siphon facet with one-way edge flow is
   compatible with persistence for every rate vector.
7. **Theorem E** (A2): present the existence part as E2's Corollary E2.5′ and the graph/manifold description as the
   addition.

## 9. What I could not check

- The general tropical certificates CM, L, S as implemented in `d2cert.py` (only their worked examples and the
  overall counts); the sampled necessary test covers the ALL verdicts only through hyperbolic vertex equilibria.
- D2's Remark on 3-cycles, part (ii).
- E's and Q1's large-family counts other than those listed (the WR families of E, D-count of drain sets, the 705-network
  family of G — G's and E's scripts were re-run, their counts not re-implemented).
- G's §9 statements about Craciun v3, OpenAI's preprints and mmislan's repository (outside the theorems in scope).
- Wiuf's Theorem 7.2 itself (only the translation of its hypotheses).

## 10. Scripts and logs (all in this folder)

| script | content | log |
|---|---|---|
| `scripts/v3lib.py` | own exact library: Γ, dual, WR, counts, siphons, LP feasibility, rank-one faces, facet types in any rank | – |
| `scripts/v3_examples.py` | 75 exact claims on the examples of G, E, Q1, E2, D2 (1.6 s) | `logs/v3_examples.log` |
| `scripts/v3_rank1_enum.py` | independent recount of Q1's eight rank-one families | `logs/rank1_*.log` |
| `scripts/v3_random.py` | RW1–RW6 (S1-WR, weak S1-WR, Theorem W, Props 2.1, 2.2, WR-Σ), seed 271027 | `logs/v3_random.log` |
| `scripts/v3_s2_sample.py` | exact Hurwitz-vertex test of D2's ALL verdicts at extreme class scales (93 s) | `logs/v3_s2_sample.log` |
| `scripts/v3_propE_split.py` | NONE verdicts resting only on Prop. E, by case and μ (32 s) | `logs/v3_propE_split.log` |
| `scripts/v3_drain_count.py` | drain-species count at s = 2 (199) | `logs/v3_drain_count.log` |
| `scripts/v3_a2.py` | A2's Theorems A, AF, E, L, O and §7, exact (3.6 s) | `logs/v3_a2.log` |
| `scripts/v3_numerics.py`, `scripts/v3_orbit_numerics.py` | numerical illustrations (not evidence): G's D, example U, X*, Theorem O's cycle | `logs/v3_numerics.log`, `logs/v3_orbit_numerics.log` |
| `scripts/certify_s2_v3copy.py` | D2's `certify_s2.py`, output redirected here | `logs/rerun_D2_certify_s2.log`, `logs/rerun_D2_undecided_s2.txt` |
| authors' scripts, run unchanged with `python -B` | G `check_gac_duality.py`, Q1 `examples.py`, D2 `examples_d2.py`, E `g_examples.py`, `e_examples.py`, A2 `certify_out_all.py`, `certify_secondary.py`, `cycle_candidate.py` | `logs/rerun_*.log` |
