# Splitting-invariant preclusters, toric steady states, and Iosif's duality redone

Agent Q2 — **Claude Opus 5.5 (effort xhigh)**, 27 October 2026, second round on result 17
(`problems/iosif-siphon-precluster-duality/`). Working notes; nothing here quotes the author's e-mail. Status tags:
**proved** (complete proof here), **checked** (exact computation on the finite families listed in §7, not a proof),
**refuted** (explicit counterexample, checked exactly), **open**.

## Summary

| # | Statement | Status |
|---|---|---|
| 1.2 | Conradi–Flockerzi clusters (Fact 2(b), Algorithm 1) and the thesis clusters (Def. 4.6) are the *proportionality classes of the coordinates of ker A_k*, where A_k is the species × source-complex rate matrix (the matrix Σ of Pérez Millán et al. restricted to source complexes); CF = meet over k | proved |
| 2.1 | Hierarchy (H2): P ⪯ P_cs ⪯ π ⪯ CF ⪯ κ_V for every positive flux V; all inclusions strict (B₃, R, B₂dup / K₃, XY4) | proved; checked |
| 2.2 | π = P_cs ∨ M(N*) (flat closure = class-span closure joined with full coupling) | proved; checked |
| 3.1 | **CF is invariant** under relabelling, duplicating a reaction, and splitting a choice (or subdividing a reaction) through a new intermediate species | proved; checked |
| 3.2 | P, P_cs and π are not: P fails duplication (N₄.₁₂) and splitting (B₃); P_cs and π fail duplication (B₂ → B₂dup) and splitting (B₂dup) | proved; checked |
| 3.3 | Any duplication-invariant extension X ⪰ P of Iosif's preclusters contains the full-coupling classes: X(N) ⪰ P(N) ∨ MIPS(N*) | proved |
| 4.1 | EC + connected + (H1) + (H2) ⇒ CF is one block (κ_V too, for every V); π can have several (K₃) | proved |
| 4.2 | Weakly reversible with deficiency 0 ⇒ CF = κ_V = linkage classes | proved; checked |
| 5.1 | dim ker A_k ≤ #κ(k), with equality iff ker A_k has a basis with disjoint supports (PMDSC Condition 3.1) iff the positive steady states are the toric set T_k = {x : x^{y_c−y_c′} = ψ_c/ψ_c′, c ~ c′} | proved; checked |
| 5.2 | **Toric criterion in terms of P̂ = CF**: N is toric for every k with one rate-independent partition iff rank A_V = m − #P̂(N) for every V ∈ C°; then the partition is P̂(N). In general rank A_V ≥ m − #P̂(N) | proved |
| 5.3 | "Internal cycles refine P̂" (Remark 23 for P̂) is the isolation property; it implies 5.2's condition (rate-free form of Iosif's thesis Thm 4.17) | proved; checked |
| 6.1 | **Theorem 7 redone for P̂ ∈ {π, CF, κ_V}**: MIPS(N*) ⪯ P̂(N) always; P̂(N) = MIPS(N*) ⇔ P̂(N) ⪯ MIPS(N*) ⇔ #blocks equal ⇔ every doubling pair is fully coupled | proved; checked |
| 6.2 | Species side: MIPS(N) ⪯ P̂(N*) always; equality ⇔ every twin pair is conservation-coupled; **holds for every network with (H1), (H3) whose product complexes are sources (all WR networks of Conjecture 22)** | proved; checked |
| 6.3 | WR networks: P̂(N) = MIPS(N*) ⇔ every linkage class is a directed cycle (contrast: manuscript Cor. 10(2), P ≠ MIPS(N*) at deficiency 0) | proved |
| 6.4 | Question 24 for P̂: transversality persists; every block of P̂(N*) lies in a minimal siphon ⇔ no twin pair and every MIPS lies in a minimal siphon | proved; checked |
| 6.5 | Iosif's diagram for P̂: "MIPS refine preclusters" true always; "preclusters dual to MIPS" ⇔ coupled twin (doubling) pairs; "internal cycles refine preclusters" ⇔ isolation ⇒ toric; on the species side of a twin-free network ⇔ the extreme conservation laws have disjoint supports | proved; counts in logs |
| 7 | CF is computable in practice: a sound lower bound (flat closure + two rules) and an upper bound (κ_V at generic and special V) coincide on 63 717 of 63 722 network sides tested; the 5 others are settled by Theorem 4.2 | checked |
| open | An exact finite algorithm for CF in general; a closed combinatorial description of CF for WR networks with δ ≥ 1 (Lemma 1.2 and Prop. 5.1 were not found in Pérez Millán et al. 2012, Conradi–Flockerzi 2012 or the thesis) | open |

**Answer to the brief in short.** The right splitting-invariant notion is the Conradi–Flockerzi cluster partition
P̂ = CF, read as the proportionality classes of the kernel of the complex rate matrix. It is forced: every
duplication-invariant extension of Definition 19 contains the full-coupling classes (3.3), which is exactly why the
manuscript's counterexamples R and T to "MIPS refine preclusters" disappear. For P̂, the refinement MIPS ⪯ P̂ holds on
both sides of the duality for every network, and equality holds exactly when every doubling pair (twin pair) is
coupled — so the duality of Conjecture 22 becomes an exact theorem with the same exceptional class as
Theorem 7(1) of the manuscript, and holds identically on the species side of every weakly reversible network
satisfying (H1)–(H3). The price is that this equality is built in: for those networks P̂(N*) *is* MIPS(N). The
content of P̂ is dynamical, not combinatorial: "internal cycles refine P̂" is Conradi–Flockerzi's isolation property
and gives toric steady states for all rate constants (5.2, 5.3), and P̂ is a single block for every connected
equal-conflict network (4.1), so the Petri-class structure of P's Theorem D(e) (the profile rule) concerns the
non-invariant Definition 19 only. On the siphon side (Question 24) nothing new appears: P̂ is transversal to the
minimal siphons wherever it differs from the MIPS (6.4).

## 0. Notation

Manuscript notation: N has species X₁..Xₙ, reactions R_t: y_t → y′_t (t = 1..r), Ye, Yp, Γ = Yp − Ye, s = rank Γ;
N* the dual; (H1) conservative, (H2) consistent, (H3) rows of Γ pairwise distinct. C = ker Γ ∩ ℝʳ_{≥0} (flux cone),
E its extreme-ray matrix with rows n_t. Under (H2), C° := C ∩ ℝʳ_{>0} ≠ ∅, and C° is the relative interior of C, i.e.
C° = {Eλ : λ > 0} (the only facets of C are coordinate hyperplanes and no coordinate vanishes on C); the rows n_t
represent the coordinate functionals of ker Γ (E spans ker Γ), so "n_s ∈ span{n_t : t ∈ B}" iff "every v ∈ ker Γ
vanishing on B vanishes at s".

- **Source classes** 𝒮(N): reactions grouped by source complex (= the components of the doubling graph G of
  Definition 19, singletons included). m = number of distinct source complexes; Q ∈ {0,1}^{r×m} the indicator.
- For k ∈ ℝʳ_{>0}: **A_k := Γ diag(k) Q** (n × m). Column c is Σ_{t: y_t = c} k_t (y′_t − c). The mass-action field
  is f(x) = A_k φ(x) with φ_c(x) = x^{c}; x > 0 is a steady state iff φ(x) ∈ ker A_k. (A_k is the matrix Σ of Pérez
  Millán–Dickenstein–Shiu–Conradi 2012 restricted to the source complexes; for EC networks it is agent P's A_k.)
  For V ∈ C° we use A_V (k = V; then 1 ∈ ker A_V because A_V 1 = ΓV = 0).
- M(N) = ray MIPS of N; M* := MIPS(N*) = classes of fully coupled reactions (manuscript Prop. 9).
- Partitions are compared by refinement: A ⪯ B iff every block of A lies in a block of B. A ∨ B is the join.

## 1. Definitions

**Definition 1.1.**
1. P(N): Iosif's preclusters (Definition 19; the three variants agree, manuscript Lemma 5).
2. P_cs(N): agent P's class-span closure (`span_closure_preclusters`): start from 𝒮(N); a block C is merged into a
   block B with |B| ≥ 2 when some s ∈ C has n_s ∈ span{n_t : t ∈ B}; iterate.
3. **π(N), the flat closure**: the finest partition coarser than 𝒮(N) all of whose blocks are *flats* of the dual
   matroid M*(Γ), i.e. sets B with [n_s ∈ span{n_t : t ∈ B} ⇒ s ∈ B]. Equivalently, every block B is all of [r] or
   the complement of the support of some v ∈ ker Γ. It exists and is computed by iterating the merge "C meets the
   closure of B ⇒ merge" (each merge is forced in every flat partition coarser than 𝒮, by monotonicity of closure).
4. **κ(k)** (for k with ker A_k ∩ ℝᵐ_{>0} ≠ ∅) and **κ_V** (V ∈ C°): source complexes c, c′ are equivalent iff the
   coordinate functionals u ↦ u_c and u ↦ u_c′ are proportional on ker A_k (for κ_V: equal on ker A_V, since
   1 ∈ ker A_V); lifted to reactions through Q.
5. **CF(N) := ∧_{V∈C°} κ_V** (meet: t ≈ t′ iff t ~ t′ in κ_V for every V ∈ C°). We call P̂ := CF the
   *corrected preclusters*. On the species side, P̂(N*) = CF(N*) (needs (H1) for N).

**Lemma 1.2 (what Conradi–Flockerzi and the thesis compute).** Assume (H2).
(a) For V, W ∈ C°: W/V (componentwise) is constant on source classes iff W = diag(V) Q ω with ω ∈ ker A_V ∩ ℝᵐ_{>0}.
(b) The clusters of Conradi–Flockerzi 2012 (Fact 2(b), Algorithm 1: rows j, l are clustered iff
n_jλ/n_jμ = n_lλ/n_lμ for all λ, μ ∈ Λ(E) = {λ ≥ 0 : Eλ > 0} satisfying the doubling relations) form the partition CF(N).
(c) For fixed k, the thesis clusters (Def. 4.6: equality of n_iλ/n_iμ for all λ, μ ∈ D, D = {λ ∈ Λ(E) : Eλ/k
constant on doubling sets}) form κ(k), and κ(k) = κ_V for every V ∈ C° of the form k∘Qψ (e.g. V = k∘a^{Ye} for a
positive steady state a). Every V ∈ C° is of this form (k = V, ψ = 1). Hence CF(N) = ∧_k κ(k).

*Proof.* (a) W ∈ ker Γ and W = diag(V) w with w = W/V; w is constant on source classes iff w = Qω, and then ΓW = 0
iff A_V ω = 0. (b) V = Eλ and W = Eμ range over C° (§0). The doubling relations say that W/V is equal on each doubling
pair, i.e. constant on source classes. By (a) the admissible W for fixed V form the set diag(V)Q(ker A_V ∩ ℝᵐ_{>0}),
a nonempty (it contains V) relatively open subset of the linear space diag(V)Q ker A_V; two linear functionals agree
on it iff they agree on ker A_V. So rows j, l are clustered iff ω_{c(j)} = ω_{c(l)} on ker A_V for all V. (c) For
V = k∘Qψ ∈ C° we have A_V = A_k diag(ψ), so ker A_V = diag(ψ)⁻¹ ker A_k, and ω_c = ω_c′ on ker A_V iff
u_c/ψ_c = u_c′/ψ_c′ on ker A_k. Since ψ ∈ ker A_k has positive entries, this says exactly that the two coordinate
functionals are proportional on ker A_k (the factor is then forced to be ψ_c/ψ_c′); it does not depend on ψ. The set D
of the thesis is {λ : Eλ = k∘Qψ, ψ > 0}, and the same argument as (b) applies. ∎

So the clusters are not a combinatorial refinement of preclusters but a property of the kernel of the complex rate
matrix — the object whose disjoint-support bases give toric steady states (§5). Pérez Millán et al. 2012 work with
partitions of the complexes such that ker Σ has a basis with these supports; Lemma 1.2 identifies the
Conradi–Flockerzi/thesis clusters with the coarsest candidate for such a partition (prior-art check: §8).

## 2. The hierarchy

**Proposition 2.1.** Under (H2): P ⪯ P_cs ⪯ π ⪯ CF ⪯ κ_V for every V ∈ C°, and CF ⪯ κ(k) for every k as in 1.1(4).

*Proof.* P ⪯ P_cs: an edge of G⁺ joins a doubling pair {i, j} (one source class) or adjoins s with
n_s ∈ span{n_i, n_j}; the P_cs block containing the class of i has ≥ 2 elements and is closed under its own span,
so it contains s. P_cs ⪯ π: every P_cs merge "s ∈ C, n_s ∈ span(B)" is a closure merge. π ⪯ κ_V (ratio argument of
Conradi–Flockerzi Fact 2(a), made class-wise): let ω ∈ ker A_V, w = Qω, W = diag(V)w ∈ ker Γ. Suppose w ≡ ρ on a
set B (true initially on source classes). Then u = W − ρV ∈ ker Γ vanishes on B, hence at every s with
n_s ∈ span(B), so w_s = ρ. If s lies in a block B′ on which w ≡ ρ′, then ρ′ = ρ. By induction along the merges
defining π, w is constant on every π-block, for every ω ∈ ker A_V; hence π ⪯ κ_V. CF ⪯ κ_V is the definition of a
meet, and CF ⪯ κ(k) by Lemma 1.2(c). ∎

**Proposition 2.2.** π = P_cs ∨ M*. Explicitly: the blocks of π are the P_cs-blocks with ≥ 2 elements and the
full-coupling classes formed by the P_cs-singletons.

*Proof.* The closure of a singleton {s} is its parallel class {t : n_t ∈ ℝn_s} = its full-coupling class (two nonzero
nonnegative proportional rows are positively proportional). A final P_cs block B with |B| ≥ 2 is closed (no s outside
has n_s ∈ span(B)); in particular it contains the full-coupling classes of its elements. So the full-coupling class of a
P_cs-singleton consists of P_cs-singletons. The partition described is therefore a partition into flats coarser than
𝒮, so π refines it; conversely every flat partition coarser than 𝒮 contains the P_cs merges (2.1) and the parallel
classes. ∎

**All inclusions are strict** (exact, `logs/check_examples.log`):

| | network | witness |
|---|---|---|
| P ≺ P_cs | B₃: X → Y_a (a = 1,2,3), Y_a → Z, Z → X (agent P's Remark 2.5) | P = {all but Z→X}, {Z→X}; P_cs = one block |
| P_cs ≺ π | R: X₁ ⇄ X₂ (no doubling pair) | P_cs = singletons; π = M* = one block |
| π ≺ CF | B₂dup = B₂ with Z → X duplicated (6 reactions); K₃ = complete digraph on X, Y, Z (WF-FC, WR, δ = 0); smallest in the families: n = 2, {2X₂→X₁+X₂, 2X₂→2X₁, X₁+X₂→2X₂, X₁+X₂→2X₁} | B₂dup: π = {R1–R4}, {R5, R6}, CF = one block; K₃: π = 𝒮 (3 blocks), CF = one block (Thm 4.1, sympy certificate) |
| CF ≺ κ_V (generic V) | XY4 = {X+Y→2X, X+Y→2Y, X+2Y→2X+Y, X+2Y→3Y} (n = 2, r = 4, s = 1, (H1)–(H3)) | κ_V is one block for generic V; CF = π = {R1,R2}, {R3,R4}, separated at V with V₁ = V₂, V₃ = V₄ (each source class then balances on its own) |

XY4 in detail: A_V has the two columns (V₁ − V₂)(1, −1) and (V₃ − V₄)(1, −1), whose sum is ΓV = 0. For generic
V ∈ C° the column is nonzero, ker A_V = ℝ1 and κ_V is one block; at V = (1, 1, 2, 2) ∈ C°, A_V = 0 and κ_V has two
blocks. So the thesis clusters κ(k) are one block for generic k and two blocks for k = (1, 1, 1, 1) (every point is
then a steady state), and CF = two blocks. The opposite jump also occurs: for Z₄ = {2X₂→2X₁, 2X₁→2X₂, X₁+X₂→2X₁,
X₁+X₂→2X₂}, κ_V = CF = {R1}, {R2}, {R3,R4} at generic V, but κ_V = {R1,R2}, {R3,R4} (coarser) at V = (1,1,2,2): a
special κ(k) can be coarser than CF. (Claim X6.)

## 3. Invariance

**Operations.** (O1) relabelling species or reactions. (O2) *duplication*: add a copy t′ of a reaction t (rate
constants k_t + k_t′ give the same dynamics). (O3) *splitting a choice*: for a source class c with source y and a
nonempty A ⊆ c, add a new species W, replace each t ∈ A by W → y′_t, and add t_A: y → W. With |A| ≥ 2 and A ≠ c this
splits a choice into two nested choices; every way of splitting a k-way choice into a tree of binary choices is a
composition of such steps; with |A| = 1 it subdivides a reaction. In each case there is a correspondence
φ: reactions(N′) → reactions(N) (t′ ↦ t; t_A ↦ an element of A), an isomorphism ker Γ(N′) ≅ ker Γ(N)
(v′ ↦ v with v_t = v′_t + v′_{t′}, resp. v ↦ (v, Σ_{t∈A} v_t)) mapping C°(N′) onto C°(N), and (H1), (H2) are
preserved (for (O3) set z_W = z·y).

**Theorem 3.1.** CF and the family (κ_V)_V are invariant: CF(N′) = φ⁻¹(CF(N)) (and t_A lies in the block of A).

*Proof.* (O1) clear. (O2): t′ joins the source class of t, m′ = m, and the column of that class in A′_{V′} contains
V′_tΓ_t + V′_{t′}Γ_t = (V′_t + V′_{t′})Γ_t; so A′_{V′} = A_{V} with V the image of V′. Hence κ′_{V′} = φ⁻¹(κ_V); as V′
ranges over C°(N′), V ranges over C°(N). (O3): the source classes of N′ are c′ = (c∖A) ∪ {t_A} (source y), A (source
W), and the others. Write S = Σ_{t∈A} V_t > 0. In A′_{V′} the row W has entries S (column c′), −S (column A), 0
elsewhere; on the old species the columns c′ and A add up to the column c of A_V, the other columns are unchanged. So
ω′ ∈ ker A′_{V′} iff ω′_{c′} = ω′_A and the vector ω with ω_c := ω′_{c′} lies in ker A_V. This isomorphism of kernels
shows κ′_{V′} = φ⁻¹(κ_V) with c′ and A in one class. ∎

**Proposition 3.2 (the combinatorial notions are not invariant).**
(i) P fails (O2): N₄.₁₂ = 2X₁ ⇄ X₁+X₂ has P = {R1}, {R2}; with R1 duplicated, P is one block (n_{R2} = n_{R1} + n_{R1′}).
P fails (O3): P(B₃) has two blocks; splitting the choice of X into X → Y₁, X → W, W → Y₂, W → Y₃ gives one block.
(ii) P_cs and π fail (O2): B₂ (one block) versus B₂dup (two blocks: n_{Z→X} + n_{Z→X′} = n_{X→Y₁} + n_{X→Y₂} is in the
span of the class {X→Y₁, X→Y₂} but neither duplicate is). They also fail (O3): splitting the choice {X→Y₁, X→Y₂} of
B₂dup behind W creates t_A with n_{t_A} = n_{X→Y₁} + n_{X→Y₂} = n_{Z→X} + n_{Z→X′}, which lies in the spans of both
classes, so π and P_cs become one block. (Checked: X2, X3, X4, X9.)

**Proposition 3.3 (full coupling is forced).** Let X be any rule assigning to each network with (H2) a partition of
its reactions with X(N) ⪰ P(N), invariant under (O2). Then X(N) ⪰ P(N) ∨ M*(N).
*Proof.* Let s, t be fully coupled in N: v_s = μ v_t on ker Γ. Duplicate t into t, t′. In N′ the coordinate functionals
satisfy f_s = μ(f_t + f_{t′}), so n′_s ∈ span{n′_t, n′_{t′}}; {t, t′} is a doubling pair, so s is adjoined to it in
P(N′) ⪯ X(N′) = φ⁻¹(X(N)). Hence s and t lie in one block of X(N). ∎

So any invariant repair of Definition 19 contains the MIPS of the dual (reaction side) and, applied to N*, the MIPS
of N (species side). CF is the canonical such repair: it is exactly the relation "the ratio W/V of two positive
steady-state fluxes is forced to be equal" (Lemma 1.2), of which Definition 19, P_cs and π are sound but incomplete
approximations. Note also that loops (y → y) at a source, harmless dynamically, change P but not π or CF.

## 4. Structure of CF in two classes

**Theorem 4.1 (equal conflict).** Let N be EC (two reactions sharing a reactant species have the same source),
connected, with (H1) and (H2). Then ker A_V = ℝ1 for every V ∈ C°; hence κ_V is one block for every V and CF(N) is
one block. Also κ(k) is one block for every k admitting a positive steady state.
*Proof.* N is strongly connected (agent P's Lemma SC) and every species is consumed ((H2)). For V ∈ C°, A_V is
agent P's A_k with k = V, and 1 ∈ ker A_V. P's Lemma CF (a nonzero nonnegative kernel vector is positive and spans
the kernel) gives ker A_V = ℝ1. ∎
So in the Petri class the corrected preclusters are trivial; the profile rule of P's Theorem D(e) and the
three-block π(K₃) describe artefacts of the non-invariant definitions. With Theorem T (agent P): for a connected EC
network with (H1), (H2), CF is one block and the positive steady states (when they exist) are toric with one class;
they exist for every k iff s = m − 1.

**Theorem 4.2 (weakly reversible, deficiency 0).** If N is weakly reversible with δ = 0 then for every V ∈ C°,
ker A_V = span{1_L : L linkage class}; hence κ_V = CF(N) = the partition of the reactions by linkage class.
*Proof.* All complexes are sources. Write Γ = Y I (complex matrix Y, incidence I) and L_V = I diag(V) Q (the
weighted Laplacian with out-flows V), so A_V = Y L_V. If A_V ω = 0 then L_V ω ∈ ker Y ∩ im I = 0 (δ = 0). So
ker A_V = ker L_V, which for a weakly reversible graph with positive weights is spanned by ℓ vectors u^L, positive on
the linkage class L and zero elsewhere. 1 ∈ ker A_V forces u^L ∝ 1_L. ∎

Corollary: in the GAC setting with δ = 0, P̂ is the linkage-class partition, independent of the rate constants. For
δ ≥ 1, A_V can have 1_L ∉ ker A_V and linkage classes can merge (question Q2-12 for E2).

**Proposition 4.3 (coupled doubling pairs).** If every doubling pair is fully coupled (e.g. no doubling pair), then
Q ker A_V = diag(V)⁻¹ ker Γ, and κ_V = π = CF = M* for every V ∈ C°.
*Proof.* For u ∈ ker Γ and t, t′ in one source class, u_t/V_t = u_{t′}/V_{t′} (full coupling), so
diag(V)⁻¹ ker Γ ⊆ Q(ℝᵐ) and it equals Q ker A_V. Two coordinates agree on diag(V)⁻¹ ker Γ iff the functionals are
proportional on ker Γ, i.e. full coupling. ∎

## 5. Toric steady states

**Proposition 5.1.** Let k admit ψ ∈ ker A_k ∩ ℝᵐ_{>0} (e.g. k has a positive steady state a, ψ = φ(a)). Then
dim ker A_k ≤ #κ(k) (number of classes of source complexes), with equality iff ker A_k = ⊕_{J∈κ(k)} ℝ ψ|_J iff ker A_k
has a basis with pairwise disjoint supports (Condition 3.1 of Pérez Millán et al. 2012 for the source complexes). In
that case the positive steady states for k are exactly
T_k = {x > 0 : x^{y_c − y_c′} = ψ_c/ψ_c′ for all c ~ c′ in κ(k)},
a toric (log-affine) set of dimension n − rank{y_c − y_c′ : c ~ c′} (nonempty iff the log-linear system is solvable,
PMDSC Condition 3.6; automatic when a positive steady state is known to exist).

*Proof.* K := ker A_k. In each class J the coordinate functionals are proportional with factors ψ_c/ψ_c′ (forced by
ψ ∈ K). Choosing one c_J per class, u ↦ (u_{c_J})_J is injective on K, so dim K ≤ #κ(k), and equality holds iff it is
onto iff each ψ|_J (ψ on J, zero elsewhere) lies in K, i.e. K = ⊕ ℝψ|_J. Conversely, if K has a basis b¹..b^d with
disjoint supports, each coordinate is a multiple of one basis coefficient, so each support lies in a class, different
supports give independent functionals, and the supports cover all source complexes (ψ > 0); hence the supports are
the classes and d = #κ(k). Finally φ(x) ∈ K ∩ ℝᵐ_{>0} iff φ(x)|_J ∝ ψ|_J with a positive factor for every J, i.e.
x^{y_c}/x^{y_c′} = ψ_c/ψ_c′ for c ~ c′. ∎

**Corollary 5.2 (toric criterion in terms of P̂).** Assume (H2). For every V ∈ C°, rank A_V ≥ m − #P̂(N). The
following are equivalent:
(i) rank A_V = m − #P̂(N) for every V ∈ C°;
(ii) for every k with ker A_k ∩ ℝᵐ_{>0} ≠ ∅, ker A_k has a basis with disjoint supports, with the same partition for
all such k.
Then κ(k) = κ_V = P̂(N) for all k, V, and every k with a positive steady state has the toric steady-state set T_k with
the partition P̂(N). Necessary condition: #P̂(N) ≥ m − s.
*Proof.* m − rank A_V = dim ker A_V ≤ #κ_V ≤ #P̂ (P̂ ⪯ κ_V). Under (i) equality holds throughout, so κ_V = P̂ (a
refinement with the same number of blocks) and 5.1 applies. Under (ii) the common partition is κ(k) for every k
(5.1), so P̂ = ∧κ(k) is that partition and dim ker A_V = #P̂. rank A_V ≤ s gives the necessary condition. ∎
Examples: EC networks (Thm 4.1: P̂ one block, rank A_V = m − 1 for all V) — this is the toric half of agent P's
Theorem T; WR δ = 0 (Thm 4.2: rank A_V = m − ℓ). XY4 satisfies Condition 3.1 for every k (rank A_k ≤ 1, m = 2) but the
partition jumps at special k (one block generically, two blocks when both classes balance), so (ii) fails there.

**Corollary 5.3 (isolation = "internal cycles refine P̂").** Every internal cycle (support of an extreme flux ray) of
N lies in a block of P̂(N) iff ker Γ = ⊕_{B∈P̂} (ker Γ ∩ ℝ^B) (each block is a separator of the matroid of Γ). This is
the isolation property of Conradi–Flockerzi (and of the thesis, Def. 4.8, with P̂ in place of κ(k)), and it implies
(i) of 5.2: steady states toric for every k with partition P̂ (rate-free form of the thesis' Theorem 4.17).
*Proof.* The first equivalence holds because the extreme rays span ker Γ. Given the splitting, for ω ∈ ker A_V the
flux W = diag(V)Qω splits into its parts W|_B ∈ ker Γ, so ω|_B ∈ ker A_V for each block B; each part contains
1|_B ≠ 0. Hence dim ker A_V ≥ #P̂ ≥ #κ_V ≥ dim ker A_V. ∎
Checked (K6): every network side with isolation satisfies dim ker A_V = #κ_V at a random V. The converse fails (XY4).
On the species side of a twin-free network (N* without doubling pairs) the two are equivalent: both say that the
conservation space has a basis with disjoint supports, i.e. that **the extreme conservation laws of N have pairwise
disjoint supports** (then these supports are the MIPS); in the logs the counts "isolation (D)" and "Condition 3.1 (D)"
coincide in every family.

## 6. The duality question redone

**Theorem 6.1 (manuscript Theorem 7 for P̂).** Let N satisfy (H2) and let P̂ be any of π, CF, κ_V (V ∈ C°). Then
M* = MIPS(N*) ⪯ P̂(N), and the following are equivalent: P̂(N) = MIPS(N*); P̂(N) ⪯ MIPS(N*); the two partitions have
the same number of blocks; every doubling pair is fully coupled (⇔ P(N) ⪯ MIPS(N*), manuscript Thm 7(1)).
*Proof.* M* ⪯ π (2.2) ⪯ CF ⪯ κ_V (2.1). Refinement plus equal cardinality is equality. If every doubling pair is
fully coupled, 4.3 gives κ_V = π = CF = M*. If P̂(N) ⪯ M*, the doubling pairs (which lie in one block) are fully
coupled. ∎

**Corollary 6.2 (species side, Conjecture 22 for P̂).** Let N satisfy (H1). Then MIPS(N) ⪯ P̂(N*), and
P̂(N*) = MIPS(N) iff every twin pair of N is conservation-coupled. If moreover every product complex of N is a source
complex and (H3) holds (in particular N weakly reversible with (H1)–(H3)), N has no twins (manuscript Cor. 9) and
P̂(N*) = MIPS(N).
*Proof.* Theorem 6.1 for N* (consistent by (H1)), whose doubling pairs are the twin pairs and whose full-coupling
classes are the ray MIPS of N (manuscript Prop. 9 with N** = N). ∎

**Corollary 6.3 (weakly reversible, reaction side).** Let N be weakly reversible (no loops). Then P̂(N) = MIPS(N*)
iff no complex is the source of two reactions (every linkage class is a directed cycle).
*Proof.* Manuscript Cor. 10(1): in a weakly reversible network a doubling pair is never fully coupled. Without
doubling pairs Theorem 6.1 applies. ∎

**Table 1 of the manuscript for P̂** (any of π, CF, κ_V; ray MIPS):

| formulation | for P̂ | counterexamples that survive |
|---|---|---|
| A3 M(N) ⪯ P̂(N*); B3 M* ⪯ P̂(N) | **always true** | none (R, T, N₄.₁₂ are repaired) |
| A1, A2, A4 (species side) | ⇔ twin pairs conservation-coupled; always for WR + (H1)–(H3) | U, F, W* (twins not coupled) |
| B1, B2, B4 (reaction side) | ⇔ doubling pairs fully coupled; for WR ⇔ cycles only | W, N₃, F (smallest: N₃, n = 2, r = 3) |
| R23 internal cycles ⪯ P̂(N) | ⇔ isolation (5.3) ⇒ toric for all k | XY4, Z₄; the 3-cycle T now satisfies it |

The manuscript's Corollary 10(2) (P(N) ≠ MIPS(N*) for every WR network of deficiency 0) is an artefact of
Definition 19: for P̂, a WR deficiency-0 network satisfies B1 iff its linkage classes are cycles, and then
P̂(N) = M* = linkage classes (Thm 4.2).

**Proposition 6.4 (Question 24 for P̂).** Let N satisfy (H1).
(a) A block of P̂(N*) containing a twin pair lies in no minimal siphon of N (agent P's Theorem A).
(b) Every block of P̂(N*) lies in a minimal siphon of N iff N has no twin pair and every ray-MIPS of N lies in a
minimal siphon. In particular this holds for twin-free networks without critical siphons (agent P's Prop. G), e.g.
every WR network with (H1)–(H3) and no critical siphon.
(c) Dually (under (H2)): every block of P̂(N) lies in a minimal T-trap (= minimal siphon of N*) iff N has no doubling
pair and every full-coupling class lies in a minimal T-trap.
*Proof.* (a) Theorem A of agent P. (b) With a twin pair, (a) applies to its block. Without, P̂(N*) = M(N) (6.2).
(c) is (b) for N*. ∎

Counts for the other siphon relations with P̂ = CF (exact where the bounds coincide; `logs/fam_*.log`):
"every minimal siphon of N lies in one block of P̂(N*)" (Q2) holds for 11 352 of the 11 409 general networks
(n = 3, molecularity ≤ 3, r ≤ 3) but fails for 2 628 of 7 076 WR networks (n = 4, molecularity ≤ 3) and 940 of 2 352
(n = 5, molecularity ≤ 2); for WR networks with (H1)–(H3) it says that every minimal siphon lies inside one MIPS.
Equality with the minimal siphons (Q1) holds for 1 449 of the 11 409. So Question 24 gains nothing from the
correction: wherever P̂(N*) differs from M(N) it is transversal to the minimal siphons, and where it does not, the
question is the MIPS–siphon relation of agent P's Prop. G.

**Remark 6.5 (Iosif's closing diagram, faithful reading for P̂).** Left column (N*): internal cycles of N*
→"Refine"→ preclusters of N*; middle arrow: preclusters of N* ↔ MIPS of N; right column: MIPS of N →"Refine"→ siphons.
For P̂:
1. "MIPS refine preclusters" is a theorem for every network (6.1, 6.2) — the arrow goes the other way from the one the
   manuscript tested in A2, and it is the one that duplication invariance forces (3.3).
2. "Preclusters dual to MIPS" holds iff twin pairs are coupled; always in the GAC setting (WR + (H1)–(H3)); on the
   reaction side iff doubling pairs are coupled (for WR: cycles only).
3. "Internal cycles refine preclusters" (Remark 23) is the isolation property; for N* of a twin-free N it says that
   the extreme conservation laws of N have disjoint supports. It holds for 11 347 of 11 409 general networks
   (n = 3) on the species side, but for 4 451 of 7 076 and 1 582 of 2 352 in the WR families with n = 4, 5. Where it
   holds, N* has toric steady states for all rate constants with partition P̂(N*) = M(N) (5.3) — this is the only
   "consequence for toric dynamical systems" we can attach to the duality of preclusters.
4. "MIPS refine siphons": unchanged (agent P's Prop. G: true when there are no critical siphons, e.g. WF-EC).

## 7. Computation

**Lemma 7.1 (two-block species rule; sound).** Let σ be a partition with σ ⪯ κ_V for every V ∈ C° (e.g. σ = π).
Suppose a species p is touched (Γ_{pt} ≠ 0) by reactions of exactly two blocks B₁, B₂ whose functionals
f_i(V) = Σ_{t∈B_i} Γ_{pt}V_t are not identically zero on ker Γ, and f₁ is sign-definite on the extreme rays. Then B₁ ∪ B₂
lies in one block of κ_V for every V ∈ C°.
*Proof.* ω constant on σ-blocks; row p of A_Vω = 0 reads ω₁f₁(V) + ω₂f₂(V) = 0 and, with ω = 1, f₁ + f₂ = 0
(other blocks contribute f ≡ 0 on ker Γ). f₁(V) ≠ 0 on C° = {Eλ : λ > 0} since f₁ ≥ 0 (≤ 0) on the rays and ≢ 0. ∎

**Lemma 7.2 (max-principle rule; sound; generalises agent P's Lemma CF).** Let σ be as in 7.1, with blocks unions of
source classes. Call p exclusive if all its consumers lie in one block β(p). Let S be a strongly connected set of ≥ 2
blocks of the graph B → β(p) (B produces the exclusive species p) such that every species consumed by a block of S is
exclusive and all its producers lie in S. Then S lies in one block of κ_V for every V ∈ C°.
*Proof.* For an exclusive p, row p of A_Vω = 0 gives ω_{β(p)}cons(p) = Σ_B ω_B prod_B(p), and with ω = 1,
Σ_B prod_B(p) = cons(p) > 0: ω_{β(p)} is a convex combination of the values on the producer blocks. Take B* ∈ S with
maximal ω over S; it consumes some p (every reaction has a nonzero source under (H1)), so every producer block of p has
the maximal value. The maximal set is closed under predecessors inside S, hence equals S. ∎

**Algorithm.** Lower bound `pi_plus`: flat closure, then Lemmas 7.1 and 7.2, re-closing after every merge (sound:
lower ⪯ CF). Upper bound `cf_upper`: meet of κ_V for three random V ∈ C° (positive integer combinations of the extreme
rays) and for the special V = V_U ⊕ V_{U^c} whenever a union U of lower-bound blocks and its complement are both
consistent subnetworks (then each part balances alone and U separates). CF is determined when the bounds coincide.
For a one-block upper bound, `one_block_certificate` proves CF = one block when some (|σ|−1)-minor of the block matrix
X(V) = [Γ_B V_B]_B, as a polynomial in the ray weights λ, has all coefficients of one sign.

**Results** (`scripts/check_families.py`, logs `logs/fam_*.log`, on agent E's records; every network has (H1),
(H2); both N and N*; 31 861 networks, 63 722 sides; random V with weights in 1..10⁶):

| family | networks | CF determined N / N* | P ≠ P_cs (N) | P_cs ≠ π (N) | π ≠ CF N / N* | κ_gen ≠ CF (N) | isolation N / N* | species-side Q2 (all; (H3)) |
|---|---|---|---|---|---|---|---|---|
| all (2,2,5) | 100 | all / all | 0 | 10 | 2 / 0 | 0 | 20 / 100 | 100; 100 of 100 |
| all (2,2,6) | 116 | all / all | 0 | 10 | 2 / 0 | 0 | 20 / 116 | 116; 116 of 116 |
| all (2,3,4) | 2 784 | all / all | 27 | 63 | 38 / 0 | 5 (smallest XY4) | 224 / 2 784 | 2 784; all |
| all (3,2,4) | 3 883 | all / all | 0 | 1 776 | 0 / 0 | 0 | 1 894 / 3 883 | 3 883; all |
| all (3,3,3) | 11 409 | all / all | 0 | 9 424 | 0 / 0 | 0 | 10 277 / 11 347 | 11 352; 11 329 of 11 332 |
| all (4,1,12) | 84 | all / all | 39 | 2 | 37 / 0 | 0 | 84 / 84 | 84; all |
| all (4,2,3) | 151 | all / all | 0 | 136 | 0 / 1 | 0 | 142 / 141 | 137; 106 of 116 |
| FC (3,6) / (4,5) | 8 / 76 | all / all | 1 / 7 | 3 / 12 | 1 / 0, 0 / 1 | 0 | all / all | all; all |
| WR (2,2,5) zero cplx | 16 | all / all | 0 | 4 | 0 / 0 | 0 | 4 / 16 | all; all |
| WR (3,2,5) zero cplx | 479 | all / all | 30 | 145 | 0 / 0 | 0 | 123 / 479 | all; all |
| WR (3,2,6) | 2 587 | 2 582 / all | 49 | 256 | 10 / 0 | 0 | 127 / 2 587 | all; all |
| WR (4,3,4) | 7 076 | all / all | 0 | 6 133 | 0 / 1 | 0 | 6 875 / 4 451 | 4 448; 3 743 of 6 139 |
| WR (5,2,5) | 2 352 | all / all | 9 | 1 567 | 0 / 0 | 0 | 2 247 / 1 582 | 1 412; 678 of 1 564 |
| WR (6,2,5) homogeneous | 740 | all / all | 1 | 560 | 0 / 7 | 0 | 730 / 348 | 283; 15 of 251 |

CF is determined by the bounds on 63 717 of the 63 722 sides; the 5 others (WR (3,2,6)) are weakly reversible with
δ = 0 (complete digraphs on three bimolecular complexes, e.g. 2X₃ ⇄ X₂ ⇄ X₁+X₃ ⇄ 2X₃), where CF is the single
linkage class by Theorem 4.2. Claims K1–K7 (Props 2.1, 2.2, Thm 6.1 on both sides, Thms 4.1, 4.2, Cor 5.3,
Prop 6.4) PASS on every family. Two corrections made during the work: K5 first tested the lower bound instead of the
theorem's statement (it flagged exactly those 5 sides); and with random weights in 1..30 the generic V was
sometimes special (Z₄ shows special V can give a *coarser* κ_V), which inflated the κ_gen ≠ CF counts; the table is
from the corrected rerun. isolation N* equals Condition 3.1 at generic V for N* in every family (Cor 5.3);
on the N side Condition 3.1 at generic V is slightly more frequent (229 vs 224 in all (2,3,4), 132 vs 127 in
WR (3,2,6)): XY4-type networks.

## 8. What is left

- **Exact algorithm for CF.** CF is the meet over a semialgebraic family; the two bounds settle 63 717 of 63 722
  tested sides and theory settles the rest, but no finite exact procedure is proved in general (the separating V can
  be special: Z₄, XY4). A real-algebraic decision (rank strata of A_V on C°) would close it.
- **WR with δ ≥ 1**: describe CF (when do linkage classes merge?). Candidate Q2-12 for E2.
- **Prior art.** Lemma 1.2 (clusters = proportionality classes of ker A_k) and Corollary 5.2 were not found in
  Conradi–Flockerzi 2012, the thesis, or agent L's notes. Pérez Millán et al. 2012 (arXiv:1102.1590, text in the
  git-ignored `certificates/sources/pmdsc/`) define Condition 3.1 (a basis of ker Σ with disjoint supports, §3) and
  the doubling matrix D (§5), but do not mention the Conradi–Flockerzi clusters, the isolation property, or the
  proportionality classes of the coordinates of ker Σ (searched for flockerzi, cluster, isolation, doubling,
  proportional: only reference-list hits and the doubling matrix). Proposition 5.1 is elementary linear algebra once
  stated; present it as an observation linking the two literatures. The invariance results (§3) and the duality
  theorems (§6) concern Iosif's notions and were not found anywhere (agent L §18). A scan of the AI campaigns of
  PROTOCOL.md §1 was not repeated here (agent L did it on the same day for the same key phrases).
- Lemmas 7.1/7.2 only certify merges; a complete combinatorial closure rule equal to CF is open.

## Files

- `scripts/pclust.py` — library (imports `crn.py`, `petri.py`): source classes, flat closure π, full coupling, A_V,
  κ_V, CF bounds (Lemmas 7.1, 7.2, special V), one-block certificate (sympy), operations `duplicate`,
  `split_choice`.
- `scripts/check_examples.py` — claims X1–X9 on R, T, W, U, D, F, N₃, N₄.₁₂, B₂, B₃, K₃, XY4, Z₄, Y5, B₃split,
  B₂dup, N₄.₁₂dup (`logs/check_examples.log`, all PASS, ~20 s, needs sympy for X5).
- `scripts/check_families.py`, `scripts/run_families.sh` — claims K1–K7 and reports on agent E's families
  (`logs/fam_*.log`, all PASS; ~2–3 min per family on WSL).
- `candidates.md` (for E2), `PROGRESS.md`. No `NEEDS.md`.
