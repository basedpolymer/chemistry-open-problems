# The Petri-net route to Iosif's duality: reverse-duals, clusters, siphons and absolute consistency

Agent P — **Claude Opus 5.5 (effort xhigh)**, 27 October 2026. Working notes for result 17
(`problems/iosif-siphon-precluster-duality/`), questions Q-A and Q-B of `../BRIEF-exploration.md`. Not a public
document yet: no quotation of the author's e-mail is used or needed here.

## Summary

| # | Statement | Status |
|---|---|---|
| 1.1 | Iosif's dual N* is the weighted reverse-dual rd(N) = (T, S, F⁻¹) (Pre and Post transposed separately) | proved (definitional); checked on 11 networks |
| 1.2 | Dictionary: conservation laws ↔ flux vectors, siphons ↔ T-traps, traps ↔ T-siphons, S-components ↔ T-components, clusters ↔ clusters, FC ↔ FC, EC+balanced joins ↔ itself, (H1) ↔ (H2), (H3) ↔ distinct columns, rank and #clusters preserved | proved (short proofs or Desel–Esparza Prop. 6.20); checked |
| 2.2 | The doubling classes (components of Iosif's initial precluster graph G) are the reaction parts of the Petri clusters **iff N is equal conflict (EC)**; twin classes are the species parts iff N* is EC | proved |
| 2.3 | For EC networks with balanced joins the initial preclusters of N and of N* are the two halves of one self-dual partition (the clusters) | proved |
| 2.5 | Iosif's pairwise span rule is not stable under splitting a choice: B₂ (2 branches) has one precluster, B₃ (3 branches) two; the class-span (Conradi–Flockerzi ratio) closure gives one block in both | proved, checked |
| A | **Transversality** (all networks): a minimal siphon contains at most one species of each twin class; a minimal T-trap at most one reaction of each doubling class | proved; checked on 11 748 EC networks and on the general families |
| B | **Exact criterion for Q3**: every precluster of N* lies in a minimal siphon of N iff N has no twins and every species lies in a minimal siphon (dually for N) | proved; checked |
| C | **Well-formed EC networks have no critical siphons**; minimal siphons = supports of extreme conservation laws (weighted case; pure nets) | proved from Teruel–Silva (Thm 28, Thm 15, Cor 13) + Gordan; checked |
| C′ | **Their duals have no critical siphons either**: every T-trap of a WF-EC network contains an internal cycle (no purity, no balanced joins). With C: for pure WF-EC networks, minimal siphons of N = internal cycles of N*, and minimal siphons of N* = internal cycles of N (fails in general: D) | proved (from Theorem T by a limit of rate constants); checked |
| D | **Structure of well-formed free-choice (WF-FC) networks**: N* is WF-FC; extreme conservation laws = indicators of minimal siphons (S-components); extreme flux rays = indicators of minimal T-traps; ray MIPS = atoms of minimal siphons; P(N*) is computed from twins + minimal siphons by an explicit rule; no critical siphons on either side | proved from Desel–Esparza (Prop. 5.4, Thm 6.14, 6.21); checked on 75 WF-FC networks |
| K8 | Outside WF-FC, P(N*) is not determined by the twin classes and the minimal siphons, nor P(N) by the doubling classes and minimal T-traps (explicit pairs with (H1)–(H3)) | checked (exact) |
| E | **Obstruction**: every formulation of Table 1 of the manuscript fails on a WF-FC network satisfying (H1)–(H3) (R, T, W, W*, F); in every WF-FC network with n, r ≥ 2, A1 and B1 fail, A2 ⇔ JF, B2 ⇔ CF, Q3 ⇔ JF/CF (preclusters trivial) | proved; checked |
| K5 | For EC, connected, (H1), (H2): s ≥ m − 1, m = number of source complexes | proved (elementary) |
| T | **Absolute consistency of EC networks**: EC + connected + (H1) + (H2). Then [positive steady states for every k] ⇔ s = m − 1 ⇔ Petri well-formed; and then for every k the positive steady states are {x : x^{y_c} = λ u_c(k)} with a unique positive u(k): a nonempty toric set of dimension n − s | proved (self-contained: Lemma CF + continuity); credits: the toric conclusion is a special case of Pérez Millán–Dickenstein–Shiu–Conradi 2012 (Thm 3.8, Cor. 3.12) once the kernel condition is known, and the kernel condition for ordinary free-choice nets with rational rates is Campos–Chiola–Silva 1991 (Lemma 5.1, Thms 5.3–5.4); new: weighted EC, all real k, mass-action converse, duality T*; checked symbolically (W*, F, D) and on 8 757 WF-EC networks with sampled k |
| P | **Persistence transfers within WF-EC networks**: N* (and N, if pure or balanced) is persistent for all rate constants | proved (C, C′ + Angeli–De Leenheer–Sontag Thm 2) |
| T* | **A duality with dynamical content**: for EC networks with balanced joins (e.g. ordinary free choice), N has positive steady states for all rate constants iff N* does; both steady-state sets are then toric (dimensions n − s and r − s) | proved |
| T-cex | Outside EC this fails: U has a positive steady state for every k, U* has none when κ₂ ≥ κ₁ | proved, checked (sympy) |
| F | WR + EC, no zero complex ⇒ δ = 0 (Mairesse–Nguyen 2009 for ordinary nets), minimal siphons = supports of extreme conservation laws, no critical siphon (so the GAC holds by classical results; the GAC is now claimed in full anyway); WR + ordinary FC + (H3) ⇒ every complex is a single species | proved; checked |
| H | **MIPS refine preclusters away from singletons** (all networks with (H1)): every precluster of N* with ≥ 2 elements is a union of MIPS of N; dually for N | proved; checked on the general families |
| G | "MIPS 'refine' siphons": (iii) minimal siphons = supports of extreme conservation laws implies (i) every minimal siphon is a union of MIPS and (ii) every MIPS lies in a minimal siphon; (iii) holds for all well-formed EC networks; D refutes (i)–(iii) | proved; counts in logs |
| open | Theorem C without purity (unbalanced joins); a characterisation of P(N) = cluster partition beyond the definitional one; Birch-type uniqueness for Theorem T; absolute consistency of N* for unbalanced EC N | open |

**Independent checks by agent E** (`../computation/logs/tests_final.log`, `tests_K14.log`, 27 Oct 2026, on 434 984
networks of E's families, both sides of the duality): Theorem A (K1) on 869 968 network sides, Theorem B (K2) on
223 668, Theorem D (K3a–g) on 492 WF-FC sides, Proposition H on 546 818, Theorems C and C′ on 5 582 WF-EC sides,
Proposition G on 197 776, Theorem T (K10) on 6 208 and Corollary T* (K11) on 4 982, the siphon/cycle bijection (K14)
on 5 582: all PASS. E also found that Theorem C held on all 2 116 impure WF-EC sides tested (so purity may be
unnecessary), and that K7 as first stated failed on 26 networks with the zero complex (now restricted, Prop. F).

The answer to Q-B in short: **the Petri-net reverse-dual makes Iosif's duality exact for clusters, siphons
(↔ T-traps) and invariants; preclusters and siphons are related by *transversality*, not refinement — every minimal
siphon meets every precluster seed in at most one element (Theorem A), so no refinement-type relation can hold
except where the preclusters are trivial (Theorems B, E), even in the class where Hack's duality holds; and the
duality statements with content that survive live in the well-formed equal-conflict class: there, the absence of
critical siphons and persistence hold on both sides (C, C′, P), and positive toric steady states for all rate
constants are invariant under N ↦ N* when joins are balanced (Theorem T*).**

## 0. Notation

Manuscript notation: N has species X₁..Xₙ, reactions Rⱼ: yⱼ → y′ⱼ, Ye, Yp, Γ = Yp − Ye, s = rank Γ, N* the dual
(species k₁..k_r, reaction ρᵢ: k^{(Ye)ᵢ·} → k^{(Yp)ᵢ·}). (H1) conservative: ∃ z > 0, zᵀΓ = 0. (H2) consistent:
∃ v > 0, Γv = 0. (H3) rows of Γ pairwise distinct. P(N) preclusters (Def. 19, all variants agree, manuscript
Lemma 5), M(N) ray MIPS, S(N) minimal siphons, n_t the rows of the matrix of extreme flux rays.

Petri reading. Places = species, transitions = reactions, Pre(Xᵢ, Rⱼ) = (Ye)ᵢⱼ, Post(Xᵢ, Rⱼ) = (Yp)ᵢⱼ.
•Rⱼ = supp yⱼ (reactant species), Rⱼ• = supp y′ⱼ; •Xᵢ = reactions producing Xᵢ, Xᵢ• = reactions consuming Xᵢ.

- *ordinary*: Ye, Yp ∈ {0,1}^{n×r}.
- *EC (equal conflict, Teruel–Silva Def. 6.3)*: •Rⱼ ∩ •Rₖ ≠ ∅ ⇒ yⱼ = yₖ. Equivalently: distinct source complexes
  have disjoint supports. *FC (free choice, Desel–Esparza)* = ordinary + EC.
- *balanced joins*: every source complex has all its nonzero coefficients equal.
- *JF (join-free)*: every source complex is a·Xᵢ. *CF (choice-free)*: every species is consumed by at most one
  reaction.
- *doubling class*: reactions with the same nonzero source complex (Teruel–Silva's equal-conflict set). The
  doubling pairs of the manuscript are the pairs inside a doubling class. *twin class*: species with the same nonzero
  row of Ye (= doubling classes of N*).
- *cluster* (Desel–Esparza Def. 4.4, "[x]"): connected component of the undirected bipartite graph with an edge
  {Xᵢ, Rⱼ} iff (Ye)ᵢⱼ > 0. Write c = (S_c, T_c). m = number of distinct (nonzero) source complexes.
- *connected*: the bipartite graph with an edge for each nonzero entry of Ye or Yp is connected.
- *T-trap* of N: a nonempty set U of reactions with U• ⊆ •U (every species produced by U is consumed by U).
  *T-siphon*: •U ⊆ U•.
- *rank equation*: s = (number of clusters) − 1; for EC networks with every species consumed this is s = m − 1.
- *WF-EC*: EC + connected + (H1) + (H2) + rank equation; *WF-FC*: WF-EC and ordinary. By the Rank Theorem
  (Desel–Esparza Thm 6.14 for FC, Teruel–Silva Thm 28 for EC) these are exactly the EC/FC networks whose Petri net is
  *well-formed* (admits a live and bounded marking).
- *critical siphon*: a siphon containing the support of no nonzero z ≥ 0 with zᵀΓ = 0.

Under (H2) every species that occurs is consumed by some reaction: if (Ye)ᵢ· = 0 then Γᵢ· = (Yp)ᵢ· ≥ 0 and
Γᵢ·v = 0 for some v > 0 forces (Yp)ᵢ· = 0. So under (H2) the clusters are in bijection with the source complexes
when N is EC, and m = number of clusters.

## 1. The dictionary

**Proposition 1.1.** N* = rd(N) as weighted nets: the places of N* are the reactions of N, its transitions are the
species of N, and Pre_{N*}(Rⱼ, Xᵢ) = Pre_N(Xᵢ, Rⱼ), Post_{N*}(Rⱼ, Xᵢ) = Post_N(Xᵢ, Rⱼ). In particular every arc is
reversed and its weight is kept, and catalytic self-loops are kept as self-loops.

*Proof.* The reaction ρᵢ of N* consumes kⱼ with coefficient (Ye)ᵢⱼ = Pre_N(Xᵢ, Rⱼ) and produces kⱼ with coefficient
(Yp)ᵢⱼ = Post_N(Xᵢ, Rⱼ). So the arc Xᵢ → Rⱼ of N (weight Pre_N(Xᵢ, Rⱼ)) becomes the arc kⱼ → ρᵢ of N*, from the
place kⱼ = Rⱼ to the transition ρᵢ = Xᵢ, with the same weight; and the arc Rⱼ → Xᵢ of N (weight Post_N(Xᵢ, Rⱼ))
becomes the arc ρᵢ → kⱼ. Every arc is reversed and places and transitions are exchanged: this is the reverse-dual
(Desel–Esparza Def. 6.19, "dual net"; Teruel–Silva §2, "reverse-dual"), with weights. ∎

(`scripts/check_dictionary_examples.py`, claim D1, checks Pre and Post entrywise on R, T, W, W*, U, D, F, N₃, N₄.₁₂,
B₂, B₃.)

**Proposition 1.2 (dictionary).** For every network N:

| N | N* = rd(N) | proof |
|---|---|---|
| species / reactions | reactions / species | 1.1 |
| conservation law z (zᵀΓ = 0) | flux vector (Γᵀz = 0) | Γ(N*) = Γᵀ |
| (H1) conservative | (H2) consistent | same |
| (H3): rows of Γ distinct | columns of Γ distinct | same |
| rank s | rank s | same |
| siphon Σ (•Σ ⊆ Σ•) | T-trap Σ (Σ• ⊆ •Σ in N*) | below |
| trap Σ (Σ• ⊆ •Σ) | T-siphon | below |
| critical siphon | T-trap of N* containing no internal cycle of N* | from the two lines above |
| S-component (ordinary) | T-component | Desel–Esparza Prop. 6.20(6) |
| cluster (S_c, T_c) | cluster (T_c, S_c) | below |
| twin class | doubling class | definitions |
| doubling class | twin class | definitions |
| FC | FC | Desel–Esparza Prop. 6.20(2) |
| EC with balanced joins | EC with balanced joins | below |
| EC | "twin-closed": species sharing a consumer have equal Ye rows | definitions |
| connected | connected | same underlying graph |
| MIPS(N) | fully coupled classes of N* | manuscript Prop. 9 |
| supports of extreme conservation laws | internal cycles | manuscript Lemma 6 (Iosif Thm 21) |
| WF-FC / WF-EC with balanced joins | same | rows above (Hack's Duality Theorem; Teruel–Silva Cor. 32) |

*Proofs.* Siphons: for Σ ⊆ species, the reactions ρᵢ (i ∈ Σ) of N* produce exactly the kⱼ with (Yp)ᵢⱼ > 0, i.e.
the reactions of N producing a species of Σ (•Σ in N), and consume exactly the kⱼ with (Ye)ᵢⱼ > 0, i.e. Σ• in N.
So Σ is a T-trap of N* iff •Σ ⊆ Σ•, i.e. iff Σ is a siphon of N (Definition 13 of the source: every reaction with a
product in Σ has a reactant in Σ). Traps are the same with the inclusion reversed. Clusters: the cluster graph of N*
has an edge {kⱼ, ρᵢ} iff (Ye(N*))ⱼᵢ = (Ye)ᵢⱼ > 0, the same bipartite graph with the two vertex classes exchanged.
EC with balanced joins: N is EC iff on each cluster block Ye[S_c, T_c] all columns are equal (= y_c restricted to
S_c), and has balanced joins iff each such column is constant; so "EC + balanced joins" says Ye[S_c, T_c] = a_c·J
(all-ones matrix), a condition invariant under transposition. Also outside the blocks Ye vanishes by definition of
clusters. ∎ (Checked: D2–D4.)

EC alone is *not* self-dual: N = 2X₁ + X₂ ⇄ X₃ is EC, but its dual has the reactions ρ₁: 2K₁ → 2K₂ and
ρ₂: K₁ → K₂, which share the reactant K₁ and have different sources.

## 2. Preclusters and Petri clusters

**Proposition 2.1.** The components of the graph G of Definition 19 (edges: doubling pairs) are the doubling
classes. Under (H1) no reaction has the zero source (z > 0, zᵀ(y′ − 0) = 0 forces y′ = 0), so they are exactly the
equal-conflict sets of Teruel–Silva. The components of G(N*) are the twin classes of N. (Definitions.)

**Proposition 2.2 (task 2).** Assume every reaction has a nonzero source. Then the doubling classes refine the
transition parts T_c of the clusters, with equality iff N is EC. Assume every species is consumed (e.g. (H2)).
Then the twin classes refine the species parts S_c, with equality iff N* is EC, iff any two species occurring in a
common source complex have identical rows of Ye.

*Proof.* Reactions with the same source share all reactants, so lie in one cluster. If N is EC, two reactions of a
cluster are joined by a chain of reactions consecutively sharing a reactant, and EC makes consecutive sources
equal, so the whole chain has one source. Conversely, if doubling classes = T_c, two reactions sharing a reactant
lie in one cluster, hence have the same source: EC. The species statement is the same for N* (Prop. 1.2). ∎

**Corollary 2.3 (the initial preclusters of N and N* are dual).** If N is EC with balanced joins, the maps
T ↦ •T (reactant species) and S ↦ S• (consuming reactions) are mutually inverse bijections between the doubling
classes of N and the twin classes of N, and both are the halves (T_c, S_c) of the self-dual cluster partition. So
before the span rule, "preclusters of N" and "preclusters of N*" correspond exactly. (Prop. 2.2 + Prop. 1.2.)

**Proposition 2.4 (preclusters as an enlargement of clusters).** For EC N with (H2), P(N) is obtained from the
cluster partition (transition parts) by the span rule, and P(N) = {T_c} iff for every doubling pair {i, j} and every
reaction s outside the class of i, n_s ∉ span{n_i, n_j}. In a WF-FC network this reads (Theorem D(e) below):
Θ(s) ∉ {Θ(i), Θ(j), Θ(i) ∪ Θ(j)}, where Θ(t) is the set of minimal T-traps (= T-components) containing t. In
particular P(N) = {T_c} whenever N is CF (no doubling pair). (Definitions; Theorem 7 of the manuscript gives the
coupling version.)

**Remark 2.5 (the pairwise span rule is not natural from the Petri side).** Let B_k: X → Y_a (a = 1..k),
Y_a → Z, Z → X. B_k is WF-FC and satisfies (H1)–(H3). Its T-components are {X→Y_a, Y_a→Z, Z→X}, so
n_{X→Y_a} = n_{Y_a→Z} = e_a and n_{Z→X} = e₁ + … + e_k. For k = 2, n_{Z→X} ∈ span{e₁, e₂}: P(B₂) is one block. For
k = 3, e₁ + e₂ + e₃ lies in no span{e_a, e_b}: P(B₃) = {all but Z→X}, {Z→X}. The three variants of the manuscript
agree on both. The ratio argument of Conradi–Flockerzi (Fact 2(a)) applies to whole doubling classes: if all
n_t (t in a class B) take the same ratio n_t·λ/n_t·μ, so does every n_s ∈ span{n_t : t ∈ B}. The resulting
"class-span closure" (implemented as `span_closure_preclusters`) gives one block for both B₂ and B₃. So Definition 19
depends on how a choice is split into binary choices; a Petri-invariant notion of precluster should use the span
of a whole cluster. (Checked: P1.)

## 3. What the classical theorems say for reaction networks

From Desel–Esparza, *Free Choice Petri Nets* (1995), read in `certificates/sources/fcbook/` (git-ignored):

- **Rank Theorem (Thm 6.14).** A free-choice net is well-formed iff (a) connected with a place and a transition,
  (b) a positive S-invariant, (c) a positive T-invariant, (d) rank = #clusters − 1. *CRN form*: an ordinary FC
  network is well-formed iff it is connected and satisfies (H1), (H2) and s = m − 1. Prop. 6.15: if (a)–(c) hold and
  N is not well-formed, then s ≥ #clusters.
- **Duality Theorem (Thm 6.21, Hack 1972).** N is a well-formed FC net iff N^d = rd(N) is. *CRN form*: N is WF-FC
  iff N* is. (With our definitions this is visible directly: every condition of WF-FC is self-dual by Prop. 1.2.)
- **Prop. 5.4 and 5.2(3).** In a well-formed FC net the minimal siphons are exactly the place sets of the
  S-components. **Thm 5.17.** Minimal T-invariants induce T-components. **Thm 4.30(a).** A minimal siphon of an FC net
  contains at most one place of each cluster. **Prop. 6.1(1) / Commoner.** A structurally live FC net is live from
  the marking with one token per place.

From Teruel–Silva 1996 (on disk, `certificates/sources/` of the problem):

- **Thm 28 (rank theorem for EC nets).** An EC net is well-formed iff it is well-structured (conservative and
  consistent) and rank = |𝓔| − 1 (𝓔 = equal-conflict sets), iff strongly connected and P-allocatable, iff strongly
  connected and T-allocatable. *CRN form*: a pure EC network is well-formed iff connected + (H1) + (H2) + s = m − 1.
- **Cor. 32 (duality).** EC nets are not closed under rd, but EC nets with balanced joins are, and well-formedness
  is preserved. *CRN form*: N EC with balanced joins and WF ⇔ N* EC with balanced joins and WF.
- **Thm 12, Cor. 13** (potential reachability of live EC systems: every marking satisfying the state equation is
  live), **Thm 15** (liveness monotonicity for bounded EC systems), **Thm 7** (strongly connected CF nets: a
  T-semiflow implies consistency, and then the minimal T-semiflow is unique).
- The authors stress (pp. 272, 296) that **siphons and traps have no satisfactory theory for weighted nets**, and
  Commoner's theorem is not generalised. Theorem C below shows that the part needed here (no critical siphon) does
  extend.

Jančar–Best–Devillers–Ošťádal (arXiv:2602.02352v2, in `certificates/sources/jbdo/`) restate Hack's duality as
Theorem 3.23 (a)⇔(d) (the exploration brief's "Thm 3.27" is a different proposition), and their §4 names the
reverse-dual of minimal siphons "minimal T-traps" (T_Y• ⊆ •T_Y) for ordinary free-choice nets; nothing beyond the
above is used here. Hack's 1972 report (MIT Project MAC TR-94, open on DSpace; read by agent L) contains the
original Well-Formedness Theorem (Thm 9, (a)⇔(c): an FC net is well-formed iff its reverse-dual is) and remarks
that he had no interpretation of the reverse-dual of a trap. Note also (agent L): the Rank Theorem 6.14 concerns
nets; the condition "every proper siphon is marked" belongs to Desel–Esparza Thm 6.17 (live and bounded systems).

**Scope (task 3).** For reaction networks the classical theory therefore covers exactly the EC networks (distinct
source complexes with disjoint supports): the Rank Theorem holds for all pure EC networks, the Duality Theorem for
EC networks with balanced joins. It does not extend to general networks: a network with the sources X and X + Y is
not EC, and Desel–Esparza Exercise 6.6 asks for a well-formed non-free-choice net whose dual is not well-formed. In
the GAC setting the EC class is small: a weakly reversible EC network has pairwise disjoint complexes (Prop. F).

Token-game confirmation of the CRN form of the Rank Theorem: on every EC network of the enumerated families with
at most 4 species (11 748 networks), the rank equation holds iff the Petri net is live from the Teruel–Silva marking
M₀(p) = Pre(p, t) (t ∈ p•) — claim "RT" in `logs/enum_ec_*.log`.

## 4. Transversality: why preclusters and siphons cannot be related by refinement

**Theorem A (transversality; no hypothesis).** Let R be a minimal siphon of N and let i ≠ j be species with
i• = j• (the same set of consuming reactions; in particular twins). Then not both i and j lie in R. Dually, a
minimal T-trap U of N (= minimal siphon of N*) contains no two reactions with the same reactant set (in particular
no doubling pair).

*Proof.* Suppose i, j ∈ R. We show that R′ = R \ {i} is a siphon, contradicting minimality (R′ ∋ j is nonempty).
Let a reaction t produce a species of R′ ⊆ R. Since R is a siphon, t consumes a species of R. If that species is i,
then t ∈ i• = j•, so t consumes j ∈ R′. Otherwise it is in R′. The dual statement is the same proof for U \ {t}:
the species produced by U \ {t} are produced by U, hence consumed by U, hence consumed by U \ {t} because the
reactant set of t equals that of the remaining t′ ∈ U. ∎

For FC nets this is the first half of Desel–Esparza Thm 4.30. In an EC network all species of a cluster have the
same consumers T_c and all reactions of a cluster have the same reactant set, so: **in an EC network every minimal
siphon meets every cluster in at most one species and every minimal T-trap meets every cluster in at most one
reaction.**

**Theorem B (exact criterion for Q3; no hypothesis).** Every precluster of N* lies in a minimal siphon of N iff N has
no twin pair and every species of N lies in some minimal siphon. Every precluster of N lies in a minimal siphon of
N* iff N has no doubling pair and every reaction lies in some minimal T-trap.

*Proof.* Every edge of the precluster graph G⁺(N*) is a twin pair or an edge {i, s}, {j, s} adjoined through a twin
pair {i, j}, so a precluster with ≥ 2 elements contains a twin pair, and then lies in no minimal siphon by Theorem A.
If there is no twin pair, G(N*) and hence G⁺(N*) have no edge, P(N*) is discrete, and the condition is that every
species lies in a minimal siphon. The second statement is the first for N*. ∎

This explains every failure of Q3 in the manuscript (D, F) and shows that Q3 can hold only when the preclusters are
trivial. Checked against direct computation on all 11 748 enumerated EC networks (claim K2) and on the general families of
the manuscript (`logs/general_*.log`, claim B); E tests it on larger families (candidate K2).

*Interpretation.* The seeds of the preclusters (twin classes of N, doubling classes of N) are the cluster halves;
minimal siphons and T-traps (and, in well-formed nets, all extreme conservation laws and flux modes) cross every
cluster at most once. The relation between preclusters and siphons is therefore one of **orthogonality** — like
the relation between the blocks of a partition and its transversals — and never one of inclusion.

## 5. Siphons of well-formed networks

**Theorem C (well-formed EC networks have no critical siphons).** Let N be a pure WF-EC network (no species on both
sides of a reaction; EC, connected, (H1), (H2), s = m − 1). Then every siphon of N contains the support of a nonzero
semi-positive conservation law. Consequently the minimal siphons of N are exactly the supports of the extreme
conservation laws (the semiflow-reading MIPS).

*Proof.* By Teruel–Silva Thm 28, N is well-formed: some marking M₀ is live, and every marking is bounded (H1).
Let R be a siphon and suppose it contains no such support. By Gordan's alternative applied to the |R| × r matrix
Γ[R, ·] (either y ≥ 0, y ≠ 0 with yᵀΓ[R,·] = 0, which extended by 0 is a conservation law supported in R; or some
x with Γ[R,·]x < 0 componentwise), there is x ∈ ℚʳ with Γ[R,·]x < 0. Adding a large multiple of an integer v > 0 in
ker Γ (H2) and clearing denominators gives σ ∈ ℤʳ_{≥0} with Γ[R,·]σ < 0. Choose an integer λ ≥ 1 and a marking
M₀′ ≥ M₀ with M₀′[R] = −λΓ[R,·]σ (≥ M₀[R] for λ large) and M₀′ large outside R, so that M := M₀′ + λΓσ ≥ 0 and
M[R] = 0. By Thm 15 (liveness monotonicity for bounded EC systems) (N, M₀′) is live. M satisfies the state equation
from M₀′, so by Cor. 13 (N, M) is live. But R is an empty siphon at M: a reaction putting tokens in R must take a
token from R, so R stays empty and the reactions of R• (nonempty: N is strongly connected, Lemma SC below) are dead.
Contradiction.
For the consequence: a minimal siphon R contains the support Z of a semi-positive conservation law; Z is a siphon
(a reaction producing i ∈ Z has z·y′ > 0, hence z·y > 0); so Z = R, and R is the support of an extreme law (a law
with minimal support; any law with support ⊆ R is again a siphon). Conversely the support Z of an extreme law is a
siphon, contains a minimal siphon R, which is the support of a law; minimality of Z gives R = Z. ∎

Status: proved from Teruel–Silva (Thm 28 for well-formedness, Thm 15, Cor. 13), which assume pure nets. Checked on
all 8 757 WF-EC networks of the enumerated families (claim "K4 min siphons = supports", class WF, 100 %), which also
contain impure ones; and in the same families no WF-EC network had a dual with a critical siphon, which is
Theorem C′ below.

**Lemma SC.** A connected network with (H1) and (H2) is strongly connected as a Petri net.
*Proof.* Otherwise take a terminal strong component K (no arc leaves K), K ≠ everything; let P_K, T_K be its species
and reactions. With z > 0, v > 0 as in (H1), (H2): 0 = Σ_{p∈P_K} z_p (Γv)_p
= Σ_{t∈T_K} v_t Σ_{p∈P_K} z_p(Post − Pre)(p,t) + Σ_{t∉T_K} v_t Σ_{p∈P_K} z_p Post(p,t) (reactions outside T_K do
not consume from P_K). For t ∈ T_K all products lie in P_K, so the inner sum is ≥ Σ_p z_p (Post − Pre)(p,t) = 0.
Both sums are ≥ 0, hence 0: no reaction outside K produces into P_K, and no species outside K is a reactant of a
reaction in K. So no arc enters K either, and K is a union of components: contradiction with connectedness. ∎

**Theorem D (structure of WF-FC networks).** Let N be WF-FC (ordinary, FC, connected, (H1), (H2), s = #clusters − 1).
Then:
- (a) N* is WF-FC.
- (b) The extreme rays of the conservation cone are exactly the 0/1 indicator vectors of the minimal siphons, which
  are the place sets of the S-components.
- (c) The extreme flux rays are exactly the indicator vectors of the minimal T-traps, which are the transition sets
  of the T-components (the internal cycles).
- (d) With Σ(i) the set of minimal siphons containing species i: i, j are in one ray MIPS iff Σ(i) = Σ(j). So M(N)
  is the partition into atoms of the family of minimal siphons, every minimal siphon is a union of MIPS, every MIPS
  lies in a minimal siphon. Dually, the fully coupled classes M(N*) are the atoms of the internal cycles.
- (e) (*profile rule*) Twins i ≠ j have disjoint nonempty Σ(i), Σ(j), and P(N*) is the component partition of the
  graph with edges {i, j} (twins) and {i, s}, {j, s} whenever Σ(s) ∈ {Σ(i), Σ(j), Σ(i) ∪ Σ(j)}. Dually P(N) is
  given by the same rule with doubling pairs and Θ(t) = set of minimal T-traps containing t.
- (f) Neither N nor N* has a critical siphon.
- (g) No twin pair is conservation-coupled and no doubling pair is fully coupled.
- (h) N is live and bounded from the marking (1, …, 1).

*Proof.* (a) Each condition is self-dual (Prop. 1.2); this is Hack's Duality Theorem. (b) By the Rank Theorem N is
well-formed; Desel–Esparza Prop. 5.4 and 5.2(3) give minimal siphons = S-component place sets. For an S-component
R, every reaction touching R has exactly one reactant and one product in R, with coefficient 1, so 1_R is a
conservation law; its support is a minimal siphon, so (as in the proof of Theorem C) it is extreme, and an extreme
ray is determined by its support up to scaling. Conversely an extreme law has a minimal support, which is a siphon
containing a minimal siphon R; 1_R has support ⊆ it, so they are equal. (c) Apply (b) to N* and use Prop. 1.2.
(d) The column w_i of the matrix of extreme laws is the indicator vector of Σ(i), nonzero by (H1); two nonzero 0/1
vectors are positively proportional iff equal. (e) Twins have the same consumers, so Σ(i) ∩ Σ(j) = ∅ by Theorem A;
nonempty by (H1). If w_s = a w_i + b w_j, then w_s equals a on Σ(i), b on Σ(j) and 0 elsewhere; being a nonzero
0/1 vector, Σ(s) ∈ {Σ(i), Σ(j), Σ(i) ∪ Σ(j)}. The converse is clear. Apply to G⁺(N*); the dual statement follows from
(a). (f) Every siphon contains a minimal one, which supports a law by (b); same for N* by (a). (g) From (e) and its
dual. (h) Desel–Esparza Thm 6.17 (or Prop. 6.1): the marking with one token per place marks every proper siphon;
bounded by (H1). ∎

Checked on all 75 WF-FC networks of the families (claims K3a–K3h): `logs/enum_ec_n3_m3_r6_ord1.log`,
`logs/enum_ec_n4_m2_r5_ord1.log`, and the WF-FC members of the weighted families.

So, in the Petri class, **the preclusters of N* are a function of the clusters and the minimal siphons of N**
(rule (e)), and the preclusters of N are a function of the clusters and the minimal T-traps; the two statements are
exchanged by Hack's duality. This is the most precise positive form of "a duality relation between preclusters
and siphons" that we found. Outside WF-FC the extreme conservation laws are no longer determined by the siphons,
and neither are the preclusters:

**Obstruction K8 (outside WF-FC the preclusters of N* are not determined by the twins and the siphons).**
N₁: X₃ → X₄, X₂+X₃ → X₁+X₂, X₁+X₄ → 2X₃ and N₂: 2X₃ → 2X₄, X₂+X₃ → 2X₂, X₂+X₃ → 2X₁, X₁+X₄ → 2X₃ both satisfy
(H1)–(H3), have the same twin classes ({X₁,X₄}), the same species parts of clusters and the same minimal siphons
({X₁,X₃}, {X₂}, {X₃,X₄}), but P(N₁*) = {X₁,X₃,X₄}, {X₂} and P(N₂*) = {X₁,X₂,X₃,X₄}. (`scripts/check_k8.py`, all PASS;
found by `scripts/k8_search.py`.) Caveat: X₂ is a pure catalyst in N₁ (zero row of Γ, allowed by (H3) as
stated); neither network is EC. A pair without zero rows, on the reaction side: M₁: X₃ → X₁+X₂, X₂+X₃ → X₁+X₃,
X₁+X₂ → X₃, X₁+X₂ → 2X₂ and M₂: 2X₃ → X₂+X₃, X₁+X₃ → X₁+X₂, X₁+X₂ → 2X₃, X₁+X₂ → 2X₁ both satisfy (H1)–(H3),
have the same doubling class {R₃, R₄} and the same minimal T-traps {R₄}, {R₁,R₃}, {R₂,R₃}, but P(M₁) is one block
and P(M₂) = {R₁}, {R₂}, {R₃,R₄} (found by `check_general.py`, family n = 3, molecularity ≤ 2, r ≤ 4, where 14 such
groups occur; one species-side group occurs in the family n = 3, molecularity ≤ 3, r ≤ 3). So the profile rule of
Theorem D(e) has no analogue for general networks.

## 6. The obstruction: the Petri class does not rescue the precise formulations

**Theorem E.** (i) Every formulation of Table 1 of the manuscript (A1–A4, B1–B4, C, Q1–Q6 on both sides, R23) fails
on one of R, T, W, W*, F, which are WF-FC networks satisfying (H1)–(H3) whose duals are WF-FC (claims O1, O2:
A1, A4, B1, B4, Q1, Q6 fail on all five; A2 on W*, F; A3 on R, T, W; B2 on W, F; B3 on R, T, W*; C on R, T, W*, F;
Q2, Q4, Q5 on R, T and W or W*; Q3 on W*/W and F; R23 on R, T, W*).
(ii) Let N be WF-FC with n ≥ 2 species and r ≥ 2 reactions. Then P(N*) ≠ M(N) and P(N) ≠ M(N*);
P(N*) ⪯ M(N) iff N is JF; P(N) ⪯ M(N*) iff N is CF; Q3 holds on the species side iff N is JF and on the reaction side
iff N is CF.

*Proof of (ii).* By Theorem D(g) and the manuscript's Theorem 7(1)/Cor. 8(1), P(N) ⪯ M(N*) iff there is no doubling
pair, i.e. N is CF (for an ordinary FC network a doubling pair exists iff some species is consumed twice). If N is
CF, then (Lemma SC) N is strongly connected, and every species has exactly one consumer c(p). For a semiflow
v ≥ 0 with v_t > 0 and p ∈ t•, the balance at p gives Pre(p, c(p)) v_{c(p)} ≥ Post(p, t) v_t > 0; so supp v is closed
under successors and equals all reactions. Hence every nonzero semiflow is positive, the flux cone has a single
extreme ray, all reactions are fully coupled, and M(N*) is one block of size r ≥ 2 while P(N) is discrete: B1
fails. If N is not CF, B2 and hence B1 fail. Applying this to N* (WF-FC by D(a)) gives the species side; "CF for N*"
is "JF for N". For Q3: Theorem B, and every species (reaction) lies in a minimal siphon (T-trap) by D(b), D(c) and
(H1), (H2). ∎

Compare with Cor. 10 of the manuscript (B1 fails for WR deficiency-zero networks): here B1 fails for every WF-FC
network, a class that contains non-WR networks (W*, F, B_k) and is closed under duality.

**The manuscript's counterexamples in Petri terms** (exact; "crit" = number of critical siphons; "abs. cons." =
positive steady state for every rate vector, decided by Theorem T for EC networks and by hand for U, N₃):

| network | N: EC / balanced / ordinary / WF / crit / abs. cons. | N*: EC / balanced / ordinary / WF / crit / abs. cons. |
|---|---|---|
| R, T, W, W*, F | y / y / y / y / 0 / yes | y / y / y / y / 0 / yes |
| D (D* ≅ D) | y / y / n / n / 2 / no (only k₁ = k₂) | same |
| U | n / y / n / n / 3 / yes | n / y / n / n / 1 / no (κ₂ ≥ κ₁) |
| N₃ | n / y / n / n / 1 / no (needs k₂ > k₃) | n / n / n / n / 3 / not decided |

## 7. Absolute consistency: a duality with dynamical content

For an EC network with source complexes y_c (c = 1..m, pairwise disjoint supports) and rate constants k > 0, the
mass-action vector field is f(x) = A_k φ(x) with φ_c(x) = x^{y_c} and
A_k[·, c] = Σ_{t∈T_c} k_t (y′_t − y_c), i.e. A_k = Γ R_k with (R_k u)_t = k_t u_{c(t)}. A point x > 0 is a steady
state iff φ(x) ∈ ker A_k. Call N *absolutely consistent* if it has a positive steady state for every k > 0.

**Lemma CF.** Let N be EC, every species consumed, and strongly connected. If u ≥ 0, u ≠ 0 and A_k u = 0, then
u > 0 and ker A_k = ℝu.
*Proof.* For a species p in the source support of cluster c = c(p), row p of A_k u = 0 reads
K_c (y_c)_p u_c = Σ_{c′} (Σ_{t∈T_{c′}} k_t (y′_t)_p) u_{c′} (K_c = Σ_{t∈T_c} k_t; the y_{c′} with c′ ≠ c vanish at p by
disjointness). If u_{c′} > 0 and p is a product of a reaction of c′, the right side is > 0, so u_{c(p)} > 0. Thus
supp u is closed under the successor relation c′ → p → c(p) of the network obtained by merging each cluster into one
transition, which is strongly connected because N is. So u > 0. If w ∈ ker A_k is not proportional to u, then
u − εw for the largest ε keeping it ≥ 0 is a nonzero nonnegative kernel vector with a zero entry: contradiction. ∎

**Theorem K5.** If N is EC, connected, (H1), (H2), then s ≥ m − 1.
*Proof.* Lemma SC gives strong connectivity. Take v > 0 in ker Γ and k = v: then u = (1, …, 1) is in ker A_k, so by
Lemma CF rank A_k = m − 1, and rank A_k ≤ rank Γ = s. ∎
(Desel–Esparza Prop. 6.15 is the FC case; this proof is elementary. Checked: all 11 748 EC networks of the families.)

**Theorem T.** Let N be EC, connected, conservative (H1) and consistent (H2), with m distinct source complexes.
(1) If s = m − 1, then for every k > 0, ker A_k = ℝ u(k) with u(k) > 0, and the positive steady states are
E_k = {x > 0 : x^{y_c} = λ u_c(k) for all c, some λ > 0},
a nonempty set whose logarithm is an affine subspace of dimension n − s (a toric set).
(2) Conversely, if N is absolutely consistent then s = m − 1.
So, for such networks: absolutely consistent ⇔ s = m − 1 ⇔ Petri well-formed (Teruel–Silva Thm 28).

*Proof.* (1) rank A_k ≤ s = m − 1 < m, so ker A_k ≠ 0 for every k. Let 𝒢 be the set of k > 0 for which ker A_k
contains a nonzero u ≥ 0. 𝒢 ∋ v (with u = 1). 𝒢 is closed in ℝʳ_{>0}: normalise u to the simplex and pass to a
limit. 𝒢 is open: at k ∈ 𝒢, Lemma CF gives ker A_k = ℝu, u > 0, so rank A_k = m − 1; nearby, rank A_{k′} ≥ m − 1
(semicontinuity) and ≤ m − 1, and the kernel line (given by the maximal minors of a fixed nonsingular
(m − 1)-row submatrix) moves continuously, so it stays spanned by a positive vector. ℝʳ_{>0} is connected, so
𝒢 = ℝʳ_{>0}, and Lemma CF gives ker A_k = ℝ u(k), u(k) > 0. Since φ(x) > 0, x is a steady state iff φ(x) = λu(k).
Taking logarithms, (y_c − y_{c₀})·log x = log(u_c/u_{c₀}) for c ≠ c₀; the y_c have disjoint nonzero supports, so
they are linearly independent, these m − 1 equations are independent and solvable, and the solution set is an
affine space of dimension n − m + 1 = n − s.
(2) If s ≥ m, the k admitting a positive steady state are k_t = v_t/u_{c(t)} with v ∈ ker Γ ∩ ℝʳ_{>0} and
u ∈ ℝᵐ_{>0}, a semialgebraic image of a set of dimension (r − s) + m − 1 < r (one dimension is lost to the common
scaling of v and u); it has measure zero, so some k has no positive steady state. With Theorem K5, s = m − 1. ∎

Checked: symbolically for W*, F (kernel spanned by positive rational functions of k) and D (not well-formed: steady
states only for k₁ = k₂), `logs/check_absolute_consistency.log`; and on all 8 757 WF-EC networks of the six families
with k = 1 and three seeded random k (`logs/toric_*.log`): kernel one-dimensional and positive in all cases; for
each of the 2 991 non-WF EC networks the property failed for some sampled k, as (2) predicts.

Remarks. (i) No (H3) and no weak reversibility are needed. W*, F and B_k are absolutely consistent and toric but not
complex balanced. (ii) **Prior art (credit).** The toric conclusion of (1) is a special case of Pérez Millán,
Dickenstein, Shiu, Conradi, *Chemical reaction systems with toric steady states*, Bull. Math. Biol. 74 (2012)
1027–1065 (arXiv:1102.1590), Theorem 3.8 / Corollary 3.12: their matrix Σ (complex-to-species rate matrix) has
columns A_k[·, c] for the source complexes and zero columns for the other complexes, so ker Σ = ker A_k ⊕ ℝ^{non-sources};
once ker A_k = ℝu(k) with u(k) > 0, their Condition 3.1 holds with the partition {sources}, {each non-source
complex} (for every k, as in their Cor. 3.12), the sign condition 3.4 holds, and the log-linear system 3.6 is solvable
because the sources are linearly independent. What the Petri side adds (not found in the literature by agent L,
§16, §18): the rank equation s = m − 1 guarantees ker A_k = ℝu(k), u(k) > 0, *for every k at once* (Lemma CF plus
the continuity argument), the converse (2), and the self-duality of Corollary T* below. Equivalence with Petri
well-formedness uses Teruel–Silva Thm 28, stated for pure nets; (1) and (2) themselves need no purity. (iii) E_k meets every stoichiometric class in exactly one point when
span{y_c − y_{c′}} equals the stoichiometric subspace (Birch), e.g. when every product complex is a source complex;
for F the two spaces differ (F is known to have toric steady states with one per class by other means). The
general uniqueness question is left open (a sign-vector criterion of Müller–Regensburger type would be the tool;
not attempted). (iv) **Prior art (credit), checked.** J. Campos, G. Chiola, M. Silva, *Properties and performance
bounds for closed free choice synchronized monoclass queueing networks*, IEEE Trans. Automat. Control 36 (1991)
1368–1382 (author's copy https://webdiis.unizar.es/~jcampos/web/pdf/ccs-tac-91.pdf, in
`certificates/sources/ccs/`, git-ignored), Lemma 5.1, Theorem 5.3 and the proof of Theorem 5.4: for a live and
structurally bounded *ordinary free-choice* net and any rational conflict-resolution rates, the vector of relative
firing frequencies (C F = 0 with the routing ratios imposed at every conflict) exists, is positive and unique, and
the rank condition rank C = |T| − 1 − (a − n) (equivalent to the rank equation) is what makes this solvable "for all
conflict resolution policies". This is the flux part of Theorem T(1) for ordinary FC nets with rational rates, and
the Rank Theorem originates there. Theorem T adds: weighted EC networks, real rate vectors, the mass-action reading
(toric steady states of dimension n − s), the converse (2) stated for mass action, and the duality T*. Lemma CF is a
weighted, rate-dependent version of Teruel–Silva Thm 7 (strongly connected CF nets have a unique minimal
T-semiflow).

**Corollary T\* (self-duality).** Let N be EC with balanced joins (for instance ordinary free choice), connected,
conservative and consistent. Then N is absolutely consistent iff N* is absolutely consistent (with dual rate
constants κ ∈ ℝⁿ_{>0}) iff s = #clusters − 1. In that case the positive steady states of N and of N* are toric, of
dimensions n − s and r − s.
*Proof.* N* is EC with balanced joins, connected, conservative (N consistent) and consistent (N conservative), has
rank s and the same number of clusters (Prop. 1.2); apply Theorem T to N and N*. ∎

**The hypothesis cannot be dropped (U).** U: X₄ → X₂+X₃, X₃+X₄ → X₁, X₁+X₂ → 2X₄ satisfies (H1)–(H3) and has the
positive steady state (k₁/k₃, 1, k₁/k₂, 1) for every k. Its dual U*: K₃ → K₂, K₃ → K₁, K₂ → K₁, K₁+K₂ → 2K₃ is
conservative and consistent, but f_{K₁} − f_{K₂} = (κ₂ − κ₁)K₃ + 2κ₃K₂ > 0 on the positive orthant when κ₂ ≥ κ₁,
so U* has no positive steady state for those κ. Neither U nor U* is EC. (Checked: U1–U3.)

What T* transfers and what it does not: the duality transfers *existence and toric structure of positive steady
states for all rate constants* within balanced EC networks; it does not transfer complex balance (W is complex
balanced for all k, W* for none, manuscript §5); persistence is handled separately (Corollary P, §9).

**Theorem C′ (well-formed EC networks have no critical T-traps; the dual of Theorem C).** Let N be EC, connected,
(H1), (H2), s = m − 1 (no purity, no balanced joins). Then every T-trap U of N contains the support of a nonzero
nonnegative flux vector. Equivalently, N* has no critical siphon. Consequently the minimal T-traps of N are exactly
its internal cycles (the supports of a nonnegative flux vector are T-traps, by the balance at each produced
species), i.e. the minimal siphons of N* are exactly the supports of the extreme conservation laws of N*.

*Proof.* Shrinking U we may assume U is a minimal T-trap; by Theorem A it contains at most one reaction t_c of each
cluster c ∈ C_U (the clusters it meets), and its products lie in R_U := ∪_{c∈C_U} S_c. For ε > 0 let k(ε) be 1 on
U and on every reaction of a cluster outside C_U, and ε on the other reactions of the clusters of C_U. By
Theorem T there is u(ε) > 0 with A_{k(ε)} u(ε) = 0; normalise Σ_{c∈C_U} u_c = 1. Multiply the balance equation of
each species p ∈ R_U by z_p (z > 0 a conservation law) and sum. The consumption side is Σ_{c∈C_U} K_c u_c (z·y_c)
with K_c = 1 + ε(|T_c| − 1). On the production side, a reaction t_c ∈ U contributes u_c z·y′_{t_c} = u_c z·y_c (its
products lie in R_U and z·y′ = z·y). Hence the remaining production into R_U — by the ε-reactions and by the
reactions of outside clusters c′ — equals ε Σ_{c∈C_U} (|T_c| − 1) u_c (z·y_c) ≤ C ε. So every outside cluster that
produces into R_U has u_{c′} ≤ C′ε. If an outside cluster c′ produces a species of an outside cluster c″, the balance
at that species gives u_{c′} ≤ C″ u_{c″}. Every outside cluster has a path to R_U (Lemma SC), so u_{c′} = O(ε) for
all of them. The flux v(ε) = R_{k(ε)} u(ε) ∈ ker Γ is bounded; a limit point v* as ε → 0 satisfies Γv* = 0, v* ≥ 0,
Σ_{c∈C_U} v*_{t_c} = 1, and v* vanishes off U (the ε-reactions and the outside clusters carry O(ε) flux). ∎

Checked: claim "K6 no critical siphon in N*", 100 % of the WF networks of every family. Together with Theorem C:
**a (pure) well-formed EC network and its dual both have no critical siphons.** In the language of the problem:
for a pure WF-EC network, **the minimal siphons of N are exactly the internal cycles of N*** (Theorem C with Iosif's
Theorem 21), and the minimal siphons of N* are exactly the internal cycles of N (Theorem C′) — a bijective
siphon/cycle duality, which fails in general (D: minimal siphons {X₁,X₂}, {X₁,X₃}; the only internal cycle of D* is
{X₁,X₂,X₃}).
If N has balanced joins, N* is again WF-EC (Prop. 1.2), and Theorem C′ applied to N* gives Theorem C for N without
the purity assumption.

## 8. Loose refinement (task 5)

**Proposition H (MIPS refine preclusters away from singletons; no hypothesis beyond (H1)).** Every precluster of N*
with at least two elements is a union of ray-MIPS blocks of N. Under (H2), every precluster of N with at least two
elements is a union of fully coupled classes (= MIPS of N*). Equivalently: M(N) restricted to the union of the
non-singleton preclusters of N* refines P(N*).

*Proof.* Let B be a precluster of N* with |B| ≥ 2, s ∈ B, and s′ with w_{s′} = μ w_s, μ > 0. Since |B| ≥ 2, s is in
a twin pair {s, j}, or was adjoined to a twin pair {i, j} because w_s ∈ span{w_i, w_j}; in the first case put i = s.
In both cases w_{s′} ∈ span{w_i, w_j}, so s′ ∈ {i, j} or s′ is adjoined to {i, j} by Definition 19; either way
s′ ∈ B. The second statement is the first for N*. ∎

So the reverse refinement A3 (M(N) ⪯ P(N*)), which the manuscript refutes with R and T, fails only through MIPS
blocks all of whose species are singleton preclusters (no twin of theirs and no span adjunction); Cor. 8(2) of the
manuscript is the exact criterion. Combined with Theorems A and B: **where P(N*) is nontrivial it is coarser than
M(N) and transversal to the minimal siphons.** This is the precise content we can give to Iosif's diagram on the
species side. Checked on 18 509 networks of the general families (claim H, `logs/general_*.log`).

### "MIPS 'refine' siphons"

Readings, for a conservative N: (i) every minimal siphon is a union of (ray) MIPS; (ii) every MIPS lies in a
minimal siphon; (iii) the minimal siphons are exactly the supports of the extreme conservation laws (semiflow
MIPS); (iv) every support of an extreme conservation law is a union of minimal siphons.

**Proposition G.** Under (H1), (iii) ⇒ (i) and (ii). (iii) holds for every pure WF-EC network (Theorem C), and in
WF-FC networks moreover M(N) is exactly the partition into atoms of the minimal siphons (Theorem D(d)). D refutes
(i), (ii) and (iii): M(D) is one block, S(D) = {{X₁,X₂}, {X₁,X₃}}.
*Proof.* If i lies in the support of an extreme law z and j ~ i (w_j = μw_i), then z_j = μz_i > 0: supports of
extreme laws are unions of MIPS, which gives (i) under (iii). Under (H1) every i lies in the support of some extreme
law, which contains the MIPS of i, which gives (ii). ∎

Counts on the EC families (`logs/enum_ec_*.log`, "K9"): (i)–(iii) hold for 100 % of WF networks and fail for
many non-WF ones. (iv) held on every EC network tested but fails in general: N₄.₁₂ = 2X₁ ⇄ X₁+X₂ (not EC) has
the conservation law X₁ + X₂ but the only minimal siphon {X₁} (`logs/general_*.log`, report K9iv).
In the weighted WF-EC case M(N) can be strictly finer than the atoms of the minimal siphons (extreme laws with the
same supports need not be proportional).

## 9. Consequences for the GAC (what the Petri route can and cannot give)

**Proposition F.** Let N be weakly reversible and EC, without the zero complex (automatic under (H1)). Then all
complexes have pairwise disjoint supports, δ = 0, each linkage class has exactly the minimal siphons "one species
from each of its complexes", these are the supports of the extreme conservation laws, and N has no critical siphon.
If moreover N is ordinary and (H3) holds, every complex is a single species (N is monomolecular). (The hypothesis on
the zero complex is needed: agent E reports 26 non-conservative WR networks with the zero complex, the smallest
0 → X₂ → X₁ → 0, for which the siphon part and the last clause fail as stated.)
*Proof.* In a WR network every complex is a source, so EC makes all complexes support-disjoint, hence linearly
independent, and s = m − ℓ: δ = 0. If a siphon Σ meets a complex c, it meets every complex c′ with c′ → c, hence by
weak reversibility every complex of the linkage class L of c; choosing one species i_c ∈ Σ ∩ supp y_c per c ∈ L,
z = Σ_c e_{i_c}/(y_c)_{i_c} satisfies z·(y′ − y) = 1 − 1 = 0 on L and vanishes on the other classes. Under (H3) and
ordinary, two species of one complex would be twins with equal Yp rows (every product is a source), contradicting
(H3) (manuscript Cor. 9). ∎ (Checked: K7, all 147 WR EC networks of the families.)

**Credit.** "WR + free choice ⇒ complexes pairwise disjoint, δ = 0" is Mairesse–Nguyen 2009 (Lemma 3.2,
Cor. 3.2; ordinary nets); the weighted case is the same argument. The siphon part is ours (not found by agent L).
Agent L also notes that a WR EC network is automatically pure and conservative, so by Teruel–Silva Thm 28 it is
well-formed iff it has one linkage class (the hypothesis of Anderson's 2011 GAC case).

By Angeli–De Leenheer–Sontag (2007) a conservative network without critical siphons is persistent, and persistent
complex-balanced systems are globally stable; so **the GAC holds for every WR EC network by classical results**, and
the equal-conflict class (where the Petri duality theory works) contained no open case of the GAC even before the
GAC was claimed in full (Craciun, arXiv:1501.02860 v3, 23 Sep 2026; OpenAI, boundedness and persistence of WR
mass-action systems, 25 Sep 2026 — agent L, §13). The relevant question is therefore which structure the duality
transfers. Theorems C, C′ and D(f) show that WF-EC networks and their duals have no critical siphons, which gives:

**Corollary P (persistence transfers within well-formed EC networks).** Let N be WF-EC (EC, connected, (H1), (H2),
s = m − 1). Then N* is persistent (no trajectory starting in ℝʳ_{>0} has an ω-limit point on the boundary) for every
choice of rate constants, under mass action and under the general kinetics of Angeli–De Leenheer–Sontag. If N is
pure or has balanced joins, the same holds for N.
*Proof.* N* is conservative because N is consistent, and has no critical siphon by Theorem C′; Angeli–De
Leenheer–Sontag, Math. Biosci. 210 (2007), Theorem 2, p. 607 (§5; journal copy in the problem's
`certificates/sources/`): if the species–reaction Petri net is conservative and each siphon contains the support of a
P-semiflow, the network is persistent, for every kinetics satisfying their §2 assumptions (real-analytic rates
depending only on the reactants, nondecreasing in them, zero exactly when a reactant is absent), which mass action
satisfies. Under (H1) for N every reaction of N* has a nonzero source, so N* is a network of that kind.
For N use Theorem C (pure) or Theorem C′ applied to N* (balanced joins, N* WF-EC). ∎
Agent G's examples show that the EC hypothesis matters: 3A ⇄ 2A+B is persistent while its dual is not, and
2A ⇄ 2B, 2A ⇄ A+B satisfies the Angeli–De Leenheer–Sontag criterion while its dual has a critical siphon
(`../gac/notes.md`). Neither network is EC (2A and 2A+B, resp. 2A and A+B, share A).

**What the duality transfers in the Petri class.** Well-formedness (Hack; Teruel–Silva for balanced joins), the
absence of critical siphons and hence persistence (C, C′, P), and — for balanced joins — absolute consistency with
toric steady states (T*). It does not transfer weak reversibility or complex balance (W, W*), and outside EC it
transfers neither persistence (agent G) nor absolute consistency (U).

## 10. What is left

- Theorem C for impure networks (catalysts): the proof uses Teruel–Silva, who assume pure nets; the computations
  include impure networks and found no exception. (Theorem C′, the dual statement, needs no purity.)
- Theorem T: uniqueness per stoichiometric class (Birch-type) for WF-EC networks;
  whether absolute consistency is self-dual beyond balanced EC (U shows it fails for non-EC networks; unbalanced EC
  networks are open).
- K8 is answered negatively outside WF-FC (pairs in §5); a pair with rows and columns of Γ pairwise linearly
  independent on both networks was not found (random search, n = 4, 570 s).
- A Petri-invariant precluster (class-span closure, Remark 2.5): its relation to Conradi–Flockerzi clusters was not
  worked out.

## Files

- `scripts/petri.py` — library (imports `../../scripts/crn.py`): reverse-dual check, clusters, EC sets, twins,
  FC/EC/balanced, rank equations, T-traps, S/T-components, token game (liveness of bounded nets), profile rule,
  class-span closure.
- `scripts/check_dictionary_examples.py` — D1–D4, O1–O3, P1 on the manuscript's networks and B₂, B₃
  (`logs/check_dictionary_examples.log`, all PASS, < 1 s).
- `scripts/enumerate_ec.py` + `scripts/run_enum.sh` — exhaustive EC/FC families, claims K1–K3, K7, RT, reports
  K4–K6, K8, K9 (`logs/enum_ec_*.log`).
- `scripts/check_toric.py` + `scripts/run_toric.sh` — Theorem T on the families with sampled k (`logs/toric_*.log`).
- `scripts/check_absolute_consistency.py` — sympy: U, U*, W*, F, D (`logs/check_absolute_consistency.log`).
- `scripts/check_general.py` (+ `scripts/run_rest.sh`, which also runs `run_toric.sh`) — Theorems A, B, Propositions H
  and G, reports K6, K8, K9(iv), on general enumeration families with (H1), (H2): (2,2,5), (2,2,6), (2,3,4),
  (3,2,4), (3,3,3), (4,1,8), (4,2,3), 18 509 networks, all PASS (`logs/general_*.log`).
- `scripts/check_k8.py` (pairs for K8, `logs/check_k8.log`, all PASS) and `scripts/k8_search.py` (random search).
- Not completed: the ordinary EC family n = 5, molecularity ≤ 2, r ≤ 5 (`logs/enum_ec_n5_m2_r5_ord1.log`) hit its
  30-minute timeout (one partition of the species carries almost all the work); n ≤ 4 families are complete.
- `candidates.md` — candidates for agent E. `PROGRESS.md`. No `NEEDS.md`: every source used is open.
- Sources (git-ignored): `certificates/sources/fcbook/` (Desel–Esparza, free PDF from the author's page),
  `certificates/sources/jbdo/` (arXiv:2602.02352v2), `certificates/sources/ccs/` (Campos–Chiola–Silva 1991, author's
  copy).
