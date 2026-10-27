# When is the dual persistent? Rank one solved, drain sets in every rank

Agent Q1, **Claude Opus 5.5 (effort xhigh)**, 27 October 2026. Second round of the gac-duality work on result 17
(`problems/iosif-siphon-precluster-duality/`). This file does not quote or use the author's e-mail. Status tags:
**proved** (complete proof here), **checked** (exact computation on the finite families listed, not a proof),
**conjectured**, **refuted** (explicit counterexample, checked exactly), **open**.

## Summary

| # | Statement | Status |
|---|---|---|
| S1 | **Rank-one face criterion.** M a conservative network of rank 1. Every face of a class (a set Z of species that can vanish together at an endpoint) has a type IN / OUT / ESC / DOM computed from the sources of M. M is persistent for **all** rate constants iff every face is IN; if some face is OUT or ESC, M is non-persistent for **every** choice of rate constants; if all non-IN faces are DOM and lie at the same end, persistence holds on a nonempty open set of rate constants and fails on another | proved; exact endpoint expansions checked on every face of every dual of §2 (3 seeded random points each) |
| S1-N | Translation to N (N with s = 1, reaction vectors b_j a): N* is persistent for all rate constants iff for every nonempty set Z of reactions of N all oriented the same way along a, the species (with a_i ≠ 0) of smallest total coefficient in the sources of Z are all net products of Z | proved |
| S1-WR | **Dichotomy.** If N has s = 1 and is weakly reversible (more generally: in each linkage class the two complexes that are extreme along a are sources), no face of N* is DOM. So N* is persistent either for all rate constants or for none | proved; no DOM face in any WR dual of §2 (e.g. 3 544 for (2,4,6), 2 742 for (2,5,5)); G5 on the 761 of 18 802 general rank-one networks satisfying the weaker hypothesis |
| C-pure | s = 1 and no reaction of N contains one of its net products in its source (e.g. N pure) ⇒ N* persistent for all rate constants | proved; checked (G3) |
| C-2 | n = 2, s = 1 (any N with (H1)–(H3)): N* persistent for all rate constants iff in every reaction the source contains strictly more of the consumed species than of the produced one | proved; checked (G4, all 7 265 two-species networks of the general families) |
| X | **Smallest example: N = 2X₂ ⇄ X₁+X₂** (reversible, one linkage class, δ = 0, (H1)–(H3)): N* = k₂ → k₁, 2k₁+k₂ → k₁+2k₂ satisfies k₂' = k₂(κ₂k₁² − κ₁); for every κ the classes k₁ + k₂ = T < (κ₁/κ₂)^{1/2} lose k₂. N* has no drain species and no drain set | proved; checked (sympy) |
| D | **Drain sets** (any s): Z, weights w > 0, degrees d with (D1) a pure-source draining reaction d_i X_i → … for each i ∈ Z and (D2) every reaction increasing ⟨w, x_Z⟩ has Σ_{i∈Z} y(i)/d_i > 1. Then for every rate constants x_Z → 0 from every point of every class with ⟨w, x_Z⟩ small; Z contains no conservation-law support and lies in a critical siphon. Generalises E's drain species (|Z| = 1) | proved |
| D-count | On E's WR families with (H1)–(H3): duals with a drain set 1 735 (s = 2: 347, s = 3: 894, s = 4: 493) versus 303 with a drain species; drain sets occur with molecularity 2 and at s = 4, where no drain species occurs | checked (|Z| ≤ 4) |
| C | **Smallest drain set without drain species: N = 2X₃ → X₂+X₃ → X₁+X₃ → 2X₃** (one 3-cycle, s = 2, δ = 0, molecularity 2): V = k₂ + 2k₃ satisfies V' = −κ₁k₃ − κ₂k₂ + 2κ₃k₁²k₂k₃ in N*, so N* loses k₂ and k₃ for every rate choice | proved; checked (sympy) |
| A-D | "No drain set (or drain species) ⇒ N* persistent for all rate constants" | **refuted**, already for s = 1, one linkage class, δ = 0 (example X). The missing condition at s = 1 is "no ESC face" |
| SOME | A network with (H1)–(H3) whose dual is persistent for some rate constants and not for others: impossible for WR N with s = 1 (S1-WR); exists outside WR: N = X₁+X₂ → 2X₂+2X₃, X₁+X₂+2X₃ → 2X₁, dual persistent iff κ₁ > κ₂; and N = X₁+X₂ → 2X₂, X₁+X₂ → 2X₁, dual persistent iff κ₁ = κ₂ | proved; checked |
| P | A face with two vanishing species can decide while all singletons pass: N = 2X₁+X₃ ⇄ 3X₁+2X₂, 2X₂+X₃ ⇄ X₁+4X₂ (WR, s = 1) | proved; checked |
| Q | C-pure fails at s = 2 outside WR: N = X₁ → X₂, X₁ → X₃, X₂+X₃ → 2X₁ (pure, (H1)–(H3)) has a pure dual with (k₂ − k₁)' = (κ₃ − κ₂)k₃, non-persistent unless κ₂ = κ₃ | proved; checked |
| open | s ≥ 2: an exact criterion; whether the all/none dichotomy persists for WR N; a rigorous facet version of S1 (Conjecture F) | open (§6) |

Counts at s = 1 (WR, (H1)–(H3), `logs/enum_s1_*.log`; E's records give the same numbers on the overlapping
families): N* is persistent for all rate constants in 160 of the 587 rank-one WR networks of E's record files (summed
over files, overlaps counted twice; E2's independent, deduplicated decision in `../computation/logs/s1_persistence.log`
gives 150 ALL / 411 NONE / 0 SOME among 561 WR networks), 381 of 3 544 in our family (2,4,6), 662 of 2 742 in (2,5,5),
245 of 829 in (3,4,5); never "some but not all".

## 0. Setting

N: species X₁..Xₙ, reactions R_j: y_j → y'_j, Γ = Yp − Ye of rank s. N*: species k₁..k_r (one per reaction of N),
reactions ρ_i (one per species X_i of N) with source Σ_j (Ye)_{ij} k_j and product Σ_j (Yp)_{ij} k_j; Γ(N*) = Γᵀ.
(H1) N conservative ⇔ N* consistent; (H2) N consistent ⇔ N* conservative; (H3) rows of Γ distinct ⇔ the reaction
vectors of N* are distinct. Mass action for N*: ρ_i has rate κ_i k^{(Ye)_{i·}}. In the source's language
(arXiv:2603.08767, §3) the dual system is the one in which the *rate constants k_j of N evolve* and the *old
concentrations play the role of rate constants κ_i*.

A conservative mass-action system is **persistent** (for a given rate vector) if no positive solution has an
ω-limit point on the boundary of the orthant, i.e. liminf x_i(t) > 0 for all i (solutions are bounded). "Persistent
for all / some / no rate constants" refers to N* with κ ∈ ℝⁿ_{>0}, and to every positive initial condition (every
compatibility class).

## 1. Rank one: the face criterion (Theorem S1)

Let M be a mass-action network with p species, rank Γ(M) = 1 and a positive conservation law. Write the reaction
vectors as c_j g with g ∈ ℤᵖ primitive; conservativity means g has entries of both signs. Put N₋ = {i : g_i < 0},
N₊ = {i : g_i > 0}. A positive class is a segment x(ξ) = x₀ + ξg, ξ ∈ [ξ₋, ξ₊], and ξ' = F(ξ) = Σ_j c_j κ_j x^{y_j}
(y_j the source of reaction j).

**Faces.** An *upper face* is a nonempty Z ⊆ N₋ (σ = +1), a *lower face* a nonempty Z ⊆ N₊ (σ = −1). The upper
endpoint x⁺ = x(ξ₊) of a class vanishes exactly on some upper face Z and is positive elsewhere; conversely every
x⁺ ≥ 0 with zero set exactly an upper face Z is the upper endpoint of the class through x⁺ − εg (ε > 0 small).
Likewise for lower endpoints. For the face Z let o_Z(j) = Σ_{i∈Z} (y_j)_i (the Z-order of the source of j),
μ = min{o_Z(j) : c_j ≠ 0}, J = {j : c_j ≠ 0, o_Z(j) = μ}; j ∈ J is *inward* if σc_j < 0 (it increases the species of
Z) and *outward* if σc_j > 0. Let β_j ∈ ℤ^{Z^c} be the source y_j restricted to the species outside Z. The face is

- **IN** if every reaction of J is inward;
- **OUT** if every reaction of J is outward;
- **ESC** if J has both kinds and some outward β_j is not in the convex hull of the inward β's;
- **DOM** if J has both kinds and every outward β_j lies in the convex hull of the inward β's.

**Lemma 1.1 (endpoint expansion).** At an upper endpoint x⁺ with zero set Z, along x = x⁺ − εg (ε ↓ 0):
F = ε^μ L(κ, x⁺) + O(ε^{μ+1}) with
L(κ, x⁺) = Σ_{j∈J} c_j κ_j Π_{i∈Z} |g_i|^{(y_j)_i} Π_{i∉Z} (x⁺_i)^{(y_j)_i}.
The same holds at a lower endpoint along x = x⁻ + εg. The endpoint attracts the nearby points of its class if σL > 0
and repels them if σL < 0.
*Proof.* x_i = ε|g_i| for i ∈ Z and x_i = x⁺_i − σεg_i > 0 otherwise, so κ_j x^{y_j} = ε^{o_Z(j)} κ_j
Π_{i∈Z}|g_i|^{(y_j)_i} Π_{i∉Z}(x⁺_i)^{(y_j)_i} (1 + O(ε)). At the upper end ξ increases toward the endpoint, so
F > 0 near it (σL > 0) moves interior points to it; at the lower end the signs are reversed. ∎

**Lemma 1.2.** For a fixed κ, M is persistent iff no endpoint of any class is attracting, i.e. iff for every face
and every point x⁺ of it, the lowest-order nonzero term of F along the class has the inward sign (or F ≡ 0 on the
class).
*Proof.* A scalar autonomous ODE on a compact segment: every solution is monotone and converges. A solution
starting near an attracting endpoint converges to it (a boundary point). Conversely, if a positive solution
converges to an endpoint, F has the outward sign on a punctured one-sided neighbourhood (F is a polynomial in ε, so
its sign is eventually constant, and it cannot be inward), so the endpoint is attracting. ∎

**Theorem S1 (proved).** Let M be conservative of rank one.
1. M is persistent for every κ ∈ ℝ^{reactions}_{>0} iff every face is IN.
2. If some face is OUT, then for every κ every class whose endpoint has that zero set has an attracting endpoint.
   If some face is ESC, then for every κ there is a class (depending on κ) with an attracting endpoint. In both
   cases M is non-persistent for every κ.
3. If every face is IN or DOM, some face is DOM, and all DOM faces are upper faces (or all lower faces), then the
   set of κ for which M is persistent contains a nonempty open set, and so does its complement.
4. (DOM faces at both ends.) Persistence for some κ is then decided by finitely many polynomial inequalities; it
   can hold on a hypersurface only (example B, §4).

*Proof.* (1) If every face is IN, σL < 0 for every κ and every x⁺ (a sum of negative terms), so no endpoint is
attracting (Lemma 1.2). If a face is not IN, pick an outward j₀ ∈ J and the point x⁺ = (0 on Z, 1 off Z); for κ_{j₀}
large and the other κ fixed, σL > 0, so that class has an attracting endpoint.
(2) OUT: σL > 0 for every κ and x⁺. ESC: let j₀ be outward with β₀ ∉ conv{β_q : q inward}; a separating vector
w ∈ ℝ^{Z^c} gives ⟨w, β₀⟩ > max_q ⟨w, β_q⟩. At x⁺ = (0 on Z, e^{tw} off Z) every outward term of σL is positive and
σL ≥ C₀ κ_{j₀} e^{t⟨w,β₀⟩} − Σ_q C_q κ_q e^{t⟨w,β_q⟩} → +∞ as t → ∞, so for every κ some class has an attracting
endpoint.
(3) Say all DOM faces are upper (outward = {c_j > 0}). Put κ_j = η for c_j > 0 and κ_j = 1 for c_j < 0. Lower faces
are IN, hence repelling for all κ. On an upper DOM face write each outward β_j = Σ_q λ_{jq} β_q (λ_j in the simplex,
q inward); weighted AM–GM gives u^{β_j} ≤ Σ_q λ_{jq} u^{β_q} for all u > 0, hence
σL ≤ Σ_q (η Σ_j C_j λ_{jq} − C_q) u^{β_q} < 0 for all u once η is below a constant depending on the (finitely many)
faces. More generally, after the AM–GM reduction: σL ≤ Σ_q (Σ_j κ_j C_j λ_{jq} − κ_q C_q) u^{β_q}, so the finitely many strict linear
inequalities Σ_j κ_j C_j λ_{jq} < κ_q C_q (all DOM faces, all inward q) define an open set of κ on which M is
persistent, nonempty because it contains the κ above. The complement contains an open set by (1) (the inequality
σL > 0 at a fixed point is open in κ).
(4) Example B below. ∎

So **at s = 1 the boundary between "persistent for all rates" and "extinction for every rate" is the face test**,
which is the parallel-sweep test of Craciun–Nazarov–Pantea / Pantea (lower-endotactic networks, Pantea 2012,
Def. 3.5) restricted to the vectors v = 1_Z, Z ⊆ N₋ or Z ⊆ N₊ — the only directions in which a rank-one class can
approach the boundary (all vanishing coordinates go to zero at the same rate). Lower-endotacticity in Pantea's sense
(all v ≥ 0) is therefore sufficient but not necessary at s = 1.

**Remark (which classes are affected).** An endpoint has zero set Z with |Z| ≥ 2 only when the ratios x₀,i/|g_i|
(i ∈ Z) tie, a codimension-one set of classes. So OUT/ESC faces with |Z| ≥ 2 produce extinction only in a thin set
of classes, but they do destroy persistence as defined (every positive initial condition). Example P (§3) shows such
faces can be decisive.

Checks: `scripts/s1lib.py` (classification, exact endpoint expansion); claim E3 of `scripts/enum_s1.py` compares,
for every face of every enumerated dual and three seeded random (κ, x⁺), the lowest nonzero coefficient of the exact
expansion of F with ε^μ L and the predicted sign (IN repelling, OUT attracting): PASS on all families.

## 2. Duals of rank-one networks

Let N have s = 1: every reaction vector is b_j a with a ∈ ℤⁿ primitive. Then Γ(N*) = b aᵀ: for M = N*, g = b
(indexed by the reactions of N) and c = a (indexed by the species of N). Under (H2) for N, b has both signs, so N* is
conservative; (H1) for N means a has both signs; (H3) for N means the a_i are pairwise distinct (the rows of Γ are
a_i bᵀ). A face of N* is a nonempty set Z of reactions of N **all oriented the same way along a** (all b_j > 0 or
all b_j < 0); the Z-order of the dual reaction ρ_i is the **load** ℓ_Z(i) = Σ_{j∈Z} (y_j)_i, the total coefficient of
X_i in the sources of Z; ρ_i is active iff a_i ≠ 0, and inward iff X_i is a net product of the reactions of Z
(sign a_i = sign b_Z).

**Proposition S1-N (proved; Theorem S1 read on N).** For N with s = 1 and (H2), N* is persistent for every choice of
rate constants iff
> **(C)** for every nonempty set Z of reactions of N with a common orientation along a, every species X_i with
> a_i ≠ 0 that minimises ℓ_Z(i) is a net product of the reactions of Z.

If (C) fails at Z and the minimisers are all net reactants of Z (OUT), every rate choice gives extinction of
{k_j : j ∈ Z} in every class reaching that face; if mixed (ESC/DOM), Theorem S1 (2)–(4) applies.

**Theorem S1-WR (dichotomy; proved).** Let N have s = 1 and suppose that in every linkage class the two complexes
that are extreme along a are sources of reactions — e.g. N weakly reversible. Then N* has no DOM face. Consequently
(if N is consistent) N* is persistent either for every choice of rate constants (iff (C)) or for none.

*Proof.* The complexes of a linkage class L lie on a line y_L + ta; write φ_i(t) = (y_L)_i + t a_i for the coefficient
of X_i in the complex at position t. Suppose a face Z, say of forward reactions (b_j > 0), were DOM: there are an
outward species o (a_o < 0, net consumed by Z) and inward species q (a_q > 0) in J, and weights λ_q ≥ 0, Σλ_q = 1, with
β_o = Σ_q λ_q β_q. Let v denote the formal species Σ_q λ_q X_q, so (y)_v := Σ_q λ_q (y)_q and a_v = Σ_q λ_q a_q > 0.
Then (y_j)_o = (y_j)_v for every reaction j ∉ Z, and Σ_{j∈Z}(y_j)_o = μ = Σ_{j∈Z}(y_j)_v (all of o and the q are in J).
In a linkage class L the function φ_o − φ_v is affine in t with slope a_o − a_v < 0, so it vanishes at no more than
one position; hence at most one position of L is the source of a reaction outside Z. The extreme forward position
t_max(L) is a source (hypothesis) and only of backward reactions, which are not in Z; so t_max(L) is that position,
φ_o(t_max) = φ_v(t_max), and every other position carries only reactions of Z. If L has no reaction in Z, both of its
extreme positions are sources of reactions outside Z, a contradiction; so every reaction of Z leaves a position
t < t_max(L), where φ_o(t) − φ_v(t) = (t − t_max)(a_o − a_v) > 0. Summing over Z ≠ ∅ gives Σ_Z (y_j)_o > Σ_Z (y_j)_v,
contradicting the equality of the orders. Backward faces are symmetric. Theorem S1 then leaves only the cases "all
faces IN" and "some face OUT or ESC". ∎

Neither (H1) nor (H3) is used. The hypothesis cannot be dropped: examples S and B (§4) have DOM faces.

**Corollary C-pure (proved).** If s = 1 and no reaction of N has one of its net products in its source (in
particular if N is pure), then N* is persistent for every choice of rate constants.
*Proof.* For a face Z of forward reactions, the species with a_i > 0 have load 0 (forward reactions net produce them),
while those with a_i < 0 occur in every forward source with coefficient ≥ |a_i| b_j ≥ 1. So every face is IN. ∎

**Corollary C-2 (proved).** If n = 2 and s = 1 (any N with (H1)–(H3)), N* is persistent for every choice of rate
constants iff every reaction's source contains strictly more of the species it consumes than of the species it
produces. Otherwise, if N is weakly reversible, N* is non-persistent for every choice.
*Proof.* With one species in each sign class of a, the load comparison over Z is a sum of the comparisons over its
reactions, so (C) reduces to singletons. The last sentence is S1-WR. ∎

**Examples (all checked exactly, `logs/examples.log`).**
- **X: N = 2X₂ ⇄ X₁+X₂** (reversible, one linkage class, δ = 0, (H1)–(H3)). N* = k₂ → k₁ (κ₁),
  2k₁+k₂ → k₁+2k₂ (κ₂), k₁ + k₂ = T, and k₂' = k₂(κ₂k₁² − κ₁) ≤ −(κ₁ − κ₂T²)k₂ on every class with κ₂T² < κ₁: k₂ → 0
  exponentially, for every rate vector (in the classes with κ₂T² > κ₁ the interior equilibrium k₁ = (κ₁/κ₂)^{1/2} is
  globally attracting). The face {k₂} is ESC; N* has neither a drain species nor a drain set. This is the
  **smallest weakly reversible network with (H1)–(H3) whose dual is non-persistent for every choice of rate
  constants** (two species, two reactions, molecularity 2), smaller than G's 3A ⇄ 2A+B. By C-2 the failing reaction
  is X₁+X₂ → 2X₂: it consumes X₁ and produces X₂, each with coefficient 1 in its source.
- **M: 3A ⇄ 2A+B** (agent G): the face of the reaction 2A+B → 3A is OUT (E's drain species).
- **P: N = 2X₁+X₃ ⇄ 3X₁+2X₂, 2X₂+X₃ ⇄ X₁+4X₂** (WR, s = 1, (H1)–(H3), a = (1, 2, −1)): every singleton face of N*
  is IN, but the face formed by the two forward reactions (loads X₁: 2, X₂: 2, X₃: 2 — a tie between the net
  reactant X₃ and the net products X₁, X₂) is ESC, so N* is non-persistent for every rate choice.
  No such network occurs in the enumerated families (molecularity ≤ 5 needed here).

**Counts (checked, `logs/enum_s1_*.log`; generator: all WR networks with s = 1 and (H1)–(H3), complexes of
molecularity ≤ m, ≤ R reactions, up to species permutation; the numbers coincide with E's records on (2,3,6): 460,
(3,3,5): 93, (3,2,6): 12, (4,3,4): 1).**

| (n, m, R) | networks | N* persistent for all κ | for none | of which OUT face | some-not-all | singleton faces decide |
|---|---|---|---|---|---|---|
| (2,2,4) | 12 | 4 | 8 | 0 | 0 | yes |
| (2,3,6) | 460 | 115 | 345 | 267 | 0 | yes |
| (2,4,6) | 3 544 | 381 | 3 163 | 1 833 | 0 | yes |
| (2,5,5) | 2 742 | 662 | 2 080 | 1 601 | 0 | yes |
| (3,2,6) | 12 | 4 | 8 | 0 | 0 | yes |
| (3,3,6) | 426 | 101 | 325 | 130 | 0 | yes |
| (3,4,5) | 829 | 245 | 584 | 291 | 0 | yes |
| (4,3,4) | 1 | 1 | 0 | 0 | 0 | yes |

(Families overlap; (2,3,5), (3,3,5), (3,4,4) are sub-families also logged.) Persistence of N* for all rate constants
holds for a minority (17–35 %) of rank-one WR networks; drain species (E.4) explain almost none of the failures
(1 of 345 in (2,3,6)): the failures are OUT faces whose minimal sources contain other species, and ESC faces.

**General (non-WR) rank-one networks** (`scripts/enum_s1_general.py`, `logs/enum_s1_general_*.log`, 18 802 networks
with (H1)–(H3) in the families (2,2,4), (2,3,4), (2,4,3), (3,2,4), (3,3,4)): G2 (no DOM in WR duals), G3 (C-pure),
G4 (C-2), G5 (no DOM under the weaker hypothesis of S1-WR, 761 networks) all PASS. DOM faces occur only outside that
hypothesis; in these families they occur only together with OUT/ESC faces or at both ends (example B); the
one-sided example S needs molecularity 4.

## 3. Drain sets in every rank (Theorem D)

**Definition.** Let M be a mass-action network. A *drain set* is a nonempty set Z of species with weights
w ∈ ℝ^Z_{>0} and integers d ∈ ℤ^Z_{≥1} such that, with V = ⟨w, x_Z⟩:
- (D1) for every i ∈ Z some reaction has source exactly d_i X_i (no other species) and strictly decreases V
  (⟨w, Δ|_Z⟩ < 0);
- (D2) every reaction that strictly increases V has a source y with Σ_{i∈Z} y_i/d_i > 1.

E's drain species (Prop. E.4) are the drain sets with |Z| = 1.

**Theorem D (proved).** Let M be conservative and Z a drain set. For every choice of rate constants and every
compatibility class P there is δ > 0 such that every solution in P with V(x(0)) < δ has x_Z(t) → 0. In particular M
is non-persistent for every choice of rate constants (start at x_Z = η·1, x_{Z^c} = 1 with η small). Moreover Z
contains the support of no nonzero semi-positive conservation law, and the zero set of any ω-limit point of such a
solution is a critical siphon containing Z.

*Proof.* Let T bound the coordinates on P. A reaction increasing V contributes at most C x_Z^{β} with β = y|_Z,
θ = Σ_i β_i/d_i > 1 (bounding the other factors by T). For i ∈ Z the reaction of (D1) contributes at most
−c_i x_i^{d_i}, c_i > 0; reactions with ⟨w, Δ|_Z⟩ ≤ 0 contribute ≤ 0. With u_i = x_i^{d_i} and U = max_i u_i,
x_Z^{β} = Π_i u_i^{β_i/d_i} ≤ U^θ ≤ U^{θ−1} Σ_i u_i, and U^{θ−1} → 0 as V → 0. Hence for V ≤ δ (δ small),
V' ≤ −(c/2) Σ_i x_i^{d_i} ≤ 0, with equality only if x_Z = 0. So {V ≤ δ} is forward invariant and V decreases to a
limit L; if L > 0, Σ_i x_i^{d_i} is bounded below on {L ≤ V ≤ δ} and V' ≤ −c' < 0 forever, impossible. So V → 0.
If z ≥ 0, z ≠ 0, were a conservation law supported in Z, z·x would be conserved and positive along a positive
solution, but z·x ≤ (max z / min w)V → 0. The last statement is the siphon property of ω-limit zero sets
(Angeli–De Leenheer–Sontag 2007, Prop. 1 / Deshpande–Gopalkrishnan 2014) and the same conservation argument. ∎

**In terms of N (for M = N*).** Z is a set of reactions of N with weights w and degrees d such that (D1) for each
j ∈ Z some species X_i of N occurs in the source of R_j with coefficient d_j and in no other source, and is net
consumed by the w-weighted reactions of Z (Σ_{j'∈Z} w_{j'} Γ_{ij'} < 0); (D2) every species of N net produced by the
w-weighted reactions of Z satisfies Σ_{j∈Z} (y_j)_i / d_j > 1. For |Z| = 1 this says: R_j has a "private" reactant of
multiplicity d₀ and all its net products already sit in its source with multiplicity > d₀ — an autocatalytic
reaction (2A+B → 3A in G's example).

**Example C (checked).** N = 2X₃ → X₂+X₃ → X₁+X₃ → 2X₃ (one 3-cycle, s = 2, δ = 0, molecularity 2, (H1)–(H3)),
N* = k₃ → k₂ (κ₁), k₂ → k₁ (κ₂), 2k₁+k₂+k₃ → k₁+k₂+2k₃ (κ₃). With V = k₂ + 2k₃,
V' = −κ₁k₃ − κ₂k₂ + 2κ₃k₁²k₂k₃ (identity), so {k₂, k₃} is a drain set (d = (1, 1), w = (1, 2)) and N* loses k₂ and k₃
for every rate choice; neither k₂ nor k₃ alone is a drain species. This is the smallest WR network in E's families
whose dual has a drain set but no drain species.

**Counts (checked, `logs/scan_records_wr.log`, all of E's WR record files, networks with (H1)–(H3), |Z| ≤ 4).**

| s | WR networks | dual with a drain species (E.4) | dual with a drain set | drain set, no drain species |
|---|---|---|---|---|
| 1 | 587 | 1 | 1 | 0 |
| 2 | 19 900 | 199 | 347 | 148 |
| 3 | 31 885 | 103 | 894 | 791 |
| 4 | 17 422 | 0 | 493 | 493 |

Drain sets occur in the molecularity-2 families ((3,2,5) with zero complex, (3,2,6), (4,2,6), (5,2,5), (5,2,6) and
(6,2,5) homogeneous), where E found no drain species, and at s = 4. Checks R1 (every drain species is a drain set),
R2 (at s = 1 a drain set forces the verdict "none") and R3 (no DOM face at s = 1) PASS.

## 4. Some rates versus all rates (task 3)

- **s = 1, N weakly reversible: impossible** (Theorem S1-WR): the dual is persistent for all rate constants or for
  none, and the borderline cases L = 0 do not arise (IN, OUT and ESC are strict).
- **s = 1, N not weakly reversible: occurs.** Example S: N = X₁+X₂ → 2X₂+2X₃, X₁+X₂+2X₃ → 2X₁ ((H1)–(H3), s = 1,
  two linkage classes, not WR). N* = k₁+k₂ → 2k₂ (κ₁), k₁+k₂ → 2k₁ (κ₂), 2k₂ → 2k₁ (κ₃), and
  k₁' = k₂((κ₂ − κ₁)k₁ + 2κ₃k₂). If κ₁ > κ₂ every positive solution converges to the interior equilibrium
  (κ₁ − κ₂)k₁ = 2κ₃k₂; if κ₁ ≤ κ₂, k₁' ≥ 2κ₃k₂² > 0 and k₂ → 0. So **N* is persistent exactly for κ₁ > κ₂** (an open
  set). Example B: N = X₁+X₂ → 2X₂, X₁+X₂ → 2X₁ ((H1)–(H3)); N* = k₁+k₂ → 2k₂, k₁+k₂ → 2k₁ with k₁' = (κ₂ − κ₁)k₁k₂:
  persistent **only on the hypersurface κ₁ = κ₂** (DOM faces at both ends).
- **s ≥ 2, N weakly reversible: open.** No exact decision procedure is available; see §6 for the facet heuristic.

## 5. What this means for Iosif's question

In the source's language the dual system is the one in which the old rate constants k_j evolve and the old
concentrations are parameters. Persistence of N* means: **no rate constant of N is driven to zero, whatever the
concentrations used as parameters.** The results above say:

1. **Small codimension does not help, small dimension does.** For s = 1 the answer is complete and combinatorial
   (condition (C) on N), and for weakly reversible N it is a dichotomy: either every choice of the parameters
   keeps all k_j positive, or every choice lets some k_j die in some class. The obstruction is local
   autocatalysis in N: a set Z of same-direction reactions whose least-loaded source species include a net reactant.
2. **Persistence of N (now known for all WR networks) does not transfer to N***, even in the smallest cases:
   2X₂ ⇄ X₁+X₂ (s = 1) and the 3-cycle 2X₃ → X₂+X₃ → X₁+X₃ → 2X₃ (s = 2) are complex balanced with one linkage class
   and deficiency 0, and their duals lose a species for every choice of rate constants. Together with P's
   Corollary P (WF-EC ⇒ N* persistent for all rate constants), this places the boundary: EC-type structure (distinct
   sources with disjoint supports, no autocatalysis) on the persistent side; autocatalytic reactions with
   low-multiplicity net reactants on the extinct side.
3. The right certificates for N* are **local** (faces, minimal loads, drain sets), not the siphon/conservation-law
   certificate of Angeli–De Leenheer–Sontag: at s = 1, 345 of 460 duals in the family (2,3,6) are non-persistent,
   almost all without a drain species.

## 6. Open points and a conjecture for s ≥ 2

- **Facets.** For general s, a set Z of species of N* whose face in a class is a facet (rank Γ(N*)[Z, ·] = 1, i.e. by
  G's Prop. 6.2 the reactions Z of N span a line) is approached along a single direction, so Lemma 1.1 applies near
  relative-interior points of the facet with the same IN/OUT/ESC/DOM types. An IN facet repels locally; an OUT facet
  attracts locally, but turning this into non-persistence needs control of the motion along the facet (it may carry
  the solution toward lower-dimensional faces). *Conjecture F:* for WR N and any s, if some facet of N* is OUT (resp.
  ESC) then N* is non-persistent for every (resp. for some) choice of rate constants. Not proved.
- **Sufficient conditions at s = 2.** N* persistent for all rate constants if N* has no critical siphon (ADLS; this is
  P's WF-EC case) or if its stoichiometric subnetworks are lower-endotactic (Pantea 2012, Thm 5.1, with bounded
  trajectories, which conservativity of N* gives); neither is necessary. An exact criterion at s = 2 is open.
- **Dichotomy for s ≥ 2.** Whether the dual of a WR network can be persistent for some but not all rate constants
  when s ≥ 2 is open.
- **One linkage class.** Example X (s = 1) and example C (s = 2) have one linkage class and deficiency 0: one linkage
  class gives no extra persistence of the dual.
- **C-pure does not extend to s = 2 without weak reversibility (example Q, proved, checked).** N = X₁ → X₂,
  X₁ → X₃, X₂+X₃ → 2X₁ is pure with (H1)–(H3) and s = 2 (not WR); N* = k₁+k₂ → 2k₃ (κ₁), k₃ → k₁ (κ₂), k₃ → k₂ (κ₃)
  is pure too, and (k₂ − k₁)' = (κ₃ − κ₂)k₃ (identity). If κ₃ ≠ κ₂, k₂ − k₁ is monotone and bounded, so
  ∫₀^∞ k₃ dt < ∞; k₃' is bounded, hence k₃ → 0: N* is non-persistent for every rate vector with κ₂ ≠ κ₃ (the hypersurface κ₂ = κ₃
  was not analysed). Whether "N weakly
  reversible and pure ⇒ N* persistent for all rate constants" holds for s ≥ 2 is open (candidate Q1-6).

## 7. Relation to the literature (novelty)

- The endpoint expansion for scalar rank-one systems is elementary; the face criterion is the restriction of the
  parallel-sweep (lower-endotactic) test of Craciun–Nazarov–Pantea 2013 and Pantea 2012 (arXiv:1103.0603, Def. 3.5,
  read) to the directions 1_Z. For the stochastic rank-one case Wiuf–Xu (arXiv:2012.07954) classify
  one-dimensional networks by similar boundary order comparisons; Pantea–Voitiuk (arXiv:2208.06310) classify
  multistationarity for rank-one networks. No deterministic classification of persistence for all / some / no rate
  constants in rank one was found (web search, 27 Oct 2026), nor anything on duals.
- Deshpande–Gopalkrishnan (arXiv:1309.3957, read: Def. 3.1, Thm 6.2) prove that networks without *drainable*
  siphons are persistent. Their "drainable" is a reachability notion; at s = 1 every face of a conservative
  consistent network admits a draining reaction, so their criterion almost never applies to N*, while S1 is exact.
  Drain sets (Theorem D) are a rate-independent extinction certificate in the opposite direction; I found no prior
  statement of it (the single-species case is E's Prop. E.4).
- AI-campaign scan (`gh api users/mmislan/repos`, `users/ipitchford/repos`, `gh search repos`): nothing on duals or
  rank-one persistence (mmislan's endotactic permanence concerns endotactic networks, which N* almost never is).

## Files

- `scripts/s1lib.py` — rank-one face analysis, exact endpoint expansion, dual.
- `scripts/enum_s1.py` — WR rank-one families; claims E1–E5; logs `logs/enum_s1_<n>_<m>_<R>.log`.
- `scripts/enum_s1_general.py` — all rank-one networks with (H1)–(H3); claims G1–G5; logs `logs/enum_s1_general_*.log`.
- `scripts/drainsets.py` — drain-set decision (exact Fourier–Motzkin, imports E's `endo.py`).
- `scripts/scan_records.py` — E's WR records (WSL `~/gacE-data`): drain species vs drain sets, s = 1 verdicts;
  claims R1–R3; log `logs/scan_records_wr.log` (77 s).
- `scripts/examples.py` — sympy identities and classifications for X, M, P, S, B, C; `logs/examples.log` (22 PASS).
- Sources (git-ignored): `certificates/sources/` (Pantea arXiv:1103.0603, Deshpande–Gopalkrishnan arXiv:1309.3957).
- `candidates.md` (for E2), `PROGRESS.md`.
