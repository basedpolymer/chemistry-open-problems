# Is persistence of the dual all-or-none in dimension two?

Agent D2, **Claude Opus 5.5 (effort xhigh)**, 27 October 2026. Third round of the gac-duality work on result 17
(`problems/iosif-siphon-precluster-duality/`), brief `../../BRIEF-D2.md`. This file does not quote or use the author's
e-mail. Status tags: **proved** (complete proof here), **checked** (exact computation on the finite family, not a
proof), **conjectured**, **refuted**, **open**. Paths are relative to `gac-duality/round3/dichotomy/`.

## Summary

| # | Statement | Status |
|---|---|---|
| WR-Σ | N weakly reversible ⇒ Yp = Ye·Σ (Σ a 0/1 matrix sending each reaction to one starting at its product); the dual reaction vector of ρ_i is s_i(Σ − I), a linear function of its source s_i; equivalently k_j' = G_{y'_j}(P) − G_{y_j}(P) with P_c = Π_{j: y_j = c} k_j and G_c = P_c ∂Ψ/∂P_c, Ψ(P) = Σ_i κ_i P^{Y_i}. Under (H3) the dual sources are pairwise distinct | proved |
| W | **N weakly reversible ⇒ no siphon facet of any class of N* is DOM, in every rank s** (generalises Q1's S1-WR, which is the case s = 1) | proved; checked on all 18 861 duals (claim W) |
| A | Rank 2: a boundary point of an ω-limit set of an interior point is an equilibrium or lies in an open edge whose zero set is a siphon; every vertex in it is an equilibrium | proved |
| P2 | **ALL certificate (rank 2).** If (A1) every siphon facet has a one-way edge flow (or is IN with a two-way flow that cannot vanish identically), (A2) every vertex equilibrium has a repelling function (Prop. R), (A3) the directions carrying a siphon facet do not positively span the plane, then N* is persistent for every κ | proved |
| E, H, CM, L, S | **NONE certificates (rank 2)**: edge equilibrium on an OUT/ESC facet (E); tropical Hurwitz vertex (H); center-manifold vertex (CM); attracting linear Lyapunov function (L); attracting sector (S). Each gives, for every κ, a class with an interior trajectory converging to the boundary | proved |
| B | N* splitting into two rank-one blocks: persistent iff both blocks are (Q1's Theorem S1 decides each) | proved |
| count | On E2's **14 825 undecided** WR networks with (H1)–(H3), s = 2: **ALL 1 885, NONE 12 769, SOME 0, undecided 171** (list in `logs/undecided_s2.txt`); on all 18 861: ALL 4 873, NONE 13 817; agrees with every E2 certificate; no network meets both an ALL and a NONE certificate | checked (exact, 19 s) |
| Q1-7 | Dichotomy at s = 2 (no SOME for WR N) | **open in general; checked** on 14 654 of the 14 825 networks and on example U (proved by hand). Not refuted: no SOME candidate survived |
| U | **Vertex ties exist.** N = 3X₃ → 2X₁+X₂ → X₂+2X₃ → 3X₃: at the vertex (T,0,0) the decisive coefficient is 2κ₁ − κ₂ (no class scaling removes it); N* is non-persistent for every κ, by different mechanisms for κ₂ > 2κ₁, κ₂ < 2κ₁, κ₂ = 2κ₁ | proved |
| F | Conjecture F (a realizable OUT facet forces NONE) | **checked**: 2 515 of the 18 861 duals have an OUT siphon facet: 2 508 certified NONE, 7 undecided, 0 ALL. **Proved** when the edge flow of the OUT facet is two-coloured, or zero with μ = 1 (Prop. E); open when it is one-way |

## 0. Setting

N: species X_1..X_n (index i), reactions R_j: y_j → y'_j (index j = 1..r), Ye, Yp, Γ = Yp − Ye of rank s = 2. N*: species
k_1..k_r, reactions ρ_i: Σ_j (y_j)_i k_j → Σ_j (y'_j)_i k_j, rate κ_i k^{s_i}, s_i = row i of Ye; Γ(N*) = Γᵀ. (H1)–(H3) for N
make N* conservative and consistent, so its positive classes are compact convex polygons.

**Planar coordinates.** Factor Γ = Π Ĝ with Π (n × 2, rows π_i) and Ĝ (2 × r, columns g_j), so Γ_ij = ⟨π_i, g_j⟩. The class
of k⁰ is the polygon P(k⁰) = {X ∈ ℝ² : k_j(X) := k⁰_j + ⟨g_j, X⟩ ≥ 0}, and N* reads X' = Σ_i κ_i k(X)^{s_i} π_i. The
edges of P have inward normals among the g_j (reaction vectors of N), the dual reactions move X along the π_i (species
of N).

**Faces.** A *facet* zero set is a nonempty set Z of reactions of N with positively parallel reaction vectors (a subset of
a *direction class* D); a *vertex* zero set is a set Z_v of reactions of rank 2 whose g_j lie in an open half-plane.
Z is a siphon of N* iff every species changed by the reactions of Z occurs in one of their sources.

**Lemma 2.1 (proved; free boundary positions).** (i) Every p ≥ 0 with zero set exactly a facet set Z lies in the relative
interior of the Z-edge of some positive class. (ii) Every p ≥ 0 with zero set exactly a vertex set Z_v is a vertex of
some positive class. In particular the coordinates of an edge point off Z, and of a vertex off Z_v, can be any positive
vector.
*Proof.* (i) Choose z ∈ ℝⁿ with ⟨γ_j, z⟩ < 0 for j ∈ Z (possible: the γ_j, j ∈ Z, are positively parallel) and small; then
k⁰ = p − Γᵀz is positive and p lies in its class; the face {k_Z = 0} of that class has dimension 2 − rank Γ[:, Z] = 1.
(ii) Same with ⟨g_j, X⟩ < 0 for all j ∈ Z_v (the open half-plane condition). ∎

## 1. Weak reversibility in the dual

**Proposition WR-Σ (proved).** If N is weakly reversible, every product complex y'_j is the source y_{σ(j)} of some
reaction σ(j); with Σ_{σ(j), j} = 1 (else 0), Yp = Ye Σ, hence Γ = Ye(Σ − I) and the reaction vector of ρ_i is
s_i(Σ − I): the product of ρ_i has coefficient s_{i,σ(j)} on k_j. Writing P_c = Π_{j: y_j = c} k_j and
Ψ(P) = Σ_i κ_i Π_c P_c^{Y_ic} (Y the complex matrix), the dual system is k_j' = G_{y'_j} − G_{y_j} with
G_c = P_c ∂Ψ/∂P_c = Σ_i κ_i Y_ic P^{Y_i}: each k_j gains the "potential" of the head complex of R_j and loses that of its
tail; in log-coordinates p = log P, G = ∇_p ψ with ψ convex. Under (H3) two dual reactions with equal sources would
have equal reaction vectors, so **the dual sources are pairwise distinct**.
*Proof.* k_j' = Σ_i κ_i k^{s_i}Γ_ij = Σ_i κ_i P^{Y_i}(Y_{i,y'_j} − Y_{i,y_j}) because k^{s_i} = Π_j k_j^{Y_{i,y_j}} = P^{Y_i}. ∎

**Theorem W (proved; any rank).** Let N be weakly reversible and Z a siphon facet set of N* in any rank (reaction vectors
γ_j = b_j a, b_j > 0, j ∈ Z). Then Z is not DOM: there are no outward species o (net consumed by Z) and inward species q
(net produced) of minimal load μ = Σ_{j∈Z}(y_j)_o = Σ_{j∈Z}(y_j)_q with (y_l)_o = Σ_q λ_q (y_l)_q for every l ∉ Z
(λ in the simplex).
*Proof.* Put φ(y) = y_o − Σ_q λ_q y_q on complexes. Then φ(y_l) = 0 for every source y_l, l ∉ Z, and each reaction j ∈ Z
changes φ by b_j(a_o − Σλ_q a_q) < 0. Let j ∈ Z have source y with φ(y) ≤ 0. Its product y + b_j a has φ < 0, so it is
not the source of a reaction outside Z; by weak reversibility it is the source of some reaction, which lies in Z, and its
product has φ smaller still. This gives an infinite sequence of distinct complexes, impossible. So φ(y_j) > 0 for all
j ∈ Z, and Σ_{j∈Z} φ(y_j) = ℓ_Z(o) − Σ_q λ_q ℓ_Z(q) = μ − μ = 0, a contradiction. ∎

For s = 1 this is Q1's Theorem S1-WR (and by Theorem S1, the dichotomy). For s = 2 it removes κ-dependence at the
level of facets; Example U (§5) shows that κ-dependence can still arise at vertices.

## 2. Where an ω-limit set can meet the boundary (rank 2)

Fix κ and a positive class P (a polygon). Let Ω be the ω-limit set of an interior point.

**Lemma A (proved).** If p ∈ Ω ∩ ∂P then either p is an equilibrium, or p lies in an open edge whose zero set Z is a
siphon of N*, the field is tangent to that edge, and the total trajectory of p stays in the open edge. Every vertex in Ω is
an equilibrium.
*Proof.* Ω is compact and invariant: through p there is a total trajectory in Ω ⊆ P (unique, the field is polynomial).
So f(p) and −f(p) lie in the tangent cone of P at p. At a vertex this cone is pointed, so f(p) = 0. In an open edge
it is a half-plane, so f(p) is parallel to the edge. If Z is not a siphon, some dual reaction whose source avoids Z
produces a species of Z; at p every reaction with source avoiding Z fires, none of them consumes Z, and all move k_Z along
the same positive vector (rank one), so the normal speed at p is > 0, a contradiction. On a siphon edge the plane
{k_Z = 0} is invariant under the polynomial flow, so the total trajectory of p stays in the edge (it cannot reach a vertex,
which would be an equilibrium). ∎

On a siphon edge only the dual reactions with source avoiding Z act, and they move along the edge: the **edge flow** is
F(p) = Σ_{i : ℓ_Z(i) = 0} c_i κ_i p^{s_i} with π_i = c_i d. It is *one-way* if all c_i have one sign (then no
equilibrium in the open edge, and the downstream vertex is an equilibrium), *two-way* if both signs occur, *zero* if there
is no such reaction (the edge consists of equilibria).

**Lemma BM (Butler–McGehee for the semiflow on P; proved).** Let e ∈ Ω be an equilibrium with a neighbourhood U whose
closure contains no compact invariant set (made of total trajectories in P) other than {e}, and Ω ≠ {e}. Then Ω contains
a point q ≠ e whose forward orbit converges to e.
*Proof.* x(t) enters and leaves U infinitely often (otherwise Ω ⊆ cl U would be {e}). Take t_n → ∞ with x(t_n) → e and
s_n < t_n the last time before t_n with x(s_n) ∈ ∂U; t_n − s_n → ∞ by continuity of the flow at the equilibrium e. A
limit point q of x(s_n) lies in ∂U ∩ Ω and its forward orbit stays in cl U, so its ω-limit set is invariant in cl U, i.e.
equals {e}. ∎

**Theorem P2 (ALL certificate; proved).** Suppose that for every κ and every class:
(A1) every siphon facet Z has a one-way edge flow, or is IN with a two-way edge flow that cannot vanish identically on an
edge; (A2) every vertex equilibrium e admits a neighbourhood U_e and a C¹ function V_e ≥ 0 with V_e > 0 on int P and
V_e' > 0 on U_e ∩ int P; (A3) some edge of P is not a siphon edge. Then no interior trajectory of P has an ω-limit point on
∂P. Condition (A3) holds for every class if the direction classes carrying a siphon facet do not positively span ℝ².
*Proof.* Suppose Ω ∩ ∂P ≠ ∅. By Lemma A and the 1-D dynamics on siphon edges, Ω contains a boundary equilibrium e (a
vertex equilibrium, or an equilibrium in an open IN edge, finitely many by (A1)). Every boundary equilibrium e is
interior-repelling: V = V_e (vertex) or V = ξ, the normal coordinate of an IN edge (ξ'/ξ > 0 near the open edge) increases
along interior trajectories near e, so no interior point converges to e; and e is isolated in the sense of Lemma BM (a
total interior trajectory in a small cl U would have its ω-limit in {V constant} ∩ ∂P, which contains no invariant set
besides e; boundary orbits near e move monotonically along the edges). Since x(0) is interior, Ω ≠ {e}; Lemma BM gives
q ∈ Ω converging to e, q ≠ e, necessarily on ∂P, hence (Lemma A) on an open siphon edge flowing into e; its α-limit e' is
the next boundary equilibrium upstream on the boundary circle, and e' ∈ Ω. Repeating, the chain moves monotonically
around ∂P through siphon edges only, and since there are finitely many boundary equilibria it closes up: every edge of P
is a siphon edge, contradicting (A3). The last sentence: the inward normals of the edges of a bounded polygon positively
span ℝ². ∎

**Proposition R (repelling functions; proved).** At a vertex equilibrium with zero set Z_v let a = k_{e1}, b = k_{e2} be
the coordinates of the two extreme normals (local coordinates in which P is the quadrant), so k_l = α_l a + β_l b
(α_l, β_l ≥ 0) for l ∈ Z_v, and every changed species contributes C_i(k) Π_{l∈Z_v}(α_l a + β_l b)^{(s_i)_l} with
C_i(k) → κ_i v^{w_i} > 0. If, for V = a^λ b^μ (V'/V = λa'/a + μb'/b) or V = λa + μb, every negative monomial of V'/V
(resp. V') is *strictly Newton-dominated* by the positive ones — some point p of the convex hull of the positive
exponents has p ≤ e componentwise, p ≠ e — then V' > 0 on U ∩ int P for every κ and every class.
*Proof.* Weighted AM–GM: x^p ≤ Σ λ_k x^{e_k}, and x^e = x^p x^{e−p} with x^{e−p} → 0 in the open quadrant; finitely many
terms, continuous positive coefficients. ∎

## 3. Certificates of extinction for every κ (rank 2)

In each case "tropical" means: the off-face coordinates are taken as exp(tθ) (allowed by Lemma 2.1), t → ∞, and a sign
condition holds for the terms of largest ⟨θ, w⟩ with margin ≥ 1 (an exact LP for θ, solved by Fourier–Motzkin). Every κ
then gets a class (t large depending on κ) with an interior trajectory converging to the boundary.

- **Prop. E (edge equilibrium; proved).** Z an OUT or ESC siphon facet, θ with the minimal-order transverse terms
  dominated by outward ones, and either (two-coloured) two edge-flow terms of opposite signs tied for the largest
  ⟨θ, s_i⟩, all other edge-flow terms below; or (zero flow) no edge flow and μ = 1. Then for every κ some edge point p*
  with F(p*) = 0 has an attracting normal speed. Two-coloured: along exp(tθ + u), u in a segment of the tied face, F
  changes sign while the outward term dominates the normal speed; the zero p* is an equilibrium with ξ' = ξ^μ(L + O(ξ)),
  L < 0 near p*; a Ważewski box argument (top side entering, side exits open and disjoint) gives an interior trajectory
  staying in the box, along which ξ decreases to 0. Zero flow, μ = 1: normally attracting line of equilibria; the
  along-edge drift is O(ξ) with ξ exponentially decreasing.
- **Prop. H (tropical Hurwitz vertex; proved).** Vertex equilibrium; the linearisation is Σ_i C_i u_i v_iᵀ over the
  species of local order 1 (rank-one, Metzler), trace Σ C_i ⟨u_i, v_i⟩, determinant Σ_{i<i'} C_iC_{i'} (u_i × u_{i'})(v_i × v_{i'}).
  θ makes the leading trace term negative and the leading nonzero determinant term positive ⇒ Hurwitz ⇒ the vertex
  attracts a neighbourhood in P. (E2's Prop. E2.7 is the case θ = ±(1,…,1).)
- **Prop. CM (center manifold; proved).** The a-column of the linearisation vanishes identically (no local monomial a),
  δ = ∂b'/∂b has a negative leading term, and with m, k the smallest a-powers on the axis b = 0 in a' and in b': if
  m < k (or b' ≡ 0 on the axis) the leading coefficient of a' on the axis is negative; if m = k,
  Σ C_iC_{i'}(u_i × u_{i'}) > 0 over (species with a^k) × (species with b) — which is the sign of −(f_k|δ| + βg_k), the
  leading coefficient of the reduced equation a' = c a^m on the center manifold b = h(a) ≈ (g_k/|δ|)a^k ≥ 0. Then every
  point of P near the vertex converges to it (stable fibres project the quadrant onto the branch a > 0).
- **Prop. L (attracting linear function; proved).** V = λa + μb (λ, μ > 0) with every positive monomial of V' Newton-
  dominated by negative ones, strictly or on the boundary of the Newton polygon with the dominating coefficients
  tropically larger. Then V' < 0 on U ∩ int P; the small triangles {V ≤ c} are forward invariant and every interior point
  in them has its ω-limit set in ∂P (LaSalle).
- **Prop. S (attracting sector; proved).** With b = a u, take u ∈ [e^{−t(ω+ε)}, e^{−t(ω−ε)}] (or u ∈ [0, e^{−t(ω−ε)}]):
  a' = a^{d_a}(…) has a negative leading term with margin, and u' = a^{d_u}(…) is positive on the lower ray and negative on
  the upper one (among the tied leading terms in (θ, ω), the lowest u-power is positive and the highest is negative). The
  sector near the vertex is forward invariant, a decreases to 0, and b = au → 0.
- **Prop. B (products; proved).** If N* splits into two rank-one blocks, classes are products of segments and N* is
  persistent iff both blocks are; Q1's Theorem S1 (`s1lib.classify`) decides each block.

Drain species, drain sets and E2.7 are special cases or are subsumed on this family (claims E2a, E2b).

## 4. Results on E2's networks (`scripts/certify_s2.py`, `logs/certify_s2.log`, 19 s)

All 18 861 WR networks with (H1)–(H3) and s = 2 of the first-round records (rebuilt exactly by
`scripts/extract_undecided.py`, which reproduces E2's counts, PASS). Claims, all **PASS**:

- W: no DOM siphon facet in any of the 18 861 duals (Theorem W).
- K: no network meets a NONE certificate and the ALL conditions at once.
- E2a / E2b: the 2 988 networks certified ALL by E2 (no critical siphon) are certified ALL here; the 1 048 certified NONE
  by E2 are certified NONE here.
- F: 2 515 duals have a realizable OUT siphon facet: 2 508 NONE, 7 undecided, 0 ALL.
- D: on E2's 14 825 undecided networks: **ALL 1 885, NONE 12 769, SOME 0, undecided 171**. First applicable NONE
  certificate: E 3 513, H 2 897, CM 5 289, L 535, S 449, B 86.

Smallest examples (exact witnesses in the log):
- ALL: N = 2X₃ → X₂ → 2X₁+X₃ → 2X₃ (E2's smallest undecided), N* = 2k₃ → 2k₂, k₂ → k₁, 2k₁+k₃ → k₂+2k₃. By hand: the
  only invariant edge is k₃ = 0 (IN, flow toward (T,0,0)); the vertex (T,0,0) is a saddle with eigenvalues −κ₂ (along the
  edge) and κ₃T² (inward), V = k₃ with k₃'/k₃ = κ₃k₁² − 2κ₁k₃ > 0 near it; the edges k₁ = 0, k₂ = 0 are not siphons.
- NONE [E]: N = 2X₃ → 2X₂+X₃ → X₁+X₂+X₃ → 2X₃; the facet {k₃} is ESC with zero edge flow and μ = 1.
- NONE [H]: N = X₃ → X₂ (R₁), 2X₃ → X₁+X₂ (R₂), X₂ → X₃ (R₃), X₁+X₂ → 2X₃ (R₄), N* = k₄ → k₂, k₃+k₄ → k₁+k₂,
  k₁+2k₂ → k₃+2k₄: the vertex {k₁, k₄} is Hurwitz in the classes where k₃ is small at the vertex (θ = (0, −1) on (k₂, k₃)).
- NONE [CM]: N = 2X₃ → X₂+X₃ → X₁+X₂ → 2X₃, N* = k₃ → k₂, k₂+k₃ → k₁+k₂, 2k₁+k₂ → k₁+2k₃. At (0,T,0), a = k₁, b = k₃:
  a' = κ₂Tb − κ₃Ta², b' = −(κ₁+κ₂T)b + 2κ₃Ta², center manifold b ≈ 2κ₃Ta²/(κ₁+κ₂T), reduced a' ≈ κ₃Ta²(κ₂T − κ₁)/(κ₁+κ₂T):
  attracting for every class with κ₂T < κ₁.
- NONE [L]: N = 2X₃ → 2X₂+X₃ → 2X₁+X₃ → 2X₃: at (T,0,0), V = k₂ + 2k₃ has V' = −2κ₁k₃² − 2κ₂k₂² + 2κ₃T²k₂k₃ + …, negative
  definite when κ₃²T⁴ < 4κ₁κ₂.
- NONE [S]: N = 3X₃ → X₂+2X₃ → 2X₁+X₂ → 3X₃: at (T,0,0) a slow curve k₃ ≈ (3κ₃T³/κ₂)k₂ attracts and a decreases along it
  (θ = −1/3 for T, ω = 1).
- NONE [B]: N = X₄ ⇄ X₃, 2X₂ ⇄ X₁+X₂: N* is the product of k₃ ⇄ k₄-type block 2k₃+k₄ → k₃+2k₄, k₄ → k₃ (Q1's example X:
  NONE) and k₁ ⇄ k₂ (ALL).
- Undecided (smallest): example U below (decided by hand).

The 171 undecided networks (`logs/undecided_s2.txt`): 136 with one vertex equilibrium that none of R, H, CM, L, S
decides (κ-ties of the kind in example U, or degenerate corners), 18 with an IN facet whose two-way flow fails the
non-vanishing test, 17 with an IN facet of zero edge flow.

## 5. The dichotomy question Q1-7

**Status: open in general; no counterexample.** Every decided network is ALL or NONE, rate-independently, and the
proofs are uniform in κ. Theorem W removes the s = 1 mechanism (DOM facets) in every rank. What remains possible at
s = 2 is κ-dependence at vertices, and it does occur:

**Example U (proved).** N = 3X₃ → 2X₁+X₂ (R₁), X₂+2X₃ → 3X₃ (R₂), 2X₁+X₂ → X₂+2X₃ (R₃): one 3-cycle, WR, (H1)–(H3), s = 2.
N* = 2k₃ → 2k₁ (κ₁), k₂+k₃ → k₁+k₃ (κ₂), 3k₁+2k₂ → 3k₂+2k₃ (κ₃), classes k₁+k₂+k₃ = T, and exactly
k₂' = −κ₂k₂k₃ + κ₃k₁³k₂², k₃' = −2κ₁k₃² + 2κ₃k₁³k₂². At the vertex (T,0,0) put a = k₂, b = k₃, K = κ₃k₁³ ≤ κ₃T³:
- z = a/b satisfies the exact identity **z' = a[(2κ₁ − κ₂) + Kz(1 − 2z)]**: the sign of 2κ₁ − κ₂ is decided by κ alone,
  no class scaling changes it (both rates are ∝ b at the vertex). The edge k₂ = 0 is an OUT siphon facet flowing into the
  vertex (k₃' = −2κ₁k₃² there).
- κ₂ > 2κ₁: on {0 < z ≤ z₀, b small} with z₀ = min(1/2, (κ₂−2κ₁)/(2κ₃T³), (κ₁/(2κ₃T³))^{1/2}) we have z' < 0 and
  b' = b²(−2κ₁ + 2Kz²) < 0, so every such interior point converges to (T,0,0): non-persistent, in every class.
- κ₂ < 2κ₁: with u = b/a, exactly u' = a[(κ₂ − 2κ₁)u² − Ku + 2K] and a' = a²(K − κ₂u). The positive root u* of
  (2κ₁−κ₂)u² + Ku − 2K exceeds K/κ₂ iff K < κ₂²/κ₁; for κ₃T³ < κ₂²/(2κ₁) a sector u ∈ [u₁, u₂] around u* with u₁ > K/κ₂
  is forward invariant and a' ≤ −ca² there: non-persistent for small T.
- κ₂ = 2κ₁: u' = aK(2 − u), and a' = a²(K − κ₂u) < 0 near u = 2 for κ₃T³ < κ₂: non-persistent for small T.
So N* is NONE, but by a mechanism that switches with the sign of 2κ₁ − κ₂: no single tropical direction certifies it,
which is why the automatic certificates leave it undecided. A SOME network would need such a vertex tie whose
"repelling" sign is not rescued by another mechanism in some class. In every case examined the terms that do not depend
on the class scale push toward the vertex and win in classes of small total mass:

**Remark (3-cycles, proved).** Let N be a 3-cycle c₁ → c₂ → c₃ → c₁ and give species i the triple (p_i, q_i, r_i) of its
coefficients in (c₁, c₂, c₃); ρ_i is (p,q,r) → (q,r,p) and the classes are triangles k₁+k₂+k₃ = T. At the vertex (0,0,T)
with a = k₁, b = k₂, the local effective rates are K_i = κ_i T^{r_i}. (i) The linearisation has determinant
Σ K_iK_{i'}(1 − r_i − r_{i'}) over species (1,0,r_i) × (0,1,r_{i'}): it is positive only through the pair r = r' = 0, which
then wins as T → 0 (E2.7); so hyperbolic sinks at 3-cycle vertices never depend on κ. (ii) If some species (1,0,·)
and no species (0,1,·) exist (center direction b), the center-manifold coefficient has the sign of
Σ K_m K_X (r_m + b₀ r_X − b₀) over species (0,b₀,r_m) × (1,0,r_X) (b₀ the least power of b, ≥ 2, among species with
p = 0); the bracket is increasing in r_m and r_X, so its smallest value is that of the pair of least T-exponent, which
dominates as T → 0: if that value is negative the vertex attracts for T small (every κ), otherwise every term is ≥ 0 and
the leading reduced coefficient is ≥ 0 for every T and κ — again no dependence on κ. Ties in κ therefore need two species with the same T-weight acting on
different coordinates, as in example U.

This is evidence for the dichotomy, not a proof.

## 6. Conjecture F

Checked: no dual with a realizable OUT siphon facet is certified ALL; 2 508 of 2 515 are certified NONE, 7 undecided.
Proved by Prop. E when the OUT facet's edge flow is two-coloured (its Newton polytope has leading terms of both signs in
some direction) or zero with μ = 1. **Open** when the edge flow is one-way: then trajectories slide along the facet into
its downstream vertex (an equilibrium), and the outcome is decided at that vertex (example U: the OUT edge k₂ = 0 attracts
into the vertex iff κ₂ > 2κ₁, and a second mechanism gives extinction otherwise).

## 7. What is left

- A proof (or a counterexample) of the dichotomy at s = 2: the open part is κ-ties at vertex equilibria (example U type);
  a proof would show that whenever a tie can be repelling, the class-independent corner terms attract in some class.
- The 171 networks of `logs/undecided_s2.txt`: case splits on the sign of tied κ-combinations (as in example U) would
  decide most of the 136 vertex cases; the 35 IN facets with two-way or zero flow need Theorem P2 for continua of edge
  equilibria (the tube argument works for μ = 1).
- Conjecture F for one-way edge flows; s ≥ 3.

## 8. Relation to the literature

The tools are classical: Butler–McGehee lemma and acyclicity theorems (Butler–Waltman 1986, Hale–Waltman 1989, here
re-proved in the needed form), center manifold reduction, quasi-homogeneous (Newton polygon) analysis of planar
singular points, the parallel-sweep idea of Craciun–Nazarov–Pantea and Pantea 2012. Theorem W and the formula
k_j' = G_{y'_j} − G_{y_j} for duals of WR networks were not found in the literature or the AI campaigns scanned
(`gh api users/mmislan/repos`, `users/ipitchford/repos`, `gh search repos` on "dual reaction network persistence",
"reverse-dual mass-action", 27 October 2026); see `NEEDS.md` for the two persistence references not on disk.

## Files

- `scripts/d2lib.py` — planar factorisation Γ = ΠĜ, dual, generic edge/vertex data.
- `scripts/d2cert.py` — facets and vertices of all classes, Theorem P2 + Prop. R (ALL), Props E, H, CM, L, S, B (NONE).
- `scripts/extract_undecided.py` — rebuilds E2's list (imports `s2_persistence.py` unchanged), PASS, 47 s; data in WSL
  `~/gacD2-data/s2_all.pkl.gz`.
- `scripts/certify_s2.py` — claims W, K, E2a, E2b, F, D (all PASS, 19 s); `logs/certify_s2.log`, `logs/undecided_s2.txt`.
- `scripts/try1.py` … `try4.py`, `survey.py`, `find_e2_mixed.py` — exploratory runs (no claims).
