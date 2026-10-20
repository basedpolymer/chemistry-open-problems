# The operating specification of the generalized Petlyuk arrangement and its infinite-stage limit configuration

Agent D (design and limit configuration), Claude Opus 5.5 (effort xhigh). **Version 1.1, 2026-10-20 18:46** (v1 18:36).
Status: complete proofs of everything stated as a lemma, proposition or theorem below (not yet independently
verified); numerical evidence is collected separately in §8 and is never used in a proof. Interface summary:
`hypotheses.md`. Scripts: `scripts/` (high-precision floating point, mpmath; evidence only).

Sources used: the preferred split of `../../paper.tex` (Theorem `achieve`, version 1, independently verified),
GPT-6-Astra's slack deformation and Lemmas 1–2 (`../../total/attainment-astra/attainment.md` §2–3), Agent S's
limit relation (`../../total/sections/sections.md` Lemma 11) and Agent N's forced impurities
(`../../total/network/network.md` §11.2). Everything that is used is re-proved here, in the form needed.

---

## 0. Result

Let N ≥ 3, α_1 > α_2 > ⋯ > α_N > 0, z ∈ (0,1)^N with Σ_p z_p = 1, q ∈ ℝ; θ_1 > ⋯ > θ_{N−1} the roots of
Σ_p α_p z_p/(α_p − θ) = 1 − q (θ_j ∈ (α_{j+1}, α_j)), H_j = Σ_{p>j} α_p z_p/(θ_j − α_p), H* = max_j H_j,
e_j = H* − H_j ≥ 0.

**Theorem D.** Fix a forcing pattern π ∈ {B, T}^{N−1}. There are η_1 > 0 and C_1 ≥ 0, depending only on (α, z, q),
with the following property. Let δ > 0, η̄(δ) := min(η_1, δ/(2C_1)) (η̄ = η_1 if C_1 = 0), and choose
- η_c ∈ (0, η̄(δ)] for every non-binary column c = [i, m] (m − i ≥ 2),
- η_p ∈ (0, η̂_p) for every binary column [p, p+1], where η̂_p > 0 is the explicit bound (4.3) (π_p = B) or (4.6)
  (π_p = T), computed from the flows of F.

Then the specification F = F(δ, η, π) of Definition 5.1 (all stream flows of the generalized Petlyuk arrangement)
is well defined, every stream flow is positive, every junction and side-draw stage is adiabatic, and the only heat
input is the reboiler duty R = H* + δ. Its limit configuration 𝓛 = 𝓛(δ, η, π) (Definition 6.1: limit net flows,
junction and product compositions) satisfies, for every section, the infinite-stage limit relations with the pinch
modes of Table 1 and strictly positive pivots, and hypotheses (H1)–(H8) of `hypotheses.md` hold. Every product p has
limit impurity (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p, so the purities tend to 1 as η_bin → 0.
More generally, all conclusions hold for every *admissible* parameter (§5: validity of the non-binary columns,
positive binary margins, η_p < η̂_p), a finite list of explicit inequalities on the flows of F that can be checked by
evaluating F; the choice above is one admissible choice for every δ > 0.

The proof occupies §§2–6: Theorem D = Proposition 5.4 + Theorem 5.5 + Theorem 6.2.

**What is new with respect to the plan (BRIEF-common §3, Step 1).** (a) Astra's conditions reduce to positivity of
the deformed split (Proposition 3.3): root interlacing, Ψ_T(α_m) > 0, Ψ_B(α_i) > 0 and the four cone weights follow.
(b) The binary margins satisfy κ_p(η_nb; δ) = κ_p(η_nb; 0) + δ exactly, κ_p(0; 0) = e_p (Proposition 5.4), so
η_nb ≤ η̄(δ), linear in δ, suffices; **no hierarchy between the η's of parents and children is needed.**
(c) The forcing of a binary column only *increases* its margin (Propositions 4.2, 4.3): η_p need not be small
relative to κ_p, only below explicit bounds that keep the forced root away from the next volatility. (d) In every
section the pinch is the smallest major root (top sections) or the largest (bottom sections); every trace mode lies
outside the major spectrum; this gives the trace signs (H7) uniformly (Theorem 6.2). (e) Any forcing pattern works;
default: all B (product 1 pure, product p ≥ 2 receives η_{p−1} of p − 1).

---

## 1. Setting and notation

**1.1 The arrangement** (HS 2003 Part 3, Fig. 1; version 1's FTC). Columns c = [i, m], 1 ≤ i < m ≤ N; c is
*non-binary* if m − i ≥ 2 and *binary* otherwise. Column [i, m] has a top section T[i,m], a bottom section B[i,m]
and a junction (feed) stage J[i,m], an adiabatic equilibrium stage. The top end of T[i,m] is the stage J[i,m−1]
(if m − 1 > i) or the product stage P_i (if m − 1 = i); the bottom end of B[i,m] is J[i+1,m] or P_m. P_1 is a
total condenser, P_N a partial reboiler (an equilibrium stage; the only heat input), and P_p (2 ≤ p ≤ N−1) an
adiabatic side-draw equilibrium stage, between B[p−1,p] above and T[p,p+1] below. Thermal couplings: the vapor
leaving the top stage of T[i,m] enters its top-end stage and the liquid of that stage returns to T[i,m]; the liquid
leaving the bottom stage of B[i,m] enters its bottom-end stage and the vapor of that stage returns. The process feed
(flow 1, composition z, quality q) enters J[1,N]. In the finite-stage model each section has n_σ ≥ 0 equilibrium
stages; the results here concern the flows and the infinite-stage limit and do not depend on the n_σ.

**1.2 Sections.** A section σ has vapor flow V_σ, liquid flow L_σ and constant upward net component flow
w_σ ∈ ℝ^N (w = d in a top section, w = −b in a bottom section), Σ_i w_{σ,i} = V_σ − L_σ. Put

  Ψ_σ(t) = V_σ − Σ_i α_i w_{σ,i}/(α_i − t),  M_σ = V_σ diag(α) − w_σ αᵀ,  G_x(t) = Σ_i α_i x_i/(α_i − t),

and Ψ^{(k)}_σ, G^{(k)}_x the same sums without the k-th term (terms with w_i = 0, resp. x_i = 0, are omitted, so
these are continuous at t = α_i for such i). The stage relation of an n-stage section is x_top ∝ M_σ^{n+1} x_bot,
x_bot the liquid in equilibrium with the vapor entering the bottom, x_top the liquid entering the top
(`sections.md` Lemma 11; attainment.md (1)). For a junction or product stage the composition used at a section end
is that stage's liquid composition (all of them are equilibrium stages; for the total condenser it is the reflux
composition = distillate composition).

**1.3 Columns and feeds.** Column c has top net product d_c, bottom net product b_c, feed f_c = d_c + b_c, vapors
V_{T_c}, V_{B_c} with V_{T_c} − V_{B_c} = U_c, D_c = Σ d_c, B_c = Σ b_c, L_{T_c} = V_{T_c} − D_c,
L_{B_c} = V_{B_c} + B_c. Feeds: f_{[1,N]} = z, U_{[1,N]} = 1 − q; otherwise

  f_{[i,m]} = d_{[i,m+1]} + b_{[i−1,m]},  U_{[i,m]} = V_{T[i,m+1]} − V_{B[i−1,m]}           (1.1)

(terms of absent parents omitted). The feed function Φ_c(t) := U_c − Σ_i α_i f_{c,i}/(α_i − t) then satisfies the
exact bookkeeping identity

  Φ_{[i,m]} = Ψ_{T[i,m+1]} − Ψ_{B[i−1,m]},  Φ_c = Ψ_{T_c} − Ψ_{B_c}.                          (1.2)

---

## 2. Algebra of one section

**Lemma 2.1 (spectrum).** Let S = {i : w_{σ,i} ≠ 0}. Then det(Vφ − M_σ) = V^{N−1} Π_{k∉S}(φ − α_k) · P_S(φ), where
P_S(φ) := Ψ_σ(φ) Π_{i∈S}(φ − α_i) is a polynomial of degree |S| with leading coefficient V whose zeros are exactly
the roots of Ψ_σ. Hence the eigenvalues of M_σ are Vφ, φ ranging over the roots of Ψ_σ (*major modes*) and over the
α_k with w_{σ,k} = 0 (*trace modes*).

*Proof.* With A = V diag(α), det(λ − A + wαᵀ) = det(λ − A)(1 + αᵀ(λ − A)^{−1}w) = Π_i(λ − Vα_i)(1 + Σ_i α_i w_i/
(λ − Vα_i)); put λ = Vφ. P_S(α_i) = α_i w_i Π_{j∈S, j≠i}(α_i − α_j) ≠ 0 for i ∈ S, so the zeros of P_S are roots
of Ψ_σ. ∎

**Lemma 2.2 (one-signed net flows).** Let V > 0, w ≠ 0 with all nonzero entries of one sign, S as above, and let
a_1 > ⋯ > a_s be the α_i, i ∈ S.
(a) (*top type*, w ≥ 0, L = V − Σw > 0.) Ψ is strictly decreasing on each interval of ℝ ∖ {a_j}; Ψ > V on
(a_1, ∞); Ψ has exactly one root in each (a_{j+1}, a_j) and one root in (0, a_s), and no other root; Ψ > 0 on
(−∞, ρ_min) and on (a_1, ∞), where ρ_min ∈ (0, a_s) is the smallest root.
(b) (*bottom type*, w ≤ 0.) Ψ is strictly increasing on each interval of ℝ ∖ {a_j}; Ψ > V on (−∞, a_s); Ψ has
exactly one root in each (a_{j+1}, a_j) and one root in (a_1, ∞), and no other root; Ψ > 0 on (−∞, a_s) and on
(ρ_max, ∞), ρ_max the largest root. (L = V + Σ|w| > 0 automatically.)
In both cases all roots are simple and positive.

*Proof.* Ψ'(t) = −Σ_{i∈S} α_i w_i/(α_i − t)², of constant sign. (a) Near a_j, −α_j w_j/(a_j − t) → −∞ as t ↑ a_j
and +∞ as t ↓ a_j; Ψ(t) → V as t → −∞, Ψ(0) = V − Σw = L > 0, and on (a_1, ∞) every term −α_i w_i/(α_i − t) is
positive. Monotonicity gives one root per interval (a_{j+1}, a_j), one in (0, a_s), none in (−∞, 0] or (a_1, ∞);
these are s roots, all zeros of P_S. (b) Same with the signs reversed; on (−∞, a_s) every term α_i|w_i|/(α_i − t)
is positive, and Ψ → V > 0 as t → +∞. ∎

**Lemma 2.3 (eigenvectors; pinch identity).** For a mode φ of σ define e_φ ∈ ℝ^N and ℓ_φ ∈ ℝ^N by
- φ a root of Ψ_σ: (e_φ)_i = w_i/(α_i − φ), ℓ_φ = (α_i/(α_i − φ))_i (so ℓ_φ·x = G_x(φ));
- φ = α_k, w_k = 0: (e_φ)_i = w_i/(α_i − α_k) (i ≠ k), (e_φ)_k = Ψ^{(k)}_σ(α_k)/α_k, ℓ_φ = e_k (so ℓ_φ·x = x_k).

Then: (i) M e_φ = Vφ e_φ and ℓ_φᵀM = Vφ ℓ_φᵀ; (ii) 1·e_φ = L/φ; (iii) G_{e_φ}(t) = −Ψ_σ(t)/(t − φ) for all
t ∉ {α_i : i ∈ S} ∪ {φ}; (iv) ℓ_φ·e_ψ = 0 for modes φ ≠ ψ, and ℓ_φ·e_φ = −Ψ_σ'(φ) (root) resp. Ψ^{(k)}_σ(α_k)/α_k
(trace mode). Put ê_φ := (φ/L) e_φ (the normalized *pinch composition* of mode φ, 1·ê_φ = 1). If the N modes are
distinct, the ê_φ form a basis and every x ∈ ℝ^N has the expansion x = Σ_φ c_φ(x) ê_φ, c_φ(x) = ℓ_φ·x/ℓ_φ·ê_φ.

*Proof.* (i) For a root, α·e_φ = Σ α_i w_i/(α_i − φ) = V − Ψ(φ) = V, so (Me_φ)_i = Vα_i(e_φ)_i − Vw_i = Vφ(e_φ)_i;
for a trace mode α·e_φ = (V − Ψ^{(k)}(α_k)) + Ψ^{(k)}(α_k) = V, and the same computation works for i ≠ k, while
(Me_φ)_k = Vα_k(e_φ)_k − w_kV = Vα_k(e_φ)_k. For the left vectors, (ℓᵀM)_j = Vℓ_jα_j − (ℓ·w)α_j with ℓ_φ·w = V − Ψ(φ)
= V (root; then Vα_j(α_j/(α_j − φ) − 1) = Vφℓ_j) or ℓ·w = w_k = 0 (trace). (ii) For a root,
Σ_i w_i/(α_i − φ) = φ^{−1}Σ_i w_i(α_i/(α_i − φ) − 1) = φ^{−1}(V − Ψ(φ) − Σw) = L/φ; for a trace mode the same
computation with the k-term Ψ^{(k)}(α_k)/α_k gives φ^{−1}(V − Ψ^{(k)}(α_k) − Σw + Ψ^{(k)}(α_k)) = L/α_k.
(iii) Using 1/((α − φ)(α − t)) = (1/(α − t) − 1/(α − φ))/(t − φ): for a root, G_{e_φ}(t) =
((V − Ψ(t)) − (V − Ψ(φ)))/(t − φ) = −Ψ(t)/(t − φ); for a trace mode the sum over i ≠ k gives
−(Ψ^{(k)}(t) − Ψ^{(k)}(α_k))/(t − α_k) and the k-term adds Ψ^{(k)}(α_k)/(α_k − t), total −Ψ^{(k)}(t)/(t − α_k), and
Ψ^{(k)} = Ψ since w_k = 0. (iv) If φ is a root, ℓ_φ·e_ψ = G_{e_ψ}(φ) = −Ψ(φ)/(φ − ψ) = 0, and ℓ_φ·e_φ =
lim_{t→φ} G_{e_φ}(t) = −Ψ'(φ). If φ = α_k is a trace mode, ℓ_φ·e_ψ = (e_ψ)_k = w_k/(α_k − ψ) = 0 for ψ ≠ α_k. The
last statement is linear algebra (distinct eigenvalues; ℓ_φ·ê_φ ≠ 0 by (iv): simple roots resp. Ψ^{(k)}(α_k) ≠ 0,
which is exactly distinctness of α_k from the roots of Ψ). ∎

(iii) is Underwood's pinch identity (attainment.md (13), `sections.md` §5) in a normalization-free form.

**Lemma 2.4 (positive directions).** (a) If w ≤ 0 then M ≥ 0 entrywise. (b) If w ≥ 0 and L > 0, then M is invertible
and M^{−1} = V^{−1}diag(α)^{−1} + (VL)^{−1}diag(α)^{−1} w 1ᵀ ≥ 0. (c) If x = Σ_φ c_φ ê_φ with all c_φ ≥ 0 and all
ê_φ ≥ 0 appearing, then M^k x ≥ 0 for every k ∈ ℤ (M invertible for k < 0).

*Proof.* (a) M = V diag(α) + |w|αᵀ. (b) Sherman–Morrison: with A = V diag(α), 1 − αᵀA^{−1}w = 1 − Σw/V = L/V ≠ 0,
and A^{−1} + A^{−1}wαᵀA^{−1}/(L/V) is the stated matrix since αᵀA^{−1} = 1ᵀ/V. (c) M^k ê_φ = (Vφ)^k ê_φ, φ > 0. ∎

**Lemma 2.5 (regularized trace mode).** Let σ be fixed except for w_k, and suppose Ψ^{(k)}_σ(α_k) ≠ 0 at w_k = 0.
Then for |w_k| small Ψ_σ has a unique root φ_k(w) near α_k, real-analytic in w, with
α_k − φ_k = α_k w_k/Ψ^{(k)}_σ(φ_k), and the rescaled left eigenvector ℓ̃_k := ((α_k − φ_k)/α_k)ℓ_{φ_k} (= e_k at
w_k = 0) satisfies exactly

  ℓ̃_k·x = x_k + (w_k/Ψ^{(k)}_σ(φ_k)) G^{(k)}_x(φ_k),

so that d(ℓ̃_k·x) = dx_k + (G^{(k)}_x(α_k)/Ψ^{(k)}_σ(α_k)) dw_k at w_k = 0, x_k = 0. The right eigenvector e_{φ_k}
is analytic too and tends to the trace-mode vector of Lemma 2.3 (its k-entry w_k/(α_k − φ_k) = Ψ^{(k)}(φ_k)/α_k).

*Proof.* F(φ, w) := (α_k − φ)Ψ^{(k)}_σ(φ) − α_k w_k vanishes at (α_k, w_k = 0) with ∂_φF = −Ψ^{(k)}(α_k) ≠ 0:
implicit function theorem. On F = 0 with w_k ≠ 0, Ψ(φ_k) = Ψ^{(k)}(φ_k) − α_kw_k/(α_k − φ_k) = 0. Multiply
G_x(φ_k) = α_kx_k/(α_k − φ_k) + G^{(k)}_x(φ_k) by (α_k − φ_k)/α_k = w_k/Ψ^{(k)}(φ_k). ∎

---

## 3. Non-binary columns: the preferred split and Astra's slack direction

Fix a non-binary column c = [i, m] (k_c := m − i + 1 ≥ 3 components) with feed f (f_p > 0 for i ≤ p ≤ m, f_p = 0
otherwise) and U ∈ ℝ. Φ(t) = U − Σ_{p=i}^m α_p f_p/(α_p − t) is strictly decreasing between its poles, from +∞ to
−∞ on each (α_{p+1}, α_p), so it has exactly one root θ^c_p in each (α_{p+1}, α_p), i ≤ p ≤ m − 1 (the *feed
roots* of c); they are simple and real-analytic functions of (f, U).

**Lemma 3.1 (preferred split relative to the actual feed; paper.tex Theorem `achieve`).** Put
ρ_p := Π_{r=i}^{m−1}(α_p − θ^c_r)/Π_{p'∈[i,m−1], p'≠p}(α_p − α_{p'}) (i ≤ p ≤ m−1), V⁰_T := α_i f_i/ρ_i,
d⁰_p := V⁰_T ρ_p/α_p (i ≤ p ≤ m − 1), d⁰_p := 0 otherwise, b⁰ := f − d⁰, V⁰_B := V⁰_T − U. Then d⁰_i = f_i,
d⁰_p > 0 (i ≤ p ≤ m−1), b⁰_i = 0, b⁰_p > 0 (i+1 ≤ p ≤ m), V⁰_T, V⁰_B > 0, L⁰_T = V⁰_T − Σd⁰ > 0,
L⁰_B = V⁰_B + Σb⁰ > 0, and

  Ψ⁰_T(t) = V⁰_T Π_{r=i}^{m−1}(t − θ^c_r)/Π_{p=i}^{m−1}(t − α_p),  Ψ⁰_B(t) = V⁰_B Π_{r=i}^{m−1}(t − θ^c_r)/Π_{p=i+1}^{m}(t − α_p).

*Proof.* Partial fractions of R(t) = Π_r(t − θ^c_r)/Π_{p=i}^{m−1}(t − α_p) (equal degrees) give R = 1 +
Σ_p ρ_p/(t − α_p), so V⁰_T R(t) = V⁰_T − Σ_p α_p d⁰_p/(α_p − t) = Ψ⁰_T(t). In ρ_p, the numerator and the
denominator each have p − i negative factors (θ^c_r > α_p ⇔ r < p), so ρ_p > 0; hence V⁰_T > 0 and d⁰_p > 0. Since
d⁰ + b⁰ = f and V⁰_T − V⁰_B = U, Ψ⁰_T − Ψ⁰_B = Φ, so Ψ⁰_B vanishes at every θ^c_r; b⁰_i = 0, so Ψ⁰_B has poles only
at α_{i+1}, …, α_m and Ψ⁰_B·Π_{p=i+1}^m(t − α_p) is a polynomial of degree ≤ m − i with leading coefficient V⁰_B
and the m − i roots θ^c_r. If V⁰_B = 0 it would vanish identically, b⁰ = 0, contradicting b⁰_m = f_m > 0. So
Ψ⁰_B has the stated form; its residue at α_p gives −α_p b⁰_p = V⁰_B Π_r(α_p − θ^c_r)/Π_{p'∈[i+1,m], p'≠p}(α_p − α_{p'}),
whose right side has the sign of −V⁰_B (p − i negative factors above, p − i − 1 below), so every b⁰_p has the sign of
V⁰_B, and b⁰_m = f_m > 0 gives V⁰_B > 0. Finally L⁰_T = Ψ⁰_T(0) = V⁰_T Π_r θ^c_r/Π_{p<m}α_p > 0 and
L⁰_B = Ψ⁰_B(0) = V⁰_B Π_r θ^c_r/Π_{p>i}α_p > 0. ∎

**Lemma 3.2 (Astra's direction; attainment.md (3)–(5)).** If k_c = 3 put K ≡ 1 and h = 0. If k_c ≥ 4 put
β := Π_{p=i+1}^{m−1}α_p / Π_{r=i+1}^{m−2}θ^c_r and

  K(t) := (t − β) Π_{r=i+1}^{m−2}(t − θ^c_r) / Π_{p=i+1}^{m−1}(t − α_p),
  h_p := (α_p − β)Π_{r=i+1}^{m−2}(α_p − θ^c_r) / (α_p Π_{p'∈[i+1,m−1], p'≠p}(α_p − α_{p'}))  (i+1 ≤ p ≤ m−1),  h_p := 0 otherwise.

Then α_{m−1} < β < α_{i+1}, K(t) = 1 − Σ_p α_p h_p/(α_p − t), Σ_p h_p = 0, K(θ^c_r) = 0 for i+1 ≤ r ≤ m−2, and
K(θ^c_i), K(θ^c_{m−1}), K(α_i), K(α_m) > 0.

*Proof.* β = α_{i+1}Π_{r=i+1}^{m−2}(α_{r+1}/θ^c_r) < α_{i+1} and β = α_{m−1}Π_{r=i+1}^{m−2}(α_r/θ^c_r) > α_{m−1}.
K has equal degrees and leading coefficient 1 and simple poles α_{i+1..m−1} with residues α_p h_p, which is the
partial-fraction identity; K(0) = (−β)Π(−θ^c_r)/Π(−α_p) = βΠθ^c_r/Πα_p = 1 gives Σh = 0. At t ∈ {θ^c_i, α_i}, t
exceeds β, every θ^c_r (r ≥ i+1) and every α_p (p ≥ i+1): all factors are positive. At t ∈ {θ^c_{m−1}, α_m}, t is
below β, every θ^c_r (r ≤ m−2) and every α_p (p ≤ m−1): the numerator has m − i − 1 negative factors and so does
the denominator. ∎

**Definition (deformed column).** For η ≥ 0 put d := d⁰ + ηh, b := f − d = b⁰ − ηh, V_T := V⁰_T + η,
V_B := V⁰_B + η (so V_T − V_B = U), L_T := V_T − Σd, L_B := V_B + Σb. Note d_i = f_i, d_m = 0, b_i = 0, b_m = f_m
for every η (h_i = h_m = 0). The column is **valid** if

  (V)  d_p > 0 for i ≤ p ≤ m − 1  and  b_p > 0 for i + 1 ≤ p ≤ m.

(V) holds at η = 0 (Lemma 3.1), is an open condition, and holds for every η ≥ 0 when k_c = 3. Explicitly, (V) ⇔
η < η_0(c) := min({d⁰_p/|h_p| : h_p < 0} ∪ {b⁰_p/h_p : h_p > 0}) (min ∅ = +∞).

**Proposition 3.3 (a valid deformed column satisfies Astra's conditions).** Let η ≥ 0 satisfy (V). Then (a) holds,
and if η > 0, (b)–(e) hold:
(a) V_T = V⁰_T + η, L_T = L⁰_T + η, V_B = V⁰_B + η, L_B = L⁰_B + η are positive and Σd = Σd⁰;
Ψ_{T_c} = Ψ⁰_T + ηK and Ψ_{B_c} = Ψ⁰_B + ηK.
(b) Ψ_{T_c} has exactly the simple roots r_c, θ^c_{i+1}, …, θ^c_{m−2}, u_c and Ψ_{B_c} exactly the simple roots
v_c, θ^c_{i+1}, …, θ^c_{m−2}, s_c, where

  α_i > r_c > θ^c_i > v_c > α_{i+1},   α_{m−1} > u_c > θ^c_{m−1} > s_c > α_m.

In particular the middle roots are common, r_c ≠ v_c, u_c ≠ s_c, and Ψ_{T_c}(α_m) > 0, Ψ_{B_c}(α_i) > 0.
(c) Modes: T_c has the major modes r_c > θ^c_{i+1} > ⋯ > θ^c_{m−2} > u_c and the trace modes α_k, k ∉ [i, m−1];
B_c has the major modes v_c > θ^c_{i+1} > ⋯ > s_c and the trace modes α_k, k ∉ [i+1, m]. All modes of each section
are positive and pairwise distinct; the trace modes of T_c are > r_c (k < i) or < u_c (k ≥ m), those of B_c are
> v_c (k ≤ i) or < s_c (k > m).
(d) The pinch compositions ê^T_{(m)} := ê^{T_c}_{α_m}, ê^T_u := ê^{T_c}_{u_c}, ê^B_v := ê^{B_c}_{v_c},
ê^B_{(i)} := ê^{B_c}_{α_i} satisfy: ê^T_{(m)} > 0 and ê^B_{(i)} > 0 on [i, m]; ê^T_u > 0 on [i, m−1], = 0 on m;
ê^B_v > 0 on [i+1, m], = 0 on i; all vanish outside [i, m].
(e) Let X_c be the normalized vector with G_{X_c}(χ) = 0 for χ ∈ {r_c, θ^c_{i+1}, …, θ^c_{m−2}, s_c}
(attainment.md Lemma 1: since α_i > r_c > α_{i+1} > θ^c_{i+1} > ⋯ > α_{m−1} > s_c > α_m interlace, X_c exists, is
unique, X_c > 0 on [i, m], = 0 outside, with X_{c,p} ∝ Π_χ(α_p − χ)/(α_p Π_{p'≠p}(α_p − α_{p'}))). Then

  X_c = A_c ê^T_{(m)} + B_c ê^T_u = C_c ê^B_v + D_c ê^B_{(i)},  A_c + B_c = C_c + D_c = 1,  A_c, B_c, C_c, D_c > 0,

with B_c = −G_{ê^T_{(m)}}(s_c)/(G_{ê^T_u}(s_c) − G_{ê^T_{(m)}}(s_c)), C_c = G_{ê^B_{(i)}}(r_c)/(G_{ê^B_{(i)}}(r_c) − G_{ê^B_v}(r_c)).
In T_c's eigenbasis X_c has nonzero coefficients only on ê^T_{(m)} and ê^T_u; in B_c's only on ê^B_v and ê^B_{(i)}.

*Proof.* (a) Σh = 0 and the definition; Ψ_T = V⁰_T + η − Σα_p(d⁰_p + ηh_p)/(α_p − t) = Ψ⁰_T + ηK, and the same
for Ψ_B with V⁰_B + η and b⁰ − ηh. (b) T_c has w = d ≥ 0 with support [i, m−1] and L_T > 0: by Lemma 2.2(a) Ψ_T has
one root in each (α_{p+1}, α_p), i ≤ p ≤ m − 2, one in (0, α_{m−1}) and no other. By (a) and Lemmas 3.1–3.2,
Ψ_T(θ^c_r) = 0 for i+1 ≤ r ≤ m−2, so these are the roots in their intervals; Ψ_T(θ^c_i) = ηK(θ^c_i) > 0 and Ψ_T
decreases on (α_{i+1}, α_i), so its root r_c there exceeds θ^c_i; Ψ_T(θ^c_{m−1}) = ηK(θ^c_{m−1}) > 0 and Ψ_T
decreases on (−∞, α_{m−1}), so its root u_c there exceeds θ^c_{m−1} > α_m, and Ψ_T(α_m) > 0. B_c has w = −b ≤ 0
with support [i+1, m]: by Lemma 2.2(b) Ψ_B has one root in each (α_{p+1}, α_p), i+1 ≤ p ≤ m−1, one in (α_{i+1}, ∞)
and no other; the middle ones are θ^c_{i+1..m−2}; Ψ_B(θ^c_i) = ηK(θ^c_i) > 0 and Ψ_B increases on (α_{i+1}, ∞), so
its root v_c there is < θ^c_i < α_i and Ψ_B(α_i) > 0; Ψ_B(θ^c_{m−1}) > 0 with Ψ_B increasing on (α_m, α_{m−1}), so
s_c < θ^c_{m−1}. (c) Lemma 2.1 and (b); the trace modes are the α_k with zero net flow; for T_c the support of d is
[i, m−1] (d_i = f_i > 0, d_m = 0), for B_c the support of b is [i+1, m]. (d) Lemma 2.3 with L > 0: (ê^T_u)_p =
(u/L_T) d_p/(α_p − u) > 0 for p ≤ m−1 (u < α_{m−1}); (ê^T_{(m)})_p = (α_m/L_T) d_p/(α_p − α_m) > 0 for p < m and
(ê^T_{(m)})_m = Ψ_T(α_m)/L_T > 0; (ê^B_v)_p = (v/L_B) b_p/(v − α_p) > 0 for p ≥ i+1 (v > α_{i+1});
(ê^B_{(i)})_p = (α_i/L_B) b_p/(α_i − α_p) > 0 for p > i and (ê^B_{(i)})_i = Ψ_B(α_i)/L_B > 0. (e) *Interpolation*
(attainment.md Lemma 1): for X supported on [i, m], G_X(t)Π_{p=i}^m(α_p − t) is a polynomial of degree ≤ m − i;
prescribing its m − i distinct zeros χ fixes it up to a scalar, and its values at t = α_p (residues) give
X_p ∝ Π_χ(α_p − χ)/(α_p Π_{p'≠p}(α_p − α_{p'})); by interlacing, numerator and denominator have p − i negative factors
each, so X_c > 0 on [i, m] after normalization (the normalizing sum is a positive combination); the homogeneous version
shows uniqueness (the matrix with rows 1ᵀ and (α_p/(α_p − χ))_p is nonsingular). In T_c's
eigenbasis (Lemma 2.3) the coefficient of X_c on a major mode χ ∈ {r_c, θ^c_{i+1..m−2}} is ∝ G_{X_c}(χ) = 0 and on a
trace mode α_k, k ∉ [i, m], it is ∝ X_{c,k} = 0; so X_c = A ê^T_{(m)} + B ê^T_u with A + B = 1 (normalizations).
By Lemma 2.3(iii), G_{ê^T_{(m)}}(s_c) = −(α_m/L_T)Ψ_T(s_c)/(s_c − α_m) < 0 (s_c < θ^c_{m−1} < u_c, so Ψ_T(s_c) > 0)
and G_{ê^T_u}(s_c) = −(u_c/L_T)Ψ_T(s_c)/(s_c − u_c) > 0; G_{X_c}(s_c) = 0 gives the formula for B, and A, B > 0.
Similarly in B_c's eigenbasis X_c = Cê^B_v + Dê^B_{(i)}; Ψ_B(r_c) > 0 (r_c > θ^c_i > v_c), G_{ê^B_v}(r_c) =
−(v_c/L_B)Ψ_B(r_c)/(r_c − v_c) < 0 and G_{ê^B_{(i)}}(r_c) = −(α_i/L_B)Ψ_B(r_c)/(r_c − α_i) > 0; G_{X_c}(r_c) = 0
gives C, D > 0. ∎

*Remark.* This is Astra's Lemma 2 with the hypothesis "η small" replaced by the explicit condition (V). For three
components every η > 0 is valid. At η = 0, B_c = C_c = 0 (X_c is then the pinch of both sections: the doubly
pinched feed stage), which is why every non-binary column needs η_c > 0.

---

## 4. Binary columns: margin and forcing

Fix a binary column c = [p, p+1] with feed f = f_p e_p + f_{p+1} e_{p+1}, f_p, f_{p+1} > 0, net vapor U and vapor
V_T, V_B = V_T − U. Let θ' ∈ (α_{p+1}, α_p) be its feed root and

  Ψ⁰_T(t) := V_T − α_p f_p/(α_p − t),  Ψ⁰_B(t) := V_B + α_{p+1} f_{p+1}/(α_{p+1} − t),  κ := Ψ⁰_T(θ') = Ψ⁰_B(θ')     (4.1)

(the sharp section functions; Ψ⁰_T − Ψ⁰_B = Φ_c vanishes at θ'). κ = V_T − α_p f_p/(α_p − θ') is the vapor of the
column above the Underwood minimum computed with its actual feed: the **binary margin**.

**Proposition 4.1 (sharp binary column).** If κ > 0 then V_T > α_pf_p/(α_p − θ') > f_p, V_B > α_{p+1}f_{p+1}/
(θ' − α_{p+1}) > 0, and the roots φ⁰_T := α_p(V_T − f_p)/V_T of Ψ⁰_T and φ⁰_B := α_{p+1}(V_B + f_{p+1})/V_B of Ψ⁰_B
satisfy α_{p+1} < φ⁰_B < θ' < φ⁰_T < α_p. Moreover Ψ⁰_B(t) > κ for t > θ' and Ψ⁰_T(t) > κ for t < θ' (t outside
the poles).

*Proof.* V_T = κ + α_pf_p/(α_p − θ'), V_B = κ + α_{p+1}f_{p+1}/(θ' − α_{p+1}); α_p/(α_p − θ') > 1. Ψ⁰_T is
decreasing on (−∞, α_p) with Ψ⁰_T(θ') = κ > 0, Ψ⁰_B increasing on (α_{p+1}, ∞) with Ψ⁰_B(θ') = κ > 0; their roots
solve α_p − t = α_p f_p/V_T and t − α_{p+1} = α_{p+1} f_{p+1}/V_B. ∎

(For f_p, f_{p+1} the global split this is the classical binary condition x̂_T < x̂_B of `sections.md` Prop. 10:
κ > 0 ⇔ φ⁰_T > φ⁰_B ⇔ the top pinch lies below the bottom pinch; cf. `check_finishing_pinch.py`.)

**B-pinning** (π_p = B). For η > 0 put

  d = (f_p − η) e_p,  b = η e_p + f_{p+1} e_{p+1},  D = f_p − η,  L_T = V_T − f_p + η,  L_B = V_B + f_{p+1} + η.   (4.2)

**Proposition 4.2.** Let κ > 0 and 0 < η < η̂_p with

  η̂_p := f_p  (p = 1),   η̂_p := min(f_p, (α_{p−1} − α_p) Ψ⁰_B(α_{p−1})/α_p)  (p ≥ 2),                        (4.3)

(η̂_p > 0 since Ψ⁰_B(α_{p−1}) > κ). Then:
(a) all four phase flows and D, B are positive;
(b) T: w = d ≥ 0, the single root φ^T = α_p L_T/V_T ∈ [φ⁰_T, α_p), modes φ^T and α_k (k ≠ p); the trace modes are
> φ^T for k < p and < φ^T for k ≥ p+1;
(c) B: w = −b ≤ 0, Ψ_B(t) = Ψ⁰_B(t) + α_pη/(α_p − t) has exactly the roots φ'' ∈ (α_{p+1}, φ⁰_B) and
φ' ∈ (α_p, α_{p−1}) (φ' ∈ (α_1, ∞) if p = 1); modes: α_k (k < p) > φ' > φ'' > α_k (k > p+1);
(d) **margin**: Ψ_B(φ^T) > κ > 0; equivalently φ'' < φ^T;
(e) the forced pinch X := ê^B_{φ'} = (φ'/L_B)(η/(φ' − α_p) e_p + f_{p+1}/(φ' − α_{p+1}) e_{p+1}) is positive on
{p, p+1}, and X = A' ê^T_{(p+1)} + B' e_p with A', B' > 0, A' + B' = 1, where e_p = ê^T_{φ^T} and
ê^T_{(p+1)} = ê^T_{α_{p+1}} = (α_{p+1}/L_T)((f_p − η)/(α_p − α_{p+1}) e_p + (Ψ_T(α_{p+1})/α_{p+1}) e_{p+1}) > 0 on {p,p+1};
explicitly B' = (α_p − φ^T)G_X(φ^T)/α_p = (α_p − φ^T)φ'Ψ_B(φ^T)/(α_p L_B (φ' − φ^T)) > 0;
(f) as η → 0, φ' ↓ α_p, η/(φ' − α_p) → Ψ⁰_B(α_p)/α_p and X → ê^{B,0}_{(p)}, the pinch of the unforced bottom
section at its trace mode α_p, whose p-entry is Ψ⁰_B(α_p)/(V_B + f_{p+1}) ∈ (0, 1).

*Proof.* (a) D = f_p − η > 0, B = f_{p+1} + η, L_T = (V_T − f_p) + η > 0 (Prop. 4.1), L_B > 0. (b) Ψ_T(t) =
V_T − α_p(f_p − η)/(α_p − t) vanishes at α_p − t = α_p(f_p − η)/V_T, i.e. t = α_pL_T/V_T, which increases with η.
(c) Lemma 2.2(b) with poles α_p, α_{p+1}: one root in (α_{p+1}, α_p) and one in (α_p, ∞). On (α_{p+1}, α_p) the added
term α_pη/(α_p − t) is positive, so Ψ_B(φ⁰_B) > 0 and φ'' < φ⁰_B. φ' < α_{p−1} ⇔ Ψ_B(α_{p−1}) > 0 (Ψ_B increases on
(α_p, ∞)) ⇔ Ψ⁰_B(α_{p−1}) − α_pη/(α_{p−1} − α_p) > 0 ⇔ η < (α_{p−1} − α_p)Ψ⁰_B(α_{p−1})/α_p; and Ψ⁰_B(α_{p−1}) >
Ψ⁰_B(θ') = κ by Prop. 4.1. (d) φ^T ≥ φ⁰_T > θ' and φ^T < α_p, so Ψ_B(φ^T) = Ψ⁰_B(φ^T) + α_pη/(α_p − φ^T) >
Ψ⁰_B(φ^T) > κ. As Ψ_B increases on (α_{p+1}, α_p), Ψ_B(φ^T) > 0 ⇔ φ^T > φ''. (e) Positivity: φ' > α_p > α_{p+1}.
On the face {p, p+1}, T's modes are φ^T (vector e_p, since d = (f_p − η)e_p) and α_{p+1} (ℓ = e_{p+1}); Lemma 2.3
gives the entries of ê^T_{(p+1)}, positive because Ψ_T(α_{p+1}) > 0 (α_{p+1} < φ^T, Lemma 2.2(a)). The coefficient
of ê^T_{(p+1)} in X is X_{p+1}/(ê^T_{(p+1)})_{p+1} > 0 and that of e_p is G_X(φ^T)/G_{e_p}(φ^T) =
G_X(φ^T)(α_p − φ^T)/α_p; by Lemma 2.3(iii) applied to B, G_X(φ^T) = −(φ'/L_B)Ψ_B(φ^T)/(φ^T − φ') > 0 by (d) and
φ^T < φ'. (f) Lemma 2.5 with k = p: φ' − α_p = α_pη/Ψ^{(p)}_B(φ'), Ψ^{(p)}_B = Ψ⁰_B. ∎

**T-pinning** (π_p = T). For η > 0 put

  d = f_p e_p + η e_{p+1},  b = (f_{p+1} − η) e_{p+1},  D = f_p + η,  L_T = V_T − f_p − η,  L_B = V_B + f_{p+1} − η.   (4.4)

**Proposition 4.3.** Let κ > 0 and 0 < η < η̂_p with

  η̂_p := min(f_{p+1}, V_T − f_p)  (p+1 = N),   η̂_p := min(f_{p+1}, (α_{p+1} − α_{p+2}) Ψ⁰_T(α_{p+2})/α_{p+1})  (p+1 ≤ N−1)   (4.6)

(η̂_p > 0 since V_T > f_p and Ψ⁰_T(α_{p+2}) > κ). Then:
(a) all phase flows and D, B are positive (for p+1 ≤ N−1, L_T = Ψ_T(0) > Ψ_T(α_{p+2}) > 0);
(b) B: w = −b ≤ 0, single root φ^B = α_{p+1}L_B/V_B ∈ (α_{p+1}, φ⁰_B], modes φ^B and α_k (k ≠ p+1); trace modes
> φ^B for k ≤ p, < φ^B for k > p+1;
(c) T: Ψ_T(t) = Ψ⁰_T(t) − α_{p+1}η/(α_{p+1} − t) has exactly the roots φ^mid ∈ (φ⁰_T, α_p) and
φ^low ∈ (α_{p+2}, α_{p+1}) (φ^low ∈ (0, α_N) if p+1 = N); modes α_k (k < p) > φ^mid > φ^low > α_k (k > p+1);
(d) **margin**: Ψ_T(φ^B) > κ > 0; equivalently φ^B < φ^mid;
(e) the forced pinch X := ê^T_{φ^low} = (φ^low/L_T)(f_p/(α_p − φ^low) e_p + η/(α_{p+1} − φ^low) e_{p+1}) is
positive on {p, p+1}, and X = C' e_{p+1} + D' ê^B_{(p)} with C', D' > 0, C' + D' = 1, where e_{p+1} = ê^B_{φ^B} and
ê^B_{(p)} = (α_p/L_B)((Ψ_B(α_p)/α_p) e_p + (f_{p+1} − η)/(α_p − α_{p+1}) e_{p+1}) > 0 on {p,p+1};
(f) as η → 0, X → ê^{T,0}_{(p+1)}, the pinch of the unforced top section at its trace mode α_{p+1}.

*Proof.* Mirror image of Proposition 4.2. (a) B = f_{p+1} − η > 0; L_T > 0 by the bound (for p+1 ≤ N−1 via
Lemma 2.2(a): Ψ_T > 0 on (−∞, φ^low) ∋ 0, α_{p+2}). (b) φ^B solves t − α_{p+1} = α_{p+1}(f_{p+1} − η)/V_B; it
exceeds α_{p+1} because f_{p+1} > η, and decreases with η. (c) Lemma 2.2(a): one root in (α_{p+1}, α_p), one in
(0, α_{p+1}); on (α_{p+1}, α_p) the added term α_{p+1}η/(t − α_{p+1}) is positive, so Ψ_T(φ⁰_T) > 0 and
φ^mid > φ⁰_T; φ^low > α_{p+2} ⇔ Ψ_T(α_{p+2}) > 0 ⇔ η < (α_{p+1} − α_{p+2})Ψ⁰_T(α_{p+2})/α_{p+1}, and Ψ⁰_T(α_{p+2}) >
Ψ⁰_T(θ') = κ (Prop. 4.1). (d) φ^B ≤ φ⁰_B < θ', φ^B > α_{p+1}: Ψ_T(φ^B) > Ψ⁰_T(φ^B) > κ; Ψ_T decreases on
(α_{p+1}, α_p). (e) G_X(φ^B) = −(φ^low/L_T)Ψ_T(φ^B)/(φ^B − φ^low) < 0; the coefficient of e_{p+1} = ê^B_{φ^B} is
G_X(φ^B)/G_{e_{p+1}}(φ^B) = G_X(φ^B)(α_{p+1} − φ^B)/α_{p+1} > 0, that of ê^B_{(p)} is X_p/(ê^B_{(p)})_p > 0;
Ψ_B(α_p) > 0 since α_p > φ^B. (f) Lemma 2.5 with k = p+1. ∎

---

## 5. The specification F(δ, η, π) and the parameter region

**Definition 5.1 (specification).** Data (α, z, q), δ > 0, η = ((η_c)_{c non-binary}, (η_p)_{1≤p≤N−1}),
π ∈ {B, T}^{N−1}.
(F1) Process the non-binary columns by decreasing length m − i. For c = [i, m]: its feed (f_c, U_c) by (1.1) (for
[1,N]: z, 1 − q); its preferred split (Lemma 3.1) relative to (f_c, U_c); its deformation with η_c (Lemma 3.2,
Definition in §3): d_c, b_c, V_{T_c}, V_{B_c}, L_{T_c}, L_{B_c}. (The recursion is well defined as long as the
already processed columns are valid, Lemma 5.2(a).)
(F2) Binary columns [p, p+1]: feed and U by (1.1). Vapors, from the reboiler up:

  V_{B[N−1,N]} := R := H* + δ,  V_{T[p,p+1]} := V_{B[p,p+1]} + U_{[p,p+1]},  V_{B[p−1,p]} := V_{T[p,p+1]}.     (5.1)

Then d, b, D, L_T, L_B by (4.2) (π_p = B) or (4.4) (π_p = T) with η = η_p.
(F3) Products: P_1 := D_{[1,2]}, P_p := B_{[p−1,p]} + D_{[p,p+1]} (2 ≤ p ≤ N−1), P_N := B_{[N−1,N]}.
The **operating specification** handed to a finite-stage computation consists of the stream flows only: V, L of
every section (with the coupling streams equal to the section end flows) and P_1, …, P_N, plus the feed (1, z, q)
to J[1,N]. The net component flows d, b are not prescribed; they belong to the limit configuration (§6).

**Lemma 5.2 (balances).** Suppose every non-binary column is valid (V).
(a) Every feed f_c is > 0 on the range of c and 0 outside; (F1) is well defined.
(b) Every junction stage J_c and every side-draw stage P_p (2 ≤ p ≤ N−1) satisfies the total material and CMO vapor
balances with zero duty; the condenser removes V_{T[1,2]}; the reboiler receives R = H* + δ. Hence R is the total
heat supplied (class.md §3), for every number of stages per section.
(c) V_{T[p,p+1]} = R + Σ_{r=p}^{N−1} U_{[r,r+1]} and V_{B[p,p+1]} = R + Σ_{r=p+1}^{N−1} U_{[r,r+1]}; the feeds and
U's of the binary columns depend only on the non-binary data (and not on δ, η_bin, π).
(d) Σ_p P_p = 1, and the limit net flows satisfy every component balance: d_c + b_c = f_c, P_p x_{P_p} = n_p :=
b_{[p−1,p]} + d_{[p,p+1]}, Σ_p n_p = z.
(e) (*Flux bookkeeping.*) For t ∉ {α_i}, let Ψ_σ(t) be computed from the limit net flows and φ_{P}(t) :=
−Σ_k α_k n_{p,k}/(α_k − t) be the Underwood flux of the liquid product p. At every unit u (junction, side-draw stage,
condenser, reboiler) Σ_{σ: bottom end u} Ψ_σ(t) − Σ_{σ: top end u} Ψ_σ(t) + φ_{P}(t)[u = P_p] = Q_u + φ_F(t)[u = J[1,N]],
with Q_u the duty (0, −V_{T[1,2]} at the condenser, R at the reboiler) and φ_F(t) = 1 − q − Σ_k α_k z_k/(α_k − t).

*Proof.* (a) Induction in the order of (F1): f_{[i,m]} = d_{[i,m+1]} + b_{[i−1,m]}, where d_{[i,m+1]} (if present) is
> 0 exactly on [i, m] and b_{[i−1,m]} (if present) is > 0 exactly on [i, m], by validity of the (non-binary) parents;
a binary column is never the parent of a column. (b) At J[i,m] the vapor outflows are V_{T[i,m]}
(up) and V_{B[i−1,m]} (back to the parent's bottom section), the inflows V_{B[i,m]} and V_{T[i,m+1]}; their
difference is V_{T_c} − V_{B_c} − U_c = 0 (for [1,N]: V_T − V_B − (1 − q) = 0, the feed bringing vapor 1 − q).
The liquid balance: inflow L_{T_c} + L_{B[i−1,m]} (+ q for the feed) equals outflow L_{B_c} + L_{T[i,m+1]}, since
L_{B_c} − L_{T_c} = D_c + B_c − U_c = Σf_c − U_c and L_{B[i−1,m]} − L_{T[i,m+1]} = B_{[i−1,m]} + D_{[i,m+1]} − U_c.
At P_p (side draw): vapor V_{B[p−1,p]} out, V_{T[p,p+1]} in, equal by (5.1); liquid L_{B[p−1,p]} in,
L_{T[p,p+1]} + P_p out: L_{B[p−1,p]} − L_{T[p,p+1]} = B_{[p−1,p]} + D_{[p,p+1]} = P_p. Condenser: V_{T[1,2]} in,
L_{T[1,2]} + P_1 out. Reboiler: L_{B[N−1,N]} in, V_{B[N−1,N]} = R + P_N out. Stages inside sections carry no duty
(CMO). (c) Telescoping (5.1). (d) Summation of the component balances d_c + b_c = f_c over the tree of columns.
(e) Insert the definitions: at J_c, Ψ_{T_c} + Ψ_{B[i−1,m]} − Ψ_{B_c} − Ψ_{T[i,m+1]} = (V_{T_c} − V_{B_c} − U_c) −
Σ_k α_k(d_{c,k} + b_{c,k} − f_{c,k})/(α_k − t) = 0 by (1.1) (for J[1,N]: = (1 − q) − Σ_kα_kz_k/(α_k − t) = φ_F); at a
side draw, Ψ_{B[p−1,p]} − Ψ_{T[p,p+1]} + φ_{P_p} = (V_{B[p−1,p]} − V_{T[p,p+1]}) + Σ_kα_k(b_{[p−1,p],k} + d_{[p,p+1],k} −
n_{p,k})/(α_k − t) = 0; at the condenser, −Ψ_{T[1,2]} + φ_{P_1} = −V_{T[1,2]}; at the reboiler, Ψ_{B[N−1,N]} + φ_{P_N} =
V_{B[N−1,N]} = R. (This is class.md Lemma 5.1 summed over the stages of each section.) ∎

**Lemma 5.3 (the middle roots are the global feed roots).** If every non-binary column is valid, then for every
non-binary c = [i, m] and every r with i+1 ≤ r ≤ m−2, θ^c_r = θ_r. Its extreme feed roots θ^c_i, θ^c_{m−1} are
real-analytic in η and equal θ_i, θ_{m−1} at η = 0.

*Proof.* Induction in the order of (F1). For [1,N], Φ = 1 − q − Σα z/(α − t), whose roots are the θ_r. For
c = [i,m] ≠ [1,N], Φ_c = Ψ_{T[i,m+1]} − Ψ_{B[i−1,m]} by (1.2). If [i,m+1] is present, it is non-binary, and by
Proposition 3.3(a), Lemma 3.2 and induction Ψ_{T[i,m+1]} = Ψ⁰ + ηK vanishes at its middle roots θ_{i+1}, …, θ_{m−1}; if [i−1,m] is
present, Ψ_{B[i−1,m]} vanishes at θ_i, …, θ_{m−2}. Both families contain θ_{i+1}, …, θ_{m−2}, so Φ_c vanishes there,
and since Φ_c has a single root in each (α_{r+1}, α_r), θ^c_r = θ_r. At η = 0 (no deformation) Ψ_{T_c} = Ψ⁰_{T_c}
and Ψ_{B_c} = Ψ⁰_{B_c} vanish at *all* feed roots of c (Lemma 3.1), and the same induction shows that every feed root
of every column, binary columns included, is global: θ^c_r = θ_r for i ≤ r ≤ m−1 (HS's carry-over, paper.tex proof
of Theorem `achieve`). Analyticity: simple roots of Φ_c, whose coefficients are analytic in η. ∎

**Proposition 5.4 (the binary margins).** Write η_nb = (η_c)_{c non-binary} and κ_p(η_nb; δ) for the margin (4.1)
of the binary column [p, p+1] in F.
(a) The non-binary part of F and the feeds and U's of the binary columns do not depend on (δ, η_bin, π).
(b) κ_p(η_nb; δ) = κ_p(η_nb; 0) + δ.
(c) At η_nb = 0 the non-binary columns are at HS's preferred-split point (paper.tex Theorem `achieve`), the feed root
of [p,p+1] is θ_p, and κ_p(0; 0) = e_p = H* − H_p.
(d) There are η_1 > 0 and C_1 ≥ 0, depending only on (α, z, q), such that for every η_nb ∈ [0, η_1]^{n_nb}
(n_nb = (N−1)(N−2)/2) every non-binary column is valid and |κ_p(η_nb; 0) − e_p| ≤ C_1 max_c η_c for every p.

*Proof.* (a) Definition 5.1: (F1) does not involve R, η_bin or π. (b) By Lemma 5.2(c), V_{T[p,p+1]} = R + const with
const and the binary feed independent of R; κ_p = V_{T[p,p+1]} − α_pf_p/(α_p − θ'). (c) With η_nb = 0, (F1) is
exactly the recursive construction of Theorem `achieve` (every non-binary column at the preferred split of its actual
feed), and by Lemma 5.3 (case η = 0) every feed root of every column is global: θ_p is a common root of Ψ_{T[i,p+1]}
and Ψ_{B[i,p+1]} for every non-binary [i, p+1] (Lemma 3.1), and θ_p is the feed root of [p, p+1]. κ_p does not
depend on η_bin, so we may use the sharp binary flows (η_bin = 0), for which every product is pure with flow z_r:
n_r = z_r e_r (Lemma 5.2(d); n_r is supported on {r}). We sum the identities of Lemma 5.2(e) at t = θ_p over the
*light* units: J[i,m] with m ≤ p, the condenser P_1 and the side draws P_2, …, P_p. A section with both ends among
them contributes +Ψ − Ψ = 0. The sections with exactly one end among them are T[i, p+1], 1 ≤ i ≤ p (top end J[i,p]
or P_p light, bottom end J[i,p+1] not): a bottom section B[i,m] with an end among the light units has m ≤ p (its
bottom end J[i+1,m] or P_m), hence both ends light, and a top section T[i,m] with a light end has m ≤ p + 1. The
feed enters J[1,N], not a light unit (p ≤ N − 1). Hence

  −Σ_{i=1}^p Ψ_{T[i,p+1]}(θ_p) + Σ_{r=1}^p φ_{P_r}(θ_p) = −V_{T[1,2]}.

Summing (e) over *all* units, every section cancels and Σ_r φ_{P_r}(t) = R − V_{T[1,2]} + φ_F(t) for every t;
letting t → ∞ (φ_P → 0, φ_F → 1 − q) gives V_{T[1,2]} = R + 1 − q. With Ψ_{T[i,p+1]}(θ_p) = 0 for i < p and
φ_{P_r}(θ_p) = −α_rz_r/(α_r − θ_p):

  κ_p(0; δ) = Ψ_{T[p,p+1]}(θ_p) = R + (1 − q) − Σ_{r≤p} α_r z_r/(α_r − θ_p) = R − H_p = e_p + δ,

by the feed equation Σ_r α_rz_r/(α_r − θ_p) = 1 − q. (This is the light-side cut identity of paper.tex, Proposition
`cut`, and re-proves the side-draw balances of Theorem `achieve`: at δ = 0 the binary top vapor is
α_pf_p/(α_p − θ_p) + e_p.) (d) By induction in the order of (F1), all data of every column
(feed, U, feed roots, preferred split, h, deformed split and flows) are real-analytic functions of η_nb on a
neighborhood of 0 in ℝ^{n_nb} as long as the feeds remain positive: simple roots of Φ_c depend analytically on
(f_c, U_c) (implicit function theorem), the other quantities are rational in them with nonvanishing denominators.
At η_nb = 0 every inequality in (V) is strict (Lemma 3.1) and the binary feeds are positive; so there is η_1 > 0 such
that on the compact cube Q = [0, η_1]^{n_nb} all columns are valid, all feeds positive, and κ_p(·; 0) is C¹ (its
binary feed root is simple). Put C_1 := max_p max_Q ‖∇κ_p(·; 0)‖_1 and use the mean value theorem along the segment
from 0 to η_nb, with (c). ∎

**Definition (admissible parameters).** (δ, η, π) with δ > 0 and all η's > 0 is **admissible** if
(i) every non-binary column is valid (V) along the recursion (F1);
(ii) κ_p > 0 for every binary column (margin (4.1) computed from F);
(iii) 0 < η_p < η̂_p for every binary column (bounds (4.3), (4.6), computed from F; positive by (ii)).
These are finitely many explicit inequalities on the flows of F, checkable by evaluating F; Propositions 3.3 (for every
non-binary column) and 4.2 (π_p = B) or 4.3 (π_p = T) (for every binary column) then apply. Only (i)–(iii) are used
in §6.

**Theorem 5.5 (the admissible set is nonempty for every δ, with η̄(δ) linear in δ).** Let η_1, C_1 be as in
Proposition 5.4, δ > 0, η̄(δ) := min(η_1, δ/(2C_1)). If 0 < η_c ≤ η̄(δ) for every non-binary column, then (i) holds
and κ_p ≥ e_p + δ/2 ≥ δ/2 > 0 for every binary column, so (ii) holds; choosing then any η_p ∈ (0, η̂_p) gives an
admissible parameter.

*Proof.* Proposition 5.4(b), (d): κ_p(η_nb; δ) = κ_p(η_nb; 0) + δ ≥ e_p − C_1η̄(δ) + δ ≥ e_p + δ/2. ∎

The admissible set is open in the η variables; every quantity below is a continuous (indeed real-analytic) function
on it. Numerically (§8.3) the largest admissible common value of η_nb is ≈ (0.12–0.26)·δ for the test data, in
agreement with the linear bound.

---

## 6. The limit configuration and hypotheses (H1)–(H8)

**Definition 6.1 (limit configuration 𝓛(δ, η, π)).** For admissible parameters:
- *net flows*: d_c, b_c as in F (non-binary: (F1); binary: (4.2) or (4.4)); the net flow of every component outside
  a column's range is 0 in both of its sections;
- *junctions*: for non-binary c, X_c of Proposition 3.3(e); for binary c, X_c := ê^{B_c}_{φ'} (π_p = B,
  Proposition 4.2(e)) or X_c := ê^{T_c}_{φ^low} (π_p = T, Proposition 4.3(e));
- *products*: x_{P_p} := n_p/P_p, n_p = b_{[p−1,p]} + d_{[p,p+1]} (absent terms omitted).

Node supports. supp X_c = [i, m] for c = [i, m] (binary: {p, p+1}); supp x_{P_p} = Maj(p) :=
{p} ∪ {p−1 : π_{p−1} = B} ∪ {p+1 : π_p = T} (with π_0 := T, π_N := B, i.e. P_1 ⊆ {1, 2}, P_N ⊆ {N−1, N}). All
entries on the supports are positive.

**Table 1 (modes and pinch of every section; "above"/"below" the pinch r*).**

| type | section | modes (decreasing) | r* | closed at bottom (above r*) | closed at top (below r*) |
|---|---|---|---|---|---|
| 1 | T_c, c = [i,m] non-binary | α_k (k<i) > r_c > θ_{i+1} > ⋯ > θ_{m−2} > u_c > α_m > α_k (k>m) | u_c | α_k (k<i), r_c, θ's | α_m, α_k (k>m) |
| 2 | B_c, c = [i,m] non-binary | α_k (k<i) > α_i > v_c > θ's > s_c > α_k (k>m) | v_c | α_k (k≤i) | θ's, s_c, α_k (k>m) |
| 3 | T_c, binary B-pinned | α_k (k<p) > φ^T > α_k (k≥p+1) | φ^T | α_k (k<p) | α_k (k≥p+1) |
| 4 | B_c, binary B-pinned | α_k (k<p) > φ' > φ'' > α_k (k>p+1) | φ' | α_k (k<p) | φ'', α_k (k>p+1) |
| 5 | T_c, binary T-pinned | α_k (k<p) > φ^mid > φ^low > α_k (k>p+1) | φ^low | α_k (k<p), φ^mid | α_k (k>p+1) |
| 6 | B_c, binary T-pinned | α_k (k≤p) > φ^B > α_k (k>p+1) | φ^B | α_k (k≤p) | α_k (k>p+1) |

(Ends: types 1, 3, 5 have bottom X_c and top J[i,m−1] resp. P_p; types 2, 4, 6 have bottom J[i+1,m] resp. P_{p+1}
and top X_c.) In every top section the pinch is the **smallest** major root, in every bottom section the
**largest**; every trace mode is above all major modes or below all of them.

**Theorem 6.2.** For admissible parameters the limit configuration has the following properties.

(H1) *Flows.* Every V_σ, L_σ, D_c, B_c, P_p is positive; every feed is positive on its column's range.

(H2) *Net flows.* Top sections have w ≥ 0 and bottom sections w ≤ 0 componentwise. Non-binary c = [i, m]:
supp d_c = [i, m−1], d_{c,i} = f_{c,i}, supp b_c = [i+1, m], b_{c,m} = f_{c,m}; binary B-pinned: d = (f_p − η_p)e_p,
b = η_p e_p + f_{p+1}e_{p+1}; T-pinned: d = f_pe_p + η_pe_{p+1}, b = (f_{p+1} − η_p)e_{p+1}; the forced flows are
never zero ("never balanced"). Middle roots of every non-binary column are the global θ's (Lemma 5.3).

(H3) *Non-binary columns*: all conclusions of Proposition 3.3 (interlacing r_c > θ^c_i > v_c, u_c > θ^c_{m−1} > s_c,
common simple middle roots, Ψ_{T_c}(α_m) > 0, Ψ_{B_c}(α_i) > 0, positive pinch vectors, A_c, B_c, C_c, D_c > 0).

(H4) *Binary columns*: κ_p > 0 (≥ e_p + δ/2 for the choice of Theorem 5.5), and all conclusions of Proposition 4.2 or 4.3 (forced root located, margin
Ψ_B(φ^T) > κ_p resp. Ψ_T(φ^B) > κ_p, positive cone decomposition of X_c in the non-pinning section's basis).

(H5) *Spectra.* Table 1 holds; in every section all N modes are positive and pairwise distinct, the pinch is a major
mode, ρ_σ := max over consecutive modes of (smaller/larger) < 1.

*Limit relations.* For every section σ with bottom-end composition x_b and top-end composition x_t:
ℓ_φ·x_b = 0 for every mode φ above r*_σ, ℓ_φ·x_t = 0 for every mode φ below r*_σ, and the pivots
c_{r*}(x_b) = ℓ_{r*}·x_b/ℓ_{r*}·ê_{r*} and c_{r*}(x_t) are > 0.

(H6) *Half-profiles.* For every section, the sequences x^b_k := M_σ^k x_b/1·(M_σ^k x_b) and
x^t_k := M_σ^{−k} x_t/1·(M_σ^{−k} x_t), k ≥ 0, are well defined, lie in the simplex, and converge to ê_{r*}
geometrically (ratio ≤ ρ_σ); the chart vectors v^b_k := M^k x_b/((Vr*)^k c_{r*}(x_b)) = ê_{r*} + Σ_{φ<r*}
(c_φ(x_b)/c_{r*}(x_b))(φ/r*)^k ê_φ and v^t_k (same with φ > r*, (r*/φ)^k) satisfy inf_k 1·v_k > 0.

(H7) *Traces.* For every section σ and every trace mode α_k of σ, closed at the end J (Table 1):
(i) Ψ^{(k)}_σ(α_k) = Ψ_σ(α_k) > 0; (ii) k is lighter than every component of supp x_J if J is the bottom end,
heavier if J is the top end; hence G^{(k)}_{x_J}(α_k) < 0 resp. > 0; (iii) x_{J,k} = 0. Consequently the linearized
closed-end condition (Lemma 2.5) reads x_{J,k} = κ_σ^{(k)} q_σ with q_σ the net flow of k from J into σ and
κ_σ^{(k)} := |G^{(k)}_{x_J}(α_k)|/Ψ^{(k)}_σ(α_k) > 0. Moreover (iv) every junction J_c at which k is a trace is the
closed end, for k, of T_c (k lighter than the range of c) or of B_c (k heavier); (v) every component with nonzero
net flow in σ belongs to the supports of both end compositions of σ.

(H8) *Duty, purities, pivots.* R = H* + δ is the only heat input (Lemma 5.2). The impurity of product p is
1 − x_{P_p,p} = (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p. The pivots of (H5) are > 0 and, in Underwood form:
at a product end G_{x_P}(r*) > 0 for a top section, < 0 for a bottom section; at a child junction J[i,m],
G_{X[i,m]}(u_{[i,m+1]}) > 0 and G_{X[i,m]}(v_{[i−1,m]}) < 0; at a non-binary junction the pivots are B_c and C_c;
at a binary junction 1 for the pinning section and B' resp. C' for the other. All quantities in (H1)–(H8) are
continuous in the parameters on the admissible set, hence bounded away from 0 on its compact subsets.

*Proof.* (H1)–(H4) collect Lemma 5.2, Propositions 3.3, 4.2, 4.3 and Theorem 5.5.

*(H5) and the limit relations, type by type.* The mode lists come from Proposition 3.3(c) (types 1, 2) and
Propositions 4.2(b, c), 4.3(b, c) (types 3–6); distinctness and positivity are stated there, so ρ_σ < 1. The pinch is
the indicated major mode by definition; we check that the closed modes annihilate the corresponding ends, using
ℓ_φ·x = G_x(φ) for a major mode and ℓ_φ·x = x_k for a trace mode α_k (Lemma 2.3).
1. T_c non-binary: at x_b = X_c, G_{X_c}(r_c) = G_{X_c}(θ_r) = 0 (definition of X_c) and X_{c,k} = 0 for k < i; at
   x_t = X_{[i,m−1]}, supported on [i, m−1], x_{t,m} = x_{t,k} = 0 (k > m).
2. B_c non-binary: at x_b = X_{[i+1,m]}, supported on [i+1, m], x_{b,k} = 0 for k ≤ i; at x_t = X_c,
   G_{X_c}(θ_r) = G_{X_c}(s_c) = 0 and X_{c,k} = 0 for k > m.
3. T_c binary B-pinned: X_c is supported on {p, p+1}, so X_{c,k} = 0 for k < p; x_t = x_{P_p} with
   Maj(p) ⊆ {p−1, p} (π_p = B), so x_{t,k} = 0 for k ≥ p+1.
4. B_c binary B-pinned: x_b = x_{P_{p+1}} with Maj(p+1) ⊆ {p, p+1, p+2}, so x_{b,k} = 0 for k < p; x_t = X_c =
   ê_{φ'} is annihilated by ℓ_{φ''} (Lemma 2.3(iv)) and X_{c,k} = 0 for k > p+1.
5. T_c binary T-pinned: x_b = X_c = ê_{φ^low} is annihilated by ℓ_{φ^mid}, X_{c,k} = 0 for k < p; x_t = x_{P_p} with
   Maj(p) ⊆ {p−1, p, p+1}, x_{t,k} = 0 for k > p+1.
6. B_c binary T-pinned: x_b = x_{P_{p+1}} with Maj(p+1) ⊆ {p+1, p+2} (π_p = T), x_{b,k} = 0 for k ≤ p; X_{c,k} = 0
   for k > p+1.
Pivots. At a junction X_c the pivot of the section whose pinch vector appears in the cone decomposition is that
decomposition's coefficient: B_c (type 1), C_c (type 2), B' (type 3), 1 (type 4), 1 (type 5), C' (type 6), all
> 0. At the other end x (a product, or a child's junction for types 1–2), supp x lies entirely on one side of
r*: for a top section every component of supp x is lighter than r* (type 1: supp ⊆ [i, m−1], u_c < α_{m−1};
type 3: ⊆ {p−1, p}, φ^T < α_p; type 5: ⊆ {p−1, p, p+1}, φ^low < α_{p+1}), so G_x(r*) > 0, and
ℓ_{r*}·ê_{r*} = (r*/L)(−Ψ'(r*)) > 0 (Lemma 2.3(iv), Ψ decreasing); for a bottom section every component of supp x
is heavier than r* (type 2: ⊆ [i+1, m], v_c > α_{i+1}; type 4: ⊆ {p, p+1, p+2}, φ' > α_p; type 6: ⊆ {p+1, p+2},
φ^B > α_{p+1}), so G_x(r*) < 0 and ℓ_{r*}·ê_{r*} < 0. In both cases c_{r*}(x) > 0. This also proves the Underwood
form of the pivots in (H8) (for types 1–2 the other end is the child's junction).

*(H6).* By the limit relations x_b = c_{r*}(x_b)(ê_{r*} + Σ_{φ<r*} a_φ ê_φ) and x_t = c_{r*}(x_t)(ê_{r*} +
Σ_{φ>r*} b_φ ê_φ), so M^k x_b = c_{r*}(x_b)(Vr*)^k v^b_k and M^{−k}x_t = c_{r*}(x_t)(Vr*)^{−k} v^t_k with the
stated chart vectors, which converge to ê_{r*} at the rate ρ_σ^k. Nonnegativity: going *up* a bottom section (types
2, 4, 6) M ≥ 0 (Lemma 2.4(a)); going *down* a top section (types 1, 3, 5) M^{−1} ≥ 0 (Lemma 2.4(b), L > 0). The two
remaining directions start at a junction X_c and use the cone decompositions with nonnegative pinch vectors (Lemma
2.4(c)): going up T_c from X_c: type 1, X_c = A_cê^T_{(m)} + B_cê^T_u; type 3, X_c = A'ê^T_{(p+1)} + B'e_p; type 5,
X_c = ê_{φ^low} (constant). Going down B_c from X_c: type 2, X_c = C_cê^B_v + D_cê^B_{(i)}; type 4, X_c = ê_{φ'}
(constant); type 6, X_c = C'e_{p+1} + D'ê^B_{(p)}. Hence M^{±k}x ≥ 0, ≠ 0, so 1·v_k > 0; since 1·v_k → 1·ê_{r*} = 1,
inf_k 1·v_k > 0. The normalized sequences are therefore in the simplex.

*(H7).* (iii) is part of the limit relations. (ii): the closed-at-bottom trace modes of Table 1 are α_k with k < i
(type 1; supp X_c = [i, m]), k ≤ i (type 2; supp = [i+1, m]), k < p (types 3, 5; supp X_c = {p, p+1}), k < p (type 4;
supp x_{P_{p+1}} ⊆ {p, p+1, p+2}), k ≤ p (type 6; supp ⊆ {p+1, p+2}); the closed-at-top ones are k ≥ m (type 1;
supp ⊆ [i, m−1]), k > m (type 2; supp X_c = [i, m]), k ≥ p+1 (type 3; supp ⊆ {p−1, p}), k > p+1 (types 4, 6;
supp X_c = {p, p+1}), k > p+1 (type 5; supp ⊆ {p−1, p, p+1}). In each case k is lighter (bottom) resp. heavier (top)
than the whole support, so every term of G^{(k)}_{x_J}(α_k) = Σ_{i∈supp} α_ix_i/(α_i − α_k) is negative resp.
positive. (i): a closed-at-bottom trace of a top section is lighter than every component with nonzero net flow in it
(those are in supp x_b by (v)), so α_k exceeds every pole and Ψ_σ(α_k) > V > 0 (Lemma 2.2(a)); a closed-at-top trace
of a top section lies below its smallest root, where Ψ_σ > 0 (Lemma 2.2(a)); a closed-at-bottom trace of a bottom
section lies above its largest root, where Ψ_σ > 0, and a closed-at-top trace of a bottom section lies below every
pole, where Ψ_σ > V (Lemma 2.2(b)). The identity x_{J,k} = κ q_σ: at w_k = 0 the linearization of ℓ̃_k·x_J = 0 is
dx_k + (G^{(k)}/Ψ^{(k)})dw_k = 0 (Lemma 2.5); q_σ = w_k at a bottom end and −w_k at a top end, and the signs in (i),
(ii) give dx_k = κ q_σ with κ > 0 in all four cases. (iv): for a junction X_c, c = [i, m], and k ∉ [i, m]: if k < i,
α_k is a trace mode of T_c above its pinch (types 1, 3, 5), closed at X_c; if k > m, α_k is a trace mode of B_c below
its pinch (types 2, 4, 6), closed at X_c. (v): the components with nonzero net flow are, by type: [i, m−1] (1), [i+1, m]
(2), {p} (3), {p, p+1} (4, 5), {p+1} (6); compare with the supports of the end compositions listed above.

*(H8).* R: Lemma 5.2(b). Purities: n_p = b_{[p−1,p]} + d_{[p,p+1]}, whose entries other than p are η_{p−1} (on p − 1,
if π_{p−1} = B) and η_p (on p + 1, if π_p = T). Pivots: shown above. Continuity: all quantities are built from simple
roots and rational expressions with nonvanishing denominators. ∎

---

## 7. Data for the limit system (Agent L) and the analysis (Agent A)

1. **Majors and traces at the nodes.** X_c: majors [i, m]; P_p: majors Maj(p); every other component is a trace
   (limit composition 0). In each section the components with nonzero limit net flow are listed in (H7)(v); all others
   are traces with zero limit net flow.
2. **Per column, the major conditions at X_c** (limit relations at the junction): non-binary c: G_{X}(r_c) = 0,
   G_X(θ^T_r) = G_X(θ^B_r) = 0 (i+1 ≤ r ≤ m−2), G_X(s_c) = 0 — 2(m − i) − 2 conditions for the m − i unknowns of X_c
   and the m − i − 2 free split coordinates (d_{i+1}, …, d_{m−1} with Σd fixed); binary B-pinned: G_X(φ'') = 0;
   binary T-pinned: G_X(φ^mid) = 0 (one condition, one unknown ratio). Parents impose only trace conditions on a
   child's junction (types 1, 2 at their open-pinch ends: Table 1).
3. **Trace conditions**: one per (section, trace mode), at the closed end J of Table 1, with κ_σ^{(k)} of (H7).
   Explicitly κ_σ^{(k)} = Σ_{i∈supp x_J} α_i x_{J,i}/|α_i − α_k| / Ψ_σ(α_k).
4. **Nondegeneracy constants** for fixed admissible (δ, η): ρ_σ < 1; pivots (H8); cone weights A_c, B_c, C_c, D_c,
   A', B', C', D'; κ_p; |r_c − v_c|, |u_c − s_c|, φ^T − φ'', φ^mid − φ^B. They degenerate as η → 0 (B_c, C_c = O(η_c);
   the pivot of a forced section at its product end is O(η_p)), so the limit n → ∞ must be taken before η → 0.
5. **Exports.** `scripts/export_config.py` writes F and 𝓛 (all of the above, modes, eigenvectors, pivots, κ's) as
   JSON, format in its header; examples in `runs/config_*.json`.

---

## 8. Numerical checks (evidence only)

All in `scripts/` (Python 3 + mpmath, 50–60 digits; WSL Fedora sci environment; logs in `runs/`). From the project
root (Git Bash):

```
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/design/scripts && python check_design.py 12 50 1"   # ~15 min, ALL PASS
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/design/scripts && python crosscheck_N3.py"        # < 1 s, PASS
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/design/scripts && python eta_scaling.py"          # ~1 min, ALL PASS
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/design/scripts && python export_config.py n4_1"    # seconds, JSON + checks
```

- `check_design.py [n_per_N] [dps] [seed]`: random data N = 3..7 (log-uniform volatility ratios in [1.25, 4];
  kinds: generic, one nearly equal pair of volatilities (ratio 1.03–1.1), one tiny feed fraction (10⁻³), q = −1,
  q = 2, two equal highest peaks (q solved for)), random forcing patterns (and all-B, and (B,T) for N = 3), relative
  margins δ/H* ∈ {0.05, 0.001}, binary forcing at 25% and 90% of the explicit bound η̂_p; η_nb found by halving from
  0.1δ until (V) and κ_p ≥ δ/2 hold. Checks every family of (H1)–(H8), the limit relations, the eigen-residuals of
  M, Lemma 5.3 (middle roots global), Proposition 5.4(b, c) (κ_p(δ) − κ_p(0) = δ; κ_p(0; δ) = e_p + δ), the
  bookkeeping identity (1.2), and the stage-level bound R ≥ Σ_P max(0, φ_P(θ_j)) at the limit configuration.
  Result: see `runs/check_design_12_50_1.log` (summary in §8.1 below).
- `crosscheck_N3.py`: Agent N's finite-stage runs (network.md §11.3: α = (4, 2, 1), z = 1/3, q = 1, δ = 0.05,
  prefractionator slack 0.025, η_A = η_C = 0.003, i.e. F(0.05, (0.025; 0.003, 0.003), BT)) converge to junction
  ratios ξ → 0.696154, ζ → 0.553931 and side-product impurity 0.01768; the limit configuration of Definition 6.1
  gives ξ = 0.6961542454, ζ = 0.5539306858, 0.017681729: **agreement to all printed digits** (PASS). The script
  also rebuilds N's total flows from N's own code (HS point + circulation + forcing) and checks that they coincide
  with F(0.05, (0.025; 0.003, 0.003), BT) (max deviation 9·10⁻¹⁷, N's data being binary floats). This is independent
  evidence that 𝓛 is the limit of the finite-stage steady states at these flows.
- `export_config.py`: JSON export; every export re-runs all checks (ALL PASS for the six exported cases).

### 8.1 Summary of the main run
`python check_design.py 12 50 1` (50 digits, 12 data sets per N for N = 3..7, two δ, two forcing levels: 240
instances; 15 min, one thread): **ALL PASS** — e.g. H1 balances 7440/7440, H3 cone weights 1680/1680, H4 margins
960/960, H5 limit relations 5280/5280, H6 half-profiles in the simplex 5280/5280, H7 trace signs 17760/17760,
H8 pivots 5280/5280, κ_p(δ) − κ_p(0) = δ and κ_p(0; δ) = e_p + δ 240/240, stage-level bound at the limit 960/960.
Extremes over the run: min κ_p/δ = 0.554 (> 1/2 as enforced), max ρ_σ = 0.9988 (a nearly-equal pair of
volatilities, ratio 1.03–1.1), smallest η_nb/δ needed = 2·10⁻⁴ (q = 2, where some non-binary column has a small
d⁰_p or b⁰_p and the explicit validity bound η_0(c) is small). A second run with the final version of the scripts
and another seed (`check_design.py 4 50 2`, 80 instances, `runs/check_design_4_50_2.log`): ALL PASS.

### 8.2 Independent confirmation (Agent C, `../computation/runs/compare_D.log`)
Agent C recomputed twelve of these exports (N = 3..6, several patterns, q < 0, wide α) from the data and parameters
only, with its own code: flows, limit net flows, node compositions, modes, pivots, κ's and the cone weights agree to
≤ 1.2·10⁻⁵⁴ relative; the limit points solve C's limit system to ≤ 1.2·10⁻⁵⁴ and its Jacobian is invertible there
(σ_min between 8·10⁻⁷ and 10⁻²). C's finite-stage solver converges to these limit configurations as n grows
(`../computation/PROGRESS.md`).

### 8.3 Linear scaling of η̄(δ) (`scripts/eta_scaling.py`, `runs/eta_scaling.log`)
Six data sets (N = 3..6, q ∈ {0, 0.5, 0.6, 0.8, 1}): κ_p(0; 0) − e_p ≤ 4·10⁻⁴⁰ (Prop. 5.4(c));
max_p |κ_p(η·1; 0) − e_p|/η converges as η → 0 to 1.9–4.2 (a directional derivative, so C_1 ≳ 2–4 here); the
largest common η_nb with (i) and κ_p ≥ e_p + δ/2 is (0.12–0.26)·δ for δ = 10⁻², …, 10⁻⁵, constant in δ (PASS).
The explicit validity bounds η_0 = min_c η_0(c) are 0.02–0.19, so for small δ the binding constraint is the margin.

---

## 9. Scope, corrections to the plan, open points

- **Proved here:** Theorem D (§§2–6): the design, its admissible region, and every property (H1)–(H8) of its limit
  configuration, for every N ≥ 3, every α, z, q and every forcing pattern. The proofs are self-contained: they use
  only Underwood algebra (re-proved in §2–3), Astra's interpolation lemma (re-proved in Prop. 3.3(e)) and the flux
  bookkeeping of Lemma 5.2(e); in particular κ_p(0; δ) = e_p + δ is proved directly by a light-side cut
  (Prop. 5.4(c)), which also re-proves the side-draw balances of paper.tex Theorem `achieve`. Not proved here: that finite-stage steady
  states converge to 𝓛 (the limit system's invertibility is Agent L's; the perturbation argument Agent A's). The
  N = 3 cross-check of §8 is numerical evidence for it.
- **Corrections / precisions to BRIEF-common §3, Step 1:** (1) no hierarchy "parents ≪ children" is needed; the only
  smallness is η_c ≤ η̄(δ) = min(η_1, δ/(2C_1)) for the non-binary columns (Theorem 5.5); (2) binary forcing needs no
  smallness relative to the margin κ_p (it increases the margin), only the explicit bounds (4.3), (4.6); (3) for
  T-pinning of the last binary column, L_T = V_T − f_p − η_p must stay positive (explicit in (4.6)); (4) the margin δ
  must enter through R only; every non-binary column gets exactly Astra's direction relative to its *actual* feed
  (equal increments η_c of V_T, L_T, V_B, L_B over the preferred split of that feed, with limit split d⁰ + η_ch).
  A different circulation through a column with ≥ 4 components is a different total-flow specification, for which
  none of the above is claimed (cf. Agent L's caution). Agent N's §10.3/§11.3 flows for N = 3 are in our family
  (they add δ/2 to the four phase flows of the prefractionator: η_{[1,3]} = δ/2, h = 0); for N ≥ 4 we have not
  checked whether N's circulation is; (5) any forcing pattern works; "all B" is the simplest uniform choice
  (product 1 pure, every other product receives one designed impurity).
- **Constants.** η_1, C_1 come from a compactness/analyticity argument and are not explicit; η̄(δ) is linear in δ.
  The explicit validity bound η_0(c) and the binary bounds η̂_p are computable from F. Numerically η_nb = 0.1δ already
  works in most random instances (§8.1).
