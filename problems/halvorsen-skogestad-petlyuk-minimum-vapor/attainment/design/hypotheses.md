# Interface file: the specification F(δ, η, π), the limit configuration, and hypotheses (H1)–(H8)

Agent D (design), Claude Opus 5.5 (effort xhigh). **Version 2.1, 2026-10-20 18:47.** Proofs: `design.md` (
v1, complete; each item below says where it is proved). Statements here are kept in sync with `design.md`; a change is
logged at the bottom with the time.

## 0. Notation

Components 1..N, α_1 > … > α_N > 0, feed z > 0 (Σz = 1), quality q ∈ ℝ, global feed roots θ_1 > … > θ_{N−1}
(θ_j ∈ (α_{j+1}, α_j)), peaks H_j, H* := max_j H_j, e_j := H* − H_j.
Columns c = [i,m], 1 ≤ i < m ≤ N; *non-binary* if m − i ≥ 2, *binary* [p, p+1] otherwise. Sections T_c (top),
B_c (bottom), junction stage J_c; product stages P_p = J[p,p] (P_1 total condenser, P_N partial reboiler,
P_p side-draw stage for 2 ≤ p ≤ N−1). T_c's top end is J[i,m−1] (P_i if m = i+1), B_c's bottom end is J[i+1,m]
(P_m if m = i+1).
Section σ: vapor V_σ, liquid L_σ, net upward component flow w_σ (w_{T_c} = d_c, w_{B_c} = −b_c),
Ψ_σ(φ) = V_σ − Σ_i α_i w_{σ,i}/(α_i − φ); Ψ^{(k)}_σ = Ψ_σ without the k-term; G_x(φ) = Σ_i α_i x_i/(α_i − φ),
G^{(k)} the same without the k-term. Eigenvalues of M_σ = V diag(α) − w αᵀ are V_σ·(mode) where a *mode* is a root
of Ψ_σ ("major mode") or α_k for a component with w_{σ,k} = 0 ("trace mode").
Feed of c: f_{[1,N]} = z, U_{[1,N]} = 1 − q; otherwise f_{[i,m]} = d_{[i,m+1]} + b_{[i−1,m]},
U_{[i,m]} = V_{T[i,m+1]} − V_{B[i−1,m]} (absent parents contribute 0). Φ_c(θ) := U_c − Σ α_i f_{c,i}/(α_i − θ)
= Ψ_{T[i,m+1]}(θ) − Ψ_{B[i−1,m]}(θ) (exact flux bookkeeping).

Parameters: δ > 0; η_c > 0 for every non-binary column c; η_p > 0 for every binary column [p,p+1];
a **forcing pattern** π ∈ {B, T}^{N−1} (π_p = B: "B-pinned", component p is sent down B[p,p+1] into product p+1;
π_p = T: "T-pinned", component p+1 is sent up T[p,p+1] into product p). **Default pattern: all B**
(π_p = B for every p). Every statement below holds for every pattern. (N = 3 with π = (B, T) is network.md §11.2.)

## 1. The specification F(δ, η, π) (total flows only; this is what a finite-stage run prescribes)

(F1) Non-binary columns, by decreasing length (top-down). Given the actual feed (f_c > 0 on [i,m], 0 outside;
U_c ∈ ℝ): roots θ^c_i > … > θ^c_{m−1} of Φ_c, one in each (α_{p+1}, α_p); preferred split (V⁰_{T_c}, d⁰_c)
of paper.tex Theorem `achieve` relative to (f_c, U_c); Astra's vector h^c (attainment.md (3)–(4), local roots
θ^c; h^c = 0 for 3 components). Then
  d_c = d⁰_c + η_c h^c, b_c = f_c − d_c, V_{T_c} = V⁰_{T_c} + η_c, V_{B_c} = V_{T_c} − U_c,
  L_{T_c} = V_{T_c} − Σ d_c, L_{B_c} = V_{B_c} + Σ b_c.
(F2) Binary columns [p,p+1]: feed f (supported on {p,p+1}) and U from the parents.
  Vapors: V_{B[N−1,N]} = R := H* + δ; V_{T[p,p+1]} = V_{B[p,p+1]} + U_{[p,p+1]}; V_{B[p−1,p]} = V_{T[p,p+1]}
  (side draws carry no duty). Distillate D_{[p,p+1]} = f_p − η_p (π_p = B) or f_p + η_p (π_p = T);
  L_T = V_T − D, L_B = V_B + (f_p + f_{p+1} − D).
(F3) Products (liquid): P_1 = D_{[1,2]}, P_p = B_{[p−1,p]} + D_{[p,p+1]} (2 ≤ p ≤ N−1), P_N = B_{[N−1,N]}
(B_c := L_{B_c} − V_{B_c}). Feed 1 (quality q) to J[1,N]. The only heat input is the reboiler: R = H* + δ.

## 2. The limit configuration 𝓛(δ, η, π) (infinite-stage sharp state at these flows)

- Net flows: non-binary c as in (F1) (d_{c,m} = 0, d_{c,i} = f_{c,i}; b_{c,i} = 0, b_{c,m} = f_{c,m});
  binary B-pinned: d = (f_p − η_p) e_p, b = η_p e_p + f_{p+1} e_{p+1};
  binary T-pinned: d = f_p e_p + η_p e_{p+1}, b = (f_{p+1} − η_p) e_{p+1};
  every component outside a column's range has zero net flow in its sections ("trace flows").
- Junctions: non-binary X_c = Astra's Lemma 1 vector for (r_c, θ_{i+1}, …, θ_{m−2}, s_c) (attainment.md (9)),
  r_c = root of Ψ_{T_c} in (α_{i+1}, α_i), s_c = root of Ψ_{B_c} in (α_m, α_{m−1});
  binary B-pinned: X_c = normalized (b_k/(φ'_c − α_k))_k, φ'_c = root of Ψ_{B_c} in (α_p, α_{p−1}) (in (α_1, ∞) if p = 1);
  binary T-pinned: X_c = normalized (d_k/(α_k − φ^low_c))_k, φ^low_c = root of Ψ_{T_c} in (α_{p+2}, α_{p+1}) (in (0, α_N) if p+1 = N).
- Products: x_{P_p} = (b_{[p−1,p]} + d_{[p,p+1]})/P_p. Majors at P_p: {p} ∪ {p−1 if π_{p−1} = B} ∪ {p+1 if π_p = T}.

**Mode tables (limit).** "above"/"below" = above/below the pinch mode r*; modes above r* are closed at the
bottom end, modes below r* closed at the top end.

| section | major modes | trace modes | pinch r* | closed at bottom | closed at top |
|---|---|---|---|---|---|
| T_c, c=[i,m] non-binary | r_c, θ_{i+1..m−2}, u_c ∈ (α_m, α_{m−1}) | α_k, k ∉ [i, m−1] | u_c (smallest root) | r_c, θ's; α_k (k<i) | α_k (k ≥ m) |
| B_c, c=[i,m] non-binary | v_c ∈ (α_{i+1}, α_i), θ_{i+1..m−2}, s_c | α_k, k ∉ [i+1, m] | v_c (largest root) | α_k (k ≤ i) | θ's, s_c; α_k (k>m) |
| T_c binary, π_p = B | φ^T_c = α_p L_T/V_T ∈ (α_{p+1}, α_p) | α_k, k ≠ p | φ^T_c | α_k (k<p) | α_k (k>p) |
| B_c binary, π_p = B | φ'_c ∈ (α_p, α_{p−1}), φ''_c ∈ (α_{p+1}, α_p) | α_k, k ∉ {p,p+1} | φ'_c (largest) | α_k (k<p) | φ''_c; α_k (k>p+1) |
| T_c binary, π_p = T | φ^mid_c ∈ (α_{p+1}, α_p), φ^low_c ∈ (α_{p+2}, α_{p+1}) | α_k, k ∉ {p,p+1} | φ^low_c (smallest) | φ^mid_c; α_k (k<p) | α_k (k>p+1) |
| B_c binary, π_p = T | φ^B_c = α_{p+1} L_B/V_B ∈ (α_{p+1}, α_p) | α_k, k ≠ p+1 | φ^B_c | α_k (k ≤ p) | α_k (k>p+1) |

General principle (proved): **every T section pinches at its smallest major root and every B section at its
largest**; every trace mode lies outside the interval spanned by the major modes, so it is closed at exactly one end.

## 3. Hypotheses (statements proved in design.md for the parameter region of §4)

**(H1) Flows.** V_σ > 0, L_σ > 0 for every section; D_c > 0, B_c > 0 for every column; P_p > 0; f_c > 0 on the
range of c. (design.md Lemma 5.2, Props. 3.3, 4.2, 4.3.)

**(H2) Net flows of the limit.** Every T section has w ≥ 0 and every B section w ≤ 0 componentwise, with the
supports listed in §2 (majors strictly nonzero). Non-binary: d_{c,i} = f_{c,i}, d_{c,m} = 0, b_{c,i} = 0,
b_{c,m} = f_{c,m}. Binary: forced as in §2 (b_{[p,p+1],p} = η_p > 0 if π_p = B; d_{[p,p+1],p+1} = η_p > 0 if
π_p = T); never balanced. All trace flows are 0. The middle roots θ^c_{i+1..m−2} of every non-binary column's
feed equal the global θ_{i+1..m−2} exactly (design.md Lemma 5.3), its extreme feed roots are within O(η) of θ_i, θ_{m−1}.

**(H3) Non-binary columns** (Astra's conditions, in an explicit form). For c = [i,m]:
α_i > r_c > θ^c_i > v_c > α_{i+1}, α_{m−1} > u_c > θ^c_{m−1} > s_c > α_m; the middle roots are common simple roots
of Ψ_{T_c} and Ψ_{B_c}; Ψ_{T_c}(α_m) > 0, Ψ_{B_c}(α_i) > 0; the pinch vectors e_c (mode α_m, full support on
[i,m]), e_u (face [i,m−1]), e_v (face [i+1,m]), e_a (mode α_i, full support) are positive on their supports, and
X_c = A_c e_c + B_c e_u = C_c e_v + D_c e_a with A_c, B_c, C_c, D_c > 0 (A+B = C+D = 1).
**Validity condition (explicit):** d⁰_c + η_c h^c > 0 on {i..m−1} and b⁰_c − η_c h^c > 0 on {i+1..m}; nothing
else is needed (design.md Prop. 3.3). For 3-component columns every η_c > 0 is valid.

**(H4) Binary columns.** Margin κ_p := Ψ_{T[p,p+1]}(θ^{[p,p+1]}_p) > 0 at η_p = 0 (θ^{[p,p+1]}_p the binary's
feed root; sharp flows). Then for every admissible η_p (see §4):
π_p = B: Ψ_{B_c}(φ^T_c) > 0, equivalently X_c = A' e^T_{α_{p+1}} + B' e_p with A', B' > 0, where
e^T_{α_{p+1}} = T_c's pinch at the trace mode α_{p+1} (positive on {p,p+1}); X_c > 0 on {p,p+1} with
X_{c,p} → x̂_B ∈ (0,1) as η_p → 0.
π_p = T: Ψ_{T_c}(φ^B_c) > 0, equivalently X_c = C' e_{p+1} + D' e^B_{α_p} with C', D' > 0.
(The forcing only *increases* this margin: design.md Props. 4.2(d), 4.3(d).)

**(H5) Spectra.** For every section all modes are real, positive and pairwise distinct (major roots lie in open
intervals between consecutive α's, trace modes are α's of zero-flow components), the pinch mode is a major mode
(table §2), ρ_σ := max ratio of consecutive modes < 1.

**(H6) Half-profiles.** For every section, the normalized orbits x_k^bot = M^k x_bot/‖·‖₁ (k ≥ 0, from the bottom
end, x_bot without components on modes above r*) and x_k^top = M^{−k} x_top/‖·‖₁ (from the top end) lie in the
closed simplex and converge to e_{r*}; the unnormalized chart vectors e_{r*} + Σ_{above/below} … have
1·(·) bounded below by a positive constant. Proof tools: M_B ≥ 0 (w ≤ 0); M_T^{−1} ≥ 0 (w ≥ 0, L > 0); positive
cone decompositions (H3), (H4).

**(H7) Trace signs.** For every section σ and every trace mode α_k of σ closed at an end J:
(i) Ψ^{(k)}_σ(α_k) > 0; (ii) k is lighter than every major component at J if J is the bottom end (so
G^{(k)}_{x_J}(α_k) < 0) and heavier than every major at J if J is the top end (G^{(k)}_{x_J}(α_k) > 0);
(iii) k is a trace at J (x_{J,k} = 0 in the limit). Hence the linearized closed-end condition
ℓ̃_k·x_J = x_{J,k} + (w_{σ,k}/Ψ^{(k)}_σ(α_k)) G^{(k)}_{x_J}(α_k) = 0 reads x_{J,k} = κ_σ q_σ with
κ_σ = |G^{(k)}_{x_J}(α_k)|/Ψ^{(k)}_σ(α_k) > 0, q_σ = net flow of k from J into σ. Moreover: (iv) every node J
at which k is a trace is the closed end of at least one section for k, or a product node; (v) every component with
nonzero limit net flow in a section is a major at both ends of that section.

**(H8) Duty, purities, pivots.** R = H* + δ exactly (the only duty). Limit product compositions: P_p's impurity
is (η_{p−1}[π_{p−1}=B] + η_p[π_p=T])/P_p = O(η); in the default pattern product 1 is pure and product p ≥ 2 has
impurity η_{p−1}/P_p of component p−1. Nonzero pivots ℓ_{r*}·x at both ends of every section: at a product end
G_{x_P}(r*) ≠ 0 (sign: + at a T top end, − at a B bottom end); at a child's junction, G_{X_child}(u_parent) > 0 and
G_{X_child}(v_parent) < 0; at a non-binary junction, ℓ_u·X = B_c > 0, ℓ_v·X = C_c > 0; at a binary junction the
pinning section has pivot 1 (X is its pinch) and the other section B' > 0 resp. C' > 0. All these quantities are
continuous in (δ, η) on the admissible region, hence bounded away from 0 on compact subsets of it.

## 4. Parameter region (proved: design.md Proposition 5.4, Theorem 5.5, Theorem 6.2 = Theorem D)

There are η_1 > 0, C_1 ≥ 0 depending only on (α, z, q) such that, for every pattern π and every δ > 0, with
η̄(δ) := min(η_1, δ/(2C_1)), (H1)–(H8) hold whenever
- 0 < η_c ≤ η̄(δ) for every non-binary column c (this gives the explicit validity condition of (H3) along the
  whole top-down cascade and κ_p ≥ e_p + δ/2 for every binary column; exactly: κ_p(η_nb; δ) = κ_p(η_nb; 0) + δ and
  κ_p(0; 0) = e_p), and
- for every binary column 0 < η_p < η̂_p (design.md (4.3), (4.6)):
  π_p = B: η̂_p = f_p (p = 1), min(f_p, (α_{p−1} − α_p) Ψ⁰_B(α_{p−1})/α_p) (p ≥ 2) (keeps φ' < α_{p−1});
  π_p = T: η̂_p = min(f_{p+1}, V_T − f_p) (p+1 = N), min(f_{p+1}, (α_{p+1} − α_{p+2}) Ψ⁰_T(α_{p+2})/α_{p+1})
  (p+1 ≤ N−1) (keeps φ^low > α_{p+2} and L_T > 0);
  here Ψ⁰_T(t) = V_T − α_p f_p/(α_p − t), Ψ⁰_B(t) = V_B + α_{p+1} f_{p+1}/(α_{p+1} − t) are the unforced
  functions of the binary column.
No hierarchy between parents and children is needed; the binary η_p need not be small relative to κ_p (the
forcing increases the margin: Ψ_B(φ^T) > κ_p resp. Ψ_T(φ^B) > κ_p).
Numbering in design.md: (H1)–(H8) = Theorem 6.2; Astra's conditions = Proposition 3.3; binary = Propositions
4.1–4.3; flows/balances = Lemma 5.2; middle roots global = Lemma 5.3.

## 5. Export for Agents C, L, A

`design/scripts/export_config.py CASE [out.json] [dps]` builds F(δ, η, π) and 𝓛 (mpmath) and writes JSON
(format "hs-petlyuk-limit-config v1", documented in the script header; 1-based indices). Examples:
`design/runs/config_*.json`. Independent recomputation by C: `computation/runs/compare_D.log` (agreement ≤ 1.2e-54).

## Change log
- 18:07 v1.
- 18:36 v2: references to design.md numbering; precise parameter region (η̄(δ) linear in δ; κ identity);
  T-pinning bound for the last binary column now includes L_T > 0 (η_p < V_T − f_p); export script name.
  No change of substance in (H1)–(H8).
- 18:47 v2.1: wording only (binary section ends); design.md v1.1 defines "admissible" as the explicit, checkable
  conditions (i) validity of every non-binary column, (ii) κ_p > 0, (iii) 0 < η_p < η̂_p — (H1)–(H8) hold for every
  admissible parameter; §4 above is one admissible choice for every δ. Self-contained proof of κ_p(0; δ) = e_p + δ
  (design.md Prop. 5.4(c)).
