# Independent verification V2: the analysis and the integrated proof (result 12, version 3)

Verifier: Claude Opus 5.5 (effort xhigh), fresh context. Started 2026-10-20 19:17, verdicts written 19:44, final 20:17
(clock times in `PROGRESS.md`).
Scope (BRIEF-V2): `attainment/analysis/analysis.md` (v1.2) and the integrated proof `attainment/theorem.md` (v2.2) end
to end, with priority on §6 (perturbation, reconstruction, positivity) and §7 (proof of Theorem 1), Corollary 2, the N = 2
case (`total/sections/sections.md` §8.1, Proposition 10), the interface with the lower bound (`total/foundations/class.md`
§§4–5, Theorem 5.3, Corollary 5.5; `total/verification-early/verification-lower-bound.md`) and with
`total/network/network.md` Theorem 10.1 (i). Paths are relative to the problem folder.

**Read:** PROTOCOL.md; BRIEF-verification-common.md; BRIEF-V2.md; theorem.md; analysis.md; class.md; network.md §10;
sections.md §§1, 8.1–8.4; paper.tex Theorem `achieve`; attainment-astra/attainment.md (Lemmas 1–2, located);
limit.md §7 (Remark 7.3, cited by theorem.md §7); the early verification of the lower bound; HS 2003 Part 3
(`certificates/sources/HS2003_part3.txt` and p. 616 of the PDF, Figure 1; §6.3); the scripts and logs that I re-ran.
**Not opened:** any PROGRESS.md other than mine, any FINAL-REPORT.md, any BRIEF file other than the two above (the early
verifier's BRIEF included), V1's folder before my verdicts were written (§8 below says when I opened it).

---

## 1. Verdicts

| # | Claim | Verdict |
|---|---|---|
| 1 | Exact reformulation of the finite-stage steady states (theorem.md §2.4–2.5, Proposition 2.1; analysis.md §1, Lemma A.1), including the positivity of every mole fraction by the M-matrix argument (network.md Theorem 10.1 (i)) at every kind of unit | **CONFIRMED** |
| 2 | Section algebra used by the analysis (Lemmas 3.1, 3.3, 3.4, 3.5 = analysis.md Lemmas A.2, A.3): spectrum, Underwood and trace left eigenvectors, product form, C¹ estimate | **CONFIRMED** |
| 3 | Perturbation step (Lemma 5.2, Lemma 6.1, Theorem 6.2 = analysis.md Lemma A.4, Theorem A.5): unique zero v_n near 𝓛 for all section lengths ≥ n_0, \|v_n − v_∞\| ≤ C ρ_∞^{n+1} | **CONFIRMED** |
| 4 | From v_n to a physical steady state (Lemmas 6.3, 6.4, Theorem 6.5 = analysis.md A.6–A.8): exact profiles, S > 0 at every stage uniformly, positivity, purity, R | **CONFIRMED WITH CORRECTIONS** (one wording fix in the proof of Theorem 6.5; no mathematical change) |
| 5 | Proof of Theorem 1 (§7): choice of parameters, purity, R = H* + δ, order of limits n → ∞, η → 0, δ → 0, dependence of the constants | **CONFIRMED** (given §§4–5, which are V1's scope; see §8) |
| 6 | The arrangement of §2.2 is the generalized Petlyuk arrangement of HS 2003 Part 3, Fig. 1 | **CONFIRMED** |
| 7 | N = 2 (sections.md Proposition 10, checked line by line) | **CONFIRMED WITH CORRECTIONS** (Proposition 10 does not state the positivity claimed by Theorem 1; it holds, one sentence added to theorem.md §1) |
| 8 | Corollary 2: definition of R* (class.md §4), use of the lower bound (class.md Theorem 5.3, Corollary 5.5), effect of the early verifier's corrections | **CONFIRMED WITH CORRECTIONS** (editorial: theorem.md §10 called the early check "provisional"; it is final, and its corrections do not touch the bound) |
| 9 | analysis.md as a whole (including Corollary A.9 and the §8 interface with D and L) | **CONFIRMED** (Remark 7.6 is an unproved side remark, not used) |
| 10 | Consistency of §§4–5 with what §§6–7 use ((H1)–(H8), Lemma 5.1, Theorem 5.9) | consistent; spot re-derivations of Proposition 4.10 (c), Proposition 5.7, Proposition 5.8 and Lemma 5.3 found no gap; the full check is V1's |
| **Overall** | **Theorem 1 and Corollary 2, i.e. the Halvorsen–Skogestad conjecture in the class 𝒜 in the finite-stage-limit sense of class.md §4, for every N ≥ 2** | **CONFIRMED WITH CORRECTIONS** (minor and editorial only), jointly with V1's verdict on §§2–5 (§8) |

No gap of substance was found in my scope. Numerical checks (§5) are evidence only; every verdict rests on the
re-derivations of §2.

---

## 2. Re-derivations (by hand; then checked numerically where useful)

### 2.1 The unknowns and the exact reformulation (Proposition 2.1, Lemma A.1)

*Total balances.* For v ∈ 𝒱, 1ᵀd_c = Δ_{T_c} by definition. By induction over the columns (parents first),
1ᵀf_c = Δ_{T[i,m+1]} − Δ_{B[i−1,m]} (+1 for the feed column), and the total balance of F at J_c is
Δ_{T_c} − Δ_{B_c} = 1ᵀf_c, so 1ᵀw_{B_c} = −1ᵀb_c = Δ_{T_c} − 1ᵀf_c = Δ_{B_c}; at J[p, p] the total balance gives
1ᵀ(b_{[p−1,p]} + d_{[p,p+1]}) = P_p, so 1ᵀx_p = 1. Then 1ᵀM_σ(w_σ) = (V_σ − Δ_σ)αᵀ = L_σαᵀ.

*(a)* At a steady state every composition is in Δ, S ≥ α_N > 0; the net upward flow V y(x_k) − L x_{k+1} is the same at
every cross-section of a section (difference of two cross-sections = the stage balance), equal to d_c in T_c and, by the
junction balance, to −b_c in B_c; L S(x_k) x_{k+1} = M x_k, hence x_top = M^{n+1}x_bot/Π_k(L S(x_k)), c_σ > 0.

*(b)* Induction gives 1ᵀM^k x_bot = L^k Π_{j<k} S(x_j) (using 1ᵀx_bot = 1), so the hypothesis "iterates defined" is
automatic once S(x_j) > 0 for j < k; x_{k+1} = M x_k/(L S(x_k)) is exactly the cross-section balance; 1ᵀx_top = 1 forces
c = 1/(L^{n+1}Π S(x_k)) > 0 and x_{n+1} = x_top. Stage balances are differences of consecutive cross-section balances.
End units: the vapor entering stage 1 comes from the bottom-end unit, always an (E) stage (junction, side draw or the
reboiler; the condenser is never a bottom end), so its composition is y(x_bot); the liquid entering stage n is the top
unit's liquid x_top (for T[1, 2] the condenser's liquid). Node balances (out − in), with the end cross-section
balances: Σ_{bot(σ)=J} w_σ − Σ_{top(σ)=J} w_σ − feed + product. I re-derived the attachment list of every node
(J[i, m], i < m: bottom end of T[i, m] and B[i−1, m], top end of B[i, m] and T[i, m+1]; J[p, p]: bottom end of
B[p−1, p], top end of T[p, p+1]; nothing else) and obtained d_c − b_{[i−1,m]} + b_c − d_{[i,m+1]} − z[c = [1,N]] =
d_c + b_c − f_c = 0 at junctions and −b_{[p−1,p]} − d_{[p,p+1]} + P_p x_p = 0 at product nodes. Condenser: the top
cross-section of T[1, 2] gives V y(x_n) = L x_1 + d_{[1,2]} = (L + P_1)x_1 = V x_1, which is the composition rule of a
total condenser. Total and CMO balances are those of F (Lemma 4.8, re-checked numerically in §5), so the duties are
those of F: 0 at junctions, side draws and stages, −V_{T[1,2]} at the condenser, V_{B[N−1,N]} = R at the reboiler.

*Positivity.* With K_{u,i} = α_i/S(x_u) > 0 at every unit with a vapor outflow (all units but the condenser), the
component-i balances of all units (stages included) form 𝔄_i X_{·,i} = e_{J[1,N]} z_i with diagonal L_u + V_u K_{u,i} > 0
(total outflows, products included) and off-diagonal −(ℓ_{u'u} + v_{u'u}K_{u',i}) ≤ 0. The excess of column u' is the
flow that u' sends to product outlets (liquid products here). 𝔄_iᵀ is weakly chained diagonally dominant: every unit
reaches a product unit along positive flows (liquid runs down every section into a junction, side draw or the
reboiler, and down the grid of columns to some J[m, m]); hence 𝔄_i is a nonsingular M-matrix, 𝔄_i^{−1} = Σ_j
(D^{−1}O)^j D^{−1} ≥ 0 with (u, J[1,N]) entry > 0 iff u is reached from J[1, N] along positive flows, which holds for
every unit (vapor up and liquid down every section, both couplings at every junction). So x_{u,i} > 0 for every unit
and component, and y = Kx > 0. This uses only S > 0, positive flows and the normalization Σx = 1 (which holds by
construction: x_k = u_k/1ᵀu_k at stages, 1ᵀX_c = 1 and 1ᵀx_p = 1 at nodes). It does not need sign information on trace
flows. Note: network.md Theorem 10.1 is stated for solutions with positive x, but its part (i) only uses K > 0; theorem.md
Proposition 2.1 re-proves exactly what is needed. Every unit type is covered: stages, junctions (two thermal couplings
each, with both liquid and vapor in both directions), side-draw stages, the total condenser (an (E)-type row with no
vapor outflow), the partial reboiler.

### 2.2 The section algebra (Lemmas 3.1, 3.3, 3.4, 3.5; analysis.md A.2, A.3)

- det(Vt − M) = Π(Vt − Vα_i)(1 − V^{−1}Σα_iw_i/(α_i − t)) = (−1)^N V^{N−1}P_w(t) (rank-one determinant). Re-derived.
- 1ᵀe_φ = L/φ also for a trace mode: Ψ^{(k)}(α_k) = L − α_k Σ_{i≠k} w_i/(α_i − α_k) (from α_i/(α_i − α_k) = 1 +
  α_k/(α_i − α_k) and w_k = 0), so 1ᵀe_{α_k} = Σ_{i≠k} w_i/(α_i − α_k) + Ψ^{(k)}(α_k)/α_k = L/α_k. This is what allows the
  normalization 1ᵀê_r = 1 of **every** mode near v_∞ in Lemma 6.3 (needed there, not only for the pinch).
- Lemma 3.3 (b): Sherman–Morrison with 1 − αᵀ(V diag α)^{−1}w = L/V. Re-derived.
- Lemma 3.4: P_w(t) = (α_k − t)Π_{j≠k}(α_j − t)Ψ^{(k)}(t) − α_kw_kΠ_{j≠k}(α_j − t), so at a root t_r near α_k,
  (α_k − t_r)Ψ^{(k)}(t_r) = α_kw_k; if w_k = 0 the root near α_k is exactly α_k (uniqueness of the simple root), so
  ℓ̃ = e_k there; if w_k ≠ 0, ℓ̃ = ((α_k − t_r)/α_k)ℓ_{t_r}. ℓ̃ is analytic (t_r analytic, away from α_j, j ≠ k, and
  Ψ^{(k)}(t_r) ≠ 0 near w⁰). Its derivative at w_k = 0 = x_k is δx_k + δw_k G^{(k)}_x(α_k)/Ψ^{(k)}(α_k). Re-derived, and
  checked numerically: my network map uses exactly these formulas (not generic eigenvectors) and reproduces G_∞(v_∞) = 0
  and the zeros found with generic eigen-data by A's script (§5).
- Lemma 3.5 (a): with the dual basis, x_top = cM^{n+1}x_bot ⟺ ℓ_r·x_top = cλ_r^{n+1}ℓ_r·x_bot for all r; the r*-relation
  fixes c = π_top/(λ_{r*}^{n+1}π_bot) (needs π_bot ≠ 0), and the others become E_r = 0 (r ∈ U, multiply by
  π_botλ_{r*}^{n+1}/λ_r^{n+1}) and E_s = 0 (s ∈ D, multiply by π_bot). Conversely E = 0 with π_botπ_top ≠ 0 gives all
  relations. The zero set does not depend on the analytic normalization of the ℓ's (factors h_rh_{r*}). (c): the
  difference E^{(n)} − E^{(∞)} is a sum of μ^{n+1}× (analytic), μ ≤ ρ̄, ∇μ^{n+1} = (n+1)μ^{n+1}∇log μ: |·| ≤ Cρ̄^{n+1},
  ‖D·‖ ≤ C(n+2)ρ̄^{n+1}, with C = sup over the compact set of the eigen-data and their first derivatives, **independent of
  n**. At v_∞ the only μ's are the ratios of λ_{r*} to the other modes, the largest being the ratio to an adjacent mode:
  |G_n(v_∞)| ≤ Cρ_∞^{n+1}. Re-derived.

### 2.3 The perturbation step (Lemma 5.2, Lemma 6.1, Theorem 6.2)

*Domain.* At v_∞ every section has N simple positive modes (H5) and nonzero pivots (H8). Simple real roots stay real,
simple and ordered on a neighborhood; the eigen-data (Underwood for major modes, (3.1) for trace modes) are analytic
there; the pivots stay ≥ π_0/2. So there is a compact ball 𝒩 = B̄(v_∞, r_0) on which G_∞ and every G_n are analytic,
ρ̄ := sup max consecutive ratio < 1, and the C¹ bounds of Lemma 3.5 (c) hold with constants independent of the n_σ.

*Lemma 6.1* (Newton–Kantorovich type): Φ(v) = v − A^{−1}F(v), DΦ = A^{−1}(A − DF_∞(v)) + A^{−1}(DF_∞(v) − DF(v)), so
‖DΦ‖ ≤ 1/4 + 1/4 on the convex ball B̄(v_∞, r_1); |Φ(v) − v_∞| ≤ |Φ(v) − Φ(v_∞)| + ‖A^{−1}‖|F(v_∞)| ≤ r_1/2 + r_1/2; Banach
gives exactly one fixed point (= zero of F) in the ball, and |v* − v_∞| ≤ |v* − v_∞|/2 + ‖A^{−1}‖|F(v_∞)|. Correct.

*Are its hypotheses checked in Theorem 6.2?* Yes, all five: F_∞(v_∞) = G_∞(v_∞) = 0 (Lemma 5.1); A = DG_∞(v_∞)
invertible (Theorem 5.9); ‖A^{−1}‖ω(r_1) ≤ 1/4 from ω(r) ≤ C_3r (G_∞ is C² on 𝒩) and r_1 := min(r_0, 1/(4C_3‖A^{−1}‖));
sup_{B̄}‖DG_n − DG_∞‖ ≤ C_1 max_σ(n_σ + 2)ρ̄^{n_σ+1} ≤ C_1(n + 2)ρ̄^{n+1} ≤ 1/(4‖A^{−1}‖) for n ≥ n_0 (the function
x ↦ (x + 2)ρ̄^{x+1} decreases for x > 1/ln(1/ρ̄) − 2, so n ≥ 1/log(1/ρ̄) suffices for the max over unequal lengths);
‖A^{−1}‖|G_n(v_∞)| ≤ ‖A^{−1}‖C_2ρ_∞^{n+1} ≤ r_1/2 for n ≥ n_0. Constants C_1, C_2, C_3, ‖A^{−1}‖, r_0, r_1 depend on
(α, z, q, δ, η, π) only — never on (n_σ) — so the conclusion holds simultaneously for **all** families of section
lengths with min n_σ ≥ n_0, equal or not. Uniqueness: exactly one zero of G_n in B̄(v_∞, r_1); since every steady state
whose (d_c, X_c) lies in that ball is a zero of G_n (Proposition 2.1 (a) + Lemma 3.5 (a), pivots nonzero on 𝒩), the
steady state is locally unique. Displacement: |v_n − v_∞| ≤ 2‖A^{−1}‖C_2ρ_∞^{n+1} — the rate ρ_∞ (pinch-adjacent ratio
at v_∞), not nρ^n, because only the residual at v_∞ enters the displacement. All correct.

*Order of the choices.* r_0 (domain) → r_1 (needs only G_∞) → n_0 (needs r_1 and the C¹ bound). No circularity.

### 2.4 Reconstruction, S > 0, purity (Lemmas 6.3, 6.4, Theorem 6.5)

*Lemma 6.3.* With 1ᵀê_r = 1 for every r (possible by §2.2) and ℓ_r·ê_r = 1, M^kx_bot = Σ_r(ℓ_r·x_bot)λ_r^kê_r; for r ∈ U,
E_r = 0 gives ℓ_r·x_bot/π_bot = μ_r^{n_σ+1}b_r, whence (6.1). Exact for every 0 ≤ k ≤ n_σ + 1. Re-derived.

*Lemma 6.4.* c_0 > 0 because v^bot_k ≥ 0 (H6), v^bot_k ≠ 0 (M invertible, x_bot ≠ 0), and 1ᵀv^bot_k → 1ᵀê_{r*} = 1 (a
positive sequence with a positive limit has a positive infimum). Termwise comparison: the coefficients (ê_r, μ_r, a_s,
b_r) are analytic on the ball and Lipschitz, so they move by ≤ C|v_n − v_∞| ≤ Cρ_∞^{n+1}; |μ^k − μ_∞^k| ≤
(sup_j jρ̄^{j−1})|μ − μ_∞| uniformly in k (both ≤ ρ̄); the U-terms are ≤ Cρ̄^{n_σ+1−k} (b_r bounded because |π_top| ≥
π_0/2), the D-terms ≤ Cρ̄^k. Both comparisons are made with the **same** u_k, whose pinch coefficient is exactly 1 —
this is the point of the product form. For every k, min(ρ̄^k, ρ̄^{n_σ+1−k}) ≤ ρ̄^{(n+1)/2}, so u_k is within
ε_n = C(ρ_∞^{n+1} + ρ̄^{(n+1)/2}) of a vector ξ ≥ 0 with 1ᵀξ ≥ c_0; then 1ᵀu_k ≥ c_0 − Nε_n and αᵀu_k ≥ α_N1ᵀξ −
(Σα)ε_n (ξ ≥ 0 is used here), so S(x_k) > 0 for n ≥ n_1. **Uniform in k and in n_σ** (the constants do not depend on
n_σ, and a long section only makes the far-end terms smaller). Correct. (Observation: what is proved is the sharper
dist(u_k, {v^bot_k, v^top_{n_σ+1−k}}) ≤ C(ρ_∞^{n+1} + ρ̄^{max(k, n_σ+1−k)}); Theorem 6.5 (iii) states the weaker
ρ̄^{min(k, n_σ+1−k)}, which is true. No change made.)

*Theorem 6.5.* S > 0 at the node compositions: they are affine in v and within Cρ_∞^{n+1} of points of Δ, where S ≥ α_N.
The proof says "shrink r_1 so that S ≥ α_N/2 on the ball", which, read literally after n_0 has been fixed, would require
re-choosing n_0; the clean statement is to enlarge n_1 (correction C3 below; no mathematical change). Then Proposition 2.1
(b) applies: a steady state of the arrangement with flows F, all mole fractions > 0, R = R(F) = V_{B[N−1,N]} = H* + δ
(Lemma 4.8). Purity: x_p is affine in v, so x_{p,p} ≥ x^∞_{p,p} − Cρ_∞^{n+1}, and 1 − x^∞_{p,p} = (η_{p−1}[π_{p−1} = B] +
η_p[π_p = T])/P_p — I re-derived this from the pinning formulas (B-pinning of [p−1, p] sends η_{p−1} of component p−1
into product p; T-pinning of [p, p+1] sends η_p of p+1 into product p; the other two pinnings add no impurity).

### 2.5 Proof of Theorem 1 (§7), order of limits and constants

Given δ, ε: η_c := η̄(δ) for every non-binary column (Theorem 4.11: valid columns, κ_p ≥ e_p + δ/2 > 0). The bounds η̂_p
depend on the binary feeds and vapors, which depend on (δ, η_nb) only (Proposition 4.10 (a), Lemma 4.8). P_p = z_p ±
η_{p−1} ± η_p (re-derived: the binary feeds carry all of component p into the two columns adjacent to product p), so
P_p → z_p > 0 as η_bin → 0 and the designed impurity ≤ ε/2 for η_bin small. The parameters are admissible; Theorem 6.5
gives n_1(δ, η, π) and C(δ, η, π); n_0 ≥ n_1 with Cρ_∞^{n_0+1} ≤ ε/2. For all n_σ ≥ n_0 the steady state has purity ≥ 1 −
ε and R = H* + δ, so R*(𝒮_Petlyuk; ε) ≤ H* + δ for all δ, ε > 0, i.e. R*(𝒮_Petlyuk; ε) ≤ H* for every ε and
R*(𝒮_Petlyuk) ≤ H*. The order is: fix (δ, ε) → η → n_0. No uniformity of n_0 in (δ, η) is needed because R*(𝒮; ε) is an
infimum over all members (all stage numbers and all flows) of 𝒮. The specification F(δ, η, π) is independent of the
n_σ, as Theorem 1 states. The non-uniformity (n_0 ≳ log(1/η)/log(1/ρ_∞), pivots O(η)) is stated correctly in §7 and
is visible in my numerics: the neighborhood of real simple modes has radius < 10⁻³ (sup norm) at η_bin = 0.003 for N = 4,
and the contraction of Lemma 6.1 is observed only once 2‖A^{−1}‖|G_n(v_∞)| is below that radius (§5).

### 2.6 The arrangement is HS's generalized Petlyuk arrangement

HS 2003 Part 3, Fig. 1 (p. 616, examined in the PDF) shows, for four products, the columns C1 (ABCD), C21 (ABC), C22
(BCD), C31 (AB), C32 (BC), C33 (CD), every column a two-product column whose top and bottom are directly (two-way
liquid/vapor) coupled to the feed points of its children, C32 fed at one junction from the bottom of C21 and the top of
C22 (text p. 617: "The mixing of flows from C21 and C22 at the feed junction to C32 does not give any problems"), one
total condenser (top of C31, product A), one reboiler (bottom of C33, product D), side products B and C drawn between
C31/C32 and C32/C33, and the caption "Vapor and liquid flow rates can be set individually in each internal two-product
column". theorem.md §2.2 is the same graph for general N (column [i, m] ↔ the column separating components i..m;
J[i, m] its feed junction, receiving T[i, m+1] from above and B[i−1, m] from below; J[p, p] the side-draw stage). In
the figure the side products are drawn at the link between consecutive last-row columns; modelling the draw point as a
duty-free equilibrium stage J[p, p] is the same family 𝒮_Petlyuk (that stage can be counted as the bottom stage of
B[p−1, p]), since all stage numbers are allowed. HS's remark that "some of the heat may be supplied or removed at the
product outlets" is a freedom not used: here all heat enters at the reboiler. **CONFIRMED.**

### 2.7 N = 2: sections.md Proposition 10, line by line

Data α_1 = α > 1 = α_2, light fraction z, θ ∈ (1, α), H_1 = (1 − z)/(θ − 1); flows D = z, B = 1 − z, V_B = s,
V_T = s + 1 − q, L_T = V_T − D, L_B = s + 1 − z (L_B − L_T = q: feed-stage liquid balance ✓; V_T − V_B = 1 − q ✓).
(i) H_1 + 1 − q = αz/(α − θ) > z (feed equation), so L_T > 0 for s ≥ H_1; H_1 > (1 − z)/(α − 1) since θ < α. ✓
(ii) β_0(x) = (s/L_B)e(x) has fixed points 0 and x̂_B = (αs/L_B − 1)/(α − 1) ∈ (0, 1) (s(α − 1) > 1 − z, s < L_B), and
β_0(x) > x on (0, x̂_B) (e(x)/x decreasing). Top map τ_0(x) = e^{−1}((L_Tx + D)/V_T): 1 − e(x) = (1 − x)/(1 + (α − 1)x)
against 1 − line = (L_T/V_T)(1 − x) gives the second fixed point x̂_T = (V_T/L_T − 1)/(α − 1) and τ_0(x) < x on (x̂_T, 1)
(strict concavity: no other fixed point). x̂_T < x̂_B ⟺ αsL_T − V_TL_B = f(s) > 0; f(s) = (α − 1)(s − H_1)(s − s_−)
identically (sympy, my script) and H_1 − s_− = zα(θ − 1)/((α − θ)(α − 1)) + (1 − z)/(θ − 1) − (1 − z)/(α − 1) > 0, so
s > H_1 ⇒ f(s) > 0. ✓
(iii) Shooting in η = x_B ∈ [0, η_1], η_1 max(1, (1 − z)/z) ≤ ε (so x_D(η) = 1 − (1 − z)η/z ∈ [1 − ε, 1]; the proof says
"η ∈ [0, 1]", but only [0, η_1] is used, where x_D ∈ [0, 1] — harmless). Bottom iterates: n_B + 1 applications of
β_η from the reboiler liquid; top: e^{−1}(x_D) on the top stage and n_T applications of the top step down to the feed
stage. Every balance re-derived: reboiler = first bottom cross-section; feed stage = difference of the two adjacent
cross-section balances + overall light balance Dx_D + Bx_B = z; condenser: V_T = L_T + D and the top-stage vapor is x_D.
Comparison: β_{η_1} ≥ β_0, both increasing ⇒ x_f^bot(η_1) ≥ min(x̂_B, β_0^{n_B+1}(η_1)) → x̂_B (both cases η_1 < x̂_B
and η_1 ≥ x̂_B); the top step with x_D ≤ 1 is ≤ τ_0, both increasing, start e^{−1}(x_D(η_1)) < 1 ⇒ x_f^top(η_1) →
≤ x̂_T. IVT on [0, η_1] (mismatch −1 at η = 0). ✓ **Positivity (needed by Theorem 1, not stated in Proposition 10):**
x_B = η* > 0, x_D ∈ (0, 1); β_{η*} maps (0, 1) into (0, 1) (s e(x) + Bη* < s + B = L_B); the top step is a convex
combination of points of (0, 1) followed by e^{−1}, which maps (0, 1) onto (0, 1): every mole fraction is in (0, 1).
Hence Theorem 1 holds for N = 2 as stated (with s = H_1 + δ, R = s, the same flows for all n_T, n_B ≥ n_0).
Correction C1 adds this to theorem.md §1. Numerically confirmed (§5, `v2_binary_prop10.py`, six data sets including
q = −0.3, q = 1.3, q = 1.6, unequal n_T ≠ n_B, every balance ≤ 1e-40, every x in (0, 1)).

### 2.8 Corollary 2 and the lower bound

class.md §4: R*(𝒮; ε) = inf{R(s) : s a steady state of a member of 𝒮, every product stream of purity ≥ 1 − ε},
nondecreasing as ε ↓ 0, R*(𝒮) its limit; inf ∅ = +∞. Theorem 1 gives R*(𝒮_Petlyuk) ≤ H*.
I re-derived class.md (v3) §5 completely: Lemma 5.1 (φ linear in (v_e, n_e); subtract Σc_i × component balances from the
CMO balance), Lemma 5.2 (liquid φ = −lG_x; vapor φ = −θvΣy_i/(α_i − θ) = −(θv/S)G_x; θ_j > 0, S > 0 ⇒ every outflow of a
unit has the sign of −G, or one sign for one-phase units), Theorem 5.3 (sum over 𝒰₊; internal streams cancel; feed
terms vanish since φ_F(θ_j) = 0; crossing terms ≥ 0; R = R(𝒰₊) + R(𝒰₋)), Corollary 5.4, and Corollary 5.5 (vapor flux
= liquid flux + total flow ≥ liquid flux; overall component balance n^h_i = z_i − m_i; the m_i (i > j) and n^h_i (i ≤ j)
are disjoint parts of the impurity, ≤ ε; coefficients ≤ c_j by monotonicity). So R ≥ H_j − c_jε at every steady state of
every arrangement of 𝒜 with purity ≥ 1 − ε, R*(𝒮) ≥ H* for every 𝒮 ⊆ 𝒜. The Petlyuk steady states of Theorem 1 are
steady states of 𝒜 in the sense of class.md §3 (nonnegative, indeed positive, stream flows; (E) stages, one (C) unit,
liquid products; R = Σ max(Q_u, 0) = V_{B[N−1,N]}). Corollary 2 follows.

*The early verification* (`total/verification-early/verification-lower-bound.md`, read in full) confirmed the bound
(claims 1–2) and confirmed claim 3 with corrections. Its corrections are: renaming of 𝒰± in network.md Prop. 6.3;
"constant" in class.md §3 only for fixed product phases; same-composition two-phase units, cross-flow trays with
E_MV > 1 and scalar-coefficient rate models are covered (the team files were too pessimistic); add non-steady
(batch/cyclic) operation and component-specific mass-transfer rates to the exclusions; "lowest" means "lower or equal";
"conjecture ⇔ attainment" holds for the conjecture including HS's formula, attainment being sufficient for part (A).
**None affects the bound or Corollary 2**: they widen or clarify the scope, and the two new exclusions are outside HS's
setting (continuous steady state, equilibrium stages). The early verifier read class.md v2 and sections.md/network.md v1;
class.md v3 §5 restates the same proof, which I re-derived above. theorem.md §10 calls the early check "provisional"; the
file says "final for this early round" (correction C2).

### 2.9 Spot checks of §§4–5 used by the chain (V1's scope; done only to secure the interface)

- Proposition 4.10 (c) (κ_p(0; δ) = e_p + δ by the light-side cut): re-derived the classification of sections with
  exactly one light end (only T[i, p+1], 1 ≤ i ≤ p), the identity Σ_rφ_{P_r}(t) = R − V_{T[1,2]} + φ_F(t) and its
  t → ∞ limit V_{T[1,2]} = R + 1 − q, and the final evaluation with the feed equation. ✓ (Numerically κ_p − e_p − δ =
  O(η_nb), as Proposition 4.10 (d) says: −0.017, −0.022 at η_nb = 0.025, N = 3; −0.008..−0.010 at η_nb = 0.005, N = 4.)
- Lemma 5.3 (counts), Lemma 5.6 (node identities, identically on 𝒱), Proposition 5.7 (the linearized trace rows ξ_J =
  κδq at the closed end and the node identity give A_kξ = 0; the signs of (H7) give −κδq at both kinds of ends),
  Proposition 5.8 (two rows of a common middle root give R(θ_j)(1/Ψ'_T − 1/Ψ'_B) = 0, Ψ'_T < 0 < Ψ'_B; with R(0) = Σδd = 0,
  a polynomial of degree ≤ m − i − 2 with m − i − 1 zeros; then the interpolation matrix): ✓.
- My own implementation of Definition 4.7 / Definition 4.12 for general N reproduces D's exported configurations to
  ≤ 5e-55 (N = 3 two patterns, N = 4, N = 5), and Lemma 4.9 (§5).

---

## 3. Gaps found

| # | Where | Severity | Description | Fix |
|---|---|---|---|---|
| G1 | theorem.md §1 (N = 2) | minor | Theorem 1 claims all mole fractions positive; Proposition 10 (cited for N = 2) proves purity and R but does not state positivity | holds (§2.7); sentence added (C1) |
| G2 | theorem.md §10 | editorial | "early independent check CONFIRMED, provisional": the early verification is final; its corrections are editorial/scope only | reworded (C2) |
| G3 | theorem.md Theorem 6.5, proof | minor (wording) | "shrink r_1 so that S ≥ α_N/2 on the ball" after n_0 is fixed would need a new n_0 | replaced by "enlarge n_1" (C3); analysis.md Theorem A.8's proof has the same phrase (C4) |
| G4 | theorem.md Theorem 6.5 (iii), analysis.md A.8 (iii) | observation | the stated rate ρ̄^{min(k, n_σ+1−k)} is weaker than what Lemma 6.4 proves (ρ̄^{max(…)} for the nearer half-profile); the statement is true | none |
| G5 | analysis.md Remark 7.6 | observation | "a partial condenser or a vapor side product changes only the affine maps …" is not proved (the design and Theorem L would also have to be redone) and is not used | none (remark only) |
| G6 | sections.md Prop. 10 (iii) | observation | "for η ∈ [0, 1]" while x_D(η) < 0 for η > z/(1 − z); only η ∈ [0, η_1] is used | none (not in my files) |
| G7 | network.md Thm 10.1 | observation | stated for solutions with positive x; Proposition 2.1 uses only part (i) with K > 0 and re-proves it | none |
| G8 | Theorem 6.2 / 6.5, regime n ≥ n_0 | observation (numerical) | below n_0 the zero of G_n near 𝓛 need not be physical: in case n3TT (N = 3, q = −0.5, pattern TT) at n = 16, where 2‖A^{−1}‖\|G_n(v_∞)\| = 519 (Lemma 6.1's hypotheses fail by far), damped Newton found a zero of G_n at distance 0.066 from v_∞ whose re-shot profile has negative mole fractions (min x = −5.3); at n = 32, 48, 64 and mixed lengths the zero is physical. Likewise n3BT at n = 8 (2‖A^{−1}‖\|G_n(v_∞)\| = 289; zero at distance 0.53, min x = −2.6). This is what theorem.md §8.5 says (S > 0 is a genuine hypothesis of Proposition 2.1 (b), delivered only for n ≥ n_1); no contradiction | none |

No gap of the kinds listed in the common brief (missing case, non-uniform constant used as uniform, wrong order of
limits, purity/R* definition, positivity, mismatch with the physical equations) was found in my scope.

---

## 4. Corrections applied

### 4.1 V2's corrections

**C1 — theorem.md §1, paragraph after the proof of Corollary 2.**
Old: "For N = 2, Theorem 1 is Agent S's Proposition 10 (sections.md §8.1). The proof below is written for N ≥ 3; it
applies verbatim to N = 2 (one binary column, no non-binary column), but we do not need this."
New: "For N = 2, Theorem 1 is Agent S's Proposition 10 (sections.md §8.1), applied with V_B = s = H_1 + δ. Proposition 10
does not state the positivity claimed in Theorem 1, but its steady state has every mole fraction in (0, 1): x_B = η* > 0
and x_D = 1 − (1 − z)η*/z ∈ (0, 1), the bottom step x ↦ (s e(x) + Bη*)/L_B maps (0, 1) into (0, 1), and the top step is a
convex combination of points of (0, 1) followed by e^{−1} (V2, `verification-A/verification.md` §2.7). The proof below is
written for N ≥ 3; it applies verbatim to N = 2 (one binary column, no non-binary column), but we do not need this."

**C2 — theorem.md §10, "Results used from earlier work".**
Old: "the stage-level lower bound (class.md Theorem 5.3, Corollary 5.5; early independent check CONFIRMED, provisional) —
for Corollary 2 only;"
New: "the stage-level lower bound (class.md Theorem 5.3, Corollary 5.5; early independent verification
`total/verification-early/verification-lower-bound.md`: CONFIRMED, final for that round, its corrections being editorial
or about the stated scope, none affecting the bound; re-derived for class.md v3 by V2) — for Corollary 2 only;"

**C3 — theorem.md Theorem 6.5, proof.**
Old: "S > 0 at the node compositions because they are within Cρ_∞^{n+1} of points of Δ (shrink r_1 so that S ≥ α_N/2 on
the ball)."
New: "S > 0 at the node compositions because they are affine in v and within Cρ_∞^{n+1} of points of Δ, where S ≥ α_N
(enlarge n_1 so that S ≥ α_N/2 there)."

**C4 — analysis.md Theorem A.8, proof** (same point as C3).
Old: "S > 0 at the junction and product compositions because they are within C ρ_∞^{n+1} of points of Δ, where S ≥ α_N
(shrink r_1 if necessary so that S ≥ α_N/2 on the ball)."
New: "S > 0 at the junction and product compositions because they are affine in v and within C ρ_∞^{n+1} of points of Δ,
where S ≥ α_N (enlarge n_1 if necessary so that S ≥ α_N/2 there)."

Each correction is marked in the file with "[V2]" and the version line of theorem.md / analysis.md records it.

### 4.2 V1's corrections for theorem.md §§2–5

V1's `verification-DL/verification.md` §7 ("Corrections for theorem.md") reads: **None** — §§2–5 transcribe D and L
faithfully; V1's two precisions (limit.md Proposition 4.5 stated under (L1)–(L4); design.md Lemma 2.5 phrased through the
regularized root equation) concern statements that theorem.md already writes correctly (Proposition 5.7 is stated under
(L1)–(L4) through §5.1; Lemma 3.4 states the regularized root via P_w). I checked both points in theorem.md and agree:
nothing to integrate. V1's dated line, appended to this folder's `PROGRESS.md` at 19:55, says the same ("NO corrections for
theorem.md §§2–5"); handled: no change.

**C5 — theorem.md, status sentence of the header** (after both verdicts were available).
Old: "Status: **complete proof, not yet independently verified** (PROTOCOL: "solved (preliminary)" once integrated in the
paper; independent verification pending)."
New: "Status: **complete proof, independently verified** (PROTOCOL §4) by V1 (`verification-DL/verification.md`: design,
limit system and §§2–5, CONFIRMED, no correction to this file) and V2 (`verification-A/verification.md`: analysis,
integrated proof, Corollary 2, N = 2, CONFIRMED WITH CORRECTIONS, minor, marked [V2])."

**C6 — theorem.md §10, first line.**
Old: "**Proved (complete proofs above, not yet independently verified):**"
New: "**Proved (complete proofs above; independently verified by V1 and V2 [V2]):**"

---

## 5. Scripts run (WSL Fedora, one thread each; numerics are evidence, never proof)

| Script | What | Result |
|---|---|---|
| `analysis/scripts/check_analysis.py` (A's, re-run unchanged) | (A1)–(A6) at 𝓛, Newton zeros of G_n with generic eigen-data, reconstruction and every balance; cases n3, n3q, n4, n4T | **ALL PASS** (45 PASS, 0 FAIL; 7 min; `runs/rerun_check_analysis.log`), numbers identical to A's log |
| same, case n5 | N = 5, α = (8, 5, 3, 1.8, 1), pattern BBBB, 80 unknowns | **ALL PASS** (20 min; `runs/rerun_check_analysis_n5.log`): σ_min(DG_∞) = 0.0026, zeros physical at n = 12, 24, 36 and mixed lengths, \|v_n − v_∞\|/ρ_∞^{n+1} ∈ [0.77, 1.8] |
| `scripts/v2_binary_prop10.py` (V2, new) | sympy: f(s) = (α−1)(s−H_1)(s−s_−), H_1 + 1 − q = αz/(α−θ), H_1 − s_− identity, αsL_T − V_TL_B = f; shooting for 6 data sets (q from −0.3 to 1.6, n_T, n_B up to 160, unequal), every balance re-computed, every x in (0, 1) | **ALL PASS** (43 s; `runs/v2_binary_prop10.log`) |
| `scripts/v2_perturbation.py` (V2, new; imports nothing from D, L, A, C) | own design (Def. 4.7, general N) vs D's exports; G_n with explicit Underwood/trace eigen-data; C¹ estimate; Newton; Lemma 6.1 bound; contraction of Φ; physical re-shooting and every unit balance | see §5.1 (`runs/v2_perturbation.log`, `runs/v2_perturbation_n3TT.log`) |

### 5.1 Results of `v2_perturbation.py`

Precision 102–175 digits. Logs: `runs/v2_perturbation.log` (cases n3B, n3BT, n3q, n3TT, n4D, n4Dlong: **ALL PASS**, 79
PASS, 0 FAIL, 9 min) and `runs/v2_perturbation_n5.log` (case n5D: **ALL PASS**, 24 min).
r_dom is the largest tested radius (3e-4, 1e-3, 1e-2, sup norm) at which all random sample points of the C¹ test have real
simple modes, and r_check := min(10⁻³, r_dom) (10⁻³ for the N = 3 cases, 3e-4 for N = 4). "Perturbative regime" means
2‖A^{−1}‖\|G_n(v_∞)\| ≤ r_check, i.e. Lemma 6.1's hypothesis ‖A^{−1}‖\|G_n(v_∞)\| ≤ r_1/2 is plausible with r_1 = r_check;
only there are (N), (S) and contraction PASS/FAIL claims (a random start whose orbit leaves the domain of real simple
modes was not in the ball of Lemma 6.1 and is reported, not counted); outside the regime the script prints INFO lines
(the theorem says nothing there).

| case | data | design vs D's export | σ_min(DG_∞), ‖A^{−1}‖_∞, ρ_∞ | C¹ ratios at n = 8, 16, 32, 64, 128 (radius 3e-4): \|G_n − G_∞\|/ρ̄^{n+1}; ‖D(G_n − G_∞)‖/((n+2)ρ̄^{n+1}) | zeros of G_n, \|v_n − v_∞\|/ρ_∞^{n+1} | physical re-shooting (worst balance; min x) |
|---|---|---|---|---|---|---|
| n3B | α = (4,2,1), z = 1/3, q = 1, δ = 0.05, η_nb = 0.025, η_bin = 0.003, BB | 4.7e-55 | 0.044, 14.8, 0.778 | 114, 17, 1.0, 0.97, 0.89; 8e3, 719, 21, 21, 19 | n = 64, 48, mixed (min 38) in regime: 0.45, 0.45, 0.12; Lemma 6.1 bound holds; Φ contracts from 3 random starts (radius 1e-3) in 9–13 steps, step ratios ≤ 0.09 | ≤ 3e-66; 6.6e-40 (n = 64); also physical at n = 8, 16, 32 (INFO) |
| n3BT | same data, pattern BT | 4.7e-55 | 0.038, 21.9, 0.777 | 63, 54, 39, 20, 5.5; 5e3, 3e3, 1.5e3, 600, 138 (decreasing) | n = 64 in regime: 2.9; Φ contracts (9 steps, step ratios ≤ 0.094) | ≤ 4e-60; 1.7e-43; physical at 16, 32, 48, mixed; at n = 8 (outside the regime) damped Newton found a non-physical zero (\|v − v_∞\| = 0.53) |
| n3q | α = (5, 2.2, 1), z = (0.2, 0.5, 0.3), q = 1.4, δ = 0.02, TB | (no export) | 0.021, 31.0, 0.840 | 52, 1.3, 0.58, 0.58, 0.58; 5.5e3, 80, 26, 20, 16 | n = 64 in regime: 0.23; Φ contracts (16–17 steps, step ratios ≤ 0.24; 1 of 3 starts at radius 1e-3 left the domain of real simple modes after one step) | ≤ 6e-71; 2.6e-36; physical at 16, 32, 48, mixed (INFO) |
| n3TT | α = (3, 1.8, 1), z = (0.3, 0.3, 0.4), q = −0.5, δ = 0.03, TT (branch φ^low ∈ (0, α_N)) | (no export) | 0.031, 31.3, 0.866 | 699, 105, 2.4, 0.41, 0.36; 3.9e4, 3.3e3, 41, 8.9, 8.0 | (all INFO: 2‖A^{−1}‖\|G_n(v_∞)\| = 1.9e-3 at n = 64) ratios 0.41, 0.35, 0.20, 0.26 at n = 64, 48, 32, mixed | physical at 32, 48, 64, mixed; **non-physical zero at n = 16** (gap table G8) |
| n4D | α = (8,4,2,1), z = 1/4, q = 1, δ = 0.05, η_nb = 0.005, η_bin = 0.003, BBB | 2.2e-55 | 0.0037, 205, 0.817 | 69, 4.7, 0.44, 0.44, 0.44 (n ≤ 64); 8e3, 331, 73, 78, 78 | (INFO) 0.56, 0.85, 0.47, 0.03 at n = 48, 32, 16, mixed (min 18) | physical at all n tested; min x 7e-65 |
| n4Dlong | same as n4D, n = 64, 80, 96 and mixed (min 42), 175 digits | 2.2e-55 | same | same; 79 at n = 96 | n = 96, 80 in regime (r_check = 3e-4): 0.33, 0.34; Lemma 6.1 bound holds; Φ contracts from 3 random starts at radius 3e-4 (8–9 steps, step ratios ≤ 0.09); n = 64, mixed (INFO): 0.38, 0.009 | ≤ 7e-111; min x 7.7e-128 (n = 96); physical also at n = 64, mixed |
| n5D | α = (16,8,4,2,1), z = 0.2, q = 0.5, δ = 0.05, η_nb = η_bin = 0.002, BBBB, 80 unknowns | 4.9e-55 | 0.0011, 715, 0.790 | at radius 3e-4: 63, 7.1, 0.58, 0.58, 0.58 (n ≤ 36); 2.9e4, 1.5e3, 365, 385, 369 | (all INFO: 2‖A^{−1}‖\|G_n(v_∞)\| = 0.14 at n = 36) 1.5, 1.5, 0.32 at n = 36, 24, mixed (min 13) | physical at n = 24, 36 and mixed (min 13); min x 1.1e-87 (n = 36) |

In every case: Lemma 4.9 holds (middle feed roots global; deviation 0, or 1e-103 = rounding for N = 5); κ_p − e_p − δ = O(η_nb) (Proposition
4.10); total balances, zero duty at junctions and side draws, reboiler duty = R = H* + δ, condenser duty = −V_{T[1,2]};
G_∞(v_∞) = 0 to working precision with the explicit Underwood/trace formulas; pivots c_{r*} > 0 at both ends of every
section (min 0.0076–0.021); half-profiles (H6) ≥ 0 for k < 300 with c_0 ∈ [0.72, 1.0]. The C¹ ratios are bounded (they
decrease) as n grows, as Lemma 3.5 (c) / Lemma 5.2 state; for N = 4 the neighborhood of real simple modes is smaller
than 10⁻³ in the sup norm (consistent with r_1 = O(η)). Where the regime holds, \|v_n − v_∞\| ≤ 2‖A^{−1}‖\|G_n(v_∞)\|
(Lemma 6.1's displacement bound) and the chord map Φ contracts to v_n from random starting points (observed step ratios
≤ 0.24 ≤ 1/2). Below the regime (small n, large ‖A^{−1}‖: n3TT, n4D, n5D) the zeros found are still physical except n3TT at
n = 16 and n3BT at n = 8 (G8). With the constants measured for n3B (‖A^{−1}‖_∞ = 14.8, C¹ ratio ≈ 20 at ρ̄ ≈ 0.78, r_1 ≈ 3e-4),
the two inequalities of Theorem 6.2 hold from n_0 ≈ 45–50 — an illustration, not a certificate.

---

## 6. Novelty scan (PROTOCOL §1; quick, AI campaigns)

- `gh api "users/ipitchford/repos?per_page=100" --paginate`: 176 repositories (latest updates 2026-10-20). Keyword
  filter (petlyuk, distill, skogestad, halvorsen, underwood, vapor, column, separation, reboil, thermally coupled,
  dividing): only "three-copy-werner-undistillability" (quantum entanglement distillation), "polydegree … column"
  (algebra), "tree-network-separation", "unique-equilibrium-lsi-separation" — none about distillation columns.
- `gh api "users/mmislan/repos?per_page=100"`: 25 repositories, none related.
- `gh search repos` "Petlyuk": boblox/Petlyuk (a 2014–2018 JavaScript web site, no proof), two personal profile repos;
  "Halvorsen Skogestad", "minimum vapor distillation", "Underwood minimum vapor", "thermally coupled distillation minimum
  energy": no hits.
- Zenodo API: "Petlyuk" (1 hit, power systems), "Halvorsen AND Skogestad AND distillation" (0), "minimum vapor/vapour AND
  distillation" (0), "Petlyuk arrangement" (0), "Underwood AND distillation AND minimum" (0), "thermally coupled AND
  distillation" (5 design/simulation records, latest "Optimal design of heat pump-assisted thermally coupled distillation
  sequences", 2026-10-05, not a proof), "dividing wall column" (2, applications).
- Two web searches: only HS's own pages (Parts 2–3, the 2002 conference paper), which state the claim as a conjecture.

**Result: nothing found that proves the conjecture or its attainment part.** (The full literature search is Agent W's.)

---

## 7. What was not checked

- §§4–5 of theorem.md in full (design, limit system, Theorem L) and the files `attainment/design/`, `attainment/limit/`:
  V1's scope; only the spot checks of §2.9.
- The scripts of `attainment/computation/` (Agent C) were not re-run.

## 8. Comparison with V1 and V1's corrections

My verdicts (§1), re-derivations, gap table and corrections C1–C4 were written before I opened anything in
`verification-DL/` (PROGRESS.md, 19:44). I then read V1's `verification.md` (in the state of about 19:51, with its
verdicts and §7 written and its script tables still to be filled in; re-read at 20:16 with its tables filled in:
verdicts unchanged); I did not open V1's PROGRESS.md.

- **Agreement.** V1 CONFIRMS Theorem L, Theorem 1.4 (= Proposition 2.1), Lemmas 2.1–2.4 (= Lemmas 3.1–3.5), Theorem D
  and the transcription in theorem.md §§2–5, with no correction to theorem.md. My spot checks of §§4–5 (Proposition 4.10 (c),
  Lemmas 5.3, 5.6, Propositions 5.7, 5.8) and my own implementation of the design (agreement with D's exports to 5e-55,
  N = 3, 4, 5) agree with V1. V1 re-derived the positivity argument of Proposition 2.1 independently with the same
  conclusion (only K = α/S > 0 is used), and notes, as I do, that the node hypothesis S > 0 is redundant except at the
  condenser, which needs none.
- **Same observation from both sides.** V1 (remark 3) found finite-stage zeros not converging at n = 35–50 in two of its
  18 cases; I found a non-physical zero at n = 16 in case n3TT (G8). Both are below the (non-explicit) n_0 of Theorem 6.2
  and are consistent with theorem.md §7–§8.5: n_0 depends on (δ, η, π) and is not small when ‖A^{−1}‖ is large.
- **Independence of the numerics.** V1 built the limit system from a general eigen-solver (mpmath `eig`); I built it from
  the explicit Underwood/trace formulas and polynomial roots; A's script uses generic eigenvectors too. All three agree.
- **No disagreement.** Nothing in V1's report changes a verdict of mine. Combined verdict on the claim of
  BRIEF-verification-common: **Theorem 1 (finite-stage attainment, every N ≥ 2) and Corollary 2 (the
  Halvorsen–Skogestad conjecture in the class 𝒜, in the finite-stage-limit sense of class.md §4) are CONFIRMED WITH
  CORRECTIONS (minor/editorial only).**
