# Independent verification V1: the forced design (Theorem D) and the limit system (Theorem 1.4, Lemmas 2.x, Theorem L)

Verifier: Claude Opus 5.5 (effort xhigh), fresh context. 2026-10-20 (times in `PROGRESS.md`). Paths are relative to
`attainment/` unless stated otherwise.

**Scope** (BRIEF-V1, BRIEF-verification-common): the design and its limit configuration (`design/design.md` v1.1,
`design/hypotheses.md` v2.1, Theorem D), the limit system and its nondegeneracy (`limit/limit.md`: Theorem 1.4,
Lemmas 1.1–1.2, 2.1–2.4, Proposition 3.3, Lemmas 3.4–3.7, Lemmas 4.1–4.4, Propositions 4.5–4.6, Theorem L,
Corollary 4.7), and their transcription in `theorem.md` §§2–5.

**Read:** `chemistry-open-problems/PROTOCOL.md`; the two briefs; `design/design.md`, `design/hypotheses.md`,
`limit/limit.md` in full; `theorem.md` in full (§§2–5 line by line, §§1, 6–10 for context); the cited sources
`total/foundations/class.md` §§2–4, `total/attainment-astra/attainment.md` §§2–3 (Lemmas 1–2, eqs. (2)–(14)),
`paper.tex` Theorem `achieve` and Proposition `cut`; the scripts that I re-ran. **Not opened:** any `PROGRESS.md` other
than mine, any `FINAL-REPORT.md`, any other `BRIEF-*.md`, the folder `verification-A/` (V2) — see §9.

---

## 0. Verdicts

| Claim | Verdict |
|---|---|
| **Theorem L** (`limit.md` §4: Lemmas 4.1–4.4, Propositions 4.5 (trace block), 4.6 (major block), Theorem L), with its inputs Proposition 3.3 (spectral structure from (L1)–(L4)), Lemmas 3.4–3.7 and Corollary 4.7 | **CONFIRMED** |
| **Theorem 1.4** (exact reformulation in (d, X)) and Lemmas 1.1–1.2 | **CONFIRMED** |
| **Lemmas 2.1–2.4** (spectrum, analytic eigen-data through zero trace flows of either sign, exact regularized trace eigenvector (2.1)–(2.2), spectral form of a section) | **CONFIRMED** |
| **Theorem D** (`design.md` §§2–6: Lemmas 2.1–2.5, 3.1–3.2, Propositions 3.3, 4.1–4.3, Lemmas 5.2–5.3, Proposition 5.4, Theorem 5.5, Theorem 6.2) and `hypotheses.md` (H1)–(H8), §4 | **CONFIRMED** (two cosmetic precisions applied, §6) |
| **Transcription** of D and L in `theorem.md` §§2–5 | **CONFIRMED** (faithful; no correction needed, §7) |

No gap of substance was found. Every step listed in BRIEF-V1 was re-derived by hand (§§1–4). An independent
re-implementation (own construction of the design by linear solves of the defining conditions, own limit system built
from a general eigen-solver instead of Underwood's formulas, finite-difference Jacobian, direct check of every physical
balance of the reconstructed finite-stage steady states) agrees with every claim on 23 cases, N = 3..6, all kinds of
forcing patterns, forcing from 5 % to 99.9 % of the explicit bound, nearly equal volatilities (§8). All the solvers'
certificate scripts in my scope were re-run (all PASS), and the exhaustive combinatorial script was made to fail on ten
deliberately broken models (§8). Two low-severity precisions (applied) and three observations are recorded in §5; none
affects a statement.

---

## 1. Theorem L: re-derivation

Notation as in `limit.md`. Throughout, (L1)–(L4) are the four hypotheses of `limit.md` §3.2.

### 1.1 Spectral structure at the limit (Proposition 3.3), checked case by case

For each of the six section types I recomputed the support of w^∞, the major roots (Lemma 2.3), the pinch, the
position of every trace mode relative to the major spectrum, its closed end, the support of the limit composition at
that end, and the signs of G^{(k)}_{x_J}(α_k) and Ψ^{(k)}_σ(α_k):

| type | supp w | lighter traces k (closed at bottom) | supp at bottom | heavier traces k (closed at top) | supp at top | where (L3) enters |
|---|---|---|---|---|---|---|
| T-nb [i,m] | [i, m−1] | k < i (above r_c, Ψ > V) | X_c: [i, m] | k ≥ m (below u_c) | X_{[i,m−1]}: [i, m−1] | u_c > α_m |
| B-nb [i,m] | [i+1, m] | k ≤ i (above v_c) | X_{[i+1,m]}: [i+1, m] | k > m (below poles, Ψ > V) | X_c: [i, m] | v_c < α_i |
| T-np [p,p+1] | {p} | k < p (Ψ > V) | X_c: {p, p+1} | k ≥ p+1 (below φ^T) | x_p: S_p ⊆ {p−1, p} | φ^T > α_{p+1} |
| B-f [p,p+1] | {p, p+1} | k < p (above φ') | x_{p+1}: S_{p+1} ⊆ {p, p+1, p+2} | k > p+1 (Ψ > V) | X_c: {p, p+1} | φ' < α_{p−1} |
| T-f [p,p+1] | {p, p+1} | k < p (Ψ > V) | X_c: {p, p+1} | k ≥ p+2 (below φ^low) | x_p: S_p ⊆ {p−1, p, p+1} | φ^low > α_{p+2} |
| B-np [p,p+1] | {p+1} | k ≤ p (above φ^B) | x_{p+1}: S_{p+1} ⊆ {p+1, p+2} | k > p+1 (Ψ > V) | X_c: {p, p+1} | φ^B < α_p |

In every row the closed-end support lies strictly on the correct side of k (π_p = B excludes p+1 from S_p; π_p = T
excludes p from S_{p+1}), so x^∞_{J,k} = 0 and G^{(k)}_{x^∞_J}(α_k) is a nonempty sum of terms of one sign (< 0 at a
bottom end, > 0 at a top end), and Ψ^{(k)}_σ(α_k) = Ψ_σ(α_k) > 0 (either α_k lies beyond all poles on the side where
Ψ > V, or between ψ_low/ψ_high and the poles on the side where Ψ > 0 — exactly where (L3) is used). Hence κ > 0 for every
trace row. Pivots (Proposition 3.3 (d)): own non-binary junction — G_{X_c}(t)Π_{k∈[i,m]}(α_k − t) has degree exactly
m − i and the m − i simple zeros r_c, θ's, s_c, so G_{X_c}(u_c) ≠ 0 ⇔ u_c ≠ s_c and G_{X_c}(v_c) ≠ 0 ⇔ v_c ≠ r_c (L4);
own binary junction — pinning section: ℓ_{r*}·e_{r*} = −(S(e)/V)Ψ'(r*) ≠ 0; non-pinning section: by the pinch identity
G_{X_c}(φ^T) ∝ Ψ_{B_c}(φ^T) ≠ 0 ⇔ φ^T ≠ φ'' (L4) (mirror for T-pinning); child junctions and product nodes: one-signed
sums. All confirmed. The equivalence stated after (L4) (r_c ≠ v_c ⇔ Ψ_{T_c}(θ^c_i) ≠ 0, etc.) is also correct: a common
root in (α_{i+1}, α_i) is a root of Φ_c, hence equals θ^c_i.

### 1.2 Coordinates, counts, zero blocks (Lemmas 4.1–4.4)

Table 2 / `limit.md` §4.1: per column, 2(N−1) coordinates; I recounted: non-binary
2(N − m + i − 1) + 2 trace + (m − i − 2) + (m − i) major = 2(N − 1); B-pinned (N − 1) + (N − 2) + 1; T-pinned
(N − 2) + 1 + (N − 2) + 1. Per column and component, #trace rows = #trace coordinates (k ∉ range: 2 = 2; non-binary
k = i: B_c only, b_{c,i}; k = m: T_c only, d_{c,m}; B-pinned k = p+1: T_c only, d_{c,p+1}; T-pinned k = p: B_c only,
b_{c,p}); #major rows = 2(m − i − 1) (r_c, θ^T's, θ^B's, s_c) = #major coordinates; binary 1 = 1. Totals N(N − 1)².

Lemma 4.2 (trace quantities are linear forms in the trace coordinates of their own component): checked by induction,
parents first. The only non-obvious entries are b_{c,k} for k ∉ range(c), equal to d_{[i,m+1],k} + b_{[i−1,m],k} −
d_{c,k}, where both parental limit flows vanish because a non-binary parent's limit flows are supported on [i, m]
(both parents of any column are non-binary). Product contents x_{p,k}, k ∉ S_p: S_p is exactly the union of the
supports of b^∞_{[p−1,p]} and d^∞_{[p,p+1]}.

Lemma 4.3: by (2.2), the numerator of a trace row is x_{J,k} + w_{σ,k}·g(v) with g analytic; both factors vanish at
v_∞, so its derivative is δx_{J,k} + g(v_∞)δw_{σ,k} and the pivot (and every major quantity) enters only multiplied by
zero. Signs: q = w at a bottom end where G^{(k)} < 0, q = −w at a top end where G^{(k)} > 0, so g δw = −κ δq in both
cases. Hence the lower-left block (trace rows × major coordinates) is 0 and the trace part is block diagonal over the
components. Lemma 4.4 (node identities) holds identically on 𝒱 (it is the definition of b_c and of x_p), hence also
for variations.

### 1.3 Trace block (Proposition 4.5)

Reduction: with δt ∈ ker 𝕋_k extended by 0, each trace row gives ξ_{c(σ)} = κ_σ δq_σ. At J ∈ 𝒥_k every attached section
carries no k in the limit (Lemma 3.4: the four attached sections of J[i,m] have limit supports in [i,m]; those of
J[p,p] have supports whose union is S_p), so each attached σ has the trace mode α_k, and J is either its closed end
(contribution ξ_J/κ_σ) or its open end (contribution −ξ_{c(σ)}/κ_σ with c(σ) ∈ 𝒥_k by 1.1). This is A_k ξ = 0 with the
stated A_k. Triangularity: lighter trace of σ ⊂ column [i,m] ⇒ closed end bottom, h(top) = h(bottom) − 1 (top of
T[i,m] is J[i,m−1], bottom of B[i,m] is J[i+1,m]) and both ends, if trace nodes, have first index > k; heavier trace ⇒
closed end top, h(bottom) = h(top) + 1, both ends with second index < k. A node of 𝒥_k cannot satisfy both
(first index > k and second index < k is impossible since first ≤ second), so 𝒥_k = 𝒥^<_k ⊔ 𝒥^>_k and no entry links the
two classes; with 𝒥^<_k by increasing h and then 𝒥^>_k by decreasing h, every off-diagonal entry (row o(σ), column c(σ))
is above the diagonal. Diagonal: P_J > 0 at product nodes; J[i,m] ∈ 𝒥^<_k is the closed (bottom) end of T[i,m], whose
support starts at i > k for all three T types; J[i,m] ∈ 𝒥^>_k is the closed (top) end of B[i,m], whose support ends at
m < k for all three B types. So det A_k = Π(A_k)_{JJ} > 0. Conclusion: ξ = 0 ⇒ δq_σ = 0 ⇒ every trace coordinate of k
(each is some ξ_J or ±δq_σ) vanishes; 𝕋_k is square (1.2) and injective. **Confirmed**, for every N and every pattern
(no case distinction depends on N; the combinatorial script re-checks it mechanically for N ≤ 12, §8).

Precision (low severity, applied, §6): the statement of Proposition 4.5 says "Assume (L1)–(L3)", but the rows (2.4)
are quotients by pivots, which are nonzero only under (L4) (e.g. the trace rows of T-nb closed at J_c have pivot
G_{X_c}(u_c)). The proof only uses the numerators, whose derivative matrix is nonsingular under (L1)–(L3); Theorem L
assumes (L1)–(L4) anyway.

### 1.4 Major block (Proposition 4.6)

Triangularity: the major rows of c are closed at J_c and depend on X_c, d_c and f_c only; f_c depends on the ancestors.
(My numerical check is stronger: the major rows of c have zero derivative in the major coordinates of every column that
is not c or an ancestor of c, §8.) Non-binary diagonal block: with all other coordinates fixed, δd_c is supported on
{i+1, …, m−1} with zero sum, δw_{T_c} = δw_{B_c} = δd_c, and Ψ'(χ)δχ = R(χ) := Σ_p α_pδd_p/(α_p − χ) for every root of
either section. The two rows of a middle root θ_j vanish at v_∞ (common root), so after multiplying by their (different)
pivots they read ℓ_θ·δX + G'_X(θ)δθ^T = 0 = ℓ_θ·δX + G'_X(θ)δθ^B; subtracting, G'_X(θ_j)R(θ_j)(1/Ψ'_T − 1/Ψ'_B) = 0 with
Ψ'_T < 0 < Ψ'_B and G'_X(θ_j) ≠ 0 (simple zero, degree count of 1.1). With R(0) = Σδd = 0, R(t)Π_{p=i+1}^{m−1}(α_p − t)
has degree ≤ m − i − 2 and m − i − 1 distinct zeros (θ's > 0), so δd_c = 0; then δχ = 0, and the m − i distinct
equations ℓ_χ·δX = 0 (χ = r_c, θ's, s_c) with ΣδX = 0 force δX = 0 (homogeneous interpolation: with Σy = 0,
Σ α_k y_k/(α_k − χ) = χ Σ y_k/(α_k − χ), and Σ y_k/(α_k − t)·Π(α_k − t) has degree ≤ m − i − 1 because its leading
coefficient is ±Σy = 0; this does not even need the interlacing). Binary block: the single row (φ'' resp. φ^mid) varies
by (α_{p+1}/(α_{p+1} − φ) − α_p/(α_p − φ))δX_{c,p+1}/pivot with φ ∈ (α_{p+1}, α_p): a negative bracket. **Confirmed.**

### 1.5 Theorem L and Corollary 4.7

DG_∞(v_∞) = [[𝕄, ∗], [0, 𝕋]] with square diagonal blocks (1.2), 𝕄 block lower triangular with nonsingular diagonal
blocks, 𝕋 = ⊕𝕋_k nonsingular: det ≠ 0. G_∞(v_∞) = 0 (Lemma 3.5): non-pinch major rows vanish by the definition of X_c and
the common middle roots (non-binary) or by left/right orthogonality (binary); trace rows vanish since x^∞_{J,k} = 0 and
w^∞_{σ,k} = 0. Lemma 3.6 (domain) and Lemma 3.7 (|G_n − G_∞| ≤ Cρ̄^{n+1}, ‖D(G_n − G_∞)‖ ≤ C(n+1)ρ̄^n) follow from
continuity of the simple spectra and pivots and from D(μ^{n+1}g) = (n+1)μ^n gDμ + μ^{n+1}Dg. Corollary 4.7 is the
standard chord contraction Φ(v) = v − DG_∞(v_∞)^{−1}G_n(v): ‖DΦ‖ ≤ K·1/(2K) = 1/2 on the ball, and the ball is
mapped into itself since ‖Φ(v) − v_∞‖ ≤ ε_0/2 + K‖G_n(v_∞)‖ ≤ ε_0. All constants depend on the fixed (δ, η, π); no uniformity is claimed or
needed (the order of limits n → ∞, then η → 0, then δ → 0 is respected in `theorem.md` §7). **Confirmed.**

---

## 2. Theorem 1.4 and Lemmas 1.x, 2.x

- **Lemma 1.1.** (i) ⇔ (ii) is the cross-section balance V y(x_j) − L x_{j+1} = w solved for x_{j+1}, using
  w S(x) = wαᵀx; 1ᵀM = (V − Σw)αᵀ = Lαᵀ, so Σx_{j+1} = 1 for any x_j with S(x_j) ≠ 0; the stage balance is the
  difference of two consecutive cross-section balances. Cross-section indexing: vapor y(x_j) rises from stage j
  (stage 0 = bottom node, an equilibrium stage) and liquid x_{j+1} falls from stage j+1 (stage n+1 = top node). ✓
- **Lemma 1.2.** 1ᵀM^{k+1}x_bot = L S(x_k)·1ᵀM^k x_bot, so 1ᵀM^{n+1}x_bot = L^{n+1}ΠS(x_k) ≠ 0, and Σx_top = 1 fixes
  c = 1/(1ᵀM^{n+1}x_bot) (> 0 when every S(x_k) > 0, as in Theorem 1.4 (b)). Σx_top = 1 holds at junctions by the constraint on X_c and at product nodes because
  Σb_c = L_{B_c} − V_{B_c} for every column (consistency of F, by induction parents first) and the side-draw balance. ✓
- **Theorem 1.4.** (a) At a steady state all compositions are in the simplex, the net flow w_σ is constant along a
  section (stage component balances), b_c = f_c − d_c is the junction component balance and x_p = (b + d)/P_p the
  product balance. (b) I re-derived every unit balance: junction J_c (four attached sections with net upward flows d_c,
  −b_c, d_{[i,m+1]}, −b_{[i−1,m]}, plus z at J[1,N]), side draw, total condenser (top cross-section of T[1,2] gives
  V y_n = (L + P_1)x_1 = V x_1: the composition rule holds), partial reboiler. Positivity: for fixed K-values
  α_i/S(x_u) > 0, the component-i balances form a column diagonally dominant Z-matrix whose column excess is the flow sent
  to product outlets, every unit being linked to a product by positive flows and reached from J[1,N] (from J[i,m] vapor
  goes to J[i,m−1] and J[i−1,m], liquid to J[i+1,m] and J[i,m+1]); the unique solution is positive. The hypothesis
  S(X_c), S(x_p) > 0 at the nodes is in fact implied by S(x_0) > 0 for every section, except at the condenser, which has
  no vapor outflow and needs no K-value (harmless redundancy). The steady state is one in the sense of class.md §3
  ((E) stages with any number of outflow streams per phase, (C) total condenser, all stream flows positive). ✓
- **Lemma 2.1.** Rank-one determinant; left eigenvector equation ℓ_i(α_i − t) = (ℓ·w)α_i/V; pinch identity
  Ψ(s) = −(V/S(e))(s − t)G_e(s) from w_i S(e) = V(α_i − t)e_i (it does not use Σe = 1). ✓
- **Lemma 2.2.** Simple real roots of P_w are analytic in w; for a trace mode, P_w = (α_k − t)Q_w − α_k w_kΠ_{j≠k}(α_j − t)
  gives (2.1) exactly, Q_{w⁰}(α_k) ≠ 0 by simplicity, and ℓ̃_r, ẽ_r are analytic for either sign of w_k: for w_k ≠ 0,
  ℓ̃_r = ((α_k − t_r)/α_k)ℓ_{t_r} (k-entry 1, other entries by (2.1)); for w_k = 0, ℓ̃_r = e_k. (2.2) is the definition
  written out. Distinct moduli: all modes are positive and distinct (Proposition 3.3 (a)). ✓
- **Lemma 2.3** (one-signed flows) ✓. **Lemma 2.4**: (E) with c ≠ 0 ⇔ ℓ_r·x_top = cλ_r^{n+1}ℓ_r·x_bot for all r, and the
  ℓ_r form a basis; any analytic nonvanishing normalization of the ℓ_r is allowed (ℓ̃ has k-entry 1). ✓

---

## 3. Theorem D

- **Lemma 3.1** (preferred split relative to the actual feed): partial fractions; sign counts (p − i negative factors in
  numerator and denominator of ρ_p; p − i above and p − i − 1 below for the residues of Ψ⁰_B); V⁰_B ≠ 0 since b⁰_m = f_m;
  L⁰ = Ψ⁰(0) > 0. ✓ (It agrees with paper.tex Theorem `achieve`.)
- **Lemma 3.2** (Astra's direction): β ∈ (α_{m−1}, α_{i+1}) by comparing θ^c_r with α_{r+1}, α_r; K(0) = 1 ⇒ Σh = 0;
  K(θ^c_i), K(α_i) > 0 (all factors positive), K(θ^c_{m−1}), K(α_m) > 0 (m − i − 1 negative factors above and below). ✓
  My script obtains h independently by solving Σh = 0, K(θ^c_r) = 0 (middle r) — same vector.
- **Proposition 3.3** (Astra's conditions from validity (V) alone): Ψ_T = Ψ⁰_T + ηK, Ψ_B = Ψ⁰_B + ηK (the equal increments
  η of all four phase flows and Σh = 0); the extreme feed roots give Ψ(θ^c_i) = ηK(θ^c_i) > 0 and Ψ(θ^c_{m−1}) > 0, whence
  the interlacing by monotonicity; Ψ_T(α_m) > 0, Ψ_B(α_i) > 0; pinch vectors positive (Lemma 2.3); cone weights from the
  sign of G at s_c (resp. r_c) by the pinch identity: G_{ê_(m)}(s_c) < 0 < G_{ê_u}(s_c), G_{ê_v}(r_c) < 0 < G_{ê_(i)}(r_c).
  The closed forms of B_c and C_c are correct. **The reduction of Astra's conditions to positivity of the deformed split
  is correct; nothing else is needed.** ✓
- **Propositions 4.1–4.3** (binary): the roots φ⁰_T = α_p(V_T − f_p)/V_T, φ⁰_B = α_{p+1}(V_B + f_{p+1})/V_B; κ > 0 ⇔
  φ⁰_B < θ' < φ⁰_T; B-pinning adds α_pη/(α_p − t) to Ψ⁰_B: positive on (α_{p+1}, α_p) ⇒ φ'' < φ⁰_B and Ψ_B(φ^T) >
  Ψ⁰_B(φ^T) > κ (Ψ⁰_B increasing, φ^T ≥ φ⁰_T > θ'); negative on (α_p, ∞) ⇒ φ' < α_{p−1} ⇔ η < (α_{p−1} − α_p)Ψ⁰_B(α_{p−1})/α_p,
  with Ψ⁰_B(α_{p−1}) > κ > 0. T-pinning is the mirror image; for p + 1 ≤ N − 1, L_T = Ψ_T(0) > Ψ_T(α_{p+2}) > 0 uses only
  the monotonicity of Ψ_T (w ≥ 0), not L_T > 0, so it is not circular; for p + 1 = N the bound η < V_T − f_p is needed and
  present. Cone decompositions (e), explicit B' and the limits (f) ✓ (in (f) of Prop. 4.2 the p-entry of the unforced pinch
  is Ψ⁰_B(α_p)/(V_B + f_{p+1}) ∈ (0, 1)). **The forcing bounds and "forcing only increases the margin" are correct.**
- **Lemma 5.2** (balances, flux bookkeeping): every junction and side-draw balance (vapor and liquid) re-derived; the
  vapor cascade (5.1) telescopes; identity (e) holds at every unit (also checked numerically at random t, §8). ✓
- **Lemma 5.3** (middle roots global): induction; the parents' section functions vanish at θ_{i+1..m−1} resp.
  θ_{i..m−2}; Φ_c has one root per interval. ✓
- **Proposition 5.4.** (a), (b) ✓ (V_{T[p,p+1]} = R + const, binary feed independent of R). (c) **light-side cut**,
  re-derived: at η_nb = η_bin = 0, n_r = z_r e_r; summing the flux identities at t = θ_p over the light units
  {J[i,m] : m ≤ p} ∪ {P_1, …, P_p}, a section with exactly one light end is T[i,p+1] (i ≤ p) and nothing else (a B
  section has both ends with the same second index; a T section has second indices m and m − 1); with
  V_{T[1,2]} = R + 1 − q (sum over all units, t → ∞) and Ψ_{T[i,p+1]}(θ_p) = 0 for i < p, κ_p(0; δ) = R + 1 − q −
  Σ_{r≤p}α_rz_r/(α_r − θ_p) = R − H_p = e_p + δ. ✓ (d) compactness/analyticity: the recursion is analytic in η_nb on a
  neighborhood of 0 (simple feed roots of positive feeds, rational expressions); (V) is strict at 0; on a cube [0, η_1]^n
  all columns are valid and κ_p(·; 0) is C¹; the mean value theorem with C_1 = max‖∇κ_p‖_1 gives the bound
  ‖κ_p(η; 0) − e_p‖ ≤ C_1 max_c η_c. **The argument is complete.** What is uniform: η_1 and C_1 depend only on (α, z, q)
  (not on δ, η_bin or π); η̄(δ) = min(η_1, δ/(2C_1)) is linear in δ; η̂_p depends on (δ, η_nb) through the binary
  vapors; nothing else is uniform, and nothing else is claimed. η_1, C_1 are not explicit, but admissibility of a given
  parameter is a finite list of explicit inequalities (i)–(iii), which is how the scripts (D's and mine) choose parameters.
- **Theorem 5.5** ✓. **Theorem 6.2**: (H1)–(H4) collect the above; (H5) and the limit relations type by type (§1.1 table);
  pivots: the cone coefficients B_c, C_c, B', 1, 1, C' at the own junction, one-signed sums with ℓ_{r*}·ê_{r*} of the same
  sign at the other end; (H6): going up a B section M ≥ 0, going down a T section M^{−1} ≥ 0 (Sherman–Morrison, re-derived),
  the two other directions start at X_c with a nonnegative cone decomposition (types 1, 2, 3, 6: two positive pinch
  vectors with positive weights; types 4, 5: X_c is itself the pinch), and 1·v_k > 0 with 1·v_k → 1; (H7) (i)–(v) ✓ (table of §1.1); (H8) impurities
  (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p ✓. **The pinch rule** (top sections at the smallest major root, bottom sections at
  the largest) is consistent: with it every trace mode is outside the major spectrum on the side that closes it at a node
  where its component is absent, and all pivots are positive. ✓

---

## 4. Transcription in `theorem.md` §§2–5

Compared statement by statement with `design.md` and `limit.md`: §2 (model, arrangement, specification, unknowns,
Proposition 2.1 = Theorem 1.4 + positivity), §3 (Lemmas 3.1–3.5 = design Lemmas 2.1–2.5, limit Lemmas 2.1–2.4; product
form (3.2) = limit (2.3) multiplied by the pivots, with Lemma 3.5 (d) showing the equivalence of the two normalizations
at a zero), §4 (Lemmas 4.1–4.2, Proposition 4.3, Propositions 4.4–4.6, Definition 4.7, Lemmas 4.8–4.9, Proposition 4.10
with the cut proof, Theorem 4.11, Definition 4.12, Table 1, Theorem 4.13), §5 (hypotheses (L1)–(L4), 𝒥_k, Lemmas
5.1–5.6, Table 2, Propositions 5.7–5.8, Theorem 5.9). Every formula, index range, sign and hypothesis agrees with the
sources (checked: the node balance signs in the proof of Proposition 2.1, P_w(α_k), the root intervals in Propositions
4.5–4.6 including p = 1 and p + 1 = N, the vapor cascade, the cut sum, Table 2 row by row, the A_k entries). Harmless
differences: §2.2 takes n_σ ≥ 1 (limit.md allows n_σ ≥ 0); §5 states the derivative bound as C(n + 2)ρ̄^{n+1} instead of
C(n + 1)ρ̄^n (equivalent). **No correction needed** (§7).

---

## 5. Gaps and remarks (severity, fix)

No gap of substance. Remarks:

1. **(low; precision of a hypothesis)** `limit.md` Proposition 4.5 says "Assume (L1)–(L3)" although the rows it
   differentiates are quotients by pivots, nonzero only under (L4). Fix: applied (§6, C1). No effect on Theorem L, which
   assumes (L1)–(L4).
2. **(cosmetic)** `design.md` Lemma 2.5 says "Ψ_σ has a unique root φ_k(w) near α_k" for |w_k| small; at w_k = 0 the
   point α_k is not a root of Ψ_σ but of the regularized equation F(φ, w) = 0 of the proof (the trace mode). Fix: applied
   (§6, C2). `limit.md` Lemma 2.2 (b) and `theorem.md` Lemma 3.4 state it correctly (roots of P_w).
3. **(observation, not a gap)** The constants n_0 of Corollary 4.7 can be large. In three of my cases (forcing at 90 % of
   the explicit bound with ρ_σ up to 0.88; a case with σ_min(DG_∞) = 0.018; a nearly-equal-volatility case with
   ρ_σ = 0.97) the finite-stage check at n = 35–50 (resp. 120) did not converge from v_∞, while it does for larger n
   (100–160) (§8.3). This is consistent with the theorem (n_0 is not explicit and depends on (δ, η, π)); `theorem.md` §7
   correctly takes n_0 after fixing (δ, η). For the same reason, trace contents of order ρ^n (here down to 1e-351) need
   a working precision well beyond the forward amplification to be resolved in sign (§8.3, last paragraph).
4. **(observation)** `limit.md` §3.2 refers to `hypotheses.md` "v1, 18:07"; the current file is v2.1 (18:47), whose
   (H1)–(H8) are unchanged in substance (its change log says so, and I compared the statements). No fix needed.
5. **(observation)** L's own test bed (`limit/scripts/limitsys.py`, used by `check_limit_jacobian.py`) is a different
   design (its own R), as `limit.md` §8 states; some of its instances have R < max_j H_j with impure products. This is
   legitimate for testing Theorem L, which assumes only (L1)–(L4), and is not D's design; D's exported configurations are
   checked separately (`check_D_configs.py`, re-run: 6/6 PASS) and by my own construction.

---

## 6. Corrections applied (in `limit/` and `design/`)

**C1 — `limit/limit.md`, Proposition 4.5 (trace block), first sentence.**
Old: "**Proposition 4.5 (trace block).** Assume (L1)–(L3). For every k, 𝕋_k is nonsingular."
New: "**Proposition 4.5 (trace block).** Assume (L1)–(L4) ((L4) only makes the pivots, hence the rows (2.4), well
defined; the argument below uses only the numerators of the rows and (L1)–(L3)) [precision added by the independent
verifier V1, 2026-10-20]. For every k, 𝕋_k is nonsingular."

**C2 — `design/design.md`, Lemma 2.5, first sentence of the statement.**
Old: "Then for |w_k| small Ψ_σ has a unique root φ_k(w) near α_k, real-analytic in w, with
α_k − φ_k = α_k w_k/Ψ^{(k)}_σ(φ_k), and"
New: "Then for |w_k| small the equation (α_k − φ)Ψ^{(k)}_σ(φ) = α_k w_k has a unique solution φ_k(w) near α_k,
real-analytic in w (φ_k = α_k, the trace mode, at w_k = 0; for w_k ≠ 0 it is the root of Ψ_σ near α_k) [precision
added by the independent verifier V1, 2026-10-20], with α_k − φ_k = α_k w_k/Ψ^{(k)}_σ(φ_k), and"

Neither correction changes a statement used elsewhere (Theorem L assumes (L1)–(L4); Lemma 2.5 is used only through its
formulas, which are unchanged). For the manuscript (agent W): no change of substance is needed; if Proposition 4.5 /
`theorem.md` Proposition 5.7 is reproduced with explicit hypotheses, state them as (L1)–(L4).

---

## 7. Corrections for theorem.md

None. `theorem.md` §§2–5 transcribe D and L faithfully (§4); the two precisions of §6 concern statements that
`theorem.md` already writes correctly (Proposition 5.7 is stated under (L1)–(L4) through §5.1; Lemma 3.4 states the
regularized root via P_w). A dated line saying so was appended to V2's `PROGRESS.md` (see §9).

---

## 8. Scripts run and results

All runs under WSL Fedora (`outils/wsl` from the project root), Python 3.12, mpmath 1.3, one thread per process, at most
three processes. Logs in `verification-DL/runs/`.

### 8.1 Solvers' certificate scripts (re-run)

| script (folder) | arguments | result | time | log |
|---|---|---|---|---|
| `check_combinatorics.py` (limit) | `12` | ALL PASS: N = 3..12, every pattern (4092 configurations, up to 1452 rows) | 10 s | `rerun_check_combinatorics.log` |
| `check_limit_jacobian.py` (limit) | `4 6 2` | 16/16 instances (N = 3..6, random patterns) C1–C6 PASS: G_∞(v_∞) ≤ 1e-57, zero blocks exactly 0, trace rows as predicted, A_k M-matrices, σ_min(DG_∞) from 7.9e-5 to 9.1e-2 | 252 s | `rerun_check_limit_jacobian.log` |
| `check_D_configs.py` (limit) | — | 6/6 PASS on D's exports (n3_N, n3_B, n4_1, n5_1, N = 6, N = 7 with 252 unknowns); σ_min(DG_∞) = 9.5e-3, 7.8e-3, 3.5e-3, 2.0e-3, 2.3e-4, 9.0e-5 | 195 s | `rerun_check_D_configs.log` |
| `check_unforced.py` (limit) | — | ALL PASS: σ_min(DG_∞) = 0 exactly at η_p = 0 (null vector on a junction-ratio coordinate); ∝ η_p in the plain normalization; ≈ 0.038 in the regularized one | 21 s | `rerun_check_unforced.log` |
| `check_finite_n.py` (limit) | `N3 80,60,40,30,25` | ALL PASS: ‖v_n − v_∞‖ = 4.2e-7 … 1.3e-2, unit balances ≤ 1.4e-71, all x > 0, side impurity → 0.01768 | 11 s | `rerun_check_finite_n_N3.log` |
| `check_design.py` (design) | `4 50 2` | ALL PASS, 80 instances N = 3..7 (generic, small z, low/high q, nearly equal α, equal peaks; random patterns): every family of (H1)–(H8), flux bookkeeping, κ identities 80/80, stage-level bound at the limit 320/320; min κ_p/δ = 0.5003, max ρ = 0.9982 | 157 s | `rerun_check_design_4_50_2.log` |
| `crosscheck_N3.py` (design) | — | PASS: N's flows = F(0.05, (0.025; 0.003, 0.003), BT) to 9.4e-17; limit ξ = 0.6961542454, ζ = 0.5539306858, side impurity 0.017681729 = N's finite-stage limits | < 1 s | `rerun_crosscheck_N3.log` |
| `eta_scaling.py` (design) | — | ALL PASS: κ_p(0; 0) − e_p ≤ 3.7e-40; ‖κ_p(η·1; 0) − e_p‖/η → 3.2–4.2; largest admissible common η_nb = (0.12–0.57)·δ, constant for δ = 1e-2 … 1e-5 (linear η̄(δ)) | 65 s | `rerun_eta_scaling.log` |
| `export_config.py` (design) | `n4_1` (output to my folder) | ALL PASS (JSON written to `runs/rerun_config_n4_1.json`) | < 1 s | `rerun_export_n4_1.log` |

The solvers' scripts were run unmodified; their outputs were redirected to my folder (no file of theirs was
overwritten).

### 8.2 Mutation test of the combinatorial checker (`scripts/mutate_combinatorics.py`, own)

Imports the unmodified `limit/scripts/check_combinatorics.py`, replaces its model functions by broken versions and runs
all patterns for N = 3..7 (124 configurations). Unmodified: 0 flagged. Every mutation is detected:
product supports without designed impurities (124/124 flagged); unforced binary columns with forced bookkeeping
(124/124); non-binary T sections carrying the heavy key (124/124); non-binary B sections carrying the light key
(124/124); product supports from the reversed pattern (124/124); closed-end rule reversed (124/124); wrong top end of T
sections (crash, KeyError); **unforced binary columns with pure products — the singular design of `limit.md` §7.1**
(124/124, caught by the per-column count K3: one surplus trace row); closed-end rule reversed in B sections only
(124/124); long T sections skipping a level (120/124; the four N = 3 configurations have no T section of length ≥ 3).
Log `runs/mutate_combinatorics.log`: ALL PASS.

### 8.3 Independent re-implementation (`scripts/v1_independent.py`, own; batch `scripts/run_batch.sh`)

Nothing is imported from the solvers' code. Design: the preferred split by solving Ψ_T(θ^c_r) = 0 (all feed roots,
d_i = f_i), Astra's h by solving Σh = 0, K(θ^c_r) = 0 (middle roots), binary columns by (F2), products by (F3). Limit
configuration: X_c by solving the interpolation system; binary X_c = eigenvector of the forced pinch from a general
eigen-solver. Limit system: for every section, eigenvalues and left/right eigenvectors of M_σ = V diag(α) − wαᵀ from
`mpmath.eig` (not Underwood's formulas), left eigenvectors normalized by a fixed entry (the k-entry for a trace mode, i.e.
the regularized normalization), pinch by the rule; rows ℓ_r·x_c/ℓ_{r*}·x_c (− τ-terms for G_n); coordinates of Table 2;
Jacobian by central differences (step 10^{−dps/3}). Checks per case: D1–D9 (Theorem D: flows and balances of the
specification with zero duty at junctions and side draws, R = H* + δ; supports; middle roots global; interlacing,
Ψ_T(α_m), Ψ_B(α_i) > 0; κ_p ≥ e_p + δ/2, forcing increases the margin, (L3), (L4); κ_p(η; δ) − κ_p(η; 0) = δ and
κ_p(0; δ) = e_p + δ; flux bookkeeping at every unit; simple positive spectra, pinch rule, limit relations, pivots > 0 in the
sum-normalized eigenbasis, trace signs and κ > 0; half-profiles (chart vectors) ≥ 0 for k ≤ 400; product impurities),
L0–L5 (Theorem L: coordinates reproduce 𝓛; G_∞(v_∞) = 0; zero blocks; every trace row equal to Lemma 4.3's prediction;
major rows of c independent of non-ancestor major coordinates; diagonal blocks nonsingular; every 𝕋_k square and
nonsingular; A_k upper triangular Z-matrix with positive diagonal in the h-order; σ_min(DG_∞) > 0), and F1 (finite n: zero
of G_n by the chord iteration of Corollary 4.7 (Jacobian refreshed if slow), working precision raised by
(n+1)log10(λ_max/λ_min) digits; reconstruction by the recursion; **every component balance of every stage, junction, side
draw, condenser and reboiler, the condenser composition rule, x_top reached, all mole fractions > 0**; |v_n − v_∞|
increasing as n decreases).

Cases (`run_batch.sh`, `run_batch2.sh`; logs `runs/v1_<case>.log`, summaries `runs/batch_summary.log`,
`runs/batch2_summary.log`). "n3" = α (4, 2, 1), z = 1/3, q = 1, δ = 0.05, η_[1,3] = 0.025, η_bin = 0.003 (network.md
§11.3 data); "n4" = α (8, 4, 2, 1), z = 1/4, q = 1, δ = 0.05, η_nb = 0.005, η_bin = 0.003; "rand:N:s" = random α
(consecutive ratios 1.3–2.5), z, q ∈ [−0.5, 1.5], δ = 0.05 H*, η_nb halved from δ/10 until admissible, η_p = frac·η̂_p
(frac 0.25 unless stated); "near" = rand with one consecutive volatility ratio 1.05.

| case | pattern | frac | D1–D9, L0–L5 | max ρ_σ | σ_min(DG_∞) | F1: n → ‖v_n − v_∞‖ (all balances ≤ 6e-86, all x > 0) |
|---|---|---|---|---|---|---|
| n3 | BB | — | PASS | 0.778 | 0.81 | 60 → 1.0e-7; 40 → 1.6e-5; 30 → 1.8e-4 |
| n3 | BT | — | PASS | 0.777 | 0.31 | 60 → 6.4e-7; 40 → 1.2e-4; 30 → 1.7e-3 (side impurity 0.0177 = network.md §11.3) |
| n3 | TB | — | PASS | 0.778 | 1.25 | 60 → 1.1e-7; 40 → 1.6e-5; 30 → 1.9e-4 |
| n3 | TT | — | PASS | 0.777 | 0.31 | 60 → 6.6e-7; 40 → 1.2e-4; 30 → 1.7e-3 |
| n4 | BBB | — | PASS | 0.817 | 0.17 | 50 → 1.8e-5; 35 → 5.6e-4 |
| n4 | TTT | — | PASS | 0.817 | 0.15 | 50 → 3.9e-5; 35 → 9.1e-4 |
| n4 | TBT | — | PASS | 0.817 | 0.13 | 50 → 3.7e-5; 35 → 8.8e-4 |
| rand:4:1 | BBT | 0.25 | PASS | 0.816 | 0.11 | 50 → 1.7e-4; 35 → 3.5e-3 |
| rand:4:2 | BBB | 0.25 | PASS | 0.898 | 0.13 | 50 → 1.0e-3; 35 → 4.9e-3 |
| rand:4:3 | TTB | 0.9 | PASS | 0.911 | 0.10 | 50 → 1.4e-3; 35 → 4.3e-3 |
| rand:4:4 | BBT | 0.9 | PASS | 0.880 | 0.14 | 150 → 1.3e-8; 100 → 8.1e-6 (n = 50, 35: no convergence from v_∞, first batch) |
| rand:4:5 | BBB | 0.999 | PASS | 0.905 | 0.17 | 120 → 4.8e-7 |
| rand:4:6 | TTT | 0.999 | PASS | 0.778 | 0.094 | 120 → 2.9e-14 |
| near:3:3 | TB | 0.25 | PASS | 0.970 | 1.12 | 160 → 2.3e-3 (8 Newton refreshes); 120: iteration leaves the domain (n < n_0) |
| near:4:1 | BBT | 0.25 | PASS | 0.976 | 0.11 | 120 → 8.4e-3; 90 → 1.3e-2 |
| rand:5:1 | BBTB | 0.25 | PASS | 0.937 | 0.18 | 45 → 5.7e-4 |
| rand:5:2 | TBTB | 0.25 | PASS | 0.878 | 0.018 | 150 → 2.0e-8 (at 449 digits, see below); 110 → 3.5e-6 (n = 45: no convergence, first batch) |
| rand:5:3 | BBBB | 0.9 | PASS | 0.905 | 0.064 | — |
| rand:5:4 | TTTT | 0.05 | PASS | 0.935 | 0.089 | — |
| near:5:2 | BBBT | 0.25 | PASS | 0.980 | 0.013 | — |
| rand:6:1 | BBTBT | 0.25 | PASS | 0.874 | 0.025 | — |
| rand:6:2 | TTTTT | 0.25 | PASS | 0.902 | 0.040 | — |
| rand:6:3 | BTBTB | 0.9 | PASS | 0.901 | 0.037 | — |

In every case (60 digits; the three N = 6 cases at 50 digits): G_∞(v_∞) ≤ 1e-58 (N = 6: ≤ 1.1e-48); the zero blocks of
Theorem L vanish to the finite-difference noise level (≤ 4e-39 relative to entries up to 3e3), so they are exactly
zero; every trace row equals the prediction of Lemma 4.3 to ≤ 4e-30 (N = 6: ≤ 1.4e-23), the finite-difference accuracy;
the major rows of a column have exactly zero derivative in the major coordinates of non-ancestor columns; every 𝕋_k is
square and nonsingular and every A_k is an upper triangular Z-matrix with positive diagonal in the order of `limit.md`
Proposition 4.5 (for N = 3 with pattern TB, component 2 has no trace row and no trace coordinate at all: an empty block,
consistent with Lemma 4.1); the half-profiles are ≥ 0; the pivots are positive; κ_p(η; δ) − κ_p(η; 0) − δ ≤ 6e-61 and
κ_p(0; δ) − e_p − δ ≤ 4e-59 (N = 6: ≤ 2.5e-50). In the finite-stage runs the reconstructed profiles satisfy every unit
balance to ≤ 6e-86 (precision raised to 108–449 digits to absorb the forward amplification), reach x_top to ≤ 1.5e-85
(top gap), have all mole fractions > 0 and S ≥ α_N, and ‖v_n − v_∞‖ grows as n decreases, at the rate allowed by ρ.

Non-convergence at small n (rand:4:4 at n = 50, 35; rand:5:2 at n = 45; near:3:3 at n = 120) is reported as INFO,
not as a failure: Corollary 4.7 asserts a zero near v_∞ only for n ≥ n_0, with n_0 not explicit and depending on
(δ, η, π); in these cases ρ_σ is 0.88–0.97 or σ_min(DG_∞) is small, and the zero is found for larger n.

**One flagged run, resolved as rounding.** rand:5:2 at n = 150 with the default precision (249 digits; zero of G_n found
to residual 5e-241) reported a smallest mole fraction of −1.4e-321: the distillate content of component 5,
x_{1,5} = d_{[1,2],5}/P_1, with a reverse trace net flow d_{[1,2],5} ≈ −4e-322 — a number 80 orders of magnitude below the
accuracy to which the zero v_n was computed, so its sign was not resolved (`runs/v1_rand_5_2_TBTB.log`,
`runs/v1_rand_5_2_TBTB_n150_extra30.log`). Rerun with 230 extra digits (449 digits, residual 2e-444,
`runs/v1_rand_5_2_TBTB_n150_extra230.log`): the same zero (|v_n − v_∞| = 1.98e-8), every unit balance ≤ 1e-318, smallest
mole fraction +5.5e-351, ALL PASS. This is consistent with Theorem 1.4 (b) (S ≥ 1.15 > 0 at every stage forces positive
contents by the M-matrix argument) and is not a counterexample.

---

## 9. The other verifier

My verdicts above were written before opening anything in `verification-A/`. After writing them I appended one dated
line to `verification-A/PROGRESS.md` (as BRIEF-verification-common item 4 requires: "no corrections for theorem.md
§§2–5", with a pointer to §§6–7 here), without reading that folder. I did not read V2's files and make no comparison
with V2's verdicts.

## 10. Credits and reproduction

Verifier: Claude Opus 5.5 (effort xhigh), fresh context. Own scripts (deterministic; Python 3.12 + mpmath 1.3, standard
library otherwise; one thread): `scripts/mutate_combinatorics.py` (< 5 s), `scripts/v1_independent.py CASE [dps] [nlist]
[extra]` (N = 3: seconds; N = 4: 10–80 s; N = 5: 0.5–3 min; N = 6: 2.5–3 min), `scripts/run_batch.sh`,
`scripts/run_batch2.sh` (from `verification-DL/scripts/` in WSL: `bash run_batch.sh`). Each check prints PASS/FAIL and
the script exits non-zero on a FAIL; INFO lines (non-convergence for n below the case's n_0) are not failures.
