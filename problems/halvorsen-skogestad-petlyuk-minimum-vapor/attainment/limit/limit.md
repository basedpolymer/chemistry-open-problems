# The limit system of the generalized Petlyuk arrangement and its nondegeneracy

Agent L (limit system), Claude Opus 5.5 (effort xhigh). Attainment wave, 2026-10-20. Paths are relative to the
problem folder `problems/halvorsen-skogestad-petlyuk-minimum-vapor/`.

**Status of this file.** Complete proofs (no numerical step) of:

- **Theorem 1.4** (exact reformulation): with all stream flows prescribed, the steady states of the generalized
  Petlyuk arrangement with n_σ stages in section σ are exactly the zeros of a square system of N(N−1)² equations in
  the net splits d_c and the junction compositions X_c (N(N−1)² unknowns), at which the reconstructed stage profiles
  have S(x) > 0; at every such zero all mole fractions are positive.
- **Lemma 2.2** (analytic eigen-data, including through zero trace flows of either sign) and **Lemma 2.4**
  (spectral form of the section relation), which turn the exact system, near a configuration with simple spectra
  and nonzero pivots, into G_n(v) = 0 with G_n = G_∞ − (exponentially small, C¹-small terms) (**Lemma 3.7**).
- **Proposition 3.3** (spectral structure of the limit configuration): under four hypotheses (L1)–(L4), all of
  which are among Agent D's (H1)–(H4), every section has N simple positive modes, the pinch is a major mode, every
  trace mode is closed at exactly one end, at a node where its component is absent, with κ_σ > 0, and every pivot is
  nonzero. In other words D's (H5), (H7) and the pivot part of (H8) are *consequences* of (L1)–(L4).
- **Theorem L** (nondegeneracy): DG_∞(v_∞) is invertible; explicitly it is block upper triangular
  [[𝕄, ∗], [0, 𝕋]] (majors / traces), 𝕋 = ⊕_k 𝕋_k is block diagonal over the components and each 𝕋_k is equivalent
  to a triangular nonsingular M-matrix A_k (no chain or cycle argument needed: the closed→open moves are strictly
  monotone in i + m), and 𝕄 is block lower triangular over the columns (parents first) with invertible diagonal
  blocks (Astra's Lemma-1 matrix ⊕ a Cauchy-type matrix for a column with ≥ 3 components, a nonzero scalar for a
  binary column). This holds **for every N ≥ 3 and every forcing pattern** π ∈ {B, T}^{N−1}.
- Sections 5–6: N = 3 and N = 4 written out in full. Section 7: why the unforced system is singular (one junction
  ratio without an equation and one surplus trace equation per unforced binary column), and the role of each
  hypothesis.

What is **not** in this file: the design itself and the proof that it satisfies (L1)–(L4) (Agent D,
`../design/design.md`), the positivity of the half-profiles (D's (H6)) and the passage from v_∞ to finite-stage
steady states with S > 0 (Agent A, `../analysis/analysis.md`; Corollary 4.7 below is the implicit-function step
in the form A needs). Numerical checks (evidence only) are listed in Section 8.

---

## 0. Setting and notation

Components 1, …, N (N ≥ 3), α_1 > α_2 > … > α_N > 0; S(x) = Σ_i α_i x_i; y(x) = (α_i x_i / S(x))_i (equilibrium
vapor of a liquid x); G_x(t) = Σ_i α_i x_i/(α_i − t) (Underwood function); G^{(k)}_x(t) = the same sum without the
term i = k. The feed (flow 1, composition z > 0, quality q) enters the junction stage of column [1, N]. Molar heat of
vaporization 1; constant molar overflow (CMO).

**The arrangement** (BRIEF-common §2; HS 2003 Part 3, Fig. 1; `paper.tex` FTC). Columns c = [i, m],
1 ≤ i < m ≤ N; a column is *non-binary* if m − i ≥ 2, *binary* ([p, p+1]) otherwise. Column c has a top section
T_c, a bottom section B_c and a junction stage J_c (an equilibrium stage without duty) between them. The top end of
T_{[i,m]} is the node J[i, m−1] and the bottom end of B_{[i,m]} is the node J[i+1, m], where J[p, p] denotes the
*product node* of component p: J[1,1] is a total condenser (reflux and liquid distillate, composition x_1),
J[N,N] a partial reboiler (an equilibrium stage with duty R, liquid product x_N), and J[p,p] (2 ≤ p ≤ N−1) a
liquid side-draw equilibrium stage without duty, where B[p−1,p] ends from above and T[p,p+1] ends from below.
Thermal couplings: T_{[i,m+1]} ends at J[i,m] from below (vapor up into J[i,m], liquid X_{[i,m]} back), B_{[i−1,m]}
ends at J[i,m] from above (liquid down into J[i,m], vapor y(X_{[i,m]}) back). Hence a non-product node J[i,m] is
attached to at most four sections: T_{[i,m]} (J is its bottom end), B_{[i,m]} (top end), T_{[i,m+1]} (top end) and
B_{[i−1,m]} (bottom end). Every node is an equilibrium stage except the condenser.

We write h(J[i,m]) := i + m (so h(J[p,p]) = 2p). For a section σ, b(σ) and t(σ) are its bottom and top end nodes.

**A section.** σ has n_σ ≥ 0 equilibrium stages numbered 1, …, n_σ from the bottom, constant vapor V_σ > 0 and
liquid L_σ > 0, and net upward component flow w_σ ∈ ℝ^N (the same at every cross-section, by the component
balances), Σ_i w_{σ,i} = V_σ − L_σ. For a top section w_{T_c} =: d_c (net split to the top), for a bottom section
w_{B_c} = −b_c (b_c = net flow down). Put

  M_σ(w) = V_σ diag(α) − w αᵀ,   Ψ_σ(t; w) = V_σ − Σ_i α_i w_i/(α_i − t),   Ψ^{(k)}_σ = Ψ_σ without the term i = k.

**Operating specification.** All stream flows are prescribed (BRIEF-common, D's (F1)–(F3)): V_σ, L_σ for every
section, the product flows P_1, …, P_N > 0, the feed. They satisfy the total-flow balance and the CMO energy
balance at every unit (junction stages and side draws without duty; the reboiler duty is R). Consequences used
below: V_{T_c} − L_{T_c} = Σ_i d_{c,i} is prescribed (= D_c), P_1 = V_{T[1,2]} − L_{T[1,2]}, and P_p equals the total
net flow into J[p,p] from its two sections. The *feed* of a column is f_{[1,N]} := z and, for c ≠ [1,N],

  f_{[i,m]} := d_{[i,m+1]} + b_{[i−1,m]}   (absent parents contribute 0).

---

## 1. The exact finite-stage system

### 1.1 Unknowns

For every column c the unknowns are the net split d_c ∈ ℝ^N with Σ_i d_{c,i} = V_{T_c} − L_{T_c} and the junction
composition X_c ∈ ℝ^N with Σ_i X_{c,i} = 1. So v = (d_c, X_c)_c ranges over an affine space 𝒱 of dimension

  dim 𝒱 = (N(N−1)/2) · 2(N−1) = N(N−1)².

**Derived quantities** (affine in v, computed parents first, i.e. by decreasing m − i): f_c by the formula above,
b_c := f_c − d_c, and the product compositions

  x_p := (b_{[p−1,p]} + d_{[p,p+1]}) / P_p    (absent terms 0; x_1 = d_{[1,2]}/P_1, x_N = b_{[N−1,N]}/P_N).

Σ_i x_{p,i} = 1 because the total net flow into J[p,p] is P_p (specification). The **end compositions** of a
section are x_bot(σ) := the composition of b(σ) and x_top(σ) := the composition of t(σ) (X_c for a junction, x_p for
a product node). The vapor entering σ from below is y(x_bot(σ)) (every bottom end is an equilibrium stage), and the
liquid entering σ from above has composition x_top(σ) (for T[1,2] this is the reflux, of composition x_1).

### 1.2 Sections

**Lemma 1.1 (cross-sections and the recursion).** Let σ have n stages, flows V, L > 0, and w ∈ ℝ^N with
Σw = V − L. For a sequence x_0, x_1, …, x_{n+1} with S(x_j) ≠ 0 (0 ≤ j ≤ n) the following are equivalent:
(i) the cross-section balances V y(x_j) − L x_{j+1} = w hold for j = 0, …, n;
(ii) x_{j+1} = M(w) x_j / (L S(x_j)) for j = 0, …, n.
Moreover 1ᵀ M(w) = L αᵀ, so in (ii) Σ x_{j+1} = 1 for every x_j with S(x_j) ≠ 0; and if (i) holds then every stage
j ∈ {1, …, n} satisfies its component balance V y(x_{j−1}) + L x_{j+1} = V y(x_j) + L x_j.

*Proof.* (ii) is (i) solved for x_{j+1}: V y(x_j) − w = (V diag(α) x_j − w S(x_j))/S(x_j) = M x_j / S(x_j).
1ᵀ M = V αᵀ − (Σ w) αᵀ = L αᵀ. The stage balance is the difference of the cross-section balances j − 1 and j. ∎

**Lemma 1.2 (projective form).** In the situation of Lemma 1.1 (ii), x_{n+1} = M^{n+1} x_0 / Π_{j=0}^{n} L S(x_j).
Conversely, let x_bot, x_top ∈ ℝ^N with Σ x_top = 1, suppose x_top = c M^{n+1} x_bot for some c ∈ ℝ, and suppose the
recursion (ii) started at x_0 = x_bot has S(x_j) ≠ 0 for j = 0, …, n. Then x_{n+1} = x_top.

*Proof.* The first claim is induction. For the converse, x_{n+1} = c' M^{n+1} x_bot with c' = 1/Π L S(x_j) and
Σ x_{n+1} = 1 (Lemma 1.1), so c' 1ᵀM^{n+1}x_bot = 1; then 1ᵀM^{n+1}x_bot ≠ 0 and Σ x_top = 1 give c = c'. ∎

**Definition 1.3 (exact system).** For stage numbers (n_σ)_σ, the *exact system* is: for every section σ,

  (E_σ)  x_top(σ) and M_σ(w_σ)^{n_σ+1} x_bot(σ) are linearly dependent, and M_σ^{n_σ+1} x_bot(σ) ≠ 0.

Near a configuration where M_σ has simple real spectrum and the pivots below are nonzero, (E_σ) is equivalent to
N − 1 analytic equations (Lemma 2.4), so the system has N(N−1) · (N−1) = N(N−1)² equations, as many as unknowns.
Per column: 2(N−1) unknowns (d_c, X_c) and 2(N−1) equations (T_c, B_c), as in BRIEF-common §3.

**Theorem 1.4 (exact reformulation).** Fix a consistent operating specification and stage numbers n_σ.
(a) Every steady state of the arrangement with these flows defines, through its net flows d_c and junction
compositions X_c, a zero v of the exact system, and its stage compositions are the recursions of Lemma 1.1 started
at x_bot(σ).
(b) Conversely, let v ∈ 𝒱 satisfy (E_σ) for every σ, and suppose that for every σ the recursion of Lemma 1.1 (ii)
started at x_bot(σ) has S(x_j) > 0 for j = 0, …, n_σ, and that S(X_c) > 0, S(x_p) > 0 at every node. Assign to every
stage of σ the composition x_j (1 ≤ j ≤ n_σ), to J_c the composition X_c and to J[p,p] the composition x_p. Then
this is a steady state of the arrangement in the sense of `total/foundations/class.md` §3, with the prescribed flows,
and every mole fraction of every unit is strictly positive.

*Proof.* (a) At a steady state all compositions lie in the simplex, so S ≥ α_N > 0 everywhere; the cross-section
balances of every section hold with w_σ = the net flow (component balances of the stages between the cross-section
and an end), so Lemmas 1.1–1.2 give (E_σ), with x_top(σ) = x_{n+1}. The net flows satisfy the junction balances, i.e.
b_c = f_c − d_c, and the product balances, i.e. x_p as above.

(b) By Lemma 1.2 the recursion reaches x_{n_σ+1} = x_top(σ), so by Lemma 1.1 every cross-section balance of σ
holds with w_σ and every stage of σ satisfies its component balance; all stage compositions have Σ = 1. The
remaining units:
- *Junction J_c, c = [i,m].* Its outflows: vapor y(X_c) into T_c, liquid X_c into B_c, liquid X_c into T_{[i,m+1]}
  (if present), vapor y(X_c) into B_{[i−1,m]} (if present); its inflows: the liquid of the first stage of T_c, the
  vapor of the last stage of B_c, the vapor of the last stage of T_{[i,m+1]}, the liquid of the first stage of
  B_{[i−1,m]}, and the feed if c = [1,N]. The four cross-sections adjacent to J_c have net upward flows d_c, −b_c,
  d_{[i,m+1]}, −b_{[i−1,m]}; hence (outflow − inflow) of each component = d_c + b_c − d_{[i,m+1]} − b_{[i−1,m]}
  (− z for c = [1,N]) = d_c + b_c − f_c = 0 by the definition of b_c.
- *Side draw J[p,p].* (outflow − inflow) = (net flow from J[p,p] into B[p−1,p] and into T[p,p+1]) + P_p x_p
  = −b_{[p−1,p]} − d_{[p,p+1]} + P_p x_p = 0 by the definition of x_p. Same for the reboiler (only B[N−1,N]).
- *Condenser J[1,1]* (a total condenser: liquid outflows of the composition of the vapor inflow). The top
  cross-section of T[1,2] gives V y_n − L x_1 = d_{[1,2]} = P_1 x_1 with V = L + P_1 (specification), so the vapor
  y_n entering the condenser equals x_1: the composition rule holds and so does the balance.
- *Total flows and energy:* prescribed by the specification.
So every unit satisfies its component balances with K-values K_{u,i} = α_i / S(x_u) > 0 (equilibrium units; the
condenser is an equilibrium-type unit without vapor outflow). By `total/network/network.md` Theorem 10.1 (i), for
fixed positive K-values the component-i balances of all units form a linear system whose matrix is a nonsingular
M-matrix with nonnegative inverse, positive in the entry (u, u') when u is reachable from u' along streams with
positive flow; every unit is reachable from the feed unit J[1,N] (all flows positive), so the unique solution, which
is our x_{·,i}, is strictly positive. Normalization holds (all compositions sum to 1). Hence a steady state. ∎

*Remark 1.5.* The exact system needs no a-priori sign information: positivity is a *consequence* of S > 0 along the
profiles (which Agent A proves near the limit configuration from D's (H6)). Conversely, S > 0 cannot be dropped:
zeros of the exact system with S < 0 at some stages exist (Agent C, `../computation/PROGRESS.md`, unforced case).
Theorem 1.4 identifies the steady states with the given flows with the zeros of the exact system whose reconstructed
profiles have S > 0; Corollary 4.7 shows that for n large there is exactly one zero in a fixed ball around v_∞.

---

## 2. Spectral data and their analyticity

Fix a section (V > 0 fixed) and drop σ. Define the polynomial (in t, with coefficients polynomial in w)

  P_w(t) := Π_{i=1}^N (α_i − t) · Ψ(t; w) = V Π_i (α_i − t) − Σ_i α_i w_i Π_{j≠i} (α_j − t).

**Lemma 2.1 (spectrum of M(w)).** (i) det(V t I − M(w)) = (−1)^N V^{N−1} P_w(t): the eigenvalues of M(w) are V t
for the roots t of P_w, with multiplicities. (ii) A vector ℓ is a left eigenvector for V t iff
ℓ_i (α_i − t) = (ℓ·w) α_i / V for all i. If t ∉ {α_i}, the left eigenvectors are the multiples of
ℓ_t := (α_i/(α_i − t))_i and t is a root of Ψ; then ℓ_t · x = G_x(t). If t = α_k, then w_k = 0 and e_k is a left
eigenvector. (iii) If t ∉ {α_i} is a root of Ψ, e_t := (w_i/(α_i − t))_i satisfies M e_t = V t e_t. (iv) (Pinch
identity.) If e is an eigenvector for V t with Σ e = 1 and S(e) ≠ 0, then w_i = V (α_i − t) e_i / S(e) and

  Ψ(s; w) = −(V/S(e)) (s − t) G_e(s)   for all s ∉ {α_i}.

*Proof.* (i) det(λ − V diag α + w αᵀ) = Π(λ − V α_i) · (1 + Σ α_i w_i/(λ − V α_i)) (rank-one update); put λ = V t.
(ii) (ℓᵀM)_i = V α_i ℓ_i − (ℓ·w) α_i. If t ∉ {α_i} and ℓ ≠ 0 then ℓ·w ≠ 0 and ℓ ∝ ℓ_t, and ℓ_t·w = V is Ψ(t) = 0.
If t = α_k, the k-th equation gives (ℓ·w) α_k = 0, so ℓ·w = 0 and ℓ ∝ e_k; then e_k · w = w_k = 0. (iii) Direct:
(M e_t)_i = V α_i e_{t,i} − w_i (α·e_t) and α·e_t = V − Ψ(t) = V. (iv) From M e = V t e: w_i S(e) = V(α_i − t) e_i.
Then Σ_i α_i w_i/(α_i − s) = (V/S(e)) Σ_i α_i e_i (α_i − t)/(α_i − s) = (V/S(e)) (S(e) + (s − t) G_e(s)), using
α_i − t = (α_i − s) + (s − t). ∎

(Part (iv) is attainment-astra eq. (13); it holds also when t = α_h for a zero-flow component in e.)

**Lemma 2.2 (analytic eigen-data through zero trace flows).** Let w⁰ ∈ ℝ^N be such that P_{w⁰} has N simple real
roots t_1⁰ > … > t_N⁰. Then there is a neighborhood W of w⁰ and real-analytic functions t_r : W → ℝ, the roots of
P_w, simple and in the same order. For each r:
(a) *major mode* (t_r⁰ ∉ {α_i}): ℓ_r(w) := ℓ_{t_r(w)} = (α_i/(α_i − t_r(w)))_i is an analytic left eigenvector
and e_r(w) := e_{t_r(w)} an analytic right eigenvector, after shrinking W;
(b) *trace mode* (t_r⁰ = α_k for some k; then w⁰_k = 0): Ψ^{(k)}(α_k; w⁰) ≠ 0, and for w ∈ W

  (α_k − t_r(w)) · Ψ^{(k)}(t_r(w); w) = α_k w_k          (exact),                       (2.1)

so t_r(w) = α_k ⟺ w_k = 0, and t_r(w) − α_k has the sign of −w_k/Ψ^{(k)}(α_k; w⁰); the regularized vectors

  ℓ̃_r(w)_i = (w_k / Ψ^{(k)}(t_r; w)) · α_i/(α_i − t_r)  (i ≠ k),   ℓ̃_r(w)_k = 1,
  ẽ_r(w)_i = w_i/(α_i − t_r)  (i ≠ k),   ẽ_r(w)_k = Ψ^{(k)}(t_r; w)/α_k,

are real-analytic on W (for either sign of w_k), are left and right eigenvectors of M(w) for V t_r(w), and

  ℓ̃_r(w) · x = x_k + (w_k / Ψ^{(k)}(t_r; w)) · G^{(k)}_x(t_r)        (exact).            (2.2)

In particular ℓ̃_r = e_k on {w_k = 0}, whatever the other components of w.

*Proof.* Simple roots of a polynomial whose coefficients are analytic in w are analytic in w (implicit function
theorem, P_w'(t_r) ≠ 0), real roots stay real, and the order persists on a small W. (a) Shrink W so that t_r(w)
stays away from every α_i; Lemma 2.1 (ii)–(iii). (b) Write Q_w(t) := Π_{j≠k}(α_j − t) Ψ^{(k)}(t; w), a polynomial in
(t, w). Then P_w(t) = (α_k − t) Q_w(t) − α_k w_k Π_{j≠k}(α_j − t). At w⁰ (w⁰_k = 0), P_{w⁰} = (α_k − t) Q_{w⁰}, and
simplicity of the root α_k means Q_{w⁰}(α_k) ≠ 0, i.e. Ψ^{(k)}(α_k; w⁰) ≠ 0. Shrink W so that t_r(w) stays away from
α_j (j ≠ k) and Ψ^{(k)}(t_r(w); w) ≠ 0; dividing P_w(t_r) = 0 by Π_{j≠k}(α_j − t_r) gives (2.1). Analyticity of ℓ̃, ẽ
is then clear. If w_k ≠ 0 then t_r ≠ α_k by (2.1) and ℓ̃_r = ((α_k − t_r)/α_k) ℓ_{t_r}, ẽ_r = e_{t_r} (its k-entry
w_k/(α_k − t_r) equals Ψ^{(k)}(t_r)/α_k by (2.1)), eigenvectors by Lemma 2.1. If w_k = 0 then t_r = α_k, ℓ̃_r = e_k,
which is a left eigenvector by Lemma 2.1 (ii); for ẽ_r: α·ẽ_r = Σ_{i≠k} α_i w_i/(α_i − α_k) + Ψ^{(k)}(α_k) = V, so
(M ẽ_r)_i = V α_i ẽ_{r,i} − w_i V = V α_k ẽ_{r,i} for i ≠ k (using w_i = (α_i − α_k) ẽ_{r,i}), and
(M ẽ_r)_k = V α_k ẽ_{r,k} − w_k V = V α_k ẽ_{r,k}. (2.2): ℓ̃_r · x = x_k + (w_k/Ψ^{(k)}(t_r)) Σ_{i≠k} α_i x_i/(α_i − t_r).
The sign statement follows from (2.1) with Ψ^{(k)}(t_r) close to Ψ^{(k)}(α_k; w⁰). ∎

*Remark.* For a major mode one may equally use the normalization ℓ̃ (k-entry 1) relative to the nearest α; for fixed
(δ, η) nothing depends on this choice. For uniformity in the forcing parameters η_p → 0 (the forced root φ'_c tends to
α_p) the regularized normalization is the natural one; Theorem L is stated at fixed (δ, η).

**Lemma 2.3 (one-signed net flows).** Let V > 0 and w ∈ ℝ^N with support the interval of components [a, b].
(i) If w ≥ 0 and L := V − Σw > 0 ("top type"), Ψ(·; w) has exactly one root in each interval (α_{p+1}, α_p),
a ≤ p < b, exactly one root ψ_low ∈ (0, α_b), and no other real root; Ψ is decreasing between consecutive poles;
Ψ(t) > V for t > α_a.
(ii) If w ≤ 0 ("bottom type"), Ψ has exactly one root in each (α_{p+1}, α_p), a ≤ p < b, exactly one root
ψ_high ∈ (α_a, ∞), and no other; Ψ is increasing between consecutive poles; Ψ(t) > V for t < α_b.
All these roots are simple; the remaining N − (b − a + 1) roots of P_w are the trace modes α_k, k ∉ [a, b].

*Proof.* Ψ'(t) = −Σ_i α_i w_i/(α_i − t)², of constant sign. (i) Near a pole α_p from below Ψ → −∞, from above +∞;
on (α_a, ∞) every term −α_i w_i/(α_i − t) is positive; Ψ(−∞) = V, Ψ(0) = L > 0. (ii) Mirror image, with Ψ(+∞) = V.
P_w = Π_{i∉[a,b]}(α_i − t) · [Π_{i∈[a,b]}(α_i − t) Ψ(t)]: the bracket is a polynomial of degree b − a + 1 with the
b − a + 1 simple roots just listed. ∎

**Lemma 2.4 (spectral form of a section relation).** Let M = M(w) have N simple positive eigenvalues
λ_1 > … > λ_N with left eigenvectors ℓ_1, …, ℓ_N (any analytic normalization). Fix r* and x_bot, x_top with
ℓ_{r*}·x_bot ≠ 0 ≠ ℓ_{r*}·x_top ("pivots"). Then (E): x_top = c M^{n+1} x_bot for some c ≠ 0 holds iff for every
r ≠ r*

  r < r* (mode above the pinch, *closed at the bottom*):  ℓ_r·x_bot/ℓ_{r*}·x_bot − (λ_{r*}/λ_r)^{n+1} ℓ_r·x_top/ℓ_{r*}·x_top = 0,
  r > r* (mode below the pinch, *closed at the top*):     ℓ_r·x_top/ℓ_{r*}·x_top − (λ_r/λ_{r*})^{n+1} ℓ_r·x_bot/ℓ_{r*}·x_bot = 0,

and then c = ℓ_{r*}·x_top / (λ_{r*}^{n+1} ℓ_{r*}·x_bot).

*Proof.* (E) gives ℓ_r·x_top = c λ_r^{n+1} ℓ_r·x_bot for all r; divide by the same identity for r*. Conversely, with
c as stated, the displayed equations say ℓ_r·(x_top − c M^{n+1} x_bot) = 0 for every r (for r* by the choice of c);
the ℓ_r form a basis. ∎

So, near a point where Lemma 2.4 applies to every section (with the pinch r*_σ fixed), the exact system is

  G_n(v) = 0,   G_{n,σ,r}(v) := ℓ_r·x_c/ℓ_{r*}·x_c − τ_{σ,r}(w_σ) · ℓ_r·x_o/ℓ_{r*}·x_o,            (2.3)

where (x_c, x_o) = (x_bot, x_top) and τ_{σ,r} = (λ_{r*}/λ_r)^{n_σ+1} for r above r*, (x_c, x_o) = (x_top, x_bot)
and τ_{σ,r} = (λ_r/λ_{r*})^{n_σ+1} for r below r*; ℓ_r is the analytic left eigenvector of Lemma 2.2 (ℓ_t for a
major mode, ℓ̃ for a trace mode). The **limit system** G_∞ is (2.3) with τ ≡ 0:

  G_{∞,σ,r}(v) = ℓ_r(w_σ)·x_c(σ,r) / ℓ_{r*}(w_σ)·x_c(σ,r),                                              (2.4)

i.e. each non-pinch mode is *closed* at one end: modes above the pinch at the bottom end, modes below it at the top
end (`total/sections/sections.md` Lemma 11, with the regularized trace modes).

---

## 3. The limit configuration and the limit system

### 3.1 Section types and roots

Fix a forcing pattern π ∈ {B, T}^{N−1} (π_p = B: component p is sent down B_{[p,p+1]} into product p+1; π_p = T:
component p+1 is sent up T_{[p,p+1]} into product p; Agent D's default is all B). For limit net flows as in (L2)
below, the six section types, their supports [a_σ, b_σ] (an interval of components) and their roots (Lemma 2.3) are:

| type | section | support [a_σ, b_σ] | sign | major modes (roots of Ψ_σ) | pinch r*_σ |
|---|---|---|---|---|---|
| T-nb | T_c, c = [i,m], m−i ≥ 2 | [i, m−1] | w ≥ 0 | r_c ∈ (α_{i+1}, α_i), θ^T_j ∈ (α_{j+1}, α_j) (i < j < m−1), u_c ∈ (0, α_{m−1}) | u_c (= ψ_low) |
| B-nb | B_c, c = [i,m], m−i ≥ 2 | [i+1, m] | w ≤ 0 | v_c ∈ (α_{i+1}, ∞), θ^B_j ∈ (α_{j+1}, α_j) (i < j < m−1), s_c ∈ (α_m, α_{m−1}) | v_c (= ψ_high) |
| T-np | T_c, c = [p,p+1], π_p = B | {p} | w ≥ 0 | φ^T_c = α_p L_T/V_T ∈ (0, α_p) | φ^T_c |
| B-f | B_c, c = [p,p+1], π_p = B | {p, p+1} | w ≤ 0 | φ'_c ∈ (α_p, ∞), φ''_c ∈ (α_{p+1}, α_p) | φ'_c (= ψ_high) |
| T-f | T_c, c = [p,p+1], π_p = T | {p, p+1} | w ≥ 0 | φ^mid_c ∈ (α_{p+1}, α_p), φ^low_c ∈ (0, α_{p+1}) | φ^low_c (= ψ_low) |
| B-np | B_c, c = [p,p+1], π_p = T | {p+1} | w ≤ 0 | φ^B_c = α_{p+1} L_B/V_B ∈ (α_{p+1}, ∞) | φ^B_c |

(np = non-pinning, f = forced/pinning; the notation of the roots is Agent D's.) **Rule:** a top section pinches at
its smallest major root ψ_low, a bottom section at its largest ψ_high. The trace components of σ are k ∉ [a_σ, b_σ];
k is a *lighter* trace if k < a_σ and a *heavier* trace if k > b_σ.

### 3.2 Hypotheses

The limit configuration is built from the specification by Agent D (`../design/hypotheses.md`). Theorem L uses
exactly the following four properties of it.

**(L1) Flows.** V_σ > 0 and L_σ > 0 for every section, P_p > 0 for every product, and the specification is
consistent (total-flow and CMO balances at every unit; the only duty is the reboiler's).

**(L2) Limit net flows.** There are net splits d^∞_c with Σ_i d^∞_{c,i} = V_{T_c} − L_{T_c} such that, with
f^∞_c, b^∞_c := f^∞_c − d^∞_c derived as in §1.1:
- non-binary c = [i,m]: d^∞_{c,k} > 0 for i ≤ k ≤ m−1 and b^∞_{c,k} > 0 for i+1 ≤ k ≤ m; every other entry of
  d^∞_c and b^∞_c is 0;
- binary c = [p,p+1]: if π_p = B, d^∞_c = D_c e_p and b^∞_c = η_p e_p + b_{c,p+1} e_{p+1} with D_c, η_p, b_{c,p+1} > 0;
  if π_p = T, d^∞_c = d_{c,p} e_p + η_p e_{p+1} and b^∞_c = b_{c,p+1} e_{p+1} with all three > 0;
- common middle roots: for every non-binary c with m − i ≥ 3, θ^T_j = θ^B_j =: θ^c_j for i < j < m−1.

**(L3) Pinch separation.** For every section σ: if σ is a top section and b_σ < N, then ψ_low(σ) > α_{b_σ+1}; if σ
is a bottom section and a_σ > 1, then ψ_high(σ) < α_{a_σ−1}. Explicitly: u_c > α_m and v_c < α_i (non-binary);
φ^T_c > α_{p+1} and, for p ≥ 2, φ'_c < α_{p−1} (π_p = B); φ^B_c < α_p and, for p+2 ≤ N, φ^low_c > α_{p+2} (π_p = T).

**(L4) Margins.** In every column the top and bottom sections have different roots in the two outer intervals:
r_c ≠ v_c and u_c ≠ s_c (non-binary); φ^T_c ≠ φ''_c (π_p = B); φ^mid_c ≠ φ^B_c (π_p = T). (Equivalently, since the
common roots of Ψ_{T_c} and Ψ_{B_c} are roots of Ψ_{T_c} − Ψ_{B_c} = U_c − Σ_i α_i f_{c,i}/(α_i − t), i.e. feed roots of
column c: Ψ_{T_c} does not vanish at the extreme feed roots of c.)

*Correspondence with Agent D (hypotheses.md v1, 18:07).* (L1) = (H1). (L2) = (H2), with the middle-root clause of
(H2)/(H3). (L3) is contained in (H3) (α_{m−1} > u_c > … > s_c > α_m, α_i > r_c > … > v_c > α_{i+1}) and in the root
intervals of D's mode table (φ^T_c ∈ (α_{p+1}, α_p), φ'_c ∈ (α_p, α_{p−1}), φ^low_c ∈ (α_{p+2}, α_{p+1}),
φ^B_c ∈ (α_{p+1}, α_p)), proved in D's §4. (L4) is contained in (H3) (r_c > θ^c_i > v_c, u_c > θ^c_{m−1} > s_c) and
(H4) (Ψ_{B_c}(φ^T_c) > 0 ⟺ φ^T_c > φ''_c for π_p = B; Ψ_{T_c}(φ^B_c) > 0 ⟺ φ^B_c < φ^mid_c for π_p = T, by the
monotonicity of Lemma 2.3). Theorem L uses **only the inequalities ≠** in (L4), and nothing of (H6) (half-profiles)
or of the positivity of Astra's cone weights; those are needed by Agent A for positivity, not for nondegeneracy.
D's (H5), (H7) and the pivot part of (H8) are *proved* below from (L1)–(L4) (Proposition 3.3).

**Definition 3.1 (limit configuration v_∞).** d_c := d^∞_c for every c, and
- non-binary c = [i,m]: X_c := the vector of attainment-astra Lemma 1 for the interlacing points
  (r_c, θ^c_{i+1}, …, θ^c_{m−2}, s_c) on the face [i,m]: X_{c,k} ∝ Π_χ(α_k − χ)/(α_k Π_{k'∈[i,m], k'≠k}(α_k − α_{k'}))
  for k ∈ [i,m], 0 outside; the unique normalized vector supported on [i,m] with G_{X_c}(χ) = 0 at these m − i points;
- binary, π_p = B: X_c := normalized e_{φ'_c}(B_c) ∝ (b_{c,p}/(φ'_c − α_p)) e_p + (b_{c,p+1}/(φ'_c − α_{p+1})) e_{p+1};
- binary, π_p = T: X_c := normalized e_{φ^low_c}(T_c) ∝ (d_{c,p}/(α_p − φ^low_c)) e_p + (η_p/(α_{p+1} − φ^low_c)) e_{p+1}.
(The interlacing r_c ∈ (α_{i+1}, α_i), θ^c_j ∈ (α_{j+1}, α_j), s_c ∈ (α_m, α_{m−1}) holds by Lemma 2.3, (L2), (L3).)

**Lemma 3.2 (supports at the limit).** (a) f^∞_c > 0 on [i,m] and 0 outside. (b) X^∞_c > 0 on the range of c and 0
outside. (c) x^∞_p ≥ 0 with support exactly S_p := {p} ∪ {p−1 : π_{p−1} = B} ∪ {p+1 : π_p = T}.

*Proof.* (a) f_{[1,N]} = z > 0. Otherwise f_{[i,m]} = d_{[i,m+1]} + b_{[i−1,m]}; each existing parent is non-binary
(its length is m − i + 1 ≥ 2) and contributes, by (L2), a vector positive exactly on [i,m]; at least one parent
exists. (b) Astra's Lemma 1: numerator and denominator of X_{c,k} have the same number k − i of negative factors.
Binary: both coefficients are positive (φ'_c > α_p > α_{p+1}; φ^low_c < α_{p+1} < α_p). (c) x_p = (b_{[p−1,p]} +
d_{[p,p+1]})/P_p is a sum of two nonnegative vectors whose supports are read off (L2). ∎

### 3.3 Spectral structure at the limit

**Proposition 3.3.** Assume (L1)–(L3). For every section σ, at w^∞_σ:
(a) P_{w^∞_σ} has N simple positive roots: the major modes of the table and the trace modes α_k, k ∉ [a_σ, b_σ].
Every lighter trace mode exceeds every major mode and every heavier trace mode is smaller than every major mode.
(b) Modes above the pinch r*_σ: the lighter trace modes, and (top sections) the major modes other than ψ_low. Modes
below r*_σ: the heavier trace modes, and (bottom sections) the major modes other than ψ_high. Consequently, in
(2.4), a lighter trace mode is closed at the bottom end b(σ), a heavier trace mode at the top end t(σ), the
non-pinch major modes of a top section at its bottom end J_c and those of a bottom section at its top end J_c.
(c) Let α_k be a trace mode of σ and J its closed end. Then x^∞_{J,k} = 0; every component in the support of x^∞_J
is heavier than k if J = b(σ) and lighter than k if J = t(σ); Ψ^{(k)}_σ(α_k; w^∞_σ) > 0. Hence
G^{(k)}_{x^∞_J}(α_k) < 0 at a bottom closed end, > 0 at a top closed end, and

  κ^{(k)}_σ := |G^{(k)}_{x^∞_J}(α_k)| / Ψ^{(k)}_σ(α_k; w^∞_σ) > 0.

(d) If moreover (L4) holds, the pivots ℓ_{r*}·x^∞_{b(σ)} and ℓ_{r*}·x^∞_{t(σ)} are nonzero (ℓ_{r*} = ℓ_{r*_σ}).
(e) ρ̄_∞ := max_σ max{(smaller)/(larger) over pairs of consecutive modes of σ} < 1.

*Proof.* (a) Lemma 2.3 gives the major modes, simple and positive (ψ_low > 0 because L_σ > 0). For a top section all
major modes are < α_{a_σ} ≤ α_k for a lighter trace k; heavier traces satisfy α_k ≤ α_{b_σ+1} < ψ_low by (L3), and
ψ_low is the smallest major mode. For a bottom section all major modes are > α_{b_σ} ≥ α_k for a heavier trace;
lighter traces satisfy α_k ≥ α_{a_σ−1} > ψ_high by (L3). So all N modes are distinct, hence simple roots of P.
(b) Read off (a) and the pinch rule. The closed ends follow from (2.4).
(c) The closed ends and their limit supports (Lemma 3.2):
- top sections, lighter traces (k < a_σ), J = b(σ) = J_c with supp X^∞_c = range(c), whose smallest index is a_σ
  (non-binary: i; binary: p): all heavier than k;
- T-nb, heavier traces (k ≥ m), J = J[i, m−1], support [i, m−1]: all lighter than k;
- T-np, heavier traces (k ≥ p+1), J = J[p,p], support S_p ⊂ {p−1, p} (p+1 ∉ S_p since π_p = B): all lighter;
- T-f, heavier traces (k ≥ p+2), J = J[p,p], S_p ⊂ {p−1, p, p+1}: all lighter;
- B-nb, lighter traces (k ≤ i), J = J[i+1, m], support [i+1, m]: all heavier;
- B-f, lighter traces (k ≤ p−1), J = J[p+1,p+1], S_{p+1} ⊂ {p, p+1, p+2}: all heavier;
- B-np, lighter traces (k ≤ p), J = J[p+1,p+1], S_{p+1} ⊂ {p+1, p+2} (p ∉ S_{p+1} since π_p = T): all heavier;
- bottom sections, heavier traces (k > b_σ), J = t(σ) = J_c, support range(c), largest index b_σ: all lighter.
In each case k ∉ supp x^∞_J, i.e. x^∞_{J,k} = 0, and G^{(k)}_{x^∞_J}(α_k) = Σ_{i∈supp} α_i x^∞_{J,i}/(α_i − α_k) is a
nonempty sum of terms of one sign: negative at a bottom end, positive at a top end. Since w^∞_{σ,k} = 0,
Ψ^{(k)}_σ(·; w^∞) = Ψ_σ(·; w^∞) away from α_k. Top section: for a lighter trace, α_k > α_{a_σ} and Ψ_σ > V_σ > 0 there
(Lemma 2.3 (i)); for a heavier trace, α_k < ψ_low and Ψ_σ decreases on (−∞, α_{b_σ}) from V_σ to −∞ through its only
root ψ_low there, so Ψ_σ(α_k) > 0. Bottom section: heavier trace, α_k < α_{b_σ}, Ψ_σ > V_σ > 0 (Lemma 2.3 (ii));
lighter trace, α_k > ψ_high and Ψ_σ increases on (α_{a_σ}, ∞) through its only root ψ_high there to V_σ, so
Ψ_σ(α_k) > 0.
(d) ℓ_{r*}·x = G_x(r*). Case by case:
- *own junction of a non-binary column* (bottom of T-nb, top of B-nb): G_{X_c}(t) Π_{k∈[i,m]}(α_k − t) is a polynomial
  of degree ≤ m − i with the m − i simple zeros r_c, θ^c_j, s_c (Definition 3.1), so G_{X_c} vanishes nowhere else.
  u_c ∈ (α_m, α_{m−1}) (by (L3) and the table), where the only zero is s_c ≠ u_c; v_c ∈ (α_{i+1}, α_i) (by (L3)),
  where the only zero is r_c ≠ v_c. So G_{X_c}(u_c) ≠ 0 ≠ G_{X_c}(v_c);
- *own junction of a binary column, pinning section* (top of B-f, bottom of T-f): X_c is the eigenvector of the pinch
  mode, and ℓ_{r*}·e_{r*} ≠ 0 for a simple eigenvalue (indeed, letting s → r* in Lemma 2.1 (iv),
  G_e(r*) = −(S(e)/V) Ψ'_σ(r*) ≠ 0);
- *own junction of a binary column, non-pinning section*: for π_p = B, Lemma 2.1 (iv) applied to B_c and its pinch
  X_c = e_{φ'_c} gives G_{X_c}(φ^T_c) = −(S(X_c)/V_{B_c}) Ψ_{B_c}(φ^T_c)/(φ^T_c − φ'_c), and Ψ_{B_c}(φ^T_c) ≠ 0 because
  φ^T_c ∈ (α_{p+1}, α_p) ((L3) and the table), where the only root of Ψ_{B_c} is φ''_c ≠ φ^T_c (L4). For π_p = T the same
  with T_c, its pinch X_c = e_{φ^low_c}, the point φ^B_c ∈ (α_{p+1}, α_p) and the root φ^mid_c ≠ φ^B_c;
- *child junction* (top of T-nb: X_{[i,m−1]} > 0 on [i, m−1], with u_c < α_{m−1}; bottom of B-nb: X_{[i+1,m]} > 0 on
  [i+1, m], with v_c > α_{i+1}): G is a nonempty sum of terms of one sign;
- *product node*: top of T-np, S_p ⊂ {p−1, p}, φ^T_c < α_p; top of T-f, S_p ⊂ {p−1, p, p+1}, φ^low_c < α_{p+1};
  bottom of B-f, S_{p+1} ⊂ {p, p+1, p+2}, φ'_c > α_p; bottom of B-np, S_{p+1} ⊂ {p+1, p+2}, φ^B_c > α_{p+1}. In each
  case all α's of the support lie on one side of r*, so G_{x^∞_p}(r*) ≠ 0 (positive at a top end, negative at a
  bottom end).
(e) Finitely many ratios, each < 1 by (a). ∎

*Remark (signs, for A).* The pivots have the signs stated in D's (H8) when (L4) holds with the strict orientation of
(H3)–(H4) (r_c > v_c, u_c > s_c, φ^T_c > φ''_c, φ^B_c < φ^mid_c); Theorem L does not need the signs.

### 3.4 The limit system near v_∞

**Lemma 3.4 (trace nodes).** For a component k let 𝒥_k be the set of nodes J with x^∞_{J,k} = 0:
𝒥_k = {J[i,m] : i < m, k ∉ [i,m]} ∪ {J[p,p] : k ∉ S_p}. Every section attached to a node of 𝒥_k has w^∞_{σ,k} = 0.
Every J ∈ 𝒥_k satisfies exactly one of: (i) its first index is > k ("k lighter than J"), (ii) its second index is
< k ("k heavier than J").

*Proof.* For J = J[i,m] with m > i: the attached sections T_{[i,m]}, B_{[i,m]}, T_{[i,m+1]}, B_{[i−1,m]} have limit
supports contained in [i,m] (L2) (the parents are non-binary, with d^∞_{[i,m+1]} supported on [i,m] and b^∞_{[i−1,m]}
on [i,m]). For J = J[p,p]: S_p is the union of the supports of b^∞_{[p−1,p]} and d^∞_{[p,p+1]}, which are the supports
of w^∞ of the two attached sections. The dichotomy: k ∉ [i,m] (resp. k ∉ S_p ∋ p). ∎

**Lemma 3.5.** G_∞(v_∞) = 0.

*Proof.* Rows of non-pinch major modes (Proposition 3.3 (b)): for a non-binary column, the rows of T_c are
G_{X_c}(r_c)/G_{X_c}(u_c) and G_{X_c}(θ^c_j)/G_{X_c}(u_c), those of B_c are G_{X_c}(θ^c_j)/G_{X_c}(v_c) and
G_{X_c}(s_c)/G_{X_c}(v_c); all numerators vanish by Definition 3.1 and (L2) (common middle roots). For a binary column
the row of φ''_c (π_p = B) is ℓ_{φ''_c}·e_{φ'_c}/ℓ_{φ'_c}·e_{φ'_c} = 0 (left and right eigenvectors of different
eigenvalues are orthogonal); likewise φ^mid_c for π_p = T. Rows of trace modes: by (2.2), ℓ̃_k·x_J = x_{J,k} +
w_{σ,k}(…) = 0 at v_∞ since x^∞_{J,k} = 0 (Proposition 3.3 (c)) and w^∞_{σ,k} = 0. ∎

**Lemma 3.6 (domain).** Assume (L1)–(L4). There is a compact neighborhood 𝒩 of v_∞ in 𝒱 and constants ρ̄ < 1,
π_0 > 0 such that for v ∈ 𝒩 and every σ: the modes of M_σ(w_σ(v)) are simple, positive, ordered as at v_∞, analytic
in v (Lemma 2.2; w_σ is affine in v), with consecutive ratios ≤ ρ̄; and |ℓ_{r*}·x_{b(σ)}|, |ℓ_{r*}·x_{t(σ)}| ≥ π_0.
On 𝒩 the exact system is equivalent to G_n(v) = 0 (2.3), and G_n (for all stage numbers) and G_∞ are real-analytic.

*Proof.* Proposition 3.3, Lemma 2.2 and continuity; Lemma 2.4 gives the equivalence (the nonvanishing of
M^{n+1}x_bot follows from ℓ_{r*}·M^{n+1}x_bot = λ_{r*}^{n+1} ℓ_{r*}·x_bot ≠ 0). ∎

**Lemma 3.7 (C¹ convergence).** With n := min_σ n_σ, there is C such that on 𝒩
|G_n − G_∞| ≤ C ρ̄^{n+1} and ‖D(G_n − G_∞)‖ ≤ C (n+1) ρ̄^n.

*Proof.* Each entry of G_n − G_∞ is −μ(v)^{n_σ+1} g(v) with μ the ratio of two modes of σ (analytic, 0 < μ ≤ ρ̄) and
g = ℓ_r·x_o/ℓ_{r*}·x_o (analytic on 𝒩, hence bounded with bounded derivative); D(μ^{n_σ+1} g) =
(n_σ+1) μ^{n_σ} g Dμ + μ^{n_σ+1} Dg. ∎

(Sections of different lengths change nothing: only min_σ n_σ enters. No section of the forced design is doubly
pinched: every pivot is nonzero at both ends, Proposition 3.3 (d); this is exactly what the forcing buys, Section 7.)

---

## 4. Theorem L: the derivative of the limit system is invertible

### 4.1 Coordinates, majors and traces

We use the following affine coordinates on the factor (d_c, X_c) of 𝒱 (all other entries are *dependent*,
determined by Σ_k d_{c,k} = V_{T_c} − L_{T_c}, Σ_k X_{c,k} = 1 and b_c = f_c − d_c, given the parents):

| column | trace coordinates (limit value 0) | major coordinates | dependent |
|---|---|---|---|
| non-binary [i,m] | d_{c,k} (k ∉ [i,m]), d_{c,m}, b_{c,i}; X_{c,k} (k ∉ [i,m]) | d_{c,p} (i < p < m−1); X_{c,p} (i < p ≤ m) | d_{c,i} = f_{c,i} − b_{c,i}, d_{c,m−1}, X_{c,i} |
| binary [p,p+1], π_p = B | d_{c,k} (k ≠ p); X_{c,k} (k ∉ {p,p+1}) | X_{c,p+1} | d_{c,p}, X_{c,p} |
| binary [p,p+1], π_p = T | d_{c,k} (k ∉ {p,p+1}), b_{c,p}; X_{c,k} (k ∉ {p,p+1}) | X_{c,p+1} | d_{c,p} = f_{c,p} − b_{c,p}, d_{c,p+1}, X_{c,p} |

Each row of the table lists 2(N−1) coordinates. The change of variables (d_c, X_c)_c ↦ (coordinates)_c is affine,
triangular with respect to the order of the columns (parents first) and invertible, so these are global affine
coordinates on 𝒱. A trace coordinate carries the index of its component (d_{c,k}, b_{c,k}, X_{c,k} → k); all trace
coordinates vanish at v_∞ (L2, Lemma 3.2). The rows of G_∞ are classified as **major rows** (non-pinch major modes)
and **trace rows of component k** (the trace mode α_k of a section σ with w^∞_{σ,k} = 0).

**Lemma 4.1 (counts).** For every column c and every component k, the number of trace rows of component k carried by
the two sections of c equals the number of trace coordinates of component k of c; the number of major rows of c
equals the number of major coordinates of c, namely 2(m − i − 1) for a non-binary column and 1 for a binary column.

*Proof.* k ∉ range(c): both T_c and B_c have the trace mode α_k (2 rows); coordinates d_{c,k}, X_{c,k} (2).
Non-binary c, k = i: B_c has the trace mode α_i, T_c has i as a major (1 row); coordinate b_{c,i} (1). k = m: T_c has
the trace mode α_m (1 row); coordinate d_{c,m} (1). i < k < m: no trace mode, no coordinate. Binary, π_p = B: k = p:
major in both sections, no coordinate; k = p+1: trace mode of T_c only; coordinate d_{c,p+1}. π_p = T: k = p: trace
mode of B_c only, coordinate b_{c,p}; k = p+1: major in both, none. Major rows: a non-binary column has r_c and the
m − i − 2 middle modes in T_c, the m − i − 2 middle modes and s_c in B_c (Proposition 3.3 (b)): 2(m − i − 1); major
coordinates m − i (X) + m − i − 2 (d). A binary column has one non-pinch major mode (φ''_c or φ^mid_c) in its pinning
section and none in the other; one major coordinate. ∎

**Lemma 4.2 (trace quantities are linear in the traces of their component).** Fix k. (a) For every section σ with
w^∞_{σ,k} = 0, the function v ↦ w_{σ,k}(v) is a linear form in the trace coordinates of component k only. (b) For
every node J ∈ 𝒥_k, v ↦ x_{J,k}(v) is a linear form in the trace coordinates of component k only.

*Proof.* (a) By induction over the columns, parents first, we show: if d^∞_{c,k} = 0 then d_{c,k} is such a form,
and if b^∞_{c,k} = 0 then so is b_{c,k}. The entries with d^∞_{c,k} = 0 are, by (L2), exactly the trace coordinates
d_{c,k} of the table. If b^∞_{c,k} = 0, then either b_{c,k} is itself a trace coordinate (b_{c,i} for a non-binary
column, b_{c,p} for π_p = T), or k ∉ range(c), in which case b_{c,k} = d_{[i,m+1],k} + b_{[i−1,m],k} − d_{c,k} with
d^∞_{[i,m+1],k} = b^∞_{[i−1,m],k} = 0 (the parents' limit flows are supported on [i,m]); induction. For σ = T_c,
w_{σ,k} = d_{c,k}; for σ = B_c, w_{σ,k} = −b_{c,k}. (b) J = J_c non-product: x_{J,k} = X_{c,k}, a trace coordinate
(k ∉ range(c)). J = J[p,p]: x_{p,k} = (b_{[p−1,p],k} + d_{[p,p+1],k})/P_p with both limit flows zero (k ∉ S_p);
apply (a). ∎

**Lemma 4.3 (linearized trace rows).** Let α_k be a trace mode of σ, closed at J (so J ∈ 𝒥_k, Proposition 3.3 (c)).
Then the derivative at v_∞ of its row G_{∞,σ,k} = ℓ̃_k·x_J / ℓ_{r*}·x_J is

  DG_{∞,σ,k}(v_∞)[δv] = ( δx_{J,k} + (G^{(k)}_{x^∞_J}(α_k) / Ψ^{(k)}_σ(α_k; w^∞_σ)) · δw_{σ,k} ) / ℓ_{r*}·x^∞_J,

where δx_{J,k} and δw_{σ,k} are the linear forms of Lemma 4.2 applied to δv. With q^k_{J→σ} := the net flow of k from
J into σ (q = w_{σ,k} if J = b(σ), q = −w_{σ,k} if J = t(σ)), the numerator is δx_{J,k} − κ^{(k)}_σ δq^k_{J→σ}.
In particular the trace rows of component k do not depend, at first order at v_∞, on any major coordinate or on any
trace coordinate of another component.

*Proof.* By (2.2), ℓ̃_k·x_J = x_{J,k} + w_{σ,k} g(v) with g = G^{(k)}_{x_J}(t_k)/Ψ^{(k)}_σ(t_k) analytic on 𝒩
(Lemma 2.2, Ψ^{(k)}_σ ≠ 0 near v_∞). At v_∞ both x_{J,k} and w_{σ,k} vanish, so by the product rule the derivative of
the numerator is δx_{J,k} + g(v_∞) δw_{σ,k}, and the derivative of the quotient is that divided by the pivot (the
numerator vanishes at v_∞). At v_∞, t_k = α_k (Lemma 2.2 (b), w^∞_{σ,k} = 0). For the second form: at a bottom closed
end G^{(k)} < 0 and q = w, so g δw = −κ δq; at a top closed end G^{(k)} > 0 and q = −w, so again g δw = −κ δq. The last
statement follows from Lemma 4.2. ∎

**Lemma 4.4 (node identities).** For every component k and every node J ≠ J_{[1,N]}, identically on 𝒱,

  Σ_{σ attached to J} q^k_{J→σ}(v) + P_J x_{J,k}(v) = 0     (P_J := 0 for a non-product node).

*Proof.* For J = J[i,m], m > i: the flows from J into T_{[i,m]}, B_{[i,m]}, T_{[i,m+1]} (J is its top end) and
B_{[i−1,m]} (J is its bottom end) are d_{c,k}, b_{c,k}, −d_{[i,m+1],k}, −b_{[i−1,m],k}, whose sum is d_c + b_c − f_c = 0
(component k). For J = J[p,p]: −b_{[p−1,p],k} − d_{[p,p+1],k} + P_p x_{p,k} = 0 by the definition of x_p. ∎

*Remark (product nodes).* The compositions of the product nodes are not unknowns: x_p is fixed by the balances
(§1.1). Their major contents (p itself, and the designed impurities p − 1 if π_{p−1} = B, p + 1 if π_p = T, which are
O(η) but nonzero in the limit) enter G_∞ only through pivots, which are nonzero (Proposition 3.3 (d)), and through no
numerator: no major mode is closed at a product node (Proposition 3.3 (b): non-pinch major modes are closed at the
column's own junction). Their trace contents enter only the trace rows closed there, linearly (Lemma 4.2 (b)).

### 4.2 The trace block

For each component k let 𝒯_k be its trace coordinates, Σ_k the set of sections σ with w^∞_{σ,k} = 0 (one trace row
each), c(σ) the closed end and o(σ) the other end of the mode α_k in σ, and 𝕋_k := ∂(trace rows of k)/∂(𝒯_k) at v_∞,
a square matrix by Lemma 4.1.

**Proposition 4.5 (trace block).** Assume (L1)–(L4) ((L4) only makes the pivots, hence the rows (2.4), well
defined; the argument below uses only the numerators of the rows and (L1)–(L3)) [precision added by the independent
verifier V1, 2026-10-20]. For every k, 𝕋_k is nonsingular. More precisely, if 𝕋_k δt = 0
then the node contents ξ_J := δx_{J,k} (J ∈ 𝒥_k) solve A_k ξ = 0, where

  (A_k)_{JJ} = P_J + Σ_{σ∈Σ_k: c(σ)=J} 1/κ^{(k)}_σ,   (A_k)_{JJ'} = − Σ_{σ∈Σ_k: c(σ)=J', o(σ)=J} 1/κ^{(k)}_σ  (J ≠ J'),

and A_k is, after a permutation of 𝒥_k, upper triangular with positive diagonal: det A_k = Π_J (A_k)_{JJ} > 0, and
A_k is a nonsingular M-matrix (A_k^{-1} ≥ 0). Its column sums are P_{J'} + Σ{1/κ_σ : c(σ) = J', o(σ) ∉ 𝒥_k} ≥ 0.

*Proof.* *Reduction.* Let 𝕋_k δt = 0 and extend δt by zero to δv (all major coordinates and all trace coordinates of
other components zero). By Lemma 4.3 every trace row of k vanishes on δv:

  ξ_{c(σ)} = κ_σ δq_σ,  δq_σ := δq^k_{c(σ)→σ},  for every σ ∈ Σ_k.                                       (4.1)

The node identity of Lemma 4.4 is linear, so it holds for the variations. At J ∈ 𝒥_k every attached section belongs
to Σ_k (Lemma 3.4); it has J = c(σ), with δq^k_{J→σ} = ξ_J/κ_σ, or J = o(σ), with δq^k_{J→σ} = −δq_σ = −ξ_{c(σ)}/κ_σ
(the same net flow seen from the other end), and c(σ) ∈ 𝒥_k by Proposition 3.3 (c). Hence A_k ξ = 0.

*Triangularity.* Let σ ∈ Σ_k be a section of the column [i,m] (for a binary column m = i + 1, and J[i,i], J[m,m] are
product nodes). If k is a lighter trace of σ, then c(σ) = b(σ), o(σ) = t(σ), every end of σ lying in 𝒥_k has first
index > k (T_{[i,m]}: ends J[i,m], J[i,m−1], and k < i; B_{[i,m]}: ends J[i+1,m], J[i,m], and k ≤ i, the end J[i,m]
being in 𝒥_k only if k < i), and h(o(σ)) = h(c(σ)) − 1 (h(J[i,m]) = i + m; t(T_{[i,m]}) = J[i,m−1],
b(B_{[i,m]}) = J[i+1,m]). If k is a heavier trace of σ, then c(σ) = t(σ), o(σ) = b(σ), every end of σ lying in 𝒥_k has
second index < k (T_{[i,m]}: ends J[i,m−1], J[i,m], and k ≥ m, the end J[i,m] being in 𝒥_k only if k > m;
B_{[i,m]}: ends J[i,m], J[i+1,m], and k > m), and h(o(σ)) = h(c(σ)) + 1. By Lemma 3.4, 𝒥_k = 𝒥^<_k ⊔ 𝒥^>_k (first
index > k, resp. second index < k);
off-diagonal entries (A_k)_{JJ'} ≠ 0 occur only for J = o(σ), J' = c(σ), within 𝒥^<_k with h(J) = h(J') − 1, or within
𝒥^>_k with h(J) = h(J') + 1. Ordering 𝒥^<_k by increasing h, then 𝒥^>_k by decreasing h, every nonzero off-diagonal
entry lies above the diagonal.

*Diagonal.* P_J > 0 at product nodes (L1). A non-product J = J[i,m] ∈ 𝒥^<_k is the bottom end of T_{[i,m]}, whose
support begins at i > k: α_k is a lighter trace mode of T_{[i,m]}, closed at J, and contributes 1/κ > 0. A non-product
J = J[i,m] ∈ 𝒥^>_k is the top end of B_{[i,m]}, whose support ends at m < k: α_k is a heavier trace mode closed at J.
So every diagonal entry is positive (κ > 0 by Proposition 3.3 (c)), A_k is nonsingular, and, being a Z-matrix that
is triangular with positive diagonal, A_k^{-1} = (D(I − D^{-1}E))^{-1} = Σ_j (D^{-1}E)^j D^{-1} ≥ 0 (E ≥ 0 strictly
triangular, so the sum is finite). The column sums: each σ with c(σ) = J' contributes +1/κ_σ to (A_k)_{J'J'} and
−1/κ_σ to row o(σ) if o(σ) ∈ 𝒥_k.

*Conclusion.* ξ = 0, so δq_σ = 0 for every σ ∈ Σ_k by (4.1). Every trace coordinate of k is either X_{c,k} = ξ_{J_c}
or, up to sign, the net flow w_{σ,k} of a section σ ∈ Σ_k (d_{c,k} = w_{T_c,k}, b_{c,k} = −w_{B_c,k}, with
w^∞ = 0), and w_{σ,k} = ±δq_σ. Hence δt = 0. ∎

*Remarks.* (1) No weak-chain or path argument is needed: the closed→open moves of a light trace go strictly up
the grid (h decreases by one) and those of a heavy trace strictly down. Equivalently (BRIEF-common's formulation) A_kᵀ
is weakly chained diagonally dominant: every node is linked, by moves closed→open, to a product node (light trace
k < i: J[i,m] → J[i,m−1] → … → J[i,i+1] → J[i,i]; heavy trace k > m: J[i,m] → J[i+1,m] → … → J[m−1,m] → J[m,m]) or to
a node whose closed section opens into a node where k is major. (2) Since A_k^{-1} ≥ 0 and its diagonal is explicit,
A can read off the signs of the leading exponentially small traces from the inhomogeneous version of (4.1), whose
right-hand side comes from the τ-terms of G_n.

### 4.3 The major block

Order the columns by decreasing length m − i (ties arbitrary); this puts every column after all its ancestors
(the columns [i', m'] ⊋ [i, m]). Let 𝕄 := ∂(major rows)/∂(major coordinates) at v_∞, with rows and columns grouped by
column.

**Proposition 4.6 (major block).** Assume (L1)–(L4). 𝕄 is block lower triangular, and each diagonal block 𝕄_c is
nonsingular. Hence 𝕄 is nonsingular.

*Proof.* *Triangularity.* The major rows of c are closed at J_c (Proposition 3.3 (b)) and are functions of X_c and of
w_{T_c} = d_c, w_{B_c} = d_c − f_c; d_c is a function of the coordinates of c and of f_c, and f_c of the coordinates of
the ancestors of c. So the major rows of c depend only on the coordinates of c and of its ancestors.

*Diagonal block of a non-binary column c = [i,m].* Let δ be a variation of the major coordinates of c only (all other
coordinates fixed; so δf_c = 0, every trace coordinate fixed at 0) with 𝕄_c δ = 0. Then δd_{c,i} = δf_{c,i} −
δb_{c,i} = 0 and δd_{c,m} = 0, δd_{c,k} = 0 (k ∉ [i,m]), and Σ_k δd_{c,k} = 0; so δd_c is supported on {i+1, …, m−1}
with zero sum, and δw_{T_c} = δw_{B_c} = δd_c. For a root χ of Ψ_{T_c} the implicit function theorem gives
Ψ'_{T_c}(χ) δχ = Σ_p α_p δd_{c,p}/(α_p − χ) =: R(χ); for a root of Ψ_{B_c} (= V_B − Σ α_p w_p/(α_p − t), w = d_c − f_c),
Ψ'_{B_c}(χ) δχ = R(χ) as well. At v_∞ every major row of c has the form G_{X_c}(χ)/pivot with G_{X_c}(χ) = 0, so its
variation is [ℓ_χ·δX_c + G'_{X_c}(χ) δχ]/pivot. The rows of θ^c_j (i < j < m−1) in T_c and in B_c therefore give

  ℓ_{θ_j}·δX_c + G'_{X_c}(θ_j) R(θ_j)/Ψ'_{T_c}(θ_j) = 0 = ℓ_{θ_j}·δX_c + G'_{X_c}(θ_j) R(θ_j)/Ψ'_{B_c}(θ_j).

Subtracting: G'_{X_c}(θ_j) R(θ_j) (1/Ψ'_{T_c}(θ_j) − 1/Ψ'_{B_c}(θ_j)) = 0. Here Ψ'_{T_c} < 0 < Ψ'_{B_c} (Lemma 2.3) and
G'_{X_c}(θ_j) ≠ 0 (θ_j is a simple zero of G_{X_c}, Proposition 3.3 (d) proof), so R(θ_j) = 0 for the m − i − 2
middle roots; moreover R(0) = Σ_p δd_{c,p} = 0. Now R(t) Π_{p=i+1}^{m−1}(α_p − t) is a polynomial of degree
≤ m − i − 2 with the m − i − 1 distinct zeros θ_{i+1}, …, θ_{m−2}, 0 (θ_j > α_N > 0), hence identically zero, so every
residue α_p δd_{c,p} Π_{p'≠p}(α_{p'} − α_p) vanishes: δd_c = 0. (This is BRIEF-common's Cauchy-type matrix
(α_p/(α_p − θ_j)) (1/Ψ'_T − 1/Ψ'_B), invertible on {Σ δd = 0}.) Then δχ = 0 for every mode, and the m − i major rows of
r_c, θ^c_j, s_c say ℓ_χ·δX_c = 0 for the m − i interlacing points χ, with δX_c supported on [i,m] and Σ δX_c = 0. By
attainment-astra Lemma 1 (the matrix with rows 1ᵀ and (α_k/(α_k − χ_j))_{k∈[i,m]} is nonsingular), δX_c = 0. So 𝕄_c is
injective, and square (Lemma 4.1).

*Diagonal block of a binary column c = [p,p+1].* Its only major coordinate is X_{c,p+1} (δX_{c,p} = −δX_{c,p+1}). With
everything else fixed, δd_c = 0 and δb_c = δf_c − δd_c = 0, so the modes of T_c, B_c do not move, and the variation of
the row of φ''_c (π_p = B; the row ℓ_{φ''_c}·X_c/ℓ_{φ'_c}·X_c) is

  𝕄_c = (α_{p+1}/(α_{p+1} − φ''_c) − α_p/(α_p − φ''_c)) / ℓ_{φ'_c}·X^∞_c ≠ 0,

because α_{p+1} < φ''_c < α_p makes the bracket negative. For π_p = T the same with φ^mid_c ∈ (α_{p+1}, α_p) and the
pivot ℓ_{φ^low_c}·X^∞_c. ∎

### 4.4 Theorem L

**Theorem L (nondegeneracy of the limit system).** Let N ≥ 3, π ∈ {B, T}^{N−1}, and let a limit configuration
satisfy (L1)–(L4) (in particular, every configuration satisfying Agent D's (H1)–(H4)). Then G_∞ is real-analytic on
a neighborhood of v_∞, G_∞(v_∞) = 0, and DG_∞(v_∞) is invertible. In the coordinates of §4.1, ordering the rows as
(major rows; trace rows of k = 1, …, N) and the coordinates as (major; trace of k = 1, …, N),

  DG_∞(v_∞) = [ 𝕄   ∗  ]        with 𝕋 = diag(𝕋_1, …, 𝕋_N),
              [ 0   𝕋 ]

𝕄 block lower triangular over the columns with nonsingular diagonal blocks (Lemma-1 matrix ⊕ Cauchy-type matrix for a
non-binary column, a nonzero scalar for a binary column), and each 𝕋_k nonsingular, reducing to the triangular
M-matrix A_k of Proposition 4.5. Hence det DG_∞(v_∞) = det 𝕄 · Π_k det 𝕋_k ≠ 0.

*Proof.* Analyticity and G_∞(v_∞) = 0: Lemmas 3.6, 3.5. The zero block: Lemma 4.3. The block diagonal form of 𝕋:
Lemma 4.3 (a trace row of k has zero derivative with respect to trace coordinates of other components). 𝕋_k
nonsingular: Proposition 4.5. 𝕄 nonsingular: Proposition 4.6. A block upper triangular matrix with square
nonsingular diagonal blocks (Lemma 4.1) is nonsingular. ∎

**Corollary 4.7 (finite stages; the implicit-function step).** Under the hypotheses of Theorem L there are ε_0 > 0,
n_0 and C such that, whenever min_σ n_σ ≥ n_0, the exact system has exactly one zero v_n in the closed ball
B̄(v_∞, ε_0), and |v_n − v_∞| ≤ C ρ̄^{n+1}, n = min_σ n_σ.

*Proof.* Let K = ‖DG_∞(v_∞)^{-1}‖ and choose ε_0 such that B̄(v_∞, ε_0) ⊂ 𝒩 and ‖DG_∞(v) − DG_∞(v_∞)‖ ≤ 1/(4K) on it.
By Lemma 3.7, for n large, ‖DG_n(v) − DG_∞(v_∞)‖ ≤ 1/(2K) on the ball and K|G_n(v_∞)| = K|G_n(v_∞) − G_∞(v_∞)| ≤
K C ρ̄^{n+1} ≤ ε_0/2. Then Φ(v) := v − DG_∞(v_∞)^{-1} G_n(v) is a contraction of ratio 1/2 mapping the ball into itself;
its unique fixed point is the unique zero of G_n in the ball, and |v_n − v_∞| ≤ 2|Φ(v_∞) − v_∞| ≤ 2KCρ̄^{n+1}. By
Lemma 3.6 the zeros of G_n in 𝒩 are exactly the zeros of the exact system there. ∎

To conclude that v_n is a steady state with positive compositions one needs S > 0 along the reconstructed stage
profiles (Theorem 1.4 (b)); this is Agent A's step, from D's (H6).

---

## 5. N = 3 written out (pattern π = (B, T), network.md §11.2)

Components A, B, C = 1, 2, 3, volatilities a > b > c. Columns: the prefractionator P = [1,3], U = [1,2] (main column
above the side draw), W = [2,3] (below it). Nodes: J_P (feed stage), J_U, J_W, J[1,1] (condenser), J[2,2] (side
draw), J[3,3] (reboiler). π_1 = B (A is sent down B_U into the side product), π_2 = T (C is sent up T_W into the side
product).

**Unknowns (12)** in the coordinates of §4.1, trace coordinates marked t:
- P: ε_C := d_{P,C} (t, C), ε_A := b_{P,A} (t, A), X_{P,B}, X_{P,C}; dependent d_{P,A} = z_A − ε_A,
  d_{P,B} = D_P − z_A + ε_A − ε_C, X_{P,A} = 1 − X_{P,B} − X_{P,C}; b_P = z − d_P.
- U: d_{U,B} (t, B), d_{U,C} (t, C), X_{U,B}, X_{U,C} (t, C); dependent d_{U,A} = P_1 − d_{U,B} − d_{U,C}, X_{U,A};
  f_U = d_P, b_U = d_P − d_U, so b_{U,A} = η_1 − ε_A + d_{U,B} + d_{U,C} with η_1 := z_A − P_1, b_{U,C} = ε_C − d_{U,C}.
- W: d_{W,A} (t, A), b_{W,B} (t, B), X_{W,A} (t, A), X_{W,C}; dependent d_{W,B} = b_{P,B} − b_{W,B},
  d_{W,C} = D_W − d_{W,A} − d_{W,B}, X_{W,B}; f_W = b_P, b_W = b_P − d_W, so b_{W,A} = ε_A − d_{W,A}.
- Products: x_1 = d_U/P_1, x_2 = (b_U + d_W)/P_2, x_3 = b_W/P_3.

**Limit configuration.** Traces 0; d^∞_P = (z_A, D_P − z_A, 0), b^∞_P = (0, z_B + z_A − D_P, z_C); d^∞_U = P_1 e_A,
b^∞_U = (η_1, d^∞_{P,B}, 0); d^∞_W = (0, b^∞_{P,B}, η_2) (η_2 := D_W − b^∞_{P,B}), b^∞_W = P_3 e_C (P_3 = z_C − η_2);
X^∞_P = Lemma-1 vector of (r_P, s_P); X^∞_U ∝ (η_1/(φ'_U − a), d^∞_{P,B}/(φ'_U − b), 0);
X^∞_W ∝ (0, b^∞_{P,B}/(b − φ^low_W), η_2/(c − φ^low_W)); x^∞_1 = e_A, x^∞_2 = (η_1, z_B, η_2)/P_2, x^∞_3 = e_C.

**Rows (12)** (mode → closed end; [maj] or [tr k]):

| section | type | modes (descending), pinch * | rows |
|---|---|---|---|
| T_P | T-nb, supp {A,B} | r_P, u_P*, c | r_P → J_P [maj]: G_{X_P}(r_P); c → J_U [tr C] |
| B_P | B-nb, supp {B,C} | a, v_P*, s_P | a → J_W [tr A]; s_P → J_P [maj]: G_{X_P}(s_P) |
| T_U | T-np, supp {A} | φ^T_U*, b, c | b → J[1,1] [tr B]; c → J[1,1] [tr C] |
| B_U | B-f, supp {A,B} | φ'_U*, φ''_U, c | φ''_U → J_U [maj]; c → J_U [tr C] |
| T_W | T-f, supp {B,C} | a, φ^mid_W, φ^low_W* | a → J_W [tr A]; φ^mid_W → J_W [maj] |
| B_W | B-np, supp {C} | a, b, φ^B_W* | a → J[3,3] [tr A]; b → J[3,3] [tr B] |

**Major block** (rows r_P, s_P, φ''_U, φ^mid_W; coordinates X_{P,B}, X_{P,C}, X_{U,B}, X_{W,C}). Here no major
coordinate enters any net flow (a 3-component column has no major net-flow coordinate, m − i − 2 = 0), so 𝕄 is block
diagonal: 𝕄 = diag(𝕄_P, 𝕄_U, 𝕄_W) with

  𝕄_P = [ (ℓ_r(B) − ℓ_r(A))/G_{X_P}(u_P)   (ℓ_r(C) − ℓ_r(A))/G_{X_P}(u_P) ]      ℓ_χ(k) := α_k/(α_k − χ),
        [ (ℓ_s(B) − ℓ_s(A))/G_{X_P}(v_P)   (ℓ_s(C) − ℓ_s(A))/G_{X_P}(v_P) ],     r = r_P, s = s_P,

  𝕄_U = (b/(b − φ''_U) − a/(a − φ''_U)) / G_{X_U}(φ'_U),   𝕄_W = (c/(c − φ^mid_W) − b/(b − φ^mid_W)) / G_{X_W}(φ^low_W).

det 𝕄_P = det[[1,1,1],[ℓ_r(A),ℓ_r(B),ℓ_r(C)],[ℓ_s(A),ℓ_s(B),ℓ_s(C)]] / (G_{X_P}(u_P) G_{X_P}(v_P)) ≠ 0 (Lemma 1;
pivots nonzero since u_P ≠ s_P, v_P ≠ r_P), and 𝕄_U ≠ 0, 𝕄_W ≠ 0 because their numerators are negative
(φ''_U ∈ (b, a), φ^mid_W ∈ (c, b)) and their pivots are nonzero (Proposition 3.3 (d)).

**Trace blocks.** Writing κ for κ^{(k)}_σ of Proposition 3.3 (c), all evaluated at the limit:

- *Component A* (𝒥_A = {J_W, J[3,3]}, 3 rows: B_P → J_W, T_W → J_W, B_W → J[3,3]; coordinates ε_A, d_{W,A}, X_{W,A}):
  q's: from J_W into B_P = −ε_A, into T_W = d_{W,A}; from J[3,3] into B_W = −b_{W,A}. Equations
  X_{W,A} = −κ_{B_P} ε_A = κ_{T_W} d_{W,A},  x_{3,A} = −κ_{B_W} b_{W,A},  with the node identities
  ε_A = d_{W,A} + b_{W,A} (J_W) and P_3 x_{3,A} = b_{W,A} (J[3,3]);
  A_A = [[1/κ_{B_P} + 1/κ_{T_W}, −1/κ_{B_W}], [0, P_3 + 1/κ_{B_W}]] (rows/columns J_W, J[3,3]), where
  κ_{B_P} = −G^{(A)}_{X_W}(a)/Ψ^{(A)}_{B_P}(a), κ_{T_W} = −G^{(A)}_{X_W}(a)/Ψ^{(A)}_{T_W}(a),
  κ_{B_W} = (c/(a − c))/(V_{B_W} − c P_3/(a − c)); G^{(A)}_{X_W}(a) = b X_{W,B}/(b − a) + c X_{W,C}/(c − a) < 0,
  Ψ^{(A)}_{B_P}(a) > 0 ⟺ v_P < a, Ψ^{(A)}_{T_W}(a) = V_{T_W} + b d_{W,B}/(a − b) + c η_2/(a − c) > 0,
  V_{B_W} − c P_3/(a − c) > 0 ⟺ φ^B_W < a.
- *Component B* (𝒥_B = {J[1,1], J[3,3]}; rows T_U → J[1,1], B_W → J[3,3]; coordinates d_{U,B}, b_{W,B}):
  x_{1,B} = −κ_{T_U} d_{U,B}, P_1 x_{1,B} = d_{U,B}; x_{3,B} = −κ'_{B_W} b_{W,B}, P_3 x_{3,B} = b_{W,B};
  A_B = diag(P_1 + 1/κ_{T_U}, P_3 + 1/κ'_{B_W}), κ_{T_U} = (a/(a − b))/(V_{T_U} − a P_1/(a − b)) (> 0 ⟺ φ^T_U > b),
  κ'_{B_W} = (c/(b − c))/(V_{B_W} − c P_3/(b − c)) (> 0 ⟺ φ^B_W < b).
- *Component C* (𝒥_C = {J_U, J[1,1]}; rows T_P → J_U, B_U → J_U, T_U → J[1,1]; coordinates ε_C, d_{U,C}, X_{U,C}):
  q's: from J_U into T_P = −ε_C, into B_U = b_{U,C}; from J[1,1] into T_U = −d_{U,C}. Equations
  X_{U,C} = −κ_{T_P} ε_C = κ_{B_U} b_{U,C}, x_{1,C} = −κ'_{T_U} d_{U,C}, with ε_C = d_{U,C} + b_{U,C}, P_1 x_{1,C} = d_{U,C};
  A_C = [[1/κ_{T_P} + 1/κ_{B_U}, −1/κ'_{T_U}], [0, P_1 + 1/κ'_{T_U}]] (rows/columns J_U, J[1,1]),
  κ_{T_P} = G^{(C)}_{X_U}(c)/Ψ^{(C)}_{T_P}(c) (Ψ^{(C)}_{T_P}(c) > 0 ⟺ u_P > c), κ_{B_U} = G^{(C)}_{X_U}(c)/Ψ^{(C)}_{B_U}(c),
  κ'_{T_U} = (a/(a − c))/(V_{T_U} − a P_1/(a − c)); G^{(C)}_{X_U}(c) = a X_{U,A}/(a − c) + b X_{U,B}/(b − c) > 0.

So det DG_∞(v_∞) = ± det 𝕄_P · 𝕄_U 𝕄_W · Π (trace-row pivots)^{-1} · det A_A det A_B det A_C ≠ 0, a product of
explicit nonzero factors. In this case the triangular solution of the trace blocks is immediate: J[1,1] and J[3,3]
first (x_{1,C} = x_{3,A} = x_{1,B} = x_{3,B} = 0), then J_U and J_W.

*The default pattern (B, B)* changes only W: T_W becomes T-np (major B; rows a → J_W [tr A], c → J[2,2] [tr C]),
B_W becomes B-f (majors B, C; rows a → J[3,3] [tr A], φ''_W → J_W [maj]); S_2 = {A, B}, S_3 = {B, C}; component B has
the single trace row T_U → J[1,1]; component C acquires the node J[2,2] (row T_W → J[2,2], closed→open move
J[2,2] → J_W, which is a major node of C).

---

## 6. N = 4 written out (default pattern π = (B, B, B))

Columns [1,4]; [1,3], [2,4]; [1,2], [2,3], [3,4]; nodes J14, J13, J24, J12, J23, J34 and J11 (condenser),
J22, J33 (side draws), J44 (reboiler). Limit supports: X14 {1,2,3,4}, X13 {1,2,3}, X24 {2,3,4}, X12 {1,2},
X23 {2,3}, X34 {3,4}; S_1 = {1}, S_2 = {1,2}, S_3 = {2,3}, S_4 = {3,4}.

**Coordinates (36 = 6 × 6; 11 major, 25 trace).**

| column | trace coordinates (component) | major coordinates |
|---|---|---|
| [1,4] | d_{14,4} (4), b_{14,1} (1) | d_{14,2}; X_{14,2}, X_{14,3}, X_{14,4} |
| [1,3] | d_{13,3} (3), d_{13,4} (4), b_{13,1} (1), X_{13,4} (4) | X_{13,2}, X_{13,3} |
| [2,4] | d_{24,1} (1), d_{24,4} (4), b_{24,2} (2), X_{24,1} (1) | X_{24,3}, X_{24,4} |
| [1,2] | d_{12,2} (2), d_{12,3} (3), d_{12,4} (4), X_{12,3} (3), X_{12,4} (4) | X_{12,2} |
| [2,3] | d_{23,1} (1), d_{23,3} (3), d_{23,4} (4), X_{23,1} (1), X_{23,4} (4) | X_{23,3} |
| [3,4] | d_{34,1} (1), d_{34,2} (2), d_{34,4} (4), X_{34,1} (1), X_{34,2} (2) | X_{34,4} |

Per component: 8 trace coordinates of 1, 4 of 2, 4 of 3, 9 of 4.

**Rows (36).** Modes in descending order, pinch marked *, then the rows "mode → closed end".

| section | type | modes | rows |
|---|---|---|---|
| T14 | T-nb, [1,3] | r, θ_2, u*, α_4 | r → J14, θ_2 → J14 [maj]; α_4 → J13 [tr 4] |
| B14 | B-nb, [2,4] | α_1, v*, θ_2, s | α_1 → J24 [tr 1]; θ_2 → J14, s → J14 [maj] |
| T13 | T-nb, [1,2] | r', u'*, α_3, α_4 | r' → J13 [maj]; α_3, α_4 → J12 [tr 3, 4] |
| B13 | B-nb, [2,3] | α_1, v'*, s', α_4 | α_1 → J23 [tr 1]; s' → J13 [maj]; α_4 → J13 [tr 4] |
| T24 | T-nb, [2,3] | α_1, r'', u''*, α_4 | α_1 → J24 [tr 1]; r'' → J24 [maj]; α_4 → J23 [tr 4] |
| B24 | B-nb, [3,4] | α_1, α_2, v''*, s'' | α_1, α_2 → J34 [tr 1, 2]; s'' → J24 [maj] |
| T12 | T-np, {1} | φ^T*, α_2, α_3, α_4 | α_2, α_3, α_4 → J11 [tr 2, 3, 4] |
| B12 | B-f, {1,2} | φ'*, φ'', α_3, α_4 | φ'' → J12 [maj]; α_3, α_4 → J12 [tr 3, 4] |
| T23 | T-np, {2} | α_1, φ^T*, α_3, α_4 | α_1 → J23 [tr 1]; α_3, α_4 → J22 [tr 3, 4] |
| B23 | B-f, {2,3} | α_1, φ'*, φ'', α_4 | α_1 → J33 [tr 1]; φ'' → J23 [maj]; α_4 → J23 [tr 4] |
| T34 | T-np, {3} | α_1, α_2, φ^T*, α_4 | α_1, α_2 → J34 [tr 1, 2]; α_4 → J33 [tr 4] |
| B34 | B-f, {3,4} | α_1, α_2, φ'*, φ'' | α_1, α_2 → J44 [tr 1, 2]; φ'' → J34 [maj] |

Major rows: 2 + 2 + 1 + 1 + 1 + 1 + 1 + 1 + 1 = 11; trace rows: 8 (comp. 1), 4 (2), 4 (3), 9 (4). The counts match the
coordinates column by column (Lemma 4.1).

**Major block 𝕄 (11 × 11).** Columns ordered [1,4], [1,3], [2,4], [1,2], [2,3], [3,4]. The only major net-flow
coordinate is d_{14,2} (with δd_{14,3} = −δd_{14,2}); the feeds of the children depend on it (f_{13} = d_{14},
f_{24} = b_{14}, and through b_{13} on f_{23}), so 𝕄 is block lower triangular with nonzero subdiagonal entries only in
the column of d_{14,2}. Diagonal blocks:
- [1,4] (rows r, θ^T_2, θ^B_2, s; coordinates d_{14,2}, X_{14,2}, X_{14,3}, X_{14,4}): replacing the row θ^T_2 by
  G_{X}(u)·(row θ^T_2) − G_{X}(v)·(row θ^B_2) gives the row (G'_X(θ_2) R(θ_2) (1/Ψ'_{T14}(θ_2) − 1/Ψ'_{B14}(θ_2)), 0, 0, 0)
  with R(t) = α_2/(α_2 − t) − α_3/(α_3 − t); R(θ_2) > 0 because α_3 < θ_2 < α_2, and Ψ'_{T14} < 0 < Ψ'_{B14}, so
  det 𝕄_{[1,4]} = ± (that entry) × det(Lemma-1 3×3 minor in X for the points r, θ_2, s)/(pivots) ≠ 0;
- [1,3]: 2 × 2 Lemma-1 block for (r', s') on the face {1,2,3}; [2,4]: the same for (r'', s'') on {2,3,4};
- [1,2], [2,3], [3,4]: the scalars (α_{p+1}/(α_{p+1} − φ''_c) − α_p/(α_p − φ''_c))/G_{X_c}(φ'_c) ≠ 0.

**Trace blocks A_k** (nodes, closed→open moves σ: c(σ) → o(σ); moves into a node where k is major give only a
diagonal contribution and are marked "→ maj").
- k = 1 (lighter everywhere; 𝒥_1 = {J23, J24, J33, J34, J44}, h = 5, 6, 6, 7, 8): B14: J24 → maj (J14); T24: J24 → J23;
  B13: J23 → maj (J13); T23: J23 → maj (J22, since 1 ∈ S_2); B24: J34 → J24; T34: J34 → J33; B23: J33 → J23;
  B34: J44 → J34. In the order (J23, J24, J33, J34, J44):

      A_1 = [ 1/κ_B13+1/κ_T23   −1/κ_T24          −1/κ_B23      0                  0            ]
            [ 0                 1/κ_B14+1/κ_T24    0            −1/κ_B24            0            ]
            [ 0                 0                  P_3+1/κ_B23  −1/κ_T34            0            ]
            [ 0                 0                  0            1/κ_B24+1/κ_T34    −1/κ_B34      ]
            [ 0                 0                  0            0                  P_4+1/κ_B34   ]

- k = 2 (𝒥_2 = {J34, J44} ⊔ {J11}): B24: J34 → maj (J24); T34: J34 → maj (J33); B34: J44 → J34; T12 (heavier):
  J11 → maj (J12). A_2 = [[1/κ_B24 + 1/κ_T34, −1/κ_B34, 0], [0, P_4 + 1/κ_B34, 0], [0, 0, P_1 + 1/κ_T12]] in the order
  (J34, J44, J11).
- k = 3 (heavier everywhere; 𝒥_3 = {J22, J12, J11}, h = 4, 3, 2): T13: J12 → maj (J13); B12: J12 → J22; T12: J11 → J12;
  T23: J22 → maj (J23). In the order (J22, J12, J11):
  A_3 = [[P_2 + 1/κ_T23, −1/κ_B12, 0], [0, 1/κ_T13 + 1/κ_B12, −1/κ_T12], [0, 0, P_1 + 1/κ_T12]].
- k = 4 (heavier everywhere; 𝒥_4 = {J33, J23, J13, J22, J12, J11}, h = 6, 5, 4, 4, 3, 2): T14: J13 → maj (J14);
  B13: J13 → J23; T13: J12 → J13; B12: J12 → J22; T24: J23 → maj (J24); B23: J23 → J33; T12: J11 → J12;
  T23: J22 → J23; T34: J33 → maj (J34). In the order (J33, J23, J13, J22, J12, J11), A_4 is upper triangular with
  diagonal (P_3 + 1/κ_T34, 1/κ_T24 + 1/κ_B23, 1/κ_T14 + 1/κ_B13, P_2 + 1/κ_T23, 1/κ_T13 + 1/κ_B12, P_1 + 1/κ_T12) and
  off-diagonal entries (J33, J23) = −1/κ_B23, (J23, J13) = −1/κ_B13, (J23, J22) = −1/κ_T23, (J13, J12) = −1/κ_T13,
  (J22, J12) = −1/κ_B12, (J12, J11) = −1/κ_T12.

(κ_σ means κ^{(k)}_σ for the component of the block.) The numerical check (Section 8) reproduces these node sets and
sizes (5, 3, 3, 6 nodes for k = 1..4) for every random N = 4 instance with pattern BBB.

---

## 7. Remarks

**7.1 Why the forcing is necessary (the unforced limit system is singular).** Suppose a binary column c = [p,p+1] is
not forced: η_p = 0, i.e. d^∞_c = f_p e_p, b^∞_c = f_{p+1} e_{p+1} (exactly balanced sharp split), all other data as
in (L1)–(L4). Then both sections of c are non-pinning (T_c of type T-np, B_c of type B-np, whose pinch φ^B_c lies in
(α_{p+1}, α_p) and whose mode α_p is a *lighter trace* closed at the product node J[p+1,p+1]). No major row is closed
at J_c: every row closed at J_c is a trace row (of T_c, B_c, or of the parents, Proposition 3.3 (b)), whose numerator
x_{J_c,k} + w_{σ,k}(…) does not involve the ratio X_{c,p} : X_{c,p+1} when the traces vanish. Hence

  v_∞ + t (e_{p+1} − e_p)  (in the coordinates X_c)

is a zero of G_∞ for every t for which the two pivots G_{X_c}(φ^T_c), G_{X_c}(φ^B_c) stay nonzero (the rows of all
other sections, and the rows of T_c, B_c closed at the product nodes, do not involve X_c at all): a curve of zeros, so
DG_∞(v_∞) annihilates its tangent e_{p+1} − e_p and is singular. (In the count of Lemma 4.1 the column has one major
coordinate without a major row, and component p one trace row, that of α_p in B_c, more than trace coordinates.)
This is exactly the neutral mode found
numerically in `total/network/network.md` §10.4 (an exchange of the two majors of a binary column, "the position of
the front is a free parameter in the limit") and the free face directions of Astra's rank count
(`total/attainment-astra/attainment.md` §5). With η_p > 0, p becomes a major of B_c: the mode α_p of B_c turns into the
major pinch φ'_c > α_p, its eigenvector pins X_c, and the second major root φ''_c, closed at J_c, supplies the missing
major row; simultaneously component p loses the surplus trace row. (Symmetric for π_p = T.) Numerically
(`scripts/check_unforced.py`, `runs/check_unforced.log`): at η_p = 0 the smallest singular value of DG_∞(v_∞) is 0
to 60 digits, with the null vector on a junction-ratio coordinate; for η_p > 0 it is positive.

**7.2 Relation to Astra's rank count.** Astra's frozen-net-flow system (attainment.md (20), N = 3) has rank 10 out of 12
for two reasons, and the formulation of this file removes both. (i) The C content of the coupling liquid U = J_U is
constrained twice ("U_C = 0" from the prefractionator top and from the upper middle section), and likewise the A content
of J_W; with frozen trace net flows these rows coincide. Here the trace net flows are unknowns: the two rows read
X_{U,C} = −κ_{T_P} ε_C and X_{U,C} = κ_{B_U} (ε_C − d_{U,C}) (Section 5) and are independent; the trace system of each
component is square and triangular (Proposition 4.5). (ii) The major compositions within the AB and BC faces are free;
here they are pinned by the forced sections (7.1). Both ingredients are needed.

**7.3 Dependence on the forcing size (for the order of limits; numerical observation, not used in Theorem L).**
With the plain normalization ℓ_φ = (α_i/(α_i − φ))_i of the forced pinch φ'_c (which tends to α_p as η_p → 0), the
top pivot ℓ_{φ'_c}·X_c is of order 1/η_p, so every row of B_c is of order η_p and σ_min(DG_∞(v_∞)) ∝ η_p (1.40 η for
N = 3, 0.51 η for N = 4 in `runs/check_unforced.log`). With the regularized normalization
ℓ̃_{φ'_c} = ((α_p − φ'_c)/α_p) ℓ_{φ'_c} (Lemma 2.2's normalization relative to α_p), the top pivot tends to X_{c,p} > 0 and
σ_min(DG_∞(v_∞)) stays bounded away from 0 uniformly as η_p → 0: it increases to 0.2207 (N = 3) and 0.0378 (N = 4) as
η_p decreases from 0.1 to 1e-6 (values in [0.13, 0.221] and [0.010, 0.0378]), while it is exactly 0 at η_p = 0: the
limit system is discontinuous at η_p = 0 because the pinch changes from φ'_c to φ''_c. In that normalization the
bottom pivot of B_c is, by (2.2) (which holds for any root near α_p, here (α_p − φ'_c)Ψ^{(p)}_{B_c}(φ'_c) = −α_p b_{c,p})
and x^∞_{p+1,p} = b_{c,p}/P_{p+1},

  ℓ̃_{φ'_c}·x_{p+1} = b_{c,p} ( 1/P_{p+1} + |G^{(p)}_{x_{p+1}}(φ'_c)| / Ψ^{(p)}_{B_c}(φ'_c) ) = O(η_p),

so the exponentially small terms of G_n carry a factor 1/η_p: |G_n − G_∞| ≲ ρ̄^n/η_p. Hence n_0 ≳ log(1/η_p)/log(1/ρ̄),
the same threshold as Astra's (attainment.md after (18)): the limits must be taken in the order n → ∞, then η → 0,
then δ → 0, as in the target theorem. Nothing in Theorem L depends on this remark.

**7.4 What each hypothesis does.** (L1): L_σ > 0 puts ψ_low > 0 (positive spectrum); P_J > 0 makes the product rows of
A_k strictly dominant. (L2): the supports fix the classification of modes, nodes and coordinates, and the signs of
G^{(k)} at the closed ends; common middle roots make v_∞ a zero (Lemma 3.5). (L3): every trace mode lies outside the
major spectrum on the correct side, which gives simple spectra, closed ends at trace nodes and Ψ^{(k)}_σ(α_k) > 0
(κ > 0). If (L3) failed for a top section (ψ_low < α_{b_σ+1}), the trace mode α_{b_σ+1} would be above the pinch and
closed at the bottom end J_c, where b_σ + 1 may be major; the configuration would not be a zero of G_∞ (for a
non-binary column it would contradict X_{c,m} > 0). (L4): nonzero pivots at the own junctions. If (L4) failed in a
binary column (φ^T_c = φ''_c: the column exactly at its minimum vapor, the HS "knot" δ = 0), X_c would be the pinch of
the non-pinning section too: a doubly pinched section, outside the scope of Lemma 2.4 — this is why δ > 0 is needed.

**7.5 Status of the coordinator's Step 2 (BRIEF-common §3).** Confirmed, with the following precisions.
(a) "Trace rows are linear homogeneous in traces": the rows are nonlinear; what holds is that their derivative at v_∞
involves only the trace coordinates of their own component (Lemma 4.3), because by the exact formula (2.2) the
dependence on everything else is multiplied by x_{J,k} or w_{σ,k}, both zero at v_∞. (b) The trace block needs no
chain argument: A_k is triangular in the order of i + m (Proposition 4.5). Its column sums at product nodes are
P_J + (contributions of closed sections opening into major nodes), not always exactly P_J. (c) The binary major row
has derivative 1 only in ratio coordinates with a suitable normalization; what matters is that it is nonzero
(Proposition 4.6). (d) D's (H5), (H7) and the pivot part of (H8) need not be proved separately: they follow from
(L1)–(L4) (Proposition 3.3). (e) No section is doubly pinched in the forced design (all pivots nonzero), so no
logarithmic coordinates (sections.md §8.4, warning) are needed; the formulation (2.3) uses only (d, X) as unknowns and
needs no chart coefficients. (f) No flaw found in Step 2. One caution for Step 1 (sent to D at 18:03): a margin
circulation through the sections of a column with ≥ 4 components breaks its common middle roots unless it follows
Astra's direction; D's design adds the margin δ through R only, which is consistent.

---

## 8. Numerical checks (evidence only; every statement above is proved)

All scripts are in `scripts/`, run under WSL Fedora (`outils/wsl "cd <problem folder>/attainment/limit/scripts &&
python <script>"` from the project root), Python 3.12, mpmath 1.3, numpy, one thread; logs in `runs/`.

- `limitsys.py` (module): the exact system G_n (2.3) and the limit system G_∞ (2.4) in the coordinates (d_c, X_c) of
  §1.1, with the regularized trace eigenvectors of Lemma 2.2, for any N and any pattern; its own test design (HS
  preferred split and Astra's slack in every column with ≥ 3 components relative to its actual feed, forced binary
  columns, R = smallest reboiler duty making every binary column strictly feasible plus δ). It is only a test bed;
  D's design (`../design/`) is the reference.
- `check_limit_jacobian.py 4 6 2` (`runs/check_limit_jacobian.log`, 228 s, **16/16 PASS**): N = 3, 4, 5, 6, four random
  instances each (volatility ratios 1.25–2.2, random z, q ∈ [−1, 2], random patterns among BB, BT, TB, TT, BTB, TTB,
  BBBT, TBBTB, BTTTT, …, Astra slack 1e-4–3e-3, forcing 1e-3–1.6e-2, δ 0.01–0.08). Checks C1 (G_∞(v_∞) ≤ 3e-57; (L1)–(L4);
  middle roots common to 1e-40), C2 (∂(trace rows)/∂(majors) and ∂(trace rows of k)/∂(traces of k' ≠ k) vanish:
  exactly 0 in the finite differences), C3 (every trace row equals the prediction of Lemma 4.3 to ≤ 1e-31 relative),
  C4 (κ > 0, closed ends at trace nodes, A_k nonsingular M-matrices with A_k^{-1} ≥ 0), C5 (major block lower triangular
  over the columns, diagonal blocks nonsingular), C6 (σ_min(DG_∞) from 1.1e-4 to 9.1e-2, σ_max up to 8e4: invertible).
- `check_finite_n.py N3 80,60,40,30,25,20,16` (`runs/check_finite_n_N3.log`, 10 s, PASS): α = (4,2,1), z = 1/3, q = 1,
  pattern (B,T), η_bin = 3e-3, slack 1e-3, δ = 0.05. Newton on the exact system G_n from v_∞ (continuation n = 80 → 25):
  |v_n − v_∞| = 4.2e-7, 3.3e-5, 1.3e-3, 5.3e-3, 1.3e-2; reconstructed stage profiles reach the prescribed ends to
  ≤ 1e-71; **every component balance of every unit (stages, junctions, side draw, condenser, reboiler) ≤ 1e-71**; all
  mole fractions > 0, S ≥ 1; product impurities (distillate, side, bottoms) (3.1e-16, 0.01768, 3.5e-17) at n = 80 and
  (4.4e-6, 0.0177, 1.4e-5) at n = 25; the side-product impurity tends to (η_A + η_C)/S = 0.01768 (the value of
  network.md §11.3, obtained there with a different code and design). At n = 20 the simple continuation stops (not a
  failure of a claim: Corollary 4.7 is for n ≥ n_0).
- `check_finite_n.py N4 80,60,45,35` (`runs/check_finite_n_N4.log`, 47 s, PASS): same for α = (8,4,2,1), z = 1/4, q = 0.5,
  pattern BBB, η_bin = 3e-3: |v_n − v_∞| = 1.3e-8, 2.1e-6, 8.8e-5, 7.3e-4 (n = 80, 60, 45, 35); every unit balance ≤ 3e-43;
  all mole fractions > 0; product impurities (1.3e-18, 0.012, 0.012, 0.0119) at n = 80 → designed (0, 0.012, 0.012, 0.0119).
- `check_finite_n.py N5 80,60,45` (`runs/check_finite_n_N5.log`, 177 s, PASS): α = (10, 6, 3.5, 2, 1), z = 1/5,
  q = 0.8, pattern BTBT: |v_n − v_∞| = 2.4e-6, 8.0e-5, 9.7e-4 (n = 80, 60, 45); every unit balance ≤ 3e-37 (the
  forward reconstruction at 110 digits limits the accuracy at n = 80), all mole fractions > 0; some trace net flows
  are negative (reverse flow of a trace, allowed in 𝒜; positivity of the compositions does not need flow signs).
- `check_unforced.py` (`runs/check_unforced.log`, 21 s, PASS): Remark 7.1 and 7.3 (σ_min vs η_p, plain and regularized
  normalizations; singular at η_p = 0).
- `check_combinatorics.py 12` (`runs/check_combinatorics.log`, 11 s, PASS; standard library, exact integers): the case
  analyses of Proposition 3.3 (c), Lemma 3.4, Lemma 4.1 and Proposition 4.5 redone mechanically from the supports of
  (L2) for **every N = 3..12 and every one of the 2^{N−1} patterns** (4092 configurations, up to 1452 rows): closed ends
  at trace nodes with all components on the correct side, attached sections of trace nodes carry no k, the i/m
  dichotomy, the per-column counts and the totals N(N−1)², the h-monotone closed→open moves, nonzero diagonals and
  acyclicity of every A_k. Mutation tests (swapping the closed-end rule, or the product supports) make it fail.
- `check_D_configs.py` (`runs/check_D_configs.log`, 194 s, **6/6 PASS**): Agent D's exported limit configurations
  (`../design/runs/config_*.json`: n3_N (BT, network.md §11.3 data), n3_B (BB), n4_1 (BBB), n5_1 (BBBB), N = 6 random
  (BBBBB, q = −0.73), N = 7 random (BTBTBT)). Only D's flows and limit net flows are read; L's code recomputes every
  mode, pinch, closed end and junction composition (agreement with D's export ≤ 6e-55) and runs C1–C6: G_∞(v_∞) ≤ 3e-51,
  zero blocks exactly 0, trace rows as predicted (≤ 5e-31), A_k triangular M-matrices, σ_min(DG_∞) = 9.5e-3, 7.8e-3,
  3.5e-3, 2.0e-3, 2.3e-4, 9.0e-5 (N = 7: 252 × 252, 76 major).
- Independent confirmations by other agents (their logs): Agent C (`../computation/`) found DG_∞(v_∞) invertible with
  the block structure of Theorem L for N = 3..6 with an independent solver (analytic Jacobians, N = 6: 306 × 306,
  σ_min 7.4e-5), and Agent A (`../analysis/`) found the finite-stage zeros converging at the rate ρ^{n+1} of
  Corollary 4.7.
