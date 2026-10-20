# Finite-stage attainment of the Petlyuk minimum for every number of components

**Integrated proof. Version 2.2, 2026-10-20 19:06** (v2 18:42; v2.1 18:47: D's self-contained cut proof of κ_p(0; δ) = e_p + δ, flux bookkeeping (4.3), Lemma 4.9's proof, index fix; v2.2: C's final evidence and red-team case in Remark 8.5 and §9) (owner: Agent A; integrates `design/design.md` v1.1 (Agent D),
`limit/limit.md` (Agent L), `analysis/analysis.md` v1.2 (Agent A)). Status: **complete proof, not yet independently
verified** (PROTOCOL: "solved (preliminary)" once integrated in the paper; independent verification pending). Every
statement below is proved here or cited by file and statement from earlier, independently verified or proved work
(listed in §10). Numerical computations are reported in §9 and are never used in a proof.

Credits: design and limit configuration — Claude Opus 5.5 (effort xhigh), Agent D; limit system and its
nondegeneracy — Claude Opus 5.5 (effort xhigh), Agent L; analysis and integration — Claude Opus 5.5 (effort
xhigh), Agent A; numerical red team — Claude Opus 5.5 (effort xhigh), Agent C. Earlier ingredients: the stage-level
lower bound (Agents N, S, F of the second wave), existence for prescribed flows and the M-matrix structure of the
component balances (Agent N), the limit relation of a long section (Agent S), the slack deformation and the
interpolation lemma for a prefractionator (GPT-6-Astra via Codex CLI), HS's preferred-split construction (version 1,
Theorem `achieve`). Paths are relative to the problem folder `problems/halvorsen-skogestad-petlyuk-minimum-vapor/`.

---

## 1. Statement and overview

**Data.** N ≥ 2 components with constant relative volatilities α_1 > α_2 > … > α_N > 0; a feed of unit size with
composition z (z_p > 0, Σ_p z_p = 1) and liquid fraction q ∈ ℝ; feed roots θ_j ∈ (α_{j+1}, α_j) of
Σ_p α_p z_p/(α_p − θ) = 1 − q (j = 1, …, N−1); peaks H_j := Σ_{p>j} α_p z_p/(θ_j − α_p); H* := max_j H_j, the minimum
vapor of the generalized Petlyuk arrangement computed by Halvorsen and Skogestad (HS 2003 Part 3), referred to the
reboiler; e_j := H* − H_j ≥ 0. The physical model is the finite-stage constant-molar-overflow (CMO) equilibrium-stage
model of `total/foundations/class.md` (recalled in §2). "Infinitely many stages and pure products" is the limit
(class.md §4)

  R*(𝒮) := lim_{ε→0+} R*(𝒮; ε),  R*(𝒮; ε) := inf{R(s) : s a steady state of a member of 𝒮, every product of purity ≥ 1 − ε},

where R is the total heat supplied and 𝒮 a family of arrangements; 𝒮_Petlyuk denotes the generalized Petlyuk
arrangement with all possible numbers of stages in its sections.

**Theorem 1 (finite-stage attainment).** For every δ > 0 and ε > 0 there is n_0 such that, for every choice of
section lengths n_σ ≥ n_0 (equal or not), the generalized Petlyuk arrangement with these section lengths has an
operating specification (all stream flows prescribed, the same for all such n_σ) with total heat supply
R = H* + δ, and a steady state at which every mole fraction of every stream is positive and every product has purity
≥ 1 − ε. Hence R*(𝒮_Petlyuk) ≤ H*.

**Corollary 2 (the Halvorsen–Skogestad conjecture in the class 𝒜).** R*(𝒮_Petlyuk) = H*, and every family 𝒮 of
arrangements of the class 𝒜 satisfies R*(𝒮) ≥ R*(𝒮_Petlyuk). In words: no arrangement built from CMO equilibrium
stages, single-phase splitters and mixers, total condensers and vaporizers (any topology, recycles, heat duties
anywhere, one feed split and conditioned at will, products of either phase) needs less total vaporization, in the
infinite-stage pure-product limit, than the generalized Petlyuk arrangement, and the latter needs exactly HS's value
max_j H_j.

*Proof of Corollary 2.* The stage-level bound (class.md Theorem 5.3 and Corollary 5.5, proved for every steady state
of every arrangement of 𝒜) gives R ≥ H_j − c_j ε at purity ≥ 1 − ε, hence R*(𝒮) ≥ H* for every 𝒮 ⊆ 𝒜, in particular
for 𝒮_Petlyuk. Theorem 1 gives R*(𝒮_Petlyuk) ≤ H*. ∎

For N = 2, Theorem 1 is Agent S's Proposition 10 (sections.md §8.1). The proof below is written for N ≥ 3; it applies
verbatim to N = 2 (one binary column, no non-binary column), but we do not need this.

**Overview of the proof.**
(i) *Design* (§4). A forced operating specification F(δ, η, π): every column with three or more components is at the
preferred split of its actual feed, deformed by Astra's slack η_c > 0 (which keeps the internal common Underwood roots
exactly); every binary column [p, p+1] is *forced*: a prescribed net flow η_p > 0 of one key is sent into the
adjacent side product (pattern π_p = B: component p down the bottom section into product p+1; π_p = T: component
p+1 up the top section into product p); the only heat input is the reboiler, R = H* + δ (junctions and side draws adiabatic, the condenser removes heat). Its
infinite-stage limit configuration 𝓛 is explicit (Underwood interpolation vectors, forced pinch eigenvectors), and in
𝓛 every section is *singly pinched*: its pinch mode has nonzero coefficients at both ends.
(ii) *Limit system* (§§3, 5). With the junction compositions and column splits as unknowns, the finite-stage
steady-state equations are G_n(v) = 0, where G_n differs from a limit map G_∞ by terms of size (eigenvalue
ratio)^{n+1}, small in C¹ on a fixed neighborhood of 𝓛. The derivative of G_∞ at 𝓛 is invertible (Theorem 5.9,
"Theorem L"): it is block triangular, with a top-down block for the major variables (Underwood interpolation ⊕ a
Cauchy-type matrix in every column) and one triangular M-matrix block per trace component.
(iii) *Analysis* (§6). A quantitative inverse-function step gives a unique zero v_n of G_n near 𝓛 for all section
lengths ≥ n_0, at distance O(ρ^{n+1}); the stage profiles are reconstructed exactly from v_n and are exponentially
close to the limit half-profiles, so S(x) > 0 at every stage; the M-matrix structure of the component balances then
makes every mole fraction positive.
(iv) *Limits* (§7): n → ∞ at fixed (δ, η), then η → 0 (side-product impurities O(η)), then δ → 0.

Forcing is essential: without it the binary junction compositions are free in the limit, DG_∞ has a kernel (the
neutral mode observed in `total/network/network.md` §10.4), and the sections of the binary columns become doubly
pinched as δ → 0; see Remarks 8.2–8.3.

---

## 2. Model, arrangement, specification, unknowns

**2.1 Model** (class.md §2–3). S(x) := Σ_i α_i x_i = αᵀx; the vapor in equilibrium with a liquid x is
y(x) := diag(α) x/S(x); Δ := {x ∈ ℝ^N : x ≥ 0, 1ᵀx = 1}, 1 the all-ones vector. An equilibrium (E) stage has any
inflows, an external duty Q, and outflows of saturated liquid of composition x and saturated vapor of composition
y(x); a total condenser has vapor inflows and liquid outflows of the composition of the condensed vapor. CMO: at every
unit, vapor out − vapor in = Q (a feed of liquid fraction q brings the vapor 1 − q per unit of feed). A steady state is
an assignment of compositions satisfying every component balance and composition rule; R := Σ_u max(Q_u, 0).

**2.2 The generalized Petlyuk arrangement** (HS 2003 Part 3, Fig. 1; version 1's FTC configuration). Columns
c = [i, m], 1 ≤ i < m ≤ N (the set 𝒞, |𝒞| = N(N−1)/2); c is *non-binary* if m − i ≥ 2 and *binary* ([p, p+1])
otherwise. Nodes J[i, m], 1 ≤ i ≤ m ≤ N. Column [i, m] consists of a top section T_c = T[i, m] with bottom end J[i, m]
and top end J[i, m−1], and a bottom section B_c = B[i, m] with top end J[i, m] and bottom end J[i+1, m].
J_c := J[i, m], i < m, is the *junction* (feed stage) of c, an (E) stage without duty; the process feed enters
J[1, N]. J[p, p] is the *product node* of component p: J[1, 1] is a total condenser (reflux and liquid product 1),
J[N, N] a partial reboiler (an (E) stage delivering liquid product N, the only heat input), and J[p, p],
2 ≤ p ≤ N−1, an (E) side-draw stage without duty (liquid product p) where B[p−1, p] ends from above and T[p, p+1]
from below. A section σ consists of n_σ ≥ 1 (E) stages, numbered 1, …, n_σ upward, with constant liquid and vapor
flows L_σ, V_σ > 0; it receives vapor from its bottom-end node and liquid from its top-end node and returns liquid and
vapor to them (thermal couplings). A junction J[i, m] is attached to T[i, m] (as its bottom end), B[i, m] (top end),
T[i, m+1] (top end, if m < N) and B[i−1, m] (bottom end, if i > 1). bot(σ), top(σ) are the end nodes of σ;
Δ_σ := V_σ − L_σ; h(J[i, m]) := i + m.

**2.3 Operating specification.** An operating specification F prescribes L_σ, V_σ for every section and the product
flows P_1, …, P_N, such that at every unit the total flow balance and the CMO balance hold, with zero duty at every
junction and every side-draw stage. The duties are then fixed by the flows (network.md §10), and R(F) = V_{B[N−1,N]}.

**2.4 Unknowns.** For every column c its *top split* d_c ∈ ℝ^N (the net upward component flow in T_c) with
1ᵀd_c = Δ_{T_c}, and its *junction composition* X_c ∈ ℝ^N with 1ᵀX_c = 1. The set of v = (d_c, X_c)_{c∈𝒞} is an
affine space 𝒱 of dimension 𝒩 := 2(N−1)|𝒞| = N(N−1)². Derived quantities (affine in v, computed parents first,
which is well founded because labels shrink along parent→child):

  f_{[1,N]} := z,  f_{[i,m]} := d_{[i,m+1]}·[m < N] + b_{[i−1,m]}·[i > 1]   (net feed of [i, m] from its parents),
  b_c := f_c − d_c,  w_{T_c} := d_c,  w_{B_c} := −b_c   (net upward component flows of the sections),
  x_p := ( d_{[p,p+1]}·[p < N] + b_{[p−1,p]}·[p > 1] ) / P_p   (product compositions).

The composition of a node is X_c at J_c and x_p at J[p, p]; x_bot(σ), x_top(σ) denote the compositions of the end
nodes of σ. The total balances of F give 1ᵀw_σ = Δ_σ for every σ and 1ᵀx_p = 1.

**2.5 Section matrices.** For w ∈ ℝ^N put M_σ(w) := V_σ diag(α) − w αᵀ, Ψ_σ(t; w) := V_σ − Σ_i α_i w_i/(α_i − t),
G_x(t) := Σ_i α_i x_i/(α_i − t) (Underwood's function); Ψ^{(k)}, G^{(k)}: the same sums without the term i = k (terms
with w_i = 0, resp. x_i = 0, omitted). If 1ᵀw = Δ_σ then 1ᵀM_σ(w) = L_σ αᵀ, and the cross-section balance
V y(x_k) − L x_{k+1} = w is equivalent to L S(x_k) x_{k+1} = M x_k (N's projective form).

**Proposition 2.1 (exact reformulation; limit.md Theorem 1.4, analysis.md Lemma A.1).** Fix F and stage numbers n_σ.
(a) Every steady state defines, through its column splits and junction compositions, a point v ∈ 𝒱 at which, for
every σ, x_top(σ) = c_σ M_σ(w_σ)^{n_σ+1} x_bot(σ) with c_σ > 0; its stage compositions are the iterates
x_{k+1} = M_σ x_k/(1ᵀM_σ x_k) started at x_0 = x_bot(σ).
(b) Conversely, let v ∈ 𝒱 satisfy x_top(σ) = c_σ M_σ(w_σ)^{n_σ+1} x_bot(σ) (some c_σ ∈ ℝ) for every σ; suppose that
the iterates x^σ_k := M_σ^k x_bot(σ)/(1ᵀM_σ^k x_bot(σ)) are defined and have S(x^σ_k) > 0 for 0 ≤ k ≤ n_σ, and that
S > 0 at every node composition. Assigning x^σ_k and y(x^σ_k) to stage k of σ, X_c to J_c and x_p to J[p, p] defines a
steady state of the arrangement with the flows F; every mole fraction of every unit and stream is **strictly
positive**, c_σ > 0, x^σ_{n_σ+1} = x_top(σ), and R = R(F) = V_{B[N−1,N]}.

*Proof.* (a) At a steady state all compositions lie in Δ, so S ≥ α_N > 0. The component balances of the stages between
a cross-section and the bottom end of σ give V_σ y(x_k) − L_σ x_{k+1} = w_σ at every cross-section (0 ≤ k ≤ n_σ,
x_0 := x_bot, x_{n_σ+1} := x_top), i.e. L_σ S(x_k) x_{k+1} = M_σ x_k; iterating, x_top = M_σ^{n_σ+1} x_bot /
Π_k (L_σ S(x_k)). The junction and product balances are b_c = f_c − d_c and the formula for x_p.
(b) Fix σ and drop the index. Since 1ᵀM = Lαᵀ, 1ᵀM^{k+1}x_bot = L S(x_k) 1ᵀM^k x_bot, so x_{k+1} = M x_k/(L S(x_k)), i.e.
the cross-section balances V y(x_k) − L x_{k+1} = w hold for 0 ≤ k ≤ n. The relation and 1ᵀx_top = 1 give
x_top = M^{n+1}x_bot/(1ᵀM^{n+1}x_bot) = x_{n+1} and 1/c = L^{n+1} Π_{k=0}^n S(x_k) > 0. The balance of stage k
(1 ≤ k ≤ n) is the difference of the cross-section balances k − 1 and k. The vapor entering stage 1 is that of the
bottom-end node, of composition y(x_bot) = y(x_0); the liquid entering stage n is that of the top-end node, of
composition x_top = x_{n+1}. At a node J, by the end cross-section balances, the component balance (out − in) is
Σ_{σ: bot(σ)=J} w_σ − Σ_{σ: top(σ)=J} w_σ − (feed) + (product). At J_c, c = [i, m], this is
d_c − b_{[i−1,m]} + b_c − d_{[i,m+1]} (− z if c = [1, N]) = d_c + b_c − f_c = 0. At J[p, p] it is
−d_{[p,p+1]} − b_{[p−1,p]} + P_p x_p = 0. At the condenser, the top cross-section of T[1, 2] gives
V y(x_n) = L x_1 + d_{[1,2]} = (L + P_1) x_1 = V x_1, so the condensed vapor has the composition x_1 of the reflux and
distillate (composition rule) and the balance holds. Total and CMO balances hold by the specification, so the duties
are those of F. *Positivity.* Put K_{u,i} := α_i/S(x_u) > 0 at every unit u with a vapor outflow. For each component
i the balances just verified form a linear system 𝔄_i X_{·,i} = e_{J[1,N]} z_i in the vector X_{·,i} of i-contents of
all units, where 𝔄_i has diagonal entries L_u + V_u K_{u,i} > 0 (total liquid and vapor outflows of u, products
included) and off-diagonal entries −(ℓ_{u'u} + v_{u'u} K_{u',i}) ≤ 0 (liquid and vapor flows from u' to u). The
column of u' has excess (diagonal minus the sum of the moduli of the off-diagonal entries) equal to the flow that u'
sends to product outlets, ≥ 0, and every unit is joined to a product node along positive flows (liquid flows down
every section to a product node). Hence 𝔄_i is a nonsingular M-matrix whose inverse is entrywise positive between
units joined by a path of positive flows (network.md Theorem 10.1 (i); weakly chained diagonal dominance). Every unit
is reached from J[1, N] along positive flows, so x_{u,i} > 0 for every unit and component; y = K x > 0 too. ∎

*Remark.* Positivity needs no sign information on trace net flows or trace contents, only S > 0 at every unit. The
hypothesis S > 0 is essential: G_n below has non-physical zeros far from the limit configuration (§9, Agent C).

---

## 3. Algebra of one section

Fix a section (flows L, V > 0; the index σ is dropped) and w ∈ ℝ^N; put P_w(t) := Π_i (α_i − t) Ψ(t; w)
= V Π_i (α_i − t) − Σ_i α_i w_i Π_{j≠i}(α_j − t), a polynomial in (t, w).

**Lemma 3.1 (spectrum, eigenvectors, pinch identity; design.md Lemmas 2.1, 2.3; limit.md Lemma 2.1).**
(i) det(Vt − M(w)) = (−1)^N V^{N−1} P_w(t): the eigenvalues of M(w) are Vt, t a root of P_w; these are the roots of
Ψ (*major modes*) and the α_k with w_k = 0 (*trace modes*).
(ii) For a root φ of Ψ put e_φ := (w_i/(α_i − φ))_i and ℓ_φ := (α_i/(α_i − φ))_i, so ℓ_φ·x = G_x(φ); for a trace
mode α_k put (e_{α_k})_i := w_i/(α_i − α_k) (i ≠ k), (e_{α_k})_k := Ψ^{(k)}(α_k)/α_k, and ℓ_{α_k} := e_k. Then
M e_φ = Vφ e_φ, ℓ_φᵀM = Vφ ℓ_φᵀ, 1ᵀe_φ = L/φ, and ℓ_φ·e_ψ = 0 for modes φ ≠ ψ; ℓ_φ·e_φ = −Ψ'(φ) for a root,
Ψ^{(k)}(α_k)/α_k for a trace mode. ê_φ := (φ/L) e_φ (1ᵀê_φ = 1) is the *pinch composition* of the mode φ.
(iii) (Pinch identity) G_{e_φ}(t) = −Ψ(t)/(t − φ) for all t ∉ {α_i : w_i ≠ 0} ∪ {φ}. In particular every eigenvector
of a mode φ is annihilated by ℓ_ψ for every other mode ψ, and G_{ê_φ}(t) = −(φ/L)Ψ(t)/(t − φ).
(iv) If the N modes are distinct, the ê_φ form a basis, and x = Σ_φ c_φ(x) ê_φ with c_φ(x) = ℓ_φ·x/ℓ_φ·ê_φ.

*Proof.* (i) Rank-one update: det(λ − V diag α + w αᵀ) = Π_i(λ − Vα_i)(1 + Σ_i α_i w_i/(λ − Vα_i)); put λ = Vt. At a
root t ∉ {α_i} of P_w, Ψ(t) = 0; P_w(α_k) = −α_k w_k Π_{j≠k}(α_j − α_k) vanishes iff w_k = 0. (ii) For a root,
α·e_φ = V − Ψ(φ) = V, so (Me_φ)_i = Vα_i(e_φ)_i − Vw_i = Vφ(e_φ)_i; for a trace mode α·e_φ = (V − Ψ^{(k)}(α_k)) +
Ψ^{(k)}(α_k) = V and the same computation works (the k-th row uses w_k = 0). (ℓᵀM)_j = Vα_jℓ_j − (ℓ·w)α_j with
ℓ_φ·w = V (root) or w_k = 0 (trace). 1ᵀe_φ: Σ_i w_i/(α_i − φ) = φ^{−1}Σ_i w_i(α_i/(α_i − φ) − 1) = φ^{−1}(V − Σw) =
L/φ (similarly for a trace mode). Orthogonality and the values of ℓ_φ·e_φ follow from (iii). (iii) With
1/((α − φ)(α − t)) = (1/(α − t) − 1/(α − φ))/(t − φ): G_{e_φ}(t) = ((V − Ψ(t)) − (V − Ψ(φ)))/(t − φ) for a root; for a
trace mode the sum over i ≠ k gives −(Ψ^{(k)}(t) − Ψ^{(k)}(α_k))/(t − α_k) and the k-term adds Ψ^{(k)}(α_k)/(α_k − t).
ℓ_ψ·e_φ = G_{e_φ}(ψ) = 0 for a root ψ ≠ φ; for a trace mode ψ = α_k ≠ φ, ℓ_ψ·e_φ = (e_φ)_k = w_k/(α_k − φ) = 0.
(iv) Linear algebra (distinct eigenvalues, ℓ_φ·ê_φ ≠ 0 because the modes are simple). ∎

**Lemma 3.2 (one-signed net flows; design.md Lemma 2.2, limit.md Lemma 2.3).** Let w ≠ 0 have all its nonzero entries
of one sign, support S, and let a_1 > … > a_s be the α_i, i ∈ S.
(a) (*Top type*: w ≥ 0, L = V − Σw > 0.) Ψ is strictly decreasing on each interval of ℝ ∖ {a_j}; Ψ > V on (a_1, ∞);
Ψ has exactly one root in each (a_{j+1}, a_j), one root ψ_low ∈ (0, a_s), and no other; Ψ > 0 on (−∞, ψ_low).
(b) (*Bottom type*: w ≤ 0.) Ψ is strictly increasing on each interval of ℝ ∖ {a_j}; Ψ > V on (−∞, a_s); Ψ has
exactly one root in each (a_{j+1}, a_j), one root ψ_high ∈ (a_1, ∞), and no other; Ψ > 0 on (ψ_high, ∞).
All these roots are simple and positive.

*Proof.* Ψ'(t) = −Σ_{i∈S} α_i w_i/(α_i − t)² has constant sign. (a) Near a_j the term −α_j w_j/(a_j − t) tends to −∞
from the left and +∞ from the right; Ψ(−∞) = V, Ψ(0) = L > 0, and on (a_1, ∞) every term is positive; monotonicity
gives the roots. (b) The mirror image, with Ψ(+∞) = V. ∎

**Lemma 3.3 (positive directions; design.md Lemma 2.4).** (a) If w ≤ 0, M ≥ 0 entrywise. (b) If w ≥ 0 and L > 0,
M^{−1} = V^{−1}diag(α)^{−1} + (VL)^{−1}diag(α)^{−1} w 1ᵀ ≥ 0. (c) If all modes are positive and x = Σ_φ c_φ ê_φ with
c_φ ≥ 0 and ê_φ ≥ 0 whenever c_φ ≠ 0, then M^k x ≥ 0 for every k ∈ ℤ.

*Proof.* (a) M = V diag(α) + |w|αᵀ. (b) Sherman–Morrison with 1 − αᵀ(V diag α)^{−1}w = 1 − 1ᵀw/V = L/V and
αᵀdiag(α)^{−1} = 1ᵀ. (c) M^k ê_φ = (Vφ)^k ê_φ. ∎

**Lemma 3.4 (analytic eigen-data through zero trace flows; limit.md Lemma 2.2, analysis.md Lemma A.2, design.md
Lemma 2.5).** Let M(w⁰) have N simple positive eigenvalues Vt⁰_1 > … > Vt⁰_N. There is a neighborhood W of w⁰ on which
the roots t_1(w) > … > t_N(w) > 0 of P_w are simple and real-analytic, with analytic left and right eigenvectors
ℓ_r(w), e_r(w), ℓ_r·e_s = δ_rs; any other analytic choice differs by nonvanishing analytic factors. If t⁰_r = α_k (a
trace mode at w⁰; then w⁰_k = 0), then Ψ^{(k)}(α_k; w⁰) ≠ 0,

  (α_k − t_r(w)) Ψ^{(k)}(t_r(w); w) = α_k w_k   (exactly, w ∈ W),

and the left eigenvector normalized by its k-th entry, ℓ̃_r(w)_i := (w_k/Ψ^{(k)}(t_r; w)) α_i/(α_i − t_r) (i ≠ k),
ℓ̃_r(w)_k := 1, is analytic on W for either sign of w_k and satisfies, exactly,

  ℓ̃_r(w)·x = x_k + (w_k/Ψ^{(k)}(t_r(w); w)) G^{(k)}_x(t_r(w)).                                               (3.1)

In particular ℓ̃_r = e_k on {w_k = 0}, and the derivative of (w, x) ↦ ℓ̃_r(w)·x at a point with w_k = 0, x_k = 0 is
(δw, δx) ↦ δx_k + δw_k G^{(k)}_x(α_k)/Ψ^{(k)}(α_k; w): it involves only the k-th components.

*Proof.* Simple roots of a polynomial with coefficients analytic in w are analytic (implicit function theorem), real
simple roots stay real and keep their order, and the spectral projectors Π_{s≠r}(M − λ_s)/(λ_r − λ_s) are analytic of
rank one, which gives analytic e_r, ℓ_r. For a trace mode, write P_w(t) = (α_k − t)Q_w(t) − α_k w_k Π_{j≠k}(α_j − t),
Q_w(t) := Π_{j≠k}(α_j − t)Ψ^{(k)}(t; w). At w⁰, P_{w⁰} = (α_k − t)Q_{w⁰}, and simplicity of α_k means Q_{w⁰}(α_k) ≠ 0,
i.e. Ψ^{(k)}(α_k; w⁰) ≠ 0. Shrinking W so that t_r(w) stays away from α_j (j ≠ k) and Ψ^{(k)}(t_r; w) ≠ 0, P_w(t_r) = 0
divided by Π_{j≠k}(α_j − t_r) is the displayed identity. If w_k ≠ 0 then t_r ≠ α_k and ℓ̃_r = ((α_k − t_r)/α_k) ℓ_{t_r}
by the identity, a left eigenvector by Lemma 3.1; if w_k = 0, ℓ̃_r = e_k, a left eigenvector by Lemma 3.1 (ii). (3.1)
is the definition written out; the derivative follows because t_r = α_k when w_k = 0. ∎

**Product form.** Let r* be a designated index (the *pinch*), U := {r < r*} (modes above the pinch),
D := {s > r*} (modes below), π_bot := ℓ_{r*}·x_bot, π_top := ℓ_{r*}·x_top (the *pivots*), μ_r := λ_{r*}/λ_r (r ∈ U),
μ_s := λ_s/λ_{r*} (s ∈ D), all < 1 (λ_r := V t_r). Define

  E^{(n)}_r := (ℓ_r·x_bot) π_top − μ_r^{n+1} (ℓ_r·x_top) π_bot   (r ∈ U),
  E^{(n)}_s := (ℓ_s·x_top) π_bot − μ_s^{n+1} (ℓ_s·x_bot) π_top   (s ∈ D),                                    (3.2)

and E^{(∞)} the same with μ^{n+1} replaced by 0.

**Lemma 3.5 (product form of the section relation; analysis.md Lemma A.3, limit.md Lemma 2.4).** On W, with
π_bot π_top ≠ 0:
(a) E^{(n)}(x_bot, x_top, w) = 0 ⟺ x_top = c M(w)^{n+1} x_bot for some c ≠ 0, and then c = π_top/(λ_{r*}^{n+1}π_bot);
(b) E^{(∞)} = 0 ⟺ ℓ_r·x_bot = 0 (r ∈ U) and ℓ_s·x_top = 0 (s ∈ D): every non-pinch mode is **closed** at one end,
modes above the pinch at the bottom end and modes below it at the top end (the limit relation of sections.md Lemma 11,
with the trace modes regularized by (3.1));
(c) on any compact set of (x_bot, x_top) and a compact neighborhood W' ⊂ W of w⁰ with ρ_{W'} := sup_{W'} max_r
λ_{r+1}/λ_r < 1: |E^{(n)} − E^{(∞)}| ≤ C ρ_{W'}^{n+1} and ‖D(E^{(n)} − E^{(∞)})‖ ≤ C(n + 2) ρ_{W'}^{n+1}; at w = w⁰,
|E^{(n)} − E^{(∞)}| ≤ C ρ_0^{n+1}, ρ_0 := max(λ⁰_{r*+1}/λ⁰_{r*}, λ⁰_{r*}/λ⁰_{r*−1});
(d) changing the analytic normalization of the ℓ_r multiplies the components of E^{(n)}, E^{(∞)} by nonvanishing
analytic factors; at a zero of E^{(∞)}, DE^{(∞)} = diag(π_top (r ∈ U), π_bot (s ∈ D))·DΛ with
Λ := ((ℓ_r·x_bot)_{r∈U}, (ℓ_s·x_top)_{s∈D}). So the invertibility of the derivative of the network limit map does not
depend on these normalizations, nor on whether the rows are written as in (3.2) or divided by the pivots
(limit.md (2.3)).

*Proof.* (a) x_top = cM^{n+1}x_bot ⟺ ℓ_r·x_top = cλ_r^{n+1}ℓ_r·x_bot for all r. Given this with c ≠ 0, r = r* fixes
c; multiplying the r-relation by π_bot λ_{r*}^{n+1}/λ_r^{n+1} (r ∈ U), resp. by π_bot (s ∈ D), gives E_r = 0, E_s = 0.
Conversely, with c := π_top/(λ_{r*}^{n+1}π_bot), E = 0 gives back every relation. (b) Immediate. (c) The difference
consists of terms μ^{n+1}×(analytic), and ∇(μ^{n+1}) = (n + 1)μ^{n+1}∇ log μ. (d) Product rule at a zero. ∎

*Reading.* With a_s := ℓ_s·x_bot/π_bot and b_r := ℓ_r·x_top/π_top, (a) says x_bot ∝ e_{r*} + Σ_D a_s e_s +
Σ_U μ_r^{n+1} b_r e_r and x_top ∝ e_{r*} + Σ_U b_r e_r + Σ_D μ_s^{n+1} a_s e_s: S's chart (sections.md Lemma 12) with
t = μ^{n+1}. The coefficients (a, b) are not unknowns; they stay bounded while the pivots stay away from 0 (*singly
pinched* sections). In a *doubly pinched* section a pivot vanishes in the limit and (3.2) degenerates into a product
of two small factors equal to an exponentially small number (sections.md §8.4; attainment-astra §5): the case the
forced design avoids.

**The network maps.** For section lengths (n_σ) and designated pinches r*_σ, define, on a neighborhood of a point
v_∞ ∈ 𝒱 at which Lemma 3.4 applies to every section and all pivots are nonzero,

  G_n(v) := ( E^{(n_σ)}_σ( x_bot(σ)(v), x_top(σ)(v), w_σ(v) ) )_σ ∈ ℝ^𝒩,   G_∞(v) := ( E^{(∞)}_σ(…) )_σ,        (3.3)

N − 1 rows per section, N(N−1) sections: a square system on the 𝒩-dimensional space 𝒱. By Lemma 3.5 (a) and
Proposition 2.1, a zero of G_n near v_∞ at which the reconstructed profiles have S > 0 is a steady state, and every
steady state near v_∞ is a zero of G_n.

---

## 4. The forced design and its limit configuration (Agent D, design.md)

### 4.1 Non-binary columns: preferred split and Astra's slack direction

Fix a non-binary column c = [i, m] with feed f (f_p > 0 for i ≤ p ≤ m, 0 otherwise) and net vapor U ∈ ℝ. The feed
function Φ_c(t) := U − Σ_{p=i}^m α_p f_p/(α_p − t) decreases from +∞ to −∞ on each (α_{p+1}, α_p), so it has exactly one
root θ^c_p there, i ≤ p ≤ m−1 (the *feed roots* of c), simple and analytic in (f, U). Bookkeeping: with
U_{[i,m]} := V_{T[i,m+1]} − V_{B[i−1,m]} (absent terms 0; U_{[1,N]} := 1 − q),

  Φ_{[i,m]} = Ψ_{T[i,m+1]} − Ψ_{B[i−1,m]},   Φ_c = Ψ_{T_c} − Ψ_{B_c}.                                            (4.1)

**Lemma 4.1 (preferred split relative to the actual feed; paper.tex Theorem `achieve`; design.md Lemma 3.1).** Put
ρ_p := Π_{r=i}^{m−1}(α_p − θ^c_r)/Π_{p'∈[i,m−1], p'≠p}(α_p − α_{p'}), V⁰_T := α_i f_i/ρ_i, d⁰_p := V⁰_T ρ_p/α_p
(i ≤ p ≤ m−1), d⁰_p := 0 otherwise, b⁰ := f − d⁰, V⁰_B := V⁰_T − U. Then d⁰_i = f_i, d⁰_p > 0 (i ≤ p ≤ m−1),
b⁰_i = 0, b⁰_p > 0 (i+1 ≤ p ≤ m), V⁰_T, V⁰_B, L⁰_T, L⁰_B > 0, and
Ψ⁰_T(t) = V⁰_T Π_r(t − θ^c_r)/Π_{p=i}^{m−1}(t − α_p), Ψ⁰_B(t) = V⁰_B Π_r(t − θ^c_r)/Π_{p=i+1}^{m}(t − α_p): every feed root
is a common (active) root of both sections.

*Proof.* Partial fractions of Π_r(t − θ^c_r)/Π_{p=i}^{m−1}(t − α_p) = 1 + Σ_p ρ_p/(t − α_p) give the formula for Ψ⁰_T.
Numerator and denominator of ρ_p have p − i negative factors each, so ρ_p > 0, V⁰_T > 0, d⁰ > 0. Ψ⁰_T − Ψ⁰_B = Φ_c, so
Ψ⁰_B vanishes at every θ^c_r; since b⁰_i = 0, Ψ⁰_B Π_{p=i+1}^m(t − α_p) is a polynomial of degree ≤ m − i with leading
coefficient V⁰_B and the m − i roots θ^c_r; V⁰_B = 0 would give b⁰ = 0, contradicting b⁰_m = f_m > 0. The residues give
−α_p b⁰_p = V⁰_B Π_r(α_p − θ^c_r)/Π_{p'≠p}(α_p − α_{p'}), of the sign of −V⁰_B, and b⁰_m > 0 forces V⁰_B > 0.
L⁰_T = Ψ⁰_T(0) and L⁰_B = Ψ⁰_B(0) are positive by the product forms. ∎

**Lemma 4.2 (Astra's direction; attainment-astra (3)–(5); design.md Lemma 3.2).** If m − i = 2 put K ≡ 1, h := 0. If
m − i ≥ 3 put β := Π_{p=i+1}^{m−1}α_p/Π_{r=i+1}^{m−2}θ^c_r, K(t) := (t − β)Π_{r=i+1}^{m−2}(t − θ^c_r)/Π_{p=i+1}^{m−1}(t − α_p),
and h_p := (α_p − β)Π_{r=i+1}^{m−2}(α_p − θ^c_r)/(α_p Π_{p'∈[i+1,m−1], p'≠p}(α_p − α_{p'})) for i+1 ≤ p ≤ m−1,
h_p := 0 otherwise. Then K(t) = 1 − Σ_p α_p h_p/(α_p − t), α_{m−1} < β < α_{i+1},
Σ_p h_p = 0, K(θ^c_r) = 0 for i+1 ≤ r ≤ m−2, and K(θ^c_i), K(θ^c_{m−1}), K(α_i), K(α_m) > 0.

*Proof.* β < α_{i+1} and β > α_{m−1} by comparing the factors θ^c_r with α_{r+1} and α_r. K has equal degrees,
leading coefficient 1 and simple poles at α_{i+1..m−1} with residues α_p h_p; K(0) = 1 gives Σh = 0. At θ^c_i and α_i
all factors are positive; at θ^c_{m−1} and α_m numerator and denominator each have m − i − 1 negative factors. ∎

**Definition (deformed column).** For η ≥ 0: d := d⁰ + ηh, b := f − d, V_T := V⁰_T + η, V_B := V⁰_B + η, L_T := V_T − Σd,
L_B := V_B + Σb. Then d_i = f_i, d_m = 0, b_i = 0, b_m = f_m, Σd = Σd⁰, and (by Lemma 4.2) Ψ_{T_c} = Ψ⁰_T + ηK,
Ψ_{B_c} = Ψ⁰_B + ηK. The column is **valid** if (V): d_p > 0 (i ≤ p ≤ m−1) and b_p > 0 (i+1 ≤ p ≤ m); (V) holds for all
η ≥ 0 if m − i = 2, and in general for η < η_0(c) := min({d⁰_p/|h_p| : h_p < 0} ∪ {b⁰_p/h_p : h_p > 0}).

**Proposition 4.3 (a valid deformed column; design.md Proposition 3.3; Astra's Lemma 2).** Let η > 0 satisfy (V).
(a) V_T, L_T, V_B, L_B > 0.
(b) Ψ_{T_c} has exactly the simple roots r_c, θ^c_{i+1}, …, θ^c_{m−2}, u_c and Ψ_{B_c} exactly v_c, θ^c_{i+1}, …,
θ^c_{m−2}, s_c, with α_i > r_c > θ^c_i > v_c > α_{i+1} and α_{m−1} > u_c > θ^c_{m−1} > s_c > α_m; the middle roots are
common; Ψ_{T_c}(α_m) > 0 and Ψ_{B_c}(α_i) > 0.
(c) T_c has the trace modes α_k, k ∉ [i, m−1], and B_c the trace modes α_k, k ∉ [i+1, m]; the trace modes of T_c are
> r_c (k < i) or < u_c (k ≥ m), those of B_c are > v_c (k ≤ i) or < s_c (k > m); all modes positive and distinct.
(d) The pinch compositions ê^T_{(m)} := ê^{T_c}_{α_m}, ê^T_u := ê^{T_c}_{u_c}, ê^B_v := ê^{B_c}_{v_c},
ê^B_{(i)} := ê^{B_c}_{α_i} satisfy ê^T_{(m)}, ê^B_{(i)} > 0 on [i, m]; ê^T_u > 0 on [i, m−1], ê^B_v > 0 on [i+1, m]; all
vanish outside their stated supports.
(e) (Interpolation; attainment-astra Lemma 1.) The vector X_c with X_{c,p} ∝ Π_χ(α_p − χ)/(α_p Π_{p'∈[i,m], p'≠p}
(α_p − α_{p'})) (p ∈ [i, m]), 0 outside, χ ranging over the interlacing points r_c, θ^c_{i+1}, …, θ^c_{m−2}, s_c, is the
unique normalized vector supported on [i, m] with G_{X_c}(χ) = 0 at these m − i points; X_c > 0 on [i, m]; the matrix
with rows 1ᵀ and (α_p/(α_p − χ))_{p∈[i,m]} is nonsingular. Moreover

  X_c = A_c ê^T_{(m)} + B_c ê^T_u = C_c ê^B_v + D_c ê^B_{(i)},   A_c + B_c = C_c + D_c = 1,   A_c, B_c, C_c, D_c > 0.

*Proof.* (a) Lemma 4.1 and Σh = 0. (b) T_c has w = d ≥ 0, support [i, m−1], L_T > 0: by Lemma 3.2 (a) one root in each
(α_{p+1}, α_p), i ≤ p ≤ m−2, and one in (0, α_{m−1}). Ψ_{T_c}(θ^c_r) = ηK(θ^c_r) = 0 for the middle r; Ψ_{T_c}(θ^c_i) =
ηK(θ^c_i) > 0 with Ψ_{T_c} decreasing on (α_{i+1}, α_i) gives r_c > θ^c_i; Ψ_{T_c}(θ^c_{m−1}) > 0 with Ψ_{T_c} decreasing on
(−∞, α_{m−1}) gives u_c > θ^c_{m−1} > α_m and Ψ_{T_c}(α_m) > 0. B_c is symmetric (Lemma 3.2 (b), Ψ_{B_c}(θ^c_i) > 0 and
increasing). (c) Lemma 3.1 and (b). (d) Lemma 3.1 (ii): (ê^T_u)_p = (u_c/L_T)d_p/(α_p − u_c) > 0 for p ≤ m−1;
(ê^T_{(m)})_p = (α_m/L_T)d_p/(α_p − α_m) > 0 (p < m), (ê^T_{(m)})_m = Ψ_{T_c}(α_m)/L_T > 0; symmetrically for B_c.
(e) G_{X}(t)Π_{p∈[i,m]}(α_p − t) is a polynomial of degree ≤ m − i; prescribing its m − i zeros fixes it up to a scalar,
residues give the formula, numerator and denominator of X_{c,p} have p − i negative factors each; the homogeneous
version gives the nonsingularity. In T_c's eigenbasis the coefficient of X_c on a major mode χ ∈ {r_c, θ^c_{middle}} is
∝ G_{X_c}(χ) = 0 and on a trace mode α_k (k ∉ [i, m]) it is ∝ X_{c,k} = 0, so X_c = Aê^T_{(m)} + Bê^T_u, A + B = 1. By the
pinch identity, G_{ê^T_{(m)}}(s_c) = −(α_m/L_T)Ψ_{T_c}(s_c)/(s_c − α_m) < 0 (Ψ_{T_c}(s_c) > 0 since s_c < θ^c_{m−1} < u_c) and
G_{ê^T_u}(s_c) = −(u_c/L_T)Ψ_{T_c}(s_c)/(s_c − u_c) > 0; G_{X_c}(s_c) = 0 then forces A, B > 0. The bottom decomposition is
the same at r_c. ∎

(At η = 0, B_c = C_c = 0: X_c is the pinch of both sections — a doubly pinched feed stage — which is why every
non-binary column needs η_c > 0.)

### 4.2 Binary columns: margin and forcing

Fix a binary column c = [p, p+1] with feed f = f_p e_p + f_{p+1} e_{p+1} (f_p, f_{p+1} > 0), net vapor U, vapors V_T,
V_B = V_T − U; let θ' ∈ (α_{p+1}, α_p) be its feed root and

  Ψ⁰_T(t) := V_T − α_p f_p/(α_p − t),  Ψ⁰_B(t) := V_B + α_{p+1}f_{p+1}/(α_{p+1} − t),  κ := Ψ⁰_T(θ') = Ψ⁰_B(θ'),       (4.2)

the sharp section functions and the **binary margin** (the vapor above the column's Underwood minimum).

**Proposition 4.4 (sharp binary column; design.md Proposition 4.1).** If κ > 0, then V_T > f_p, V_B > 0, the roots
φ⁰_T := α_p(V_T − f_p)/V_T and φ⁰_B := α_{p+1}(V_B + f_{p+1})/V_B satisfy α_{p+1} < φ⁰_B < θ' < φ⁰_T < α_p, and
Ψ⁰_B > κ on (θ', ∞), Ψ⁰_T > κ on (−∞, θ'). (κ > 0 ⟺ φ⁰_T > φ⁰_B ⟺ x̂_T < x̂_B of sections.md Prop. 10.)

*Proof.* V_T = κ + α_pf_p/(α_p − θ'), V_B = κ + α_{p+1}f_{p+1}/(θ' − α_{p+1}); monotonicity of Ψ⁰_T, Ψ⁰_B. ∎

**B-pinning** (π_p = B), η > 0: d = (f_p − η)e_p, b = ηe_p + f_{p+1}e_{p+1}, L_T = V_T − f_p + η, L_B = V_B + f_{p+1} + η.

**Proposition 4.5 (design.md Proposition 4.2).** Let κ > 0 and 0 < η < η̂_p, η̂_p := f_p (p = 1),
η̂_p := min(f_p, (α_{p−1} − α_p)Ψ⁰_B(α_{p−1})/α_p) (p ≥ 2). Then:
(a) all phase flows are positive; (b) T_c has the single major mode φ^T = α_pL_T/V_T ∈ [φ⁰_T, α_p) and the trace
modes α_k, k ≠ p; (c) B_c has exactly the major modes φ'' ∈ (α_{p+1}, φ⁰_B) and φ' ∈ (α_p, α_{p−1}) (φ' > α_1 if p = 1),
and the trace modes α_k (k < p) > φ' > φ'' > α_k (k > p+1); (d) margin: Ψ_{B_c}(φ^T) > κ > 0, equivalently φ'' < φ^T;
(e) the forced pinch X := ê^{B_c}_{φ'} ∝ (η/(φ' − α_p))e_p + (f_{p+1}/(φ' − α_{p+1}))e_{p+1} is positive on {p, p+1} and
X = A'ê^T_{(p+1)} + B'e_p, A', B' > 0, A' + B' = 1, where e_p = ê^{T_c}_{φ^T} and ê^T_{(p+1)} := ê^{T_c}_{α_{p+1}} > 0 on
{p, p+1}; (f) as η → 0, X tends to the pinch of the unforced bottom section at its trace mode α_p.

*Proof.* (a) Proposition 4.4. (b) The root of V_T − α_p(f_p − η)/(α_p − t). (c) Lemma 3.2 (b); the added term
α_pη/(α_p − t) > 0 on (α_{p+1}, α_p) gives φ'' < φ⁰_B; φ' < α_{p−1} ⟺ Ψ_{B_c}(α_{p−1}) > 0 ⟺ the bound on η, and
Ψ⁰_B(α_{p−1}) > κ. (d) φ^T > θ' and Ψ_{B_c}(φ^T) = Ψ⁰_B(φ^T) + α_pη/(α_p − φ^T) > κ; Ψ_{B_c} increases on (α_{p+1}, α_p).
(e) On the face {p, p+1}, T_c's modes are φ^T (vector e_p) and α_{p+1}; the coefficient of e_p in X is
G_X(φ^T)(α_p − φ^T)/α_p with G_X(φ^T) = −(φ'/L_B)Ψ_{B_c}(φ^T)/(φ^T − φ') > 0 by (d) and the pinch identity; the other
coefficient is X_{p+1}/(ê^T_{(p+1)})_{p+1} > 0, and (ê^T_{(p+1)})_{p+1} = Ψ_{T_c}(α_{p+1})/L_T > 0 since α_{p+1} < φ^T.
(f) Lemma 3.4 with k = p. ∎

**T-pinning** (π_p = T), η > 0: d = f_pe_p + ηe_{p+1}, b = (f_{p+1} − η)e_{p+1}, L_T = V_T − f_p − η, L_B = V_B + f_{p+1} − η.

**Proposition 4.6 (design.md Proposition 4.3; mirror image of 4.5).** Let κ > 0 and 0 < η < η̂_p,
η̂_p := min(f_{p+1}, V_T − f_p) (p+1 = N), η̂_p := min(f_{p+1}, (α_{p+1} − α_{p+2})Ψ⁰_T(α_{p+2})/α_{p+1}) (p+1 < N).
Then all phase flows are positive; B_c has the single major mode φ^B = α_{p+1}L_B/V_B ∈ (α_{p+1}, φ⁰_B]; T_c has exactly
the major modes φ^mid ∈ (φ⁰_T, α_p) and φ^low ∈ (α_{p+2}, α_{p+1}) (φ^low ∈ (0, α_N) if p+1 = N); margin
Ψ_{T_c}(φ^B) > κ > 0, i.e. φ^B < φ^mid; the forced pinch X := ê^{T_c}_{φ^low} > 0 on {p, p+1} and
X = C'e_{p+1} + D'ê^B_{(p)}, C', D' > 0, C' + D' = 1, e_{p+1} = ê^{B_c}_{φ^B}, ê^B_{(p)} := ê^{B_c}_{α_p} > 0 on {p, p+1}.

*Proof.* The mirror image of Proposition 4.5 (design.md Proposition 4.3 writes it out). ∎

### 4.3 The specification and the parameter region

**Definition 4.7 (specification F(δ, η, π); design.md Definition 5.1).** Data (α, z, q); δ > 0; η = ((η_c)_{c
non-binary}, (η_p)_{1≤p≤N−1}); a pattern π ∈ {B, T}^{N−1}.
(F1) Non-binary columns, by decreasing length: feed (f_c, U_c) by (4.1)'s definitions (for [1, N]: z, 1 − q); the
preferred split relative to (f_c, U_c) (Lemma 4.1), deformed with η_c (Definition after Lemma 4.2).
(F2) Binary columns: feeds and U's likewise; vapors from the reboiler up, V_{B[N−1,N]} := R := H* + δ,
V_{T[p,p+1]} := V_{B[p,p+1]} + U_{[p,p+1]}, V_{B[p−1,p]} := V_{T[p,p+1]} (zero-duty side draws); then d, b, L_T, L_B by
B- or T-pinning with η_p.
(F3) Products: P_1 := D_{[1,2]}, P_p := B_{[p−1,p]} + D_{[p,p+1]} (2 ≤ p ≤ N−1), P_N := B_{[N−1,N]}, where D_c := Σd_c,
B_c := Σb_c.
The **operating specification** consists of the stream flows only (V_σ, L_σ, P_p, the feed); the net component flows
d, b belong to the limit configuration (§4.4).

**Lemma 4.8 (balances; design.md Lemma 5.2).** If every non-binary column is valid: every feed f_c is > 0 on the range
of c and 0 outside; every junction and side-draw stage satisfies the total and CMO balances with zero duty; the
condenser removes V_{T[1,2]}, the reboiler receives R = H* + δ, the only heat input; V_{T[p,p+1]} = R + Σ_{r≥p}U_{[r,r+1]},
V_{B[p,p+1]} = R + Σ_{r>p}U_{[r,r+1]}; Σ_p P_p = 1 and the limit net flows satisfy every component balance
d_c + b_c = f_c, P_p x_p = n_p := b_{[p−1,p]} + d_{[p,p+1]}, Σ_p n_p = z. Moreover (*flux bookkeeping*), with Ψ_σ computed
from the limit net flows and φ_{P_p}(t) := −Σ_k α_k n_{p,k}/(α_k − t) the Underwood flux of the liquid product p, at every
unit u and every t ∉ {α_i}:

  Σ_{σ: bot(σ)=u} Ψ_σ(t) − Σ_{σ: top(σ)=u} Ψ_σ(t) + φ_{P_p}(t)[u = J[p,p]] = Q_u + φ_F(t)[u = J[1,N]],             (4.3)

Q_u the duty (0 at junctions and side draws, −V_{T[1,2]} at the condenser, R at the reboiler) and
φ_F(t) := 1 − q − Σ_k α_k z_k/(α_k − t).

*Proof.* Induction for the feeds (each parent contributes a vector positive exactly on [i, m]). Vapor at J_c:
V_{T_c} + V_{B[i−1,m]} − V_{B_c} − V_{T[i,m+1]} = V_{T_c} − V_{B_c} − U_c = 0 (feed: 1 − q). Liquid at J_c: L_{B_c} − L_{T_c} =
D_c + B_c − U_c = Σf_c − U_c = L_{B[i−1,m]} + q·[c = [1,N]] − L_{T[i,m+1]}. Side draws: V_{B[p−1,p]} = V_{T[p,p+1]} by (F2) and
L_{B[p−1,p]} − L_{T[p,p+1]} = B_{[p−1,p]} + D_{[p,p+1]} = P_p. Telescoping gives the vapors; summing d_c + b_c = f_c over the
tree of columns gives Σ_p n_p = z. (4.3): insert the definitions; at J_c the left side is (V_{T_c} − V_{B_c} − U_c) −
Σ_k α_k (d_{c,k} + b_{c,k} − f_{c,k})/(α_k − t) = 0 (for J[1, N]: φ_F(t)); at a side draw (V_{B[p−1,p]} − V_{T[p,p+1]}) +
Σ_k α_k(b_{[p−1,p],k} + d_{[p,p+1],k} − n_{p,k})/(α_k − t) = 0; at the condenser −Ψ_{T[1,2]} + φ_{P_1} = −V_{T[1,2]}; at the
reboiler Ψ_{B[N−1,N]} + φ_{P_N} = R. (This is class.md Lemma 5.1 summed over the stages of each section.) ∎

**Lemma 4.9 (middle roots are global; design.md Lemma 5.3).** If every non-binary column is valid, then θ^c_r = θ_r for
every non-binary c = [i, m] and i+1 ≤ r ≤ m−2; at η_nb := (η_c)_c = 0 every feed root of every column (binary
included) is global.

*Proof.* Induction top-down with (4.1). For a valid parent with deformation η ≥ 0, Ψ_T = Ψ⁰_T + ηK and Ψ_B = Ψ⁰_B + ηK
(Definition after Lemma 4.2), where Ψ⁰_T, Ψ⁰_B vanish at all the parent's feed roots (Lemma 4.1) and K at its middle ones
(Lemma 4.2); by induction these middle roots are global. So Ψ_{T[i,m+1]} vanishes at θ_{i+1}, …, θ_{m−1} and Ψ_{B[i−1,m]}
at θ_i, …, θ_{m−2}, hence Φ_{[i,m]} vanishes at θ_{i+1}, …, θ_{m−2}, and Φ has one root per interval. At η_nb = 0 both
section functions of every non-binary column vanish at all its feed roots, and the same induction makes every feed root
of every column global (HS's carry-over). ∎

**Proposition 4.10 (binary margins; design.md Proposition 5.4).** Let κ_p(η_nb; δ) be the margin (4.2) of [p, p+1] in F.
(a) The non-binary part of F and the binary feeds and U's do not depend on (δ, η_bin, π). (b) κ_p(η_nb; δ) =
κ_p(η_nb; 0) + δ. (c) κ_p(0; 0) = e_p. (d) There are η_1 > 0 and C_1 ≥ 0, depending only on (α, z, q), such that for
η_nb ∈ [0, η_1]^{(N−1)(N−2)/2} every non-binary column is valid and |κ_p(η_nb; 0) − e_p| ≤ C_1 max_c η_c.

*Proof.* (a) (F1) does not involve R, η_bin, π. (b) V_{T[p,p+1]} = R + const (Lemma 4.8), with const and the binary feed
independent of R. (c) (A light-side cut; design.md Proposition 5.4 (c).) At η_nb = 0 every feed root of every column is
global (Lemma 4.9); in particular θ_p is the feed root of [p, p+1] and a root of Ψ_{T[i,p+1]} for every non-binary
[i, p+1] (Lemma 4.1). κ_p does not depend on η_bin (by (a)), so we may evaluate the flux identities (4.3) with the sharp
binary net flows (η_bin = 0), for which n_r = z_r e_r for every product r (Lemma 4.8). Sum (4.3) at t = θ_p over the
*light* units: the junctions J[i, m] with m ≤ p, the condenser and the side draws J[2,2], …, J[p,p]. A section with both ends
light contributes Ψ − Ψ = 0; a bottom section B[i, m] with a light end has m ≤ p, so both ends are light; a top section
T[i, m] with a light end has m ≤ p + 1; so the sections with exactly one light end are T[i, p+1], 1 ≤ i ≤ p (light top end
J[i, p]). The feed enters J[1, N], which is not light (p ≤ N − 1). Hence −Σ_{i=1}^p Ψ_{T[i,p+1]}(θ_p) +
Σ_{r=1}^p φ_{P_r}(θ_p) = −V_{T[1,2]}. Summing (4.3) over all units, every section cancels and Σ_r φ_{P_r}(t) =
R − V_{T[1,2]} + φ_F(t); letting t → ∞ gives V_{T[1,2]} = R + 1 − q. With Ψ_{T[i,p+1]}(θ_p) = 0 for i < p,
φ_{P_r}(θ_p) = −α_r z_r/(α_r − θ_p) and the feed equation Σ_r α_r z_r/(α_r − θ_p) = 1 − q:

  κ_p(0; δ) = Ψ_{T[p,p+1]}(θ_p) = R + (1 − q) − Σ_{r≤p} α_r z_r/(α_r − θ_p) = R − Σ_{r>p} α_r z_r/(θ_p − α_r) = R − H_p = e_p + δ.

(This is the light-side cut identity of paper.tex Proposition `cut` and re-proves the side-draw balances of Theorem
`achieve`.) (d) All data of
every column are real-analytic in η_nb near 0 (simple roots, rational expressions with nonvanishing denominators);
(V) is strict at η_nb = 0; compactness and the mean value theorem. ∎

**Definition (admissible parameters; design.md §5).** (δ, η, π) with δ > 0 and all η's > 0 is *admissible* if
(i) every non-binary column is valid (V) along the recursion (F1); (ii) κ_p > 0 for every binary column; (iii)
0 < η_p < η̂_p for every binary column (bounds of Propositions 4.5–4.6, computed from F). These are finitely many
explicit inequalities on the flows of F; Propositions 4.3, 4.5, 4.6 then apply to every column, and only (i)–(iii) are
used below.

**Theorem 4.11 (the admissible set is nonempty for every δ; design.md Theorem 5.5).** Let δ > 0,
η̄(δ) := min(η_1, δ/(2C_1)) (η̄ = η_1 if C_1 = 0). If 0 < η_c ≤ η̄(δ) for every non-binary column, then (i) holds and
κ_p ≥ e_p + δ/2 > 0 for every binary column; any η_p ∈ (0, η̂_p) then gives an admissible parameter. The admissible set
is open in η, and every quantity below is real-analytic on it.

*Proof.* κ_p = κ_p(η_nb; 0) + δ ≥ e_p − C_1η̄(δ) + δ ≥ e_p + δ/2. ∎

### 4.4 The limit configuration

**Definition 4.12 (limit configuration 𝓛(δ, η, π); design.md Definition 6.1).** For admissible parameters: net flows
d_c, b_c as in F (non-binary: (F1); binary: B- or T-pinning); every component outside a column's range has zero net
flow in its sections. Junctions: X_c of Proposition 4.3 (e) (non-binary); X_c := ê^{B_c}_{φ'} (B-pinned) or
ê^{T_c}_{φ^low} (T-pinned). Products: x_p := n_p/P_p, n_p := b_{[p−1,p]} + d_{[p,p+1]}. Supports: supp X_c = [i, m];
supp x_p = Maj(p) := {p} ∪ {p−1 : π_{p−1} = B} ∪ {p+1 : π_p = T}; all entries on the supports are positive. The point
v_∞ := (d_c, X_c)_c ∈ 𝒱.

**Table 1 (modes, pinch, closed ends; design.md Table 1, limit.md §3.1).**

| type | section | modes (decreasing) | pinch r* | closed at bottom (above r*) | closed at top (below r*) |
|---|---|---|---|---|---|
| 1 (T-nb) | T_c, c = [i,m] non-binary | α_k (k<i) > r_c > θ_{i+1} > ⋯ > θ_{m−2} > u_c > α_m > α_k (k>m) | u_c | α_k (k<i), r_c, θ's | α_m, α_k (k>m) |
| 2 (B-nb) | B_c, c = [i,m] non-binary | α_k (k<i) > α_i > v_c > θ's > s_c > α_k (k>m) | v_c | α_k (k≤i) | θ's, s_c, α_k (k>m) |
| 3 (T-np) | T_c, binary, π_p = B | α_k (k<p) > φ^T > α_k (k≥p+1) | φ^T | α_k (k<p) | α_k (k≥p+1) |
| 4 (B-f) | B_c, binary, π_p = B | α_k (k<p) > φ' > φ'' > α_k (k>p+1) | φ' | α_k (k<p) | φ'', α_k (k>p+1) |
| 5 (T-f) | T_c, binary, π_p = T | α_k (k<p) > φ^mid > φ^low > α_k (k>p+1) | φ^low | α_k (k<p), φ^mid | α_k (k>p+1) |
| 6 (B-np) | B_c, binary, π_p = T | α_k (k≤p) > φ^B > α_k (k>p+1) | φ^B | α_k (k≤p) | α_k (k>p+1) |

Every top section pinches at its smallest major mode and every bottom section at its largest; every trace mode lies
above all major modes or below all of them.

**Theorem 4.13 (properties of the limit configuration; design.md Theorem 6.2).** For admissible parameters:
(H1) every V_σ, L_σ, D_c, B_c, P_p is positive and every feed is positive on its column's range;
(H2) top sections have w ≥ 0 and bottom sections w ≤ 0, with the supports of Definition 4.12; the forced flows are
nonzero; middle roots are global (Lemma 4.9);
(H3), (H4) the conclusions of Propositions 4.3, 4.5, 4.6, with κ_p ≥ e_p + δ/2;
(H5) Table 1 holds: in every section all N modes are positive and distinct, the pinch is a major mode, and
ρ_σ := max over consecutive modes of (smaller)/(larger) < 1; **limit relations**: for every section, ℓ_φ·x_bot = 0
for every mode above r* and ℓ_φ·x_top = 0 for every mode below r*;
(H6) **half-profiles**: for every section, the vectors M^k x_bot (k ≥ 0) and M^{−k} x_top (k ≥ 0) are ≥ 0, ≠ 0, and the
chart vectors v^bot_k := M^k x_bot/((Vr*)^k c_{r*}(x_bot)) = ê_{r*} + Σ_{φ<r*}(c_φ(x_bot)/c_{r*}(x_bot))(φ/r*)^k ê_φ and
v^top_k := M^{−k}x_top (Vr*)^k/c_{r*}(x_top) = ê_{r*} + Σ_{φ>r*}(c_φ(x_top)/c_{r*}(x_top))(r*/φ)^k ê_φ are ≥ 0 with
inf_k 1ᵀv_k > 0 (c_φ as in Lemma 3.1 (iv));
(H7) for every trace mode α_k of σ closed at J: Ψ^{(k)}_σ(α_k) > 0, k is lighter (bottom end) resp. heavier (top end)
than every component of supp x_J, so G^{(k)}_{x_J}(α_k) < 0 resp. > 0, and x_{J,k} = 0; every junction J_c at which k is
a trace is the closed end, for k, of T_c (k lighter than c's range) or B_c (k heavier);
(H8) R = H* + δ is the only heat input; the impurity of product p is 1 − x_{p,p} = (η_{p−1}[π_{p−1} = B] +
η_p[π_p = T])/P_p; the pivots c_{r*}(x_bot), c_{r*}(x_top) are > 0 for every section (Underwood form: at a product end
G_x(r*) > 0 for a top section and < 0 for a bottom one; at a child junction G_{X[i,m]}(u_{[i,m+1]}) > 0,
G_{X[i,m]}(v_{[i−1,m]}) < 0; at a non-binary junction the pivots are B_c and C_c; at a binary junction 1 for the pinning
section and B' resp. C' for the other).

*Proof.* (H1)–(H4): Lemma 4.8, Propositions 4.3–4.6, Theorem 4.11. (H5): the mode lists are in Propositions 4.3 (c),
4.5 (b, c), 4.6. The limit relations, type by type (ℓ_φ·x = G_x(φ) for a major mode, x_k for a trace mode):
type 1: G_{X_c}(r_c) = G_{X_c}(θ_r) = 0 by Proposition 4.3 (e), X_{c,k} = 0 (k < i), and the top end X_{[i,m−1]} is
supported on [i, m−1]; type 2: the bottom end X_{[i+1,m]} is supported on [i+1, m], G_{X_c}(θ_r) = G_{X_c}(s_c) = 0,
X_{c,k} = 0 (k > m); type 3: X_c supported on {p, p+1}; the top end x_p has Maj(p) ⊆ {p−1, p}; type 4: x_{p+1} has
Maj(p+1) ⊆ {p, p+1, p+2}; X_c = ê_{φ'} is annihilated by ℓ_{φ''} (Lemma 3.1 (iii)) and has no k > p+1; types 5, 6: the
mirror images (Maj(p) ⊆ {p−1, p, p+1}, Maj(p+1) ⊆ {p+1, p+2}). *Pivots.* At a junction, the pivot of the section whose
pinch vector appears in the cone decomposition is its coefficient (B_c, C_c, B', 1, 1, C', all > 0). At the other end x
(a product or a child's junction) supp x lies entirely on one side of r* (top sections: lighter; bottom: heavier), so
G_x(r*) has a sign, and ℓ_{r*}·ê_{r*} = (r*/L)(−Ψ'(r*)) has the same sign (Ψ decreasing for top, increasing for bottom
sections): c_{r*}(x) > 0. (H6): going up a bottom section (types 2, 4, 6) M ≥ 0; going down a top section (types 1, 3, 5)
M^{−1} ≥ 0 (Lemma 3.3); the two other directions start at a junction and use the cone decompositions (types 1, 3: the
decompositions of Propositions 4.3 (e), 4.5 (e) in T_c's basis; type 5: X_c = ê_{φ^low}; type 2: in B_c's basis;
type 4: X_c = ê_{φ'}; type 6: Proposition 4.6) with Lemma 3.3 (c). By the limit relations, M^k x_bot =
c_{r*}(x_bot)(Vr*)^k v^bot_k and similarly at the top; since c_{r*} > 0 and M^{±k}x ≥ 0, ≠ 0 (M invertible), v_k ≥ 0 and
1ᵀv_k > 0; 1ᵀv_k → 1ᵀê_{r*} = 1, so inf_k 1ᵀv_k > 0. (H7): from Table 1 and the supports, k is lighter than the whole
support at a bottom closed end and heavier at a top closed end; a lighter trace of a top section exceeds every pole
(Ψ > V, Lemma 3.2 (a)); a heavier trace of a top section lies below ψ_low (Ψ > 0); a lighter trace of a bottom section
lies above ψ_high (Ψ > 0); a heavier trace of a bottom section lies below every pole (Ψ > V); and Ψ^{(k)} = Ψ since
w_k = 0. The statement on junctions is read off Table 1. (H8): Lemma 4.8 and the entries of n_p. ∎

---

## 5. The limit system and its nondegeneracy (Agent L, limit.md)

**5.1 Hypotheses.** Theorem L uses only the following properties (limit.md §3.2), all contained in Theorem 4.13:
(L1) = (H1) and the consistency of F; (L2) = the limit net flows of (H2) (non-binary [i, m]: d_{c,k} > 0 for
i ≤ k ≤ m−1, b_{c,k} > 0 for i+1 ≤ k ≤ m, all other entries 0; binary as in Definition 4.12) and the common middle roots;
(L3) *pinch separation*: u_c > α_m, v_c < α_i (non-binary), φ^T > α_{p+1}, φ' < α_{p−1} (π_p = B), φ^B < α_p,
φ^low > α_{p+2} (π_p = T), whenever these α's exist; (L4) *margins*: r_c ≠ v_c, u_c ≠ s_c, φ^T ≠ φ'' (π_p = B),
φ^mid ≠ φ^B (π_p = T). From (L1)–(L3) alone L re-derives (limit.md Proposition 3.3) the spectral structure of Table 1,
the closed ends at trace nodes, κ > 0 and, with (L4), nonzero pivots; here these facts are already in Theorem 4.13.

For a component k let 𝒥_k := {nodes J : x^∞_{J,k} = 0} = {J[i,m] : i < m, k ∉ [i,m]} ∪ {J[p,p] : k ∉ Maj(p)}. Every
section attached to a node of 𝒥_k has zero limit net flow of k, and every J ∈ 𝒥_k has either first index > k ("k
lighter than J", J ∈ 𝒥^<_k) or second index < k (J ∈ 𝒥^>_k) (limit.md Lemma 3.4).

**Lemma 5.1 (limit.md Lemma 3.5).** G_∞(v_∞) = 0. *Proof.* These are the limit relations of (H5) (Lemma 3.5 (b)). ∎

**Lemma 5.2 (domain; limit.md Lemmas 3.6–3.7).** There are a compact neighborhood 𝒩 of v_∞ and ρ̄ < 1, π_0 > 0 such
that on 𝒩, for every σ, the modes of M_σ(w_σ(v)) are simple, positive, ordered as at v_∞, analytic, with consecutive
ratios ≤ ρ̄, and |π_bot|, |π_top| ≥ π_0; G_∞ and all G_n are analytic on 𝒩, the exact system is equivalent to G_n = 0
there, and |G_n − G_∞| ≤ Cρ̄^{n+1}, ‖D(G_n − G_∞)‖ ≤ C(n + 2)ρ̄^{n+1} with n = min_σ n_σ. *Proof.* Theorem 4.13 (H5),
(H8), Lemma 3.4, continuity, Lemma 3.5 (a), (c). ∎

**5.2 Coordinates.** On 𝒱 use the affine coordinates of Table 2 (limit.md §4.1); all other entries of (d_c, X_c) are
dependent (Σ_k d_{c,k} = Δ_{T_c}, Σ_k X_{c,k} = 1, b_c = f_c − d_c given the parents). The change of variables is affine,
triangular in the order of the columns (parents first) and invertible.

**Table 2.**

| column | trace coordinates (limit value 0) | major coordinates | dependent |
|---|---|---|---|
| non-binary [i,m] | d_{c,k} (k ∉ [i,m]), d_{c,m}, b_{c,i}; X_{c,k} (k ∉ [i,m]) | d_{c,p} (i < p < m−1); X_{c,p} (i < p ≤ m) | d_{c,i} = f_{c,i} − b_{c,i}, d_{c,m−1}, X_{c,i} |
| binary, π_p = B | d_{c,k} (k ≠ p); X_{c,k} (k ∉ {p,p+1}) | X_{c,p+1} | d_{c,p}, X_{c,p} |
| binary, π_p = T | d_{c,k} (k ∉ {p,p+1}), b_{c,p}; X_{c,k} (k ∉ {p,p+1}) | X_{c,p+1} | d_{c,p} = f_{c,p} − b_{c,p}, d_{c,p+1}, X_{c,p} |

A trace coordinate carries the index k of its component. In this section the rows of G_∞ are taken divided by the
nonvanishing pivot product π_bot π_top (limit.md (2.4): the row of a mode r closed at the end J is ℓ_r·x_J/ℓ_{r*}·x_J);
by Lemma 3.5 (d) this does not affect the invertibility of the derivative at v_∞. Left eigenvectors: Underwood's ℓ_φ for
major modes, the regularized ℓ̃ of Lemma 3.4 for trace modes. The rows are **major rows** (non-pinch major modes) and
**trace rows of component k** (the trace mode α_k of a section σ with w^∞_{σ,k} = 0).

**Lemma 5.3 (counts; limit.md Lemma 4.1).** For every column c and component k, the trace rows of k carried by T_c and
B_c are as many as the trace coordinates of k of c; the major rows of c are as many as its major coordinates
(2(m − i − 1) for a non-binary column, 1 for a binary one).

*Proof.* k ∉ range(c): two trace rows (T_c, B_c), two coordinates (d_{c,k}, X_{c,k}). Non-binary, k = i: one row (B_c),
one coordinate (b_{c,i}); k = m: one row (T_c), coordinate d_{c,m}; i < k < m: none. Binary π_p = B: k = p+1: one row
(T_c), coordinate d_{c,p+1}; π_p = T: k = p: one row (B_c), coordinate b_{c,p}. Major rows: r_c and the m − i − 2 middle
modes in T_c, the middle modes and s_c in B_c; binary: φ'' resp. φ^mid. Major coordinates: m − i (X) + m − i − 2 (d),
resp. 1. ∎

**Lemma 5.4 (trace quantities are linear in the traces of their component; limit.md Lemma 4.2).** Fix k. For every
section σ with w^∞_{σ,k} = 0, w_{σ,k}(v) is a linear form in the trace coordinates of component k only; for every
J ∈ 𝒥_k, x_{J,k}(v) is such a linear form.

*Proof.* Induction over the columns, parents first: d_{c,k} with d^∞_{c,k} = 0 is a trace coordinate; b_{c,k} with
b^∞_{c,k} = 0 is either a trace coordinate (b_{c,i}, b_{c,p}) or, for k ∉ range(c), d_{[i,m+1],k} + b_{[i−1,m],k} − d_{c,k}
with both parental limit flows zero (they are supported on [i, m]). x_{J_c,k} = X_{c,k}; x_{p,k} = (b_{[p−1,p],k} +
d_{[p,p+1],k})/P_p. ∎

**Lemma 5.5 (linearized trace rows; limit.md Lemma 4.3).** Let α_k be a trace mode of σ closed at J (J ∈ 𝒥_k). With
q^k_{J→σ} the net flow of k from J into σ (q = w_{σ,k} if J = bot(σ), −w_{σ,k} if J = top(σ)) and κ^{(k)}_σ :=
|G^{(k)}_{x^∞_J}(α_k)|/Ψ^{(k)}_σ(α_k) > 0 (Theorem 4.13 (H7)), the derivative at v_∞ of the row
ℓ̃_k·x_J/ℓ_{r*}·x_J is

  ( δx_{J,k} − κ^{(k)}_σ δq^k_{J→σ} ) / ℓ_{r*}·x^∞_J.

In particular the trace rows of k depend, at first order at v_∞, only on the trace coordinates of k.

*Proof.* By (3.1) the numerator is x_{J,k} + w_{σ,k} g(v) with g analytic; x_{J,k} and w_{σ,k} vanish at v_∞, so its
derivative is δx_{J,k} + g(v_∞)δw_{σ,k}, g(v_∞) = G^{(k)}_{x^∞_J}(α_k)/Ψ^{(k)}_σ(α_k); the signs of (H7) turn g δw into
−κ δq at both kinds of ends. Lemma 5.4 gives the last statement. ∎

**Lemma 5.6 (node identities; limit.md Lemma 4.4).** For every component k and every node J ≠ J[1, N], identically on 𝒱,
Σ_{σ attached to J} q^k_{J→σ}(v) + P_J x_{J,k}(v) = 0 (P_J := 0 for a junction). *Proof.* The junction balance
d_c + b_c − f_c = 0 and the definition of x_p. ∎

**Proposition 5.7 (trace block; limit.md Proposition 4.5).** For every k, the square matrix 𝕋_k := ∂(trace rows of
k)/∂(trace coordinates of k) at v_∞ is nonsingular. More precisely, if 𝕋_k δt = 0, the contents ξ_J := δx_{J,k}
(J ∈ 𝒥_k) solve A_k ξ = 0 with

  (A_k)_{JJ} = P_J + Σ_{σ: closed at J} 1/κ^{(k)}_σ,   (A_k)_{JJ'} = −Σ_{σ: closed at J', open end J} 1/κ^{(k)}_σ  (J ≠ J'),

and A_k is, after a permutation, upper triangular with positive diagonal; it is a nonsingular M-matrix (A_k^{−1} ≥ 0).

*Proof.* Extend δt by zero. By Lemma 5.5 every trace row of k gives ξ_{c(σ)} = κ_σ δq_σ (c(σ) the closed end,
δq_σ := δq^k_{c(σ)→σ}). At J ∈ 𝒥_k every attached section has zero limit net flow of k, and the node identity (Lemma
5.6, linear) gives Σ_{σ closed at J} ξ_J/κ_σ − Σ_{σ open at J} ξ_{c(σ)}/κ_σ + P_J ξ_J = 0 (an open end's flow is
−δq_σ; the closed end of such σ is in 𝒥_k by (H7)), i.e. A_kξ = 0. *Triangularity:* if k is a lighter trace of σ,
c(σ) = bot(σ), the open end is top(σ), and h(top) = h(bot) − 1 (top of T[i,m] is J[i,m−1], bottom of B[i,m] is
J[i+1,m]); if k is heavier, c(σ) = top(σ) and h(bot) = h(top) + 1; both ends lie in the same class 𝒥^<_k or 𝒥^>_k.
Ordering 𝒥^<_k by increasing h and then 𝒥^>_k by decreasing h, every nonzero off-diagonal entry is above the
diagonal. *Diagonal:* P_J > 0 at products; a junction J[i,m] ∈ 𝒥^<_k is the closed bottom end of T[i,m] for k (k < i),
and J[i,m] ∈ 𝒥^>_k the closed top end of B[i,m] (k > m) ((H7)); κ > 0. A triangular Z-matrix with positive diagonal is
a nonsingular M-matrix. *Conclusion:* ξ = 0, so δq_σ = 0 for all σ ∈ Σ_k; every trace coordinate of k is some ξ_J or
± some δq_σ; δt = 0. ∎

**Proposition 5.8 (major block; limit.md Proposition 4.6).** Order the columns by decreasing length. The matrix
𝕄 := ∂(major rows)/∂(major coordinates) at v_∞ is block lower triangular over the columns, with nonsingular diagonal
blocks 𝕄_c; hence 𝕄 is nonsingular.

*Proof.* *Triangularity:* the major rows of c are closed at J_c and depend on X_c, d_c and f_c, i.e. on the coordinates
of c and of its ancestors. *Non-binary c = [i, m]:* let δ be a variation of the major coordinates of c only with
𝕄_c δ = 0. Then δf_c = 0, δd_{c,i} = δd_{c,m} = 0 = δd_{c,k} (k ∉ [i, m]), so δd_c is supported on {i+1, …, m−1} with
zero sum and δw_{T_c} = δw_{B_c} = δd_c. For a root χ of Ψ_{T_c} or Ψ_{B_c}, Ψ'(χ)δχ = R(χ) := Σ_p α_pδd_p/(α_p − χ). Each
major row is G_{X_c}(χ)/(pivot) with G_{X_c}(χ) = 0, of variation [ℓ_χ·δX_c + G'_{X_c}(χ)δχ]/(pivot). The two rows of a
middle root θ_j (in T_c and B_c) give G'_{X_c}(θ_j)R(θ_j)(1/Ψ'_{T_c}(θ_j) − 1/Ψ'_{B_c}(θ_j)) = 0, where Ψ'_{T_c} < 0 < Ψ'_{B_c}
(Lemma 3.2) and G'_{X_c}(θ_j) ≠ 0 (simple zero of G_{X_c}, Proposition 4.3 (e)); so R(θ_j) = 0 at the m − i − 2 middle
roots, and R(0) = Σ_pδd_p = 0. R(t)Π_{p=i+1}^{m−1}(α_p − t) is a polynomial of degree ≤ m − i − 2 with m − i − 1 distinct
zeros, hence 0: δd_c = 0. (This is the Cauchy-type matrix (α_p/(α_p − θ_j))(1/Ψ'_T(θ_j) − 1/Ψ'_B(θ_j)) of BRIEF-common,
invertible on {Σδd = 0}.) Then δχ = 0 for every mode, and the rows of r_c, θ_j, s_c say ℓ_χ·δX_c = 0 at the m − i
interlacing points with Σ δX_c = 0, δX_c supported on [i, m]: δX_c = 0 by Proposition 4.3 (e). *Binary c, π_p = B:* the
only major coordinate is X_{c,p+1}; with everything else fixed δd_c = δb_c = 0 and the modes do not move; the row of φ''
varies by (α_{p+1}/(α_{p+1} − φ'') − α_p/(α_p − φ''))δX_{c,p+1}/(pivot) ≠ 0 since α_{p+1} < φ'' < α_p. π_p = T: the same
with φ^mid ∈ (α_{p+1}, α_p). ∎

**Theorem 5.9 (Theorem L; limit.md §4.4).** Under (L1)–(L4) — in particular for every admissible F(δ, η, π) — G_∞ is
real-analytic near v_∞, G_∞(v_∞) = 0, and DG_∞(v_∞) is invertible: in the coordinates of Table 2, with rows ordered
(major; trace of k = 1, …, N) and coordinates likewise,

  DG_∞(v_∞) = [ 𝕄  ∗ ; 0  𝕋 ],   𝕋 = diag(𝕋_1, …, 𝕋_N),   det DG_∞(v_∞) = det 𝕄 · Π_k det 𝕋_k ≠ 0.

*Proof.* Lemma 5.1, Lemma 5.2; the zero block and the block-diagonal form of 𝕋 by Lemma 5.5; Propositions 5.7, 5.8; the
diagonal blocks are square by Lemma 5.3. ∎

(For N = 3 and N = 4 every row, coordinate and block is written out in limit.md §§5–6.)

---

## 6. From the limit system to finite-stage steady states (Agent A, analysis.md)

Let F be admissible and v_∞ = 𝓛. In the notation of analysis.md, the hypotheses used here are: (A1) F as in §2.3;
(A2) N distinct positive modes in every section and node compositions in Δ; (A3) nonzero pivots; (A4) G_∞(v_∞) = 0;
(A5) the half-profiles ξ^bot_k := v^bot_k, ξ^top_k := v^top_k of (H6) are ≥ 0 (with 1ᵀê_{r*} = 1 and c_{r*} =
ℓ_{r*}·x/ℓ_{r*}·ê_{r*}); (A6) DG_∞(v_∞) invertible. (A1) ⇐ Lemma 4.8; (A2) ⇐ (H5), Definition 4.12; (A3) ⇐ (H8);
(A4) ⇐ Lemma 5.1; (A5) ⇐ (H6); (A6) ⇐ Theorem 5.9.

**Lemma 6.1 (quantitative inverse function step; analysis.md Lemma A.4).** Let E be a finite-dimensional normed space,
F_∞, F ∈ C¹(B̄; E) on B̄ := B̄(v_∞, r_0), F_∞(v_∞) = 0, A := DF_∞(v_∞) invertible, ω(r) := sup_{|v−v_∞|≤r}‖DF_∞(v) − A‖,
and r_1 ∈ (0, r_0] with ‖A^{−1}‖ω(r_1) ≤ 1/4. If sup_{B̄}‖DF − DF_∞‖ ≤ 1/(4‖A^{−1}‖) and ‖A^{−1}‖|F(v_∞)| ≤ r_1/2, then F
has exactly one zero v* in B̄(v_∞, r_1), and |v* − v_∞| ≤ 2‖A^{−1}‖|F(v_∞)|.

*Proof.* Φ(v) := v − A^{−1}F(v) has ‖DΦ(v)‖ ≤ ‖A^{−1}‖(‖A − DF_∞(v)‖ + ‖DF_∞(v) − DF(v)‖) ≤ 1/2 on the convex ball
B̄(v_∞, r_1), which Φ maps into itself (|Φ(v) − v_∞| ≤ |Φ(v) − Φ(v_∞)| + ‖A^{−1}‖|F(v_∞)| ≤ r_1/2 + r_1/2); Banach's
fixed-point theorem; |v* − v_∞| ≤ |v* − v_∞|/2 + ‖A^{−1}‖|F(v_∞)|. ∎

**Theorem 6.2 (network perturbation; analysis.md Theorem A.5, limit.md Corollary 4.7).** There are r_1 > 0, n_0 and C
such that for all section lengths with n := min_σ n_σ ≥ n_0, G_n has exactly one zero v_n in B̄(v_∞, r_1), and

  |v_n − v_∞| ≤ C ρ_∞^{n+1},   ρ_∞ := max_σ max(λ^∞_{σ,r*+1}/λ^∞_{σ,r*}, λ^∞_{σ,r*}/λ^∞_{σ,r*−1}) < 1.

At v_n every section satisfies x_top = c_σ M_σ^{n_σ+1}x_bot with c_σ ≠ 0.

*Proof.* On 𝒩 (Lemma 5.2), sup‖DG_n − DG_∞‖ ≤ C_1 max_σ(n_σ + 2)ρ̄^{n_σ+1} ≤ C_1(n + 2)ρ̄^{n+1} (n ≥ 1/log(1/ρ̄), where
x ↦ (x + 2)ρ̄^{x+1} decreases), and |G_n(v_∞)| = |(G_n − G_∞)(v_∞)| ≤ C_2ρ_∞^{n+1} by Lemma 3.5 (c) at w = w^∞. G_∞ is C²,
so ω(r) ≤ C_3r; take r_0 with B̄(v_∞, r_0) ⊂ 𝒩, r_1 := min(r_0, 1/(4C_3‖A^{−1}‖)), n_0 with C_1(n_0 + 2)ρ̄^{n_0+1} ≤
1/(4‖A^{−1}‖) and ‖A^{−1}‖C_2ρ_∞^{n_0+1} ≤ r_1/2, and apply Lemma 6.1 to F = G_n. The last claim is Lemma 3.5 (a). ∎

**Lemma 6.3 (exact composite profile; analysis.md Lemma A.6).** At v_n, take for a section σ the eigen-data at
w_σ(v_n) with the right eigenvectors ê_r normalized by 1ᵀê_r = 1 near v_∞ (possible since 1ᵀê^∞_r = 1, Lemma 3.1 (ii)) and
the left eigenvectors rescaled so that ℓ_r·ê_r = 1 (so ℓ_r·x = c_r(x) and π = c_{r*}); put a_s := ℓ_s·x_bot/π_bot
(s ∈ D), b_r := ℓ_r·x_top/π_top (r ∈ U). For 0 ≤ k ≤ n_σ + 1:

  u_k := M^k x_bot/(π_bot λ_{r*}^k) = ê_{r*} + Σ_{s∈D} a_s μ_s^k ê_s + Σ_{r∈U} b_r μ_r^{n_σ+1−k} ê_r,               (6.1)

and the stage compositions of Proposition 2.1 are x_k = u_k/1ᵀu_k.

*Proof.* M^k x_bot = Σ_r (ℓ_r·x_bot)λ_r^k ê_r (ℓ_r·ê_r = 1), and E_r = 0 (r ∈ U) gives ℓ_r·x_bot = μ_r^{n_σ+1}b_rπ_bot. ∎

The profile is a bottom half-profile (decaying from the bottom end at the rates μ_s) plus a top half-profile
(decaying from the top end), glued at the pinch; at v_∞ these are v^bot_k and v^top_k of (H6).

**Lemma 6.4 (S > 0 at every stage; analysis.md Lemma A.7).** Let c_0 := min_σ min(inf_k 1ᵀv^bot_k, inf_k 1ᵀv^top_k) > 0
((H6)). There are n_1 ≥ n_0 and C such that for min_σ n_σ = n ≥ n_1, every σ and 0 ≤ k ≤ n_σ + 1,

  |u_k − v^bot_k|_∞ ≤ C(ρ_∞^{n+1} + ρ̄^{n_σ+1−k}),   |u_k − v^top_{n_σ+1−k}|_∞ ≤ C(ρ_∞^{n+1} + ρ̄^k),
  1ᵀu_k ≥ c_0/2,   S(x_k) ≥ α_N c_0/(2·1ᵀu_k) > 0.

*Proof.* The coefficients ê_r, μ_r, a_s, b_r of (6.1) are analytic in v on B̄(v_∞, r_1) and equal at v_∞ to those of
(H6), so they differ from them by ≤ C|v_n − v_∞| ≤ Cρ_∞^{n+1}; since all μ ≤ ρ̄, |μ^k − (μ^∞)^k| ≤ (sup_j jρ̄^{j−1})|μ − μ^∞|
uniformly in k. Termwise comparison of (6.1) with v^bot_k: the ê_{r*}- and D-terms differ by ≤ Cρ_∞^{n+1} and the U-terms
are ≤ Cρ̄^{n_σ+1−k} (b_r bounded since |π_top| ≥ π_0); with v^top_{n_σ+1−k}: the U-terms differ by ≤ Cρ_∞^{n+1} and the
D-terms are ≤ Cρ̄^k. Since min(ρ̄^k, ρ̄^{n_σ+1−k}) ≤ ρ̄^{(n+1)/2}, for every k one half-profile ξ (≥ 0, 1ᵀξ ≥ c_0) satisfies
|u_k − ξ|_∞ ≤ ε_n := C(ρ_∞^{n+1} + ρ̄^{(n+1)/2}); then 1ᵀu_k ≥ c_0 − Nε_n and αᵀu_k ≥ α_Nc_0 − (Σ_iα_i)ε_n; take n_1 with
Nε_n ≤ c_0/2 and (Σ_iα_i)ε_n ≤ α_Nc_0/2. ∎

**Theorem 6.5 (finite-stage steady states near the limit configuration; analysis.md Theorem A.8).** For every
admissible (δ, η, π) there are n_1 and C such that for all section lengths with min_σ n_σ = n ≥ n_1, the arrangement
operated with F(δ, η, π) has a steady state with
(i) every mole fraction of every unit and stream strictly positive;
(ii) every node composition and column split within Cρ_∞^{n+1} of its value in 𝓛 — in particular every product p of
purity ≥ 1 − (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p − Cρ_∞^{n+1};
(iii) stage compositions within C(ρ_∞^{n+1} + ρ̄^{min(k, n_σ+1−k)}) of the limit half-profiles;
(iv) heat supply R = H* + δ;
(v) local uniqueness: no other steady state has its (d_c, X_c) in B̄(v_∞, r_1).

*Proof.* Theorem 6.2 gives v_n; Lemma 6.4 gives S > 0 at every stage and 1ᵀu_k ≠ 0; S > 0 at the node compositions
because they are within Cρ_∞^{n+1} of points of Δ (shrink r_1 so that S ≥ α_N/2 on the ball). Proposition 2.1 (b) gives a
steady state with (i) and (iv) (Lemma 4.8); (ii) by Theorem 6.2 and (H8); (iii) by Lemma 6.4; (v) by Theorem 6.2 and
Proposition 2.1 (a). ∎

---

## 7. Proof of Theorem 1

Let N ≥ 3 (N = 2: sections.md Proposition 10), δ > 0, ε > 0, and fix a pattern π (e.g. all B). Choose η_c := η̄(δ) for
every non-binary column (Theorem 4.11). The bounds η̂_p > 0 of Propositions 4.5–4.6 involve only the binary feeds and
vapors, which depend on (δ, η_nb) but not on η_bin (Proposition 4.10 (a), Lemma 4.8); the product flows P_p are affine in
η_bin and positive on the admissible set, with P_p → n_{p,p} > 0 as η_bin → 0. So we can choose η_p ∈ (0, η̂_p) so small
that (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p ≤ ε/2 for every p. The parameters are
admissible, so Theorem 6.5 applies: take n_0 ≥ n_1(δ, η, π) with C(δ, η, π)ρ_∞^{n_0+1} ≤ ε/2. For all section lengths
n_σ ≥ n_0, the arrangement operated with F(δ, η, π) — a specification independent of the n_σ — has a steady state with all
mole fractions positive, every product of purity ≥ 1 − ε, and R = H* + δ. Hence R*(𝒮_Petlyuk; ε) ≤ H* + δ for all
δ, ε > 0, i.e. R*(𝒮_Petlyuk) ≤ H*. ∎

**Order of limits and constants.** n → ∞ at fixed (δ, η), then η → 0, then δ → 0. The constants n_0, C, r_1 depend on
(δ, η) and are not uniform as η → 0: the non-binary pivots B_c, C_c are O(η_c) and the bottom pivot of a forced section
is O(η_p) (limit.md Remark 7.3), so n_0 ≳ log(1/η)/log(1/ρ_∞), the same threshold as Astra's prefractionator theorem.
No uniformity is needed, because R*(𝒮; ε) is an infimum over all stage numbers.

---

## 8. Remarks

**8.1 Section lengths.** Theorem 1 holds for any family of section lengths with min_σ n_σ ≥ n_0. The limit system does not
involve the lengths because every section of 𝓛 is singly pinched (H8); with a doubly pinched section the position of a
composition front depends on the length ratios (network.md §11.1).

**8.2 Why forcing is needed (limit.md §7.1).** If a binary column is not forced (η_p = 0), both of its sections are
non-pinning; every row closed at J_c is a trace row, which does not involve the ratio X_{c,p} : X_{c,p+1} when the traces
vanish, so v_∞ + t(e_{p+1} − e_p) (in X_c) is a one-parameter family of zeros of G_∞: DG_∞(v_∞) is singular (one junction
ratio without an equation, one surplus trace row of component p). This is the neutral mode of network.md §10.4 and the
free face directions of Astra's rank count (attainment-astra §5). Forcing turns p into a major of the forced section, whose
pinch then fixes X_c, while its second major root supplies the missing major row.

**8.3 Why δ > 0 and η_c > 0 are needed.** At δ = 0 a binary column at the highest peak has margin κ_p = 0, (L4) fails
(φ^T = φ'' in the unforced limit), and X_c is the pinch of both of its sections (doubly pinched). At η_c = 0 the
non-binary junction X_c is the pinch of both sections of c (B_c = C_c = 0). Both are the HS "knots" of the Underwood
model; the design stays at positive distance from them for fixed (δ, η).

**8.4 No sign information on traces.** The trace contents and trace net flows of v_n are O(ρ_∞^{n+1}); positivity of all
contents comes from Proposition 2.1 (b). Trace net flows may have either sign — for N ≥ 4 Agent C observed trace net flows
against the "natural" direction, changing sign with n (§9) — and no step of the proof uses their signs.

**8.5 Uniqueness, non-physical zeros, and why (H6) is needed.** v_n is unique near 𝓛. The exact system has other zeros
far from 𝓛, some with negative stage mole fractions (Agent C, §9): the hypothesis S > 0 in Proposition 2.1 (b) cannot be
dropped, and is delivered near 𝓛 by Lemma 6.4, which uses the nonnegativity of the limit half-profiles (H6) and the
positivity of the pivots. These are not consequences of the algebraic hypotheses: in Agent C's red-team case
`n3_invalid` (N = 3, δ = 0.002 with a prefractionator slack 0.05 > η̄(δ), so that the margin of [2, 3] is negative),
G_∞(v_∞) = 0 and DG_∞(v_∞) is invertible with the full block structure of Theorem 5.9, but the non-pinning section
B[2, 3] has a negative pivot and a negative half-profile entry, and the zeros of G_n converging to v_∞ (n = 32, 64, 128)
have negative stage mole fractions. The parameter region of Theorem 4.11 is therefore genuinely needed.

**8.6 Relation with earlier partial results.** At δ = η = 0 the specification is HS's preferred-split point of paper.tex
Theorem `achieve` (every non-binary column at its preferred split, binary columns with excess e_p). Astra's
prefractionator theorem is the single-column case of Proposition 4.3 and Lemma 6.3 with prescribed returns; network.md
§11.2's forced impurities for N = 3 are the pattern (B, T) of Definition 4.7 (with η_{[1,3]} = δ/2, h = 0).

---

## 9. Numerical evidence (not part of the proof)

All runs in high precision (mpmath, 50–230 digits), WSL Fedora, scripts and logs in the agents' folders.
- *Design* (`design/scripts/check_design.py`, `design/runs/check_design_12_50_1.log`): 240 random instances, N = 3..7
  (generic, nearly equal volatilities, tiny feed fractions, q = −1 and 2, two equal highest peaks), random patterns, two
  δ, two forcing levels: every family of (H1)–(H8) holds (ALL PASS); κ_p(δ) − κ_p(0) = δ and κ_p(0; δ) = e_p + δ. The limit
  configuration for network.md §11.3's flows reproduces N's finite-stage limits (junction ratios 0.696154, 0.553931,
  side impurity 0.01768) to all printed digits (`design/scripts/crosscheck_N3.py`).
- *Limit system* (`limit/scripts/check_limit_jacobian.py`, 16 random instances N = 3..6 with random patterns): G_∞(v_∞)
  ≤ 3e-57, the block structure of Theorem 5.9 exactly (zero blocks exactly 0 in finite differences), every trace row as
  predicted by Lemma 5.5 to 1e-31, A_k nonsingular M-matrices, σ_min(DG_∞) from 1.1e-4 to 9.1e-2. Unforced columns:
  σ_min(DG_∞) = 0 to 60 digits (`check_unforced.py`).
- *Finite stages, independent solver* (Agent C, `computation/REPORT.md`, `computation/FINAL-REPORT.md`, tables in
  `computation/runs/paper_tables.md`): 62 cases, N = 3..8 (all four patterns for N = 3, all eight for N = 4; q from −10
  to 10; α from (1.05, 1.02, 1) to (100, 10, 1); z down to 1e-4; equal highest peaks; δ down to 0.001; forcing up to
  0.95 f_p; n up to 768, 4000 for nearly equal volatilities; precision up to about 1500 digits; every steady state
  re-checked stage by stage), D's 12 exports recomputed independently (agreement ≤ 1.2e-54), L's code cross-validated
  (|G_n| ≤ 1e-95 at C's steady states): in all 61 cases inside the admissible region DG_∞ is invertible with the block
  structure of Theorem 5.9, the finite-stage zeros are physical for n ≥ n_0(observed), |v_n − v_∞| decays like ρ_∞^n;
  N = 3 to n = 256 with σ_min(DG_n) → σ_min(DG_∞) = 0.00729 and product impurities 5.6e-49 / 0.0177 (designed) / 4.3e-54;
  unequal section lengths (one section at 24 stages, the others at 192; N = 4 random lengths in [24, 96]) all physical;
  without forcing σ_min(DG_n) decays geometrically (8.2e-5, 9.8e-7, 1.4e-8 at n = 16, 24, 32) while with forcing it stays
  near 0.007; σ_min(DG_∞) ∝ η in the plain normalization; non-physical zeros of G_n exist far from 𝓛.
- *Analysis* (`analysis/scripts/check_analysis.py`, `analysis/runs/check_analysis.log`): with generic numerical
  eigen-data (not Underwood's formulas) and D's design: (A1)–(A6) hold; Newton zeros of G_n (3.3) for n = 8..64 and unequal
  lengths, N = 3 (q = 1, pattern BB; q = −0.3, pattern TB) and N = 4 (patterns BBB, BTB); every stage, junction, side-draw,
  condenser and reboiler balance re-checked with the actual streams (≤ 1e-66); S > 0 and all mole fractions > 0;
  |v_n − v_∞|/ρ_∞^{n+1} between 0.16 and 5.5 (bounded, as Theorem 6.2 predicts).

---

## 10. Status and dependencies

**Proved (complete proofs above, not yet independently verified):** Theorem 1 for every N ≥ 3 (N = 2: sections.md
Proposition 10) and Corollary 2. With version 2's lower bound, **the Halvorsen–Skogestad conjecture holds in the class 𝒜
for every number of components**, in the finite-stage-limit sense of class.md §4: R*(𝒮) ≥ max_j H_j = R*(𝒮_Petlyuk) for
every family 𝒮 of arrangements of 𝒜.

**Results used from earlier work (cited by file and statement):** the stage-level lower bound (class.md Theorem 5.3,
Corollary 5.5; early independent check CONFIRMED, provisional) — for Corollary 2 only; the M-matrix structure of the
component balances at fixed K-values (network.md Theorem 10.1 (i); its short proof is recalled in Proposition 2.1).
Everything else is re-proved here: HS's preferred split and carry-over (paper.tex Theorem `achieve`, version 1,
independently verified: Lemmas 4.1, 4.9, and the side-draw balances at the HS point through the cut identity of
Proposition 4.10 (c)), Astra's Lemmas 1–2 (Proposition 4.3), S's Lemma 11 (Lemma 3.5 (b)), N's projective form (§2.5).

**Not claimed:** uniformity of n_0 in (δ, η); explicit values of η_1, C_1 (compactness); any statement about steady
states far from 𝓛; non-ideal mixtures, several pressures, several feeds (outside HS's setting).
