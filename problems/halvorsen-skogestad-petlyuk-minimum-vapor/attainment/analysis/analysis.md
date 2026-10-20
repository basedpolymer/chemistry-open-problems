# Finite-stage attainment: the analytic step (Agent A)

Claude Opus 5.5 (effort xhigh), Agent A (analysis and integration). **Version 1.2, 2026-10-20 18:48 (v1.1 18:20: §8 against hypotheses.md v1; v1.2: §9 numerical check, pointer to theorem.md; v1.3, 19:45, independent verifier V2: one wording correction marked [V2] in the proof of Theorem A.8, see `../verification-A/verification.md` §4).**
Paths are relative to the problem folder `problems/halvorsen-skogestad-petlyuk-minimum-vapor/`. References:
[class] `total/foundations/class.md`; [S] `total/sections/sections.md`; [N] `total/network/network.md`;
[Astra] `total/attainment-astra/attainment.md`; [D] `attainment/design/`; [L] `attainment/limit/`.

## 0. Status

This file proves the **analytic half** of the attainment theorem: *if* an operating specification of the
generalized Petlyuk arrangement and an infinite-stage limit point satisfy the hypotheses (A1)–(A6) of §3 (design
properties, to be supplied by Agent D, and the invertibility of the limit Jacobian, Theorem L of Agent L), *then*
for all section lengths n_σ ≥ n_0 the arrangement has a steady state with all mole fractions positive, within
C ρ^{n+1} (n = min n_σ) of the limit point, with the product purities of the limit point up to C ρ^{n+1}, and with
the heat duty R of the specification. Every statement below is proved; nothing is numerical. The hypotheses are
stated in §3 and are the exact interface with D and L (§8 maps them onto D's (H1)–(H8)).

Main results: Lemma A.1 (a zero of the network map with S > 0 on every stage is a physical steady state, all mole
fractions > 0), Lemma A.2 (analytic eigen-data through zero trace flows; exact form of the trace left eigenvector),
Lemma A.3 (product form of the section relation; exponential C¹ bounds), Lemma A.4 (quantitative inverse function
step), Theorem A.5 (network perturbation), Lemmas A.6–A.7 (stage profiles and S > 0), Theorem A.8 (physical steady
states, purities, R), Corollary A.9 (attainment with the order of limits n → ∞, η → 0, δ → 0).

Three points differ from the coordinator's plan (BRIEF-common §3, Step 3); none is a flaw of the plan, but each
simplifies or sharpens it:
1. No chart coefficients are needed as unknowns: the section relation is written in **product form**
   (Lemma A.3) directly in the junction compositions and net flows; S's Lemma 12 chart reappears only in the
   reconstruction of the stage profiles (Lemma A.6).
2. The distance to the limit point is O(ρ^{n+1}), not O(n ρ^n): the C¹ bound O(n ρ^n) is needed only to make the
   Newton map a contraction, while the displacement is controlled by the residual G_n(v_∞) = O(ρ^{n+1}).
3. The **section lengths are arbitrary**: every section may have its own n_σ, only min_σ n_σ ≥ n_0 is required. The
   limit system does not involve the lengths because, in the forced design, no section is doubly pinched.

---

## 1. The arrangement, the specification, the unknowns

**Data.** N ≥ 3 components, α_1 > … > α_N > 0, feed of size 1 with composition z (z_p > 0, Σ z_p = 1) and liquid
fraction q ∈ ℝ. S(x) = Σ_i α_i x_i = αᵀx, y(x) = diag(α) x / S(x), Δ = {x ∈ ℝ^N : x ≥ 0, 1ᵀx = 1} (1 = all-ones
vector), G_x(φ) = Σ_i α_i x_i/(α_i − φ).

**The generalized Petlyuk arrangement** (HS 2003 Part 3, Fig. 1; version 1's FTC; BRIEF-common §2). Columns
c = [i, m], 1 ≤ i < m ≤ N (the set 𝒞, |𝒞| = N(N−1)/2); units J[i, m], 1 ≤ i ≤ m ≤ N. Column [i, m] has a top section
T[i, m] with bottom end J[i, m] and top end J[i, m−1], and a bottom section B[i, m] with top end J[i, m] and bottom
end J[i+1, m]. J[i, m] (i < m) is an equilibrium (E) stage without duty, the *junction* of [i, m]; the process feed
enters J[1, N]. J[1, 1] is a total condenser (reflux and liquid product 1), J[N, N] a partial reboiler (an (E) stage,
liquid product N, the only heat input), J[p, p] (2 ≤ p ≤ N−1) an (E) side-draw stage without duty (liquid product p)
between B[p−1, p] above and T[p, p+1] below. Every section σ consists of n_σ ≥ 1 equilibrium stages with constant
liquid and vapor flows L_σ, V_σ; it receives vapor from its bottom-end unit and liquid from its top-end unit, and
returns liquid to the former and vapor to the latter (thermal couplings). We write bot(σ), top(σ) for its end units.

**Specification.** An operating specification F prescribes L_σ, V_σ for every section and the product flows P_p.
(A1) below lists what is required of it (this is Agent D's design). Put Δ_σ := V_σ − L_σ.

**Unknowns.** For each column c = [i, m] its *top split* d^c ∈ ℝ^N and its *junction composition* X^c ∈ ℝ^N, with
1ᵀd^c = Δ_{T[c]} and 1ᵀX^c = 1. The set of such v = (d^c, X^c)_{c∈𝒞} is an affine space 𝒱 of dimension
𝒩 := 2(N−1)|𝒞| = N(N−1)². Derived quantities, defined recursively from the feed column downward (labels shrink, so
the recursion is well founded):

  f^{[1,N]} := z,  f^{[i,m]} := d^{[i,m+1]} [m < N] + b^{[i−1,m]} [i > 1]  (the net feed of [i, m] from its parents),
  b^c := f^c − d^c,  w_{T[c]} := d^c,  w_{B[c]} := −b^c  (net upward component flows of the sections),
  x_p := ( d^{[p,p+1]} [p < N] + b^{[p−1,p]} [p > 1] ) / P_p  (product compositions).

The composition of a unit is X^c at J[c] (c ∈ 𝒞) and x_p at J[p, p]; x_bot(σ), x_top(σ) denote the compositions of
the end units of σ. All these maps v ↦ (w_σ, x_bot(σ), x_top(σ)) are affine. By the total balances of (A1),
1ᵀw_σ = Δ_σ for every σ and 1ᵀx_p = 1 (proof: 1ᵀb^c = 1ᵀf^c − Δ_{T[c]} = −Δ_{B[c]} by the total balance at J[c],
and 1ᵀ(d^{[p,p+1]} + b^{[p−1,p]}) = Δ_{T[p,p+1]} − Δ_{B[p−1,p]} = P_p by the total balance at J[p, p]).

**Section matrices.** For a section σ and w ∈ ℝ^N put M_σ(w) := V_σ diag(α) − w αᵀ. If 1ᵀw = Δ_σ then
1ᵀM_σ(w) = L_σ αᵀ. Its eigenvalues are V_σ φ with Ψ_{σ,w}(φ) := V_σ − Σ_i α_i w_i/(α_i − φ) = 0, together with
V_σ α_h for every h with w_h = 0; a left eigenvector for V_σ φ (φ ∉ {α_i}) is (α_i/(α_i − φ))_i, i.e. ℓ·x = G_x(φ)
([S] Lemma 11, [class] §6).

**Lemma A.1 (from a zero to a physical steady state).** Let v ∈ 𝒱 and suppose that for every section σ:
(a) x_top(σ) = c_σ M_σ(w_σ)^{n_σ+1} x_bot(σ) for some c_σ ∈ ℝ;
(b) the vectors u^σ_k := M_σ(w_σ)^k x_bot(σ), 0 ≤ k ≤ n_σ + 1, satisfy 1ᵀu^σ_k ≠ 0, and the compositions
x^σ_k := u^σ_k/(1ᵀu^σ_k) satisfy S(x^σ_k) > 0 for 0 ≤ k ≤ n_σ;
and suppose S > 0 at the composition of every unit. Assign to stage k of σ (k = 1..n_σ, counted upward) the liquid
x^σ_k and the vapor y(x^σ_k). Then this is a steady state of the arrangement with the flows F (in the sense of
[class] §3), x^σ_{n_σ+1} = x_top(σ), c_σ > 0, every mole fraction of every unit and stream is **strictly positive**,
and the total heat supplied is R(F) = V_{B[N−1,N]}.

*Proof.* Fix σ, write M = M_σ(w_σ), L = L_σ, V = V_σ, w = w_σ, n = n_σ. Since 1ᵀM = L αᵀ, 1ᵀu_{k+1} = L S(x_k) 1ᵀu_k,
so x_{k+1} = M x_k/(L S(x_k)) for 0 ≤ k ≤ n, i.e. L S(x_k) x_{k+1} = V diag(α) x_k − w S(x_k), i.e.

  V y(x_k) − L x_{k+1} = w   (0 ≤ k ≤ n)   (cross-section balances).

By (a) and 1ᵀx_top = 1, x_top = M^{n+1} x_bot/(1ᵀM^{n+1} x_bot) = x_{n+1}, and 1/c_σ = 1ᵀu_{n+1} = L^{n+1} Π_{k=0}^{n} S(x_k) > 0.
*Stages.* For 1 ≤ k ≤ n the component balance of stage k is (L x_{k+1} + V y(x_{k−1})) − (L x_k + V y(x_k)) =
(V y(x_{k−1}) − L x_k) − (V y(x_k) − L x_{k+1}) = w − w = 0. The vapor entering stage 1 is the vapor of the
bottom-end unit, of composition y(x_0) = y(x_bot); the liquid entering stage n is that of the top-end unit, of
composition x_{n+1} = x_top. (For σ = T[1, 2] the top-end unit is the total condenser; for σ = B[N−1, N] the
bottom-end unit is the reboiler, an (E) stage.)
*Junctions.* At a unit J the component balance reads, with the two end cross-section balances,
Σ_{σ: bot(σ)=J} (L_σ x^σ_1 − V_σ y(x_J)) + Σ_{σ: top(σ)=J} (V_σ y(x^σ_n) − L_σ x_J) + (feed) − (product)
= −Σ_{bot(σ)=J} w_σ + Σ_{top(σ)=J} w_σ + (feed) − (product). At J[c], c = [i, m], the sections with bottom end J[c]
are T[c] and B[i−1, m], those with top end J[c] are B[c] and T[i, m+1]; the balance is
−d^c + b^{[i−1,m]} − b^c + d^{[i,m+1]} (+ z if c = [1, N]) = f^c − d^c − b^c = 0 by definition of b^c. At J[p, p]
it is d^{[p,p+1]} + b^{[p−1,p]} − P_p x_p = 0 by definition of x_p; for the condenser (J[1,1], no vapor outflow)
V y(x^σ_n) − (L + P_1) x_1 = w + L x_1 − (L + P_1) x_1 = d^{[1,2]} − P_1 x_1 = 0 (and y(x^σ_n) = x_1 follows since
V = L + P_1). Total and CMO energy balances hold by (A1); hence the duties are those of (A1) and
R = Σ_u max(Q_u, 0) = V_{B[N−1,N]}.
*Positivity.* Set K_{u,i} := α_i/S(x_u) > 0 at every unit u with a vapor outflow (all units except the condenser,
by hypothesis). The steady-state equations just verified are exactly the equations F_{u,i} = 0 of [N] §10 with these
K, so for each component i the vector (x_{u,i})_u solves the linear system −A_i(K) X = (feed) z_i of [N] Theorem 10.1
(i). With all flows positive, every unit is reachable from J[1, N] and every unit reaches a product outlet along
positive flows (liquid flows down to the reboiler, which delivers product N; vapor and liquid are exchanged both
ways at every coupling), so −A_i(K) is a nonsingular M-matrix whose inverse is entrywise positive between units
joined by a path ([N] Thm 10.1 (i)); hence x_{u,i} > 0 for all u and i. Stage compositions inside the sections are
units as well. Vapor compositions y = K x are then positive too. ∎

*Remark.* Positivity needs **no** information on the signs of the trace net flows or of the trace contents: only
S > 0 at every unit (the observation recorded in BRIEF-common §1, "positivity tool").

---

## 2. Sections: analytic eigen-data and the product form

Fix a section σ (flows L, V > 0) and w^∞ ∈ ℝ^N with 1ᵀw^∞ = Δ_σ.

**Lemma A.2 (eigen-data).** Suppose M(w^∞) has N distinct positive eigenvalues λ^∞_1 > … > λ^∞_N > 0. Then there is
an open neighborhood W of w^∞ in ℝ^N such that:
(a) the eigenvalues λ_1(w) > … > λ_N(w) > 0 of M(w) are simple and real-analytic on W, and there are real-analytic
right and left eigenvectors e_r(w), ℓ_r(w) with ℓ_r·e_s = δ_rs; any other analytic choice differs by
e_r ↦ h_r e_r, ℓ_r ↦ ℓ_r/h_r with h_r analytic and nonvanishing;
(b) if 1ᵀe_{r*}(w^∞) ≠ 0, one may (and we do) normalize 1ᵀe_{r*} = 1 on W (after shrinking W);
(c) (trace modes) if w^∞_k = 0 and λ^∞_ϱ = V α_k for some index ϱ = ϱ(k), then Ψ^{(k)}_{w^∞}(α_k) ≠ 0, where
Ψ^{(k)}_w(φ) := V − Σ_{i≠k} α_i w_i/(α_i − φ); the left eigenvector of the branch λ_ϱ may be normalized by
(ℓ_ϱ)_k = 1, and then, with φ(w) := λ_ϱ(w)/V,

  ℓ_ϱ(w)·x = x_k + w_k · G^{(k)}_x(φ(w)) / Ψ^{(k)}_w(φ(w)),   G^{(k)}_x(φ) := Σ_{i≠k} α_i x_i/(α_i − φ),   (2.1)

exactly, for all w ∈ W and x ∈ ℝ^N; in particular ℓ_ϱ(w) = e_k whenever w_k = 0, and the derivative of
(w, x) ↦ ℓ_ϱ(w)·x at a point (w, x) with w_k = 0 and x_k = 0 is (δw, δx) ↦ δx_k + δw_k G^{(k)}_x(α_k)/Ψ^{(k)}_w(α_k);
(d) τ_{rs}(w) := (λ_s(w)/λ_r(w))^{n+1} (r < s) satisfies 0 < τ_{rs} ≤ ρ_W^{n+1} and
|∇_w τ_{rs}| ≤ (n+1) C_W ρ_W^{n+1} on W, where ρ_W := sup_W max_r λ_{r+1}/λ_r < 1 and
C_W := sup_W max_{r<s} |∇_w log(λ_s/λ_r)| < ∞ (after shrinking W to a compact neighborhood).

*Proof.* (a) The characteristic polynomial p(λ, w) = det(λ − M(w)) is polynomial in (λ, w) and ∂_λ p ≠ 0 at a
simple root, so each simple eigenvalue is analytic (implicit function theorem) and stays real, positive and simple
near w^∞; the spectral projector P_r(w) = Π_{s≠r} (M(w) − λ_s(w))/(λ_r(w) − λ_s(w)) is analytic of rank one. Take
e_r := P_r(w) g_r / (h_rᵀ P_r(w) g_r) and ℓ_r := P_r(w)ᵀ h_r for fixed vectors g_r, h_r with h_rᵀ P_r(w^∞) g_r ≠ 0.
Uniqueness up to scalars: the eigenspaces are one-dimensional. (b) Rescale e_r by 1/(1ᵀe_r) and ℓ_r by 1ᵀe_r.
(c) Up to the factor V^{N−1}, the characteristic polynomial in φ (λ = Vφ) is p̃(φ) := Π_i (φ − α_i)·Ψ_w(φ), a
polynomial (the poles cancel): p̃(φ) = V Π_i(φ − α_i) + Σ_i α_i w_i Π_{j≠i}(φ − α_j). At φ = α_k only the i = k term of
the sum survives, so p̃(α_k) = α_k w_k Π_{j≠k}(α_k − α_j) vanishes iff w_k = 0. When w_k = 0,
p̃(φ) = (φ − α_k) Π_{i≠k}(φ − α_i) Ψ^{(k)}_w(φ), so α_k is a simple root iff Ψ^{(k)}_w(α_k) ≠ 0; this proves
Ψ^{(k)}_{w^∞}(α_k) ≠ 0. For w ∈ W with w_k ≠ 0, λ_ϱ(w) = Vφ with φ ∉ {α_i} (φ ≠ α_k because w_k ≠ 0; φ ≠ α_i for
i ≠ k because φ is close to α_k), and ℓᵀM = λℓᵀ reads ℓ_i (Vα_i − λ) = (ℓ·w) α_i, so
ℓ_i = (ℓ·w) α_i / (V(α_i − φ)) and ℓ·w ≠ 0 (otherwise ℓ would be a multiple of a unit vector e_j with λ = Vα_j,
impossible). Normalizing ℓ_k = 1 gives ℓ_i = α_i (α_k − φ)/(α_k (α_i − φ)). The root equation Ψ_w(φ) = 0 reads
Ψ^{(k)}_w(φ) = α_k w_k/(α_k − φ), hence α_k − φ = α_k w_k/Ψ^{(k)}_w(φ) and ℓ_i = α_i w_k / ((α_i − φ) Ψ^{(k)}_w(φ)) for
i ≠ k, which is (2.1). The right side of (2.1) is analytic on W (shrink W so that Ψ^{(k)}_w(φ(w)) ≠ 0 and
φ(w) ≠ α_i for i ≠ k); both sides are analytic and agree on the dense set {w_k ≠ 0}, hence everywhere on W. The
derivative statement follows from (2.1), since φ(w) = α_k when w_k = 0.
(d) τ_{rs} = exp((n+1) log(λ_s/λ_r)) and ∇τ_{rs} = (n+1) τ_{rs} ∇ log(λ_s/λ_r). ∎

*Reading of (2.1).* A closed-end trace condition ℓ_ϱ·x_J = 0 at an end unit J of σ says
x_{J,k} = −w_k G^{(k)}_{x_J}(φ)/Ψ^{(k)}(φ): the trace content of J is proportional to the trace net flow, with the
factor κ_σ = |G^{(k)}_{x_J}(α_k)|/Ψ^{(k)}_σ(α_k) of BRIEF-common §3, and it has the sign of the net flow *from J into
σ* whenever Ψ^{(k)}_σ(α_k) > 0 (G^{(k)}_{x_J}(α_k) > 0 at a top end where k is heavier than every component present,
< 0 at a bottom end where k is lighter). These sign facts are used by L (trace block); the analysis below does not
use them.

**Pinch index, pivots, product form.** Let r* ∈ {1..N} be a designated index (the *pinch*), U := {r < r*} (modes
above the pinch), D := {s > r*} (modes below). For (x_bot, x_top, w) ∈ ℝ^N × ℝ^N × W put
π_bot := ℓ_{r*}(w)·x_bot, π_top := ℓ_{r*}(w)·x_top, μ_r := λ_{r*}/λ_r (r ∈ U), μ_s := λ_s/λ_{r*} (s ∈ D), and

  E^{(n)}_r := (ℓ_r·x_bot) π_top − μ_r^{n+1} (ℓ_r·x_top) π_bot   (r ∈ U),
  E^{(n)}_s := (ℓ_s·x_top) π_bot − μ_s^{n+1} (ℓ_s·x_bot) π_top   (s ∈ D),                (2.2)

and E^{(∞)} the same expressions with μ^{n+1} replaced by 0.

**Lemma A.3 (product form).** (a) If π_bot π_top ≠ 0, then E^{(n)}(x_bot, x_top, w) = 0 (N − 1 equations) holds
iff x_top = c M(w)^{n+1} x_bot for some c ≠ 0, and then c = π_top/(λ_{r*}^{n+1} π_bot).
(b) E^{(∞)} = 0 with π_bot π_top ≠ 0 iff ℓ_r·x_bot = 0 for all r ∈ U and ℓ_s·x_top = 0 for all s ∈ D: the
**limit relation** of [S] Lemma 11 with pinch r*.
(c) On any set where |x_bot|, |x_top| are bounded and w ∈ W: |E^{(n)} − E^{(∞)}| ≤ C ρ_W^{n+1} and
‖D(E^{(n)} − E^{(∞)})‖ ≤ C (n + 2) ρ_W^{n+1}, with C depending on the bounds and on sup_W of the eigen-data and of
their first derivatives. At a fixed point (x_bot, x_top, w^∞): |E^{(n)} − E^{(∞)}| ≤ C ρ_∞^{n+1}, where
ρ_∞ := max(λ^∞_{r*+1}/λ^∞_{r*}, λ^∞_{r*}/λ^∞_{r*−1}) (terms absent if r* = N or r* = 1).
(d) Replacing the ℓ_r by h_r ℓ_r (h_r analytic, nonvanishing) multiplies each component of E^{(n)} and E^{(∞)} by a
nonvanishing analytic function; at a zero of E^{(∞)} the derivative is multiplied on the left by an invertible
diagonal matrix. Likewise DE^{(∞)} = diag(π_top (r ∈ U), π_bot (s ∈ D))·DΛ at a zero, where
Λ := ((ℓ_r·x_bot)_{r∈U}, (ℓ_s·x_top)_{s∈D}).

*Proof.* (a) Expand in the eigenbasis: x_top = c M^{n+1} x_bot iff ℓ_r·x_top = c λ_r^{n+1} ℓ_r·x_bot for every r. If
this holds with c ≠ 0, then π_top = c λ_{r*}^{n+1} π_bot fixes c; for r ∈ U, multiplying ℓ_r·x_top = c λ_r^{n+1} ℓ_r·x_bot
by π_bot λ_{r*}^{n+1}/λ_r^{n+1} gives μ_r^{n+1}(ℓ_r·x_top)π_bot = (ℓ_r·x_bot)π_top, i.e. E_r = 0; for s ∈ D,
multiplying ℓ_s·x_top = c λ_s^{n+1} ℓ_s·x_bot by π_bot gives E_s = 0. Conversely, define c := π_top/(λ_{r*}^{n+1}π_bot)
≠ 0; the r*-relation holds by definition, and E_r = 0 (r ∈ U) gives
ℓ_r·x_top = (λ_r/λ_{r*})^{n+1} (π_top/π_bot) ℓ_r·x_bot = c λ_r^{n+1} ℓ_r·x_bot; E_s = 0 (s ∈ D) likewise.
(b) Immediate. (c) E^{(n)} − E^{(∞)} = −Σ μ^{n+1}·(products of analytic factors); use Lemma A.2 (d) (μ_r and μ_s
are among the ratios λ_s/λ_r, r < s) and the product rule. At w^∞ the largest μ is the ratio of λ^∞_{r*} to an
adjacent eigenvalue. (d) The factors h_r cancel in μ and appear as h_r h_{r*} in E_r; at a zero of a product f g with
f = 0 the derivative is g Df. ∎

*Remark (relation with [S] Lemma 12).* Writing a_s := ℓ_s·x_bot/π_bot (s ∈ D) and b_r := ℓ_r·x_top/π_top
(r ∈ U), (a) says x_bot ∝ e_{r*} + Σ_D a_s e_s + Σ_U μ_r^{n+1} b_r e_r and x_top ∝ e_{r*} + Σ_U b_r e_r +
Σ_D μ_s^{n+1} a_s e_s: this is S's chart (u_X, u_Y) with t = μ^{n+1}. We do not use (a, b) as unknowns; they are
functions of (x_bot, x_top, w), bounded as long as the pivots stay away from 0.

---

## 3. The network map, its limit, and the hypotheses

For section lengths n = (n_σ)_σ define G_n on a neighborhood of v_∞ in 𝒱 by collecting, over all 2|𝒞| = N(N−1)
sections,

  G_n(v) := ( E^{(n_σ)}_σ( x_bot(σ)(v), x_top(σ)(v), w_σ(v) ) )_σ ,   G_∞ := same with E^{(∞)}_σ,

each block having N − 1 components, so G_n maps the 𝒩-dimensional affine space 𝒱 to ℝ^𝒩. The pinch index r*_σ of
each section is part of the data supplied by the design. By Lemma A.3 (a), a zero of G_n at which all pivots are
nonzero satisfies hypothesis (a) of Lemma A.1.

**Hypotheses** (to be supplied by D, (A1)–(A5), and by L, (A6); §8 gives the correspondence with D's (H1)–(H8)).

- **(A1) Specification.** All L_σ, V_σ, P_p are > 0; at every unit the total and the CMO (energy) balances hold;
  the duty is zero at every junction J[i, m] (i < m) and every side-draw stage J[p, p] (2 ≤ p ≤ N−1); J[1, 1] is a
  total condenser and J[N, N] a partial reboiler. Then R(F) = V_{B[N−1,N]}.
- **(A2) Limit point and spectra.** A point v_∞ ∈ 𝒱 (same flows F) such that, for every section σ, M_σ(w^∞_σ) has N
  distinct positive eigenvalues, a pinch index r*_σ is designated, and every unit composition of v_∞ (all X^c and
  x_p) lies in Δ.
- **(A3) Pivots.** For every σ, π^∞_bot := ℓ^∞_{r*}·x^∞_bot ≠ 0 and π^∞_top := ℓ^∞_{r*}·x^∞_top ≠ 0.
- **(A4) Limit relations.** G_∞(v_∞) = 0, i.e. (Lemma A.3 (b)) ℓ^∞_r·x^∞_bot = 0 for r < r*_σ and
  ℓ^∞_s·x^∞_top = 0 for s > r*_σ, for every σ.
- **(A5) Half-profiles.** For every σ, 1ᵀe^∞_{r*} ≠ 0, and with the normalization 1ᵀe_{r*} = 1 = ℓ_{r*}·e_{r*} the
  limit half-profiles

    ξ^bot_k := e^∞_{r*} + Σ_{s∈D} (ℓ^∞_s·x^∞_bot/π^∞_bot) (λ^∞_s/λ^∞_{r*})^k e^∞_s ,
    ξ^top_k := e^∞_{r*} + Σ_{r∈U} (ℓ^∞_r·x^∞_top/π^∞_top) (λ^∞_{r*}/λ^∞_r)^k e^∞_r ,   k = 0, 1, 2, …

  are componentwise ≥ 0. (By (A4), ξ^bot_k = M_∞^k x^∞_bot/(π^∞_bot (λ^∞_{r*})^k) is the forward orbit of the limit
  bottom-end composition and ξ^top_k = (λ^∞_{r*})^k M_∞^{−k} x^∞_top/π^∞_top the backward orbit of the top-end one,
  M_∞ := M_σ(w^∞_σ). (A5) then forces π^∞_bot > 0, π^∞_top > 0 (k = 0) and e^∞_{r*} ∈ Δ (k → ∞).)
- **(A6) Theorem L.** DG_∞(v_∞) : T𝒱 → ℝ^𝒩 is invertible. By Lemma A.3 (d) this is equivalent to the
  invertibility of the derivative of the limit relations Λ(v) = ((ℓ_r·x_bot)_{r<r*}, (ℓ_s·x_top)_{s>r*})_σ for any
  analytic normalization of the left eigenvectors (in particular the Underwood normalization ℓ·x = G_x(φ) for
  modes with φ away from the α_i, and the trace normalization (2.1)).

Useful sufficient conditions for (A5) (this is how D proves it): if w^∞_σ ≥ 0 (a top section with nonnegative
limit net flows) and L_σ > 0, then M_∞^{−1} = V^{−1}diag(α)^{−1} + (V L)^{−1} diag(α)^{−1} w 1ᵀ ≥ 0 ([Astra] (19)), so
ξ^top_k ≥ 0 for all k as soon as x^∞_top ∈ Δ and π^∞_top > 0; if w^∞_σ ≤ 0 (a bottom section), M_∞ ≥ 0 and
ξ^bot_k ≥ 0 as soon as x^∞_bot ∈ Δ and π^∞_bot > 0. On the other side a **cone decomposition** suffices: if
x^∞_bot = A e^∞_{r*} + Σ_D A_s e^∞_s with A > 0, A_s ≥ 0 and e^∞_s ≥ 0 for the s with A_s ≠ 0, then ξ^bot_k ≥ 0
for all k ([Astra] Lemma 2 and (16)); symmetrically at the top. (When the forward orbit is nonnegative,
π^∞_bot > 0 is automatic: M_∞^k x_bot/(λ^∞_{r*})^k → π^∞_bot e^∞_{r*} with e^∞_{r*} ∈ Δ.)

---

## 4. The perturbation theorem

**Lemma A.4 (quantitative inverse function step).** Let E be a finite-dimensional normed space, B̄ := B̄(v_∞, r_0) ⊂ E,
F_∞, F ∈ C¹(B̄; E), F_∞(v_∞) = 0, A := DF_∞(v_∞) invertible. Let ω(r) := sup_{|v−v_∞|≤r} ‖DF_∞(v) − A‖, choose
r_1 ∈ (0, r_0] with ‖A^{−1}‖ ω(r_1) ≤ 1/4, and suppose sup_{B̄} ‖DF − DF_∞‖ ≤ 1/(4‖A^{−1}‖) and
‖A^{−1}‖ |F(v_∞)| ≤ r_1/2. Then F has exactly one zero v* in B̄(v_∞, r_1), and |v* − v_∞| ≤ 2 ‖A^{−1}‖ |F(v_∞)|.

*Proof.* Φ(v) := v − A^{−1}F(v) has DΦ = A^{−1}(A − DF(v)) = A^{−1}(A − DF_∞(v)) + A^{−1}(DF_∞(v) − DF(v)), so
‖DΦ‖ ≤ 1/4 + 1/4 = 1/2 on the convex set B̄(v_∞, r_1) and Φ is 1/2-Lipschitz there. For |v − v_∞| ≤ r_1,
|Φ(v) − v_∞| ≤ |Φ(v) − Φ(v_∞)| + |A^{−1}F(v_∞)| ≤ r_1/2 + r_1/2. By Banach's fixed point theorem Φ has exactly one
fixed point v* in the ball; fixed points of Φ are the zeros of F. Finally
|v* − v_∞| ≤ |Φ(v*) − Φ(v_∞)| + |A^{−1}F(v_∞)| ≤ |v* − v_∞|/2 + ‖A^{−1}‖ |F(v_∞)|. ∎

**Theorem A.5 (network perturbation).** Assume (A1)–(A4) and (A6). There are r_1 > 0, n_0 and C such that for all
section lengths (n_σ) with n := min_σ n_σ ≥ n_0, the map G_n has exactly one zero v_n in B̄(v_∞, r_1), and

  |v_n − v_∞| ≤ C ρ_∞^{n+1},   ρ_∞ := max_σ max( λ^∞_{σ,r*+1}/λ^∞_{σ,r*}, λ^∞_{σ,r*}/λ^∞_{σ,r*−1} ) < 1.

At v_n all pivots are nonzero, so for every σ, x_top(σ) = c_σ M_σ(w_σ)^{n_σ+1} x_bot(σ) with c_σ ≠ 0.

*Proof.* By (A2) and Lemma A.2 each section has a compact neighborhood W_σ of w^∞_σ with analytic eigen-data and
ρ_σ := sup_{W_σ} max λ_{r+1}/λ_r < 1; put ρ_1 := max_σ ρ_σ < 1. By (A3) and continuity there is r_0 > 0 such that for
v ∈ B̄(v_∞, r_0): w_σ(v) ∈ W_σ and |π_bot|, |π_top| ≥ π_0/2 for every σ, where π_0 := min_σ min(|π^∞_bot|, |π^∞_top|).
G_∞ and G_n are analytic on B̄(v_∞, r_0) (compositions of affine maps with Lemma A.2's eigen-data). By Lemma A.3 (c),
applied section by section with n_σ in place of n,

  sup_{B̄(v_∞,r_0)} ‖DG_n − DG_∞‖ ≤ C_1 max_σ (n_σ + 2) ρ_1^{n_σ+1} ≤ C_1 (n + 2) ρ_1^{n+1}
                    (for n ≥ 1/log(1/ρ_1), where x ↦ (x+2)ρ_1^{x+1} is decreasing),
  |G_n(v_∞)| = |(G_n − G_∞)(v_∞)| ≤ C_2 ρ_∞^{n+1}   ((A4) and the second part of Lemma A.3 (c)).

By (A6), A := DG_∞(v_∞) is invertible; G_∞ is C² so ω(r) ≤ C_3 r, and r_1 := min(r_0, 1/(4 C_3 ‖A^{−1}‖)) works in
Lemma A.4. Choose n_0 with C_1 (n_0 + 2) ρ_1^{n_0+1} ≤ 1/(4‖A^{−1}‖) and ‖A^{−1}‖ C_2 ρ_∞^{n_0+1} ≤ r_1/2. Lemma A.4
with F = G_n gives the unique zero and |v_n − v_∞| ≤ 2 ‖A^{−1}‖ C_2 ρ_∞^{n+1}. The last statement is Lemma A.3 (a). ∎

*Constants.* r_1, n_0, C depend on the data (α, z, q) and on the specification and limit point, i.e. on (δ, η)
in D's family; they are not uniform as η → 0 (the pivots of the non-binary columns are proportional to η by
[Astra] (18), B, C ~ η, and ‖A^{−1}‖ grows), which is why n → ∞ is taken first (§6).

---

## 5. Stage profiles and positivity

Let v = v_n be the zero of Theorem A.5 and σ a section (we drop σ). With the eigen-data at w = w_σ(v_n), the
pivots π_bot, π_top ≠ 0, a_s := ℓ_s·x_bot/π_bot (s ∈ D) and b_r := ℓ_r·x_top/π_top (r ∈ U):

**Lemma A.6 (exact composite profile).** For 0 ≤ k ≤ n_σ + 1, u_k := M^k x_bot/(π_bot λ_{r*}^k) equals

  u_k = e_{r*} + Σ_{s∈D} a_s μ_s^k e_s + Σ_{r∈U} b_r μ_r^{n_σ+1−k} e_r ,                (5.1)

and the stage composition of Lemma A.1 is x_k = u_k/(1ᵀu_k) (when 1ᵀu_k ≠ 0).

*Proof.* M^k x_bot = Σ_r (ℓ_r·x_bot) λ_r^k e_r. For r ∈ U, E_r = 0 gives ℓ_r·x_bot = μ_r^{n+1} b_r π_bot, so the
r-term of u_k is b_r μ_r^{n+1} (λ_r/λ_{r*})^k e_r = b_r μ_r^{n+1−k} e_r. The others are immediate. ∎

So the profile is the superposition of a **bottom half-profile** (decaying from the bottom end at the rates μ_s)
and a **top half-profile** (decaying from the top end), glued at the pinch e_{r*}. At v_∞ these are exactly the
ξ^bot_k and ξ^top_k of (A5).

**Lemma A.7 (S > 0 at every stage).** Assume (A1)–(A6). Let c_0 := min_σ min(inf_k 1ᵀξ^bot_k, inf_k 1ᵀξ^top_k).
Then c_0 > 0, and there are n_1 ≥ n_0 and C such that for n ≥ n_1, every σ and every 0 ≤ k ≤ n_σ + 1:

  |u_k − ξ^bot_k|_∞ ≤ C (ρ_∞^{n+1} + ρ_1^{n_σ+1−k}),   |u_k − ξ^top_{n_σ+1−k}|_∞ ≤ C (ρ_∞^{n+1} + ρ_1^k),
  1ᵀu_k ≥ c_0/2,   S(x_k) ≥ α_N c_0/(2·1ᵀu_k) > 0 :

the stage composition at k is exponentially close to the bottom half-profile far from the top and to the top
half-profile far from the bottom, and S > 0 at every stage.

*Proof.* c_0 > 0: ξ^bot_k ≥ 0 by (A5), ξ^bot_k ≠ 0 (M_∞ is invertible and x^∞_bot ≠ 0), so 1ᵀξ^bot_k > 0, and
1ᵀξ^bot_k → 1ᵀe^∞_{r*} = 1 as k → ∞; a positive sequence with a positive limit has a positive infimum. Same for ξ^top.
Closeness: the quantities e_r, μ_r, a_s, b_r in (5.1) are analytic functions of v on B̄(v_∞, r_1), equal at v_∞ to
the coefficients of (A5), so they differ from them by ≤ C |v_n − v_∞| ≤ C ρ_∞^{n+1}; since all μ, μ^∞ ≤ ρ_1,
|μ^k − (μ^∞)^k| ≤ k ρ_1^{k−1} |μ − μ^∞| ≤ (sup_j j ρ_1^{j−1}) |μ − μ^∞|, uniformly in k. Comparing (5.1) termwise with
ξ^bot_k: the e_{r*}- and D-terms differ by ≤ C ρ_∞^{n+1}, and the U-terms are ≤ C ρ_1^{n_σ+1−k} (b_r is bounded,
since |π_top| ≥ π_0/2). Comparing with ξ^top_{n_σ+1−k}: the U-terms differ by ≤ C ρ_∞^{n+1} and the D-terms are
≤ C ρ_1^k. Since min(ρ_1^k, ρ_1^{n_σ+1−k}) ≤ ρ_1^{(n+1)/2}, for every k one of the two half-profiles, call it ξ
(≥ 0, 1ᵀξ ≥ c_0), satisfies |u_k − ξ|_∞ ≤ ε_n := C(ρ_∞^{n+1} + ρ_1^{(n+1)/2}). Then 1ᵀu_k ≥ c_0 − N ε_n and
αᵀu_k ≥ α_N 1ᵀξ − (Σ_i α_i) ε_n ≥ α_N c_0 − (Σ_i α_i) ε_n. For n ≥ n_1 with N ε_n ≤ c_0/2 and
(Σ_i α_i) ε_n ≤ α_N c_0/2 we get 1ᵀu_k ≥ c_0/2 and S(x_k) = αᵀu_k/1ᵀu_k ≥ α_N c_0/(2·1ᵀu_k) > 0. ∎

---

## 6. Conclusion: physical steady states, purities, heat duty, order of limits

**Theorem A.8 (finite-stage steady states near the limit point).** Assume (A1)–(A6). There are n_1 and C such that
for all section lengths with min_σ n_σ =: n ≥ n_1 the generalized Petlyuk arrangement operated with the flows F has
a steady state (in the sense of [class] §3) such that:
(i) every mole fraction of every unit and every stream is strictly positive;
(ii) every unit composition (junctions and products) and every column split is within C ρ_∞^{n+1} of its value at
v_∞; in particular the purity of product p is ≥ x^∞_{p,p} − C ρ_∞^{n+1};
(iii) the stage compositions of every section are within C(ρ_∞^{n+1} + ρ_1^{min(k, n_σ+1−k)}) of the limit
half-profiles (Lemma A.7);
(iv) the heat supplied is R = R(F) = V_{B[N−1,N]};
(v) it is the only zero of G_n in B̄(v_∞, r_1) (local uniqueness in the junction/net-flow variables).

*Proof.* Theorem A.5 gives v_n; by Lemma A.3 (a) hypothesis (a) of Lemma A.1 holds; by Lemma A.7 hypothesis (b)
holds and S > 0 at the stages; S > 0 at the junction and product compositions because they are affine in v and
within C ρ_∞^{n+1} of points of Δ, where S ≥ α_N (enlarge n_1 if necessary so that S ≥ α_N/2 there) [V2]. Lemma A.1 gives
(i) and (iv); (ii) is Theorem A.5 (the product compositions are affine in v); (iii) is Lemma A.7; (v) is
Theorem A.5. ∎

**Corollary A.9 (attainment, conditional on D and L).** Suppose that for every δ > 0 there is η_0(δ) > 0 and, for
every η ∈ (0, η_0(δ)], a specification F(δ, η) of the generalized Petlyuk arrangement with
R(F(δ, η)) ≤ max_j H_j + δ and a limit point v_∞(δ, η) satisfying (A1)–(A6), whose product compositions satisfy
x^∞_{p,p} ≥ 1 − c(δ) η for every p. Then for every δ > 0 and ε > 0 there is n_0 such that, for all section lengths
n_σ ≥ n_0, the arrangement has a steady state with every product of purity ≥ 1 − ε and R ≤ max_j H_j + δ. Hence
R*(𝒮_Petlyuk) ≤ max_j H_j, and with the lower bound ([class] Cor. 5.5) R*(𝒮_Petlyuk) = max_j H_j.

*Proof.* Given δ, ε, take η := min(η_0(δ), ε/(2 c(δ))), then n_0 := n_1(δ, η) of Theorem A.8, enlarged so that
C(δ, η) ρ_∞(δ, η)^{n_0+1} ≤ ε/2. For n_σ ≥ n_0 the steady state of Theorem A.8 has purities
≥ 1 − c(δ)η − ε/2 ≥ 1 − ε and R = R(F(δ, η)) ≤ max_j H_j + δ. So R*(𝒮_Petlyuk; ε) ≤ max_j H_j + δ for all δ > 0 and
ε > 0, i.e. R*(𝒮_Petlyuk) ≤ max_j H_j ([class] §4). The lower bound R*(𝒮) ≥ max_j H_j holds for every 𝒮 ⊆ 𝒜. ∎

**Order of limits.** n → ∞ at fixed (δ, η) (Theorem A.8: n_0 = n_0(δ, η), typically n_0 ≳ log(1/η)/log(1/ρ_∞) since
the non-binary pivots are O(η)); then η → 0 at fixed δ (limit purities 1 − O(η)); then δ → 0 (R → max_j H_j). No
uniformity in (δ, η) is claimed or needed: R*(𝒮; ε) is an infimum over all stage numbers.

---

## 7. Remarks

1. **Section lengths.** Theorem A.8 holds for *any* family of section lengths with min n_σ ≥ n_1: equal lengths,
   or lengths growing at different rates. This is possible because the limit relations (A4) and the Jacobian (A6) do
   not involve the lengths, which in turn rests on (A3): every section has nonzero pivots at both ends in the limit
   ("singly pinched"). In a doubly pinched section ([S] §8.4 warning) one pivot tends to 0, the coefficients a or b
   blow up, the product form (2.2) degenerates to a product of two small factors equal to an exponentially small
   number ([Astra] §5), and the position of the composition front depends on the length ratios ([N] §11.1).
2. **Why forcing is needed.** With the unforced flows (binary columns exactly balanced, d_p = f_p), the binary
   junction composition J[p, p+1] is unconstrained by the limit relations of both sections of [p, p+1] (the top
   section imposes only trace conditions at its bottom end, the bottom section only at its top end), so DG_∞ has a
   kernel along the p:p+1 edge: (A6) fails. This is the neutral mode of [N] §10.4 (smallest singular value of the
   finite-n Jacobian ~ ρ^n). Forcing p down the bottom section at net rate η_p > 0 makes the junction the pinch of
   that section, an O(1) condition with nonzero derivative (L's binary block).
3. **No sign information on traces is used.** The trace contents and trace net flows at v_n are O(ρ_∞^{n+1}) by
   Theorem A.5; the contents are positive by Lemma A.1, and the trace net flows have whatever signs the balances
   give them (reverse trace flows are allowed in 𝒜).
4. **Uniqueness.** v_n is unique near v_∞ in the variables (d^c, X^c); other steady states with the same flows may
   exist elsewhere ([S] §8.3 caution); they are irrelevant for attainment. The exact system also has non-physical
   zeros far from v_∞ (Agent C, `computation/scripts/spurious_zeros.py`).
5. **(A5) is independent of (A1)–(A4), (A6).** Agent C's red-team case `n3_invalid` (outside D's parameter region, a
   negative binary margin) satisfies (A1)–(A4) and (A6) but not (A5) (a negative pivot and a negative half-profile
   entry in B[2, 3]), and the zeros of G_n converging to v_∞ are non-physical. So Lemma A.7 genuinely needs (A5).
6. **Other product phases or end conventions.** A partial condenser or a vapor side product changes only the affine
   maps of §1 and the end units in Lemma A.1; Theorem A.5 and Lemmas A.6–A.7 are unchanged.

---

## 8. Correspondence with Agent D's hypotheses (H1)–(H8) and Agent L's Theorem L

Against `design/hypotheses.md` v1 (18:07) and L's interface (limit/PROGRESS.md, 18:03):

- **(A1)** ⇐ (H1) (all L_σ, V_σ, D_c, B_c, P_p > 0) and the construction (F1)–(F3): zero duty at J_c is
  V_{T_c} − V_{B_c} = U_c := V_{T[i,m+1]} − V_{B[i−1,m]} (F1, F2), zero duty at a side draw is V_{B[p−1,p]} =
  V_{T[p,p+1]} (F2), and the total balance at J_c is Δ_{T_c} − Δ_{B_c} = 1ᵀf_c (definitions of L_{T_c}, L_{B_c});
  R = V_{B[N−1,N]} = H* + δ (H8).
- **(A2)** ⇐ (H2) (limit net flows), (H5) (all modes real, positive, pairwise distinct; pinch = a designated major
  mode, D's table §2), and D's §2 compositions: X_c > 0 on its support ([Astra] Lemma 1, pinch eigenvectors), product
  compositions (b + d)/P_p ≥ 0 by (H2).
- **(A3)** ⇐ (H8) "nonzero pivots" (with (H3): ℓ_u·X_c = B_c > 0, ℓ_v·X_c = C_c > 0; (H4): A', B' > 0 resp.
  C', D' > 0; G_{X_child}(u_parent) > 0, G_{X_child}(v_parent) < 0; product ends).
- **(A4)** holds by construction. Check against D's mode table (closed-end conditions, with ℓ_ϱ = e_k for trace
  modes at zero trace flow, Lemma A.2 (c)): non-binary T_c: G_{X_c}(r_c) = G_{X_c}(θ_j) = 0 by Lemma 1, X_{c,k} = 0
  for k < i, and the top end J[i, m−1] (a junction, m − 1 > i) is supported on [i, m−1]; non-binary B_c: symmetric;
  binary B-pinned T: X_c supported on {p, p+1}, product p contains no k > p (its majors are p and p−1 if
  π_{p−1} = B); binary B-pinned B: X_c = e_{φ'} so ℓ_{φ''}·X_c = 0, X_c has no k > p+1, product p+1 contains no
  k < p; T-pinned: mirror. So every closed-end condition holds at 𝓛 for every pattern π.
- **(A5)** ⇐ (H6).
- **(A6)** ⇐ Theorem L (`limit/limit.md`), whose inputs (L's list 3(a)–(d)) are contained in (H1)–(H5), (H7).
- **Purity hypothesis of Cor. A.9** ⇐ (H8): the impurity of P_p is (η_{p−1}[π_{p−1} = B] + η_p[π_p = T])/P_p ≤
  2 max_p η_p / min_p z_p once the η's are small; **η_0(δ)** ⇐ D's §4 (parameter region): all non-binary
  η_c ≤ η̄(δ) and the binary η_p in the stated open ranges, which contain (0, η'] for some η' > 0.

(H7) (trace signs κ_σ > 0) is used by L inside (A6), not by the analysis.

The integrated proof (design, Theorem L and this analysis in one notation) is `../theorem.md`; there §6 = this file's
§§3–6, and §3 merges Lemmas A.2–A.3 with limit.md Lemmas 2.1–2.4.

---

## 9. Numerical check of the analysis (evidence only)

`scripts/check_analysis.py` (Python 3, mpmath; WSL Fedora; one thread; working precision 50 + (n_max + 2)·log10(1.5
α_1/α_N) digits, because the forward reconstruction of a section amplifies rounding by (λ_max/λ_min)^{n+1}; at a fixed
60 digits the n = 64 reconstruction lost accuracy — a precision effect, fixed, not a mathematical one). It takes F and
𝓛 from D's library (read-only) and is otherwise independent of D, L and C: generic numerical eigen-data (mpmath `eig`,
analytic normalization by the limit eigenvectors), the network map G_n of (2.2)–§3 in the unknowns (d^c, X^c) only,
damped Newton with finite-difference Jacobian, continuation from the largest n down, reconstruction of every stage
profile, and a re-check of every stage, junction, side-draw, condenser and reboiler balance with the actual streams.
Results (`runs/check_analysis.log`, `runs/check_analysis_n5.log`, **ALL PASS**):

| case | (A1)–(A6) at 𝓛 | σ_min(DG_∞) | n tested | |v_n − v_∞|/ρ_∞^{n+1} | min S, min x | worst balance |
|---|---|---|---|---|---|---|
| N = 3, α = (4,2,1), z = 1/3, q = 1, BB | PASS | 0.018 | 8, 16, 32, 64, mixed [32, 64] | 2.0, 0.83, 0.41, 0.45, 0.16 | ≥ 1.009, > 0 | ≤ 1e-66 |
| N = 3, z = (0.5, 0.2, 0.3), q = −0.3, TB | PASS | 0.026 | same | 2.2, 5.4, 5.4, 5.5, 0.17 | ≥ 1.013, > 0 | ≤ 1e-79 |
| N = 4, α = (6, 3.5, 2, 1), q = 0.6, BBB | PASS | 0.0086 | 8, 16, 32, 48, mixed [32, 64] | bounded (≤ 3.9) | ≥ 1.012, > 0 | ≤ 1e-81 |
| N = 4, same data, BTB | PASS | 0.0098 | same | bounded | > 0 | ≤ 1e-81 |
| N = 5, α = (8,5,3,1.8,1), z = 0.2, q = 0.8, BBBB | PASS | 0.0026 (80×80) | 12, 24, 36, mixed [24, 48] | 1.8, 1.07, 0.87, 0.77 | ≥ 1.012, > 0 | ≤ 1e-75 |

The ratio |v_n − v_∞|/ρ_∞^{n+1} stays bounded, as Theorem A.5 predicts (rate ρ_∞, not n ρ_∞^n); the composite-profile
distance of Lemma A.7 decreases with n (e.g. 0.25, 0.11, 4e-3, 6e-6 for N = 3, n = 8, 16, 32, 64); side-product
impurities tend to the designed η_p/P_{p+1}. Runtime: N = 3 about 20 s, N = 4 about 10 min, N = 5 about 22 min.
