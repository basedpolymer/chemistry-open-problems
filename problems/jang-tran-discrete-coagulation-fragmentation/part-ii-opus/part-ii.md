# Part (ii) of Jang–Tran Conjecture 2 — global mass-conserving solutions for every mass 0 < m ≤ 1

Credit: **Part (ii) (and part (i) as a corollary, and the case m = 1): Claude Opus 5.5 (effort xhigh).**
Refinements adopted from GPT-6-Astra's audit (effort xhigh): the x-variable form of the boundary estimate in
Lemma 5.6(b), the discrete Scheffé argument in Section 7, the optional upper curvature bound.
Status: **complete proof, preliminary** (not yet independently verified). The integrated manuscript is
`../paper.tex`; this file is the self-contained Markdown version of the time-dependent part.

Every step below is a proof unless explicitly labelled *remark* or *numerical evidence*. All algebraic identities
are checked exactly by `scripts/check_identities.py` (sympy, PASS/FAIL); floating-point sanity checks (not part of the
proof) are in `scripts/numerics_check.py`.

---

## 1. The problem

Source: J. Jang and H. V. Tran, *Discrete coagulation–fragmentation equations with multiplicative coagulation
kernel and constant fragmentation kernel*, arXiv:2409.17974v1 (2024); Adv. Contin. Discrete Models **2025**, paper 86,
doi:10.1186/s13662-025-03946-4. Conjecture 2 (arXiv v1, printed p. 7, end of Section 1.4), verbatim up to typesetting:

> **Conjecture 2.** Consider (1.1) with kernels (1.2) and initial data satisfying (1.3).
> (i) For 1/2 < m < 1, (1.1) admits a unique stationary solution ρ̃ : ℕ → [0,∞) such that m₁(ρ̃) = m. Further,
> ρ̃ satisfies the recursive formula (1.9).
> (ii) For 1/2 ≤ m₁(ρ₀) < 1, (1.1) admits a unique mass-conserving solution ρ for all time. This solution ρ is in
> fact a classical solution.

Here (1.2) is a(j,k) = jk, b(j,k) = 1, and (1.3) is m₂(0) = Σ j² ρ₀(j) < ∞. With ρ_n = ρ(n,t), system (1.1) reads

$$
\dot\rho_n=\tfrac12\sum_{i=1}^{n-1} i(n-i)\rho_i\rho_{n-i}-n\rho_n\sum_{j\ge1}j\rho_j-\tfrac{n-1}2\rho_n+\sum_{j>n}\rho_j ,
\qquad n\ge1. \tag{1.1}
$$

We write m_k(ρ) = Σ n^k ρ_n and K(ρ)_n for the right-hand side of (1.1) (with the actual first moment m₁(ρ) in the
coagulation loss term). For a bounded g : ℕ → ℝ, Jang–Tran's weak form (1.5) is

$$
\frac{d}{dt}\sum_n g(n)\rho_n=\tfrac12\sum_{j,k}\big(g(j+k)-g(j)-g(k)\big)jk\rho_j\rho_k-\tfrac12\sum_{j}\sum_{k=1}^{j-1}\big(g(j)-g(k)-g(j-k)\big)\rho_j. \tag{1.5}
$$

**Solution concepts.** Jang–Tran (Definition 1) call ρ a solution if ρ(·,0) = ρ₀ and (1.5) holds for every bounded
g; they do not specify the time regularity. We use:

* a **mass-conserving solution** on [0,∞) with datum ρ⁰ is a family of nonnegative sequences ρ(t) with
  m₁(ρ(t)) = m₁(ρ⁰) for all t ≥ 0, such that for every bounded g the function t ↦ Σ g(n)ρ_n(t) is absolutely
  continuous on compact intervals and (1.5) holds for almost every t. (This is the weakest reading of Definition 1;
  uniqueness is proved in this class.)
* it is **classical** if moreover every ρ_n is C¹ on [0,∞) and (1.1) holds at every t ≥ 0. (Existence is proved in
  this, stronger, class.)
* a **stationary solution** is a time-independent mass-conserving solution, i.e. ρ̃ ≥ 0 with m₁(ρ̃) < ∞ and
  K(ρ̃)_n = 0 for all n.

## 2. Results

**Theorem 1 (global well-posedness, 0 < m ≤ 1).** Let ρ⁰ = (ρ⁰_n) be nonnegative with
m = Σ n ρ⁰_n ∈ (0,1]. No other moment is assumed. Then (1.1) has exactly one mass-conserving solution with datum ρ⁰.
It is classical; every ρ_n belongs to C^∞([0,∞)); t ↦ ρ(t) belongs to C([0,∞); ℓ¹₁) ∩ C¹([0,∞); ℓ¹), where
ℓ¹₁ = {u : Σ n|u_n| < ∞}; (1.5) holds for every bounded g with a classical derivative at every t ≥ 0; and

$$ \textstyle\sum_n\rho_n(t)=m(1-m)+e^{-t/2}\big(\sum_n\rho^0_n-m(1-m)\big). $$

In particular Conjecture 2(ii) holds, with the hypothesis (1.3) removed and the window enlarged to 0 < m ≤ 1.
Together with Jang–Tran's Theorem 1.1 (no global mass-conserving solution when m > 1) this also settles the
discrete form of their Conjecture 1 (critical mass m = 1).

**Theorem 2 (stationary states and relaxation; part (i)).** Let 0 < m < 1. The recursion (1.9),

$$
\tilde\rho_n=\frac{2m(1-m)+\sum_{i=1}^{n-1}i(n-i)\tilde\rho_i\tilde\rho_{n-i}-2\sum_{i=1}^{n-1}\tilde\rho_i}{(2m+1)n+1},
$$

defines strictly positive numbers with Σ ρ̃_n = m(1−m) and Σ n ρ̃_n = m; ρ̃ is a stationary solution, and it is the only
stationary solution with m₁ = m. Every mass-conserving solution with mass m converges to ρ̃ in ℓ¹₁ as t → ∞.
For m = 1 every component of the solution of Theorem 1 tends to 0 although its mass stays equal to 1.

Theorem 2 is Conjecture 2(i) for every m ∈ (0,1) (the range m ≤ 1/2 was Jang–Tran's Lemma 5.7). Its proof
(Section 8) is a corollary of the machinery of Theorem 1. The asymptotics
ρ̃_n ~ c(m) R(m)^{-n} n^{-5/2} (R > 1) are proved in `../paper.tex` (singularity analysis due to GPT-6-Astra, checked
and adapted there).

## 3. Key observations (overview)

1. **κ = 1/τ.** The zeroth-order coefficient of Jang–Tran's Hamilton–Jacobi equation (1.7) is
   κ(x) = 1/2 + 1/(e^x − 1) = ½ coth(x/2) = 1/τ(x) with **τ(x) = 2 tanh(x/2)**. So (1.7) is Tran–Van's continuous
   equation with x replaced by τ(x) in the term F/x. We write k = τ' = 1 − τ²/4 ∈ (0,1], and note 0 < τ(x) ≤ x, τ < 2.
2. **Discrete half-slope invariant.** W = τ(x)(m + F_x) − 2F ≥ 0 holds for the transform of every nonnegative
   sequence of mass m, thanks to the identity (Lemma 5.2)
   τ·j(1+z^j) − 2(1−z^j) = 2(1−z)/(1+z) Σ_{i=1}^{j−1}(1−z^i)(1−z^{j−i}),  z = e^{−x},
   and W ≥ 0 is propagated along characteristics with a *nonnegative* source (Van's continuous invariant has a zero
   source). This upgrades the curvature barrier τF_xx > −1 from m < 1/2 (Jang–Tran, Lemma 3.6) to all m ≤ 1.
3. **No viscosity theory.** Characteristics of (1.7) leave the half-line through x = 0. With the barrier they never
   cross, so the method of characteristics gives a global classical solution F on (0,∞) × [0,∞); no boundary
   condition at x = 0 is needed, and the slope F(x,t)/x → m at x = 0 follows from an ODE inequality along
   characteristics.
4. **Positivity from truncation, mass from comparison.** Finite truncations of (1.1) (coagulation j + k > L
   suppressed) are positive and conserve mass m exactly; their transforms are classical *supersolutions* of (1.7).
   A transport comparison gives F^L ≥ F; any limit point therefore has m₁ ≥ lim F(x,t)/x = m: no mass is lost.
5. **Uniqueness and part (i) are elementary** once m₁ ≡ m: m₀(t) is explicit and (1.1) becomes a lower-triangular
   system ("closed system"), which has at most one solution and whose components converge as t → ∞ to the recursion
   (1.9).

## 4. The closed triangular system; uniqueness

**Proposition 4.1.** Let ρ be a mass-conserving solution with mass m > 0, and n₀(t) = Σ_n ρ_n(t). Then
(a) n₀(t) = m(1−m) + e^{−t/2}(n₀(0) − m(1−m)) for all t ≥ 0;
(b) every ρ_n is C^∞ on [0,∞) and, with c_n = mn + (n+1)/2,

$$
\dot\rho_n+c_n\rho_n=n_0(t)-\sum_{j<n}\rho_j+\tfrac12\sum_{i=1}^{n-1}i(n-i)\rho_i\rho_{n-i}\qquad(t\ge0); \tag{4.1}
$$

(c) two mass-conserving solutions with the same datum coincide.

*Proof.* All series below converge absolutely because Σ nρ_n(t) = m. With g = 1_{n} in (1.5), ρ_n is absolutely
continuous and ρ̇_n = K(ρ)_n a.e., where the coagulation loss is nmρ_n because m₁ ≡ m. With g ≡ 1 the right side of
(1.5) equals −m²/2 + ½Σ(n−1)ρ_n = ½(m − m² − n₀), so n₀ is absolutely continuous and solves this linear equation
a.e., which gives (a). Since Σ_{j>n}ρ_j = n₀ − Σ_{j≤n}ρ_j, (4.1) holds a.e. Its right side involves only n₀ (smooth)
and ρ_1, …, ρ_{n−1}; by induction on n each ρ_n equals its variation-of-constants formula with a continuous (then C^∞)
forcing, so ρ_n is C^∞ and (4.1) holds everywhere. (c): two solutions have the same n₀ (same n₀(0)), and (4.1) for
n = 1, 2, … are linear scalar equations with identical forcings, by induction. ∎

(Jang–Tran proved (c) by viscosity-solution uniqueness, their Proposition 1.2; the argument above is direct.)

## 5. The auxiliary Hamilton–Jacobi equation and its characteristics

Fix m ∈ (0,1] and ρ⁰ ≥ 0 with Σ nρ⁰_n = m. Put

$$ F_0(x)=\sum_n(1-e^{-nx})\rho^0_n,\qquad H_m(x,u,p)=\tfrac12(p-m)(p-m-1)+\frac{u}{\tau(x)}-m , $$

and consider

$$ F_t+H_m(x,F,F_x)=0\quad(x>0,\ t\ge0),\qquad F(x,0)=F_0(x). \tag{5.1} $$

F₀ is C^∞ on (0,∞) (all derivative series converge uniformly on [δ,∞)), F₀ > 0, 0 ≤ F₀' = Σ n e^{−nx}ρ⁰_n ≤ m,
F₀'' ≤ 0, and F₀(x)/x ↑ m as x ↓ 0 (monotone convergence). **We do not assume that the solution of (5.1) is a
transform of a nonnegative sequence**; it is only used as a lower barrier.

**Lemma 5.1 (initial data).** W₀ := τ(m + F₀') − 2F₀ ≥ 0 and −m/e ≤ τF₀'' ≤ 0 on (0,∞).

*Proof.* With z = e^{−x}, τ = 2(1−z)/(1+z), and for every integer j ≥ 1,
(1+z)(1+z+⋯+z^{j−1}) = 1 + 2(z+⋯+z^{j−1}) + z^j, hence

$$ j(1+z^j)-(1+z)\frac{1-z^j}{1-z}=\sum_{i=1}^{j-1}(1+z^j-z^i-z^{j-i})=\sum_{i=1}^{j-1}(1-z^i)(1-z^{j-i}), $$

i.e. τ j(1+z^j) − 2(1−z^j) = 2(1−z)(1+z)^{−1}Σ_{i=1}^{j−1}(1−z^i)(1−z^{j−i}) ≥ 0 (exactly checked for j ≤ 60, A5).
Multiply by ρ⁰_j and sum, using Σ jρ⁰_j = m: W₀ ≥ 0. Next, −τF₀'' = Σ_j jρ⁰_j·(τ j e^{−jx}) and
τ j e^{−jx} ≤ x j e^{−jx} ≤ 1/e. ∎

**Characteristic system.** For a label y > 0 let (X,P,U)(t) = (X,P,U)(t;y) solve

$$
\dot X=P-m-\tfrac12,\qquad \dot P=\frac{k(X)}{\tau(X)^2}U-\frac{P}{\tau(X)},\qquad
\dot U=\tfrac12(P^2+m-m^2)-\frac{U}{\tau(X)}, \tag{5.2}
$$

(X,P,U)(0) = (y, F₀'(y), F₀(y)). These are Ẋ = ∂_pH_m, Ṗ = −∂_xH_m − p∂_uH_m, U̇ = p∂_pH_m − H_m (checks A4a–c).
The vector field is C^∞ on {X > 0} × ℝ²; let [0, t_e(y)) be the maximal interval of existence. Write q = m − P,
W = τ(X)(m+P) − 2U, D = U/τ(X) − P, λ = U/τ(X), η = m − λ, v = q + ½, ε = τ(X)²v/4. Exact identities
(checks A4d, A4h, A4i, A4k):

$$ \dot W+\frac{2+k}{2\tau}W=q(2m-q)(1-k),\qquad 1-k=\tfrac{\tau^2}4, \tag{5.3} $$
$$ W=\tau(q-2D), \tag{5.4} $$
$$ \dot\lambda=\frac{(P-\lambda)^2+(m-\lambda)(1-m+\lambda)}{2\tau}-\frac{\varepsilon\lambda}{\tau},\qquad
   \dot\eta=-\frac{(\eta-q)^2+\eta(1-\eta)}{2\tau}+\frac{\varepsilon\lambda}{\tau}. \tag{5.5} $$

**Lemma 5.2 (invariant region).** For every y > 0 and t ∈ [0, t_e(y)): 0 ≤ P ≤ m, U > 0, W ≥ 0. Consequently
0 < U ≤ mτ(X) < 2m, 0 ≤ λ ≤ (m+P)/2 ≤ m, −(m+½) ≤ Ẋ ≤ −½, y/(m+½) ≤ t_e(y) ≤ 2y, and X(t) ↓ 0 as t ↑ t_e(y).

*Proof.* Fix T < t_e(y). On [0,T], X ≥ δ > 0 and |P|, |U| ≤ A for some δ, A (continuity). Let
V = P⁻ + (P−m)⁺ + U⁻ + W⁻ (a⁺ = max(a,0), a⁻ = max(−a,0)); V is Lipschitz and, a.e.,
V' = −Ṗ 1_{P<0} + Ṗ 1_{P>m} − U̇ 1_{U<0} − Ẇ 1_{W<0}. Put δ' = τ(δ).
- On {P < 0}: −Ṗ = −kU/τ² + P/τ ≤ U⁻/δ'².
- On {P > m}: since 2U = τ(m+P) − W ≤ τ(m+P) + W⁻ and k ≤ 1,
  Ṗ ≤ [k(m+P) − 2P]/(2τ) + W⁻/(2τ²) ≤ (m−P)/(2τ) + W⁻/(2δ'²) ≤ W⁻/(2δ'²).
- On {U < 0}: −U̇ = −½(P² + m(1−m)) + U/τ < 0 (here m ≤ 1 is used).
- On {W < 0}: by (5.3), −Ẇ ≤ −q(2m−q)(1−k). If 0 ≤ P ≤ m then q(2m−q) ≥ 0; if P < 0 then q(2m−q) = m² − (P⁻)²
  ≥ −A P⁻; if P > m then q(2m−q) = −(P−m)⁺(2m + (P−m)⁺) ≥ −(2m+A)(P−m)⁺. So −Ẇ ≤ (2m+A)(P⁻ + (P−m)⁺).
Hence V' ≤ L V a.e. with L = L(δ, A, m), and V(0) = 0 by Lemma 5.1 and 0 ≤ F₀' ≤ m, F₀ > 0. Grönwall gives V ≡ 0.
Then U ≤ τ(m+P)/2 ≤ mτ, λ ≤ (m+P)/2, and the bounds on Ẋ follow. Since U̇ ≥ −U/τ, U(t) ≥ F₀(y)exp(−∫₀^t ds/τ(X)) > 0.
X(t) ≤ y − t/2 gives t_e ≤ 2y; X(t) ≥ y − (m+½)t gives t_e ≥ y/(m+½). If X stayed ≥ δ > 0 on [0, t_e) the solution
would stay in the compact set [δ, y] × [0, m] × [0, 2m] and could be continued; since X decreases, X(t) → 0. ∎

**Lemma 5.3 (strip condition and curvature barrier).** Let J = ∂_yX, Π = ∂_yP, Υ = ∂_yU (continuous on
Ω = {(y,t) : y > 0, 0 ≤ t < t_e(y)} by smooth dependence on initial data). Then Υ = PJ and J > 0 on Ω, and
Γ := τ(X)Π/J satisfies −1 < Γ ≤ √(2m).

*Proof.* Differentiating (5.2) in y, A := Υ − PJ satisfies Ȧ = −∂_uH_m·A = −A/τ (check A4e), and A(0) = F₀'(y) −
F₀'(y)·1 = 0; so A ≡ 0. Using Υ = PJ, J̇ = Π and Π̇ = −J(κ''U + 2κ'P) − κΠ with κ = 1/τ. Let
t₁ = sup{t < t_e : J > 0 on [0,t]} (J(0) = 1). On [0,t₁), r = Π/J satisfies the Riccati equation
ṙ = −r² − κr − 2κ'P − κ''U (A4f) and Γ = τr satisfies (A4g)

$$ \tau\dot\Gamma=-\Gamma^2-\Gamma\big(1+k(q+\tfrac12)\big)-2kD. \tag{5.6} $$

By (5.4) and W ≥ 0, D ≤ q/2; by U ≥ 0, D ≥ −P ≥ −m. Initially Γ(0) = τ(y)F₀''(y) ∈ [−m/e, 0] ⊂ (−1, √(2m)].
At a first time where Γ = −1, (5.6) gives τΓ̇ = k(q + ½ − 2D) ≥ k/2 > 0 (A4j), impossible for a first downward
crossing; at a first time where Γ = √(2m), (5.6) gives τΓ̇ ≤ −2m − √(2m) + 2km < 0, impossible as well. So
−1 < Γ ≤ √(2m) on [0,t₁). If t₁ < t_e, then X ≥ X(t₁) > 0 on [0,t₁] and r ≥ −1/τ(X(t₁)), hence
J(s) = exp∫₀^s r ≥ exp(−t₁/τ(X(t₁))) > 0 on [0,t₁); by continuity J(t₁) > 0 and J > 0 slightly beyond t₁, contradicting
the definition of t₁. So J > 0 on Ω. ∎

The lower barrier is the only place where the half-slope invariant is essential: with the crude bound D ≤ q
(which is all that F ≤ m(1−e^{−x}) gives) the right side of (5.6) at Γ = −1 is k(q + ½ − 2D) ≥ k(½ − q), positive only
when q < ½; this is exactly Jang–Tran's restriction m < 1/2.

**Proposition 5.4 (global classical solution of (5.1)).** For every t ≥ 0 the set Ω_t = {y : t_e(y) > t} is an
interval (a(t), ∞), with a(0) = 0 and a(t) ≥ t/2, and X(t;·) maps it increasingly onto (0,∞). The map
Φ(y,t) = (X(t;y), t) is a C^∞ diffeomorphism from Ω onto (0,∞) × [0,∞), and F := U∘Φ^{−1} satisfies:
F ∈ C^∞((0,∞) × [0,∞)) (one-sided t-derivatives at t = 0), F(·,0) = F₀, F_x = P∘Φ^{−1}, F_t + H_m(x,F,F_x) = 0, and

$$ 0\le F_x\le m,\qquad 0<F\le m\tau(x),\qquad \tau(m+F_x)-2F\ge0,\qquad -1<\tau F_{xx}\le\sqrt{2m}. $$

*Proof.* Ω is relatively open in (0,∞) × [0,∞) (lower semicontinuity of t_e) and Φ is C^∞ with Jacobian determinant
J > 0. *Interval property:* let y < y' with y ∈ Ω_t, and let s* = sup{s ∈ [0,t] : [y,y'] × [0,s] ⊂ Ω}. For s < s*
and η ∈ [y,y'], X(s;η) ≥ X(s;y) ≥ X(t;y) > 0 (monotone in η because J > 0; decreasing in s). So on [0,s*) all
characteristics with labels in [y,y'] remain in the compact set K' = [X(t;y), y'] × [0,m] × [0,2m]; since the
vector field is smooth near K', there is a uniform extension time, so [y,y'] × [0, min(t, s* + c)] ⊂ Ω for some
c > 0; hence s* = t and y' ∈ Ω_t. So Ω_t = (a(t),∞); Ω_t is open, a(t) ∉ Ω_t when a(t) > 0, and y ≤ t/2 implies
t_e(y) ≤ 2y ≤ t, so a(t) ≥ t/2. *Onto:* X(t;y) ≥ y − (m+½)t → ∞ as y → ∞. If X(t;y) ≥ ℓ > 0 for all y ∈ Ω_t, then
for s < t_e(a(t)) ≤ t continuous dependence gives X(s;a(t)) = lim_{y↓a(t)} X(s;y) ≥ lim X(t;y) ≥ ℓ, contradicting
X(s;a(t)) → 0 as s → t_e(a(t)) (Lemma 5.2). So X(t;·) is an increasing bijection (a(t),∞) → (0,∞); hence Φ is a
bijection with positive Jacobian, i.e. a C^∞ diffeomorphism. *Equation:* differentiating U(t;y) = F(X(t;y),t) in y
gives Υ = F_x J, so F_x = Υ/J = P (Lemma 5.3); in t, U̇ = F_xẊ + F_t, so F_t = U̇ − PẊ = −H_m. Differentiating
F_x(X,t) = P in y gives F_xx = Π/J. The bounds are Lemmas 5.2–5.3. ∎

**Lemma 5.5 (comparison).** Let G ∈ C¹((0,∞) × [0,∞)) satisfy G_t + H_m(x,G,G_x) ≥ 0, 0 ≤ G_x ≤ m and
G(·,0) ≥ F₀. Then G ≥ F on (0,∞) × [0,∞).

*Proof.* w = F − G is C¹ and, since H_m(x,u,p) − H_m(x,ũ,p̃) = (p−p̃)((p+p̃)/2 − m − ½) + (u−ũ)/τ,
w_t + b w_x + w/τ ≤ 0 with b = (F_x + G_x)/2 − m − ½ ∈ [−(m+½), −½] continuous. Fix x₁ > 0, t₁ > 0, and let Y
solve Ẏ = b(Y,t), Y(t₁) = x₁, backwards in time (Peano). Backwards Y increases, Y(t) ∈ [x₁ + (t₁−t)/2,
x₁ + (m+½)(t₁−t)], so Y exists on [0,t₁] and stays in (0,∞). Then g(t) = w(Y(t),t) is C¹ with
ġ ≤ −g/τ(Y), so (g exp∫₀^t 1/τ(Y))' ≤ 0 and w(x₁,t₁) ≤ w(Y(0),0)·exp(−∫₀^{t₁}1/τ(Y)) ≤ 0. No boundary condition at
x = 0 is needed because backward transport curves move away from it. ∎

**Lemma 5.6 (slope at the boundary).**
(a) For every m ∈ (0,1] and t ≥ 0, F(x,t)/x → m as x ↓ 0.
(b) If m < 1, then for every t₀ > 0 and δ ∈ (0, min{1, t₀/2}),

$$ 0\le m-\frac{F(x,t)}{\tau(x)}\le K x^{\alpha}\qquad(0<x<\delta,\ t\ge t_0),\qquad
\alpha=\frac{1-m}{2m+1},\quad K=m\delta^{-\alpha}+\frac{2C}{2-\alpha}\delta^{2-\alpha},\quad C=\frac{m(m+\frac12)}4. $$

*Proof.* At t = 0 this is F₀(x)/x → m. Let t > 0, 0 < δ < min{1, t/2}, 0 < x < δ, and follow the characteristic
through (x,t). Its label is y ≥ x + t/2 > δ and X decreases at speed ≥ ½, so it crossed X = δ at a time
s₁ ∈ [t − 2δ, t) ⊂ (0,t), and X ∈ [x,δ] on [s₁,t]. Along this segment regard λ, η as functions of X ∈ [x,δ];
d/dX = −v^{−1} d/ds with ½ ≤ v ≤ m + ½ ≤ 3/2, and ε ≤ (m+½)τ²/4 ≤ (3/8)X², ελ/τ ≤ Cτ ≤ CX (since λ ≤ m, τ ≤ X).

(b) If m < 1: η ≤ m gives η(1−η) ≥ (1−m)η, so by (5.5) η̇ ≤ −(1−m)η/(2τ) + CX. As η ≥ 0 and τ ≤ X,
dη/dX = −η̇/v ≥ (1−m)η/(2vτ) − CX/v ≥ αη/X − 2CX. Integrating (X^{−α}η)' ≥ −2CX^{1−α} over [x,δ] and using η(δ) ≤ m
gives η(x,t) ≤ x^α(mδ^{−α} + 2Cδ^{2−α}/(2−α)). The constant does not depend on t ≥ t₀ when δ < t₀/2.

(a) Let m ≤ 1 and θ ∈ (0,m). Since 1 − m + λ ≥ λ, (5.5) gives λ̇ ≥ (λ/τ)((m−λ)/2 − ε). Choose δ also with
(3/8)δ² ≤ θ/4. Then on [s₁,t]: whenever λ ≤ m − θ, λ̇ ≥ θλ/(4τ) > 0; in particular {λ ≥ m − θ} is never left after it is
entered. While λ ≤ m − θ, dλ/dX ≤ −θλ/(4vτ) ≤ −(θ/6)λ/X, so λ(X) ≥ λ(δ)(δ/X)^{θ/6}. The entry value is
λ(δ) = F(δ,s₁)/τ(δ) ≥ λ_* := min_{s∈[0,t]} F(δ,s)/τ(δ) > 0 (F is continuous and positive). Hence for
x < δ(λ_*/(m−θ))^{6/θ} we get λ(x,t) = F(x,t)/τ(x) ≥ m − θ; also λ ≤ m. As θ is arbitrary and τ(x)/x → 1, the claim
follows. ∎

(For m = 1 the convergence in (a) is not uniform in time, and cannot be: by Theorem 2 the components tend to 0.)

## 6. Positive truncations

For L ≥ 2 consider the truncated system (coagulation events with i + j > L suppressed):

$$
\dot\rho^L_n=\tfrac12\sum_{i=1}^{n-1}i(n-i)\rho^L_i\rho^L_{n-i}-n\rho^L_n\sum_{j=1}^{L-n}j\rho^L_j-\tfrac{n-1}2\rho^L_n+\sum_{j=n+1}^{L}\rho^L_j,
\qquad 1\le n\le L, \tag{6.1}
$$

with ρ^L_n(0) = ρ⁰_n (n < L), ρ^L_L(0) = ρ⁰_L + L^{−1}Σ_{j>L} jρ⁰_j. Put F^L(x,t) = Σ_{n≤L}(1 − e^{−nx})ρ^L_n(t).

**Lemma 6.1.** (6.1) has a unique global solution; it is nonnegative and Σ_{n≤L} nρ^L_n(t) = m for all t.
Moreover F^L ∈ C^∞((0,∞) × [0,∞)) (C¹ in t), 0 ≤ F^L_x ≤ m, F^L(·,0) ≥ F₀, and

$$ F^L_t+H_m(x,F^L,F^L_x)=E_L,\qquad E_L=\tfrac12\sum_{\substack{i,j\le L\\ i+j>L}}ij\rho^L_i\rho^L_j(1-e^{-ix})(1-e^{-jx})\ \ge0. \tag{6.2} $$

*Proof.* The right side of (6.1) is polynomial; it is quasi-positive (ρ^L_n = 0 and ρ^L ≥ 0 imply ρ̇^L_n ≥ 0), so the
nonnegative orthant is invariant; multiplying by n and summing gives d/dt Σ nρ^L_n = 0 (coagulation conserves i + j,
fragmentation conserves size). Nonnegativity and conserved mass bound the solution, which is therefore global. The
initial mass is Σ_{n<L} nρ⁰_n + Lρ⁰_L + Σ_{j>L} jρ⁰_j = m. For (6.2), take g(n) = 1 − z^n (z = e^{−x}) in the weak form of
(6.1). Coagulation: g(i+j) − g(i) − g(j) = −(1−z^i)(1−z^j), so its contribution is
−½(Σ_n n(1−z^n)ρ^L_n)² + E_L = −½(m − F^L_x)² + E_L. Fragmentation: for each n,
Σ_{k=1}^{n−1}(g(n) − g(k) − g(n−k)) = −(n−1)(1+z^n) + 2(z − z^n)/(1−z); summing against −½ρ^L_n and writing
n₀ = Σρ^L_n, Σz^nρ^L_n = n₀ − F^L, Σnz^nρ^L_n = F^L_x gives ½(m + F^L_x) − F^L(1/(1−z) − ½) = ½(m + F^L_x) − F^L/τ.
Hence F^L_t = −½(m − F^L_x)² + ½(m + F^L_x) − F^L/τ + E_L, which is (6.2) (checks A3, A3b; the same computation with
E = 0 is Jang–Tran's (1.7) for the full system, A2). Finally F^L_x = Σ ne^{−nx}ρ^L_n ∈ [0,m], and
F^L(x,0) − F₀(x) = Σ_{j>L} jρ⁰_j[(1−e^{−Lx})/L − (1−e^{−jx})/j] ≥ 0 because s ↦ (1−e^{−sx})/s decreases. ∎

## 7. Proof of Theorem 1

*Uniqueness* is Proposition 4.1(c).

*Step 1 (barrier).* By Lemmas 6.1 and 5.5 (with G = F^L), F^L ≥ F on (0,∞) × [0,∞) for every L.

*Step 2 (compactness).* Extend ρ^L_n = 0 for n > L. Then 0 ≤ ρ^L_n ≤ m/n and |ρ̇^L_n| ≤ m²/2 + m² + m/2 + m. By
Arzelà–Ascoli and a diagonal argument there is L_k → ∞ such that ρ^{L_k}_n → ρ_n locally uniformly on [0,∞) for every
n; ρ_n is continuous, ρ_n ≥ 0, ρ_n(0) = ρ⁰_n, and by Fatou M(t) := Σ nρ_n(t) ≤ m. Since
Σ_{n>N}ρ^L_n(t) ≤ m/(N+1) for all L, the transforms converge: F^{L_k}(x,t) → F^ρ(x,t) := Σ(1−e^{−nx})ρ_n(t).

*Step 3 (no mass loss).* By Step 1, F^ρ ≥ F. As 1 − e^{−nx} ≤ nx, M(t) ≥ F^ρ(x,t)/x ≥ F(x,t)/x for all x > 0;
letting x ↓ 0 and using Lemma 5.6(a), M(t) ≥ m. Hence M(t) = m for every t ≥ 0.

*Step 4 (limit equation).* Fix t. Since ρ^{L_k}(t) → ρ(t) coordinatewise and both have mass m, the mass is tight:
given ϵ > 0 choose N with Σ_{n≤N} nρ_n(t) > m − ϵ; then Σ_{n>N} nρ^{L_k}_n(t) < 2ϵ for k large. Hence for fixed n,
μ^{L_k}_n(t) := Σ_{j≤L_k−n} jρ^{L_k}_j(t) → m. Write (6.1) in integral form on [0,t]. The coagulation gain is a
finite sum; the loss terms nρ^{L}_n μ^{L}_n are bounded by m² and converge pointwise to nmρ_n, so
dominated convergence applies; the fragmentation gains converge uniformly on [0,t] by the number-tail bound. Hence
ρ_n(t) = ρ⁰_n + ∫₀^t K(ρ(s))_n ds for all n, t, where the coagulation loss is nm ρ_n = n m₁(ρ(s)) ρ_n.

*Step 5 (regularity, weak form).* t ↦ ρ(t) ∈ ℓ¹₁ is continuous (discrete Scheffé: coordinatewise convergence of
nonnegative sequences with equal masses m implies Σ n|ρ_n(t_k) − ρ_n(t)| = 2Σ n(ρ_n(t) − ρ_n(t_k))⁺ → 0 by dominated
convergence). K : ℓ¹₁ → ℓ¹ is locally Lipschitz: the coagulation part is bilinear with ℓ¹-norm ≤ (3/2)‖u‖_{1,1}‖v‖_{1,1}
and the fragmentation part is linear with norm ≤ 3/2. So s ↦ K(ρ(s)) ∈ C([0,∞); ℓ¹), the integral identity holds in ℓ¹,
and ρ ∈ C¹([0,∞); ℓ¹) with ρ̇ = K(ρ). For bounded g, t ↦ Σ g(n)ρ_n(t) is therefore C¹ with derivative Σ g(n)K(ρ)_n,
which is the right side of (1.5) after rearranging absolutely convergent double series. Thus ρ is a classical
mass-conserving solution; Proposition 4.1 gives C^∞ components and the formula for n₀. (By uniqueness the whole family
ρ^L converges, not only a subsequence.) ∎

*Remark 7.1 (identification).* F^ρ = F, i.e. the auxiliary characteristic solution is exactly the transform of the
solution: F^ρ is C¹ (ρ ∈ C¹(ℓ¹)) and solves (5.1) by the computation of Lemma 6.1 with E = 0; then w = F^ρ − F solves
w_t + bw_x + w/τ = 0 with w(·,0) = 0, and the transport argument of Lemma 5.5 applied to ±w gives w ≡ 0. Not needed
above. (Numerically, `numerics_check.py` N3b shows F^L − F ↓ 0 as L grows.)

## 8. Proof of Theorem 2 (part (i) as a corollary)

Let 0 < m ≤ 1 and ρ the solution of Theorem 1 for some datum of mass m. By Proposition 4.1, (4.1) holds with
n₀(t) → m(1−m). For a scalar equation u̇ + cu = f(t) with c > 0 and f(t) → f_∞, u(t) → f_∞/c. Applying this to (4.1) for
n = 1, 2, … in turn: ρ_n(t) → ρ̃_n where ρ̃ is given by c_n ρ̃_n = m(1−m) − Σ_{j<n}ρ̃_j + ½Σ i(n−i)ρ̃_iρ̃_{n−i}, which is
exactly (1.9) (2c_n = (2m+1)n + 1). Hence ρ̃_n ≥ 0. For m = 1 this recursion gives ρ̃ ≡ 0.

Now let m < 1. *Mass:* by the number-tail bound, F^ρ(x,t) → F̃(x) := Σ(1−e^{−nx})ρ̃_n as t → ∞, and Σρ̃_n = m(1−m).
From Step 1 of Section 7 (and uniqueness), F^ρ ≥ F, and Lemma 5.6(b) with t₀ = 1 gives
F^ρ(x,t) ≥ τ(x)(m − Kx^α) for 0 < x < δ, t ≥ 1. Letting t → ∞, F̃(x) ≥ τ(x)(m − Kx^α), so
Σ nρ̃_n = lim_{x↓0} F̃(x)/x ≥ m; Fatou gives ≤ m. *Stationarity:* with these two sums, Σ_{j>n}ρ̃_j = m(1−m) − Σ_{j≤n}ρ̃_j,
so (1.9) is the equation K(ρ̃)_n = 0 with m₁(ρ̃) = m: ρ̃ is a stationary solution. *Positivity:* ρ̃₁ = m(1−m)/(m+1) > 0 and,
in (1.9), the numerator equals 2Σ_{j≥n}ρ̃_j + Σ i(n−i)ρ̃_iρ̃_{n−i} ≥ 2(n−1)ρ̃₁ρ̃_{n−1} (n ≥ 2, using the term i = 1 and,
for n ≥ 3, i = n−1; for n = 2 the term is ρ̃₁²), so ρ̃_n > 0 by induction. *Uniqueness:* a stationary solution σ with
m₁(σ) = m is a constant mass-conserving solution, so by Proposition 4.1(a) Σσ_n = m(1−m) and by (4.1) with ρ̇ = 0 it
satisfies (1.9), i.e. σ = ρ̃. *Relaxation:* ρ_n(t) → ρ̃_n for all n and both have mass m, so by Scheffé
Σ n|ρ_n(t) − ρ̃_n| → 0. ∎

## 9. Comparison with the conjecture and with the literature; what is and is not claimed

* Conjecture 2(ii) assumes (1.3) (finite second moment) and 1/2 ≤ m < 1. Theorem 1 needs only a finite first moment,
  and covers 0 < m ≤ 1. Uniqueness is among all mass-conserving solutions in the weak (absolutely continuous) sense of
  Definition 1; existence is classical. The second moment of the solution is not claimed to be finite for t > 0 when
  m₂(0) = ∞ (it is not needed); numerically it is finite and bounded when m₂(0) < ∞ (N4: m₂(2) ≈ 2.04, 4.54, 12.3 for
  m = 0.6, 0.75, 0.9, monomer data).
* Conjecture 2(i): Theorem 2 for every m ∈ (0,1).
* m = 1: global existence and uniqueness of a mass-conserving solution (new; Jang–Tran's Conjecture 1 included
  m = 1), mass escaping to infinity as t → ∞ (consistent with Jang–Tran Corollary 5.2/Remark 1).
* m > 1: no global mass-conserving solution (Jang–Tran Theorem 1.1, not reproved here). Where our argument breaks:
  the face U = 0 of the invariant region needs m ≤ 1 (U̇ = ½(P² + m − m²) at U = 0).
* Relation to Van (arXiv:2607.02293): the invariant W ≥ 0 is the discrete analogue of Van's half-slope invariant
  2M − xM_x ≥ 0 (M = mx − F), with x replaced by τ(x); in the discrete case it has a nonnegative source
  q(2m−q)τ²/4. Unlike Van (and Tran–Van, Jang–Tran), we never use viscosity solutions, vanishing viscosity, or the
  Bernstein property of the HJ solution: positivity comes from the finite truncations, mass conservation from the
  comparison with the characteristic solution.

## 10. Points an independent verifier should check first

1. Lemma 5.3 and Proposition 5.4: the no-crossing argument and the global coverage of (0,∞) by characteristics
   (interval property of Ω_t, surjectivity, C¹ regularity of F up to t = 0).
2. Lemma 5.6(a): uniformity over the family of characteristics arriving near x = 0 at time t (the reference level
   X = δ and the positivity bound λ_*).
3. Lemma 6.1: the identity (6.2), in particular the sign E_L ≥ 0 and F^L(·,0) ≥ F₀ (exact checks A3, A3b).
4. Section 7, Step 4: tightness of the mass at fixed t, which is what allows passing to the limit in the coagulation
   loss term.
5. Section 8: the uniform-in-time bound of Lemma 5.6(b) is what prevents mass from escaping as t → ∞.

## 11. Numerical evidence (not part of the proof)

`scripts/numerics_check.py` (see its header) checks, for m ∈ {0.5, 0.75, 0.9, 0.97} and monomer or heavy-tail
(m₂ = ∞) data: the invariants and the barrier along characteristics (min Γ ≈ −0.33), the bound of Lemma 5.6(b), the
comparison F^L ≥ F (violations at the level of round-off, 1e−15), the decrease of F^L − F with L, agreement of the
closed triangular system with the truncated system (L = 1600, difference ≤ 1e−12), and the long-time limit.
Small truncations (L of a few hundred) are inaccurate for m close to 1: the stationary state then has a long
n^{−5/2} tail up to sizes of order R(m)/(R(m)−1).
