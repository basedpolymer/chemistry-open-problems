# Brief C: the minimal number of complexes for ACR with multistability

Status (27 Oct 2026): **solved modulo the dimension-two theorem (D2) of the coordinator.** The minimal number of
complexes of an at-most-bimolecular network admitting, for some rate constants, absolute concentration robustness
(ACR) together with multistability is **6**. The lower bound 6 is proved below assuming (D2) (Section 2.4); it holds
unconditionally for full-dimensional networks. The bound is attained by an explicit network with 3 species,
6 reactions, 6 reactant complexes, 6 complexes and dimension 3, certified in exact rational arithmetic. The
Kaihnsa–Nguyen–Shiu bound for ACR with nondegenerate multistationarity is 5 complexes, so the conjecture is **true
for the number of complexes** (6 > 5), provided (D2) holds.

(Coordinator's note: (D2) is proved in the manuscript, `../paper.tex`, Theorem "dimension two".)

Correction to the brief (item 2): for m = 6 and dimension 3, two linkage classes are *not* excluded. The case ℓ = 2,
δ = 1 with linkage-class deficiencies δ₁ = δ₂ = 0 escapes the Deficiency One Theorem (which needs Σδ_θ = δ), and it
really occurs: the weakly reversible network N2 below has ℓ = 2, δ = 1, δ₁ = δ₂ = 0, ACR and multistability.

## 1. Setting and definitions

Species X₁, …, X_n; complexes y₁, …, y_m ∈ Z^n_{≥0} (columns of the complex matrix Y); reactions y_i → y_j with rate
constants κ_ij > 0. Laplacian A_κ: A_κ e_i = Σ_{y_i→y_j} κ_ij (e_j − e_i). Mass-action system
dx/dt = f(x) = Y A_κ Ψ(x), Ψ(x) = (x^{y_1}, …, x^{y_m}).

- S = span{y_j − y_i : y_i → y_j}, s = dim S (the dimension); classes (x⁰ + S) ∩ R^n_{>0}.
- ℓ = number of linkage classes, t = number of terminal strong linkage classes, δ = m − ℓ − s; δ_θ = m_θ − 1 − s_θ
  for the θ-th linkage class. Terminal complex = complex in a terminal strong linkage class; weakly reversible ⟺
  every complex terminal.
- Kinetic subspace K = Y(Im A_κ) ⊆ S.
- Bimolecular = at most bimolecular: every complex (reactant and product) has molecularity ≤ 2.
- Nondegenerate positive steady state (KNS §2): f(x*) = 0, x* > 0, Im(df(x*)|_S) = S.
- ACR (KNS Def. 2.6): (G, κ) has a positive steady state and x_i takes the same value at every positive steady state;
  unconditional ACR = ACR for every κ.
- Multistability: for the same κ, some class contains two positive steady states exponentially stable relative to the
  class (all eigenvalues of df|_S with negative real part); they are nondegenerate, so all KNS bounds apply. The
  question is existential in κ.

## 2. Results used from the literature (exact forms)

**2.1 Kernel of the Laplacian** (Feinberg–Horn, ARMA 66 (1977) 83–97; finite Markov chains). For every κ > 0,
dim ker A_κ = t, with a basis of nonnegative vectors whose supports are the terminal strong linkage classes. Hence
rank A_κ = m − t and every vector of ker A_κ vanishes on nonterminal complexes. *Sketch:* A_κ^T is the generator of a
continuous-time Markov chain on the complexes; ker A_κ = invariant measures, spanned by those of the closed classes
(= terminal strong linkage classes), each positive exactly on its class. ∎

**2.2 Deficiency Zero Theorem** (Horn 1972; Feinberg 1972, 1987; KNS Lemma 2.4). If δ = 0: not weakly reversible ⟹
no positive steady state for any κ; weakly reversible ⟹ exactly one positive steady state in every positive class
for every κ.

**2.3 Deficiency One Theorem** (M. Feinberg, *The existence and uniqueness of steady states for a class of chemical
reaction networks*, Arch. Rational Mech. Anal. 132 (1995) 311–370; stated in Feinberg, Chem. Eng. Sci. 42 (1987)
2229–2268). If (i) each linkage class contains exactly one terminal strong linkage class, (ii) δ_θ ≤ 1 for all θ, and
(iii) Σδ_θ = δ, then: if (G, κ) has a positive steady state, every positive class contains precisely one positive
steady state; if G is weakly reversible, a positive steady state exists for every κ. Weak reversibility is not needed
for uniqueness; (i) replaces it (same hypotheses as Joshi–Shiu, arXiv:1412.5257, Thm 3.2).

*Remark.* KNS Lemma 2.5 states the theorem only for weakly reversible networks but their proof of Thm 4.9(2) applies
it to arbitrary networks with ≤ 4 complexes, ℓ = 1, δ = 1; this is justified by Feinberg's general version plus
Lemma 3.1, which forces t = 1 there.

**2.4 KNS** (Eur. J. Appl. Math. 35 (2024) 566–600, arXiv:2307.04186). Thm 1.2: bimolecular + ACR + nondegenerate
multistationarity ⟹ n ≥ 3, ≥ 3 reactant complexes, ≥ 5 complexes, dimension ≥ 2; if full-dimensional, ≥ n + 2
reactant complexes and dimension ≥ 3. Prop. 3.10: one-dimensional networks never have ACR with nondegenerate
multistationarity.

**Hypothesis (D2)** (coordinator): no bimolecular two-dimensional network has, for some κ, ACR and multistability.
With Prop. 3.10: ACR + multistability (bimolecular) ⟹ s ≥ 3.

## 3. Kinetic-subspace lemma

**Lemma 3.1.** If (G, κ) has a nondegenerate positive steady state x* and G has a nonterminal complex, then K = S and
s ≤ m − t − 1, i.e. **δ ≥ t − ℓ + 1**.

*Proof.* K ⊆ S since columns of A_κ are combinations of e_j − e_i over reactions. f takes values in K, so
df(x*)v ∈ K for all v; nondegeneracy gives S = df(x*)(S) ⊆ K. By 2.1, dim K = (m − t) − dim(ker Y ∩ Im A_κ). The
vector A_κΨ(x*) lies in Im A_κ ∩ ker Y (since YA_κΨ(x*) = f(x*) = 0), and is nonzero, because Ψ(x*) > 0 is positive at
a nonterminal complex while every element of ker A_κ vanishes there. So s ≤ m − t − 1 and δ = m − ℓ − s ≥ t − ℓ + 1. ∎

## 4. A general lower bound: m ≥ s + 3

**Proposition 4.1.** If for some κ a positive class contains two distinct nondegenerate positive steady states (any
molecularity), then δ ≥ 1, and δ ≥ 2 when ℓ = 1. Hence **m ≥ s + 3**, with equality only for (ℓ, δ) ∈ {(1,2), (2,1)}.

*Proof.* δ ≥ 1 by 2.2. If ℓ = 1 and δ = 1: either weakly reversible (t = ℓ = 1) or, by Lemma 3.1, 1 ≥ t; so t = 1
and 2.3 applies (δ₁ = δ = 1), giving at most one positive steady state per class, a contradiction. Then
m = δ + ℓ + s with δ + ℓ ≥ 3. ∎

**Corollary 4.2.** Bimolecular with ACR and multistability for some κ: (a) m ≥ 5 (s ≥ 2, KNS Thm 1.2(3); recovers KNS
Thm 4.9(2) with the gap above closed); (b) full-dimensional ⟹ m ≥ 6 (KNS Thm 1.2(5)); (c) under (D2), m ≥ 6.

## 5. One linkage class

**Proposition 5.1.** Let ℓ = 1, reactant complexes y₁, …, y_r, d_i = y_i − y₁ ∈ S, v_i = Σ_{y_i→y'} κ(y' − y_i),
g(u) = Σ_i e^{⟨u,d_i⟩} v_i for u ∈ S, π(x) = orthogonal projection of log x on S.

- (a) f(x) = x^{y₁} g(π(x)).
- (b) For every u ∈ S and positive class P, P ∩ exp(u + S^⊥) is one point; hence every positive class contains
  exactly |Z| positive steady states, Z = g^{-1}(0).
- (c) df(x)|_S = x^{y₁} Dg(π(x)) Q_x at a steady state, with Q_x w = P_S(w/x) symmetric positive definite on S;
  nondegenerate ⟺ Dg(π(x)) invertible.
- (d) ACR in X_k ⟹ e_k ∈ S (X_k in no conservation law) and ⟨u, e_k⟩ is constant on Z; and conversely.
- (e) A nondegenerate positive steady state forces r ≥ s + 1; if r = s + 1 every class has exactly one positive steady
  state. So multistationarity with ℓ = 1 needs r ≥ s + 2.

*Proof.*

- (a) Complexes of one linkage class differ by vectors of S.
- (b) Φ(x) = Σ x_i(log x_i − 1 − u_i) is strictly convex and coercive on cl P, so it has a minimiser there. The
  minimiser is not on the boundary: the directional derivative toward P tends to −∞ there. Interior critical points
  satisfy log x − u ⊥ S, and there is exactly one of them.
- (c) Direct differentiation; ⟨w, Q_x w⟩ = Σ w_i²/x_i > 0.
- (d) log x_k = u_k + w_k with w ∈ S^⊥ arbitrary.
- (e) Dg(u) has image in span{v_i} and kernel ⊇ S ∩ span{d_i}^⊥. If r = s + 1, then d₂…d_r is a basis of S and the
  kernel of λ ↦ Σλ_i v_i is a line ω. Then e^{⟨u,d_i⟩} = ω_i/ω₁ determines u. ∎

## 6. Structure of 6-complex examples

**Proposition 6.1** (assume (D2)). If G is bimolecular with 6 complexes and has ACR and multistability for some κ,
then s = 3 and either:

- **(i)** ℓ = 1, δ = 2, t ≤ 2, r ∈ {5, 6}, and the ACR species is in no conservation law; or
- **(ii)** ℓ = 2, δ = 1, δ₁ = δ₂ = 0, t₁ = t₂ = 1, with linkage classes of sizes (2,4) or (3,3).

Both occur: N1 and N3 are type (i), N2 is type (ii) with sizes (2,4).

*Proof.* Prop. 4.1 with s ≥ 3 gives s = 3 and (ℓ,δ) ∈ {(1,2),(2,1)}. For (i): t ≤ δ + ℓ − 1 = 2 (Lemma 3.1) or t = 1
if weakly reversible; then Prop. 5.1(d),(e). For (ii): Σδ_θ ≤ δ since s ≤ s₁ + s₂; t ≤ 2 forces t_θ = 1; Σδ_θ = 1 = δ
would trigger 2.3, so δ₁ = δ₂ = 0, s_θ = m_θ − 1, s₁ + s₂ = 4, s_θ ≤ 3. ∎

*Remark.* "Σδ_i ≤ 1 so the deficiency-one theorem applies" fails when Σδ_i = 0. Joshi–Shiu (arXiv:1412.5257, §1) note
deficiency-one networks with two deficiency-zero linkage classes where neither theorem applies. KNS's minimal example
{A+B → 2C → 2B, C → A} is of this kind.

## 7. Certified examples

Species X, Y, Z; all three networks are full-dimensional, so the single class is R³_{>0}. Certificates:
`scripts/certify_examples.py`.

### 7.1 Main example N1

Reactions: Z → Y+Z (k₁ = 5), 2Z → Z (k₂ = 1), Y → Y+Z (k₃ = 1), Y+Z → 2X (k₄ = 4), X+Z → 2Z (k₅ = 10), 2X → Y (k₆ = 7).

**Counts:** 3 species, 6 reactions, 6 reactant complexes, 6 complexes (Z, 2Z, Y, Y+Z, X+Z, 2X), dimension 3; ℓ = 1,
δ = 2. Strong linkage classes: {Y, Y+Z, 2X} (terminal), {Z}, {2Z}, {X+Z}; so t = 1, and the network is not weakly
reversible.

    dx/dt = 2k₄yz − k₅xz − 2k₆x²        = 8yz − 10xz − 14x²
    dy/dt = k₁z − k₄yz + k₆x²           = 5z − 4yz + 7x²
    dz/dt = −k₂z² + k₃y − k₄yz + k₅xz   = −z² + y − 4yz + 10xz

**Theorem 7.1.**

- (a) For every κ, dx/dt + 2dy/dt = z(2k₁ − k₅x), so every positive steady state has x = 2k₁/k₅. A positive steady
  state always exists: on x = c, y = (k₁z + k₆c²)/(k₄z), and k₄z·dz/dt = −k₂k₄z³ + k₁k₄z² + (k₁k₃ − k₄k₆c²)z + k₃k₆c²,
  which is positive at 0 with negative leading coefficient. So N1 has unconditional ACR in X.
- (b) For the rates above, x = 1, y = (5z+7)/(4z), and 4z·dz/dt = −(2z−1)(z−1)(2z−7). The positive steady states are
  exactly P₁ = (1, 19/4, 1/2), P₂ = (1, 3, 1), P₃ = (1, 7/4, 7/2).
- (c) Exact Jacobians and characteristic polynomials:
  - P₁: J = [[−33,4,28],[14,−2,−14],[5,−1,−10]], λ³ + 45λ² + 206λ + 30 = (λ+5)(λ² + 40λ + 6). Eigenvalues −5 and
    −20 ± √394 ≈ −39.8494, −0.150567; Routh–Hurwitz 45·206 > 30.
  - P₂: J = [[−38,8,14],[14,−4,−7],[10,−3,−4]], λ³ + 46λ² + 47λ − 50. Eigenvalues −44.9291, −1.71846, +0.647595;
    unstable.
  - P₃: J = [[−63,28,4],[14,−14,−2],[35,−13,−4]], λ³ + 81λ² + 632λ + 1050. Eigenvalues −72.4803, −6.17294,
    −2.34681; Routh–Hurwitz 81·632 > 1050.
  - All three are nondegenerate. Eigenvalue enclosures are in `eigen_enclosures.py` (arb, 200 bits).
- (d) By the implicit function theorem plus (a), ACR and multistability hold on an open set of κ. ∎

### 7.2 N2 (two linkage classes, weakly reversible, 7 reactions)

Reactions: 2Z → X+Z (1), X+Z → 2Z (5), Y+Z → X (2), Y+Z → X+Y (5), X → 2X (3), X+Y → X (1), 2X → Y+Z (3).

Linkage classes {2Z, X+Z} and {Y+Z, X, X+Y, 2X}; t = 2, δ = 1, δ₁ = δ₂ = 0. The only failing hypothesis of the
Deficiency One Theorem is Σδ_θ = δ.

    dx/dt = z² + 7yz + 3x − 5xz − 6x²,  dy/dt = −2yz − xy + 3x²,  dz/dt = −z² − 7yz + 5xz + 3x².

- For every κ, dx/dt + dz/dt = x(k₄ − k₇x), and a positive steady state always exists (constant term k₆k₇c³ > 0,
  leading coefficient −k₁k₂).
- For the rates above, (2z+1)·dz/dt = −(z−1)(2z−1)(z−3) with y = 3/(2z+1).
- Steady states:
  - (1, 3/2, 1/2): stable, λ³ + 20λ² + (117/4)λ + 15/2.
  - (1, 1, 1): unstable, λ³ + 21λ² + 17λ − 6.
  - (1, 3/7, 3): stable, λ³ + 35λ² + 73λ + 30.

### 7.3 N3 (one linkage class, 8 reactions)

Reactions: 2Z → Y+Z (1), Y → X (6), Y+Z → 2Z (4), Y+Z → 2X (1), X → X+Z (1), X+Z → Y+Z (10), X+Z → 2X (3), 2X → X (1).

- dx/dt + dy/dt + dz/dt = x − x².
- Plane equation: −2(z−1)(z−3/2)(z−2).
- Steady states:
  - (1, 1, 1): stable, λ³ + 22λ² + 28λ + 1.
  - (1, 3/2, 2): stable, λ³ + (69/2)λ² + 77λ + 1.
  - (1, 23/18, 3/2): unstable, λ³ + (169/6)λ² + (403/9)λ − 1/2.

### 7.4 Counts

| network | species | reactions | reactant complexes | complexes | dim | ℓ | t | δ | weakly reversible |
|---|---|---|---|---|---|---|---|---|---|
| N1 | 3 | 6 | 6 | 6 | 3 | 1 | 1 | 2 | no |
| N2 | 3 | 7 | 6 | 6 | 3 | 2 | 2 | 1 (δ₁ = δ₂ = 0) | yes |
| N3 | 3 | 8 | 6 | 6 | 3 | 1 | 1 | 2 | no |

ACR mechanism in all three ("weighted integral feedback"): a weight vector w with w_X ≠ 0 and a pair of complexes
{y, y + X}. All reactions out of the other four complexes preserve w·(complex), so w·f = x^y(a − bx).

## 8. Main theorem

**Theorem 8.1.** Let μ be the minimal number of complexes of an at-most-bimolecular network with ACR and
multistability for some κ.

- (a) μ ≤ 6 (N1).
- (b) μ ≥ 5 unconditionally, and μ ≥ 6 under (D2); so **μ = 6 under (D2)**.
- (c) Among full-dimensional networks the minimum is 6 unconditionally.
- (d) Every 6-complex example has dimension 3 and is of type (i) or (ii) of Prop. 6.1; both types occur.

The only possible exception to μ = 6 is a two-dimensional, non-full-dimensional 5-complex example, which (D2)
excludes.

## 9. How the examples were found (exploration only)

- **Inverse LP:** prescribe steady states on x = c, impose f = 0 there (linear in κ), and take LP vertices over the
  reactions allowed by a structural ACR design. These vertices are sparse networks: at most 8 reactions with
  3 points, at most 6 with 2.
- **Scripts:** `explore_n3.jl` found the first hit; `search_lp.py` found 12 hits over 6864 designs (8, 7 and
  6 reactions, both types).
- **Integer rates:** chosen by matching the plane cubic to −a(z−z₁)(z−z₂)(z−z₃).
- **Not found:** no 5-reaction 6-complex network.

## 10. Scripts (`complexes/scripts/`)

- `certify_examples.py`: exact certificate, standard library only. Runtime under 0.1 s, 45 PASS.
- `eigen_enclosures.py`: python-flint enclosures, run in WSL. Runtime under 0.1 s, 9 PASS.
- `simulate.py`: floating-point illustration, 6 PASS.
- `explore_n3.jl` and `search_lp.py`: exploration only (no PASS/FAIL).

## 11. Proved / numerical / open

- **Proved:** Lemma 3.1; Prop. 4.1 and Cor. 4.2; Prop. 5.1; Prop. 6.1; the three certified examples; Thm 8.1(a), (c),
  (d).
- **Conditional on (D2):** μ ≥ 6 for networks that are not full-dimensional.
- **Numerical only:** the searches and `simulate.py`.
- **Open:**
  - (D2) itself (now proved in the manuscript).
  - Whether 5 reactions are possible with 6 complexes.
  - A full classification of 6-complex examples, e.g. whether type (ii) with linkage classes of sizes (3,3) occurs.
