# Brief E: a certified 3-species example, the minimal number of reactions, and four reactant complexes

Status: 27 October 2026. Certificate scripts are in `scripts/`. Numerical exploration only is in `scripts/explore/`.

## 0. Summary

| Question | Answer | Status |
|---|---|---|
| E1: a 3-species bimolecular network with ACR and multistability | G7 (6 reactions, 6 reactant complexes, 7 complexes, κ = (1,1,5,4,1,6)), three rational steady states, two stable | proved (exact) |
| E1: fewest reactions on 3 species | 6. No full-dimensional 3-species bimolecular network with at most 5 reactions has three nondegenerate positive steady states, with or without ACR (verifier: the restriction to full-dimensional networks is necessary, see `../verification/item3_scope_example.py`) | proved (computer-assisted, exact) |
| E1: 6 reactions with 5 reactant complexes | not found; one design family excluded by proof | open |
| E1: 6 complexes (coordinator's lower bound) | not found in this sub-study (found by the complexes sub-study: network N1) | see `../complexes/report.md` |
| E2: exactly 4 reactant complexes | impossible for n = 4 (proved). For n ≥ 5 there are exactly 9 structures (proved); 6 are excluded by proof and 3 (S3, S7, S9) by an exhaustive numerical search over all realizable stoichiometric subspaces | proved except S3, S7, S9 |
| E2: 5 reactant complexes | no example found | open. The minimum is ≥ 4 (KNS with dim 3), ≤ 6 (G7, N1), conjecturally 6 |

"Dimension ≥ 3" (the coordinator's theorem) is used in two places only:

- (i) a 3-species network with ACR and multistability is full-dimensional (Section 2);
- (ii) Section 4 assumes dim S = 3.

## 1. Conventions

- **Network.** Reactions y → y′ with complexes in Z^n_{≥0}, |y|, |y′| ≤ 2. Mass action: f(x) = Σ κ_r x^{y_r}(y′_r − y_r).
  S = span of the reaction vectors.
- **Matrix form.** f = N m(x), where m(x) lists the distinct reactant monomials. The column
  N_a = Σ_{y_r = y_a} κ_r (y′_r − y_r) lies in the cone of the reaction vectors of y_a:
  - entries are ≥ 0 outside supp(y_a);
  - the entry sum is ≤ 0 if |y_a| = 2;
  - all entries are ≥ 0 if y_a = 0.
- **Jacobian.** J = N diag(m) Y^T diag(1/x).
- **ACR** (KNS Def. 2.6). A positive steady state exists, and x_i takes the same value at all positive steady states.
- **Multistability.** For one κ, one compatibility class contains two positive steady states at which J restricted to S
  has all eigenvalues with negative real part. In dimension 3 this means det(λI − J|_S) = λ³ + a1λ² + a2λ + a3 with
  a1 > 0, a3 > 0 and a1a2 > a3. In particular det J|_S < 0.
- **Lemma 1.1 (kinetic subspace).** At a nondegenerate positive steady state, Im N = S. Hence every reaction vector
  lies in Im N.
  - *Proof.* Im J ⊆ Im N ⊆ S, and J|_S is invertible. ∎

## 2. Three species, five reactant complexes: Gale-dual reduction

Here n = 3, S = R³, there are exactly 5 distinct reactant complexes, and κ is fixed.

**Proposition 2.1.** Suppose (G, κ) has a nondegenerate positive steady state. Then rank N = 3, rank [Y; 1^T] = 4, and
C = ker N ∩ R^5_{>0} is a nonempty 2-dimensional cone. Notation:

- r1, r2 are the extreme rays of the closure of C, and μ(s) = (1 − s) r1 + s r2;
- c ∈ Z^5 is the primitive affine dependence (Σ c_i y_i = 0, Σ c_i = 0);
- g(s) = Σ c_i log μ_i(s) on (0, 1).

Then:

1. Positive steady states correspond bijectively to the zeros of g in (0, 1), via m(x) ∈ R_{>0} μ(s).
2. det J(x(s)) = σ · Q(s) · g′(s), with Q > 0 and σ = ±1 independent of s.

*Proof.*

- rank N ≤ 2 gives degeneracy (KNS Prop. 3.1(1)).
- If y_i · w = β for all i (w ≠ 0), then f(x∘e^{tw}) = e^{βt} f(x), so J(x)(x∘w) = 0 at a steady state, which is
  degenerate.
- m(x) = αμ(s) means Y^T ξ − (log α)1 = log μ(s), with ξ = log x. The matrix [Y^T | −1] has rank 4 and left kernel
  Rc, so the system is solvable iff c · log μ(s) = 0, and then uniquely. This gives 1.
- For 2, use complementary-minor (Grassmann) duality. For N with kernel basis [μ(s), v], v = r2 − r1, the 2×2 minors
  of the kernel basis do not depend on s. Apply the same duality to Y, with ker Y = span(c, c′) and c′ · 1 = 1.
  Cauchy–Binet then gives

  det(N diag(αμ) Y^T) = γη α³ Π μ_i · Σ_{a<b} (v_b/μ_b − v_a/μ_a)(c′_a c_b − c′_b c_a)
  = γη α³ Π μ_i · det([1, v/μ]^T [c′, c]) = γη α³ Π μ_i · g′(s).

  Finally det J = det(N diag(λ) Y^T) / Π x_k. ∎

(Exact check: `scripts/check_gale_formula.py`, 40 random instances, PASS.)

**Corollary 2.2.** Multistability implies that g has ≥ 3 zeros and that g′ changes sign ≥ 2 times.

- *Proof.* Both stable states have det J < 0, so g′ has the same sign at them and g has a zero in between. With at
  most one sign change of g′, g would be strictly monotone on at most two pieces and could not vanish three times. ∎

**Proposition 2.3 (degree bound).** For every κ there are at most D(c) = Σ_{c_i>0} c_i positive steady states, or a
continuum.

- *Proof.* The steady states are the roots of Π_{c_i>0} μ_i^{c_i} − Π_{c_i<0} μ_i^{−c_i}, a polynomial in s of degree
  ≤ D. ∎
- Among the 54 classes (up to relabelling species) of 5-sets of bimolecular complexes: D = 2 for 33 (excluded for any
  number of reactions); 6 are degenerate (rank [Y; 1^T] < 4); D = 3 for 13; D = 4 for 2.

**Proposition 2.4 (ACR structure).** Let P be the number of distinct zeros ("poles") of the forms λ ↦ λ_i on
P(ker N). Suppose there is ACR in X_j together with at least 3 steady states. Then P ∈ {4, 5}. If P = 4, the unique
pair {a, b} with a common pole satisfies y_a − y_b = ±e_j, and the ACR is structural.

*Proof.*

- The numerator of g′ has degree ≤ P − 2, so P ≥ 4.
- On the steady states log x_j = h(s) + const, where h = Σ d_i log μ_i, Σ d_i = 0 and Σ d_i y_i = e_j.
- If the pole-weight vectors of c and d are proportional, some e = d − βc ≠ 0 has zero weight on each pole class:
  - P = 5 forces e = 0, a contradiction;
  - P = 4 forces e = ε(e_a − e_b), so ε(y_a − y_b) = e_j.
- Otherwise, for P = 4, kill the weight of one pole: φ = ag + bh has ≤ 3 poles, so φ′ has ≤ 1 zero and φ never takes
  a value three times. But ACR forces equal values of φ at the three steady states. ∎

**Theorem 2.5 (five reactions are not enough).** No full-dimensional at-most-bimolecular 3-species network with at
most 5 reactions has, for any κ, three positive steady states, unless they are all degenerate. (Verifier's note:
the word "full-dimensional" was missing; a two-dimensional 4-reaction counterexample to the unrestricted statement is
in `../verification/item3_scope_example.py`.) So none is multistable. With dim ≥ 3,
ACR plus multistability on 3 species needs ≥ 6 reactions, and 6 is attained (Theorem 3.1).

*Proof.*

- With at most 4 distinct reactant complexes, use KNS Prop. 3.1(3) and 3.5.
- Otherwise N = Γ diag(κ), and g′ does not depend on κ.
- `scripts/enum5_reactions.py` enumerates the 54 classes of reactant 5-sets × 9^5 product choices = 3,188,646
  networks, all decided exactly:

| Stage | Networks |
|---|---|
| rank 3 | 1,251,698 |
| consistent | 628,446 |
| rank [Y; 1^T] = 4 | all consistent ones |
| poles P = 2 / 3 / 4 / 5 | 36,303 / 265,885 / 253,227 / 73,031 |
| sign changes of g′ on (0,1), P ≥ 4 (exact sympy root counting) | 0 for 298,879; 1 for 27,379; never ≥ 2 |

- For P ≤ 3 there is at most one sign change, by degree. ∎

## 3. The certified example (E1)

**Theorem 3.1.** The network G7 on X, Y, Z:

| | Reaction | Rate constant |
|---|---|---|
| R1 | X+Y → 2Y | k1 = 1 |
| R2 | Y → 2Z | k2 = 1 |
| R3 | 0 → X+Z | k3 = 5 |
| R4 | X+Z → 2Z | k4 = 4 |
| R5 | 2Z → X+Z | k5 = 1 |
| R6 | Y+Z → Y | k6 = 6 |

Counts: 3 species, 6 reactions, 6 reactant complexes {X+Y, Y, 0, X+Z, 2Z, Y+Z}, 7 complexes (adding 2Y), 2 linkage
classes, 2 terminal strong linkage classes, dimension 3, deficiency 2.

It has ACR in X for every κ. For κ = (1, 1, 5, 4, 1, 6) it has exactly three positive steady states, all in the
single class R³_{>0}:

| Steady state (x, y, z) | a1 | a2 | a3 | a1a2 − a3 | Eigenvalues | Type |
|---|---|---|---|---|---|---|
| (1, 2, 1) | 16 | 70 | 4 | 1116 | −0.0579065…, −7.97105 ± 2.35356i | exponentially stable |
| (1, 5/4, 3/2) | 55/4 | 435/8 | −5/8 | 23945/32 | +0.0114610…, −6.88073 ± 2.68109i | saddle (a3 < 0) |
| (1, 10/9, 5/3) | 124/9 | 470/9 | 20/27 | 58220/81 | −0.0142378…, −6.88177 ± 2.16044i | exponentially stable |

Eigenvalue enclosures are certified with arb, radius ≤ 1e−7. Both ACR and multistability persist on an open set of κ.

*Proof.* The system is:

- dx/dt = −k1xy + k3 − k4xz + k5z²
- dy/dt = y(k1x − k2)
- dz/dt = 2k2y + k3 + k4xz − k5z² − k6yz

Steps:

- **ACR.** At a positive steady state y > 0, so x = k2/k1 for every κ. Rank 3 means the only class meeting R³_{>0} is
  R³_{>0}.
- **Reduction to a cubic.** At x = k2/k1:
  - adding the first and third equations gives y(k6z − k2) = 2k3;
  - the first equation gives k2y = k5z² − (k4k2/k1)z + k3.
  - With ζ = k6z/k2, the steady states are the roots ζ > 1 of Pζ³ − (P+Q)ζ² + (Q+1)ζ − 3 = 0, with
    P = k5k2²/(k6²k3), Q = k4k2²/(k1k6k3) and y = 2k3/(k2(ζ − 1)).
- **Roots for this κ.** P = 1/180, Q = 2/15, and the cubic equals (ζ − 6)(ζ − 9)(ζ − 10)/180. Equivalently, the lex
  Gröbner basis at x = 1 is {y − (z² − 4z + 5), 6z³ − 25z² + 34z − 15}.
- **Stability.** J = [[−k1y − k4z, −k1x, −k4x + 2k5z], [k1y, k1x − k2, 0], [k4z, 2k2 − k6z, k4x − 2k5z − k6y]]. At the
  three rational points its characteristic coefficients are the rationals in the table. Routh–Hurwitz holds exactly
  at the first and third point. a3 < 0 at the second gives a positive real eigenvalue.
- **Open set.** The roots are simple and det J ≠ 0, so the implicit function theorem gives an open set of κ with the
  same signs. ACR is structural. ∎

Scripts:

- `scripts/certify_example.py`: exact Gröbner basis, exact root isolation, arb enclosures; 20 PASS, < 1 s.
- `scripts/simulate_example.py`: illustration only (floating point).

**Minimality of G7.**

- Reactions: 6 is optimal (Theorem 2.5).
- Reactant complexes: 5 is open (Section 5).
- Complexes: 7 here; 6 is attained by N1 of the complexes sub-study.

## 4. Exactly four reactant complexes (E2)

Setting: n species, 4 reactant complexes, dim S = 3, so k = n − 3 ≥ 1. A has rows y_a − y_4.

**Proposition 4.1 (toric structure).** At a nondegenerate positive steady state:

- rank N = 3, ker N = Rp with p > 0, and rank A = 3;
- the positive steady states form the coset x* ∘ exp(ker A);
- ACR in X_i ⟺ e_i ⊥ ker A;
- in the class of x*, the steady states are x* ∘ e^{Vt} with W^T(x* ∘ (e^{Vt} − 1)) = 0, where V is a basis of ker A
  and W a basis of S^⊥;
- all reaction vectors lie in Im N = S (Lemma 1.1).

**Theorem 4.2 (n = 4).** No κ gives ACR together with three positive steady states in one class.

*Proof.*

- In a class ω · x = T, the steady states are the zeros of φ(t) = Σ_j ω_j x*_j e^{tv_j} − T.
- Write f = x^{y_4} N_123 (e^{Aξ} − c). Then J|_S ≅ x^{y_4} diag(c) A diag(1/x) N_123.
- Cauchy–Binet, with minors of A ∝ (−1)^j v_j and minors of N_123 ∝ (−1)^j ω_j, gives
  det J|_S = σ · (positive) · φ′(t).
- So two stable states force three zeros of φ.
- Laguerre's rule of signs for exponential sums bounds the number of zeros by the number of sign changes of the
  coefficients, ordered by exponent.
- The ACR species has v_i = 0 and merges with the constant term. By Lemma 4.3 there are at most 3 distinct
  exponents, hence at most 2 zeros. ∎

**Lemma 4.3.** If v_i = 0, the other three entries of v are not distinct and nonzero.

*Proof.* Every reactant complex satisfies y · v = β.

- **β = 0.** The only candidates are 0, X_i, 2X_i and at most one X_m + X_p. Since 0, X_i, 2X_i are collinear,
  rank A < 3.
- **β ≠ 0.** The candidates are:
  - X_m and X_i + X_m, where v_m = β (one index);
  - 2X_{m′}, where v_{m′} = β/2;
  - at most one X_p + X_q, where v_p + v_q = β.

  Four complexes need all of them. Then {p, q} = {m′, l}, so v_l = v_{m′}, a contradiction. ∎

The machine check (`scripts/e2_reactant_sets.py`) finds 26 sets with an ACR species and none with three distinct
nonzero exponents.

**Proposition 4.4 (n ≥ 5).** Every species must occur in a reactant complex, so n ≤ 8. The exact list
(`scripts/e2_reactant_sets.py`, 8 s):

| | Reactant complexes | ACR species | Exponent rows of ker A |
|---|---|---|---|
| S1 | 0, E, C+D, A+B | E | A(−1,0) B(1,0) C(0,−1) D(0,1) |
| S2 | 0, 2E, C+D, A+B | E | as S1 |
| S3 | E, 2E, C+D, A+B | E | as S1 |
| S4 | E, D, C+E, A+B | C | A(−1,1) B(1,0) D(0,1) E(0,1) |
| S5 | E, D+E, C+E, A+B | C, D | A(−1,1) B(1,0) E(0,1) |
| S6 | E, D+E, C+D, A+B | D | A(−1,1) B(1,0) C(0,1) E(0,1) |
| S7 | E, D+E, 2C, A+B | D | A(−1,1) B(1,0) C(0,1/2) E(0,1) |
| S8 | E, D+E, B+C, A+C | D | A(−1,1) B(−1,1) C(1,0) E(0,1) |
| S9 (n = 6) | F, E+F, C+D, A+B | E | A(−1,0,1) B(1,0,0) C(0,−1,1) D(0,1,0) F(0,0,1) |

There are no structures for n = 7 or 8.

**Proposition 4.5 (S4, S5, S6, S8).** Each class has at most 2 positive steady states, or a continuum, for every N.

*Proof.*

- The class equation is Q1(W/U − 1) + Q2(U − 1) + Q3(W − 1) = 0 in R², with U = e^{t1}, W = e^{t2}.
- So (W/U, U, W, 1) ∈ ker[Q1, Q2, Q3, −ΣQ] = span{(1,1,1,1), (ρ, 0)}, which gives m_i = 1 + βρ_i.
- The monomial relation m1 m2 = m3 becomes (1 + βρ1)(1 + βρ2) = 1 + βρ3, which has at most two roots β. ∎

**Proposition 4.6 (S1, S2).** A class containing a nondegenerate positive steady state contains no other positive
steady state, for every N.

*Proof.* Write f = N_0 + N_E x_E + N_CD x_C x_D + N_AB x_A x_B. A second steady state is x* + z with z ∈ S, z_E = 0
and the two product equations.

- **Eliminate θ_E.** The row-E balance gives e = −N_E[E] > 0 (e = 0 makes x_E conserved, hence degenerate).
- **Sign constraints.** After eliminating θ_E, z_i = ℓ_i(σ, τ) for i ∈ {A, B, C, D}. Bimolecularity gives
  ℓ_A, ℓ_B ≥ 0 at (1, 0) and ℓ_C, ℓ_D ≥ 0 at (0, 1). The balance gives ℓ_i(π) ≤ 0 for all four, with
  π = (p_CD, p_AB) > 0.
- **Disjoint wedges.**
  - (u + ℓ)(v + ℓ′) = uv with all factors positive forces ℓℓ′ < 0 or ℓ = ℓ′ = 0.
  - {ℓ_A ℓ_B < 0} lies in cone((1, 0), π) ∪ its negative; {ℓ_C ℓ_D < 0} lies in cone(π, (0, 1)) ∪ its negative. They
    meet only on Rπ, where z = 0.
- **Degenerate leftover.** The only remaining cases are ℓ_A ≡ ℓ_B ≡ 0 or ℓ_C ≡ ℓ_D ≡ 0, which make the Jacobian
  singular. ∎

S2 is the same with N_2E in place of N_E.

**S3, S7, S9 (numerical only).**

- With unconstrained N, S7 and S9 do have multistable classes (`explore/general_direct.py`, logs gd_6.log,
  gd_8.log). These N are not realizable by Lemma 1.1.
- **Realizable subspaces.** `explore/list_subspaces.py` lists every S spanned by reaction vectors that uses all four
  reactant complexes and admits a positive flux.
- **Class counts.** `explore/count_gen.py` samples 120–300 classes per S. It then:
  - samples fluxes and requires rank N = 3;
  - tests stability whenever 3 states occur.
- The number of steady states in a class depends only on S, so each count covers every network with that
  stoichiometric subspace.

| Structure | Valid S | Kinetically degenerate | Class-degenerate | Max steady states per class (rest) |
|---|---|---|---|---|
| S3 | 242 | 43 | 6 | 2 (193 subspaces) |
| S7 | 305 | 23 | 32 | 2 (250 subspaces) |
| S9 | 859 | 97 | 144 | 2 (618 subspaces) |

No class with 3 states was found and there were no two-stable hits. For S7, a second resultant-based solver agrees
(max 2).

**Conclusion for E2.**

- 4 reactant complexes are impossible for n = 4 (Theorem 4.2) and for S1, S2, S4, S5, S6, S8 (Propositions 4.5,
  4.6). They are numerically impossible for S3, S7, S9.
- So 4 is conjecturally not attained.
- No example with 5 reactant complexes was found; the smallest known has 6 (G7, N1).

## 5. Partial results on 5 reactant complexes, 3 species

1. Only 15 reactant 5-sets can carry three steady states (Proposition 2.3).
2. Structural ACR (P = 4) needs a reduced circuit of degree 3. With x_j fixed at its ACR value, two equations in four
   reduced monomials remain. Only {1, y, z², yz} (or its coordinate swap) and {y, z, y², z²} qualify.
3. **Proved.** The design X+Y → 2Y, Y → W with reactant complexes {X+Y, Y, 0, 2Z, Y+Z} has at most 2 positive steady
   states, for any κ and any number of reactions.
   - At x = x*, the steady states are the zeros on (0, z_max) of H − K, where H = (B1 − β2z)/(α1 − A2z) and
     K = (β3z² − B4)/(A3z² + A4).
   - H is a Möbius function (decreasing, or increasing and convex) and K is nondecreasing.
   - If H decreases there is one zero. If H increases, H − K > 0 at both ends, and the zeros are those of a cubic, so
     at most two are simple.
4. Numerical searches found no multistability:
   - `explore/design5rc.py`: all 22,500 six-reaction designs with reactant complexes {X+Y, Y, c0, 2Z, Y+Z},
     c0 ∈ {0, X, 2X}, each with 180,000 values of κ.
   - `explore/design_fast.py`: 725,000 designs with at most 5 reactant complexes.

## 6. Proved / numerical / open

**Proved:**

- Theorem 3.1 (G7, refuting the species part of the conjecture: 3 species suffice).
- Propositions 2.1–2.4 and Theorem 2.5 (minimum 6 reactions on 3 species, given dim ≥ 3).
- Theorem 4.2 with Lemma 4.3; Propositions 4.4, 4.5, 4.6.
- Section 5, item 3.

**Numerical only:**

- S3, S7, S9 cannot be multistable.
- No design with 6 reactions and 5 reactant complexes was found.
- The simulation.

**Open:**

- Proofs for S3, S7, S9.
- Existence of any example with 5 reactant complexes.
- 6 reactions with 5 reactant complexes.

## 7. How to run

Run from WSL in `problems/kaihnsa-acr-multistability/examples/scripts` (or via `outils/wsl "…"`):

| Command | What it does | Runtime |
|---|---|---|
| `python certify_example.py` | Theorem 3.1 | < 1 s |
| `python enum5_reactions.py 12` | Theorem 2.5 | ≈ 2 min |
| `python e2_reactant_sets.py 8` | Lemma 4.3 check, Proposition 4.4 | 8 s |
| `python check_gale_formula.py` | Proposition 2.1(2) check | < 1 s |
| `python simulate_example.py` | illustration | 1 s |

Exploration scripts and their logs are in `scripts/explore/` (floating point, not certificates).
