# Cases 9 and 10 of Banaji–Boros–Hofbauer: no Andronov–Hopf bifurcation implies no periodic orbit

Status: **solved (preliminary)**. The argument below is complete, but no independent verifier has checked it yet.
Solver: Claude Opus 5.5 (effort xhigh), 24 October 2026.

## 0. Statement

**Source.** M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*,
Stud. Appl. Math. 152 (2024) 249–278, doi:10.1111/sapm.12639, arXiv:2304.02303. The open question is in §6
(Conclusions, p. 22 of arXiv v1, identical text in the published version):

> "However, we leave it open whether the networks in cases 9 and 10 admitting no Andronov–Hopf bifurcations can have
> a limit cycle. For example, we do not know whether the system 2X → 3X + 2Y (κ₁), X + Y → 0 (κ₂), Y → X + Y (κ₃),
> ẋ = κ₁x² − κ₂xy + κ₃y, ẏ = 2κ₁x² − κ₂xy admits a limit cycle."

"Cases 9 and 10" refer to arXiv Theorem 9, which is Theorem 3 of the published version, with the same text (checked
on the Europe PMC full text of PMC10953466). Throughout, BBH means this paper. Lemma and theorem numbers are those of
arXiv v1. The published numbers are: arXiv Lemma 2 = published Lemma 1, Lemma 3 = Lemma 2, Lemma 4 = Lemma 3,
Theorem 9 = Theorem 3, Theorem 10 = Theorem 4. Equation numbers agree.

**Theorem A (main result).** Take a three-reaction mass-action network on two species X, Y whose source complexes are
2X, X+Y, Y (case 9) or 2X, X+Y, 0 (case 10). Target complexes are arbitrary, with no molecularity bound. Suppose the
network does not satisfy the Andronov–Hopf condition of BBH Theorem 9 (published Theorem 3) for its case. Then, for
every choice of rate constants κ ∈ ℝ³₊, the mass-action system has **no periodic orbit** in the closed quadrant
ℝ²_{≥0}. In particular it has no limit cycle.

The same holds after exchanging X and Y; BBH list the cases up to this exchange.

**Corollary B (classification).** For a network in case 9 or 10, consider the statements:
- (a) it has a periodic orbit for some κ;
- (b) it has a limit cycle for some κ;
- (c) it satisfies BBH's Andronov–Hopf condition.

Then (a) ⟺ (c), and (b) ⟹ (a). Theorem A gives (a) ⟹ (c); (b) ⟹ (a) is trivial. The converse (c) ⟹ (b) is BBH's
claim that the Andronov–Hopf bifurcation is supercritical. They obtain it from a computation of the first focal value
in their Mathematica notebook; we did not re-check it. Granting that claim, (a), (b) and (c) are all equivalent.

The example of §6 is the case k = 1 of family 9a below. It was already settled in this repository, result 2,
`problems/bbh-three-reaction-no-limit-cycle/`; Theorem A re-proves it by a different route.

**Answer to the question: no.** No network of case 9 or 10 that admits no Andronov–Hopf bifurcation has a limit cycle,
or even a periodic orbit, for any rate constants.

## 1. Setting and notation

We follow BBH (13)–(18). Reaction i (i = 1, 2, 3) is

  aᵢX + bᵢY → (aᵢ + cᵢ)X + (bᵢ + dᵢ)Y,

with rate constant κᵢ > 0, where (aᵢ, bᵢ) is the source and (cᵢ, dᵢ) ≠ (0, 0) is the reaction vector. The cᵢ, dᵢ are
integers with aᵢ + cᵢ ≥ 0 and bᵢ + dᵢ ≥ 0. The sources are:
- case 9: (a, b) = (2,0), (1,1), (0,1);
- case 10: (a, b) = (2,0), (1,1), (0,0).

The admissible ranges are therefore

  c₁ ≥ −2, d₁ ≥ 0, c₂ ≥ −1, d₂ ≥ −1, and
  (case 9) c₃ ≥ 0, d₃ ≥ −1;  (case 10) c₃ ≥ 0, d₃ ≥ 0.

The mass-action system is ẋ = Σ cᵢκᵢmᵢ, ẏ = Σ dᵢκᵢmᵢ, with monomials m = (x², xy, y) in case 9 and
m = (x², xy, 1) in case 10. Set u = c × d, that is u₁ = c₂d₃ − c₃d₂, u₂ = c₃d₁ − c₁d₃, u₃ = c₁d₂ − c₂d₁, as in
BBH (16). Set A := c₁u₁ − c₃u₃; this is not to be confused with the function A(x) of §3.

**BBH's Andronov–Hopf condition (Theorem 9, verbatim).**
- Case 9: either
  c₁ > 0, c₂ = −1, c₃ > 0, d₁ > 0, d₂ ≥ −1, d₃ ≥ −1 and ½(d₃/c₃ + d₁/c₁) < d₂/c₂ < d₁/c₁,
  or
  c₁ > 0, c₂ = −1, c₃ = 0, d₁ > 0, d₂ ≥ −1, d₃ = −1 and d₂/c₂ < d₁/c₁.
- Case 10: c₁ > 0, c₂ = −1, c₃ > 0, d₁ > 0, d₂ = −1, d₃ ≥ 0 and ½(d₃/c₃ + d₁/c₁) < d₂/c₂ < d₁/c₁.

**Facts used throughout.**

(F1) *Boundary.* The open quadrant ℝ²₊ is forward invariant. If cᵢ < 0 then aᵢ ≥ 1, so every negative term of ẋ
carries the factor x, and similarly for y; hence x(t) ≥ x(0)·exp(−∫h) > 0 for some continuous h (BBH §2.4).
Consequently a periodic orbit Γ ⊂ ℝ²_{≥0} that meets ℝ²₊ lies entirely in ℝ²₊: it is the forward orbit of any of its
points. A periodic orbit contained in the boundary ∂ℝ²_{≥0} is impossible, because the boundary is homeomorphic to ℝ
and contains no Jordan curve. **So every periodic orbit in ℝ²_{≥0} lies in ℝ²₊.**

(F2) *Positive equilibria.* At a positive equilibrium the flux vector (κᵢmᵢ) is a positive vector in ker Γ, where
Γ = (c; d). If rank Γ = 2, then ker Γ = span(u), so a positive equilibrium exists only if u ∈ ℝ³₊ ∪ ℝ³₋. In that case
κᵢmᵢ(x̄, ȳ) = μuᵢ with μuᵢ > 0, and the equilibrium is unique:
- case 9: x̄ = κ₃u₂/(κ₂u₃), ȳ = κ₁κ₃u₂²/(κ₂²u₁u₃);
- case 10: x̄ = (κ₃u₁/(κ₁u₃))^{1/2}, ȳ = κ₃u₂/(κ₂u₃x̄).

(This is BBH Lemma 2; `certify_lienard.py` S1 checks the formulas.)

(F3) *Normalisation.* Assume u ∈ ℝ³₊ ∪ ℝ³₋. Put x = x̄X, y = ȳY and t = (x̄/μ)T. The system becomes

  X' = Σ cᵢuᵢMᵢ(X, Y),  Y' = ρ Σ dᵢuᵢMᵢ(X, Y),  with ρ := x̄/ȳ > 0,

where Mᵢ = mᵢ. The equilibrium is now (1, 1). This is a linear change of coordinates and a constant rescaling of time,
so it preserves orbits and periodicity. As κ ranges over ℝ³₊, ρ ranges over (0, ∞). At (1, 1),

  det J = ρu₁u₂u₃ (case 9), 2ρu₁u₂u₃ (case 10);
  tr J = A − ρd₁u₁ (case 9), A + ρd₂u₂ (case 10).

Both are checked in S1. So sgn det J = sgn u, which is BBH (17). Below we drop the capitals and write x, y for X, Y.

(F4) *Standard planar tools.* We use three standard facts.
- A periodic orbit Γ of a C¹ planar vector field is a Jordan curve. The equilibria in its bounded complementary
  component have indices summing to +1 (Poincaré index theory), and a saddle has index −1.
- If Γ lies in the convex set ℝ²₊, then its interior lies in ℝ²₊. Indeed ℝ² ∖ ℝ²₊ is connected, unbounded and disjoint
  from Γ, so it lies in the unbounded component of ℝ² ∖ Γ.
- Bendixson–Dulac: if B is C¹ and div(BF) < 0 everywhere (or > 0 everywhere) on a simply connected domain, then
  there is no periodic orbit in that domain.

## 2. Reduction to three families

**Lemma 2.1 (classes without periodic orbits).** A case-9 or case-10 network has no periodic orbit in ℝ²_{≥0}, for any
κ, if it satisfies any of (i)–(v):
- (i) rank Γ ≤ 1;
- (ii) rank Γ = 2 and u ∉ ℝ³₊ ∪ ℝ³₋ (dynamically trivial);
- (iii) u ∈ ℝ³₋;
- (iv) u ∈ ℝ³₊ and c₁ ≤ 0;
- (v) case 9, u ∈ ℝ³₊, c₁ > 0 and d₁ = 0. Such a network is necessarily 2X → (2+c₁)X, X+Y → (1+d₂)Y, Y → 0 with
  c₁, d₂ ≥ 1 ("generalised LVA type").

*Proof.*

(i) The velocity always lies in the line im Γ, so every trajectory lies in a line parallel to it. A line contains no
Jordan curve.

(ii) im Γᵀ = u^⊥. If uᵢ = 0, then eᵢ ∈ u^⊥. If uᵢ > 0 > uⱼ, then w = |uⱼ|eᵢ + uᵢeⱼ ∈ u^⊥. In both cases there is a
nonzero w ≥ 0 with w = Γᵀλ for some λ. Then d/dt(λ·(x, y)) = Σ wᵢκᵢmᵢ > 0 on ℝ²₊, so λ·(x, y) is strictly increasing
along orbits in ℝ²₊. Hence there is no periodic orbit in ℝ²₊, and by (F1) none in ℝ²_{≥0}.

(iii) By (F2), (F3) the unique positive equilibrium has det J < 0, so it is a saddle. A periodic orbit in ℝ²₊ would
enclose, inside ℝ²₊, equilibria of total index +1 (F4). The only candidate is the saddle, of index −1. Contradiction.

(iv) A direct computation (S7) gives:
- case 9: div(F/(xy)) = (−d₁κ₁x³ + c₁κ₁x²y − c₃κ₃y²)/(x²y²);
- case 10: div(F/(xy)) = (−d₁κ₁x³ + c₁κ₁x²y − c₃κ₃y − d₃κ₃x)/(x²y²).

Here d₁, c₃ ≥ 0, and in case 10 also d₃ ≥ 0. If c₁ < 0, the term c₁κ₁x²y is negative. If c₁ = 0, then d₁ ≥ 1 since
(c₁, d₁) ≠ 0, and the term −d₁κ₁x³ is negative. So the divergence is negative on the simply connected set ℝ²₊, and
Bendixson–Dulac applies. Then use (F1). This is BBH Lemma 4 (published Lemma 3) in the present setting.

(v) We have d₁ = 0, so u₃ = c₁d₂ > 0 gives d₂ ≥ 1, and u₂ = −c₁d₃ > 0 gives d₃ = −1 (d₃ ≥ −1). Then
u₁ = −c₂ − c₃d₂ > 0 forces c₂ < 0, so c₂ = −1, and then c₃d₂ < 1, so c₃ = 0. The system is
ẋ = c₁κ₁x² − κ₂xy, ẏ = d₂κ₂xy − κ₃y, and div(F/(xy)) = c₁κ₁/y > 0 on ℝ²₊ (S7). This is the argument of BBH §4.2 for
c₁ = 1. ∎

**Lemma 2.2 (the Hopf networks).** Assume u ∈ ℝ³₊, c₁ > 0, d₁ > 0 and A = c₁u₁ − c₃u₃ > 0. Then the network
satisfies BBH's Andronov–Hopf condition.

*Proof.*

Case 9. First, c₂ = −1. Suppose instead c₂ ≥ 0.
- From u₃ = c₁d₂ − c₂d₁ > 0 we get c₁d₂ > c₂d₁ ≥ 0, so d₂ ≥ 1.
- From u₁ = c₂d₃ − c₃d₂ > 0 we get c₂d₃ > c₃d₂ ≥ 0, so c₂ > 0 and d₃ > 0.
- If c₃ = 0, then u₂ = −c₁d₃ < 0, a contradiction.
- If c₃ > 0, then d₁ < c₁d₂/c₂ and d₃ > c₃d₂/c₂, so u₂ = c₃d₁ − c₁d₃ < 0, a contradiction.

So c₂ < 0, i.e. c₂ = −1. Now the conditions:
- u₃ = c₁d₂ + d₁ > 0 is d₂/c₂ = −d₂ < d₁/c₁.
- If c₃ > 0: A = −2c₁c₃d₂ − c₁d₃ − c₃d₁ > 0 is ½(d₃/c₃ + d₁/c₁) < −d₂ = d₂/c₂.
- If c₃ = 0: u₂ = −c₁d₃ > 0 gives d₃ = −1.

The sign and range conditions hold by admissibility.

Case 10. Dynamical nontriviality needs some dᵢ < 0. Since d₁ > 0 and d₃ ≥ 0, this gives d₂ = −1. Then
u₃ = −c₁ − c₂d₁ > 0 forces c₂ = −1. Next u₁ = c₃ − d₃ > 0 gives c₃ > 0. The two inequalities follow as in case 9. ∎

The converse also holds: every network satisfying BBH's condition has u ∈ ℝ³₊, c₁ > 0, d₁ > 0 and A > 0. This
follows from BBH's proof of Theorem 9 and is checked exhaustively in a box by `classify_networks.py` (C2). Theorem A
does not use it.

**Lemma 2.3 (the remaining networks).** Assume u ∈ ℝ³₊, c₁ > 0, d₁ > 0 and A ≤ 0. Then the network belongs to exactly
one of the following families. In each, c₁, c₃ ≥ 1.

| family | reactions | conditions |
|---|---|---|
| **9a** | 2X → (2+c₁)X + d₁Y, X+Y → 0, Y → c₃X + (1+d₃)Y | d₁ > c₁, −1 ≤ d₃ < c₃, d₁/c₁ + d₃/c₃ ≥ 2 |
| **9b** | 2X → (2+c₁)X + d₁Y, X+Y → Y, Y → c₃X | d₁ ≥ 1, c₃d₁ ≥ c₁ |
| **10a** | 2X → (2+c₁)X + d₁Y, X+Y → 0, 0 → c₃X + d₃Y | d₁ > c₁, 0 ≤ d₃ < c₃, d₁/c₁ + d₃/c₃ ≥ 2 |

In every family c₂ = −1 and d₂ ∈ {−1, 0}. The quantity

  k := c₃u₃/(c₁u₁)

satisfies k ≥ 1, and A = c₁u₁(1 − k).

*Proof.*

Case 9. Since c₁u₁ > 0 and A ≤ 0, we get c₃u₃ > 0, so c₃ ≥ 1. The argument of Lemma 2.2 (which only used u > 0,
c₁ > 0, d₁ > 0) gives c₂ = −1. Then:
- u₁ = −d₃ − c₃d₂ > 0, with d₃ ≥ −1 and c₃ ≥ 1, excludes d₂ ≥ 1.
- If d₂ = −1: u₁ = c₃ − d₃ > 0 and u₃ = d₁ − c₁ > 0. The condition A ≤ 0 reads 2c₁c₃ − c₁d₃ − c₃d₁ ≤ 0, i.e.
  d₁/c₁ + d₃/c₃ ≥ 2. This is family 9a.
- If d₂ = 0: u₁ = −d₃ > 0 gives d₃ = −1, and A = c₁ − c₃d₁ ≤ 0. This is family 9b.

Case 10. As in Lemma 2.2, d₂ = −1, c₂ = −1, u₁ = c₃ − d₃ > 0 and u₃ = d₁ − c₁ > 0. Then A ≤ 0 is the stated
inequality. This is family 10a.

Finally, A = c₁u₁ − c₃u₃ = c₁u₁(1 − k), and c₁u₁ > 0, so A ≤ 0 is k ≥ 1. ∎

(The BBH example is 9a with c₁ = c₃ = 1, d₁ = 2, d₃ = 0, and k = 1. In the box of all targets with coefficients at
most 20, `classify_networks.py` finds, in case 9, 29 425 networks in 9a, 6 784 in 9b, 342 of LVA type and 21 850 Hopf
networks; in case 10, 27 691 networks in 10a and 11 999 Hopf networks. Every network is classified.)

**Theorem A follows** from Lemmas 2.1–2.3 and Theorem 3.1 below. Indeed, a network that does not satisfy BBH's
condition is in a class of Lemma 2.1, or has u ∈ ℝ³₊, c₁ > 0, d₁ > 0.

(In case 10, u ∈ ℝ³₊ and c₁ > 0 already force d₁ > 0. Since u > 0 and Σuᵢdᵢ = 0, some dᵢ is negative. As d₁, d₃ ≥ 0,
this gives d₂ = −1, and then u₃ = −c₁ − c₂d₁ > 0 forces d₁ > 0.)

In the second case A ≤ 0 by Lemma 2.2 (contrapositive), and Lemma 2.3 puts the network in 9a, 9b or 10a.

## 3. The families 9a, 9b, 10a

**Theorem 3.1.** Every network of family 9a, 9b or 10a has no periodic orbit in ℝ²_{≥0}, for every κ ∈ ℝ³₊.

### 3.1 Normalised form

By (F3), with c₂ = −1, put

  α₁ = c₁u₁, α₂ = u₂, α₃ = c₃u₃, β₁ = d₁u₁, β₂ = −d₂u₂, β₃ = d₃u₃.

All α's are positive, and β₂ ≥ 0 because d₂ ≤ 0. Since Σcᵢuᵢ = Σdᵢuᵢ = 0, we have α₂ = α₁ + α₃ and
β₂ = β₁ + β₃. Also k = α₃/α₁ ≥ 1. The system is linear in y:

  ẋ = A(x) + B(x)y,  ẏ = ρ(C(x) + D(x)y),

with:

| | A(x) | B(x) | C(x) | D(x) |
|---|---|---|---|---|
| case 9 (9a, 9b) | α₁x² | α₃ − α₂x | β₁x² | β₃ − β₂x |
| case 10 (10a) | α₁x² + α₃ | −α₂x | β₁x² + β₃ | −β₂x |

In 10a, d₂ = −1, so β₂ = u₂ = α₂. All of the following is checked in `certify_lienard.py` S3 with the integer data of
each family:

- (3.1) BC − AD = u₁u₂u₃·x²(1 − x) in case 9, and u₁u₂u₃·x(1 − x²) in case 10, with u₁u₂u₃ > 0.
- (3.2) The interval I is I = (x₀, ∞) with x₀ = α₃/α₂ ∈ (0, 1) in case 9, and I = (0, ∞) in case 10. B < 0 on I,
  and D/B > 0 on I.
  - In case 9: D(x₀) = −u₁u₃ < 0 and D' = −β₂ ≤ 0.
  - In 10a: D/B = 1.
- (3.3) The x-nullcline is the graph of N(x) := −A(x)/B(x):
  N = α₁x²/(α₂x − α₃) in case 9, N = (α₁x² + α₃)/(α₂x) in case 10.
  N(1) = 1. N'(1) = (α₁ − α₃)/α₁ = 1 − k in case 9, and N'(1) = (α₁ − α₃)/α₂ in case 10. In both cases N'(1) ≤ 0.

**Lemma 3.2 (confinement).** Every periodic orbit in ℝ²_{≥0} lies in ℝ²₊ ∩ (I × ℝ).

*Proof.* Periodic orbits lie in ℝ²₊ by (F1). This settles case 10, where I = (0, ∞).

In case 9, let (x(t), y(t)) be a periodic solution, so x, y > 0. Let x attain its minimum m at time t*; then
ẋ(t*) = 0. If m ≤ x₀, then B(m) = α₃ − α₂m ≥ 0 and y(t*) > 0, so

  ẋ(t*) = α₁m² + B(m)y(t*) ≥ α₁m² > 0,

a contradiction. Hence x(t) > x₀ for all t. ∎

### 3.2 Liénard form

**Lemma 3.3 (Liénard transformation).** On the strip I × ℝ, put

  f(x) := N'(x) − ρD(x)/B(x),  F(x) := ∫₁ˣ f,  g(x) := −ρ(BC − AD)(x)/B(x)²,  G(x) := ∫₁ˣ g,

and Y := y − N(x) + F(x). Then

  (1/B(x)) ẋ = Y − F(x),  (1/B(x)) Ẏ = −g(x).

(This is a symbolic identity for arbitrary C¹ functions A, B, C, D; see S2.) The map (x, y) ↦ (x, Y) is a C^∞
diffeomorphism of I × ℝ onto itself. The time change dτ = B(x(t)) dt is strictly monotone, since B < 0 on I. So every
periodic orbit of the mass-action system in I × ℝ is mapped to a periodic orbit of the Liénard system

  (L)  x' = Y − F(x),  Y' = −g(x)  on I × ℝ.

Moreover:
- F(x) = N(x) − 1 − ρψ(x), with ψ(x) := ∫₁ˣ D/B. ψ is strictly increasing on I by (3.2).
- f(1) = N'(1) − ρD(1)/B(1) < 0, by (3.2), (3.3) and k ≥ 1. (S2 also checks that f(1) = −tr J(1,1)/B(1).)
- G = ρG̃ with G̃ independent of ρ. By (3.1), g = ρu₁u₂u₃x²(x − 1)/B² in case 9 and ρu₁u₂u₃x(x² − 1)/B² in case 10,
  so (x − 1)g(x) > 0 for x ∈ I ∖ {1}. Hence (1, 0) is the only equilibrium of (L) in I × ℝ.

### 3.3 A comparison lemma for Liénard systems

**Lemma 3.4.** Let I be an open interval containing x̄. Let F ∈ C²(I) with F(x̄) = 0 and f := F', and let g ∈ C¹(I)
with (x − x̄)g(x) > 0 for x ∈ I ∖ {x̄}. Put G(x) = ∫_{x̄}^x g. Assume:
- (H1) f(x̄) < 0;
- (H2) F(a) > F(b) for all a < x̄ < b in I with G(a) = G(b).

Then x' = Y − F(x), Y' = −g(x) has no periodic orbit in I × ℝ.

Comparison arguments of this kind, in the variable z = G(x) (Filippov's transformation), are classical for Liénard
systems. We give a self-contained proof.

*Proof.* Suppose Γ is a periodic orbit in I × ℝ, with minimal period T. Then Γ is a Jordan curve. It does not contain
the only equilibrium (x̄, 0).

*Step 1: Γ meets the line ℓ = {x = x̄} in exactly two points.*
- At a point of Γ ∩ ℓ we have x' = Y ≠ 0. So the zeros of x(t) − x̄ are simple, and there are finitely many in one
  period.
- There is at least one zero. If x(t) > x̄ for all t, then Y' = −g(x) < 0 throughout, so Y is strictly decreasing,
  which is impossible on a periodic orbit. The case x(t) < x̄ is the same.
- The zeros cut Γ into arcs that lie alternately in {x > x̄} ("right arcs") and {x < x̄} ("left arcs").
- On a right arc, Y is strictly decreasing, from a value Y > 0 at the crossing where the orbit enters to a value Y < 0
  where it leaves. So the arc is the graph x = X(Y) of a continuous function on an interval [β⁻, β⁺] ∋ 0, with
  X > x̄ inside and X = x̄ at the ends. X is C¹ in the interior, because Y' ≠ 0 there. Left arcs are similar.
- Distinct right arcs are disjoint, since Γ is simple and they do not share endpoints.
- Two disjoint right arcs over [α, β] and [α', β'] are strictly nested. Say α < α'. If also β < β', then on [α', β]
  the difference X' − X is negative at α' (X'(α') = x̄ < X(α')) and positive at β. So the arcs would meet,
  contradicting disjointness. Hence [α', β'] ⊂ (α, β). Left arcs are nested in the same way.
- Suppose there were more than one right arc. Let R* be the outermost right arc, from P = (x̄, Y_top) down to
  Q = (x̄, Y_bot). The left arc L₁ starting at Q ends at some P₁ ≠ P. P₁ is the top of an inner right arc, so it lies
  strictly between Q and P on ℓ. The left arc L' ending at P starts at some Q' ≠ Q, the bottom of an inner right arc,
  so Y_bot < Y(Q') < 0. Then L₁ lies over [Y_bot, Y(P₁)] and L' over [Y(Q'), Y_top], two intervals that are not
  nested. This contradicts the nesting property.
- Hence Γ = R ∪ L, with one right arc R and one left arc L, both graphs over the same interval [Y_bot, Y_top],
  Y_bot < 0 < Y_top.

*Step 2: comparison.* Write the arcs as x = X_R(Y) and x = X_L(Y), and put Z_R = G(X_R) and Z_L = G(X_L). On the
open interval, dX/dY = x'/Y' = (F(X) − Y)/g(X). Hence

  Z_R' = F(X_R) − Y,  Z_L' = F(X_L) − Y,  and  D(Y) := Z_L − Z_R satisfies D' = F(X_L) − F(X_R).

D is continuous on [Y_bot, Y_top], C¹ inside, and vanishes at both ends.

By (H1), F is strictly decreasing near x̄, with F(x̄) = 0. As Y tends to either end of the interval, X_L, X_R → x̄,
so F(X_L) > 0 > F(X_R) near each end. Hence D' > 0 near both ends. Thus D > 0 on some interval (Y_bot, Y_bot + ε)
and D < 0 on some interval (Y_top − ε, Y_top).

Let Y* be the smallest Y ≥ Y_bot + ε/2 with D(Y) ≤ 0. Then Y* lies in the open interval, D(Y*) = 0, and D > 0 just
to the left of Y*, so D'(Y*) ≤ 0. But D(Y*) = 0 means G(X_L(Y*)) = G(X_R(Y*)) with X_L(Y*) < x̄ < X_R(Y*). By (H2),
D'(Y*) = F(X_L(Y*)) − F(X_R(Y*)) > 0. This is a contradiction. ∎

### 3.4 The key inequalities

For the systems of §3.2, with x̄ = 1:

  F(a) − F(b) = [N(a) − N(b)] + ρ[ψ(b) − ψ(a)],

and the second bracket is positive for a < b. Moreover G(a) = G(b) ⟺ G̃(a) = G̃(b), a condition independent of ρ.
So (H2) holds **for every ρ > 0** as soon as:

  (★)  N(a) ≥ N(b) whenever a < 1 < b in I and G̃(a) = G̃(b).

**Lemma 3.5 (case 9).** For families 9a and 9b, (★) holds, with strict inequality.

*Proof.* Substitute x = (η + k)/(1 + k). This is an increasing affine bijection from η ∈ (0, ∞) onto I = (x₀, ∞),
with η = 1 ↔ x = 1. Using α₂ = α₁(1 + k) and α₃ = α₁k (S4):

  N = (η + k)²/((1 + k)²η),  and  (g/ρ) dx = K·(η + k)²(η − 1)/η² dη  with K > 0.

Hence G̃ = K·H_k(η), where

  H_k(η) := ∫₁^η (σ + k)²(σ − 1)/σ² dσ = η²/2 + (2k−1)η + (k² − 2k) log η + k²/η − (k² + 2k − ½).

H_k is strictly decreasing on (0, 1], strictly increasing on [1, ∞), and H_k(1) = 0.

*Claim.* For k ≥ 1 and 0 < a < 1: H_k(k²/a) > H_k(a).

Let P(a) := H_k(k²/a) − H_k(a). Differentiating (S5),

  P'(a) = −Q(a)/a³,  Q(a)/a² = 8k(k − 1) + (6k − 2)s + s²,  where s := (a − k)²/a.

For k ≥ 1 and a ∈ (0, 1) we have a ≠ k, so s > 0 and Q > 0. Hence P is strictly decreasing on (0, 1]. Also
P(1) = H_k(k²) ≥ 0, because k² ≥ 1 and H_k ≥ 0 on [1, ∞). So P(a) > 0 for a ∈ (0, 1).

*Conclusion.* Let η_a < 1 < η_b with H_k(η_a) = H_k(η_b). Then H_k(η_b) < H_k(k²/η_a). Both arguments are greater
than 1, and H_k is increasing there, so η_b < k²/η_a, i.e. η_aη_b < k². Finally (S5),

  N(a) − N(b) = (η_a − η_b)(1 − k²/(η_aη_b))/(1 + k)² > 0,

since both factors are negative. ∎

**Lemma 3.6 (case 10).** For family 10a, (★) holds, with strict inequality.

*Proof.* By (3.1), g/ρ = (u₁u₃/u₂)·h'(x), with h(x) = x²/2 − log x (S6). So G̃(a) = G̃(b) ⟺ h(a) = h(b).

The function t ↦ h(1/t) − h(t) vanishes at t = 1 and has derivative −(1 − t²)²/t³ < 0. So h(1/a) > h(a) = h(b) for
0 < a < 1. Since h is increasing on (1, ∞), this gives b < 1/a, i.e. ab < 1 ≤ k.

Then N(a) − N(b) = (a − b)(α₁ − α₃/(ab))/α₂ > 0, because a < b and α₃/(ab) > α₃ = kα₁ ≥ α₁. ∎

### 3.5 Proof of Theorem 3.1

Fix a network of 9a, 9b or 10a and κ ∈ ℝ³₊, and pass to the normalised system (F3) with ρ = x̄/ȳ. Suppose there were
a periodic orbit in ℝ²_{≥0}.
- By Lemma 3.2 it lies in I × ℝ.
- By Lemma 3.3 it is mapped to a periodic orbit of (L) in I × ℝ.
- The hypotheses of Lemma 3.4 hold with x̄ = 1: the sign condition on g and (H1) by Lemma 3.3, and (H2) by
  Lemmas 3.5 and 3.6 together with F(a) − F(b) = [N(a) − N(b)] + ρ[ψ(b) − ψ(a)].

So no such orbit exists. This contradiction proves Theorem 3.1, and with §2 also Theorem A. ∎

## 4. Remarks

1. **Relation to result 2.** For the subfamily of result 2 (9a with d₃ = 0), k = d₁/c₁ − 1. So k ≥ 1 is exactly
   d₁ ≥ 2c₁, the condition of result 2, which was proved there with the Dulac function y^{−b}. The present proof is
   independent of it.

   Monomial Dulac functions x^α y^β do not seem to suffice in general. This is a side computation, done by hand and
   not scripted. In x^{1−α}y^{1−β}div(x^α y^β F), the cubic and the lower-order homogeneous parts must both be ≤ 0 on
   the quadrant. Comparing them, a suitable (α, β) exists for all rate constants only when d₃ ≥ 0 and d₁ ≥ 2c₁. In
   9a with d₃ = −1, in 9b, and in 9a/10a with d₁ < 2c₁, it fails for large κ₁/κ₂. This is why a different tool is
   needed.

2. **Sharpness in case 10, and a stronger statement there.** Lemma 3.6 shows ab < 1 at every pair, whatever k is.
   For any case-10 network with u ∈ ℝ³₊ and c₁ > 0, including the Hopf networks, which have k < 1, the structure of
   §3.1 holds. Then

   F(a) − F(b) = (b − a)[ρ − (α₁ − α₃/(ab))/α₂],  and  (α₁ − α₃/(ab))/α₂ < (α₁ − α₃)/α₂ = N'(1).

   Also tr J(1, 1) < 0 ⟺ f(1) < 0 ⟺ ρ > N'(1), since ψ' = 1. So Lemma 3.4 applies whenever tr J < 0. Consequently:
   **every case-10 network with u ∈ ℝ³₊, c₁ > 0 has no periodic orbit whenever the positive equilibrium satisfies
   tr J < 0.** For the Hopf networks, periodic orbits can therefore exist only when the equilibrium is repelling or has
   tr J = 0; Lemma 3.4 does not cover tr J = 0. The computation is the same as in Lemma 3.6. It is not needed for
   Theorem A, and is not separately scripted beyond S6.

3. **Case-9 Hopf networks (numerical observation only).** For case-9 Hopf networks with c₃ > 0, we computed
   numerically the supremum over G-paired points of (N(b) − N(a))/(ψ(b) − ψ(a)) for six networks. In every one it is
   attained in the limit a, b → 1, where it equals the Hopf value ρ_H. This suggests the same conclusion, namely no
   periodic orbit whenever tr J < 0, but **we have not proved it**. It is not part of any claim.

4. **What is not claimed.**
   - The number of limit cycles of the Hopf networks: BBH's supercriticality gives at least one near the bifurcation.
   - Networks with more species, more reactions, or other source configurations.
   - The tr J = 0 case of Remark 2.

## 5. Verification map

| Claim | Where checked |
|---|---|
| Normalisation, uniqueness of equilibrium, det J and tr J at (1,1) (F2, F3) | `certify_lienard.py` S1 |
| Liénard identity (Lemma 3.3), f(1) = −tr J/B(1) | S2 |
| u > 0, (3.1), (3.2), (3.3), A = c₁u₁ − c₃u₃ = α₁ − α₃ = BBH inequality, f(1) < 0 in each family | S3 |
| η-substitution, N(η), G̃ = K·H_k with K > 0 | S4 |
| H_k formula, P' = −Q/a³, Q/a² = 8k(k−1) + (6k−2)s + s², P(1) = H_k(k²), N(a) − N(b) identity | S5 |
| Case 10: g = ρ(u₁u₃/u₂)h', the h-inequality, N(a) − N(b) identity | S6 |
| Dulac 1/(xy) for classes (iv) and (v), boundary behaviour | S7 |
| Classification (Lemmas 2.1(v), 2.2, 2.3, and the converse of 2.2) in the box T = 20, about 1.7·10⁸ networks | `classify_networks.py` (C1–C5) |
| Numerical sanity: no cycles in 128 backward integrations; Hopf controls do show cycles; (H2) by quadrature | `numerical_sanity.py` (not a proof) |

Everything in §§1–3 is either proved in the text or an exact symbolic identity checked by `certify_lienard.py`. The
finite box enumeration and the floating-point runs are cross-checks only; the proof does not use them.

## 6. Assumptions and non-rigorous steps

- The proof of Theorem A has no numerical or unverified step. It uses three standard facts, (F4): Poincaré index
  theory, the Jordan curve theorem, and the Bendixson–Dulac criterion.
- It also uses BBH's definition of "admits an Andronov–Hopf bifurcation" through the explicit conditions of their
  Theorem 9. Theorem A is stated for networks failing those conditions.
- Corollary B's direction (c) ⟹ (b) relies on BBH's computation of the first focal value, which we did not re-check.
