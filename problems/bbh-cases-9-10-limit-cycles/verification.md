Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Independent verification: Banaji–Boros–Hofbauer, cases 9 and 10 without Andronov–Hopf bifurcation

Date: 25 October 2026.

**Inputs used.**
- The claimed result, as stated in the task: "no periodic orbit for any network of case 9 or 10 that admits no
  Andronov–Hopf bifurcation".
- The deliverables in this folder: `README.md`, `result.md`, `paper.tex`, `references.*`, `scripts/` and
  `certificates/`.
- `PROTOCOL.md`, and the verification report of result 2 as a model of the format.
- The primary source, downloaded by me:
  - arXiv:2304.02303v1 (25 pages; Theorem 9 and its proof read on rendered page images, pp. 14–16);
  - the published version, from the Europe PMC full-text XML of PMC10953466.

I did not read any solver log, progress file or working notes. I derived the classification (§2 below) and wrote
`verification/classify_independent.py` before opening `result.md`, `paper.tex` or the solver's scripts.

## Verdict: **CONFIRMED WITH CORRECTIONS**

Theorem A is correct and complete. Take any network of case 9 (sources 2X, X+Y, Y) or case 10 (sources 2X, X+Y, 0),
with arbitrary targets, that fails the Andronov–Hopf condition of BBH Theorem 9 (published Theorem 3). For every
κ ∈ ℝ³₊ its mass-action system has no periodic orbit in the closed quadrant, and so no limit cycle. This answers the
question of BBH §6 in full, negatively, and contains the example of §6 (result 2) as a special case.

I checked every step of the proof by hand and re-derived its identities with my own code. I found no mathematical
error and no gap. Numerically I searched for counterexamples (16 095 return-map integrations, including the
boundary k = 1 of the Hopf region and continuation from Hopf networks) and found none. In five Hopf controls the same
detector does find the stable cycle.

All corrections are editorial; the mathematics was not changed. The one factual slip was in the README summary
table, which omitted the inequalities u > 0 that define the families. Details are in §6.

## 1. Source and statement

- **Quote.** The §6 sentence quoted in `README.md`, `result.md` and `paper.tex` matches arXiv v1, p. 22 (text
  extracted with pdftotext). In the published version (§6) the sentence is identical; only "Theorem 9" becomes
  "Theorem 3".
- **Numbering.** I confirmed the map stated in `result.md` (arXiv Lemmas 2, 3, 4 = published Lemmas 1, 2, 3;
  Theorem 9 = Theorem 3) from sentences that match between the two versions. For example, "by Lemma 2 [resp. 1],
  each corresponding mass-action system has a unique positive equilibrium", and "by Lemma 4 [resp. 3], ... of the
  Lotka–Volterra type with no diagonal term".
- **Definitions.** These were taken from BBH §4.1 and §4.5.
  - Reaction i is aᵢX + bᵢY → (aᵢ+cᵢ)X + (bᵢ+dᵢ)Y, with integer cᵢ, dᵢ and nonnegative targets.
  - The cases are defined up to exchange of X and Y, with the sources listed positively oriented.
  - "Admits an Andronov–Hopf bifurcation" is defined by the explicit conditions of Theorem 9. I transcribed them from
    the page image of p. 15, and they agree with `result.md` §1 and `paper.tex` §1.
  - The set the question is about is therefore all case-9/10 networks that fail those conditions.
- **Scope.** The question is about limit cycles. Theorem A proves more: there is no periodic orbit at all.

## 2. Independent re-derivation of the classification (before reading the solver's work)

**Trace criterion.** From (14) and (18) of BBH, at the unique positive equilibrium (when u = c × d ∈ ℝ³₊):

  tr J = μ(P/x̄ + Q/ȳ),  with P = Σ aᵢcᵢuᵢ and Q = Σ bᵢdᵢuᵢ.

- Case 9 gives P = c₁u₁ − c₃u₃ =: A and Q = −d₁u₁. Case 10 gives P = A and Q = d₂u₂.
- Every positive point is the equilibrium for suitable κ, so x̄/ȳ ranges over (0, ∞).
- Hence tr J changes sign as κ varies (with det J > 0) iff P·Q < 0.
- C1 of `classify_independent.py` checks the trace and determinant formulas symbolically. C2 checks that P·Q < 0
  coincides exactly with BBH's literal Theorem 9 conditions, in every network of both boxes below.

**Non-Hopf networks with u > 0, derived by hand.** Networks with u ∉ ℝ³₊ are dynamically trivial or have a saddle.

| verifier's family | definition | sign of tr J for all κ | solver's name |
|---|---|---|---|
| F1 (case 9) | c₁ ≥ 1, c₂ = −1, c₃ = 0, d₁ = 0, d₂ ≥ 1, d₃ = −1 (generalised LVA) | > 0 | LVA, Lemma 2.1(v) |
| F2a (case 9) | c₁ ≥ 1, c₂ = −1, c₃ ≥ 1, d₁ ≥ 1, d₂ = 0, d₃ = −1, c₁ ≤ c₃d₁ | < 0 | 9b |
| F2b (case 9) | c₁ ≥ 1, c₂ = d₂ = −1, c₃ ≥ 1, d₁ > c₁, −1 ≤ d₃ < c₃, c₃d₁ + c₁d₃ ≥ 2c₁c₃ | < 0 | 9a |
| F3 (case 9) | c₁ ∈ {−2, −1, 0}, c₂ = −1, u > 0 | < 0 | "c₁ ≤ 0", Lemma 2.1(iv) |
| F4 (case 9) | c₁ ∈ {−2, −1}, c₂ ≥ 0, d₂ = −1, u > 0 | < 0 | "c₁ ≤ 0" |
| G1 (case 10) | c₁ ≥ 1, c₂ = d₂ = −1, d₁/c₁ + d₃/c₃ ≥ 2, u > 0 | < 0 | 10a |
| G2 (case 10) | c₁ ≤ 0, c₂ = d₂ = −1, u > 0 | < 0 | "c₁ ≤ 0" |
| G3 (case 10) | c₁ < 0, c₂ ≥ 0, d₂ = −1, u > 0 | < 0 | "c₁ ≤ 0" |

Two consequences of the derivation:
- **The question is non-trivial.** Apart from F1, every non-Hopf network with u > 0 has a stable equilibrium for
  every κ. So "no Hopf" excludes periodic orbits only if no pair of cycles (or semistable cycle) can surround a
  stable focus. This is exactly what the solver's Liénard argument proves for 9a, 9b and 10a.
- **Theorem 9 contains no degenerate case.** The case P = Q = 0 (tr J ≡ 0) never occurs (C3).

`verification/classify_independent.py` checks the families against the trace signs (C4, C5) in two boxes:
- c_i, d_i ≤ 8: 1.7·10⁶ networks, 8 s;
- the solver's box (all targets ≤ 20): 1.7·10⁸ networks, 216 s.

**Comparison with `certificates/classify_networks.log`.** In the box T = 20 all counts agree exactly.

| class | case 9 | case 10 |
|---|---|---|
| dynamically trivial | 84 985 037 | 85 034 404 |
| rank < 2 | 26 475 | 27 016 |
| saddle | 65 560 | 35 910 |
| Hopf | 21 850 | 11 999 |
| LVA = F1 | 342 | – |
| 9b = F2a | 6 784 | – |
| 9a = F2b | 29 425 | – |
| 10a = G1 | – | 27 691 |
| "c₁ ≤ 0" | 48 527 = F3 15 372 + F4 33 155 | 46 980 = G2 12 980 + G3 34 000 |

## 3. Proof check, step by step

**Lemma 2.1 (classes without periodic orbits).**
- (F1), boundary: correct. A periodic orbit meeting the open quadrant lies in it by forward invariance. The boundary
  of the quadrant is homeomorphic to ℝ and contains no Jordan curve.
- (i), rank Γ ≤ 1: trajectories lie on lines.
- (ii), dynamically trivial: the Stiemke-type argument is correct. When rank Γ = 2, im Γᵀ = u^⊥, and the vector w
  constructed in the proof is ≥ 0 and nonzero, which gives a strictly increasing linear function.
- (iii), u < 0: the unique positive equilibrium has det J < 0, so it is a saddle. By convexity the interior of a
  periodic orbit stays in ℝ²₊, so index theory gives a contradiction.
- (iv), u > 0 and c₁ ≤ 0: I recomputed x²y² div(F/(xy)) (V2). It is c₁κ₁x²y − c₃κ₃y² − d₁κ₁x³ in case 9, plus the
  term −d₃κ₃x in case 10 (where c₃κ₃y replaces c₃κ₃y²). It is < 0 on ℝ²₊ when c₁ < 0, and also when c₁ = 0, since
  then d₁ ≥ 1. This disposes of F3, F4, G2 and G3 at once, whatever the sign of A.
- (v), generalised LVA: the divergence is c₁κ₁/y > 0 (V2).

**Lemmas 2.2 and 2.3.**
- Every implication was re-derived by hand. This includes c₂ = −1 whenever u > 0, c₁ > 0 and d₁ > 0, and d₂ ≥ 1
  being excluded in the families.
- The case analysis is exhaustive. Rank and u-sign classes come first; then u > 0 with c₁ ≤ 0; then, in case 9,
  c₁ > 0 with d₁ = 0. Case 10 with u > 0 and c₁ > 0 forces d₁ > 0.
- What remains, u > 0, c₁ > 0 and d₁ > 0, splits into A > 0 (the BBH Hopf networks, Lemma 2.2) and A ≤ 0
  (9a, 9b, 10a).
- k = c₃u₃/(c₁u₁) ≥ 1 is equivalent to A ≤ 0.

**Normalisation (F2, F3).**
- The equilibrium formulas were re-derived by hand.
- The scaling x = x̄X, y = ȳY, t = (x̄/μ)T gives X' = Σcᵢuᵢmᵢ and Y' = ρΣdᵢuᵢmᵢ, with ρ = x̄/ȳ (V1). The single
  parameter ρ sweeps (0, ∞).

**§3.1 structure.**
- With c₂ = −1, I re-derived α₂ = α₁ + α₃ and β₂ = β₁ + β₃ (V3).
- BC − AD = u₁u₂u₃x²(1 − x) in case 9 holds for every d₂, not only d₂ ∈ {−1, 0}, and BC − AD = u₁u₂u₃x(1 − x²) in
  case 10 (V3).
- D(x₀) = −u₁u₃ and D' = −β₂ ≤ 0 (because d₂ ≤ 0), hence D/B > 0 on I.
- N(1) = 1, and N'(1) = 1 − k (case 9) or (α₁ − α₃)/α₂ (case 10).

**Lemma 3.2 (confinement to I × ℝ).** Correct: at the minimum m ≤ x₀ of x one would have ẋ ≥ α₁m² > 0.

**Lemma 3.3 (Liénard form).**
- I checked the identity on the actual rational functions A, B, C, D of each family (V4), not only on abstract
  functions. With W = y − N(x) and f = N' − ρD/B:
  - ẋ/B = W;
  - d/dt(W + F(x))/B = −g, where g = −ρ(BC − AD)/B².
- The shear is a diffeomorphism of I × ℝ. The time change dτ = B dt reverses orientation (B < 0) but preserves
  periodic orbits.
- f(1) = −tr J(1,1)/B(1) (V4), and f(1) < 0 because N'(1) ≤ 0 and D(1)/B(1) > 0.

**Lemma 3.4 (comparison lemma).** I read the proof in `result.md` line by line. It is self-contained, cites no
external theorem, and is correct.
- *Step 1.*
  - Crossings of x = x̄ are transversal, because x' = Y ≠ 0 there: the only equilibrium is (x̄, 0).
  - On each side Y is strictly monotone, so every arc is a graph over a Y-interval containing 0.
  - Disjoint arcs on one side are strictly nested, by the intermediate value theorem.
  - The argument that there is only one right arc (the two left arcs at the outermost right arc would be distinct
    and not nested) is complete.
- *Step 2.*
  - With Z = G(X), Z' = F(X) − Y holds on the open interval, where Y' ≠ 0. Then Δ = Z_L − Z_R vanishes at both ends.
  - (H1) gives Δ > 0 near the bottom and Δ < 0 near the top.
  - At the first zero Y* of Δ (well defined: Δ > 0 on a neighbourhood of the bottom end), Δ'(Y*) ≤ 0. This
    contradicts (H2).
- Both C¹ hypotheses hold on I, because B ≠ 0 there.

**Reduction of (H2) to (★).**
- F(a) − F(b) = [N(a) − N(b)] + ρ[ψ(b) − ψ(a)], and the second bracket is > 0 for a < b because D/B > 0.
- G = ρG̃, so the pairing condition does not depend on ρ.
- Hence (★) with strict inequality gives (H2) for every ρ > 0. This is the key structural observation, and it is
  correct.

**Lemma 3.5 (case 9).**
- I re-derived the following by hand and in V5–V6:
  - x = (η + k)/(1 + k) maps (0, ∞) onto (x₀, ∞), with x₀ = k/(1 + k);
  - N = (η + k)²/((1 + k)²η);
  - (g/ρ)dx = K(η + k)²(η − 1)/η² dη, with K = u₁u₂u₃/(α₁²(1 + k)⁴) > 0;
  - the primitive H_k.
- For P(a) = H_k(k²/a) − H_k(a) I used a different factorisation from the solver's:

    P'(a) = −(a + k)²((a − 1)² + k² − 1)/a³.

  It is < 0 on (0, 1) for every k ≥ 1. With P(1) = H_k(k²) ≥ 0 this gives P > 0 on (0, 1), as claimed.
- The step from η_aη_b < k² to N(a) > N(b) is correct (V6).

**Lemma 3.6 (case 10).**
- g/ρ = (u₁u₃/u₂)h' with h = x²/2 − log x.
- h(1/t) − h(t) is decreasing, so ab < 1 at every G-paired pair.
- N(a) − N(b) = (a − b)(α₁ − α₃/(ab))/α₂ > 0 when k ≥ 1 (V7).

**Remarks.**
- *Remark 1, first paragraph:* correct (k = d₁/c₁ − 1 when d₃ = 0).
- *Remark 1, monomial Dulac side claim.* It is not used anywhere. I re-derived it for case 9 by hand:
  - in 9a, a monomial Dulac function x^α y^β works for all κ iff d₃ ≥ 0 and d₁ ≥ 2c₁;
  - with d₃ = −1, or in 9b, the x²y coefficient grows like κ₁ while the admissible discriminant only grows like
    √κ₁.
  - I did not check the case-10 part. A note to this effect was added.
- *Remark 2 (case 10: no periodic orbit whenever tr J < 0, Hopf networks included):* correct.
  - F(a) − F(b) = (b − a)[ρ − (α₁ − α₃/(ab))/α₂] and tr J(1,1) = α₂(N'(1) − ρ) (V7);
  - together with ab < 1, (H1) and (H2) hold exactly when tr J < 0;
  - the remark itself says that tr J = 0 is not covered.
- *Remark 3:* explicitly numerical and not claimed; I did not check it.

**Corollary B.** Correctly conditional on BBH's first-focal-value computation, which neither the solver nor I
re-checked. Theorem A does not depend on it.

## 4. Scripts re-run (WSL Fedora 44, Python 3.12.15, sympy 1.14.0, numpy 2.5.4, scipy 1.18.1)

| script | result | wall time | comparison |
|---|---|---|---|
| `scripts/certify_lienard.py` | 64 PASS, exit 0 | 3.8 s | output identical to `certificates/certify_lienard.log` |
| `scripts/classify_networks.py` | 9 PASS, exit 0 | 43 s | identical except the runtime line |
| `scripts/numerical_sanity.py` | 5 PASS, exit 0 | 175 s | identical except last floating digits (~1e-13) of three control lines |
| `verification/classify_independent.py 8` | 12 PASS | 7.6 s | – |
| `verification/classify_independent.py T20` | 12 PASS | 216 s | all class counts equal to the solver's |
| `verification/verify_symbolic.py` | 54 PASS | 14 s | – |
| `verification/verify_numeric.py` | 4 PASS | 265 s | – |

**What the solver's scripts actually check.**
- `certify_lienard.py` checks the exact identities and the sign facts of each family. For 9a it uses the symbolic
  parametrisation d₁ = c₁ + p, d₃ = c₃ − q with p, q > 0.
- The classification itself is proved in the text. `classify_networks.py` only cross-checks it in a finite box,
  which `result.md` states correctly.
- `numerical_sanity.py` is floating point and is labelled as such.
- I found nothing printed as PASS that is not actually tested.

## 5. Numerical search for counterexamples (`verification/verify_numeric.py`; can refute, cannot prove)

The search works on the normalised system, on the section S = {x = 1, y > 1}.
- S is transversal: there ẋ = α₁(1 − y) in case 9 and α₂(1 − y) in case 10.
- Every periodic orbit around (1, 1) crosses S.
- I computed the displacement δ(r) = y_return − 1 − r of the first-return map (LSODA, rtol 1e-10, orbit-preserving
  time change). A periodic orbit requires δ = 0 somewhere.

**R1: 29 non-Hopf family members.**
- The members:
  - 8 random members each of 9a and 10a, and 6 of 9b, with c₁, c₃ ≤ 6 and d₁ ≤ 12 (these 22 include one 9a and one
    9b member with k = 1);
  - the BBH example;
  - 9a with k = 1 and d₃ = −1, and 9a with k = 1.1;
  - two 9b members with k = 1;
  - 10a with k = 1, and 10a with k = 1.1.
- The grid: ρ ∈ [10⁻⁴, 10³] (15 values) and r ∈ [10⁻³, 10³] (37 values), so 16 095 integrations.
- Outcomes:
  - 4 320 returns, all with δ < 0;
  - 11 338 orbits entered a 10⁻⁶-neighbourhood of the attracting equilibrium without returning;
  - 437 (2.7%) were inconclusive within the time limit, mostly in 10a at the smallest ρ and largest r.
- The largest δ/r was −0.107, attained at the boundary k = 1. There was no sign change and no near-zero
  displacement.

**R2: positive controls.** Five Hopf networks (three of case 9, two of case 10) at ρ = 0.97ρ_H, where tr J > 0. The
displacement changes sign in all five: the stable cycle is detected.

**R3: continuation across the Hopf boundary in network space.**
- I used six pairs, each of a Hopf network (k < 1) and a non-Hopf neighbour (k ≥ 1) differing in one coefficient.
  For example (c₁, c₃) = (10, 1) with d₁ = 19 → 20, 21, in both cases.
- ρ ran over [0.5, 2]·ρ_H of the Hopf twin, plus small values of ρ.
- There was no cycle in any non-Hopf neighbour; the largest δ/r was −0.146.

**R4: hypothesis (H2) by direct quadrature.**
- This check uses mpmath at 30 digits and does not use the closed forms of the paper.
- At 348 G-paired points (a < 1 < b), over 29 members, ρ ∈ {0.001, 0.3, 5} and four values of a, F(a) > F(b)
  always.

## 6. Issues found and corrections applied

| # | issue | severity | correction |
|---|---|---|---|
| 1 | `README.md` summary table of families listed only "d₁/c₁ + d₃/c₃ ≥ 2" (9a, 10a). It omitted d₁ > c₁ and d₃ < c₃, i.e. u > 0, so as printed it included networks that are not in the families. `result.md` and `paper.tex` were correct. | minor (summary only) | Conditions completed, with a note that they mean u > 0 |
| 2 | `paper.tex` abstract and introduction said that only cases 9 and 10 admit an Andronov–Hopf bifurcation. BBH case 8 admits a vertical (degenerate) one. | minor (wording) | Rephrased: apart from the vertical bifurcation of case 8, only cases 9 and 10 have such networks, and there the bifurcation is supercritical |
| 3 | `paper.tex` used Γ without defining it, and cited BBH Lemma 4 / Theorem 9 without saying that the numbering is that of arXiv v1. | minor | Γ defined. Numbering note added (arXiv Lemmas 2–4 = published Lemmas 1–3, Theorem 9 = Theorem 3) |
| 4 | `paper.tex`, case 10: D = −α₂x was written without saying that β₂ = α₂, and "some dᵢ < 0 is needed" was not justified. The positivity of ψ(b) − ψ(a) was not tied to D/B > 0. Theorem A did not mention the X↔Y exchange. | minor (rigour of exposition) | One-line justifications added |
| 5 | Status and credits. README and paper said "preliminary" and "effort not recorded" for the solver (the owner confirmed effort xhigh). The paper had never been compiled. | editorial | Status set to solved (confirmed with corrections). Solver effort set to xhigh in README and paper (`result.md` already said xhigh). Verifier credit added. `paper.pdf` compiled (6 pages). Verification paragraph added to the paper's §4 |
| 6 | `result.md` Remark 1: side claim about monomial Dulac functions, not scripted. | observation | Verifier's note added: case 9 re-derived, case 10 not checked, not used in the proof |
| 7 | `README.md`: the case-10 sharpness sentence was slightly loose ("whenever tr J < 0" without saying it is per κ, and without the c₁ ≤ 0 networks). | minor | Made precise |

No change was made to any script, certificate or log of the solver.

## 7. Novelty (re-checked independently; details appended to `references.md`)

- **Citing works.**
  - OpenAlex lists 10 records and Semantic Scholar 12, the newest from 5 Oct 2026.
  - I grepped 2406.13451, 2409.06877 and 2610.06130 in full text. Only 2406.13451 mentions cases 9 and 10, and only
    to place its Hopf networks.
- **Authors.** No new paper by Banaji, Boros or Hofbauer after Banaji 2609.07265. `balazsboros/reaction_networks`
  is unchanged since 30 Sep 2026.
- **AI campaigns.**
  - Evidence Press (`ipitchford`, 193 repositories): all 149 READMEs pushed since 1 Sep 2026 were grepped for the
    BBH keywords, with no hit.
  - Palomar (`mmislan`, 27 repositories): all READMEs grepped, with no hit. `palomar-lienard-limit-cycles` concerns
    classical Liénard lower bounds and is unrelated.
- **Zenodo.** Required-term queries returned nothing relevant.
- **GitHub code search** for the DOI and the arXiv id found only this repository and the authors' notebook
  repository.
- **Web search** found nothing either.
- Gemini was not used.

**Conclusion: no prior answer was found. The result appears to be new.**

## 8. What I only read, without recomputing

- BBH's supercriticality of the Hopf bifurcation in cases 9 and 10 (first focal value, Mathematica notebook). It is
  used only in Corollary B (c) ⟹ (b), which is stated as conditional on it. The positive controls of R2 agree with
  it.
- BBH Lemma 2 (published Lemma 1, uniqueness of the positive equilibrium). It is not needed, because the solver gives
  explicit formulas, which I re-derived.
- The standard planar facts: Jordan curve theorem, Poincaré index, Bendixson–Dulac.

## 9. Reservations

- None on the correctness of Theorem A.
- The novelty search cannot exclude unindexed or very recent work.
- Corollary B (non-Hopf ⟺ no limit cycle) inherits BBH's unpublished-in-detail focal-value computation for the
  direction "Hopf ⟹ limit cycle". Theorem A, which is the answer to the question asked, does not.

## Files written or changed by the verifier

- `verification.md` (this report).
- `verification/classify_independent.py` with `classify_independent.log` (box ≤ 8) and
  `classify_independent_T20.log` (solver's box). This is my own classification.
- `verification/verify_symbolic.py` and `verify_symbolic.log`: exact sympy checks V1–V8.
- `verification/verify_numeric.py` and `verify_numeric.log`: numerical search R1–R4, not a proof.
- `README.md`, `result.md`, `paper.tex` and `references.md`: the corrections listed in §6, and the novelty re-check.
- `paper.pdf`: compiled with `templates/build_tex.py`.
