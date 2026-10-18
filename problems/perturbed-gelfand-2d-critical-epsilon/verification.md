# Independent verification: perturbed Gelfand equation on the disk, critical parameter

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

Files read: `paper.tex`/`paper.pdf`, `README.md`, `references.md`, `references.bib`, `scripts/`, `certificates/`, `PROTOCOL.md`, KLO (arXiv:1612.07722v1 and the author-hosted manuscript of the JDE version), and the full text of Du–Lou (2001), supplied by the repository owner during the verification. I did not open `PROGRESS.md`, `FINAL-REPORT.md`, the `BRIEF-*` files or `tmp/`.

## Verdict

**CONFIRMED WITH CORRECTIONS.**

The global theorem (Theorem 2 in `paper.pdf`, label `thm:global`; the template shares the theorem and equation counters) answers the conjecture as stated by Korman–Li–Ouyang. There is a unique ε₀ with |ε₀ − 0.24210616559524695| < 10⁻¹²:
- every curve with 0 < ε < ε₀ is exactly S-shaped: two simple folds with signs +,−,+, exact multiplicity 1/2/3/2/1, and ordered solutions;
- every curve with ε ≥ ε₀ is strictly increasing;
- at ε₀ there is a single ordinary cusp.

I found no mathematical error and no gap in the logical chain. The corrections are editorial: an omitted discrepancy with the published KLO estimate ε₀ ≈ 0.2385, two arguments stated too briefly, an outdated Du–Lou source statement, and status/attribution. Status: **solved**.

## 1. Statement and sources

- **Verbatim quotation.** The quotation in `README.md`/`paper.tex` matches arXiv:1612.07722v1, printed p. 3, word for word, including the bracketed reference numbers [16] and [12].
- **Published version.** The author-hosted manuscript of the published version (JDE 263 (2017), dated 27 June 2017) states the same conjecture with the same content, but its next sentence differs: it adds the Huang–Wang one-dimensional proof and the numerical estimate "(ε₀ ≈ 0.2385)". The solver noted this only in `references.md`. It must be in the paper, because the certified value contradicts it (see §3.4); added (correction C1).
- **Monotone.** In KLO, "monotone" means the curve has no turn, i.e. λ is a monotone function of α. At ε₀ the curve is strictly increasing with one degenerate point, so "for ε ≥ ε₀ the curve is monotone" holds.
- **Exactly S-shaped.** Du–Lou, pp. 213–214 (definition in the introduction, Remarks (a)–(b)): exactly 1, 2, 3 solutions according as λ ∉ [λ₁, λ₂], λ ∈ {λ₁, λ₂}, λ ∈ (λ₁, λ₂), and exactly two turning points. The paper's definition (two simple zeros of λ_α with signs +,−,+), together with λ(0⁺) = 0 and λ(∞) = ∞ (both re-derived), gives exactly this.
- **Conclusion.** Theorem 2 (`thm:global`) answers the conjecture completely, for all ε > 0 and all amplitudes.

## 2. Mathematical argument: what I re-derived myself

### 2.1 Reduction

- **Change of variables.** With q = ε + ε²u, u/(1+εu) = 1/ε − 1/q. Hence q'' + q'/r = −μ e^{−1/q} with μ = λε²e^{1/ε}, q(0) = a = ε + ε²α, q(1) = ε.
- **Rescaling.** With s = √μ·r the equation becomes parameter-free, and μ = s_ε(a)², where s_ε(a) is the unique radius with q = ε (q is strictly decreasing).
- **Derivative.** λ_α = e^{−1/ε}·2 s_ε(a) s_ε'(a), with s_ε' = −q_a/q_s and q_s < 0, so sign λ_α = sign v(s_ε(a)) with v = q_a. This gives (eq. foldlevel) exactly: λ_α = 0 ⇔ ε = 𝓔(a).
- **Sign of v at the first zero of q.** If S₀(a) is the first zero of q, η(a) = S₀(a)², so η' = −2S₀ v(S₀)/q_s(S₀): sign v(S₀) = sign η'(a), as claimed.
- **Second derivative.** At a fold, d/da[v(s_ε(a), a)] = z, because s_ε' = 0 there. Hence λ_αα has the sign of z, and the fold signs +,−,+ follow.
- **Consistency with the original variables.** q_a = w_α, q_aa = ε⁻² w_αα, and x_red = x_orig + 2 log ε + 1/ε.

### 2.2 Qualitative inputs: hypotheses checked against KLO and Du–Lou

**KLO Theorem 3.1.** Hypotheses: n = 2; f ∈ C²(R̄₊); f > 0 and f' > 0 for u > 0; f''f − f'² < 0. Conclusion: a nontrivial Dirichlet kernel satisfies w > 0 on [0,1). For f(u) = e^{−1/(u+β)}, β ≥ 0, every hypothesis holds; KLO themselves give this example.
- I checked the proof for β = 0, where h'(u) → ∞ at the boundary. Then α₂ = 0 and Case 2(b) of the proof applies; it goes through.
- **Original source.** Du–Lou **Theorem 3** (p. 216) states exactly the needed fact: for n = 2 and every ε ≥ 0, not necessarily small, the Dirichlet kernel at a degenerate positive solution of Δw + μe^{−1/(ε+w)} = 0 does not change sign. No discrepancy.

**Application: "v has at most one zero while q ≥ 0".** Suppose v had zeros S₁ < S₂ with c = q(S₂) ≥ 0. Then w(r) = q(S₂r) − c > 0 solves Δw + S₂² e^{−1/(c+w)} = 0, w = 0 on the boundary. Moreover v(S₂r) is a radial Dirichlet kernel with an interior zero at S₁/S₂, which contradicts Theorem 3. Zeros are simple by ODE uniqueness. The argument is correct, including c = 0.

**KLO Theorem 4.1.** Hypotheses for g(q) = e^{−1/q}, extended by 0:
- g ∈ C²(R̄₊), g(0) = g'(0) = 0, g > 0;
- g strictly log-concave; g'' = g(1−2q)q⁻⁴, so g is convex-concave with β = 1/2;
- g → f₀ = 1 at infinity;
- qg'/g = 1/q > 1 on (0, 1/2).

All verified, and KLO state that the theorem applies to (4.3).
- **Original source.** The theorem is Du–Lou **Theorem 2** plus **Lemma 2** (τ''(0) > 0). On p. 220 Du–Lou state explicitly that "(η₀, v_{η₀}) is the only degenerate solution of (3)". Again no discrepancy for this g.
- **Gap in the paper (C3).** The paper did not derive the signs of η'(a) in (eq. limitcurve) from "only one critical point". Argument: η(a) = S₀(a)² is smooth; η'(a) = 0 iff the solution is degenerate, i.e. only at a_m; η attains its minimum η₀ at a_m. Hence η' < 0 on (0, a_m) and η' > 0 on (a_m, ∞). Added.

**Remark on the published hypothesis list.** KLO Theorem 4.1 requires (4.5) only on (0, β), while its proof invokes Theorem 3.2, which assumes uf' > f for all u > 0. The proof applies the Theorem 3.2 argument only to solutions with v(0) ≤ β, where (4.5) suffices. For the present g this is moot, since Du–Lou prove the result directly.

### 2.3 Delicate steps

- **Tail exclusion, derivative version.** For a > a_m: η' > 0 ⇒ v(S₀) > 0. Since v has at most one simple zero on [0, S₀], it has none. Hence v > 0 at every positive stopping level, and λ_α > 0 whenever ε + ε²α ≥ a_m.
  - For a = a_m: v(S₀) = 0 at q = 0, which is not physical.
  - This is the derivative statement, stronger than the monotonicity of KLO Lemma 4.1. Correct.
- **𝓔(a) → 0 as a ↑ a_m.** The paper only said "by continuity at the limiting nondegenerate fold" (C4). Precise argument: v(·, a_m) > 0 on [0, S₀(a_m)); by continuous dependence, v(·, a) has no zero on [0, S₀(a_m) − δ] for a close to a_m. Hence T(a) → S₀(a_m) and 𝓔(a) → q(S₀(a_m), a_m) = 0. Nondegeneracy is not even needed. Added.
- **Wronskian on (0, 1/2].**
  - Identity: [s(v z_s − z v_s)]' = −s g''(q) v³, recomputed.
  - Integrating to the first zero S of v gives S z(S) v_s(S) = ∫₀^S s g'' v³ ds > 0, because g'' > 0 for q < 1/2 and v > 0 on [0, S). With v_s(S) < 0, this gives z(S) < 0.
  - Combined with 𝓔' = −pz/j, where p, j < 0 (also recomputed), this yields 𝓔' > 0 on (0, 1/2]. Correct.
- **Extension g = 0 for q ≤ 0.**
  - g is C^∞ (all derivatives tend to 0 as q → 0⁺), so the extended problem and its variational equations are classical.
  - After q = 0, v_xx = 0, so v is affine in x.
  - "At most one downward crossing" in the extended problem: a crossing before S₀ forces v(S₀) < 0, and an affine function cannot then cross downward again. Correct.
  - Consequently, a certified downward crossing at a = 3/2 with q > 0 gives η'(3/2) < 0. A certified crossing at 8/5 with q < 0 gives v(S₀) > 0, i.e. η'(8/5) > 0. Hence 3/2 < a_m < 8/5.
- **Assembly.** 𝓔' > 0 on (0, a_c), 𝓔' < 0 on (a_c, a_m), 𝓔 → 0 at both ends, and 𝓔''(a_c) = −pD/j < 0. So every level 0 < ε < ε₀ is crossed exactly twice and transversally, ε₀ exactly once at the cusp, and ε > ε₀ never.
  - For ε > ε₀, positivity of λ_α follows because the unique zero of v occurs below the stopping level, or not at all if a ≥ a_m.
  - The identification ε₀ = ε_c uses the uniqueness of the zero of z on the C-strip.
- **Multiplicity, ordering, stability, analytic exclusions.**
  - Re-derived λ(0⁺) = 0 (α ≥ λ/4 since f ≥ 1) and λ(∞) = ∞ (α ≤ λe^{1/ε}/4).
  - Ordering by continuation: openness via the positive kernel and Hopf, closedness via ODE uniqueness. Correct.
  - Theorem 5 (`thm:tails`): the Wronskian [r(uv' − vu')]' = −rvλ(uf' − f) and the identity uf'/f = 1/(4ε) − (1−εu)²/(4ε(1+εu)²). Correct.
  - Explicit Gelfand solution at ε = 0: re-verified by hand; the maximum is 2 at α = 2 log 2.
- **Original-variable cusp script (`validated_cusp.py`)**, checked by reading:
  - the complex Volterra bounds (|1 + εw| > 1, |f|, |f_w| < e⁸, contraction 10⁴e⁻⁶⁰);
  - the Cauchy radii 10⁴e⁻⁶⁰ k! j! 100^j;
  - the Jacobian layout of F = (w, v, z) in (α, ε, x);
  - the Krawczyk operator and its strict inclusion;
  - the implicit-differentiation formulas for B_α, B_αα, B_ε;
  - the rectangle argument: mean-value bracket of the first zero, B convex in α, positive at both amplitude edges.

  All correct.

## 3. Computations

### 3.1 Re-running the solver's certificates (WSL Fedora, python-flint 0.9.0, outputs in a scratch copy)

| Script | Result | Time |
|---|---|---|
| `check_analytic.py` | 20 × PASS | 1 s |
| `validated_cusp.py` | 23 × PASS, "PASS local ordinary cusp" | 37 s |
| `validated_reduced.py --full --workers 6` | 110 leaves, no bisection, all PASS (E(1/2) > 0.238, E(3/2) > 0.0078, E(8/5) < 0, overlap) | 143 s |

The regenerated `cusp-certificate.json` and `reduced-cover.json` are identical to the delivered ones, apart from runtimes.

**Code review of `validated_reduced.py`.**
- Taylor coefficients: the order-14 polynomial plus the order-15 coefficient evaluated on a Picard tube; Lagrange remainder, valid componentwise.
- Tube: strict inclusion check, with q > 0 required in the tube.
- Origin: the complex Volterra bound with Re q ≥ 0.24 and Lipschitz constant < 18, then Cauchy radii e⁻⁵⁰ m! 8^m.
- Mean-value intersections (m = 2, 1, 0 for q and p): each uses only full-family enclosures of the next derivative, and the window version uses the same time on both sides. Valid.
- Flat step: the bounds b₀–b₃ follow from the monotonicity of e^{−1/q}q^{−m} on (0, 1/m). The second-order step formula U + hV + h²F/2 is valid across q = 0.
- Event windows: v > 0 at the left edge, v < 0 at the right edge, j < 0 on the whole window.
- Cover: contiguous and exact, using decimal balls containing the rational endpoints.
- SciPy only proposes window centres.

I found no rigor defect.

### 3.2 Independent re-certification with a different implementation (`verification/`)

Written from scratch in Julia 1.13 with IntervalArithmetic.jl 1.0.12; no code is shared with the solver. Differences from the solver's method:
- **Derivatives in a** are carried automatically as truncated 3-jets q(a+δ) = Q₀ + Q₁δ + Q₂δ² + Q₃δ³. The Taylor recursion is built from jet reciprocal and jet exponential; no variational equation is written by hand.
- **Origin.** A real-variable bootstrap in t = s²/4: p = −∫₀ᵗ G, Q − A = ∫₀ᵗ p/τ, with a jet ball of radius 10⁻⁸ and t₀ = e⁻³⁰ (e⁻⁶⁰ for the 160-bit runs). There are no complex estimates.
- **Integration.** Order-16 Taylor with the Lagrange remainder on a Picard tube; a separate second-order step for q ≤ 1/40, with bounds valid across q = 0.
- **Event windows** are found automatically, with no floating-point proposals: the last step end where v > 0 for every member, and the first where v < 0.
- **Partition and time bookkeeping.** Dyadic boxes of width 1/512 (564 boxes covering [1/2, 1.6015625]). Every advance of x is checked to be exact.
- **Wrapping control.** Plain box integration wraps catastrophically near q = 0 (v reached ±186), so I added first-order mean-value recentring. It uses the jets themselves (∂ₐQ_m = (m+1)Q_{m+1}) with the midpoint solution integrated in lockstep. This is the same principle as the solver's, implemented independently.
- **Classification.** Each box is labelled N (q > 0, z < 0), C (q > 0, D = c − kz/j > 0) or P (z > 0). The checks are v(x_L) > 0 > v(x_R), j < 0 and p < 0 on the window, and the pattern N…N C…C P…P.

**Results: `cover.jl`, see `cover.log`.**

- 564 leaves of width 1/512 cover [1/2, 1.6015625]; none needed bisection. Pattern: 46 N, 1 C, 517 P, all checks PASS (227 s on 6 threads).
- The C box is [0.58984375, 0.591796875]; it contains the transformed cusp centre a_c ≈ 0.5904717.
- Margins: max z on the N boxes −8.1·10⁻⁴; min z on the P boxes 2.5·10⁻³; min D on the C box 1.72; min q on the N/C boxes 0.2380.
- Point values (refined windows):
  - 𝓔(1/2) ∈ [0.238114, 0.238115] (> 0.238);
  - 𝓔(3/2) ∈ [0.0078977, 0.0078978] (> 0.0078, hence a_m > 3/2);
  - event level at 8/5: [−0.0351864, −0.0351863] < 0, hence a_m < 8/5.
- These agree with the solver's coarse enclosures (0.238109, 0.238124), (0.007889, 0.007907) and (−0.035196, −0.035177).
- Hence z < 0 at the event on [1/2, 0.5898], z is strictly increasing on [0.5898, 0.5918], z > 0 on [0.5918, 1.6], and a_m ∈ (3/2, 8/5). Together with the analytic part on (0, 1/2], this is exactly what the global argument needs. The solver's constants (z < −0.01 on [0.5, 0.58], D > 1, z > 0.01) are not needed in this form; my cover proves the same sign pattern with its own partition.

**Critical value and cusp, entirely in the reduced variables: `cusp_reduced.jl`, see `cusp_reduced.log`.** Steps:
- bisection with rigorous Float64 point evaluations locates the sign change of ζ(a) = z(T(a), a);
- 160-bit interval evaluations at the two bracket ends a⁻ < a⁺ re-prove ζ(a⁻) < 0 < ζ(a⁺);
- a 160-bit box integration of [a⁻, a⁺] encloses ε₀ = 𝓔(a*), T(a*) and D on the bracket;
- one interval Newton step (160-bit ζ(a_mid)) narrows a*; then α_c = (a* − ε₀)/ε₀² and log(λ_c/4) = T − 2 log ε₀ − 1/ε₀.

This uses neither the solver's original-variable shooting nor its local Krawczyk certificate.

All checks PASS (54 s, one thread). Rigorous outward bounds, rounded outward:

| Quantity | Verifier's enclosure | Claim in the paper |
|---|---|---|
| a* = argmax 𝓔 | [0.59047172793910032, 0.59047172793910046] | (0.5904717279294, 0.5904717279490) ✓ |
| ε₀ = 𝓔(a*) | [0.24210616559523735, 0.24210616559523776] | \|ε₀ − 0.24210616559524695\| < 10⁻¹² ✓ (difference ≈ 9.4·10⁻¹⁵) |
| α_c = (a* − ε₀)/ε₀² | [5.94324340648535, 5.94324340648539] | \|α_c − 5.943243406485813\| < 10⁻¹⁰ ✓ (≈ 4.5·10⁻¹³) |
| log(λ_c/4) | [−0.28558378244631059, −0.28558378244629871] | within 10⁻¹² of −0.28558378244629573 ✓ |
| λ_c | [3.00630147886945, 3.00630147886949] | (3.006301478866, 3.006301478873) ✓ |
| D(a*) = ζ'(a*) | [2.4622880, 2.4622882] > 0 | ordinary cusp (B_αα > 0) ✓ |

- **Bracket.** The Float64 bisection gives [0.5904717278899625, 0.5904717280063778], with rigorous signs at both ends, re-proved in 160-bit arithmetic: ζ(a⁻) ≈ −1.21·10⁻¹⁰, ζ(a⁺) ≈ +1.66·10⁻¹⁰.
- **Newton step.** One interval Newton step, with ζ(a_mid) ∈ [2.23322, 2.23326]·10⁻¹¹ and D enclosed on the bracket, narrows a* to width 1.3·10⁻¹⁶.
- **Uniqueness.** By the cover, ζ < 0 left of the C box, D > 0 on it and ζ > 0 right of it. So a* is the unique maximiser of 𝓔 and ε₀ = 𝓔(a*).
- **Consistency checks.** D(a*)·ε₀⁴ ≈ 0.008460 agrees with the solver's B_αα = 0.0084598423 at the cusp (B_αα = ε⁴ q_aaa). In reduced variables, B_ε = j/p > 0 holds automatically at a fold.
- **Origin.** The 160-bit runs start at x₀ = −60. With x₀ = −30, the bootstrap uncertainty (≈ e⁻³⁰M) limits z to about 10⁻¹⁰, enough for the cover but not for the 10⁻¹⁰ claim on α_c.

### 3.3 Direct check on the original problem, without the reduction (`verification/direct_original.py`)

Floating-point shooting in x = log(s²/4) with the variational equation for w_α, DOP853, rtol 10⁻¹². The sign of λ_α is the sign of v at the first zero of w. Evidence only, not proof.

| ε | folds (direct) | fold α (max, min of λ) | a = ε + ε²α |
|---|---|---|---|
| 0.05 | 2 | 1.5608737, 538.56862 | 0.0539022, 1.3964216 |
| 0.15 | 2 | 2.1714079, 42.761507 | 0.1988567, 1.1121339 |
| 0.23 | 2 | 3.8988268, 10.052371 | 0.4362479, 0.7617704 |
| 0.2420 | 2 | 5.6920048, 6.2110940 | 0.5753466, 0.6057465 |
| 0.2422 | 0 | none | none |
| 0.245 | 0 | none | none |

- **Agreement with the paper.** The ε = 0.242 row agrees with the paper's table to 8 digits. The λ_α sign is positive at α = 10⁻³ and at the end of each scan.
- **Second fold for small ε.** a₊ = 1.3964, 1.4216, 1.4709, 1.4950, 1.5069 for ε = 0.05, 0.04, 0.02, 0.01, 0.005. So a₊ ↑ a_m, and a₊ > 3/2 exactly when ε < 𝓔(3/2) ≈ 0.0079. This is consistent with the paper's statement, which is restricted to ε ≤ 0.0078. Its form "α₊ ≈ a/ε² with a ∈ (3/2, 8/5)" is therefore correct only for ε below about 0.0079; for ε = 0.05 one has a₊ ≈ 1.40.
- **Cross-check (`crosscheck.jl`).** Each of the eight folds found by direct shooting satisfies 𝓔(a) = ε to the printed precision in my validated reduced integrator: for example 𝓔(1.39642155) ∈ [0.0499999, 0.0500001] and 𝓔(0.57534657) ∈ [0.241999, 0.242001]. The sign of z is negative at every first fold and positive at every second fold. This checks the fold equivalence λ_α = 0 ⇔ ε = 𝓔(a) and the fold signs on the original equation.

### 3.4 The published estimate ε₀ ≈ 0.2385 and the three-dimensional value

- **The estimate is wrong for the disk.** The JDE version of KLO prints "ε₀ ≈ 0.2385". The certified value, my independent enclosure and the direct shooting all give ε₀ ≈ 0.242106; at ε = 0.242 the curve still has two folds.
- **Dimension three.** The same reduction in dimension three, computed in floating point (`verification/dimension_check.py`), gives max 𝓔 ≈ 0.2387971 at a ≈ 0.6485. This agrees with the Minamoto–Nakao value 0.238797186 (error bound 1.2·10⁻⁴), and supports the solver's reading that their result is three-dimensional.
- **Where the discrepancy belongs.** It concerns a numerical remark of KLO, not a theorem, but a referee would ask about it. It is now discussed in the paper (C1).

## 4. Novelty

**Searches re-run.**
- OpenAlex: KLO 2017 has 7 citers; Du–Lou 2001 has 42 citers; all titles screened. The 2026 citers concern the Kernevez–Thomas model on the unit ball and a regularized MEMS model; the 2025 citer is Taira's monograph.
- Semantic Scholar: the same 7 citers of KLO; the Du–Lou query was rate-limited.
- Korman's publication list for 2017–2026, including two Korman–Schmidt continuation papers, read in their author PDFs.
- arXiv API `all:"perturbed Gelfand"`: two hits, both already known.
- Web searches: Huang, Hung, Wang (all one-dimensional or other operators); Lessard–Pugliese 2024/2026 and Blanco–Lessard 2026 (rigorous cusp validation in other models).
- GitHub: `ipitchford` (169 repositories) and `mmislan` (24 repositories), names and descriptions screened; `gh search repos` returned nothing; Zenodo returned nothing.

**Not read in full.** Minamoto–Nakao (2007): I did not read the full text, but the dimension-three computation in §3.4 corroborates the solver's reading. Taira's 2025 monograph: description only, as the solver noted.

**Conclusion.** No prior resolution of the two-dimensional global conjecture was found. This is a bounded search, not a proof of novelty.

## 5. What I only read (not recomputed)

- **The published theorems of KLO and Du–Lou.** I checked hypotheses and read the proofs of KLO Theorem 3.1 and Du–Lou Theorem 3. For Theorem 2 / Lemma 2 / KLO Theorem 4.1 I read the proofs, but I did not re-prove the bifurcation-direction computation τ''(0) > 0 or the Korman–Li–Ouyang 1996 results invoked by KLO.
- **The solver's original-variable Krawczyk certificate.** I reviewed it as code and reran it. Its numerical conclusions (ε_c, α_c, λ_c, and the ordinary-cusp signs) are reproduced in the reduced variables, but I did not reimplement the 14-dimensional original-variable system.
- **Gidas–Ni–Nirenberg radial symmetry.** Classical; not re-checked.

## 6. Corrections made

| # | File | Correction |
|---|---|---|
| C1 | `paper.tex`, `README.md`, `references.md` | Added the published KLO wording, its estimate "ε₀ ≈ 0.2385", why that estimate is incompatible with the result, and the three-dimensional comparison. |
| C2 | `paper.tex`, `references.md`, `README.md` | Du–Lou is now read in full. Added the precise original sources: Theorem 3, n = 2 and any ε ≥ 0, for exactly e^{−1/(ε+w)}; Theorem 2, Lemma 2 and the p. 220 uniqueness of the degenerate solution. Removed "full text not obtained". |
| C3 | `paper.tex` §"Two established qualitative inputs" | Added the derivation of the sign pattern of η'(a) from "only one critical point". |
| C4 | `paper.tex` §"Completion" | Replaced "continuity at the limiting nondegenerate fold gives 𝓔(a_m) = 0" with an explicit argument. |
| C5 | `paper.tex` abstract | "Convexity handles arbitrarily small transformed amplitudes" → "A Wronskian identity …" (the argument is a Wronskian identity that uses the convexity of g on (0, 1/2]). |
| C6 | `paper.tex`, `README.md` | Status and attribution: "solved", verifier credit in `\PaperCredits` and in the README header; the "not independently verified" sentence of the local Theorem 3 (`thm:local`) and the "Verification status" section rewritten; numerical section cites the verifier's direct check. |
| C7 | `README.md` | Removed the link to `PROGRESS.md`, which is git-ignored and therefore dead in the public repository; added the verification files and commands. |
| C8 | `paper.tex` §"Completion" | Removed leftovers of earlier drafts: "replaces the earlier coarse amplitude bound 10¹⁰" (no such bound appears in the paper) and "phase-2 decomposition". |
| C9 | `references.md` | Status line updated; Du–Lou entry updated; the verifier's additional searches appended. |

None of these changes alters a theorem statement, a constant or a certificate.

## 7. Reservations

1. **Trusted computing base.** Both computer parts depend on interval libraries (Arb/FLINT via python-flint; IntervalArithmetic.jl) and on the correctness of two independent implementations. The agreement of two independent implementations, on different libraries and with different origin, remainder and event-window designs, makes a shared implementation error unlikely.
2. **Published qualitative theorems.** The proof relies on Du–Lou Theorems 2 and 3 and Lemma 2 (equivalently KLO Theorems 3.1 and 4.1). These are refereed results; I checked their hypotheses and read their proofs, but did not re-prove the Crandall–Rabinowitz bifurcation-direction computation.
3. **Tolerances.** The local precision claims (|ε₀ − e₀| < 10⁻¹², |α_c − a₀| < 10⁻¹⁰) rest on the solver's 160-bit certificate and on my 160-bit reduced-variable run. The global theorem itself needs only the qualitative sign pattern, which both covers establish with wide margins (min |z| ≥ 8·10⁻⁴ on the N/P boxes of my cover, min D ≥ 1.7 on the C box).
4. **Novelty** rests on a bounded literature search.
