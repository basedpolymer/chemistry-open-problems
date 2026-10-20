# Independent verification: Jang–Tran Conjecture 2 (discrete critical coagulation–fragmentation)

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, 20 October 2026. Inputs read: `PROTOCOL.md`, the
verifier brief, `README.md`, `paper.tex`, `references.md`, `references.bib`, `scripts/`, `certificates/`
(including the source texts of Jang–Tran arXiv:2409.17974v1 and Van arXiv:2607.02293v2, and the page image of
Conjecture 2), `part-ii-opus/part-ii.md` and `part-ii-opus/scripts/`. Not opened: the `PROGRESS.md` files,
`FINAL-REPORT.md`, the other brief files, `tmp/`.

**PDF requests.** None needed. (Tran–Van, CPAM 2022, doi:10.1002/cpa.21979, and the two Escobedo et al. papers
are cited only for attribution and were not needed for the verification.)

## Verdict

**CONFIRMED WITH CORRECTIONS.**

The three theorems of the manuscript are correct as stated: for the kernels a(i,j) = ij, b(i,j) = 1,

1. for every nonnegative datum with finite mass m ∈ (0,1] there is exactly one mass-conserving solution (uniqueness
   in the weakest reading of Jang–Tran's Definition 1: absolutely continuous in time, weak form a.e.), and it is
   classical; no second moment is needed;
2. for every m ∈ (0,1) the recursion (1.9) gives strictly positive numbers with Σρ̃ₙ = m(1−m), Σnρ̃ₙ = m; this is
   the unique stationary solution of mass m and attracts every mass-conserving solution in ℓ¹₁; for m = 1 every
   component tends to 0;
3. ρ̃ₙ = b/(2√π) R^{−n} n^{−5/2}(1 + O(1/n)) with R(m) > 1 and the stated b.

This proves Conjecture 2 of Jang–Tran (both parts, the whole window, without their assumption (1.3)). I found no
mathematical gap. The corrections concern the **statement of the consequence for Conjecture 1** (the only claim
that could be read as an overclaim), one citation, and documentation; they are listed below.

## 1. Quotations and notions (checked against the sources)

* Conjecture 2: compared word by word with the extracted arXiv v1 text and with the rendered page image (printed
  p. 7, end of Section 1.4). The README/paper quotation is exact up to typesetting. The journal version
  (open-access HTML of doi:10.1186/s13662-025-03946-4, retrieved by the verifier) states the same conjecture,
  slightly reworded, at the end of its Section 1.4.
* Conjecture 1 (arXiv v1, Section 1.2): "The critical mass of the initial data should be m = m₁(ρ₀) = 1. For m > 1,
  there is gelation. For 0 < m ≤ 1, (1.1) has a unique mass-conserving solution for all time", for data satisfying
  (1.3). Journal version: "For m > 1, gelation occurs."
* Van, Section 1.2, printed p. 4: "so the discrete analogue of the window 1/2 ≤ m ≤ 1 remains open as of this
  writing." Exact.
* Notions: system (1.1)–(1.2), assumption (1.3) m₂(ρ₀) < ∞, weak form (1.5) and Definition 1, (1.6), (1.7), (1.9)
  are reproduced correctly. Jang–Tran do not define "mass-conserving solution" or the time regularity in
  Definition 1; the manuscript's reading (absolutely continuous in time, (1.5) a.e., m₁ constant) is the weakest
  reasonable one, uniqueness is proved in that class and existence in the classical class, which covers every
  stronger reading. A stationary solution is a time-independent solution with m₁ = m, as in Conjecture 2(i).
* Claims beyond the conjecture:
  - *m = 1*: correctly stated (global existence, uniqueness, mass 1 for all t, components → 0; consistent with
    Jang–Tran's Corollary 5.2 / Remark 1 and Lemma 5.6). The only place where m < 1 is needed is the uniform-in-time
    slope bound (Lemma 20(b)), which is used only for the stationary state. Checked.
  - *No second moment*: correct. Nothing in Sections 1–3 uses m₂(ρ⁰): F₀ is C^∞ on (0,∞) for any finite-mass datum,
    the initial bounds (W₀ ≥ 0, −m/e ≤ τF₀'' ≤ 0) use only Σjρⱼ = m, the truncation puts the tail mass at size L.
    The manuscript correctly does not claim m₂(ρ(t)) < ∞ for t > 0.
  - *Conjecture 1*: the 0 < m ≤ 1 half is Theorem 4 (for all finite-mass data, a fortiori under (1.3)). The m > 1
    half ("there is gelation") is supplied only in the sense of Jang–Tran's Theorem 1.1: no solution conserves mass
    for all time. If "gelation" is read as "there is a global weak solution whose mass drops in finite time", the
    existence part of that statement is proved neither by Jang–Tran nor here. The original wording "this settles
    their Conjecture 1" was therefore **qualified** (correction C1). I also checked that the m > 1 non-existence
    follows independently from the manuscript's Proposition 7: its derivation of n₀(t) = m(1−m) +
    e^{−t/2}(n₀(0) − m(1−m)) holds on any interval [0,T], and for m > 1 it forces n₀ = 0 (hence m₁ = 0) at
    T* = 2 log(1 + n₀(0)/(m(m−1))).

## 2. What I re-derived by hand

Numbering as in the compiled `paper.pdf` (one shared counter): Theorems 4 (well-posedness), 5 (stationary
states), 6 (decay); Proposition 7 (closed system); Lemma 10 (initial data), Lemma 15 (invariant region), Lemma 16
(barrier), Proposition 18 (global solution), Lemma 19 (comparison), Lemma 20 (boundary slope), Lemma 22
(truncations), Remark 24 (F^ρ = F), Lemmas 27–28 (singularity analysis). The README calls Theorems 4–6
"Theorems 1–3".

* **Transform equation.** With g(n) = 1 − zⁿ, z = e^{−x}: coagulation gives −½(m − F_x)² (since
  g(j+k) − g(j) − g(k) = −(1−z^j)(1−z^k)); fragmentation gives ½(m + F_x) − F((1−z)^{−1} − ½), and
  (1−z)^{−1} − ½ = ½ + 1/(eˣ−1) = ½coth(x/2) = 1/τ, τ = 2tanh(x/2). Hence F_t + H_m(x,F,F_x) = 0 with
  H_m = ½(p−m)(p−m−1) + u/τ − m, which is Jang–Tran's (1.7). For the truncation, the loss
  Σₙ g(n) nρₙ Σ_{j≤L−n} jρⱼ symmetrises to ½Σ_{i+j≤L}(g(i)+g(j))ijρᵢρⱼ, giving exactly
  −½(m − F^L_x)² + E_L with E_L = ½Σ_{i,j≤L, i+j>L} ij ρᵢρⱼ(1−zⁱ)(1−z^j) ≥ 0. Mass conservation and
  quasi-positivity of the truncation, and F^L(·,0) ≥ F₀ (monotonicity of s ↦ (1 − e^{−sx})/s), checked.
* **Characteristics and W.** Ẋ = H_p = P − m − ½, Ṗ = −H_x − pH_u = kU/τ² − P/τ, U̇ = pH_p − H =
  ½(P² + m − m²) − U/τ, with k = τ' = 1 − τ²/4. Substituting U = (τ(m+P) − W)/2 into Ẇ gives
  Ẇ + (2+k)W/(2τ) = (1−k)(m² − P²) = q(2m−q)τ²/4. Initial identity
  τj(1+z^j) − 2(1−z^j) = 2(1−z)/(1+z) Σ_{i<j}(1−zⁱ)(1−z^{j−i}) re-derived for general j (telescoping
  (1+z)(1+…+z^{j−1}) = 1 + 2(z+…+z^{j−1}) + z^j).
* **Invariant region (Lemma 15).** Each face of V = P⁻ + (P−m)⁺ + U⁻ + W⁻ checked, including the use of k ≤ 1
  and m + P > 0 on {P > m}, the strict sign on {U < 0} (needs m ≤ 1), and the source term on {W < 0}. Grönwall
  from V(0) = 0. Consequences (U ≤ mτ, λ ≤ (m+P)/2, speed bounds, t_e ∈ [y/(m+½), 2y], X → 0) checked.
* **Barrier and J > 0 (Lemma 16).** Strip condition Ȧ = −A/τ for A = Υ − PJ; Π̇ = −J(κ''U + 2κ'P) − κΠ with
  κ' = −k/τ², κ'' = k/(2τ) + 2k²/τ³ (using k' = −τk/2); Riccati for r = Π/J; for Γ = τr the Γ-terms are
  −Γ² − Γ(1 + k(q+½)) and the remaining terms reduce to −2kD exactly because τ²/2 + 2k = 2. At Γ = −1 the right
  side is k(q + ½ − 2D) ≥ k/2 since W = τ(q − 2D) ≥ 0; at Γ = √(2m) it is ≤ −2m − √(2m) + 2km < 0. The
  continuation argument for J (J = exp∫r with r ≥ −1/τ(X(t₁))) is correct. This is the step where the
  half-slope invariant replaces Jang–Tran's restriction m < 1/2 (their Lemma 3.6 needs 1 − 2δ > 0 with δ ≤ m;
  checked in their text).
* **Coverage (Proposition 18).** Interval property: for s < s*, X(s;η) ≥ X(s;y) ≥ X(t;y) > 0 (J > 0 and X
  decreasing), so the characteristics stay in a compact subset of the domain and a uniform extension time gives
  s* = t; s* > 0 because t_e(η) ≥ η/(m+½). Ω_t open by lower semicontinuity of t_e. Surjectivity: upper end by
  X ≥ y − (m+½)t; lower end by the limit characteristic of label a(t), which would stay above ℓ > 0, contradicting
  X(s;a) → 0. Hence Φ is a bijective local diffeomorphism; F = U∘Φ⁻¹, F_x = P, F_t = −H_m. Behaviour as x → ∞
  needs no condition (labels y → ∞ cover large x).
* **Comparison (Lemma 19).** H_m(p,u) − H_m(p̃,ũ) = (p−p̃)((p+p̃)/2 − m − ½) + (u−ũ)/τ; the backward
  transport curve moves away from x = 0 with speed in [½, m+½], so no boundary condition is needed. Applies to
  G = F^L (C¹, 0 ≤ F^L_x ≤ m).
* **Boundary slope (Lemma 20).** λ̇ formula re-derived (the ε-term is exactly −λτ²v/(4τ)); with m ≤ 1,
  1 − m + λ ≥ λ gives λ̇ ≥ (λ/τ)((m−λ)/2 − ε). Part (a): the characteristic through (x,t), x < δ < t/2, crosses
  X = δ at s₁ ∈ (t − 2δ, t), so its entry value is bounded below by λ_* = min_{[0,t]}F(δ,·)/τ(δ) > 0, a quantity
  that depends only on (δ,t), not on the characteristic; the trapping of {λ ≥ m − θ} and the power-law growth
  λ(X) ≥ λ(δ)(δ/X)^{θ/6} then give F(x,t)/τ(x) ≥ m − θ for x small. This is uniform over the characteristics
  arriving near x = 0 at the fixed time t, which is all that the no-mass-loss step needs. Part (b) (m < 1):
  dη/dX ≥ αη/X − 2CX with α = (1−m)/(2m+1), C = m(m+½)/4; integrating (X^{−α}η)' over [x,δ] gives the stated K,
  independent of t ≥ t₀ because the reference level X = δ < t₀/2 is reached at a positive time.
* **Existence (Section 3).** Bounds 0 ≤ ρ^L_n ≤ m/n, |ρ̇^L_n| ≤ 3m²/2 + 3m/2; Arzelà–Ascoli; number tail
  Σ_{n>N}ρ^L_n ≤ m/(N+1) gives F^{L_k} → F^ρ. No mass loss: M(t) ≥ F^ρ(x,t)/x ≥ F(x,t)/x → m (using
  1 − e^{−nx} ≤ nx and Lemma 20(a)); Fatou gives M ≤ m. Tightness at each fixed t follows from coordinatewise
  convergence and equal masses, so μ^{L}_n(s) → m for every s, and the loss term passes to the limit by dominated
  convergence (bound m²). Discrete Scheffé, the ℓ¹₁ → ℓ¹ Lipschitz bounds (3/2 for both parts) and the C¹(ℓ¹)
  regularity checked.
* **Uniqueness (Proposition 7).** g = 1_{n} and g ≡ 1 in the weak form; ṅ₀ = ½(m − m² − n₀); the closed system
  ρ̇ₙ + cₙρₙ = n₀ − Σ_{j<n}ρⱼ + ½Σ i(n−i)ρᵢρ_{n−i}, cₙ = mn + (n+1)/2, is lower triangular; induction. Correct,
  and simpler than Jang–Tran's viscosity-solution proof.
* **Stationary state (Section 4, proof of Theorem 5).** ρₙ(t) → ρ̃ₙ by induction (u̇ + cu = f, f → f_∞);
  Σρ̃ₙ = m(1−m) by the uniform number tail; Σnρ̃ₙ ≥ m from F^ρ ≥ F ≥ τ(m − Kx^α) uniformly for t ≥ 1 (this is
  exactly where the uniform-in-time bound is needed, and it is available because δ and K do not depend on t);
  Fatou for ≤. Stationarity, strict positivity (numerator ≥ (n−1)ρ̃₁ρ̃_{n−1}), uniqueness (a stationary σ of mass
  m has Σσₙ = m(1−m) and satisfies (1.9)) and relaxation in ℓ¹₁ checked. For m = 1 the recursion gives 0.
* **Singularity analysis (Section 5).** Generating-function identity (1−z)(aP − P²) − 2zN + (1+z)H = 0 re-derived
  from K(ρ̃) = 0; the r-equation; the blow-up exclusion near z = 1 ((1−z)r → 1/3 versus → 0); analytic
  continuation across z = 1 via the regular ODE and the identity theorem; R < ∞ via P' ≥ P/(3z); P(R⁻) = m + ½
  (otherwise regular continuation contradicts Pringsheim); the quadratic turning point of z(p) with
  b² = RD(R,p_*)/B(R) = ((R−1)²(m+½) + 3R/2)/(R²−1); uniqueness of the dominant singularity (|P(w)| < P(R) for
  |w| = R, w ≠ R, by ρ̃₁, ρ̃₂ > 0, then the implicit function theorem and H' = P/z); Δ-analyticity; transfer with
  Γ(−3/2) = 4√π/3. Correct. Note that once R > 1 is known, P(1) = m also follows from (P) at z = 1, consistent with
  Theorem 5.

## 3. What I re-ran

| Script | Result |
|---|---|
| `scripts/check_evolution_identities.py` | 24 PASS, ALL PASS (~8 s). I read the code: the checks are genuine (finitely supported generic sequences for the transform equations, symbolic X,P,U for the characteristic identities). |
| `scripts/check_identities.py` | all PASS; `certificates/algebra.json` regenerated byte-identical. |
| `scripts/check_coefficients_arb.py` | 18 PASS; `certificates/coefficients-arb.json` regenerated byte-identical. |
| `scripts/explore_recursion.py` | PASS for the six exact rational coefficients. |
| `scripts/novelty_audit.py` | re-run in a scratch copy (so as not to overwrite the recorded raw responses); same retrievals and same known failures as recorded in `references.md`. |
| `scripts/numerics_evolution.py` | As delivered: 52 PASS, then stalled (> 38 min, stopped) in the long-time check m = 0.6, N = 1600, T = 30, where LSODA switches to its stiff mode with dense 1600 × 1600 Jacobians — the same place as in the solvers' session (`verification/numerics_evolution_lsoda_stopped.log`). After switching that one integration to DOP853 (correction C7): **59 PASS, 0 FAIL in 4 min 5 s**, including the checks never completed before (long-time truncation = closed system, difference 4.7·10⁻¹³; N6 at m = 1; heavy-tail comparison, min(F^N − F) = 0) (`certificates/numerics_evolution.log`). |

## 4. Independent checks (`verification/`, written from scratch, different methods)

* `verification/verifier_stationary.py` (log `verifier_stationary.log`, ALL PASS, ~20 s):
  - exact rationals, n ≤ 40, m ∈ {1/2, 51/100, 3/4, 9/10, 99/100, 999/1000}: ρ̃ₙ > 0, number remainders > 0,
    partial mass < m; m = 1 gives ρ̃ ≡ 0;
  - Arb balls (300 bits), n ≤ 1500, m ∈ {513/1024, 3/4, 19/20}: same inequalities, rigorous but finite;
  - the radius R(m) computed independently by integrating the regularised r-equation across z = 1 and then the
    inverse function z(P) up to P = m + ½: R(0.6) = 1.0237248…, R(0.75) = 1.0048049…, R(0.9) = 1.0002643…;
    with Arb coefficients, ρ̃ₙR^n n^{5/2}/(b/(2√π)) − 1 = −2.6·10⁻², −1.3·10⁻² at n = 1500, 3000 (m = 0.6) and
    −6.0·10⁻², −3.2·10⁻² at n = 3000, 6000 (m = 0.75): an O(1/n) error, Richardson-extrapolated limit within
    7·10⁻⁴ and 4·10⁻³ of 1 (Theorem 6, constant included);
  - Σnρ̃ₙ (partial sum + asymptotic tail) = m to 10⁻¹⁶ at m = 0.6 and 0.75 and to 3·10⁻⁷ at m = 0.9 (floating
    point, n ≤ 30000); the paper's figure "mass above size 3000 ≈ 2·10⁻⁹ at m = 0.75" is confirmed (2.18·10⁻⁹).
* `verification/verifier_dynamics.py` (log `verifier_dynamics.log`): datum ρ⁰ⱼ = c j^{−5/2} (infinite second
  moment) at m = 0.9 and m = 1; fixed-step RK4 for the truncated and closed systems, implicit Radau for the
  characteristics. ALL PASS (12 checks, 7.4 min):
  - truncations L = 250, 500, 1000, 2000 stay nonnegative and conserve mass; their sizes ≤ 20 converge to the
    closed triangular system (which presupposes mass conservation of the limit): max difference at t ∈ {1, 3}
    8.2·10⁻³, 3.3·10⁻³, 9.1·10⁻⁴, 1.3·10⁻⁴ at m = 0.9, and 2.1·10⁻², 1.3·10⁻², 7.8·10⁻³, 4.1·10⁻³ at m = 1
    (slower at m = 1, where 11 % of the mass sits above size 1000 at t = 3 for L = 2000; the mass above L/2
    decreases with L in both cases);
  - along the characteristics used: 0 ≤ P ≤ m, U > 0, W ≥ 0, min Γ ≈ −0.16 > −1, J > 0.5;
  - comparison F^L ≥ F at x ∈ {0.02, 0.1, 0.5, 2}, t ∈ {1, 3}, the gap decreasing with L to ~10⁻⁶–10⁻⁵ at
    L = 2000 (consistent with F^ρ = F);
  - m = 1: F(x,t)/τ(x) = 0.555, 0.843, 0.974, 0.997 at x = 10⁻¹…10⁻⁴ (t = 1) and 0.228, 0.443, 0.724, 0.932
    (t = 3): convergence to 1 at each fixed t, slower for larger t, as Lemma 20(a) (non-uniform at m = 1)
    allows.

## 5. Corrections applied

* **C1 (claim, README and paper).** "This settles their Conjecture 1" now states that gelation for m > 1 is
  understood as in Jang–Tran's Theorem 1.1 (no solution conserves mass for all time), notes that this also follows
  from Proposition 7 (n₀ vanishes at T* = 2 log(1 + n₀(0)/(m(m−1)))), recalls that Conjecture 1 assumes (1.3),
  and says that the existence of global weak solutions losing mass for m > 1 is not addressed. Abstract
  reworded accordingly. Conjecture 1 is now quoted in the README.
* **C2 (citation, paper).** The Δ-domain is Definition VI.1 in Section VI.3 of Flajolet–Sedgewick, not
  "Section VI.1" (checked in the book text).
* **C3 (bibliography).** `references.bib`: the Flajolet–Sedgewick note listed Lemma IX.2 and Theorems VI.1–VI.4,
  which the final manuscript does not use; it now lists Theorem IV.6, Definition VI.1, Theorems VI.1 and VI.3. The
  Jang–Tran note claimed "journal version, p. 7", which was not verifiable; it now says that the journal version
  states the same conjecture at the end of its Section 1.4 (checked on the article HTML).
* **C4 (wording).** Van's continuous result is an unrefereed preprint: "Van settled the window" →
  "Van's preprint extends existence to the window" (paper and README).
* **C5 (status and credits).** README status "solved (preliminary)" → "solved" (with the verdict and "not
  refereed"); attribution and `\PaperCredits` now read "Independent verification: Claude Opus 5.5 (effort xhigh)";
  the paper subtitle and abstract no longer say verification is pending; the "verification priorities" paragraph of
  the paper is replaced by a short account of the verification, and Section 7 is retitled "Computer checks and
  independent verification".
* **C6 (references.md).** Flajolet–Sedgewick entry no longer says it is used "in the stationary continuation proof"
  without qualification; a section records what the verifier read (including the journal HTML) and the repeated
  novelty search.
* **C7 (scripts and README numerics).** `scripts/numerics_evolution.py`: the long-time truncated integration
  (m = 0.6, N = 1600, T = 30) now uses DOP853 instead of LSODA (LSODA never finished it in the solvers' session nor
  in my first rerun); header runtime corrected ("about 6 minutes (2 threads)" was not achievable with LSODA; now
  4–6 min on one thread); stale labels "Lemma 4.6", "Section 6", "Lemma 5.6(a)" from an earlier draft replaced by
  the paper's statements (the N2 check uses a τ-variable variant of the bound of Lemma 20(b), now said so). The
  complete log is `certificates/numerics_evolution.log`; README "Reproduce" and run-status paragraphs and the file
  table updated (verifier scripts, logs, `verification.md`), and the README now gives the correspondence between its
  Theorems 1–3 and Theorems 4–6 of the PDF. The paper's paragraph on the stopped run is replaced by a pointer to the
  README.
* `paper.pdf` recompiled.

No change to any proof was needed.

## 6. Remaining caveats

* This is an AI verification, not peer review. The arguments are elementary but long; a human referee should
  still read Sections 2–3.
* "Gelation for m > 1" in Conjecture 1 is settled only in the non-existence sense (C1).
* The quantitative statements are proved for the auxiliary characteristic solution F; the identification F^ρ = F
  (Remark 24) is not needed and is supported numerically (F^L − F decreases to 0 with L).
* Numerical checks (both the solvers' and mine) are evidence only. At m = 1 the convergence of truncations is slow
  (mass moves to large sizes), as expected; nothing in the proof depends on it.
* Novelty: bounded search (databases, arXiv API, Zenodo, GitHub, one web search) up to 20 October 2026; no
  resolution of the discrete window was found. The root `README.md` results table has no entry for this problem
  yet; it should be added by the repository owner/orchestrator (the brief did not list it among the files the
  verifier may edit).
