Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Independent verification: the Banaji–Boros–Hofbauer three-reaction example

Date: 18 October 2026.

Inputs used:
- the claimed result;
- the deliverables in this folder (`paper.tex`/`paper.pdf`, `README.md`, `references.md`, `scripts/`, `certificates/`);
- `PROTOCOL.md`;
- the primary source, which I downloaded myself: arXiv:2304.02303v1 (25 pages) and the published version (Stud. Appl.
  Math. 152 (2024) 249–278, open access as PMC10953466).

I did not read `PROGRESS.md`, any solver log, `travail/` or `rapports/`.

## Verdict: **CONFIRMED WITH CORRECTIONS**

Theorem 1 is correct and answers the question exactly as it was asked. The example
2X → 3X + 2Y, X + Y → 0, Y → X + Y has no periodic orbit in the closed positive quadrant for any κ > 0, so in
particular no limit cycle. This is a stronger statement than the one asked for, not a weaker variant.

Theorem 2 (the subfamily) is also correct. Periodic orbits can only occur for c₁ < d₁ < 2c₁, which is exactly the
Andronov–Hopf region of BBH. The Dulac proof is complete and short, and every step was re-derived by hand and with
independent code.

Corrections were needed, but none affects the mathematics:
- **Wrong third reaction.** The subfamily was written with the third reaction **Y → (1+c₃)X + Y** in the abstract,
  in Theorem 2, in the README and in the script docstring. With BBH's notation (target = source + reaction vector)
  and the ODE that is actually analysed (ẋ = c₁κ₁x² − κ₂xy + c₃κ₃y), the reaction is **Y → c₃X + Y**. With the
  printed version, the case c₃ = 1 would be Y → 2X + Y, which is not the BBH example.
- **Theorem numbering.** "Theorem 9" is the arXiv numbering. In the published version, which the bibliography
  cites, the same classification is **Theorem 3**.

The other corrections are wording or rigour improvements. I also added one remark (sharpness, and the extra values
c₁ ≤ 0), checked symbolically.

## 1. Source and statement

- **Verbatim quote.** arXiv v1, §6 Conclusions, p. 22, checked on a rendered page image: "However, we leave it open
  whether the networks in cases 9 and 10 admitting no Andronov–Hopf bifurcations can have a limit cycle. For
  example, we do not know whether the system [2X →κ₁ 3X + 2Y, X + Y →κ₂ 0, Y →κ₃ X + Y; ẋ = κ₁x² − κ₂xy + κ₃y,
  ẏ = 2κ₁x² − κ₂xy] admits a limit cycle."
  - The quote in `paper.tex`/`README.md` is exact.
  - The published version (§6) has the same text, and its displayed formula image is identical. The only
    difference is that it refers to "Theorem 3" instead of "Theorem 9".
- **ODE ↔ network.** I built the mass-action ODE from the source and target complexes with a generic constructor
  (V1 in `verification/verify_symbolic.py`). It is exactly the displayed ODE, with rate laws κ₁x², κ₂xy, κ₃y.
- **Case 9 and notation.** In BBH (13)–(16), reaction i is aᵢX + bᵢY → (aᵢ+cᵢ)X + (bᵢ+dᵢ)Y, with integer
  coefficients and nonnegative targets.
  - Case 9 has sources 2X, X+Y, Y, in this order, which is positively oriented (determinant 1).
  - The subfamily is c = (c₁, −1, c₃), d = (d₁, −1, 0), so its targets are (2+c₁)X + d₁Y, 0 and **c₃X + Y**.
    The example is c₁ = c₃ = 1, d₁ = 2.
- **BBH Hopf condition.** For case 9 with c₃ > 0, the condition is: c₁ > 0, c₂ = −1, d₁ > 0, d₂ ≥ −1, d₃ ≥ −1 and
  ½(d₃/c₃ + d₁/c₁) < d₂/c₂ < d₁/c₁.
  - I read this from the rendered page 15 of arXiv v1, and it is identical in the journal's Theorem 3.
  - With d₂/c₂ = 1 and d₃ = 0 it becomes c₁ < d₁ < 2c₁ (V12).
  - The example sits on the boundary d₁ = 2c₁, where the condition fails by equality.
- **Integrality.** The coefficients are integers (stoichiometric coefficients), so "positive integers" in Theorem 2
  is the right hypothesis. The proof uses neither integrality nor c₁ ≥ 1; see the Remark added in §4.

The result answers the example question in full. It does not claim the broader question (all non-Hopf networks of
cases 9 and 10), and the manuscript says so.

## 2. What I recomputed myself, and how

All in `verification/verify_symbolic.py` (sympy, exact; log `verification/verify_symbolic.log`, 34 checks, ALL PASS),
plus checks by hand.

1. **Dulac computation.**
   - For B = x^α y^β I computed x^{1−α}y^{1−β} div(BF) in general.
   - With B = y^{−b}: y^{b+1} div(BF) = −b d₁κ₁ x² + (2c₁κ₁ + (b−1)κ₂) xy − κ₂ y².
   - I obtained the same result by a second, direct computation. The paper's coefficients at b = b* agree.
2. **Optimal b.**
   - Δ(b) = (2c₁κ₁ + (b−1)κ₂)² − 4b d₁κ₁κ₂ is a convex quadratic in b, with leading coefficient κ₂².
   - Its unique minimiser is b* = 1 + 2κ₁(d₁−c₁)/κ₂, with Δ(b*) = 4d₁κ₁(κ₁(2c₁−d₁) − κ₂).
   - The paper's statement that b minimises Δ among the functions y^{−b} is therefore correct.
3. **Sign.**
   - With d₁ = 2c₁ + e (e ≥ 0), −Δ/4 = 2c₁eκ₁² + 2c₁κ₁κ₂ + e²κ₁² + eκ₁κ₂, a polynomial with positive coefficients,
     so Δ < 0.
   - Since the y² coefficient −κ₂ is negative and Δ < 0, the form is negative definite. In particular the x²
     coefficient −b*d₁κ₁ is negative (b* > 1).
   - As a spot check, 400 random exact rational samples gave Q < 0.
   - Example: Δ = −8κ₁κ₂.
4. **Domain, Dulac and the boundary** (checked by hand).
   - The open quadrant is convex, hence simply connected, and B = y^{−b} is C^∞ there. Since div(BF) < 0
     everywhere on it, Bendixson–Dulac excludes any periodic orbit contained in it.
   - Boundary: on y = 0, x > 0 we have ẏ = d₁κ₁x² > 0; on x = 0, y > 0 we have ẋ = c₃κ₃y > 0; the origin is an
     equilibrium.
   - So a periodic orbit through a positive-axis point p would be outside the closed quadrant just before reaching
     p. It therefore cannot lie in the closed quadrant, and cannot contain the origin either.
   - The original wording of this step ("it would have to reach that point from the inside") was loose; I
     rewrote it with the one-sided Taylor argument.
   - The statement covers all periodic orbits, not only limit cycles.
5. **Case d₁ ≤ c₁.**
   - The rank of Γ is 2 (since c₃ ≠ 0), and ker Γ = span(u) with u = c × d = (c₃, c₃d₁, d₁ − c₁).
   - u is neither positive nor negative when d₁ ≤ c₁: this includes d₁ = c₁, where u₃ = 0. So there is no positive
     flux vector and no positive equilibrium.
   - A periodic orbit in the open quadrant would enclose an equilibrium in its Jordan interior, which lies in the
     open quadrant. Contradiction. The boundary argument above excludes the axes.
6. **Boundary cases.**
   - d₁ = 2c₁ falls in the Dulac case, with Δ = −8c₁κ₁κ₂ < 0. d₁ = c₁ falls in the "no equilibrium" case. Both
     are handled correctly.
   - When Δ = 0, which can only happen for c₁ < d₁ < 2c₁, the form is Q = −κ₂(y − r x)² with r = d₁/(2c₁ − d₁) > 0
     (V9). It is semidefinite and vanishes only on a ray, a set of measure zero. Green's formula still gives a
     contradiction, so Dulac still applies.
   - When Δ > 0, Dulac with y^{−b} gives nothing. This is expected: see the next item.
7. **New observation: the Dulac threshold coincides with the Hopf curve** (V7, V8).
   - At the positive equilibrium (d₁ > c₁), y*/x* = d₁κ₁/κ₂ and det J > 0.
   - tr J = c₃d₁κ₃(κ₁(2c₁−d₁) − κ₂) / (κ₂(d₁−c₁)), so tr J has the sign of κ₁(2c₁−d₁) − κ₂ = Δ(b*)/(4d₁κ₁).
   - The Dulac function therefore excludes periodic orbits exactly when tr J ≤ 0. On the other side, the
     equilibrium is a repellor and BBH's supercritical Hopf bifurcation produces limit cycles.
   - So in this subfamily the y^{−b} argument is sharp, not merely sufficient. I added this as a remark.
8. **Extension to c₁ ∈ {−2, −1, 0}** (V13).
   - These values are admissible, since the target (2+c₁)X + d₁Y is still a complex. For them, with d₁ ≥ 1,
     −Δ(b*)/4 = d₁κ₁(κ₂ + κ₁(d₁ − 2c₁)) > 0, so the same proof applies verbatim.
   - These networks are in case 9 and have no Hopf bifurcation, because that requires c₁ > 0. I added them to the
     remark.
   - The case d₁ = 0 is trivial: ẏ = −κ₂xy < 0 in the open quadrant, and y = 0 is invariant.
9. **Example data** (V10). The only positive equilibrium is (2κ₃/κ₂, 4κ₁κ₃/κ₂²), with tr J = −2κ₃ and
   det J = 4κ₁κ₃²/κ₂, as stated in Remark 3.

**Numerical sanity check (not part of the proof).** File `verification/verify_numeric.py` (scipy LSODA, single
thread; log `verification/verify_numeric.log`, ALL PASS).
- **Scaling.** With x = (κ₃/κ₂)X, y = (κ₃/κ₂)Y and t = T/κ₃, the subfamily depends only on L = κ₁/κ₂.
- **Part A** (d₁ ≥ 2c₁). I ran six members: (c₁,c₃,d₁) = (1,1,2) (the example), (1,1,3), (1,2,2), (2,1,4), (2,3,5)
  and (1,3,6). Each used 21 values of L in [10^−2.5, 10^2.5] and 12 random initial conditions per value, so 1512
  runs in all. Every trajectory either converged to the equilibrium or left every bounded set (|z| > 10⁷). No
  bounded non-equilibrium limit set appeared.
  - Side remark: about 20% of the initial conditions escape, so the equilibrium is not globally attracting. This
    is irrelevant to the claim.
- **Part B** (Hopf region). I ran (2,1,3), (3,1,4) and (3,2,5), where L_H = 1/(2c₁ − d₁).
  - At L = 0.90 and 0.98 L_H, the trajectories converge to the equilibrium.
  - At L = 1.01 and 1.04 L_H there is a stable periodic attractor. It is reached from inside and from outside
    with the same amplitude, agreeing to 10⁻⁷ relative.
  - The amplitude ratio is 1.91–1.98, close to the factor 2 expected for a supercritical Hopf bifurcation
    (√4 = 2).
  - This confirms that the region c₁ < d₁ < 2c₁ in Theorem 2 cannot be shrunk, and it is consistent with point 7.

## 3. The solver's scripts: what they actually check

I re-ran `scripts/certify_dulac.py` under WSL. It exits 0 and its output is byte-identical to
`certificates/certify_dulac.log`, both before and after my comment-only edits.

Its PASS lines test what they announce:
- y^{b+1}div(BF) is a quadratic form;
- the y² coefficient is −κ₂;
- the discriminant identity, and Δ = −8κ₁κ₂ for the example;
- the axis values;
- the example's equilibrium and Jacobian;
- u = c × d.

Several things are only printed or commented, not checked:
- the sign conclusion for d₁ ≥ 2c₁;
- the optimality of b;
- the specialisation of BBH's Hopf condition;
- the correspondence between the networks and the ODEs.

The last gap is why the misprint (1+c₃)X + Y went unnoticed. My script checks all four of these (V2, V4, V5, V12).

## 4. Corrections applied

1. **Third reaction of the subfamily (substantive misprint).** "Y → (1+c₃)X + Y" became "Y → c₃X + Y" in the
   `paper.tex` abstract and Theorem 2, in the `README.md` result paragraph, and in the `scripts/certify_dulac.py`
   docstring.
2. **Theorem numbering.** `paper.tex` now states that theorem numbers refer to arXiv v1, and that the classification
   is Theorem 3 in the published version, with identical cases, notation and §6 text. The same note was added to
   `README.md`, `references.md` and a comment in `scripts/certify_dulac.py`.
3. **Theorem 2, wording.**
   - c, d were called "reaction vectors". They are now described as the X- and Y-components of the reaction vectors
     in BBH's notation (13).
   - "has a periodic orbit" was made precise as "in the closed positive quadrant" / "in ℝ²≥0", matching
     Theorem 1.
4. **Proof, boundary step.** It is now a precise one-sided argument, and the proof notes that B is smooth on the open
   quadrant.
5. **Remark 3.** It now lists the sign conditions of BBH's case 9 together with the inequality, all of which hold in
   the subfamily.
6. **New Remark 4 (added at verification).** It records (i) the sharpness result, that the Δ = 0 curve is the Hopf
   curve, including the semidefinite Δ = 0 case, and (ii) the extension to c₁ ∈ {−2, −1, 0}.
7. **"What remains open".** It now says that all of case 9 outside the family c₂ = d₂ = −1, d₃ = 0, and all of case
   10, remain open. The old list was incomplete: it omitted c₂ ≠ −1, for instance.
8. **Status and credits.**
   - `\PaperSubtitle` and `\PaperCredits`: status "independently verified", with verifier credit.
   - `README.md`: status changed to "solved", plus attribution, a reproduce section and the list of files.
   - `references.md`: the verifier's novelty re-check was appended.
   - `paper.pdf` was recompiled (3 pages).
   - The **root `README.md` results table was not edited**, because that is outside my write permission; it should
     be updated by the orchestrator.

## 5. Novelty (re-checked independently)

- **Citing works.**
  - Semantic Scholar lists 12 citing records and OpenAlex 10 (W4386927497), from 2023–2026.
  - I grepped the full text of every arXiv citing work for this question: 2406.13451 (BBH, Nonlinear Dyn. 2024),
    2312.12897 (inheritance of bifurcations), 2406.05057 (Craciun–Erban), 2610.06130 (Yeung–Erban, Oct 2026),
    2608.18509, 2608.01116, 2507.09272, 2503.00008, 2406.09514, 2402.18523 and 2307.04186. BBH appears only in
    their bibliographies or introductions. 2406.13451 mentions cases 8–10 only to classify its 104 Hopf networks.
  - The three non-arXiv citing works (Wnt pathway, delayed CRN, a thesis on 2-reaction networks) were checked by
    title.
- **Later work of the authors.**
  - arXiv author listings up to Oct 2026, including Banaji 2609.07265 (fold and cusp computations) and Boros et al.
    2510.03621 (disguised toric locus).
  - B. Boros's publication page, and his GitHub `balazsboros/reaction_networks` (`3reactions/` holds only the
    focal-value notebook; last commits Sept 2026).
  - Boros's LCC2025 abstract and Hofbauer's January 2026 slides (S-systems).
  - None of these addresses the question.
- **AI campaigns.**
  - `mmislan`: 23 Palomar repositories, none on planar three-reaction networks.
  - `ipitchford`: 166 Evidence Press repositories. I filtered names and descriptions and read the READMEs of the two
    CRN-related ones (`reverse-reaction-attracting-oscillations`, 18 Oct 2026, on Feliu–Kaihnsa;
    `quadratic-mass-action-sharpness`, on equilibrium counts). Nothing on BBH §6.
  - GitHub code search for `sapm.12639` and `2304.02303` found only the authors' repository.
- **Zenodo.** I searched with required-term queries (authors, "three-reaction", Dulac + reaction) and with
  free-text queries; there were no relevant records.
- **Web.** I searched the key phrases and the verbatim sentence, with no prior answer.

**Conclusion: no prior solution found.** The answer for the example and for the subfamily appears to be new.

## 6. What I only read, without recomputing

- BBH's supercriticality of the Hopf bifurcation (first focal value, computed in their Mathematica notebook). This
  is used only for "limit cycles do occur" in the Hopf region, and my numerics agree with it.
- BBH Lemma 2 (unique positive equilibrium for nondegenerate dynamically nontrivial networks), which is not needed by
  the proof.
- The page numbers of the quote in the journal version (only the text was compared).

## 7. Reservations

- None on the correctness of Theorems 1 and 2.
- The broader BBH question (all non-Hopf networks in cases 9 and 10) remains open, and the manuscript does not claim
  it. "Solved" refers to the example asked about explicitly, plus the subfamily.
- The novelty search cannot exclude unindexed or very recent work, such as conference talks or papers in
  preparation by the authors.

## Files written by the verifier

- `verification.md` (this report).
- `verification/verify_symbolic.py` and `verification/verify_symbolic.log`: exact symbolic checks V1–V13, a few
  seconds.
- `verification/verify_numeric.py` and `verification/verify_numeric.log`: numerical sanity check, a few minutes, not
  a proof.
