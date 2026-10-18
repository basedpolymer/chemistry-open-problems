Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Independent verification: a cubic weakly reversible system with an algebraic limit cycle (W^a(3) ≥ 1)

Date: 18 Oct 2026. Inputs used: the task statement (claim and network), the primary source, the deliverables
(`paper.tex`, `paper.pdf`, `README.md`, `references.md`, `scripts/certify_cycle.py`, `certificates/`). Not read:
`PROGRESS.md`, solver logs, `travail/`, `rapports/`.

## Verdict: **CONFIRMED WITH CORRECTIONS**

The network of ten reactions belongs to the class W_3 of Craciun and Erban, and its mass-action system has the circle
(x−2)² + (y−2)² = 1 as a hyperbolic, orbitally asymptotically stable, algebraic limit cycle. So W^a(3) ≥ 1, which
answers the published open question in the affirmative. I found no mathematical error. The corrections are editorial
or add precision (listed below); none changes the theorem or the proof strategy.

## 1. Source and statement

- I downloaded arXiv:2406.05057v2 myself and read §2 (Definitions 1–3, the standing assumptions on networks (2.1),
  condition (2.10)), §3, §4 (Definition 5, Lemmas 9–12, Table 2), §5 (Theorem 2) and §6. I also checked the journal
  version (PMC12098472, open access) for the same sentences, Table 2 and Definitions 2, 3 and 5, and they match.
- Both quoted sentences are verbatim: §4, just before Lemma 11 ("The question of whether cubic weakly reversible
  two-species systems can give rise to algebraic limit cycles remains open, leaving us with inequality W^(a)(3) ≥ 0 in
  Table 2."), and §6 ("it remains an open question whether a cubic weakly reversible system can have an algebraic limit
  cycle."). Table 2 gives W^a(3) ≥ 0.
- Definitions, checked against the source:
  - Network (2.1): reactions α X + β Y → γ X + δ Y with positive rate constants; no four-tuple repeated;
    (α, β) ≠ (γ, δ); both species take part in some reaction.
  - S_n: order max(α+β) ≤ n.
  - M_n (Definition 2, (2.10)): max(α+β, γ+δ) ≤ n, so reactants **and** products are bounded.
  - W_n (Definition 3): the subset of M_n that is weakly reversible, meaning that every edge of the E-graph lies on an
    oriented cycle.
  - W^a(n) (Definition 5): the maximum number of algebraic limit cycles of ODEs in W_n. An algebraic limit cycle (§4)
    is a closed isolated trajectory that is a closed component of an algebraic curve h = 0.
  - The definitions put no constraint on the number of reactions or on which complexes are allowed beyond
    molecularity ≤ n, and they do not require a "generic" cycle. "Cubic" in the open sentence means n = 3, i.e. W_3.
    The claim therefore answers exactly the question that was asked, not a weaker variant.

## 2. What I recomputed myself (own code: `verification/verify_independent.py`, log `verification/verify_independent.log`, ALL PASS, ~3 min)

The reaction list was retyped from the task statement; nothing was copied from the solver's script.

| Step | Method (independent of the solver's) | Result |
|---|---|---|
| Membership in M_3 | molecularity of every source **and** every product complex | all ≤ 3; order 3 (cubic) |
| Assumptions of (2.1) | no trivial reaction, no repeated reaction, both species change | hold |
| Weak reversibility | **reaction by reaction**: edge u→v lies on a cycle iff u is reachable from v | all 10 edges on cycles; 8 complexes, 1 linkage class, rank 2, deficiency 5 |
| Mass-action field | rebuilt from the reaction list, compared with the paper | identical f, g (also checked by hand) |
| Invariance | (i) multivariate division of f h_x + g h_y by h; (ii) exact substitution of the rational parametrisation x = 2+(1−t²)/(1+t²), y = 2+2t/(1+t²) | remainder 0; cofactor s identical to the paper's; identity vanishes identically on the parametrisation |
| No equilibrium on C | the tangential speed P = −(y−2)f + (x−2)g on C is a positive multiple of a polynomial in t = tan(θ/2): (1+t²)⁴P = 4(1+t²)²(130t⁴+135t³+522t²+329t+198); exact Sturm count over all of ℝ, plus P(1,2) = 520 at t = ∞ | no real root, so P > 0 on C; exact bound P > 18759/50 = 375.18 (sampled minimum 375.1907). Since F·∇h = 0 on C, F = P·(unit tangent), so F ≠ 0 on C |
| ∮ div F dt | **E1**: rigorous mpmath interval arithmetic (not Arb), composite midpoint rule on 2048 cells with the remainder h³/24·φ''(cell) enclosed by interval evaluation of the symbolic second derivative | [−3.6847574, −3.6847469] < 0 |
| ∮ div F dt | **E2**: closed form by residues: the integral of a rational function over ℝ (after t = tan(θ/2) and cancellation, 3 simple poles in the upper half-plane, isolated exactly by sympy CRootOf) | −3.684751929326756730604462… (inside E1 and inside the solver's Arb ball [−3.68475192932676 ± 3.3e−15]) |
| Period | same two methods | E1 [0.01033203957, 0.01033204871]; E2 0.01033204340829214261807… |
| Dynamics (floating point, sanity only) | scipy DOP853 on the ODE and its variational equation over one period from (3, 2) | return error 3e−14; det(monodromy) = 0.025103401369 = e^Λ; eigenvalues {1, 0.0251}; trajectories started at (2.8, 2) and (3.2, 2) reach \|h\| ≈ 4e−12 after 40 periods |
| Remark "construction" | linear system "f h_x + g h_y = s h" in the 10 rates and 6 cofactor coefficients on the N10 support | 15 equations, 16 unknowns, nullity 1; the generator is the stated integer vector |
| N17 (Remark "a second example") | all the checks above | all PASS; ∮ div F dt = −0.6365396683814778…, T = 1.915981759823322…; 10 complexes, deficiency 7 |

**Stability argument.** C is a compact invariant curve without equilibria. It is therefore a periodic orbit,
traversed counter-clockwise with θ̇ = P (the radius is 1). For a planar C¹ field, the derivative of the Poincaré map of
a periodic orbit equals exp(∮ div F dt) (Perko, *Differential Equations and Dynamical Systems*, §3.4, the planar
derivative formula). This is applied correctly, with the time integral converted by dt = dθ/P. Here the derivative is
e^Λ ≈ 0.0251 < 1, so the orbit is hyperbolic and attracting. Being hyperbolic, it is isolated among periodic orbits,
so it is a limit cycle in the sense of Craciun and Erban ("closed isolated solution"). It lies in [1, 3]², inside the
open positive quadrant, and is the whole zero set of h, hence an algebraic limit cycle. The paper also calls
∮ s dt = ∮ div F dt a consistency check. That identity is correct (h is a transversal coordinate with ḣ = s h) and is
not needed for the proof.

**Side check (`verification/milp_min_support.py`, HiGHS, 1 thread, ~20 s; floating point; not used in the proof).**
Among the 90 possible reactions between the 10 complexes of molecularity ≤ 3, the minimum number of reactions of a
weakly reversible network meets the following conditions: rates in [1, 10⁵], C invariant, and |P| ≥ 1 with constant
sign at 96 points of C. That minimum is **10** for both orientations. This is consistent with the paper's statement
that minimisation returned a 10-reaction support. The optimal support is not unique: HiGHS returned a different
10-reaction network.

## 3. The solver's script: re-run and audit

`scripts/certify_cycle.py` re-ran in about 5 s under WSL (Python 3.12.15, sympy 1.14.0, python-flint 0.9.0). Its output
is byte-identical to `certificates/certify_cycle.log`, ALL PASS. Here is what each PASS actually tests:

- "maximal molecularity 3" runs over the set of all complexes, sources and products together, so it does test M_3.
  "Strongly connected" implies weak reversibility (sufficient, and true here).
- The script does **not** check the standing assumptions of (2.1) (no repeated reaction, source ≠ product). They hold
  (my A4–A6), so this is not a defect of the result.
- The invariance test is `sp.div(num, h)` with zero remainder, which is a valid exact certificate.
- The resultant test is valid: h is monic in y, so a common real zero would give a real root in [1, 3]. The degree-4
  resultant, which looks low next to the generic degree 8, is explained by P ≡ −2(97x² − 97xy − 262x − 38y + 175)
  (mod h): P restricted to C is a trigonometric polynomial of degree 2.
- The Arb positivity loop really covers [0, 2π]: ball i + 0.5 ± 0.5, scaled by 2π/N.
- The `acb.integral` integrand ignores the `analytic` flag. This is acceptable for an integrand built from field
  operations and the entire functions cos and sin: on a ball where P may vanish, ball division returns an
  indeterminate value and the integrator subdivides. The enclosures are therefore rigorous. My independent E1 and E2
  computations agree to 13+ digits.
- The stability message is based on ∮ s dt < 0. ∮ div F dt is only checked to be nonzero and to overlap ∮ s dt.
  Together these imply ∮ div F dt < 0, so the logic is sound, if indirect. My E1 checks the sign of ∮ div F dt
  directly.
- The header claimed "Runtime ~ 1 min" while the README said ~3 s. I corrected it to the measured ~5 s.

## 4. What I only read (not recomputed)

- Perko §3.4, the derivative of the planar Poincaré map: standard, cited from memory. I did not re-check the exact
  theorem number in the 3rd edition.
- Correctness of the libraries (sympy Sturm sequences and CRootOf isolation, mpmath `iv`, Arb). I reduced the risk
  by using three independent tools for the key sign (Arb, mpmath iv, exact residues) and two for positivity of P
  (resultant, and a Sturm count on the half-angle polynomial).
- The construction process in the Remark (the LP that found N17, the SCIP run): not re-run. My HiGHS MILP gives a
  consistent result.
- Lemma 11 and Theorem 2 of the source: read, not re-proved; the result does not depend on them.

## 5. Corrections applied

1. `paper.tex` / `paper.pdf` (recompiled, 4 pages):
   - Credits and subtitle: "Independent verification: pending" / "preliminary" replaced by the verifier credit and
     "independently verified".
   - Abstract: the enclosure "[−3.684751929326 ± 10⁻¹²]" was harmonised with the body and the README
     ("[−3.68475192932676 ± 3.3×10⁻¹⁵]"), with a mention of the independent confirmation.
   - Proof, step (1): the weak-reversibility argument deduced "every reaction lies on a directed cycle" from a union
     of directed paths and cycles. The deduction was correct, but it went implicitly through strong connectivity. The
     step now lists three explicit directed cycles that cover all ten reactions. The step also now states
     explicitly that every product, not only every reactant, has molecularity ≤ 3 (condition (2.10)), and that the
     standing assumptions of (2.1) hold.
   - Proof, step (3): the explicit resultant, the reduction P ≡ −2(97x² − 97xy − 262x − 38y + 175) mod h, the
     half-angle form P = 4(130t⁴ + … + 198)/(1+t²)², and the exact bound P > 375.18 were added, so a referee can check
     positivity by hand.
   - Proof, step (4): the justification of Arb's rigour was made precise. The old wording was "analytic near the
     path; P has no zero on [0, 2π]"; the new wording explains why complex poles are handled automatically. The two
     independent evaluations (interval midpoint rule, residues) were added.
   - Remark "construction": added the verifier's bounded, sampled MILP check of 10-reaction minimality, labelled as not
     used in the proof.
   - Remark "a second example": rephrased to remove an overfull line. "Computer verification": runtime 3 s → 5 s, and
     a paragraph on `verification/verify_independent.py` was added.
2. `scripts/certify_cycle.py`: the runtime in the docstring was corrected (code unchanged; output still identical to
   the certificate log).
3. `README.md`: status changed to solved (confirmed with corrections), verifier attribution added, the reproduce
   commands and file table now include the verification files, and a history line was added.

Not changed (outside this folder, which I may not edit): the results table in the repository's root `README.md`
should now show this problem as **solved** with "Independent verification: Claude Opus 5.5 (effort xhigh)".

## 6. Novelty check (quick)

- OpenAlex (journal DOI): 2 citing works, both 2026: a Bull. Math. Biol. survey and a J. Nonlinear Sci. paper on
  switching ODEs. I read only the titles; neither title suggests a cubic weakly reversible algebraic limit cycle. The
  solver's `references.md` says the same.
- arXiv listings: the latest papers of R. Erban (up to arXiv:2610.06130, Oct 2026), G. Craciun (up to 2601.08801),
  B. Boros (up to 2510.03621) and T. Plesa (up to 2602.02374). Full-text arXiv searches for "algebraic limit cycle(s)"
  with "chemical"/"reaction", and for "weakly reversible" with "limit cycle", return only Craciun–Erban itself.
- Book of abstracts of LCC2025 (UAB, where Erban and Craciun lectured): no claim of a cubic example.
- GitHub `mmislan`: 23 Palomar repositories, none on limit cycles or W^a(n).
- Conclusion: no prior or concurrent solution found. The search is not exhaustive: Google Scholar was not queried,
  and the two citing works were not read in full.

## 7. Reservations

- The result is the lower bound W^a(3) ≥ 1. This settles the published question, which asks for existence; the exact
  value of W^a(3) remains open.
- The certificate depends on computer algebra and interval arithmetic. These are standard and were cross-checked with
  three tools, but they are not formally verified.
- The MILP minimality statement is floating-point, bounded and sampled. It is explicitly not part of the theorem.
- The Perko theorem number is cited from memory.
