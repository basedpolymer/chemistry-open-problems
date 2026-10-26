# Completion and verification of result 8: branches 1–3

Verifier and completion editor: Codex (GPT-6, effort not recorded). Date: 26 October 2026.
Original all-M solver: Claude Opus 5.5 (effort xhigh).

## Verdict: CONFIRMED WITH CORRECTIONS

The local criticality and stability question is resolved for every integer cluster size on the three
specified branches. Branch 1 is subcritical at M = 15 and supercritical with orbitally asymptotically
stable small periodic orbits on the mass hyperplane for every M ≥ 16. Branches 2 and 3 are supercritical
on the centre manifold for every M ≥ 44 and M ≥ 87, but their small periodic orbits are unstable in the
mass hyperplane. Their exact integer onsets are 15, 44, 87. For 3 ≤ M ≤ 14 there is no Hopf at any K > 0.

This report audits the prior all-M deliverables and supplies corrections and completion certificates.
I used the original statement, primary source, manuscripts, code and public certificate files; I did
not use the original solver's private reasoning or progress notes. The new onset checks are verifier
additions using the already available onset-chart machinery. They have not been reviewed by a further
agent or a human. The original finite-range report in `../verification.md` remains unchanged.

## Mathematical audit

1. **Model and restriction.** The mass row annihilates both B and S. The constant equilibrium has a simple
   zero eigenvalue transverse to the mass hyperplane; the restriction removes it. The quadratic form is
   B(x,y) = x₁Sy + y₁Sx, and the cubic form vanishes. Substitution in the standard projection formula
   gives the five-term transfer expression in AM5 of the manuscript. The zero-frequency inverse must
   be taken in the mass hyperplane.
2. **Regular formula.** The characteristic determinant and its epsilon-regular version are algebraically
   identical. The eigenvector, adjoint and resolvent atoms satisfy the bulk recurrences and boundary/mass
   equations. Divided differences remove the boundary-layer cancellation before interval evaluation.
   The geometric and mass sums follow from finite geometric series and differentiation. Exponential
   remainders are bounded uniformly, including epsilon = 0. In the mean-value evaluation the disks are
   held fixed and reused at the midpoint; differentiating their numerical bounds would be incorrect.
3. **Existence and branch identity.** The parametric Krawczyk inclusion is justified with the integral
   of the *preconditioned* parameter derivative. This avoids an invalid vector mean-value argument
   with different componentwise intermediate points. Every adjacent epsilon-box pair is linked by a
   newly certified endpoint zero contained in both boxes. The zero at epsilon = 0 is identified in
   the correct interval between consecutive tangent poles.
4. **Spectrum.** The explicit localization inequalities exclude the entire unbounded exterior.
   The region bounds and monotonicity checks reduce H2 to polygon counts. Counts are four on branch 1,
   and (6,4), (8,4) for the large and narrow polygons of branches 2, 3. Their extra zeros correspond
   to positive-real-part eigenvalues, sufficient for orbital instability; their simplicity and an
   algebraic-multiplicity identity are not assumed.
5. **Onsets.** The chart s = 1/(1+K) is regular at s = 0. Its s-parametrized part has dM/ds > 0;
   its M-parametrized part extends to 162 and meets the epsilon chart. Every stored box and every
   junction is re-proved, and every finite integer point through 161 is linked. Nonsingular local
   systems and escape at both ends make the resulting graph the entire specified connected component.
   This closes the previously open onset statement for branch 3.
6. **Limit.** Direct symbolic evaluation of the regular coefficient at epsilon = 0 yields the claimed
   rational expression in the Hopf parameter. Exact real-root isolation verifies all normalization
   and resolvent divisors, including their denominators. The numerical consistency diagnostic S1 is
   distinguished from the exact algebraic claims S2–S4. The discrete limit is therefore proved;
   the separate formal renewal-equation derivation is not used.
7. **Hopf conclusion.** Positive crossing and negative coefficient give the supercritical side and
   attraction on the centre manifold. The transverse spectrum makes branch 1 stable and branches 2–3
   unstable. At M = 15 the positive coefficient gives the subcritical exception. All claims concern
   sufficiently small positive-concentration cycles near the Hopf point.

## Corrections made

- Replaced floating-point `atan2` accumulation and tolerance-based integer rounding with Arb arguments
  of endpoint ratios and a unique enclosed integer. Every accepted edge image lies in a certified
  rotated right half-plane, so no branch-cut choice is guessed.
- Certified the geometry used by the contour argument: exact-ball supporting half-plane checks and
  an outward left chord for branch 1; increasing ordinates, exact top/bottom endpoints and containment
  in the monotonicity rectangle for the level-curve contours.
- Replaced floating-point coverage comparisons with exact rational comparisons. The stored endpoint is
  **0.078811099999999995**, slightly below the old rounded label 0.0788111, but strictly above 1/sqrt(161).
  A check requesting the larger rounded endpoint correctly failed even though all mathematical rows
  passed. Correcting the requested endpoint and rerunning all three coefficient checks gave PASS.
  No interval or mathematical claim about integer M was weakened.
- Added explicit endpoint links. A mere overlap of two boxes does not identify their enclosed roots.
  At the onset chart transition, a separate narrow point is used: the first full M-box need not fit
  inside the last s-box even when their curves meet.
- Added exact Hurwitz certificates for M = 3,…,14, exact onset certificates for branches 1–3, and the
  symbolic discrete-limit identity. The symbolic divisor checker now checks denominators as well as
  numerators; its sign threshold comparison uses exact algebraic numbers.
- Tightened the overlap-file audit to reject missing or duplicate M and non-PASS records. The
  coefficient comparisons printed from old decimal CSV values are consistency diagnostics; the
  coefficient proof comes from the fresh ball evaluations.
- Replaced stale partial-result claims in both READMEs and the manuscript. The full all-M proof is now
  integrated in `../paper.tex` through `proof.tex`, with a compiled, visually inspected PDF.

## Exhaustive rechecks

Environment: Fedora 44 under WSL; Python **3.12.15**, python-flint **0.9.0**, SymPy **1.14.0**, NumPy and
mpmath. Arb precision 192 bits (256 bits in limiting evaluations). The runner uses separate processes
because the regular-coefficient library has mutable differentiation context. No sub-agents were used.

| Certificate/check | Scope | Outcome |
|---|---|---|
| Direct reduced matrices and rigorous eigensolver | All 337 points: j=1, M=15…160; j=2, 44…160; j=3, 87…160 | PASS |
| Compact coefficient/Hopf cover | All 562 + 589 + 834 = 1,985 boxes, including epsilon=0 | PASS |
| Entire spectrum, branch 1 | All 562 boxes; count 4, localization and region | PASS |
| Entire spectrum, branches 2–3 | All 589 and 834 boxes; counts (6,4) and (8,4) | PASS |
| Epsilon branch links | All 561 + 588 + 833 = 1,982 junctions, exact coverage and limiting labels | PASS |
| Onset curves | 3,075 + 3,171 + 2,368 = 8,614 boxes; every junction | PASS |
| Integer-to-onset links | 147 + 118 + 75 = 340 points, through M=161 | PASS |
| Finite-range overlap | All 99,840 points per branch for 161 ≤ M ≤ 100,000 | PASS |
| Small-size exclusion | Every K > 0 for M=3,…,14, exact Hurwitz polynomials | PASS |
| Formula checks F1–F7, O1–O4 | Symbolic identities; dense systems and direct sums | PASS |
| Limit checks S1–S4 and A1–A4 | Exact symbolic identity/divisors/sign; certified limiting evaluations | PASS |
| Winding regression | Multiple roots, reversed orientation, no roots, boundary-zero and open-path rejection | PASS |

The complete reproducible command is `python scripts/recheck_all.py` from this folder. The archived
current logs, regenerated onset/spectral files, dense certificates, execution summary and SHA-256
manifest are in `verification/recheck-2026-10-25/` (runs on 25–26 October). The summary aggregates all
21 required successful jobs and records their run provenance, including the reruns after correcting
the coverage target; it does not claim a single invocation. Historical failed attempts and their
corrections are retained in `earlier-attempts/`. The manifest hashes text with LF line endings so that
Git checkout conventions do not change the digests. These checks prove the infinite
range with parameter intervals, not samples of large integers.

## Scope and remaining questions

The result is complete for local criticality and orbital stability on the three specified connected
branches. It does not assert that no other disconnected Hopf components exist, that the equilibrium
is stable for every 0 < K < K₁(M), or that distant periodic orbits are stable. Higher branches retain
the status of result 25. The trust base is standard finite-dimensional bifurcation theory, the
Krawczyk and argument-principle theorems, Arb and exact polynomial arithmetic; this is not a proof
formalized in a proof assistant or a human-refereed result.
