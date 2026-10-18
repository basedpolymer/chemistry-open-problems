# History

## 18 October 2026

- Repository created (private).
- Results 1–3 (reaction networks): general family and proof for Banaji–Feliu Remark 4.23; Dulac proof and
  subfamily classification for the Banaji–Boros–Hofbauer example; rigorous Arb certificate and a 10-reaction
  integer-rate network for Craciun–Erban. The n = 3 case of result 1, the Dulac function of result 2 and the
  17-reaction network of result 3 come from the exploration phase of 17 October 2026.
- Result 4 (polymer physics): DGH charged di-block polymer, solved by GPT-6-Astra (Codex CLI).
- Result 5 (reaction networks): Boros–Hofbauer four-limit-cycle conjecture proved (focal values + Krawczyk, and
  validated Taylor integration of the return map), solver sub-agent Claude Opus 5.5 (xhigh).
- Result 6 (chemical engineering): Jäschke–Skogestad temperature conjecture disproved for M_j ≥ 3 and corrected
  recursion proved, solver sub-agent Claude Opus 5.5 (xhigh).
- Independent verification started (Claude Opus 5.5, xhigh) for results 5 and 6.
- Result 6 independently verified: CONFIRMED WITH CORRECTIONS (wording of the Appendix B claim softened, M ≥ 4
  argument completed, three further exact counterexamples including M = 3 with constant hot temperatures).
- Under investigation: Becker–Döring equations with atomization (Hopf criticality); upper bound for the
  connective constant of self-avoiding walks on the cubic lattice.
- Result 5 independently verified: CONFIRMED WITH CORRECTIONS (focal values recomputed by another algorithm,
  Krawczyk redone with another interval library, uniqueness re-checked with msolve, independent validated integrator;
  corrections to the abstract and to a remark).
- Result 7 (reaction networks): Banaji–Feliu Remark 5.8 answered for every n by an explicit family (orchestrator,
  Claude Opus 5.5 xhigh); Remark 5.6 of the same paper found to misquote its source for n = 1.
- Result 8 (nucleation kinetics, partial): Becker–Döring with atomization, sign of the first Lyapunov coefficient
  certified for 15 ≤ M ≤ 10^5 (supercritical except M = 15), solver sub-agent Claude Opus 5.5 (xhigh).
- Result 3 independently verified: CONFIRMED WITH CORRECTIONS (weak reversibility edge by edge, invariance,
  absence of equilibria by Sturm and an exact lower bound, divergence integral by two other rigorous methods;
  hypotheses of the network definition (2.1) stated and proof details added).
- Result 1 independently verified: CONFIRMED WITH CORRECTIONS (equilibria recounted by Gröbner bases for
  n = 3..6; the displayed instance of N_4 in the manuscript had wrong signs and was replaced, N_5 added; check count
  corrected to 300; minimality statement made precise for n = 2).
- Result 7 independently verified: proof CONFIRMED WITH CORRECTIONS, but the verifier found prior work: an
  anonymous unrefereed manuscript (Zenodo, GitHub ipitchford, 4 Oct 2026) answering Remark 5.8 for every n >= 2,
  built on the Phillipson–Rojas (2013) extremal fewnomial system, which is also the equilibrium system of our family.
  Status changed to "already solved". The novelty check in PROTOCOL.md now includes Evidence Press and Zenodo; a scan of
  its 166 repositories found no overlap with the other results.
- Repository made public; history scrubbed of local paths; licenses added (MIT for code, CC BY 4.0 for texts and data).
- Result 2 independently verified: CONFIRMED WITH CORRECTIONS (typo in the third reaction of the subfamily fixed:
  Y -> c3 X + Y; the Dulac threshold shown to coincide with the Hopf curve, hence optimal; theorem numbering of the
  published version noted).
- Result 4 independently verified: CONFIRMED WITH CORRECTIONS (proof re-derived line by line; own transfer-matrix
  and dynamic-programming checks; the free energy has a first-order kink at beta_*; analytic proof that
  beta_* > 2 pi / 9 by Poisson summation; DGH 1992 only predicted the value by a one-line free-energy argument).
- Result 9 (polymer physics, record): rigorous upper bound mu(Z^3) <= 4.722048 for the connective constant of
  self-avoiding walks, improving 4.7387 (Pönitz–Tittmann 2000), plus new bounds for d = 4, 5, 6; optimised
  forbidden-set automata with exact Collatz–Wielandt certificates (largest: 2.06e8 states), solver sub-agent
  Claude Opus 5.5 (xhigh). Independent verification pending.
- Result 8 independently verified: CONFIRMED WITH CORRECTIONS. Signs of l1 re-certified in other coordinates
  for all M <= 200 and samples up to 1e5, subcritical M = 15 and supercritical M = 16 confirmed by periodic-orbit
  shooting; the stability claim was overstated and is now rigorous only for 16 <= M <= 160 (spectral hypothesis
  certified there, plus M = 250, 300) and conditional beyond; "first crossing" proved for 15 <= M <= 57.
