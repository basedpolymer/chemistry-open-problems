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
- Result 9 independently verified: CONFIRMED WITH CORRECTIONS. The verifier wrote its own certificate checker from
  the mathematical specification before reading the solver's; all five certificates pass (K = 18 in 151 s), with the
  same exact ratios; tampered certificates rejected; Pönitz–Tittmann table reproduced by an independent
  implementation; corrections editorial (certificate-format documentation, a file size in SHA256SUMS).
- Result 10 (combustion theory): the Korman–Li–Ouyang conjecture on the perturbed Gelfand equation in the unit
  disk is proved: a unique eps_0 = 0.2421061655952..., exactly S-shaped curve below, monotone above. Du–Lou change of
  variables (eps becomes a stopping level), analytic Wronskian argument for small amplitudes, Arb cover of a compact
  interval, Du–Lou/KLO qualitative theorems for the tail. Solver GPT-6-Astra via Codex CLI (xhigh), two runs.
  Independently verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS, with an independent Julia cover;
  the published estimate eps_0 ~ 0.2385 of KLO is wrong for the disk (close to the 3D value 0.2388).

## 20 October 2026

- Result 11 (aggregation kinetics): Jang–Tran Conjecture 2 on the discrete coagulation–fragmentation equations with
  multiplicative coagulation and constant fragmentation proved, both parts. The Hamilton–Jacobi equation of the
  Bernstein transform has zeroth-order coefficient 1/τ(x), τ = 2 tanh(x/2), so Van's half-slope invariant transfers
  to the discrete setting; global classical solutions by characteristics for every mass 0 < m ≤ 1 without a second
  moment, comparison with positive truncations, part (i) through the triangular fixed-mass system. Solver Claude
  Opus 5.5 (xhigh); novelty audit, singularity analysis and first draft GPT-6-Astra via Codex CLI (xhigh).
  Independently verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS (scope of the Conjecture 1
  consequence stated precisely, a citation fixed, the long-time numerical check repaired: 59 PASS).
- Result 12 (chemical engineering): the Halvorsen–Skogestad conjecture (2003) proved for every number of
  components in a stated class and model: networks of two-section columns with interval streams (containing
  Agrawal's regular-column space: 8 / 152 / 6128 configurations for N = 3 / 4 / 5), strict infinite-stage Underwood
  model, no reverse flows, any feed quality. Exact cut identity for the Underwood flux, column sign rules and a
  sign-graph reachability argument; the Petlyuk value (highest V_min peak, previously proved only for the extreme
  peaks) is attained by the fully thermally coupled arrangement. In the permissive (TG-literal) model: proved for
  N <= 4, for extreme peaks and for merge-free mixed networks. Solvers: three Claude Opus 5.5 (xhigh) agents
  (theory, computation, general N). Independently verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS
  (degenerate operating points and further exclusions stated, TG equation numbering, computer-verification section
  aligned with the actual runs); the verifier also found an interval-certified point of the permissive model below
  the Petlyuk value for N = 5 (q = 0), an artifact of that formulation rather than a physical counterexample.
- Result 12, version 2 (second wave, four Claude Opus 5.5 agents, plus GPT-6-Astra and a literature agent): a
  stage-level argument proves the Halvorsen–Skogestad lower bound for every adiabatic arrangement. Every unit of a
  constant-molar-overflow equilibrium-stage network conserves the Underwood flux phi(theta_j) up to its heat duty, the
  feed contributes zero at a feed root, and all streams leaving an equilibrium stage have the same sign of phi;
  summing the balances over the units of one sign gives R >= H_j for every arrangement (splitters, recycles,
  multi-feed columns, reverse flows, Murphree stages, finite stages and impure products up to O(eps)). Found
  independently by three agents; an early independent verifier (Claude Opus 5.5, xhigh) confirmed the bound and its
  scope, with editorial corrections. Conjecture fully proved for N = 2 (finite-stage attainment by a single column);
  for N >= 3, attainment by the Petlyuk is proved in Underwood's infinite-stage model and observed numerically at
  finite stages (geometric decay of impurities), and the isolated prefractionator is proved for any N (GPT-6-Astra);
  finite-stage attainment of the closed network for N >= 3 remains open. Targeted novelty search: no prior proof or
  argument found (search shortened). Full verification of version 2 pending.

## 24 October 2026

All entries of this day are preliminary: no independent verifier has checked them yet. Agents: Claude Opus 5.5
sub-agents of the orchestrator (effort not recorded); credits per problem in each folder.

- Swelling exponent of self-avoiding walks on Z^3 (owner's selection of 22 October): first wave of four agents
  (literature, lower bound, upper bound, computation lab). Literature: no rung beyond Madras 2014
  (E_N|w_N|^2 >= N^{4/9}/36) and Duminil-Copin–Hammond 2013 (o(n^2)) found for Z^3; the "ICM 2026" reference of the
  owner's report (doi:10.1137/25M180679X) is Hutchcroft's long-range percolation survey, and Guttmann–Jensen 2022 is
  2D numerics. Upper bound: unconditional Theorem A (Krachun–Panagiotis route moved to Z^3 with an arch lemma) and,
  under sum_T B_T(z_c)^2 < infinity, E_n|w_n|^2 = O(n^2/(log n)^2); exact obstruction to quantifying
  Duminil-Copin–Hammond. Lower bound: multiplicity structure of Madras's reflection and a conditional theorem (a
  sublinear multiplicity bound gives an exponent above 4/9). Lab: exact enumeration with 48-fold symmetry to n = 19;
  several candidate lemmas refuted exactly (E_16[s_1.s_16] < 0). Owner supplied Kesten 1963/1964, Madras 1988/1995,
  Chayes–Chayes 1986, Hammersley–Welsh 1962 and Madras–Slade chapters 3, 4, 7, 8 and selected pages.
- Result 13 (reaction networks), solved (preliminary): Banaji–Boros–Hofbauer's question on cases 9 and 10 answered
  (no periodic orbit without Andronov–Hopf bifurcation), generalizing result 2. Verification pending.
- Started: Becker–Döring with atomization for every M (extension of result 8), knotting probability of equilateral
  hexagons, minimum number of matches in heat-exchanger networks, skewness of Taylor dispersion in regular-polygon
  ducts (owner supplied Aminian–Camassa–McLaughlin 2018), Frank–Kasper tilings in the diblock foam model, Kirchhoff
  analogue of Šoltés's problem, spectral conjecture on fullerene duals. Current states in the root README.
- Result 14 (chemical engineering, transport), solved (preliminary): Aminian–Camassa–McLaughlin's claim on the sign
  of solute skewness in regular-polygon ducts proved for long and short times (computer-assisted); the claim "at all
  times" holds on certified windows only. Verification pending.
- Result 8 extended (preliminary): Becker–Döring with atomization, first branch supercritical with stable cycles for
  every M ≥ 16 (subcritical at M = 15); branches 2 and 3 supercritical on the centre manifold, orbits unstable in the
  full system. Interval sweeps in ε = M^{-1/2} covering every M ≥ 161 at once. Not yet verified by an agent.
- SAW computation lab finished: exact enumeration to n = 19, testbed and candidate lemmas (report written by the
  orchestrator in `computation/report.md`), pivot estimate ν = 0.5876 ± 0.0005 (numerics only).
