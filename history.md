# History

## 18 October 2026

- Repository created (private).
- Results 1–3 (reaction networks): general family and proof for Banaji–Feliu Remark 4.23; Dulac proof and
  subfamily classification for the Banaji–Boros–Hofbauer example; rigorous Arb certificate and a 10-reaction
  integer-rate network for Craciun–Erban. The n = 3 case of result 1, the Dulac function of result 2 and the
  17-reaction network of result 3 come from the exploration phase of 17 October 2026.
- Result 4 (polymer physics): DGH charged di-block polymer, solved by GPT-6-Astra.
- Result 5 (reaction networks): Boros–Hofbauer four-limit-cycle conjecture proved (focal values + Krawczyk, and
  validated Taylor integration of the return map), solver sub-agent Claude Opus 5.5 (xhigh).
- Result 6 (chemical engineering): Jäschke–Skogestad temperature conjecture disproved for M_j ≥ 3 and corrected
  recursion proved, solver sub-agent Claude Opus 5.5 (xhigh).
- Verification started (Claude Opus 5.5, xhigh) for results 5 and 6.
- Result 6 verified: CONFIRMED WITH CORRECTIONS (wording of the Appendix B claim softened, M ≥ 4
  argument completed, three further exact counterexamples including M = 3 with constant hot temperatures).
- Under investigation: Becker–Döring equations with atomization (Hopf criticality); upper bound for the
  connective constant of self-avoiding walks on the cubic lattice.
- Result 5 verified: CONFIRMED WITH CORRECTIONS (focal values recomputed by another algorithm,
  Krawczyk redone with another interval library, uniqueness re-checked with msolve, separate validated integrator;
  corrections to the abstract and to a remark).
- Result 7 (reaction networks): Banaji–Feliu Remark 5.8 answered for every n by an explicit family (orchestrator,
  Claude Opus 5.5 xhigh); Remark 5.6 of the same paper found to misquote its source for n = 1.
- Result 8 (nucleation kinetics, partial): Becker–Döring with atomization, sign of the first Lyapunov coefficient
  certified for 15 ≤ M ≤ 10^5 (supercritical except M = 15), solver sub-agent Claude Opus 5.5 (xhigh).
- Result 3 verified: CONFIRMED WITH CORRECTIONS (weak reversibility edge by edge, invariance,
  absence of equilibria by Sturm and an exact lower bound, divergence integral by two other rigorous methods;
  hypotheses of the network definition (2.1) stated and proof details added).
- Result 1 verified: CONFIRMED WITH CORRECTIONS (equilibria recounted by Gröbner bases for
  n = 3..6; the displayed instance of N_4 in the manuscript had wrong signs and was replaced, N_5 added; check count
  corrected to 300; minimality statement made precise for n = 2).
- Result 7 verified: proof CONFIRMED WITH CORRECTIONS, but the verifier found prior work: an
  anonymous unrefereed manuscript (Zenodo, GitHub ipitchford, 4 Oct 2026) answering Remark 5.8 for every n >= 2,
  built on the Phillipson–Rojas (2013) extremal fewnomial system, which is also the equilibrium system of our family.
  Status changed to "already solved". The novelty check in PROTOCOL.md now includes Evidence Press and Zenodo; a scan of
  its 166 repositories found no overlap with the other results.
- Repository made public; history scrubbed of local paths; licenses added (MIT for code, CC BY 4.0 for texts and data).
- Result 2 verified: CONFIRMED WITH CORRECTIONS (typo in the third reaction of the subfamily fixed:
  Y -> c3 X + Y; the Dulac threshold shown to coincide with the Hopf curve, hence optimal; theorem numbering of the
  published version noted).
- Result 4 verified: CONFIRMED WITH CORRECTIONS (proof re-derived line by line; own transfer-matrix
  and dynamic-programming checks; the free energy has a first-order kink at beta_*; analytic proof that
  beta_* > 2 pi / 9 by Poisson summation; DGH 1992 only predicted the value by a one-line free-energy argument).
- Result 9 (polymer physics, record): rigorous upper bound mu(Z^3) <= 4.722048 for the connective constant of
  self-avoiding walks, improving 4.7387 (Pönitz–Tittmann 2000), plus new bounds for d = 4, 5, 6; optimised
  forbidden-set automata with exact Collatz–Wielandt certificates (largest: 2.06e8 states), solver sub-agent
  Claude Opus 5.5 (xhigh). Verification pending.
- Result 8 verified: CONFIRMED WITH CORRECTIONS. Signs of l1 re-certified in other coordinates
  for all M <= 200 and samples up to 1e5, subcritical M = 15 and supercritical M = 16 confirmed by periodic-orbit
  shooting; the stability claim was overstated and is now rigorous only for 16 <= M <= 160 (spectral hypothesis
  certified there, plus M = 250, 300) and conditional beyond; "first crossing" proved for 15 <= M <= 57.
- Result 9 verified: CONFIRMED WITH CORRECTIONS. The verifier wrote its own certificate checker from
  the mathematical specification before reading the solver's; all five certificates pass (K = 18 in 151 s), with the
  same exact ratios; tampered certificates rejected; Pönitz–Tittmann table reproduced by a separate
  implementation; corrections editorial (certificate-format documentation, a file size in SHA256SUMS).
- Result 10 (combustion theory): the Korman–Li–Ouyang conjecture on the perturbed Gelfand equation in the unit
  disk is proved: a unique eps_0 = 0.2421061655952..., exactly S-shaped curve below, monotone above. Du–Lou change of
  variables (eps becomes a stopping level), analytic Wronskian argument for small amplitudes, Arb cover of a compact
  interval, Du–Lou/KLO qualitative theorems for the tail. Solver GPT-6-Astra (xhigh), two runs.
  Verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS, with a separate Julia cover;
  the published estimate eps_0 ~ 0.2385 of KLO is wrong for the disk (close to the 3D value 0.2388).

## 20 October 2026

- Result 11 (aggregation kinetics): Jang–Tran Conjecture 2 on the discrete coagulation–fragmentation equations with
  multiplicative coagulation and constant fragmentation proved, both parts. The Hamilton–Jacobi equation of the
  Bernstein transform has zeroth-order coefficient 1/τ(x), τ = 2 tanh(x/2), so Van's half-slope invariant transfers
  to the discrete setting; global classical solutions by characteristics for every mass 0 < m ≤ 1 without a second
  moment, comparison with positive truncations, part (i) through the triangular fixed-mass system. Solver Claude
  Opus 5.5 (xhigh); novelty audit, singularity analysis and first draft GPT-6-Astra (xhigh).
  Verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS (scope of the Conjecture 1
  consequence stated precisely, a citation fixed, the long-time numerical check repaired: 59 PASS).
- Result 12 (chemical engineering): the Halvorsen–Skogestad conjecture (2003) proved for every number of
  components in a stated class and model: networks of two-section columns with interval streams (containing
  Agrawal's regular-column space: 8 / 152 / 6128 configurations for N = 3 / 4 / 5), strict infinite-stage Underwood
  model, no reverse flows, any feed quality. Exact cut identity for the Underwood flux, column sign rules and a
  sign-graph reachability argument; the Petlyuk value (highest V_min peak, previously proved only for the extreme
  peaks) is attained by the fully thermally coupled arrangement. In the permissive (TG-literal) model: proved for
  N <= 4, for extreme peaks and for merge-free mixed networks. Solvers: three Claude Opus 5.5 (xhigh) agents
  (theory, computation, general N). Verified (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS
  (degenerate operating points and further exclusions stated, TG equation numbering, computer-verification section
  aligned with the actual runs); the verifier also found an interval-certified point of the permissive model below
  the Petlyuk value for N = 5 (q = 0), an artifact of that formulation rather than a physical counterexample.
- Result 12, version 2 (second wave, four Claude Opus 5.5 agents, plus GPT-6-Astra and a literature agent): a
  stage-level argument proves the Halvorsen–Skogestad lower bound for every adiabatic arrangement. Every unit of a
  constant-molar-overflow equilibrium-stage network conserves the Underwood flux phi(theta_j) up to its heat duty, the
  feed contributes zero at a feed root, and all streams leaving an equilibrium stage have the same sign of phi;
  summing the balances over the units of one sign gives R >= H_j for every arrangement (splitters, recycles,
  multi-feed columns, reverse flows, Murphree stages, finite stages and impure products up to O(eps)). Found
  separately by three agents; an early verifier (Claude Opus 5.5, xhigh) confirmed the bound and its
  scope, with editorial corrections. Conjecture fully proved for N = 2 (finite-stage attainment by a single column);
  for N >= 3, attainment by the Petlyuk is proved in Underwood's infinite-stage model and observed numerically at
  finite stages (geometric decay of impurities), and the isolated prefractionator is proved for any N (GPT-6-Astra);
  finite-stage attainment of the closed network for N >= 3 remains open. Targeted novelty search: no prior proof or
  argument found (search shortened). Full verification of version 2 pending.

## 24 October 2026

All entries of this day are preliminary: no verifier has checked them yet. Agents: Claude Opus 5.5
sub-agents of the orchestrator (effort xhigh); credits per problem in each folder.

- Swelling exponent of self-avoiding walks on Z^3 (selected on 22 October): first wave of four agents
  (literature, lower bound, upper bound, computation lab). Literature: no rung beyond Madras 2014
  (E_N|w_N|^2 >= N^{4/9}/36) and Duminil-Copin–Hammond 2013 (o(n^2)) found for Z^3; the "ICM 2026" reference of the
  scouting report (doi:10.1137/25M180679X) is Hutchcroft's long-range percolation survey, and Guttmann–Jensen 2022 is
  2D numerics. Upper bound: unconditional Theorem A (Krachun–Panagiotis route moved to Z^3 with an arch lemma) and,
  under sum_T B_T(z_c)^2 < infinity, E_n|w_n|^2 = O(n^2/(log n)^2); exact obstruction to quantifying
  Duminil-Copin–Hammond. Lower bound: multiplicity structure of Madras's reflection and a conditional theorem (a
  sublinear multiplicity bound gives an exponent above 4/9). Lab: exact enumeration with 48-fold symmetry to n = 19;
  several candidate lemmas refuted exactly (E_16[s_1.s_16] < 0). Sources read: Kesten 1963/1964, Madras 1988/1995,
  Chayes–Chayes 1986, Hammersley–Welsh 1962 and Madras–Slade chapters 3, 4, 7, 8 and selected pages.
- Result 13 (reaction networks), solved (preliminary): Banaji–Boros–Hofbauer's question on cases 9 and 10 answered
  (no periodic orbit without Andronov–Hopf bifurcation), generalizing result 2. Verification pending.
- Started: Becker–Döring with atomization for every M (extension of result 8), knotting probability of equilateral
  hexagons, minimum number of matches in heat-exchanger networks, skewness of Taylor dispersion in regular-polygon
  ducts, Frank–Kasper tilings in the diblock foam model, Kirchhoff
  analogue of Šoltés's problem, spectral conjecture on fullerene duals. Current states in the root README.
- Result 14 (chemical engineering, transport), solved (preliminary): Aminian–Camassa–McLaughlin's claim on the sign
  of solute skewness in regular-polygon ducts proved for long and short times (computer-assisted); the claim "at all
  times" holds on certified windows only. Verification pending.
- Result 8 extended (preliminary): Becker–Döring with atomization, first branch supercritical with stable cycles for
  every M ≥ 16 (subcritical at M = 15); branches 2 and 3 supercritical on the centre manifold, orbits unstable in the
  full system. Interval sweeps in ε = M^{-1/2} covering every M ≥ 161 at once. Not yet verified by an agent.
- SAW computation lab finished: exact enumeration to n = 19, testbed and candidate lemmas (report written by the
  orchestrator in `computation/report.md`), pivot estimate ν = 0.5876 ± 0.0005 (numerics only).
- Work paused (about 12:30 UTC). Final states, none verified by an agent yet:
  hexagon knotting probability certified in (3.0675e-5, 8.4945e-4), Hake's Lemmas 3.9 and 3.12 false; heat-exchanger
  instances 37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22 closed exactly, lower bound 32 for the balanced15 family,
  unbalanced10 = 25 by exact SCIP (certificate not checked); Frank–Kasper DFM screen of 124 structures with relaxed
  boxes (sigma lowest; published values of BCC, P, M, H, mu shown not minimal); Kirchhoff–Šoltés: no new graph,
  exhaustive to n <= 11, cyclomatic number <= 5 and most vertex-transitive graphs with n <= 47, class theorems;
  fullerene duals: arXiv-v1 Conjecture 1(a) already refuted by the authors, published Conjecture 1(b) proved in two
  residue classes mod 24, the authors' 2025 Conjecture 3 refuted for alpha = beta. Verification of results 13, 14 and
  of the extension of result 8 was not started.

### Frank–Kasper continuation

GPT-6-Astra (OpenAI; effort xhigh), one agent, no sub-agents. Restored the
Zenodo catalogue (MD5 verified) and read the 2010/2012 articles plus the PNAS
main article. The PNAS SI was subsequently retrieved
via Europe PMC; all twelve unconstrained DFM energies and caption ratios
agree with the preprint inputs.

- Original six certificates, catalogue check and Milestone 1 replayed successfully.
- Added a separate exact energy/topology audit, including rigorous affine
  transport to the intended boxes: 6 PASS. Saved standalone integer topology
  certificates for all 121 catalogue entries: 121 PASS. Six adversarial and
  invariance regression tests pass.
- Added two symmetry-breaking restarts for σ, E11, H and μ through levels
  0–2; all 24 stages pass the floating-point geometry check and keep the
  observed ordering. These are numerical evidence only.
- Proved the universal lower bound F ≥ 1, a limitation of isolated-cell
  graph bounds, and a certified lower bound F > 1.0699 for the one-parameter
  flat A15 formula. Its finite-energy cell-collapse limit illustrates why
  degenerating sequences need separate analysis.
- Added the proof manuscript and same-agent verification record. Corrected
  claims equating floating-point values with exact certificates or failed
  relaxations with energetic exclusions. **Global σ optimality remains
  unproved**; the catalogue infimum is enclosed in [1, 1.068575909].

### Frank–Kasper continuation, 25 October 2026

GPT-6-Astra (OpenAI; effort xhigh), without sub-agents. Derived a uniform
DFM lower bound F > 1.00000025 using a quantitative moment defect and the
Kepler packing theorem, with exact rational constants. The fully reproduced
elementary Blichfeldt argument gives F > 1.000000004. These bounds are uniform
in periods, positive volumes and box shapes, but remain far below σ.

Reduced all homogeneous affine deformations of any fixed partition to convex
minimization with a global supporting-plane certificate. Exact tensor
integration and Arb verification enclose the six stored affine-family infima;
all five competitor families lie above the σ witness. Three analytic tests
pass. This resolves these affine subproblems, not the unrestricted topology
comparison. Manuscript, source record, READMEs and progress journal updated;
fresh review remains pending and global σ optimality is unproved.

### Frank–Kasper: topology-wide bounds and auxiliary-apex witnesses, 25 October 2026

GPT-6-Astra (OpenAI; effort xhigh), without sub-agents. Constructed and
exactly checked linear divergence fields for all 121 topology classes,
covering all internal shapes and box strains; 34 improve the uniform
bound. Proved that this construction, with the separate ball-moment bound,
has a ceiling below 1.066 across the catalogue. It therefore cannot separate
competitors from the current σ witness. This ceiling limits the proof
method, not physical energies.

Completed 48 positive-cone barrier stages for the eight previously unrelaxed
entries and certified eight improved dyadic witnesses. Y04, Y05 and Y06 are
embedded despite negative cones at true volume centroids, giving exact
counterexamples to necessity of that diagnostic. Corrected its documentation.
Added four linear-calibration tests and three barrier/centroid tests; the
existing six audit tests still pass. Updated manuscript, READMEs, audit and
progress journal. Global topology ordering remains unresolved.

## 25 October 2026

- Credits corrected: the Codex model of the Frank–Kasper continuation and of the Brownian
  coagulation selection is GPT-6-Astra (previously written "GPT-6"); every GPT-6-Astra run, including the internal
  reviewer sub-agent of result 4, and the Claude Opus 5.5 agents of 24 October ran at effort xhigh (previously
  "effort not recorded"). Mentions of the Codex CLI dropped from the credits.
- Root README: the Status column now links to each manuscript (or to the result file when there is no manuscript
  yet) and gives the date of the result; new sections for the work in progress and for candidate problems.
- Started (Claude Opus 5.5, xhigh): bimolecular (2, 4, 2) networks with several stable equilibria
  (Banaji–Boros–Hofbauer, Nonlinear Dyn. 2024), VIPR check of the unbalanced10 closure, infinite families of
  β-Kirchhoff Šoltés graphs (Gottwald–Majstorović Ergotić–Došlić 2026, Problem 4); verification of
  result 13.
- Result 15 (chemical engineering, heat integration), solved (preliminary) for four instances: the minimum number of
  matches of unbalanced10 is 25. The exact SCIP proof that 24 matches are infeasible (10 233 nodes) now has a VIPR
  certificate (0.54 GB, kept outside the repository, SHA-256 recorded) accepted by VIPR's checker and by a
  separately written checker; the LP model is re-derived row by row from the instance data. VIPR's checker
  (commit 30f2951 and current master) accepts two invalid toy certificates; a hardening patch is included. With the
  closures of 24 October (37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22, lower bound 32 for balanced15), the folder now has
  a manuscript. Not yet verified by an agent.
- Result 13 verified (Claude Opus 5.5, xhigh, fresh context): CONFIRMED WITH CORRECTIONS. Hopf criterion and
  classification re-derived before reading the solver's files (class counts identical on 1.7·10⁸ networks), every
  lemma checked by hand, 16 095 return-map integrations without a cycle; corrections editorial (family conditions
  u > 0 restored in the README table, case 8's degenerate Hopf bifurcation mentioned, a symbol defined). Status
  changed to solved; manuscript compiled.
- Result 16 (reaction networks), solved (preliminary): Banaji–Boros–Hofbauer's question on bimolecular (2, 4, 2)
  networks answered: no network beyond their ten admits two asymptotically stable equilibria in the closed quadrant,
  boundary equilibria included. Exhaustive enumeration (27 195 networks, 5 716 classes), boundary lemmas, at most one
  stable positive equilibrium (exact polynomial certificates), a new index theorem for a stable origin, exact rate
  constants for the ten networks. Solver Claude Opus 5.5 (xhigh). Not yet verified.
- The word "independent" is no longer used for verification in the READMEs, manuscripts and protocol; "independent"
  in the sense of a different method became "separate".
- Started (Claude Opus 5.5, xhigh): upper bounds for the connective constants of the diamond and face-centred cubic
  lattices; Iosif's siphon–precluster duality conjecture (arXiv:2603.08767); Becker–Döring with atomization, Hopf
  branches j ≥ 4.
- Result 17 (reaction networks), solved (preliminary), negative answer: Iosif's conjectured duality between
  preclusters and maximal invariant polyhedral supports, and between preclusters and siphons (arXiv:2603.08767), is
  false in every precise formulation extracted (equality, refinement in both directions, cardinality, three readings
  of the supports): exact counterexamples satisfying all hypotheses, the smallest being X₁ ⇄ X₂ and the 3-cycle.
  Exact criteria for each relation are proved; the one-step and iterated definitions of preclusters coincide.
  Solver Claude Opus 5.5 (xhigh). Not yet verified.
- β-Kirchhoff Šoltés graphs (Gottwald–Majstorović Ergotić–Došlić 2026, Problem 4), partial result (preliminary):
  the source's Theorem 15 already settles every β < 2/3; every rational β < 2/3 is moreover the exact proportion of
  good vertices of infinitely many graphs; a barrier theorem shows that splice constructions with branches of at most
  nine vertices cannot exceed 2/3; C₅ is the only connected graph with at most 10 vertices in which two thirds of the
  vertices are good. Open for β ≥ 2/3. Solver Claude Opus 5.5 (xhigh). Not yet verified.
- Diamond and face-centred cubic lattices (polymer physics), record (preliminary): new rigorous upper bounds for the
  self-avoiding-walk connective constants, μ(diamond) ≤ 2.8958880 (previous best 2.90940, Alm 1993) and
  μ(fcc) ≤ 10.1630815 (previous best 10.240586, derived here from the 2017 exact counts; Alm's published bound
  10.27777). Forbidden-set automata with loop memory K = 26 (diamond, 348 320 130 states) and K = 11 (fcc), as in
  result 9; every certificate passes a checker written separately from the generator, the small ones also a
  brute-force checker. Solver Claude Opus 5.5 (xhigh); checks re-run by GPT-6-Astra (xhigh). Not yet verified.
- Becker–Döring with atomization, Hopf branches j ≥ 4 (follow-up of result 8), partial result (preliminary):
  for every branch j = 4, …, 13 and every integer M, computer-assisted proof (ball arithmetic) that the branch has
  exactly one Hopf point from its onset M_j = 142, 210, 289, 379, 480, 591, 713, 844, 985, 1135 on (onsets counted
  in the limit K → ∞ for M ≤ 1150), that the Hopf hypotheses hold, that the equilibrium has 2j − 2 unstable
  eigenvalues there, and that the first Lyapunov coefficient is negative except at the onsets M = 142 (j = 4) and
  M = 1135 (j = 13), where the bifurcation is subcritical; the bifurcating orbits are unstable. For every j,
  ℓ₁/M tends to result 8's closed-form limit, which is negative (supercritical for large M, threshold not
  explicit). Open: j ≥ 14 at moderate M. Solver Claude Opus 5.5 (xhigh), editing GPT-6-Astra (xhigh), branches
  9–13 run by Claude Opus 5.5 (xhigh). Not yet verified.
- Results 15, 16 and 17 verified (fresh context): 15 and 16 CONFIRMED, 17 CONFIRMED WITH CORRECTIONS; statuses set
  to solved. β-Kirchhoff Šoltés families: the 2/3 barrier is stated for splices with branches of at most 9 vertices.
- Root README: the "Under investigation" section is merged into Results as results 18–25, each with its current
  status.

## 26 October 2026

- Result 8 completed for branches 1–3 (Codex, GPT-6, effort not recorded): the all-M argument is
  audited and integrated into the manuscript. Stable supercritical cycles on branch 1 for every M ≥ 16,
  the subcritical exception M = 15, and unstable supercritical cycles on branches 2–3 are certified
  for every integer on their curves. Exact onsets 15, 44, 87, all branch junctions, and the closed-form
  limit of the discrete Lyapunov coefficient are proved. The winding counter now uses Arb arguments
  and a unique enclosed integer; contour geometry and coverage use rigorous checks. The stored upper
  endpoint is exactly 0.078811099999999995 (the former 0.0788111 was a rounded label), still above
  1/sqrt(161). All 1,985 infinite-range coefficient boxes, the full spectral covers and all 337 small-M
  points are rechecked. The original finite-range verification is retained; the completion audit is
  all-M/verification.md. Result 25 and global-in-parameter questions retain their separate scope.

## 27 October 2026

- Closed the Kaibel–Petlyuk candidate as **already solved**, now result 26: Halvorsen and Skogestad,
  AIChE Annual Meeting 2006, paper 216d, already give the strict vapour penalty and an algebraic
  calculation in the ideal infinite-stage Underwood model. Added an attributed exposition, explicit
  strict-inequality proof and exact rational certificates for seven examples. The printed `min` in
  the source's equation (7) is corrected to `max`, consistently with its text and Table 1. The source
  example is reproduced. Audit, exposition, certificates and editing: GPT-6-Astra (effort High). No
  new research priority or separate fresh-context verification is claimed.
- Credits corrected: the completion of result 8 (26 October) is by GPT-6-Astra, effort xhigh (previously
  "Codex, GPT-6, effort not recorded"); remaining "via Codex CLI" mentions dropped.
- Result 27, Kaihnsa–Nguyen–Shiu conjecture on ACR with multistability (Eur. J. Appl. Math. 2024, §6), **solved**:
  true for dimension, reactant complexes and complexes, false for species. Main theorem: a two-dimensional
  bimolecular mass-action system with ACR has at most one exponentially stable positive steady state per
  compatibility class (normal form on the ACR hyperplane, connectedness of the positive steady-state set from
  bimolecular products, global ACR). Consequences: dimension ≥ 3, ≥ 4 reactant complexes, ≥ 6 complexes
  (kinetic-subspace lemma and deficiency-one theorem). The 3-species, 6-complex, 6-reaction network N1 has ACR for
  every rate constant vector and two stable steady states (exact certificate), so dimension 3, species 3 and
  complexes 6 are attained; the minimal number of reactant complexes (4, 5 or 6) stays open. Also noted: KNS state
  the deficiency-one theorem for weakly reversible networks only but use it for general ones; the general version
  and the kinetic-subspace lemma close that gap. Solver and two sub-studies: Claude Opus 5.5 (xhigh). Verified the
  same day in a fresh context (Claude Opus 5.5, xhigh): CONFIRMED WITH CORRECTIONS (scope of the 5-reaction
  enumeration restricted to full-dimensional networks, two wording fixes, one missing sentence). The candidate entry
  is removed from the root README.
