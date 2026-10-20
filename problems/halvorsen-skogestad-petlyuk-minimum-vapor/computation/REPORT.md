# Computation report (Agent C)

Computation: Claude Opus 5.5 (effort xhigh). This report was written by the computation agent at the end of its log and
copied here verbatim by the orchestrator (the agent's tool policy did not allow it to create this file). Paths are
relative to the problem folder.

alpha_1 > ... > alpha_N > 0, F = 1, z_p > 0, any real q, CMO, constant pressure, infinitely many stages, pure
saturated-liquid products; R = total reboiler vapor; V_P = max_j P_j - (1-q) = max_j H_j. Class = Shah–Agrawal basic
regular-column configurations with all TC/HX choices, TG conventions (arXiv:2010.12113 §4 = Nallasivam 2016):
partial condenser -> saturated-vapor transfer, reboiler -> saturated-liquid transfer, side draws without
exchanger, nonnegative net flows, one feed stage per column, Underwood at each column's own roots. Interior roots:
STRICT (V_T = Underwood value; Halvorsen thesis §4.2.1 p. 97, §3.2.4) or PERMISSIVE (TG literal, V_T >= Upsilon).

## 2. Enumeration — PASS (export_configs.py): 3+5=8, 18+134=152, 203+5925=6128 = Shah–Agrawal 2010 Table 1.
Lists: certificates/configurations_N{3,4,5}.json.

## 3. Validation (numerical)
- FTC closed form = TG's 8 certified optima: 402.703, 272.485, 260.043, 896.351, 695.581, 929.098, 902.712,
  542.047 (published 402.7 ... 542.0).
- tg_cases_N5.py: min over all 6128 configurations, each of the 8 cases (strict, multistart ipopt) = FTC value to
  1e-10 relative; FTC attains it; 62–205 configurations tie per case; 4–14% of configuration solves failed.
- Ternary alpha=[4,2,1], z=1/3, q=1 = Halvorsen Table 6.1 (2.0324, 2.0718, 1.8819, 1.5556, 1.3657).
  Fidkowski–Królikowski: side stripper = side rectifier at q=1 (5e-11 on 20 random points); orderings
  Petlyuk <= SS/SR <= indirect/direct on 25 random points with q in [-0.3, 1.3].
- SCIP global vs ipopt (N=4): strict, alpha=[8,4,2,1], z=1/4, q=1, all 152 configurations: max rel difference
  2.2e-9, SCIP dual >= V_P for all non-tying configurations (tying FTC-type ones: dual reached V_P(1-8e-6) at the
  300 s limit); permissive, P2-dominant point alpha=[47.5774,25.5263,15.1071,1], z=[.1358,.5169,.0923,.255],
  first 100 configurations: max rel difference 3.7e-9, min SCIP dual/V_P = 0.99999983 (ties).

## 4. Main result: STRICT model, every N (proof by Agent C; checked by T; independent flux proof by G)
Cut identity R - H_j = FR(U) + FC(U^c) + sum_{U^c->U} Psi_s(theta_j) - sum_{U->U^c} Psi_s(theta_j) for every
node partition separating light from heavy products; sign lemmas S3, S4, S5 (strict), Root; sign-graph invariant
=> R >= H_j for all j at every feasible point of every two-section column network, every alpha, z > 0, q real.
Exact/combinatorial certificates (scripts/run_certificates.sh, ~30 s, ALL PASS):
 export_configs.py; verify_identity_exact.py (identity in Fraction arithmetic, all 8+152+6128 configurations, all
 j, 5 partitions each, arbitrary junction duties incl. extra heaters/coolers at every node);
 verify_sign_paths.py (0 sign-consistent light->heavy paths, all basic configurations N=3..6, strict model);
 check_lemmas.py 4 20 7 0 (numerical: L1 sign rule, recovery ordering, S5, sign graph — PASS).
Larger numerical lemma runs (all 0 violations): N=4 strict 204 points/612 tests; N=4 permissive 198/594;
N=5 strict 123/492 and 866/3464 (certificates/runs/check_lemmas_N5_strict*.log).
Tie condition (exact): if j* is the strict top peak, R - V_P >= FR(light-only nodes) + FC(heavy-only nodes).

## 5. PERMISSIVE model
Sign-consistent paths (S3, S4, Root only): N=3: 0; N=4: 1 (family K, j=2; settled by T Lemma 5 / G3);
N=5: 70 paths, 58 (basic, j) pairs, of which 38 merge-free (G3) and 20 with a mixed merge node = the open set
(864 (variant, j) pairs); N=6: 5795 paths in 2300 configurations (certificates/potential_paths.json).
Numerics (never below 1): K: bad path feasible for ABC:TC variants, min R/(P_2-(1-q)) in the bad-path region
1.00000005; N=5 all-TC variants of the 48 families: 56 triples, bad paths feasible, min ratio 1.0000000193;
merge families (badpath_parallel, stopped at 233/1248 tasks, 181 (variant, j) pairs): 56 feasible bad paths,
38 in-region minimisations, min ratio 1.0000000258.

## 6. Counterexample search (numerical evidence only)
N=4 strict q=1, "mixed" distribution (volatility ratios 1.02–20, half the points with trace components):
239 points x 152 configurations, min ratio 1 - 1.4e-9 (noise; one sub-1 value not reproducible -> a-posteriori
verification added to the solver). N=5: the 8 TG points x 6128 configurations. Adversarial joint (alpha, z, q,
operating point) minimisation of V_c/(P_j-(1-q)): ratio -> 1, never below beyond 1e-6. No counterexample in any
model. Not run (stopped): q in [0,1]/subcooled/superheated random searches, permissive random N=5 search.

## 7. Tight cases
N=4, q=1 (certificates/ties_by_peak_N4_q1.json): j*=1: 12 configurations tie, j*=2: 16 (incl. K with ABC:TC,
BCD:TC), j*=3: 12; FTC at every point; 2717/2720 observed ties satisfy the exact tie condition (3 near-ties from
trace components).

## 8. Certified vs numerical — see 4 (certified: analytic proof + exact checks) vs 3, 5-numerics, 6, 7 (numerical).

## Numbers for the paper (T): counts 8/152/6128; TG cases reproduced to 1e-10 by min over all 6128
configurations (strict); SCIP cross-check N=4 at 2 points (strict all 152: max diff 2.2e-9; permissive 100);
random search 239 x 152 (N=4, q=1), no ratio below 1 - 1.4e-9; verify_identity_exact / verify_sign_paths /
check_lemmas all PASS via scripts/run_certificates.sh; permissive open set N=5 = 20 (basic, j) pairs / 864
(variant, j) pairs, no counterexample found there (min ratio 1.0000000258 on 181 pairs tested).
