# Independent verification — Halvorsen–Skogestad conjecture (generalized Petlyuk minimum vapor)

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, 2026-10-20. Not opened: any `PROGRESS.md`,
`FINAL-REPORT.md`, other `BRIEF-*` files, `tmp/`. `computation/REPORT.md` (the computation agent's final report,
copied by the orchestrator as a deliverable) was read to check the computer-verification section.

PDF requests: none (all sources needed were present in `certificates/sources/`).

## Verdicts

| Claim | Verdict |
|---|---|
| 1. Strict model, every N: every two-section column network has R ≥ max_j H_j; the FTC attains it | **CONFIRMED WITH CORRECTIONS.** The theorem is correct as a statement about the model as defined (TG's nominal-split formulation plus V_T = Y_T whenever two or more components are nominally distributing, nonnegative net flows, single feed, pure saturated-liquid products, the class of Definition 1). Corrections concern how the model and class are described (degenerate points, TG equation numbering, Nallasivam's model not identical, reverse flows), and the computer-verification section. |
| 2. Permissive model, partial results (N ≤ 4; extreme peaks for every network; merge-free mixed networks under (N); the FTC itself; simple diamonds) | **CONFIRMED WITH CORRECTIONS** for every proved statement. The description of the remaining N = 5 cases as "open" is **corrected**: the permissive bound is **false** for the pair (basic configuration 161, j = 3): I found and certified (interval Krawczyk test) a permissive-feasible point with R/V_FTC ≈ 0.999455, all flows ≥ 1e-6 and all inequality slacks ≥ 1e-6. It is an artifact of the permissive formulation (excess vapor in a column with an active root), not a physical counterexample (see §6). The other 19 pairs remain undecided. |
| 3. Exclusions stated as hypotheses | Honest, but incomplete: added (a) degenerate operating points, where the strict model is more restrictive than infinite-stage physics, (b) splits in which a product keeps every component of the feed, (c) merging of streams with different labels. (H1) is a genuine restriction of the operating points (Halvorsen §4.3.3). |

Status: in PROTOCOL vocabulary, **solved** (complete proof, independently verified) for claim 1 *in the stated
class and model*; **partial result** for the permissive model, now including a negative result for N = 5. The HS
conjecture in its literal wording ("any adiabatic distillation arrangement") is not settled outside the class.

## 1. Statement and scope

**Quotations.** Checked against `HS2003_part3.txt` (text of the publisher PDF): the abstract sentence, the §6.3
sentences (page 626) and the §5.3 sentence (page 624) are verbatim. The definition of an adiabatic section (§6) is
quoted correctly. HS Part 3 also states the Petlyuk result (eq. 8) with the vapor referred to the top; H_j =
P_j^T − (1 − q) is the reboiler-referenced value, and R (total vapor generated) is the right objective for "need for
vaporization" (R − Q_total = −(1 − q) for liquid products, so minimizing R is minimizing heat input).

**Is the strict model the right formalization?** Halvorsen's thesis (§4.2.1, p. 97): roots are active when they lie
between the volatilities of the components distributed to both ends, and "If we have any active roots then V = Vmin";
§3.2.4 defines a component as distributing if it "appears in both products, or is exactly at the limit of becoming
distributing". Hence:

- At *nondegenerate* points (every nominally distributing component l..k has positive flow in both products) the
  interior roots are active, and the strict model (V_T = Y_T if l ≤ k − 1) is exactly the infinite-stage behaviour.
  The permissive model, which also allows V_T > Y_T there, is not physical (confirmed dramatically by §6).
- At *degenerate* points (a nominally distributing component absent from one product) physics is governed by the
  actual distribution: with fewer than two actually distributing components, no root is active and an infinite-stage
  column may carry excess vapor, which the strict model forbids; and TG's equalities at the nominal interior roots are
  then not physically required (e.g. ABCDE → ABCD/BCDE with d_D = 0 but B, C distributing: physics needs only
  V_T = g_T(θ_2) ≥ g_T(θ_3), TG imposes g_T(θ_2) = g_T(θ_3)). Relabelling the split does not produce a network of the
  class in general, because the delivered stream would not carry its receiver's label (merge nodes). **So Theorem
  "strict" does not cover these physical points.** The paper said the strict model is the infinite-stage behaviour
  "at any operating point"; corrected (new paragraph "Degenerate operating points" in §2.2, scope item 3, README).
  The sign-rule lemma only uses V_T ≥ g_T(θ) at a common root, valid at any infinite-stage point with nonnegative
  flows, so the sign-rule results (path lemma, extreme peaks, FTC theorem, diamond) very probably extend; the
  sign-graph proof of the strict theorem does not extend as written (it needs the interior lemma, which fails when the
  nominally distributing components do not actually distribute). I checked by hand that the degenerate points of
  family K (d_C = 0 or b_B = 0 in the feed column) still satisfy R ≥ H_2 via the cut identity with U ∋ ABC (resp. BCD).
- TG's constraints are numbered (A18)–(A23) in the arXiv v1 PDF, not (B1)–(B6) (corrected in paper and
  references.md). Nallasivam et al. 2016 (Eqs. 14, 17–18b) use V ≥ Underwood value at *every* root plus explicit
  recovery-ordering constraints, not TG's equalities; "identical to Nallasivam" was wrong (corrected).

**Competitor class.** Definition 1 is stated precisely and contains the Shah–Agrawal space (my own enumeration
reproduces 3/8, 18/152, 203/6128). Exclusions listed by the paper: several feed stages, stream splitters, recycles,
non-interval streams, reverse flows. Added: splits with k = m or l = i (a product keeping every component of the
feed — physically possible at low reflux, excluded by i < l, k < m); streams of different labels merged at one feed
stage. Remark added: pure products withdrawn as saturated vapor would only raise the bound (a heavy vapor product
contributes θ_j z_p/(θ_j − α_p) > α_p z_p/(θ_j − α_p) in the node balance), so (H2)'s phase restriction is harmless.

**Reverse flows (H1).** Halvorsen §4.3.3 shows reverse net flows in a directly coupled prefractionator operated
away from the preferred split (B entering with the top liquid and flowing down). So (H1) restricts the operating
points of thermally coupled networks, not only a modelling convenience; the paper's phrase "who meets reverse flows
only away from optimal operation" was replaced by an accurate statement. Whether reverse flows can lower R remains
open (stated).

**Overclaiming.** Title and abstract qualify the claim by the class and the Underwood model; I added "and model",
"in this model (without reverse flows)" for the FTC formula, the degenerate-point caveat, and the permissive
counterexample. The README status now reads "solved in a precisely defined class and model" and says explicitly
that the literal conjecture is settled only within it.

## 2. Re-derivations (by hand, from the definitions)

- **Underwood flux and cut identity.** Re-derived the node balance Σ_end Λ_s − Σ_start Λ_s = Q_c − Q_r − Λ^feed(x)
  − α_p z'_p/(α_p − x) from conservation of vapor and net component flows at every node type (feed stage with heat
  duties, partial condenser / reboiler on transfer streams, thermal couplings, total condenser / reboiler and liquid
  side draws of products); summing over U^c gives (eq:cut) with H_j from the heavy products and Λ^feed(θ_j) = 0.
  Correct, needs no sign assumption. Checked also in exact rationals (§4).
- **Lemma mono, sign rule, interior lemma.** Correct, including degenerate cases: a root at an interval end forces
  f = d = b = 0 for that component (no pole); f ≡ 0 forces U = 0 and A = V_T, B = V_B ≥ 0; the three cases
  (l−1 ≤ j ≤ k, k < j, l > j+1) are exhaustive and disjoint.
- **Strict theorem (sign graph).** Checked the four ways a mixed node can be entered (down its own top section, up
  its own bottom section, up a producer's top, down a producer's bottom), that light and heavy nodes are never
  adjacent, that producers of mixed nodes are mixed, and that a heavy node can only be entered down a mixed column's
  bottom section; the invariant (i) B ≥ 0, (ii) A ≤ 0 if dist is mixed is path-independent. Correct for every network
  of Definition 1 with heat duties at any node. Agent G's flux proof (G2) is also correct and does not actually need
  (N) (it uses only the sign rule, the interior lemma and superadditivity of min(0, ·)).
- **Corollary (tie condition).** Correct; wording sharpened to "with positive duty".
- **Path lemma, extreme peaks.** Correct for general networks (all mixed columns lie on the root's chain, merges are
  impossible when every mixed column has ≤ 1 mixed child).
- **FTC theorem (grid).** Correct (descending induction A ≥ X_m, telescoping Y_m = X_m − A).
- **N = 4:** the case split (path / diamond / K) is exhaustive (ABC without BC forces AB, BCD without BC forces CD,
  so K is unique with 16 variants); diamond lemma, recovery-ordering lemma (diagonal + monotonicity in d_B) and
  family-K lemma (four sign cases, X/(X+Y) + Z/(Z+W) ≤ 1, Rat decreasing) are correct.
- **Lemma G1 (zero counting).** Correct: sign changes only at r ∈ [l−1, k], each such interval carries ≥ 2 zeros
  (interior roots: a zero with equal end signs; boundary roots: nonnegative value between two negative ends),
  total ≥ n + 2 > deg P.
- **Lemmas G4/G5, Theorem G3.** Checked all four types and the three sub-cases of type (d); correct under (N)
  (needed for the zeros ζ_y and for G1). **G8 (simple diamonds):** verified the light-form and heavy-form block
  identities and the two cases; correct.
- **Attainment.** Correct: partial fractions, sign counting (every d_p has the sign of V_T, every b_p that of V_B,
  fixed by d_i = f_i > 0, b_m = f_m > 0), 0 < L_T < V_T, L_B > V_B, carry-over of θ_i..θ_{m−1} to every child,
  binary columns with excess e_r, and the side-draw imbalances δ_p = 0 from the cut identity for r = N−1, …, 1.
  The definition of V_T ("d_i = f_i") was ambiguous; clarified.

## 3. Re-runs of the solvers' certificate scripts

| Script | Result | Runtime |
|---|---|---|
| `computation/scripts/run_certificates.sh` (export_configs, verify_identity_exact, verify_sign_paths, check_lemmas 4 20 7 0) | ALL PASS (8/152/6128; exact cut identity on all configurations N ≤ 5; 0 strict paths N = 3..6; permissive paths 0/1/70/5795; lemma checks 69 points) | 31 s |
| `theory/scripts/run_theory_checks.sh` (check_algebra, check_ftc_achievability, signcheck certify4) | ALL PASS (1200 FTC instances N = 3..8; 16 = family K) | 91 s |
| `general/scripts/recovery_order_lp.py 300 1`, `classify_mixed.py 4`, `classify_mixed.py 5`, `gamma_check.py 4 20 3 1` | PASS / counts reproduced / PASS | < 10 s each |
| `theory/scripts/opencount.py` | 864 (variant, j) pairs, 720 configurations, 16 basic, 20 (basic, j) pairs; identical to Agent C's 20 path pairs (checked set equality); 58 sign-hard basic pairs = 58 path pairs; union 2160 configurations | s |
| `theory/scripts/signcheck.py 5` (full listing) | 1312 sign-hard pairs for j = 2 and j = 3; the committed `theory/data/signcheck_N5.txt` (1456) was stale and has been replaced | 6 min |

## 4. Independent implementation (`verification/`, no code shared with the solvers)

- `indep_checks.py` (ALL PASS, ~7 s): (a) max_j H_j at TG's eight five-component cases equals their certified
  optima (402.703, 272.485, 260.043, 896.351, 695.581, 929.098, 902.712, 542.047; TG print 4 flows for Case D, the 5th
  is 5); all eight have an extreme highest peak (j = 1 or 4), so they are covered by the extreme-peak theorem in both
  models; (b) Halvorsen Table 6.1 in closed form: Petlyuk 1.3657, direct split 2.0718, indirect split 2.0324 (with
  TG's convention of a saturated-vapor distillate after a partial condenser; a total condenser would give 2.3657);
  (c) exact rational check of the node-balance cut identity on 400 random networks of Definition 1 (several columns
  per label, merges, heaters/coolers at every node, random U, rational x); (d) own construction of the attainment
  point, 60 digits, 200 instances N = 3..7, q ∈ [−1, 2], every strict constraint, side draws, R = max H (residual 5e-56).
- `indep_nlp.py`: own enumeration and NLP (pyomo + ipopt, multistart, a-posteriori feasibility check at 1e-7, and an
  independent evaluation of the light-cut identity at every returned point). N = 3: all 8 configurations reproduce
  Table 6.1 (1.8819 SS/SR, 2.0324, 2.0718, 1.3657, 1.5556). N = 4: 13 runs (152 configurations each) at 8 parameter
  points (q = −0.3, 0.3, 0.4, 0.5, 1, 1.3; four with the middle peak highest), strict and/or permissive: no configuration
  below V_FTC beyond 1e-7 (minimum ratio 1 − 4.5e-8); at middle-peak points 16 configurations tie (including the
  all-coupled family-K variants), as Agent C reports; the paper's example values 1.588798 (FTC and K-TC), 1.647249 (K
  with HX on ABC, BCD, AB) and 2.481686 (indirect) are reproduced; `tie_check.py`: all 120 observed ties satisfy the
  tie corollary (PASS).
- `open_n5_probe.py`: permissive model, the 20 open N = 5 pairs, 220 random parameter points (H_j highest, other
  peaks within 10 %), all-coupled and all-exchanger variants. One point below V_FTC: basic 161, j = 3, all TC, q = 0
  (ratio 0.99944); all other values ≥ 1 − 1.5e-7.

## 5. Computer-verification section vs. the actual runs (computation/REPORT.md)

Corrected in the paper: random multistart searches were run only for N = 4, q = 1 (239 points × 152 configurations)
and, for N = 5, only at TG's eight points (× 6128); adversarial runs "never below beyond 1e-6" (not 1e-9); SCIP was
used at two points only (strict, all 152 at α = (8,4,2,1), z = 1/4, q = 1; permissive, first 100 configurations at
the example point), and for tying configurations the SCIP dual reached V_FTC(1 − 8e-6) at the time limit; random
searches with q ≠ 1 and a random permissive N = 5 search were not run. Counts 8/152/6128, 70/48, 5795/2300, 1312,
2160, 864/720/20 and 12–16 ties confirmed. The N = 5 statement "settled by the sign rule" rests on a floating-point
LP (integer data, HiGHS); now said so. `general/general-N.md` §7–8 carried stale counts (1456, 28, 1152): note added.

## 6. The permissive counterexample (N = 5)

Configuration (basic 161, all thermal couplings; B, C, D liquid side draws): ABCDE→ABCD/BCDE, ABCD→ABC/BCD,
BCDE→BCD/CDE, ABC→AB/BC, BCD→BC/D, CDE→C/DE, AB→A/B, BC→B/C, DE→D/E; q = 0;
α = (10675204/742583, 5824255/854876, 529647/171259, 1673613/789461, 1) ≈ (14.376, 6.813, 3.093, 2.120, 1); rational
feed roots θ ≈ (12.325, 5.554, 2.530, 1.950) and the exactly rational z ≈ (0.1904, 0.1883, 0.1709, 0.0389, 0.4116)
having them as roots (so H and V_FTC are exact rationals): H ≈ (0.33442, 0.32912, 0.46994, 0.43302),
V_FTC = H_3 ≈ 0.4699411. Permissive point: R ≈ 0.4696851 (ratio 0.999455).

Exact check (as requested by the orchestrator): the raw ipopt point violated some inequalities by 1e-8 and was
degenerate (B ≈ 0 in BCDE); `refine_candidate.py` re-solved it with hard bounds (no bound relaxation), flows ≥ 1e-6
and every inequality slack ≥ 1e-6, at the rational data; `krawczyk_certificate.py` then (1) polishes the 57
equalities in 60 digits (76 unknowns), (2) fixes 19 variables at exactly representable values and Newton-solves the
57 × 57 system (pivot block condition number ≈ 420), (3) proves with the Krawczyk operator in interval arithmetic
(interval Jacobian by forward AD over intervals) that an exact zero exists within 1e-30, and (4) bounds every TG
inequality (all roots strictly inside (α_{r+1}, α_r), feed equations, interior equalities, boundary inequalities,
0 ≤ Y ≤ V_T, Y − U ≥ 0, flows, vapors) and R over the box: min slack ≥ 1e-6, R < V_FTC. **PASS.** The FTC's own
minimum at these data is V_FTC (Theorem ftc), so the FTC is not optimal in TG's literal model here. In the strict
model the same configuration gives R/V_FTC = 1 to 1e-7 (`strict_at_counterexample.py`), as Theorem strict requires.

Mechanism: the type-(d) column BCDE (C and D distributing with positive flows on both sides, root θ_{BCDE,3} active)
runs 0.006 above its threshold (V_T ≈ 0.03015, Y_T ≈ 0.02420) and its distillate feeds the merge node BCD, whose
other producer ABCD is of type (b): not a simple diamond. This is exactly the gap identified by Agent G.

**Physical meaning.** I do not know an adiabatic arrangement with pure products realizing this point, and the natural
candidates fail:
- *Infinite-stage column:* impossible (active root ⇒ V = V_min, Halvorsen §4.2.1).
- *Finite stages:* excess vapor at a fixed product distribution needs fewer stages, but a section with finitely many
  stages cannot remove B from the residue or E from the distillate completely; as the impurities tend to zero both
  sections become infinite, the feed-end compositions tend to the two section pinches, and pinch matching forces
  V → V_min. So the point is at best a limit of finite-stage operations with impure intermediate streams (E in the
  stream sent to BCD, B in the stream sent to CDE), which the conjecture does not compare and which would need extra
  separation.
- *Vapor bypass / stream splitting and remixing (no heat exchange):* routing the excess vapor e from the CDE node to
  the BCD node, or adding a V = L section between them, moves components (vapor from the CDE node carries E into the
  BCD node, liquid from the BCD node carries B down; a V = L section transports C and D in opposite directions),
  whereas the permissive excess leaves every net flow unchanged — a different point, outside Definition 1, and with
  E entering the BCD column it no longer separates pure products downstream.
Conclusion: a flaw of the permissive (TG-literal) formulation, not a counterexample to a broader physical reading of
the HS conjecture; it does show that global optimization with TG's model can return configurations below the FTC
value (TG's eight published cases are unaffected: q = 1 and extreme peaks). This is an argument, not a proof that no
exotic adiabatic arrangement could emulate the point.

## 7. Novelty

Re-checked (details appended to `references.md`): OpenAlex (119 citing works; 16 from 2024–2026 with abstracts
scanned), web search, Zenodo API (4 queries), GitHub (`ipitchford`: 171 repositories, none on distillation;
`mmislan`: 24, none; `gh search repos`): no proof or refutation of the conjecture, no competing AI candidate.

## 8. Corrections applied

`paper.tex` (recompiled, 19 pages, no undefined references, no overfull lines):
1. `\PaperCredits`: "Independent verification: Claude Opus 5.5 (effort xhigh)"; subtitle "independently verified".
2. Abstract: strict model described via nominally distributing components; "in this class and model"; FTC formula
   "in this model (without reverse flows)"; the N = 5 permissive counterexample and its status; degenerate points.
3. §2.1: added excluded constructions (product keeping all feed components; merging streams of different labels).
4. §2.2: TG numbering (A18)–(A23); Nallasivam's formulation is related, not identical; strict model = infinite-stage
   behaviour only when the nominally distributing components actually distribute; new paragraph "Degenerate operating
   points"; (H1) paragraph rewritten (reverse flows do occur away from the preferred split).
5. Corollary tie: "with positive duty".
6. Remark after Theorem 1: part 2 does not extend to all configurations (Section sec:cex).
7. New subsection "A five-component counterexample in the permissive model" with data, certificate, mechanism and
   physical meaning.
8. Attainment proof: V_T fixed by d_i = f_i (clarified), d_m := 0.
9. Computer section: Agent C's searches and SCIP cross-checks described as actually run; TG cases have extreme
   peaks; signcheck LP is floating point; new bullet for the independent verification.
10. Scope: item 1 (more exclusions), item 3 (degenerate points), item 4 (pair (161, 3) fails, 19 undecided, LP
    classification), item 5 (vapor products).

`README.md`: attribution and status (solved in class and model, verified with corrections; permissive bound false for
N = 5; degenerate points and literal conjecture), strict/permissive description, permissive results and
counterexample, floating-point exceptions, verification reproduction commands, files table.
`references.md`: TG numbering, Nallasivam note, verifier's novelty re-check.
`general/general-N.md`: verifier's note (stale N = 5 counts; the open case is false for (161, 3); proofs checked).
`theory/data/signcheck_N5.txt`: replaced by the output of the current `signcheck.py 5` (1312 per j).
No solver script was modified. Note: my run of `run_certificates.sh` rewrote `computation/certificates/configurations_N*.json`
(deterministic export; same content).

## 9. Remaining caveats

- Theorem strict is about TG's nominal-split model; degenerate physical points are not covered (§1).
- (H1) no reverse flows; single feed; class of Definition 1 (no multi-feed columns, splitters, recycles).
- The N = 5 permissive classification "settled by the sign rule" is a floating-point LP; the 19 remaining permissive
  pairs are undecided (no other point below V_FTC found).
- Numerical searches (theirs and mine) are local; they are evidence, not proofs.
- The root `README.md` of the repository has no row for this problem yet (left to the orchestrator).

## Files written by the verifier

`verification.md`; `verification/indep_checks.py`, `indep_nlp.py`, `tie_check.py`, `open_n5_probe.py`,
`open_pairs_N5.json`, `analyze_candidate.py`, `certify_candidate.py` (exact fixed-root LP attempt; fails because with
the roots fixed at rounded rationals three feed rows of BCDE become dependent with inconsistent right-hand sides —
superseded by the Krawczyk certificate),
`refine_candidate.py`, `krawczyk_certificate.py`, `strict_at_counterexample.py`, `permissive_counterexample_N5.json`,
and `verification/logs/` (all run outputs).
