# Minimum number of matches on the Letsios–Kouyialis–Misener benchmarks: exact bounds

Credits: Solver (Secs. 1–3, first exact runs; 24 Oct 2026): Claude Opus 5.5 (effort xhigh). Solver (Sec. 4 items
5–7, the VIPR closure of unbalanced10 and Sec. 7; 25 Oct 2026): Claude Opus 5.5 (effort xhigh). Verification: GPT-6-Astra (effort xhigh); see `verification.md`.

Status of each claim: see the table in `README.md`. This file gives the definitions, the proofs of every valid
inequality used, the certification pipeline, and every assumption or non-rigorous step.

## 1. Problem and data conventions

An instance (file `data/mip_instances/<set>/<name>.dat` of github.com/dimletsios/min_matches_heuristics) has hot
streams i ∈ H (|H| = n, hot utilities included with their fixed loads), cold streams j ∈ C (|C| = m, cold
utilities included), and temperature intervals t = 0, …, k−1 (t = 0 hottest). QH_it ≥ 0 is the heat supplied by
hot stream i in interval t, QC_jt ≥ 0 the heat demanded by cold stream j in interval t. All numbers are read as
exact decimal rationals. Notation: H_i[a..b] = Σ_{a≤t≤b} QH_it, H_i(≤t) = H_i[0..t], D_j[a..b], D_j(≤t) likewise,
H_i = H_i(≤k−1), D_j = D_j(≤k−1), N = n + m. The heat cascade is R_t = Σ_{s<t}(Σ_i QH_is − Σ_j QC_js), t = 0..k.

*Checked exactly for all 54 data files* (`scripts/check_solutions.py` re-checks it for the instances below): the
instance is balanced (R_k = 0), R_t ≥ 0 for all t, the R values printed in the file equal the exact cascade, and
every stream has positive total heat.

**Feasible network.** Numbers q_ijt ≥ 0 (heat of hot i received by cold j in interval t) such that
- Σ_i q_ijt = QC_jt for all j, t (cold balances), and
- r_it := H_i(≤t) − Σ_{s≤t} Σ_j q_ijs ≥ 0 for all i, t, with r_{i,k−1} = 0
  (heat of hot i can only be used in its own or a colder interval; the ΔT_min shift is built into the intervals).

Its *match set* is E(q) = {(i,j) : Σ_t q_ijt > 0}. The **minimum number of matches** is min |E(q)| over feasible
networks. This is exactly the optimum of the transshipment MILP of Letsios et al. (2018), Sec. 2 (variables
y_ij ∈ {0,1}, q_ijt, r_it; Σ_t q_ijt ≤ U_ij y_ij with any valid U_ij), and of their transportation MILP.

## 2. Combinatorial lower bound: the block lemma

**Definition.** A nonempty set K ⊆ H ∪ C is a *block* if Σ_{i∈K∩H} H_i = Σ_{j∈K∩C} D_j and, for every t,
Σ_{i∈K∩H} H_i(≤t) ≥ Σ_{j∈K∩C} D_j(≤t). Let c_max be the largest number of parts in a partition of H ∪ C into
blocks (c_max ≥ 1 since H ∪ C itself is a block).

**Lemma 1.** For every feasible network q, |E(q)| ≥ N − c_max.

*Proof.* Let G be the bipartite graph on H ∪ C with edge set E(q). Every stream has positive total heat, so it is
incident to at least one edge; hence every connected component K of G has |K| ≥ 2 vertices and at least |K| − 1
edges. All heat of a hot stream i ∈ K goes to its neighbours, which lie in K, and all demand of a cold stream
j ∈ K comes from its neighbours, which lie in K; therefore Σ_{i∈K} H_i = Σ_{i∈K, j∈K, t} q_ijt = Σ_{j∈K} D_j. For
every t, Σ_{j∈K} D_j(≤t) = Σ_{j∈K} Σ_{s≤t} Σ_{i∈K} q_ijs = Σ_{i∈K} (H_i(≤t) − r_it) ≤ Σ_{i∈K} H_i(≤t). So K is a
block, the components form a partition of H ∪ C into c(G) ≤ c_max blocks, and |E(q)| ≥ Σ_K (|K| − 1) = N − c(G)
≥ N − c_max. ∎

**Lemma 2 (structure at the bound).** If |E(q)| = N − 2 and c_max = 2, then G is a forest with exactly two
components, whose vertex sets {A, B} form a 2-block partition, and the restriction of q to A (resp. B) is a feasible
network of the sub-instance A (resp. B) with |A| − 1 (resp. |B| − 1) matches (a spanning tree).

*Proof.* From the proof of Lemma 1, N − 2 = |E| ≥ N − c(G) forces c(G) = 2 and |E ∩ K| = |K| − 1 per component. ∎

Lemma 2 is used only to *search* for solutions (the upper bounds are separately certified by explicit networks).

**Computation of c_max (`scripts/block_bound.py`).** With L = lcm of the decimal denominators, every stream has an
integer weight w_v = L·H_i (hot) or −L·D_j (cold) and an integer cumulative profile. All subsets with zero weight
are enumerated exactly by meet-in-the-middle (all 2^{⌊N/2⌋} and 2^{⌈N/2⌉} subset sums of the two halves, exact
hash join on the integer sums). A zero-sum subset is a block iff its cumulative profile (sum of int64 profiles of
the two halves; all partial sums are bounded by Σ|w_v| < 2^62, asserted) is ≥ 0 everywhere. Then:
- a 2-block partition exists iff some block A with stream 0 ∈ A ≠ H ∪ C has a complement that is a block;
- a 3-block partition {A, B, C} exists iff there are blocks A ∋ 0, X = complement(A) a block (a union of blocks is a
  block: balance and the cascade inequalities are additive), and a block B ⊊ X containing the lowest-index stream
  of X with X \ B a block. The search over all such (A, B) is exhaustive.
Partitions into ≥ 4 blocks give 3-block partitions by merging, so "no 3-block partition" means c_max ≤ 2.

Results (exhaustive, exact):

| instance | N | zero-sum subsets | blocks | 2-block partitions | 3-block partition | c_max | bound N − c_max |
|---|---|---|---|---|---|---|---|
| 37sp-yfyv | 38 | 4 833 | 2 431 | 14 | none | 2 | **36** |
| 20sp1 | 21 | 15 | 12 | 4 | none | 2 | **19** |
| 23sp1 | 24 | 17 | 12 | 3 | none | 2 | **22** |
| 22sp1 | 24 | 1 | 1 | 0 | none | 1 | 23 |
| balanced12, b12_random0/1/2 | 27 | ≈ 6·10^4 each | ≈ 1.8·10^4 | 0 | none | 1 | 26 |
| balanced15 | 33 | 1 817 217 | 518 054 | 0 | none | 1 | **32** |
| balanced15_random0 / 1 / 2 | 33 | 1.75 / 1.78 / 3.42 ·10^6 | 0.63 / 0.65 / 1.00 ·10^6 | 0 | none | 1 | **32** |
| unbalanced15 | 33 | 470 203 | 94 888 | 0 | none | 1 | 32 (weaker than the MILP bound) |

For the unbalanced (17, 20) instances N − 1 is far below the known bounds; for the 162-stream large-scale instances
the enumeration is out of reach.

## 3. The strengthened transshipment MILP (all rows proved valid)

Let rc_i(t) = min(R_{t+1}, H_i(≤t)) for 0 ≤ t ≤ k − 2 and rc_i(−1) = rc_i(k−1) = 0.

**Lemma 3.** In every feasible network, 0 ≤ r_it ≤ rc_i(t).
*Proof.* Summing the definitions of r_it over i and using the cold balances,
Σ_i r_it = Σ_i H_i(≤t) − Σ_j D_j(≤t) = R_{t+1}. All r_i't ≥ 0, so r_it ≤ R_{t+1}; also r_it ≤ H_i(≤t). ∎

Given a feasible network q with match set E, put y_ij = 1 iff (i,j) ∈ E. **Proposition 4.** (q, r, y) satisfies
every row below; hence if the MILP with the extra row Σ y ≤ K is infeasible, every feasible network has ≥ K + 1
matches.

- *Balances:* hot_i_t: Σ_j q_ijt + r_it − r_{i,t−1} = QH_it; cold_j_t: Σ_i q_ijt = QC_jt (definitions).
- *Bounds:* 0 ≤ r_it ≤ rc_i(t) (Lemma 3); 0 ≤ q_ijt ≤ u_ijt := min(QC_jt, QH_it + rc_i(t−1)), because the heat
  that hot i delivers in interval t is at most what is available there, r_{i,t−1} + QH_it.
- *dl_i_j_t:* q_ijt ≤ u_ijt y_ij (same bound; q_ijt = 0 if (i,j) ∉ E).
- *link_i_j:* Σ_t q_ijt ≤ U_ij y_ij, where U_ij is the minimum, over all labellings λ: {0..k−1} → {in, out}, of
  cap(λ) = Σ_{λ(u)=in} QC_ju + Σ_{λ(u)=out} QH_iu + Σ_{λ(u)=in, λ(u+1)=out} rc_i(u).
  *Proof that Σ_t q_ijt ≤ cap(λ) for every λ:* intervals labelled `in` contribute q_iju ≤ QC_ju. Split the `out`
  intervals into maximal runs [a..b]; the heat i delivers (to anybody) in [a..b] is r_{i,a−1} + H_i[a..b] − r_ib ≤
  r_{i,a−1} + H_i[a..b], and r_{i,a−1} ≤ rc_i(a−1), which is paid by the in→out transition at a − 1 (or is 0 if
  a = 0). U_ij is computed by a two-state dynamic programme; `scripts/check_lp_model.py` reconstructs the optimal
  labelling and re-evaluates its capacity directly. If U_ij = 0 the pair can never carry heat and has no variable.
- *W1_c{j}_{s}_{T}* (window of a cold stream, QC_jT > 0):
  Σ_i min(D_j[s..T], rc_i(s−1) + H_i[s..T]) y_ij ≥ D_j[s..T].
  *Proof:* D_j[s..T] = Σ_i Σ_{s≤t≤T} q_ijt, and the i-th inner sum is ≤ D_j[s..T], ≤ r_{i,s−1} + H_i[s..T] (heat i
  delivers in the window), and 0 if y_ij = 0.
- *W2_h{i}_{s}_{T}* (window of a hot stream, QH_is > 0, right-hand side > 0):
  Σ_j min(H_i[s..T], D_j[s..T]) y_ij ≥ H_i[s..T] − rc_i(T).
  *Proof:* decompose the chain flow of hot stream i into origin→delivery amounts g(a→t) ≥ 0, a ≤ t, with
  Σ_t g(a→t) = QH_ia and Σ_a g(a→t) = Σ_j q_ijt (possible for any flow on a path, e.g. first-in-first-out), and
  split each delivery among cold streams proportionally to q_ijt. Let O_j be the heat originating in [s..T] and
  delivered to j inside [s..T]: O_j ≤ min(H_i[s..T], D_j[s..T]) y_ij. Heat originating in [s..T] and delivered after
  T is part of the carry r_iT, so Σ_j O_j ≥ H_i[s..T] − r_iT ≥ H_i[s..T] − rc_i(T).
- *P_<row>* (partner count): if Σ_e a_e y_e ≥ b is a W1/W2 row with y binary, at least p of its y_e equal 1, where p
  is the least number of coefficients (largest first) whose sum reaches b; so Σ_e y_e ≥ p.
- *maxm:* Σ y ≤ K is the hypothesis being refuted.

All rows are generated by `scripts/hen_lpfile.py` (`write_lp_strong`) with exact decimal coefficients, and
**re-derived separately** from the data by `scripts/check_lp_model.py` (exact comparison, PASS/FAIL), which also
checks that every omitted variable is forced to 0 by the bounds above.

## 4. Certification of "no network with ≤ K matches" (VIPR)

1. SCIP 10.0.2 built from source with exact rational solving (SoPlex 8.0.2 with GMP/MPFR/Boost), run with
   `exact/enable = TRUE`, **cutting planes off** (`separating/maxrounds = 0`, `separating/maxroundsroot = 0`),
   **presolving off** (`presolving/maxrounds = 0`) and **conflict analysis off** (`conflict/enable = FALSE`)
   (`certificates/vipr/scip_exact_safe.set`), single thread, on the LP file of Sec. 3 (objective Σ y, plus the row
   maxm: Σ y ≤ K). Output: a VIPR certificate (`RTP infeas`) of a pure LP-based branch-and-bound proof.
   - With cutting planes on at all nodes, the Gomory mixed-integer cut derivations written by SCIP 10.0.2 could
     not be completed and verified (`Failed to derive constraint L6380_gom430_x9` on the 10sp1 test, with or without
     presolving, before and after canonicalising all rationals). With Gomory cuts at the root node only the 10sp1
     certificate verified, but the unbalanced10 certificate did not (`Failed to derive constraint
     L5381_gom21_s1484`, a root cut). We therefore use **no cutting planes at all**.
   - **Conflict analysis must be off** (`conflict/enable = FALSE`). Our first unbalanced10 run (conflict analysis on,
     21 203 nodes, 2.19 GB certificate) failed verification at a conflict-analysis row: `AggrRow_5376` claims
     Σ_{116 of the 122 y} y ≤ 18, derived "weakly" from the row Σ y ≤ 24 with the instruction `{ 0 }` (complete
     using *global* bounds). This is valid only under the local bounds y = 1 of the six dropped variables, so the
     derivation as written is not valid; all 9 828 weak derivations of that certificate are such conflict rows. With
     conflict analysis off, the certificates contain no weak derivations at all.
2. `scripts/check_vipr_model.py` checks that every constraint of the certificate's problem section equals (up to a
   positive factor) a row or bound of the LP file, every integer variable is a binary of the LP file, and the RTP is
   `infeas`. So the certificate proves infeasibility of (a relaxation of) the LP-file MILP.
3. `viprcomp` (VIPR, github.com/scipopt/vipr, commit 30f2951) completes the "weak" derivations; it is untrusted.
4. `viprchk` checks the completed certificate in exact GMP rational arithmetic. **Patch:** GMP's C++ stream input
   does not canonicalise rationals ("5/5" is read as 5/5, and GMP compares non-canonical values incorrectly: 5/5 ≠ 1),
   while SCIP prints non-canonical fractions; unpatched `viprchk` therefore rejects even correct certificates
   ("Noninteger value for integer variable"). We add `canonicalize()` after each rational is read (4 places; diff in
   `certificates/viprchk_canonicalize.patch`). Canonicalisation does not change any value, so the patch cannot make
   an invalid certificate pass. The patched checker accepts all VIPR example certificates that the original accepts
   and rejects the two deliberately incomplete examples.

5. **Hardening patch** (25 Oct 2026, `certificates/viprchk_harden.patch`, applied after the canonicalisation patch).
   Reading `viprchk.cpp` (commit 30f2951; master 99eec19 of 24 Oct 2026 has the same code) shows two inference
   rules that are checked too weakly:
   - `Constraint::round()` (rounding step `rnd`) requires integer coefficients on integer variables but silently
     accepts nonzero coefficients on *continuous* variables; Chvátal–Gomory rounding is invalid then.
   - `canUnsplit()` (step `uns`) checks that the branching constraints are a x ≤ β and a x ≥ β + 1 with a integral on
     integer variables, but not that β is an integer; with fractional β the disjunction misses integer points.
   The toy certificates `certificates/vipr/viprchk_tests/bad_rnd_continuous.vipr` (x ≥ 1/2, x ≤ 1/2, x continuous,
   "rounded" to x ≥ 1) and `bad_uns_fractional.vipr` (x = 1 integer, branching x ≤ 1/2 | x ≥ 3/2) "prove" feasible
   problems infeasible and are **accepted by the unpatched viprchk** (also with the canonicalisation patch only). The
   patch adds both checks, plus bounds checks on the referenced constraint indices (done by the original only in
   debug builds, otherwise out-of-range indices are undefined behaviour). It only adds rejections; the patched
   checker rejects both toy certificates, still accepts the VIPR examples the original accepts, and rejects the
   incomplete ones.
6. **Second checker** `scripts/vipr_check_independent.py`, written from the VIPR 1.1 specification without using
   viprchk's code (Python, exact rationals via gmpy2): asm / lin / rnd (support must be integer variables with
   integer coefficients) / uns (integral β required) / no sol, weak or incomplete steps; dominance as in the
   specification; anything after the reference index on a derivation line (SCIP writes a `global` flag) is ignored,
   i.e. never used to drop assumptions; success iff an assumption-free falsehood is derived.
7. **Negative tests** `scripts/vipr_mutation_tests.sh`: four corrupted copies of the unbalanced10 certificate
   (maxm rhs 24 → 25; one wrong multiplier; swapped branching assumptions; an over-strong rounding) must be rejected
   by both checkers.

Trusted base for a VIPR-certified claim: the proofs of Secs. 1 and 3, `check_lp_model.py`, `check_vipr_model.py`,
and either `viprchk` (+ both patches) with GMP or `vipr_check_independent.py` with gmpy2, and the C++/Python
toolchains. SCIP, SoPlex, viprcomp and the LP generator are not trusted.

**Result for unbalanced10 (25 Oct 2026).** Certificate of "no network with ≤ 24 matches": 540 587 890 bytes, sha256
`17967e67e607271c2fc60fc77a185da2badcbb84459d7deaeaddd5a4601ec352`, 191 697 derivations (10 232 asm, 170 552 lin,
5 797 rnd, 5 116 uns), no weak/incomplete step; SCIP: infeasible after 10 233 nodes (same count as on 24 Oct with a
different build). `verify_vipr.sh`: steps 1–5 PASS (viprchk 68 s on the raw file, 64 s after viprcomp; separate
checker 199 s; 7 min 40 s in total). Also accepted by the unpatched and canonicalise-only viprchk builds. Mutation
tests: see `certificates/vipr/mutation_tests.log`. With the 25-match network (`certificates/solutions/unbalanced10.json`,
checked by `check_solutions.py` and `check_network_lkm18.py`), **unbalanced10 = 25** (solved; see `verification.md`).

## 5. Floating-point pitfalls observed (why exact checks matter)

- 37sp-yfyv: the 36-match network returned by floating-point SCIP 10 (pyscipopt, default tolerances) is infeasible
  in exact arithmetic (max-flow deficit exactly 0.1 heat units, the resolution of the data). All six CPLEX/Gurobi
  match lists in the repository for this instance contain 37 pairs with y > 1e−10 and none of the 37 can be removed
  exactly. An exactly feasible 36-match network exists nevertheless (found through Lemma 2: two spanning trees on
  the 2-block partition of sizes 11 and 27).
- The repository's solution writer prints `y=1` for every y > 1e−10, so its match lists may contain pairs carrying
  negligible heat; we always recompute exact flows by integer max-flow and prune unnecessary pairs.

## 6. Assumptions and non-rigorous steps

- The instance is the data file as printed (decimal numbers taken as exact). The underlying HEN design data
  (temperatures, FCp, minimum-utility LP) were not re-derived; the min-utility-cost step of Letsios et al. is taken
  as given.
- Lower bounds of earlier work quoted in the README come from floating-point MILP solvers with tolerances; they are
  not certified (we reproduce some and certify others).
- Floating-point SCIP/HiGHS runs are used only to find solutions and to explore; every upper bound is an explicit
  network checked exactly, every lower bound is Lemma 1 (exact enumeration) or a VIPR certificate.

## 7. State of the work (25 Oct 2026)

- **Certified (exact):** 37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22 (Lemma 1 enumeration + exact networks);
  **unbalanced10 = 25** (VIPR certificate checked by two checkers + exact network); LB 32 for balanced15 and
  balanced15_random0/1/2; LB 26 (balanced12 family), 23 (22sp1), 32 (unbalanced15) by Lemma 1; all 22 upper bounds of
  `certificates/solutions/` checked exactly by two separate programs (`check_solutions.py`,
  `check_network_lkm18.py`, the second against the original 2018 MILP with its three big-M files).
- **History of unbalanced10 (24 Oct 2026).** Three exact-mode SCIP runs all reported "infeasible": run 1 (conflict
  analysis on, 21 203 nodes) and run 2 (root Gomory cuts, 10 129 nodes) produced certificates that FAIL `viprchk` at a
  conflict row and a Gomory row respectively (logs `certificates/vipr/run1_*`, `run2_*`; steps 1–3 PASSED for both).
  The cut-free, conflict-free run (10 233 nodes) was VIPR-checked on 25 Oct after re-running it with a fresh build
  (Sec. 4, item 7 and the result paragraph). Floating point SCIP also reports infeasibility (7 505 nodes).
- **Tried on 25 Oct, not closed (exploration only, no claim):** floating-point SCIP 10.0.2 with default settings
  (cuts on) on the strengthened model, 1800 s limit: unbalanced15 ≤ 35, 22sp1 ≤ 24, balanced12 family ≤ 27
  (balanced12_random1 ≤ 28) undecided (see README). Cut-free root LP bounds of the strengthened model: unbalanced15
  28.9 (36 needed), 22sp1 18.2 (25 needed), balanced12 22.2 (28 needed), against 21.1 (25 needed) for unbalanced10;
  the certified cut-free pipeline was therefore not run on them.
- **Not attempted:** unbalanced17/20 families beyond the upper bounds, large-scale instances.
