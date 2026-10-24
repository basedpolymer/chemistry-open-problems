# Minimum number of matches: exact values on the Letsios–Kouyialis–Misener heat-recovery benchmarks

**Attribution.** Solver: Claude Opus 5.5 (effort not recorded), sub-agent of the orchestrator. No helper agents.
Verification: pending.

**Status: partial result (not yet independently verified; work stopped on 24 Oct 2026 at the owner's request).**
Of the 25 benchmark instances left open in 2018, three are closed with exact certificates (an exactly checked
network plus a combinatorial lower bound proved by exhaustive exact enumeration): **37sp-yfyv = 36, 20sp1 = 19,
23sp1 = 22** (the last with a new upper bound, previously 23). For **unbalanced10** (25 vs 24), exact-arithmetic
SCIP proves that no network with ≤ 24 matches exists, but the VIPR check of that run was not completed (two
earlier certificates failed VIPR verification because of SCIP certificate defects, see below), so unbalanced10 = 25
is *not yet certified*. The lower bound of the four balanced15 instances improves from 31 to **32** (certified).
All other lower bounds of the 2018 study are floating-point solver output; several are certified here (equal or
lower values), none is contradicted.

## Problem

D. Letsios, G. Kouyialis, R. Misener, *Heuristics with performance guarantees for the minimum number of matches
problem in heat recovery network design*, Comput. Chem. Eng. 113 (2018) 57–85, doi:10.1016/j.compchemeng.2018.03.002,
arXiv:1709.04688. Instances, models and logs: github.com/dimletsios/min_matches_heuristics.

Given hot and cold streams (utilities included, with loads fixed by the minimum-utility-cost step) and their heat in
each shifted temperature interval, find a feasible heat exchange (heat goes only to the same or colder intervals)
using the fewest hot–cold stream pairs ("matches"). This is the NP-hard MILP of the sequential HEN synthesis method
(transshipment model of Papoulias–Grossmann 1983, transportation model of Cerdá–Westerberg 1983).

**Open statement (arXiv v2, Sec. 9.2, p. 30), verbatim:** *"Observe that state-of-the-art approaches cannot, in
general, solve moderately-sized problems with 30-40 streams to global optimality. For example, none of the test cases
in the Grossmann (2017) or large scale test sets is solved to global optimality within the specified time limit."*
The runs used *"absolute gap 0.99, relative gap 4%"* (same section; the repository code sets only the 4 % relative
gap), so instances reported as solved may stop with a one-unit gap; from the repository logs, 37sp-yfyv (36 vs bound
35) and unbalanced10 (25 vs bound 24) are such cases.

## Results

Previous values: best upper bound / best lower bound over the 6 runs per instance (CPLEX 12.6.3 and Gurobi 6.5.2 on
three formulations) in `data/mip_solutions/` of the repository ("LKM18"), re-derived by us from the `.sol` files
(lower bounds rounded up); FS04 = Furman–Sahinidis 2004 (CPLEX 7.0, 7 h), quoted from LKM18 Table G.9. These lower
bounds are floating-point solver output, not certificates.

| instance | previous UB / LB | new UB / LB | status | UB certificate | LB certificate |
|---|---|---|---|---|---|
| 37sp-yfyv | 36 / 35 (LKM18; FS04 36) | **36 / 36** | **optimal = 36** | exact network (36) | Lemma 1, c_max = 2 (exact enumeration) |
| 20sp1 | 19 / 16 (LKM18; FS04 19) | **19 / 19** | **optimal = 19** | exact network (19) | Lemma 1, c_max = 2 |
| 23sp1 | 23 / 17 (LKM18; FS04 23) | **22 / 22** | **optimal = 22** (new UB) | exact network (22) | Lemma 1, c_max = 2 |
| unbalanced10 | 25 / 24 (LKM18) | 25 / 24 certified-by-nothing-new; ≤ 24 infeasible per exact SCIP | open (preliminary 25) | exact network (25) | exact SCIP 10.0.2 run: infeasible (10 233 nodes); VIPR check **not completed** |
| unbalanced15 | 36 / 35 (LKM18) | 36 / 32 certified | open | exact network (36) | Lemma 1 |
| 22sp1 | 25 / 23 (LKM18; FS04 25) | 25 / 23 certified | open | exact network | Lemma 1, c_max = 1 |
| balanced12 | 28 / 26 | 28 / 26 certified | open | exact network | Lemma 1, c_max = 1 |
| balanced12_random0 / 1 / 2 | 28 / 26, 29 / 27, 28 / 26 | 28 / 26, 29 / 26, 28 / 26 certified | open | exact networks | Lemma 1, c_max = 1 |
| balanced15 | 35 / 31 | **35 / 32** | open, **record LB** | exact network | Lemma 1, c_max = 1 |
| balanced15_random0 / 1 / 2 | 35 / 31 each | **35 / 32** each | open, **record LB** | exact networks | Lemma 1, c_max = 1 |
| unbalanced17, _random0 / 1 / 2 | 43 / 39, 43 / 39, 43 / 40, 43 / 40 | UB reproduced exactly | open | exact networks | – |
| unbalanced20, _random0 / 1 / 2 | 50 / 45, 51 / 45, 51 / 45, 50 / 45 | UB reproduced exactly | open | exact networks | – |
| large_scale0 / 1 / 2 | 175 / 135, 219 / 123, 239 / 131 | not attempted | open | – | – |

"certified" means: proved by an exact argument here (for balanced12_random1 the certified bound 26 is below the
uncertified solver bound 27). Lemma 1 (block lemma, `result.md` Sec. 2): every connected component of a feasible
match graph is heat-balanced and cascade-feasible, so #matches ≥ (number of streams) − c_max, where c_max, the
largest number of such blocks partitioning the streams, is computed by exhaustive exact enumeration.

Side findings:
- The 36-match network found by floating-point SCIP for 37sp-yfyv is infeasible in exact arithmetic (deficit 0.1),
  and every CPLEX/Gurobi match list in the repository for it needs 37 pairs exactly; an exact 36-match network
  exists (two spanning trees on a 2-block partition, `certificates/solutions/37sp-yfyv.json`).
- SCIP 10.0.2 exact-mode certificates failed VIPR verification when Gomory cuts were used (also root-only) and when
  conflict analysis was on (a conflict row valid only under local bounds was marked as globally derivable); only
  cut-free, conflict-free runs are usable. `viprchk` needed a canonicalisation patch (GMP reads "5/5"
  non-canonically and then compares it wrongly).

## Reproduce

```bash
cd problems/hen-minimum-matches
git clone https://github.com/dimletsios/min_matches_heuristics certificates/sources/min_matches_heuristics
python3 scripts/check_solutions.py      # exact check of all 22 networks (upper bounds), ~10 s, ALL PASS
bash scripts/run_block_bounds.sh        # Lemma 1 lower bounds (exact enumeration), ~1 min, ALL PASS
python3 scripts/side_findings.py        # 37sp-yfyv floating-point findings, ALL PASS
```

## How to resume (unbalanced10 and further instances)

- Build: SCIP 10.0.2 (`-DEXACTSOLVE=on`) + SoPlex 8.0.2 with GMP/MPFR/Boost; VIPR tools from github.com/scipopt/vipr
  (commit 30f2951) with `certificates/viprchk_canonicalize.patch` applied to `viprchk.cpp`; `viprcomp` needs TBB.
- Run: `scip -c "set load certificates/vipr/scip_exact_safe.set" -c "set certificate filename u10.vipr"
  -c "read certificates/vipr/unbalanced10_le24.lp" -c optimize -c quit` (≈ 5 min; certificate ≈ 0.54 GB; the
  solver log of the completed run is `certificates/vipr/unbalanced10_le24_scip.log`).
- Check: `VIPRCOMP=... VIPRCHK=... bash scripts/verify_vipr.sh unbalanced10 24 u10.vipr` (LP regeneration,
  independent row re-derivation, certificate-problem match, completion, patched viprchk). Steps 1–3 passed on the
  earlier runs; step 4 has not been run on the cut-free, conflict-free certificate.
- Next targets: unbalanced15 (≤ 35), 22sp1 (≤ 23), balanced12/15 families (≤ 26 / ≤ 32) with the same pipeline;
  LP files are generated by `scripts/hen_lpfile.py:write_lp_strong` with `hen_models.strong_rows`.

## Files

- `result.md` — definitions, Lemmas 1–2, validity proofs of every MILP row, certification pipeline, assumptions.
- `scripts/hen_instance.py` (exact parser, consistency checks), `hen_exact.py` (exact max-flow, exact solution
  checker), `block_bound.py` (Lemma 1 enumeration, PASS/FAIL), `check_solutions.py` (PASS/FAIL per network),
  `hen_cuts.py`, `hen_models.py`, `hen_lpfile.py` (strengthened MILP; exact LP files), `check_lp_model.py`
  (independent re-derivation of every LP row), `check_vipr_model.py` (certificate problem = LP model),
  `build_solution.py`.
- `certificates/solutions/*.json` — exact rational networks; `claimed_upper_bounds.json`.
- `certificates/block_bounds.log`; `certificates/vipr/` — LP model, SCIP settings and logs of the unbalanced10 runs
  (run1: conflict analysis on, VIPR FAIL; run2: root Gomory cuts, VIPR FAIL; final: cut/conflict-free, solver says
  infeasible, VIPR not checked). Certificates themselves (0.5–2.2 GB) are not stored.
- `certificates/side/` — 37sp-yfyv floating-point side findings and log.
- `certificates/viprchk_canonicalize.patch` — the patch to VIPR's checker.
- `references.md` (sources, novelty check), `NEEDS.md`.
