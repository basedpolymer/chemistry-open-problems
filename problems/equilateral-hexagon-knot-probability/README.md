# Knotting probability of random equilateral hexagons: a certified two-sided enclosure

**Status: partial result (preliminary, not yet independently verified)** — a certified enclosure
L ≤ p₆ ≤ U of the knotting probability, with a first positive lower bound and an upper bound
**[RUN]** times smaller than the best published one; the exact value remains open.

**Attribution.** Solver (novelty check, lemmas and proofs, certifier, branch and bound, write-up): Claude Opus 5.5
(effort not recorded), sub-agent of the orchestrator, 24 Oct 2026. Verifier and editor: not yet assigned.

## Problem

Let p₆ be the probability that a random closed equilateral hexagon in R³ (unit edges, modulo translation and
rotation, with the standard = symplectic measure of Cantarella–Shonkwiler, Ann. Appl. Probab. 26 (2016) 549–596)
is knotted (necessarily a trefoil). The best rigorous bounds were 0 ≤ p₆ ≤ (14 − 3π)/192 ≈ 0.0238

> "Theorem 4.1. The probability that an equilateral hexagon is knotted is at most (14−3π)/192."
> — K. Hake, *Knotting probability of equilateral hexagons*, J. Knot Theory Ramifications 28 (2019) 1950040,
> arXiv:1810.11882, §4

while Monte Carlo gives about 1.3–1.4·10⁻⁴ (Cantarella–Shonkwiler §4.1: "1.3 × 10⁻⁴ ± 0.2 × 10⁻⁴"; Hake §4:
1.370402·10⁻⁴). The most recent paper on the subject leaves exact hexagon knotting probabilities open:

> "Returning to the original question about knot probabilities, it would be interesting to determine the exact
> probability of knotting, even for hexagons, in any other random knot model with a continuous state space."
> — J. Geisz, C. Peterson, C. Shonkwiler, *Trefoil probabilities and polyhedra*, arXiv:2609.24860 (Sept. 2026), §5
> (who solve the six-points-on-the-sphere model, p = 1/(16π²), and cite Hake's bound as the known result for
> equilateral hexagons).

Novelty check: `references.md` (no improvement of Hake's bound and no positive lower bound found).

## Result

> **Theorem (certified).** **[RUN]** ≤ p₆ ≤ **[RUN]**.

(Monte Carlo, 10⁸ samples, numerics only: p₆ ≈ **[RUN-MC]**.)

**Method in one paragraph.** In the action-angle coordinates of the triangulation T135 (diagonal lengths
d ∈ P, dihedral angles θ ∈ T³; uniform measure = Lebesgue × Haar), elementary arguments (Δ-moves, a
shared-vertex lemma, half-spaces; Section 2 of `result.md`) show that, up to the symmetries
θ ↦ −θ and S₃, every knotted hexagon with θ ∈ (0,π)³ satisfies one of two explicit *pierce patterns*
R = A12 ∧ A23 ∧ A13 or L = B12 ∧ B23 ∧ B13, where each condition involves only **two** of the three angles
(e.g. A12: "edge v4v5 pierces the triangle v1v2v3"). Conversely, wherever the R (or L) sign conditions hold
strictly the hexagon is embedded, so knottedness is constant on connected regions and is certified at one point
by an exact Fox 3-colouring of a certified diagram. The volume of the patterns is enclosed by an adaptive
branch and bound over boxes of the 3-dimensional d-chamber: for each box, interval arithmetic (outward rounding,
centred forms in d and θ, certified trigonometric table) classifies cells of three 2-dimensional angle grids, and
the number of 3-dimensional cells in the intersection of the three "cylinders" is the exact trace
Σ P12[a,b]·|P23[b,·] ∧ P13[a,·]|. The bounds are exact dyadic rationals.

**Side result (Hake's lemmas).** Hake's Lemmas 3.9 and 3.12 (used in her proof of Theorem 4.1) are false:
`scripts/hake_lemma_check.py` certifies knotted hexagons of curl +1 with d1 > d2 > d3 and θ2 > π/2, and with
d1² > d2² + d3² and θ1 < π/2 (about 0.6 % and 1.4 % of all knots in the chamber, all left-handed in the Monte
Carlo sample). Her proof is therefore incomplete; her stated bound is nevertheless true (it exceeds our U). Her
remark "(14−3π)/192 < 1/42" is also false (0.023829 > 0.023810).

## Reproduce

```bash
bash scripts/build.sh                       # g++ (C++20), python3 + mpmath
python3 scripts/test_interval.py            # interval arithmetic vs exact rationals      PASS/FAIL
python3 scripts/check_geometry.py           # genericity, symmetries, chart (ball arithmetic, sympy)
python3 scripts/hake_lemma_check.py         # certified counterexamples to Hake's Lemmas 3.9/3.12
python3 scripts/check_boxes_mc.py build/hexcert build/trig14.bin     # certifier vs Monte Carlo (numerics)
python3 scripts/mc_validate.py 100000000 2026 L U                    # 10^8-sample Monte Carlo (numerics)
bash scripts/reproduce_small.sh             # short B&B run (5 rounds) + full re-verification
python3 scripts/verify_certificate.py certificates/hexknot.leaves.gz certificates/hexknot.summary.json \
        build/hexcert build/trig14.bin [--sample 5000]                # re-check the main certificate
```

The main certificate was produced by `scripts/bnb_driver.py` (resumable: the state file is rewritten after every
round, the same command resumes). **[RUN details]**. With k cores, the leaf evaluations are independent:
split the leaf list into k parts and run `verify_certificate.py`-style recomputation (or the driver's certifier
calls) in parallel; the exact sums are additive.

## Files

* `result.md` — statement, lemmas with proofs, the certified integration scheme and its error control, results,
  numerical validation, assumptions and limitations.
* `src/interval.hpp` (interval arithmetic, forward AD), `src/hexcert.cpp` (certifier for one d-box),
  `src/test_interval.cpp`.
* `scripts/` — `build.sh`, `make_trig_table.py`, `bnb_driver.py`, `export_certificate.py`,
  `verify_certificate.py`, `reproduce_small.sh`, `test_interval.py`, `check_geometry.py`, `check_boxes_mc.py`,
  `mc_validate.py`, `hake_lemma_check.py`.
* `certificates/` — `hexknot.leaves.gz` (all leaves: box, θ-level, counts), `hexknot.summary.json` (exact bounds,
  SHA-256), logs of all checks and of the run, `SHA256SUMS.txt`.
* `references.md`, `NEEDS.md` (paywalled items, none needed for the proof), `PROGRESS.md` (git-ignored log).
