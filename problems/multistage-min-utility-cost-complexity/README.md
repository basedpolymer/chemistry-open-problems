# Complexity of the multistage minimum utility cost problem

<!-- plain-words:start -->
**The question.** In a heat-exchanger network, hot streams pass through a fixed number of stages and give heat to cold
streams; streams may split and remix in each stage. How hard is it to choose the splits so that as little external
heating and cooling as possible is needed? Letsios and co-authors (2020) call this an intriguing open question.

**The answer.** As printed, the paper's model forgets that heat only flows from hot to cold, and is then trivial. Once
that physical constraint is added, in either of its two usual forms, the problem is strongly NP-hard, already with two
stages and already for the question "can the utility cost be zero?". With one stage, one form is a linear program and
the other is still NP-hard.

**How.** An exact identity on squared temperatures shows that a zero-utility design with balanced energy and balanced
"squared-temperature budget" can have no mixing loss and no temperature-driving-force slack. The temperatures are then
forced onto a grid, and each hot stream must pick one "bin": a 3-Partition instance.

<!-- plain-words:end -->

**Status: solved (preliminary).** Verification pending.

**Attribution.** Solver and manuscript: GPT-6-Astra (effort xhigh), two runs on 5 November 2026 (the second one stopped at
its usage limit after the manuscript and certificates were written). Deliverables completed (README, credits, compilation,
re-run of the certificates): Claude Opus 5.5 (effort xhigh). Verification: pending.

## Problem

D. Letsios, R. Baltean-Lugojan, F. Ceccon, M. Mistry, J. Wiebe, R. Misener, *Approximation algorithms for process systems
engineering*, Computers & Chemical Engineering 132 (2020) 106599, doi:10.1016/j.compchemeng.2019.106599
(arXiv:1909.12328), §6.4, p. 10:

> "The complexity of the multistage minimum utility cost problem is an intriguing open question."

and §7, future direction 5, p. 11: "Determining the computational complexity of the multistage minimum utility cost
problem is an intriguing future direction. Because of stream mixing, the problem exhibits commonalities with pooling.
However, no hardness reduction formalizes this insight of domain experts."

The model is the QCQP (5) of the paper: hot streams cross stages 1, …, ℓ, cold streams cross them in the opposite
direction, every stream splits among all streams of the other kind in each stage and remixes non-isothermally, and the
objective is the hot plus cold utility cost. The journal version (checked on a rendering of p. 10) corrects the
monotonicity signs of arXiv v1 but, like arXiv v1, has **no temperature-approach constraint** relating hot and cold
temperatures; (5m) also indexes cold streams by the hot set, read here as the cold set.

## Result

**Theorem** (paper.tex, Theorems 1–2 and §§5–6).
1. *Literal model.* The optimum of (5) as printed is c_CU·max(E_H − E_C, 0) + c_HU·max(E_C − E_H, 0), with an optimal
   rational solution computable in polynomial time. The model allows heat to flow from a colder hot branch to a hotter
   cold branch, so it cannot be the intended problem.
2. *Physical completions.* Add, for every exchanger that carries heat, either the counter-current endpoint approaches
   (model C: hot inlet ≥ cold branch outlet + δ, hot branch outlet ≥ cold inlet + δ) or the stage-boundary approaches of
   Yee–Grossmann's SYNHEAT (model B). For every fixed number of stages ℓ ≥ 2, deciding whether the minimum utility cost is
   zero is **strongly NP-hard** in both models, with all hot–cold pairs allowed, unit utility prices, integer data and any
   fixed approach δ ≥ 0. Hence both are strongly NP-hard when ℓ is part of the input (unary or binary).
3. *One stage.* Model C is a linear program (polynomial). Model B is strongly NP-complete.
4. *Membership.* Model B's decision version is in NP (so strongly NP-complete); model C is in ∃ℝ, and its NP membership
   is not settled. Any number of stages can be compressed to at most 2(n + s) − 1 without changing the attainable
   targets, so these memberships also hold for a binary stage count.

Key step: for each exchanger, f(h² − a²) − g(b² − c²) = q[(h − b) + (a − c)]; summed over a network, the difference of
the hot and cold "squared-temperature budgets" equals the approach slacks plus the mixing variances, all nonnegative. In
the reduction both budgets balance exactly, so every slack and variance vanishes, active matches share their stage-boundary
temperatures, and each hot item stream must drop to a single integer level j, i.e. join bin j.

**Scope.** The answer concerns these two explicit completions, which are the standard ones (SYNHEAT and the counter-current
exchanger), not every conceivable engineering superstructure. Strong NP-hardness of C does not depend on its NP
membership, which stays open.

## Reproduce

From the repository root (WSL, `~/envs/sci`; each script prints PASS/FAIL per claim and exits non-zero on failure;
runtime a few seconds):

```bash
wsl -e bash -lc 'source ~/envs/sci/bin/activate && cd /mnt/c/Users/demor/chemistry-open-problems/problems/multistage-min-utility-cost-complexity && python scripts/verify_reduction.py && python scripts/verify_model_audit.py'
py -3 templates/build_tex.py problems/multistage-min-utility-cost-complexity/paper.tex
```

The proofs do not depend on these computations; the scripts check the polynomial identities symbolically, 36 rational
zero-utility witnesses (2–5 stages, δ = 0, 1, 3/2) against every equation and both approach versions, the equality-case
enumeration on small 3-Partition instances, and the one-stage multiplier identities. `scripts/explore_physical.py` is a
floating-point exploration with SCIP; its log contains unresolved runs (time limits) and is not a certificate.

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: model, literal theorem, quadratic identity, reductions, one stage, membership |
| `scripts/verify_reduction.py` | Exact checks of the reduction (SymPy, fractions); log `certificates/reduction-checks.log`, data `reduction-checks.json` |
| `scripts/verify_model_audit.py` | Exact audit of the literal model and the one-stage LP example; `certificates/exact-checks.json` |
| `scripts/explore_physical.py` | Floating-point exploration (not a certificate); `certificates/physical-exploration.*` |
| `scripts/novelty_*.py`, `inspect_audit.py` | Read-only novelty queries (network-dependent) |
| `references.md` | Sources read, novelty audit and its limits |
| `NEEDS.md` | Source request (resolved) |
