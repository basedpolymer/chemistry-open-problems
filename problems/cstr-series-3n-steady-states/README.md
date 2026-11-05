# n cooled CSTRs in series: 3ⁿ steady states, 2ⁿ of them stable

<!-- plain-words:start -->
**The question.** A first-order exothermic reaction runs in n stirred tank reactors in series. Each reactor has its
own cooling jacket. One such reactor can have three steady states. In 1987 Bisbis showed that two reactors can have
nine. He conjectured that n reactors can have at most 3ⁿ, of which at most 2ⁿ are stable.

**The answer.** True, and both bounds are reached for every n. For n identical reactors (same volume, same jacket,
same coolant temperature) there are explicit operating parameters with exactly 3ⁿ steady states, exactly 2ⁿ of them
stable.

**How.** The upper bounds follow from the series structure: each reactor sees only its predecessor, has at most three
states for any feed, and its middle state is a saddle. For the construction, the Arrhenius rate cannot exceed k₀.
With a small k₀V/F, an ignited reactor still converts only a small fraction of its feed, so every downstream reactor
gets almost the same amount of reactant whatever happened upstream. Strong cooling makes each reactor almost
insensitive to the temperature of its feed. Five elementary inequalities finish the proof. The ignited states are
far hotter than any real reactor (above 3E/R).

<!-- plain-words:end -->

**Status: solved (preliminary).** Verification pending.

**Attribution.** Solver: Claude Opus 5.5 (effort xhigh). Verification: pending.

## Problem

B. Y. Bisbis, *On the steady states and dynamic behavior of two CSTR's in series*, PhD thesis, University of
Missouri–Columbia, December 1987 (supervisor P. C.-H. Chan), doi:[10.32469/10355/98283](https://doi.org/10.32469/10355/98283)
(open access, read in full). Chapter II, Concluding remarks, p. 27:

> "For the case of n - CSTR's in series, one may conjecture that reiteration of the same methods used here, along
> with an induction reasoning, would lead to the following : The maximum possible number of steady states is 3ⁿ , of
> which at most only 2ⁿ are stable, the rest unstable.."

The model is eqs. (5)–(8) of the thesis (p. 4), extended to n reactors in series: one feed at T₀ and one flow rate.
The kinetic parameters B, Da₀ and ε = RT₀/E are common to all reactors ("chosen only once", p. 29), and so is τ. The
volume ratio α_k = V_k/V₁, the jacket coefficient β_k and the coolant temperature η_k are chosen per reactor.

## Result

**Theorem.**
1. For every n and all parameter values (also in the exponential limit ε = 0), the cascade has at most 3ⁿ steady
   states. At most 2ⁿ of them are stable. Any steady state in which some reactor sits on its middle branch has a
   positive real eigenvalue.
2. For every n ≥ 1 and every ε > 0, put m = max(1, ε) and take
   B = 6400·n·m/ε², Da₀ = e^(−1/ε)/(4n), τ = 1, α_k = 1, β_k = 200m − 1, η_k = (1 − 20ε)/(20ε²) for every k.
   Then there are exactly 3ⁿ steady states, all hyperbolic. Exactly 2ⁿ are asymptotically stable (sinks); the other
   3ⁿ − 2ⁿ are unstable. In physical terms: identical reactors with k₀V/F = 1/(4n), UA/(FρC_p) = 200m − 1, coolant
   at E/(20R), adiabatic rise 6400·n·m·E/R, and any feed temperature.
3. The outer states need not be stable. Two certified one-reactor examples have three steady states each: in the
   first only one state is stable (the upper one is an unstable focus), in the second none is.

The proof of part 2 is analytic and holds for every n. Interval certificates confirm it independently for n ≤ 6 with
ε = 1/20 (up to 729 states, 64 stable) and for n ≤ 4 with ε = 1 and ε = 5. They use a Krawczyk test on the thesis's
own equations. With ε = 1/20 and n = 3 the thesis parameters are B = 7 680 000, Da₀ = e^(−20)/12, τ = 1, α_k = 1,
β_k = 199, η_k = 0. This gives 27 states, 8 of them stable.

**Caveats.** The attaining parameters are far outside practical operation: ignited reactors run between 3E/R and
9E/R. Whether 3ⁿ is reached within realistic ranges is not addressed. The construction needs ε > 0, because it relies
on the bounded Arrhenius factor. The exponential limit ε = 0, used in the thesis's numerical example for n = 2, is
covered by the upper bounds but not by the construction. Varma (1980) and Li (1994) are paywalled and were checked
through abstracts and citing works only, as were Kubíček et al. (1980) and Retzloff et al. (1992), which treat two
reactors (see [references.md](references.md), [NEEDS.md](NEEDS.md)).

## Reproduce

From the repository root (WSL, `~/envs/sci`):

```
wsl -e bash -lc 'source ~/envs/sci/bin/activate && cd /mnt/c/Users/demor/chemistry-open-problems/problems/cstr-series-3n-steady-states && python scripts/check_model_algebra.py && python scripts/check_construction_inequalities.py && python scripts/certify_cascades.py && python scripts/certify_unstable_outer_state.py'
py -3 templates/build_tex.py problems/cstr-series-3n-steady-states/paper.tex
```

Each script prints PASS/FAIL per claim and exits non-zero on failure. Total runtime about one minute.

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: model, upper bounds, window lemma, explicit family, unstable outer states |
| `scripts/check_model_algebra.py` | SymPy: change of variables, scalar reduction, convex–concave factor, Jacobian determinant and trace identities, block triangular structure |
| `scripts/check_construction_inequalities.py` | Arb: the five numerical inequalities behind the explicit family |
| `scripts/certify_cascades.py` | Arb: 3ⁿ steady states and 2ⁿ stable ones, Krawczyk-verified in the thesis variables (n ≤ 6 for ε = 1/20; n ≤ 4 for ε = 1, 5) |
| `scripts/certify_unstable_outer_state.py` | Arb: the two one-reactor examples with unstable outer states |
| `certificates/cascade_n3_eps1_20.json`, `cascade_n4_eps1_20.json` | Every steady state for n = 3, 4 (ε = 1/20): word, enclosures, block traces and determinants |
| `certificates/cascade_summary.json`, `unstable_outer_state.json`, `*.log` | Summaries and logs of all runs |
| `references.md` | Sources and novelty check |
| `NEEDS.md` | Paywalled papers useful for the priority check |
