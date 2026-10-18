# Banaji–Feliu, Remark 5.8: can quadratic (n, n+2, n) networks reach n + 1 equilibria?

**Status: solved (preliminary): complete proof for every n, awaiting independent verification.**

**Attribution.** Construction, proof and manuscript: Claude Opus 5.5 (effort xhigh), orchestrating agent.
Independent verification: pending.

## Problem

M. Banaji and E. Feliu, *Positive equilibria in mass action networks: geometry and bounds*, J. Math. Biol. 93
(2026), doi:10.1007/s00285-026-02429-8, arXiv:2409.06877v5. Theorem 5.5: quadratic (n, n+k, n) networks with n + 2
distinct sources have at most n + 1 positive nondegenerate equilibria. Remark 5.8 (p. 41):

> "However, we were not able to find a quadratic (3, 5, 3) network admitting four positive nondegenerate equilibria,
> and thus achieving the bound in Theorem 5.5. We leave open the question of whether the bound in Theorem 5.5 can be
> achieved by quadratic (n, n + 2, n) networks for n ≥ 3."

## Result

**Yes, for every n ≥ 1.** The network 𝒩ₙ with sources 0, 2X₁, X₁+X₂, …, Xₙ₋₁+Xₙ, Xₙ,

  0 → Σ αᵢ₊₁Xᵢ, 2X₁ → 2X₁ + Σ βᵢ₊₁Xᵢ, Xᵢ + Xᵢ₊₁ → Xᵢ₊₁, Xₙ → 0,

has kernel exactly span(α, β), so the kernel can be prescribed. Equilibria then correspond to the roots of an
explicit one-variable function, and geometrically separated ratios αⱼ/βⱼ = 64^(j−1) give exactly n + 1 positive
nondegenerate equilibria, by a short estimate that is uniform in n.

Small example (n = 3, products of molecularity ≤ 20): 0 → X₁+3X₂+16X₃, 2X₁ → 16X₁+3X₂+X₃, X₁+X₂ → X₂,
X₂+X₃ → X₃, X₃ → 0 has four positive nondegenerate equilibria for a suitable rate constant of the second reaction.

**Side remark.** For n = 1 the family gives 0 → X, X → 0, 2X → 3X, with two positive nondegenerate equilibria. This
contradicts Remark 5.6 of the same paper, which misquotes Lemma 20 of Banaji–Boros–Hofbauer (Nonlinear Dyn. 2024): the
lemma lists exactly this network as the basic quadratic rank-one network with multiple equilibria. The negative
statement holds for *bimolecular* rank-one networks (their Corollary 19).

## Reproduce

```bash
python scripts/certify_family.py     # exact checks for the n = 3 example, n = 1, and the R = 64 family n = 1..8 (~10 s)
python scripts/other_example_n3.py   # an unrelated (3,5,3) example found by computer search
```

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: theorem for all n, proof, examples |
| `scripts/certify_family.py` | Exact (sympy) + 80-digit checks of the family |
| `scripts/other_example_n3.py` | Exact check of the computer-found example |
| `certificates/*.log` | Outputs (all PASS) |
| `references.md` | Sources and novelty check |
