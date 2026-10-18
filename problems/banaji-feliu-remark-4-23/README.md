# Banaji–Feliu, Remark 4.23: are the equilibrium bounds for (n, n, n−1) networks sharp?

**Status: solved (preliminary): complete proof, awaiting independent verification.**

## Problem

M. Banaji and E. Feliu, *Positive equilibria in mass action networks: geometry and bounds*, J. Math. Biol. 93
(2026), doi:10.1007/s00285-026-02429-8, arXiv:2409.06877v5. Corollary 4.21(i): (n, n, n−1) networks admit no more
than n positive nondegenerate equilibria on any stoichiometric class. Remark 4.23 (p. 30 of v5):

> "We leave open whether the bound on positive nondegenerate equilibria in Corollary 4.21(i) is achieved for
> general n and, if so, the minimum source molecularity of a network needed to achieve the bound."

## Result

**Yes for every n, and the minimum source molecularity is exactly 3.** For each n ≥ 3 the network 𝒩ₙ with
sources X₁ and 2Xᵢ₋₁ + Xᵢ (i = 2, …, n) and explicit products (product molecularity ≤ 5) has, for suitable rate
constants, exactly n positive nondegenerate equilibria on one stoichiometric class. Since quadratic (n, n, n−1)
networks have at most 2 (Banaji–Feliu, Cor. 5.16), molecularity 3 is minimal; n = 2 is covered by the authors'
own example. The proof is short and computer-free: the reaction vectors form a "path" through the species, the
equilibria lie on the monomial curve x = C∘λ^q with q = (1, −1, 3, −5, 11, …), and the conservation law
x₁ + x₂ − x₃ − x₄ + x₅ + x₆ − … = T becomes a signomial whose coefficient signs alternate, so any n prescribed
roots can be realised (Descartes–Laguerre).

Smallest instance (n = 3): X₁ → 2X₁ + X₃, 2X₁ + X₂ → X₁ + 2X₂, 2X₂ + X₃ → X₂ with rates (1/25, 1/22500, 1/1296)
has exactly the equilibria (25, 36, 1), (50, 18, 8), (75, 12, 27) on x₁ + x₂ − x₃ = 60 (two stable, one saddle).

## Reproduce

```bash
python scripts/certify_family.py 3 12   # exact rational checks of the whole construction, < 1 s
python scripts/example_n3.py            # the n = 3 instance with sympy
```

## Files

| File | Content |
|---|---|
| `notes.md` | Theorem, complete proof, explicit example |
| `scripts/certify_family.py` | Exact (python-flint) verification for n = 3…12 |
| `scripts/example_n3.py` | Exact (sympy) verification of the n = 3 instance |
| `certificates/certify_family_n3-12.log` | Output: 301 checks, all PASS |
| `references.md` | Sources and novelty check |

## History

The n = 3 case was found during the exploration phase (17 Oct 2026), together with computer-found networks for
n ≤ 9; the uniform family 𝒩ₙ and the general proof were obtained on 18 Oct 2026.
