# Verification status

**Self-checks passed; independent verdict pending.**

Author of both the argument and these checks: GPT-6-Astra (effort xhigh), 31 October 2026.
This document is not a fresh-context verification report and does not use the verdict CONFIRMED.
The repository's appropriate status is **solved (preliminary)**: the mathematical argument has both
bounds, but no separate verifier has assessed the new work.

## What is complete

- Lower bound: every network with at most four reactant complexes is excluded under the stated
  bimolecularity, ACR and exponential-stability hypotheses. Rank at most two is covered by the
  cited earlier result; rank three is covered by the new proof and exhaustive exact certificates.
- Upper bound: a six-reaction network on four species has exactly five reactant complexes,
  global ACR and two exponentially stable positive equilibria in the class `a+b=9`.
- Thus the exact minimum is five. No unresolved five-versus-six alternative remains in this argument.

## Checks performed

`scripts/run_all.py` runs the four certificates and records their full logs, environment, exit codes
and source hashes in [certificates/execution.json](certificates/execution.json).

1. `certify_sources.py`: every four-element source set covering four through eight species, exact
   rank and ACR tests, canonicalization, all 26 four-species exponent bounds and all nine larger types.
2. `certify_subspaces.py`: all 40,385 rank-three spans for S1, S2, S3, S7 and S9; exact membership,
   dual-cone test, Cauchy-Binet injectivity/degeneracy tests, symbolic kinetic-rank exclusions,
   and membership of every remaining space in the nineteen proved normal forms.
3. `certify_normal_forms.py`: every scalar curve preserves its conservation laws and one binomial;
   the other binomial equals the stated scalar equation times a positive factor. The augmented
   determinant identities and concavity derivatives are checked symbolically.
4. `certify_five_example.py`: reconstructs the vector field from the reaction list, counts sources
   and complexes, computes the stoichiometric rank and conservation law, proves global ACR,
   factors the complete scalar equilibrium equation, and checks all stability determinants and
   traces exactly. The coupled variant has the same restricted characteristic polynomial.

The handwritten parts were reviewed for the maximal-network reduction, the zero dual cone criterion,
the determinant orientation, removal of inert species, and connectedness of every positive parameter
domain. The S4/S5/S6/S8 quadratic argument uses a linear three-group embedding and an affine-line
intersection, with a fixed determinant orientation. The rank-three classification allows arbitrary
numbers of reactions from a source.

## Corrections during self-checking

- Some scalar determinant identities differed from initial guesses by a fixed orientation sign.
  The proof only compares signs at two roots within the same system; the exact factors are now
  recorded and every form is checked up to its fixed orientation.
- The five-source factorization check originally compared two differently expanded SymPy expressions
  structurally. It now checks that their rational-function difference is exactly zero.
- These were development corrections, not an independent verification process.

## Remaining review, not an unresolved mathematical subcase

A separate verifier should examine the hand arguments and reimplement or audit the exhaustive
enumeration. The lower bound also depends on the parent manuscript's dimension-two theorem, whose
existing verification is in [../verification.md](../verification.md). No human review or claim of
research priority is implied. Numerical exploration files are not inputs to any certificate.
