# Kaihnsa–Nguyen–Shiu: does ACR with multistability need larger networks than ACR with multistationarity?

<!-- plain-words:start -->
**The question.** A species has “absolute concentration robustness” if it takes the same value in every steady state, an ideal regulation mechanism. Kaihnsa, Nguyen and Shiu conjectured that combining it with two stable states requires larger networks than combining it with several steady states.

**The answer.** True for dimension (≥ 3), reactant complexes (≥ 4) and complexes (≥ 6); false for the number of species (3 are enough).

**How.** A normal form shows that in dimension 2 robustness allows only one stable state. A 3-species network with two stable states is then built and certified exactly.

<p align="center"><img src="../../figures/acr-multistability.svg" alt="Robust concentration with two stable steady states" width="720"></p>

<p align="center"><sub>One combination of the equations forces x = 1 at every steady state, yet two stable states coexist.</sub></p>

More detail, with every step explained: [guided reading, result 27](https://basedpolymer.eu/chemistry-open-problems#result-27) · [en français](https://basedpolymer.eu/chemistry-open-problems/fr#d27).

<!-- plain-words:end -->

**Status: solved; verified (verdict CONFIRMED WITH CORRECTIONS, see [verification.md](verification.md)).** The
conjecture holds for three of the four quantities (dimension, reactant complexes, complexes) and fails for the number
of species.

**Attribution.** Solver: Claude Opus 5.5 (effort xhigh), 27 Oct 2026: novelty check, dimension-two theorem,
consequences for reactant complexes and complexes, manuscript. Sub-studies by Claude Opus 5.5 (effort xhigh): certified
examples, reactions and reactant complexes (`examples/`); minimal number of complexes (`complexes/`). Verification:
Claude Opus 5.5 (effort xhigh), fresh context, 27 Oct 2026 ([verification.md](verification.md), `verification/`).

## Problem

N. Kaihnsa, T. Nguyen, A. Shiu, *Absolute concentration robustness and multistationarity in reaction networks:
conditions for coexistence*, Eur. J. Appl. Math. 35 (2024) 566–600, doi:10.1017/S0956792523000335, arXiv:2307.04186.
Section 6 (Discussion), p. 34 of arXiv v1:

> "In the present work, our interest in multistationarity comes from the fact that it is a necessary condition for
> multistability. Another interesting direction, therefore, is to investigate conditions for coexistence of ACR and
> multistability, rather than multistationarity. The "minimal" networks in the current work admit only two positive
> steady states and are not multistable. Hence, we conjecture that the lower bounds (on dimension and the numbers of
> species, reactant complexes, and all complexes) for the coexistence of ACR and multistability are strictly larger
> than the bounds proven here for ACR and multistationarity."

The bounds proven there (Theorem 1.2) for at-most-bimolecular networks with ACR and nondegenerate multistationarity
are 3 species, 3 reactant complexes, 5 complexes and dimension 2, all attained. Definitions: ACR as in their
Definition 2.6; multistability means that, for the same rate constants, one stoichiometric compatibility class
contains two exponentially stable positive steady states (the usual notion, e.g. Jiao–Tang–Zeng, J. Math. Biol. 2025).

## Result

| | ACR + multistationarity (KNS) | ACR + multistability: lower bound | best example | conjecture |
|---|---|---|---|---|
| dimension | 2 | **3** | 3 | true |
| species | 3 | 3 | **3** | **false** |
| reactant complexes | 3 | **4** | 6 | true (minimum 4, 5 or 6: open) |
| complexes | 5 | **6** | 6 | true |

- **Dimension two (main theorem).** A bimolecular mass-action system of dimension 2 with ACR has at most one
  exponentially stable positive steady state in each compatibility class. Proof: on the ACR hyperplane the vector field
  is h·α for a fixed vector α, which gives a normal form f = ψα + (x₁−c)Pℓ and det J = −P·(∇ψ·ℓ); two stable states
  force P to change sign along the positive steady-state set; bimolecularity of the products makes that set connected
  and smooth; global ACR forces every zero of P on it to be a tangency, so P cannot change sign. With KNS Prop. 3.10
  (dimension 1) this gives dimension ≥ 3, and with KNS Thm 5.1/5.5 at least 4 reactant complexes.
- **Complexes.** If a network has a nondegenerate positive steady state and a nonterminal complex, then δ ≥ t − ℓ + 1
  (kinetic-subspace argument). With dimension ≥ 3, five complexes force one linkage class of deficiency one with one
  terminal class, excluded by Feinberg's deficiency-one theorem. Hence at least 6 complexes. (This also closes a small
  gap in KNS: their Lemma 2.5 states that theorem for weakly reversible networks only but is applied to general ones.)
- **Three species suffice (counterexample to the species part).** The network
  N1 = {Z → Y+Z, 2Z → Z, Y → Y+Z, Y+Z → 2X, X+Z → 2Z, 2X → Y} has ACR in X for every rate constant vector
  (ẋ + 2ẏ = z(2k₁ − k₅x)); for κ = (5, 1, 1, 4, 10, 7) its positive steady states are exactly (1, 19/4, 1/2) and
  (1, 7/4, 7/2), exponentially stable, and (1, 3, 1), a saddle (exact Routh–Hurwitz). It has 3 species, 6 reactions,
  6 complexes and dimension 3, so three of the four minima are attained by one network.
- **Further results** (manuscript §6): a second certified example G7 with rational steady states; a weakly reversible
  six-complex example with two linkage classes; with 3 species at least 6 reactions are needed (such a network is
  full-dimensional, and an exact enumeration of all 3 188 646 networks with 5 reactions and 5 reactant complexes
  shows that no full-dimensional one has two stable states); with exactly 4 reactant complexes, 4 species are
  impossible and 6 of the 9 possible structures for more species are excluded by proof (the other 3 numerically).

**What is computer-assisted.** The examples (exact rational arithmetic, ball-arithmetic eigenvalue enclosures), the
5-reaction enumeration, and the list of 4-reactant structures. The dimension-two theorem and the bounds on dimension,
reactant complexes and complexes are proved by hand. **Numerical only:** the exclusion of the reactant structures
S3, S7, S9 (not needed for the conjecture) and the searches that found the examples.

**Open.** The exact minimal number of reactant complexes (4, 5 or 6); whether the dimension-two theorem extends to
networks with bimolecular reactants but larger products.

## Reproduce

From this folder, under WSL Fedora 44 with Python 3.12, sympy 1.14, python-flint 0.9:

```bash
python scripts/check_dimension_two_example.py          # sanity check of the determinant formula (8 s)
python complexes/scripts/certify_examples.py           # N1 and two further 6-complex examples, exact (< 1 s)
python complexes/scripts/eigen_enclosures.py           # eigenvalue enclosures, arb (< 1 s)
python examples/scripts/certify_example.py             # G7, exact (< 1 s)
python examples/scripts/check_gale_formula.py          # determinant formula of the 5-reactant reduction (< 1 s)
python examples/scripts/e2_reactant_sets.py 8          # 4-reactant structures (8 s)
python examples/scripts/enum5_reactions.py 12          # all 3 188 646 five-reaction networks (2 min, 12 cores)
python verification/certify_examples_own.py            # verifier: N1, N2, G7 from the reaction lists (< 1 s)
python verification/global_step_example.py             # verifier: why the global ACR step is needed (< 2 s)
python verification/dim2_relaxed_lp.py                 # verifier: LP search against the dimension-two proof (14 min)
python verification/dim2_network_scan.py 12            # verifier: dimension-two networks on 3 species (19 min)
python verification/item3_scope_example.py             # verifier: scope of the 5-reaction enumeration (< 1 s)
```

Each script prints PASS/FAIL per claim and exits non-zero on failure.

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: dimension-two theorem with full proof, consequences, certified examples |
| `scripts/` | Sanity check of the dimension-two determinant formula |
| `complexes/report.md`, `complexes/scripts/` | Minimal number of complexes: lemma, structure of 6-complex examples, certified N1–N3 |
| `examples/report.md`, `examples/scripts/` | Certified G7, minimal number of reactions on 3 species, 4 reactant complexes |
| `references.md` | Sources and novelty check |
| `NEEDS.md` | Sources obtained and not obtained |
| `verification.md`, `verification/` | Verification report and the verifier's own scripts and logs |
