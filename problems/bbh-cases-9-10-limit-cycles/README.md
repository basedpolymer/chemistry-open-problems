# Banaji–Boros–Hofbauer, cases 9 and 10: can a network without Andronov–Hopf bifurcation have a limit cycle?

**Status: solved (preliminary).** The proof is complete but has not yet been checked by an independent verifier. The
answer is **no**.

**Attribution.**
- Solver (classification, Liénard argument, scripts, write-up): Claude Opus 5.5 (effort xhigh), 24 October 2026.
- Builds on result 2 of this repository (`problems/bbh-three-reaction-no-limit-cycle/`), which settled the example
  and one subfamily.
- Independent verification: pending.

## Problem

M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*, Stud. Appl. Math.
152 (2024) 249–278, doi:10.1111/sapm.12639, arXiv:2304.02303. Conclusions, §6, p. 22 of arXiv v1 (the published
version has the same text):

> "However, we leave it open whether the networks in cases 9 and 10 admitting no Andronov–Hopf bifurcations can have
> a limit cycle. For example, we do not know whether the system 2X → 3X + 2Y, X + Y → 0, Y → X + Y,
> ẋ = κ₁x² − κ₂xy + κ₃y, ẏ = 2κ₁x² − κ₂xy admits a limit cycle."

Cases 9 and 10 are the two-species, three-reaction networks with source complexes 2X, X+Y, Y and 2X, X+Y, 0. Target
complexes are arbitrary. Which of them admit an Andronov–Hopf bifurcation is decided by Theorem 9 of arXiv v1, which
is Theorem 3 of the published version.

## Result

**Theorem A.** Let a case-9 or case-10 network fail the Andronov–Hopf condition of BBH Theorem 9. Then for every
choice of rate constants its mass-action system has **no periodic orbit** in the closed positive quadrant, and in
particular no limit cycle.

With BBH's result that the Hopf bifurcation in these cases is supercritical, this gives a classification: in cases
9 and 10, a network has a limit cycle (equivalently, a periodic orbit) for some rate constants if and only if it
admits an Andronov–Hopf bifurcation.

**How the proof goes.**
1. Elementary arguments, or BBH's own Dulac function 1/(xy), dispose of every network except three explicit
   families, all with the reaction 2X → (2+c₁)X + d₁Y, c₁ ≥ 1:

   | family | other two reactions | condition |
   |---|---|---|
   | 9a | X+Y → 0, Y → c₃X + (1+d₃)Y | d₁/c₁ + d₃/c₃ ≥ 2 |
   | 9b | X+Y → Y, Y → c₃X | c₃d₁ ≥ c₁ |
   | 10a | X+Y → 0, 0 → c₃X + d₃Y | d₁/c₁ + d₃/c₃ ≥ 2 |

   The BBH example is in 9a, with equality in the condition.
2. These systems are linear in y. A shear along the x-nullcline and a time change turn them into a Liénard system
   x' = Y − F(x), Y' = −g(x). Its damping splits as F = (nullcline) − ρ·(increasing function). Here ρ is the only
   essential rate parameter, and G = ∫g is ρ times a ρ-independent function.
3. A comparison lemma, proved in full, rules out closed orbits whenever F(a) > F(b) at all points a < 1 < b with
   G(a) = G(b). Because of the splitting in step 2, this condition reduces to one inequality between the nullcline
   and G that does not involve ρ. That inequality is proved in closed form, H_k(k²/a) > H_k(a) for k ≥ 1, and is
   exactly where "no Hopf" (k ≥ 1) enters.

For case 10, the same computation also shows that **every** network with a positive non-saddle equilibrium,
including the Hopf ones, has no periodic orbit whenever tr J < 0 (Remark 2 of `result.md`).

Full proofs: [`result.md`](result.md).

## Reproduce

```bash
python3 scripts/certify_lienard.py      # exact sympy identities and sign claims (64 checks), < 10 s
python3 scripts/classify_networks.py    # exhaustive classification in the box of targets <= 20 (1.7e8 networks), ~1 min
python3 scripts/numerical_sanity.py     # floating-point sanity checks, not part of the proof, ~3 min
```

Each script prints PASS/FAIL per claim and exits with status 1 on failure. The logs are in `certificates/`. The
environment was Python 3.13, sympy 1.14, numpy 2.5, scipy 1.18, single thread.

## Files

| File | Content |
|---|---|
| `result.md` | Statement, complete proofs (classification §2, Liénard argument §3), remarks, verification map |
| `paper.tex`, `references.bib` | Draft manuscript written by the solver; not compiled, since there is no LaTeX in the container. The editor agent will finalise it |
| `scripts/certify_lienard.py` | Symbolic certificate: normalisation, Liénard identity, family data, key inequalities, Dulac classes |
| `scripts/classify_networks.py` | Exhaustive check of the classification (Lemmas 2.1–2.3) for all targets with coefficients ≤ 20 |
| `scripts/numerical_sanity.py` | Backward integration (no cycles), Hopf controls (cycles found), quadrature check of the comparison, sharpness observation |
| `certificates/*.log` | Outputs of the three scripts |
| `references.md` | Sources read and novelty check |
| `NEEDS.md` | Paywalled sources needed (none) |
