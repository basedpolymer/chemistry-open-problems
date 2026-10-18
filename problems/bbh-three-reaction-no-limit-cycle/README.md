# Banaji–Boros–Hofbauer: does the three-reaction example have a limit cycle?

**Status: solved (preliminary): complete proof, awaiting independent verification.** The broader question
(all case-9 and case-10 networks without Andronov–Hopf bifurcation) remains open beyond the subfamily below.

**Attribution.** Dulac argument (scout) and subfamily (verifier): Claude Opus 5.5 exploration agents (effort high), 17 Oct 2026. Certificate script and manuscript: Claude Opus 5.5 (effort xhigh). Independent verification: pending.

## Problem

M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*, Stud. Appl. Math.
152 (2024) 249–278, doi:10.1111/sapm.12639, arXiv:2304.02303. Conclusions (§6):

> "we leave it open whether the networks in cases 9 and 10 admitting no Andronov–Hopf bifurcations can have a limit
> cycle. For example, we do not know whether the system 2X → 3X + 2Y, X + Y → 0, Y → X + Y,
> ẋ = κ₁x² − κ₂xy + κ₃y, ẏ = 2κ₁x² − κ₂xy admits a limit cycle."

## Result

**No: the example has no periodic orbit for any κ > 0.** The Dulac function y^(−1−2κ₁/κ₂) turns the divergence into
a negative definite quadratic form (discriminant −8κ₁κ₂), and the axes are crossed inwards.

The same argument classifies the whole subfamily 2X → (2+c₁)X + d₁Y, X + Y → 0, Y → (1+c₃)X + Y (c₁, c₃, d₁ ≥ 1):
periodic orbits can only exist for c₁ < d₁ < 2c₁, which is exactly the Andronov–Hopf region of the authors'
Theorem 9. In this subfamily, "no Andronov–Hopf bifurcation" implies "no limit cycle".

## Reproduce

```bash
python scripts/certify_dulac.py   # symbolic identities (sympy), < 10 s
```

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: Theorems 1–2 and their proofs |
| `scripts/certify_dulac.py` | Symbolic verification of every identity |
| `certificates/certify_dulac.log` | Output, all PASS |
| `references.md` | Sources and novelty check |

## History

Found during the exploration phase (17 Oct 2026, scout and verifier); subfamily classification and write-up on
18 Oct 2026.
