# The Kaibel–Petlyuk vapour penalty

<!-- plain-words:start -->
**The question.** The Kaibel column separates four components in one dividing-wall shell. Does it always need strictly more vapour than the full Petlyuk arrangement?

**The answer.** Yes, in the ideal infinite-stage model, but this was already shown by Halvorsen and Skogestad in 2006. The folder is an audit: an explicit proof, exact examples, and a misprint (“min” for “max”) corrected in the 2006 paper.

**How.** Underwood’s method; the Kaibel requirement is the maximum of two binary requirements.

**Overlaps.** [Result 12](../halvorsen-skogestad-petlyuk-minimum-vapor/) already implies V(Kaibel) ≥ V(Petlyuk); this folder adds the strict inequality and a formula, both known since 2006.

<p align="center"><img src="../../figures/petlyuk-vmin-diagram.svg" alt="Minimum-vapour diagram with two peaks" width="720"></p>

<p align="center"><sub>The Petlyuk pays only the highest peak of this diagram; a Kaibel column must meet two binary separations at once.</sub></p>

More detail, with every step explained: [guided reading, result 26](https://basedpolymer.eu/chemistry-open-problems#result-26) · [en français](https://basedpolymer.eu/chemistry-open-problems/fr#d26).

<!-- plain-words:end -->

**Status: already solved by Halvorsen and Skogestad (2006); candidate closed by literature audit, 27 October 2026.**

**Attribution.** Original result: Ivar J. Halvorsen and Sigurd Skogestad.
Literature audit, explicit proof exposition, certificates and editing: **GPT-6-Astra (effort High)**.
No separate fresh-context verifier was used.

The candidate asked whether the four-product Kaibel arrangement always needs strictly more energy than
the fully extended Petlyuk arrangement, and requested a closed form. In the ideal infinite-stage
Underwood model, **yes**: with four positive feed fractions and distinct positive relative volatilities,
the minimum vapour-flow difference is strictly positive. Both the comparison and an algebraic
calculation were already given in 2006. This is a documented prior answer, not a new research result.

## Primary source and correction

I. J. Halvorsen and S. Skogestad, *Minimum Energy for the Four-Product Kaibel Column*,
AIChE Annual Meeting 2006, paper 216d:
[full paper](https://skoge.folk.ntnu.no/publications/2006/halvorsen_kaibel_aiche06/paper216d.pdf),
[conference record](https://proceedings.aiche.org/conferences/aiche-annual-meeting/2006/proceeding/paper/216d-minimum-energy-four-product-kaibel-distillation-column).

The full sentence settling the question is quoted verbatim in Section 1 of [the note](paper.pdf),
from Section 2, PDF page 4 of the source. There is no source declaring this candidate open:
the repository had explicitly left its literature check pending.

**Equation (7) of the source prints `min` twice where `max` is required.** Both binary sections
must meet their vapour requirements; their maximum is binding. This agrees with the source's
Sections 2 and 4 and its numerical Table 1. The original PDF page was inspected visually.
The certificates reproduce the published rounded values 1.38 (Petlyuk) and 1.83 (Kaibel).

## Complete answer in the stated model

Let the relative volatilities be `a > b > c > d > 0`, the four positive mole fractions be
`z_A, z_B, z_C, z_D` with sum one, and `F > 0` the feed rate. Set `s = 1 - q`, where `q`
is the feed thermal parameter. The algebra works for every finite real `q`, including the usual
two-phase range `[0,1]`. Define

\[
h(x)=\frac{az_A}{a-x}+\frac{bz_B}{b-x},\qquad
g(x)=\frac{cz_C}{c-x}+\frac{dz_D}{d-x}.
\]

Let `theta_A`, `t`, `theta_C` be the unique roots of `h(x)+g(x)=s` in `(b,a)`, `(c,b)`,
`(d,c)`, respectively. Put `u=h(t)>0`, `v=u-s=-g(t)>0`. Then

\[
P=\max\left\{\frac{az_A}{a-\theta_A},\ u,\ s+\frac{dz_D}{\theta_C-d}\right\},
\qquad
K=\max\left\{\frac{u(a-t)}{a-b},\ s+\frac{v(t-d)}{c-d}\right\}.
\]

The minimum overhead vapour flows are `F P` and `F K`, respectively, and their gap is
**`F (K-P) > 0`**. The bottom-flow gap is identical. The note proves all strict comparisons,
including ties between the three Petlyuk requirements, and derives the simplified Kaibel formula
by eliminating the propagated binary roots. All required roots solve an explicit polynomial of
degree at most four; the formula is algebraic, not asserted to be rational in the feed data.

For the published example `F=1`, `alpha=(6,4,2,1)`, `z=(0.3,0.2,0.2,0.3)`, `q=1`, the
certified intervals imply:

| Quantity | Strict enclosing interval |
|---|---|
| Petlyuk minimum | `(1.37983442383, 1.37983442384)` |
| Kaibel minimum | `(1.82607549321, 1.82607549322)` |
| Difference | `(0.44624106937, 0.44624106939)` |

These are normalized vapour flows. Equating their comparison with heat-duty comparison requires
the same latent-heat convention. Constant relative volatilities, constant pressure and constant
molar overflow, pure products and the infinite-stage limit are part of the model. No claim is made
for finite-stage attainment, absent components, coincident volatilities, arbitrary nonideal mixtures
or arbitrary product purities. No uniform positive gap over all feeds is claimed.

## Reproduce

Python 3.10 or later, standard library only; the checks take less than one second in the tested environment.
From the repository root:

```bash
python problems/kaibel-petlyuk-energy-gap/scripts/verify.py
python problems/kaibel-petlyuk-energy-gap/scripts/check_rejections.py
python templates/build_tex.py problems/kaibel-petlyuk-energy-gap/paper.tex
```

On Windows, `py -3.12` can replace `python`. The PDF build requires the repository's LaTeX toolchain.
To regenerate the exact data before checking it:

```bash
python problems/kaibel-petlyuk-energy-gap/scripts/generate.py
```

The checker uses rational endpoint signs to certify each unique root, then exact interval arithmetic
to certify the flows. Its seven examples cover both bottleneck sections, `q=0`, `q=1`, interior
values, `q=-1/2`, `q=3/2`, and a rational middle root. These are computational checks of examples;
the general result is proved analytically in the note.

## Files

- [paper.tex](paper.tex), [references.bib](references.bib), [paper.pdf](paper.pdf): attributed exposition and full comparison proof.
- [scripts/verify.py](scripts/verify.py): exact checker, with `PASS`/`FAIL` output and nonzero failure exit.
- [scripts/generate.py](scripts/generate.py): deterministic optional certificate generator.
- [scripts/check_rejections.py](scripts/check_rejections.py): rejection checks on deliberately corrupted certificates.
- [certificates/examples.json](certificates/examples.json): exact inputs, root intervals and output bounds.
- [certificates/verification.log](certificates/verification.log): successful check output.
- [certificates/rejections.log](certificates/rejections.log): successful rejection-check output.
- [references.md](references.md): primary-source audit and source hash.
- [verification.md](verification.md): checks and their limits; no separate verifier claim.

The source PDF and rendered source page are retained locally under ignored `certificates/sources/`.
