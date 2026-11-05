# Becker–Döring bubblelator: a rigorous Hopf bifurcation in the limit model

<!-- plain-words:start -->
**The question.** Monomers are injected at a constant rate, clusters grow by absorbing them, and large clusters
are removed (a "bubblelator", or gas evolution oscillator). Niethammer, Pego, Schlichting and Velázquez reduced this
Becker–Döring model to a limit model and found numerically that it starts to oscillate (a Hopf bifurcation) when the
removal rate drops below a critical value. They wrote that a rigorous Hopf theorem had not yet been proven for a model
of this type, which has an infinitely long memory.

**The answer.** It does oscillate. Below an explicit critical removal rate, a branch of small periodic
solutions is born from the steady state. It is unique, its direction is the one the authors predicted, and at
that rate the steady state first becomes linearly unstable. This is proved for size-independent removal at every
cluster-growth exponent, and for size-dependent removal at every point of the grid on which the authors computed
(1111 parameter pairs).

**How.** The periodic solutions are found directly in a space of periodic functions, which avoids the infinite
memory (Lyapunov–Schmidt reduction). The three numerical conditions this needs, plus a fourth for the stability
change, become statements about one explicit integral. These are proved by hand when removal does not depend on
size, and in certified interval arithmetic otherwise.
<!-- plain-words:end -->

**Status: solved (preliminary)**, not yet independently verified.

**Attribution.** Solver: Claude Opus 5.5 (effort xhigh).

## Problem and source

B. Niethammer, R. L. Pego, A. Schlichting, J. J. L. Velázquez, *Oscillations in a Becker–Döring model with injection
and depletion*, SIAM J. Appl. Math. 82(4) (2022) 1194–1219, [doi:10.1137/20M1398664](https://doi.org/10.1137/20M1398664),
[arXiv:2102.06751](https://arxiv.org/abs/2102.06751) (v2). The limit model (1.12)–(1.14) reduces, for solutions defined
for all times, to the equation (3.7) with infinite delay

  u'(τ) = 1 − ∫₀^∞ e^{u(τ−z)} exp(−η z^{β+1}/(β+1)) z^ν dz,   ν = α/(1−α), β = r/(1−α),

where η is the rescaled removal rate. §1.4 (p. 4):

> "To our knowledge, a mathematically rigorous Hopf bifurcation theorem has been proven only for finite delay horizons."

§3 (p. 13):

> "To our knowledge, a rigorous Hopf bifurcation theorem has not yet been proven for a model of this type, though it
> seems plausible that one might extend existing methods for retarded functional differential equations with finite
> delays (RFDE) which are based upon rescaling to fix the temporal period and Lyapunov-Schmidt reduction [18, 12]."

The source argues formally. It proves transversality only for β = 0 (§4.1). It checks transversality (3.15) and
nonresonance (3.28) numerically on a grid ("verified numerically in the Appendix C", Tables 2–3). Its direction
coefficients (3.22)–(3.23) come from a formal expansion. The source's own wording is weaker than "open problem". General
Hopf theorems with fading memory exist (Staffans, JDE 1987); what was missing is a theorem for this model with its
hypotheses checked.

## Result

With G(λ) = ∫₀^∞ exp(−z^{β+1}/(β+1) − λz) z^ν dz and H(t) = Re G(it), the conditions at a frequency t₀ > 0 are:
(H1) H(t₀) = 0 and Im G(it₀) < 0; (H2) H′(t₀) ≠ 0; (H3) ϑ₀·ikt₀ + G(ikt₀) ≠ 0 for every k ≥ 2, where ϑ₀ = −Im G(it₀)/t₀;
and (H4): every other zero t of H has −Im G(it) < ϑ₀t.

* **Theorem A (Hopf theorem).** Under (H1)–(H3), at the explicit critical rate η₀ = (ϑ₀/G(0))^{β+1} and frequency
  κ₀ = t₀ϑ₀/G(0):
  - ±iκ₀ are simple characteristic roots crossing transversally.
  - A unique C^∞ branch of non-constant periodic solutions u_ε = u₀ + ε cos(κ(ε)τ) + O(ε²) bifurcates. It is unique up
    to time translation among solutions C¹-close to u₀ with period close to 2π/κ₀.
  - η(ε) = η₀(1 + η₂ε² + O(ε⁴)) and κ(ε) = κ₀(1 + κ₂ε² + O(ε⁴)), with the source's formulas (3.22)–(3.23), now proved.
  - Under (H4) and H′(t₀) < 0, all characteristic roots lie in Re λ < 0 for η > η₀, and exactly two lie in Re λ > 0
    just below η₀.
* **Theorem B (β = 0, removal independent of size).** For every ν > 0 the Hopf points are
  η_j = cos^{ν+2}ω_j / sin ω_j, with ω_j = π(1+4j)/(2(ν+1)) and 0 ≤ j < ν/4. All of (H1)–(H3) hold there, and (H4)
  holds at j = 0. The direction coefficient is η₂ < 0 for every ν > 0, with limits −1/36 as ν → 0 and
  1/(10π) − 3/20 as ν → ∞.
* **Theorem C (β ≥ 0, certified grid).** The grid is ν ∈ {1/9, 1/4, 3/7, 1/2, 2/3, 1, 3/2, 2, 7/3, 4, 9}, that is
  α ∈ {0.1, …, 0.9} together with α = 1/3 and α = 2/3, and β ∈ {0, 0.1, …, 10}. At all 1111 points (H1), (H2) with
  H′(t₀) < 0, (H4) and η₂ < 0 are certified in Arb ball arithmetic.
* **Theorem D (ν = 0, β = 3, the source's Table 1).** H has exactly 17 zeros in (0, 27.77]. The first eleven
  crossings are certified, with alternating crossing directions and η₂ < 0 at each. The number of unstable roots
  alternates between 0 and 2 on successive intervals of η.

In the source's sense the bifurcation is **supercritical** (η₂ < 0) everywhere it was checked. The source's
floating-point values agree with our enclosures to 10⁻⁶ relative. This covers η₀, κ₀, η₂, κ₂, Re ĝ₁(1) and |L̂(2)| in
its data files for Figures 3–4 and Tables 2–3 (909 points), and everything in its Table 1 except the κ₂ column (remark 2).

Two remarks on the source:

1. Its change of variables (3.1) drops the transport speed 1 − α and the factor 1/(1 − α) in the boundary value. Its
   equation (3.7), on which all its analysis rests, is exact for η̃ = (1 − α)^β η and u + ν log(1 − α). In the
   variables of (1.12)–(1.14) the critical rate is therefore η₀/(1 − α)^β.
2. The κ₂ column of its Table 1 equals η₂/4 in every row, which contradicts its own formula (3.23). Certified value
   at the first crossing: κ₂ = −0.10842900, against −0.10887 in the table. A Galerkin continuation agrees with the
   certified value.

**Scope.** The theorem is about the limit model, which is where the source observes the oscillations. It is not about
the infinite Becker–Döring system with injection and depletion, whose reduction to the limit model is formal in the
source. Orbital stability of the bifurcating orbits is not proved. Nor is nonlinear stability of the steady state:
part (5) of Theorem A is about characteristic roots. For β > 0 the theorem holds at the 1111 grid points, not
in between. The grid contains every β-value of the source's figures and tables, but continuity between the points is
expected, not proved. A centred-form certificate over β-boxes was tried: it works for small ν, but for large ν the ball
overestimation needs boxes of width about 2·10⁻⁴, so it was abandoned.

The theorems are numbered 5–8 in the paper; they are Theorems A–D here.

## Reproduce

Python 3 with python-flint 0.9 (Arb); NumPy and SciPy only for script 06. The certificates were run under WSL Fedora
(`search_grid.log`, `grid_certify.log`, `odd_beta3.log`) and again under Windows with Python 3.12 and python-flint 0.9.0
(`grid_certify_windows.log`, `beta0_direction.log`, `odd_beta3_windows.log`). From `scripts/`:

```bash
python 01_beta0_direction.py          # Theorem B, direction for all nu > 0          (~1 min)
python 02_search_grid.py              # float search of the t-intervals (exploration) (~15 min, 8 cores)
python 03_certify_grid.py             # Theorem C, 1111 points                        (~3-20 min, 8 cores)
python 04_odd_beta3.py                # Theorem D                                      (~2 min)
python 05_compare_source.py           # comparison with the source's data files       (seconds)
python 06_galerkin_check.py           # floating-point check of eta_2, kappa_2         (~1 min)
```

The manuscript is compiled from the repository root with
`python templates/build_tex.py problems/becker-doring-bubblelator-hopf/paper.tex`.

Every certificate script prints PASS/FAIL per claim and exits non-zero on failure. Logs are in `certificates/`.

## Files

- [paper.tex](paper.tex), [paper.pdf](paper.pdf): manuscript with complete proofs; [references.bib](references.bib).
- `scripts/bdhopf.py`: rigorous enclosures of G and of the Hopf quantities. It also evaluates H on the saddle line
  for β = 3.
- `scripts/certify_lib.py`: the checks C1–C5 at a parameter point.
- `scripts/01`–`06`: as above.
- `certificates/grid_tintervals.csv`: the t-intervals that `03` checks.
- `certificates/grid_results.csv`: enclosures of η₀, κ₀, ϑ₀, η₂, κ₂, Re ĝ₁(1) and |L̂(2)| at all 1111 points.
- `certificates/table_grid.tex`: the excerpt printed in the paper.
- `certificates/odd_beta3.csv`: data for Theorem D.
- `certificates/*.log`: the run logs.
- [references.md](references.md): sources and the novelty check.
