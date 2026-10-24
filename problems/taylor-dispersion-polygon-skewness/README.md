# Sign of the skewness of Taylor dispersion in regular polygonal ducts

**Solver: Claude Opus 5.5 (effort not recorded)**  
**Status: solved (preliminary)** for the long-time and short-time signs (computer-assisted proofs, not yet
independently verified); **partial result** for the "all times" claim.

## Original statement

M. Aminian, R. Camassa, R. M. McLaughlin, *Mass distribution and skewness for passive scalar transport in pipes with
polygonal and smooth cross sections*, Stud. Appl. Math. **141**(3) (2018) 399–417,
[doi:10.1111/sapm.12230](https://doi.org/10.1111/sapm.12230); [arXiv:1708.05097](https://arxiv.org/abs/1708.05097)
(v1–v3). A tracer released uniformly across a straight duct in pressure-driven flow; the sign of the skewness of its
cross-sectionally averaged distribution says whether it is front- or back-loaded.

arXiv abstract (v1, v2, v3):

> Interestingly, the equilateral triangle is the only regular polygon with this property, all other polygons possess
> positive skewness at all times, although this cannot cannot be proved on finite times due to the lack of closed form
> flow solutions for such geometries.

Published abstract (p. 399), weakened:

> The equilateral triangle appears to be the only regular polygon with this property—all others possess positive
> skewness at all times.

Introduction (p. 401) and Section 4 (p. 406), where the long-time claim is based on Mathematica finite elements for
n = 3, 4, 5, 6, 8, 20 (arXiv also n = 7, 100):

> We apply this analysis, using finite element solvers, and establish remarkably that the only polygonal domain
> yielding negative long-time solute negative skewness is the triangle.

> We use a finite element package to calculate the asymptotic coefficients, and show that in fact the triangle is the
> only one of the regular polygons for which this occurs.

("This property" = positive short-time and negative long-time skewness, proved for the triangle in the paper except
the long-time value, which is a finite-element number, eq. (18): −0.07076.)

## Result

Let `u` be the mean-zero Poiseuille profile (`−Δũ = 2`, `ũ = 0` on the wall, `u = ũ − ⟨ũ⟩`) and `g1` the mean-zero
Neumann cell problem `−Δg1 = u`. We re-derive the paper's long-time formula
`Sk(τ) ~ 3Pe³⟨u g2⟩/(2 + 2Pe²⟨u g1⟩)^{3/2} τ^{−1/2}` (eq. (17)) and prove `⟨u g2⟩ = 2⟨u g1²⟩`, so for every Péclet
number the long-time sign is the sign of `∫ u g1²`.

* **n = 3 (exact).** `g1` is a polynomial of degree 6; `⟨u g1⟩ = 9/1540`, `⟨u g1²⟩ = −7101/476476000 < 0`, and the
  long-time coefficient is `3⟨u g2⟩/(2⟨u g1⟩)^{3/2} = −789√770/309400 = −0.07076…` (inradius 1), confirming the paper's
  number exactly. Geometric skewness `4√21/99 > 0`.
* **Every n ≥ 4: long-time skewness positive** (for every Pe). Computer-assisted: a-posteriori error bounds
  (maximum principle, an L² bound for harmonic functions through the torsion function, Payne–Weinberger and a trace
  identity for the Neumann problem) evaluated in ball arithmetic for 4 ≤ n ≤ 40, and a perturbation argument around
  the disk, uniform in n, for all n ≥ 30 (all quantities are explicit functions of `t = tan(π/n)`).
* **Every n ≥ 3: short-time (geometric) skewness `⟨u³⟩/⟨u²⟩^{3/2}` positive.** Computer-assisted for 4 ≤ n ≤ 150
  (Green identity plus a Rellich bound for the first-order error) and uniformly for n ≥ 140, through the exact
  expansion `⟨u³⟩ = t⁴/360 + O(t⁵)` (explicit two-sided bounds; the main part is exactly `t⁴/360 + 19t⁶/7560`); hence `S^G_n ~ (8√3/15) tan⁴(π/n) ≈ 90/n⁴`.
* **Finite times (partial).** For symmetric cross-sectionally uniform releases, the sign of `Sk(τ)` does not depend on
  Pe. Certified windows: n = 3: `Sk > 0` for `τ ≤ 3.3·10⁻³`, `Sk < 0` for `τ ≥ 371`; n = 4: `Sk > 0` for
  `τ ≤ 5.3·10⁻⁴` and for `τ ≥ 32`; similar for n = 5, 6, 8 (time unit: apothem²/κ). The sign at intermediate times
  (the claim "at all times") remains open.
* Rigorous enclosures of Table 1 of the paper (n = 4, 5, 6, 7, 8, 20, 100) agree with all printed values except
  the arXiv-only `S^G_100 = 8.6·10⁻⁷`; the certified value lies in `[8.348, 8.397]·10⁻⁷`.

Full statements and proofs: [result.md](result.md). Novelty check: [references.md](references.md) (no earlier proof
found; the only later work on regular polygonal ducts, arXiv:2605.22982, concerns the Taylor coefficient of rods).

## Reproduction

From `scripts/` (Python 3.13, python-flint 0.9.0, sympy 1.14, numpy 2.x; one core):

```sh
sh run_all.sh
```

| script | claim | runtime |
|---|---|---|
| `triangle_exact.py` | n = 3: exact `g1`, `⟨u g1²⟩ < 0`, coefficient −789√770/309400, `S^G = 4√21/99` | 5 s |
| `lt_certify.py 4 40` | `∫u g1² > 0` for 4 ≤ n ≤ 40 | 1 min |
| `lt_uniform.py 30` | `∫u g1² > 0` for all n ≥ 30 | 2 s |
| `sg_certify.py 4 150` | `⟨u³⟩ > 0` for 4 ≤ n ≤ 150 | 15 min |
| `sg_uniform.py 140 6` | `⟨u³⟩ > 0` for all n ≥ 140 | 5 s |
| `alltime_windows.py` | certified time windows (n = 3, 4, 5, 6, 8) | 1 min |
| `table1.py` | enclosures of Table 1 of the paper | 5 min |
| `numerics_check.py` | non-rigorous sanity checks (not a certificate) | 30 s |

Each certificate script prints PASS/FAIL per claim and exits non-zero on failure. Floating point is used only to
choose approximation coefficients, which are then exact binary numbers recorded in the certificates.

## Files

* `result.md` — derivations, lemmas, theorems, certificate method, limitations.
* `references.md` — sources (what was read) and novelty check. `NEEDS.md` — inaccessible sources.
* `scripts/polylib.py` — exact/ball calculus of D_n-invariant polynomials on the n-gon (homogeneous components on the
  side `x = 1`). Other scripts as above.
* `certificates/lt_certificates_4_40.json`, `certificates/sg_certificates_4_150.json` — coefficients (hex floats),
  error terms and bounds per n; `certificates/*.log` — run logs.
* `certificates/sources/` — third-party PDFs (git-ignored).
