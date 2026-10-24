# Diblock foam model screen of the Frank–Kasper catalogue

**Status: partial result — first exhaustive DFM screen of the Frank–Kasper catalogue (numerical, with exact upper-bound
certificates); not independently verified.** σ remains the minimiser of the diblock foam model over all 121 catalogue
structures (84 of Dutour Sikirić–Delgado-Friedrichs–Deza 2010 + 37 of part II, 2012), numerically, with a gap of
0.022 %. Reproduction of Reddy et al. (2018): five values within 1.6×10⁻⁶; five published values are *above* our
exactly certified upper bounds.

**Attribution.** Method, code, runs, certificates and write-up: Claude Opus 5.5 (effort not recorded), solver
sub-agent of the orchestrator, 24 Oct 2026. Verification: pending.

## Problem

Reddy, Buckley, Arora, Bates, Dorfman, Grason, *Stable Frank–Kasper phases of self-assembled, soft matter spheres*,
PNAS 115 (2018) 10233 (arXiv:1806.01726) introduce the **diblock foam model (DFM)**: minimise
𝓕 = (𝒜² ℐ)^{1/3} over curved periodic partitions with a fixed cell topology and free cell volumes (𝒜: mean cell area,
ℐ: mean second moment of each cell about its centroid, both normalised by the sphere of the mean cell volume). For 11
Frank–Kasper (FK) phases and BCC they find:

> "Among these nearly degenerate, fully unconstrained DFM structures, the σ phase overtakes A15 (minimal for fixed,
> equal volume) as the minimal energy phase (with next lowest energy for P)" — [1], Sec. II (arXiv v1, p. 7).

The catalogue of FK structures of M. Dutour Sikirić, O. Delgado-Friedrichs, M. Deza (*Space fullerenes: a computer
search for new Frank–Kasper structures*, Acta Cryst. A 66 (2010) 602; 84 structures, 9 of them in [1]) and its part II
(Dutour Sikirić & Deza, Struct. Chem. 23 (2012) 1103; 37 structures) had never been evaluated in the DFM.
**No source states this as an open problem**: "is σ the DFM minimiser over the whole catalogue?" is our framing.
Novelty check (`references.md`): no prior DFM evaluation of catalogue structures found.

## Result in one paragraph

We wrote a curved-cell DFM minimiser in which every candidate is an explicit triangulated partition whose 𝒜, cell
volumes and second moments are evaluated **exactly** (polyhedral formulas), so every computed value is an upper bound
on the infimum for its topology; meshes are refined to 2.6×10⁵–6.7×10⁵ triangles and extrapolated (observed O(h²)
convergence, ratio 3.93). It reproduces [1] for A15, C15, Z, C14 and δ within 1.6×10⁻⁶ (A15 and C15, cubic, within
4.3×10⁻⁷), the Kelvin (1.097251) and Weaire–Phelan (1.093541) foams, and the paper's flat Voronoi values; for BCC, P,
M, H and μ exact integer-arithmetic certificates show that the published values are not minimal (gaps 3×10⁻⁶ to
1.2×10⁻⁴). For the screen we relaxed the box shape of every structure (the paper keeps metallic aspect ratios fixed,
which penalises e.g. H by 3.9×10⁻⁴ and μ by ~7×10⁻⁴). Result: **σ is the lowest of all 124 structures** (121 catalogue
+ P, M, BCC), 𝓕(σ) = 1.0685748 (certified upper bound 1.068575909), followed by a tie between the new catalogue entry
**E11** (I4₁/amd) and **H**, both at +0.0219 % (2.34×10⁻⁴, about 2000 times our extrapolation uncertainty), then μ and
zra-d (K₇Cs₆ type, new to the DFM) at +0.025 %, the new **N56** at +0.027 % and P at +0.028 %. Sixteen catalogue entries
cannot be realised as Voronoi diagrams; they are geometrically frustrated and relax to values ≥ 1.5 % above σ, most
with a vanishing foam edge. The ranking is numerical (no lower bounds): see the caveats in `result.md` §5.

## Main table (box-relaxed DFM, all structures)

`data/dfm_screen_boxrelaxed.csv` has every structure with all level values. 𝓕_L is the exact energy of the relaxed
level-L partition (an upper bound for that topology in the relaxed box); "extrapolated" is the Richardson value
(heuristic error). Cubic A15, C15, BCC: see also the fixed-box level-4 values (identical boxes).

<!-- TABLE -->

## Milestone 1 — reproduction of [1] (`scripts/check_milestone1.py`: 0 FAIL)

Fixed paper boxes, levels 0–4 (up to 4.9×10⁵ triangles), Richardson-extrapolated (error ≤ 1×10⁻⁷ except BCC 7×10⁻⁷):

| structure | ours (extrapolated) | published [1] | difference | status |
|---|---|---|---|---|
| A15 | 1.0690374 | 1.069037 | +3.8×10⁻⁷ | reproduced |
| C15 | 1.0692434 | 1.069243 | +4.3×10⁻⁷ | reproduced |
| Z | 1.0692501 | 1.069250 | +9.5×10⁻⁸ | reproduced |
| C14 | 1.0690704 | 1.069071 | −5.6×10⁻⁷ | reproduced |
| δ | 1.0696755 | 1.069677 | −1.5×10⁻⁶ | reproduced |
| pσ | 1.0690656 | 1.069068 | −2.4×10⁻⁶ | lower |
| σ (a/c = 1.9) | 1.0685789 | 1.068583 | −4.1×10⁻⁶ | lower; box ambiguous (paper's CVP row fits a/c ≈ 1.93) |
| BCC | 1.0711023 | 1.071111 | −8.7×10⁻⁶ | certified: published value not minimal |
| H (Z12 at 2a+4f) | 1.0692035 | 1.069219 | −1.6×10⁻⁵ | certified: not minimal |
| P | 1.0688810 | 1.068909 | −2.8×10⁻⁵ | certified: not minimal |
| M | 1.0689983 | 1.069038 | −4.0×10⁻⁵ | certified: not minimal |
| μ | 1.0695330 | 1.069662 | −1.3×10⁻⁴ | certified: not minimal |

Also: Kelvin foam (equal-volume minimal area, BCC) 1.0972512 vs 1.097251; Weaire–Phelan 1.0935413 vs 1.093541;
A15 equal-volume DFM 1.069742 vs 1.069743; flat "Initial Voronoi" (𝒜, ℐ) rows of [1] reproduced to 6 decimals for
BCC, A15, C15, H, P, δ, pσ, M. Details and diagnostics: `result.md` §3.

## How to reproduce

Environment: Python 3.12+ with numpy, scipy, numba, gemmi (`pip install gemmi`); one thread (enforced by
`scripts/dfm/__init__.py`). Catalogue data: Zenodo 10.5281/zenodo.13691342 (`SpaceFullerene.tar.gz`, md5
582371cee665c62e01a2f47a73dd2c9d) unpacked into `certificates/sources/dutour2010/` (git-ignored).

```bash
cd problems/frank-kasper-dfm-screening/scripts
python3 check_certificates.py        # 6 exact certificates, ~1 min: PASS x6
python3 check_catalogue.py           # 121 catalogue files: valid embedded FK tilings, ~5 min: 121 PASS
python3 check_milestone1.py          # Milestone 1 (uses data/dfm_fixedbox_paper_structures.csv), ~5 min: 0 FAIL
# regenerate the runs (hours on a loaded machine; resumable, results in ../work/runs)
python3 run_dfm.py BCC A15 C15 Z C14 sigma psigma H mu M P delta --box fixed --maxlev 4
python3 screen_pipeline.py ../data/screen_order.txt --win3 5e-4 --win4 3e-4
python3 analyze.py --box fixed --out ../data/dfm_fixedbox_paper_structures.csv
python3 make_tables.py
```

## Files

| path | content |
|---|---|
| `result.md` | method, validation (Milestone 1), screen, caveats |
| `references.md` | sources (what was read) and novelty check |
| `NEEDS.md` | paywalled sources requested |
| `scripts/dfm/` | engine: `topology.py` (periodic Delaunay, foam), `mesh.py` (periodic triangulated foam, refinement), `energy.py` (exact geometry, 𝓕, gradient, virial), `relax.py` (restricted motion, smoothing, box strain, embeddedness check), `lbfgs.py`, `constrained.py` (equal-volume problems), `cvp.py` (Lloyd), `cgd.py` (catalogue reader), `embed.py` (spring re-embedding), `registry.py` (structures), `structures.py` |
| `scripts/run_dfm.py` | resumable relaxation driver (levels, checkpoints, JSON results) |
| `scripts/screen_pipeline.py` | adaptive screen (L2 all, L3/L4 near σ) |
| `scripts/analyze.py`, `scripts/make_tables.py` | extrapolation and tables |
| `scripts/certify_upper_bound.py`, `scripts/check_certificates.py`, `scripts/make_certificates.py` | exact certificates (integer arithmetic) |
| `scripts/check_milestone1.py`, `scripts/check_catalogue.py` | PASS/FAIL checks |
| `scripts/diagnostics/se_quadrature_bias.py` | one-point quadrature diagnostic for Eq. (B2) of [1] |
| `data/published_reddy2018.csv` | tables of Figs. 7–18 of [1] |
| `data/dfm_fixedbox_paper_structures.csv` | fixed-box runs (levels 0–4) of the 12 structures of [1] |
| `data/dfm_screen_boxrelaxed.csv` | screen: every structure, all levels, extrapolation |
| `data/catalogue_checks.csv` | per-file topology checks of the catalogue |
| `data/screen_order.txt` | the list of runs of the screen |
| `certificates/meshes/*.npz`, `certificates/claims.json` | integer partitions (2⁻⁴⁴ grid) certifying upper bounds |
| `certificates/logs/` | logs of the checks |
| `PROGRESS.md` | work log (git-ignored) |
