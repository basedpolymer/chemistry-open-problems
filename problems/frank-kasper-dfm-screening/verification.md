# Verification record — 24–25 October 2026

**Global verdict: INCONCLUSIVE.** The requested global ordering of topology
infima has not been proved or refuted. The numerical screen favours σ.

This is a **same-agent audit**, by Codex (GPT-6; effort not recorded), with
separate integration and geometric verification code. It is not the fresh,
independent review defined in the repository protocol. The owner explicitly
requested no sub-agents, and none were used during the continuation.

## Verified statements

| Claim | Verification | Outcome |
|---|---|---|
| Six stored partitions are embedded and have the declared topology | Spherical cell boundaries, disk faces, links and incidences; periodic dual identity; positive integer cone determinants and exact volume sum | 6 PASS |
| Their energy upper bounds hold in the intended exact boxes | Independent rational simplex integration, integer radical enclosures, Arb affine transport | 6 PASS |
| All catalogue contact complexes have embedded FK realizations | Standalone integer triangulations, oriented face pairing, spherical vertex links, degree one, FK coordination and isolated hexagons | 121 PASS |
| Original catalogue sanity checks reproduce | Restored source archive, original floating-point checker | 121 PASS, 0 FAIL |
| Original validation script passes | Recomputed flat inputs and equal-volume examples; checks of stored fixed-box data | 0 FAIL |
| Published SI supports the five source comparisons | Visual check of all twelve unconstrained DFM energies and caption ratios; source hash and transcription retained | All agree with preprint inputs |
| Invalid test cases are rejected | Missing tetrahedron, reversed triangle, duplicated face, open boundary, pinched link and wrong box; invariance checks | 6 tests, OK |
| One-parameter flat A15 expression obeys F(c) > 1.0699 for every c ≥ 0 | Positive Bernstein coefficients on [0,2], positive derivative for c ≥ 2; 192-bit outward arithmetic | PASS |
| Uniform F > 1.00000025 for all admissible periodic partitions | Moment defect, volume-weighted concentration and equal-ball packing; exact rational constants; invokes the established Kepler theorem | PASS for constants; analytical proof in manuscript; no fresh independent review |
| Uniform elementary F > 1.000000004 | Same argument with a fully reproduced Blichfeldt bound | PASS for exact rational constants |
| Global bounds over all affine deformations of six stored meshes | Exact full moment tensors, convex supporting planes, 192-bit positive-definiteness and endpoint checks | 6 PASS; repeated using saved metrics without optimization |
| Affine integration and certificate guards | Analytic cube minimum, off-diagonal tensor under shear, translation, invalid metric and false dual bound | 3 tests, OK |
| New symmetry-breaking starts give the same observed ordering | Two seeds for σ, E11, H and μ through levels 0–2 | 24 stages completed; numerical evidence only |

The mathematical justifications are in [paper.tex](paper.tex), compiled as
[paper.pdf](paper.pdf). The degree proof requires a manifold, not merely a
volume identity. The new checks establish the manifold hypotheses and the
proof covers points on the simplex skeleton as well.

The six certified upper bounds are 1.071107810 (BCC), 1.068897166 (P),
1.069011740 (M), 1.069208084 (H), 1.069546495 (μ) and 1.068575909 (σ).
The first five use exact caption boxes; σ uses its stored relaxed dyadic
box. The first five also beat the lower rounding endpoints of the printed
arXiv-v1 values, also confirmed in the published SI. This does not assert an error in a proved source theorem.

## Mathematical scope and remaining gap

The analytic universal bound is now F > 1.00000025. Together with the σ witness it gives
**1.00000025 ≤ inf F ≤ 1.068575909** for the catalogue, with or without the P, M and
BCC comparators. No competing topology has a certified global lower bound
above the σ witness. The restricted flat A15 result does not cover arbitrary
curved A15 partitions.

The new uniform bound excludes sequences with F tending to 1 without imposing
a lower cell-volume cutoff or a box aspect-ratio cutoff. Its numerical gain is
tiny: 0.00000025 compared with the remaining gap of about 0.0686. It does not
exclude degenerations near σ. The geometric argument is written out in full;
the script checks its constants, not the argument in a formal proof assistant.
The Kepler theorem is used as an established external theorem, not reproved.

The affine certificates enclose the infimum over all nonsingular affine maps
of each stored mesh: BCC [1.071107807, 1.071107811], P [1.068892909, 1.068892956],
M [1.068962667, 1.068962682], H [1.068948492, 1.068948499],
μ [1.069099093, 1.069099122], σ [1.068575891, 1.068575910].
The first five lower endpoints exceed the σ witness. These are global bounds
on six affine families, not on topology classes. The lower inherited value
for box-relaxed H uses internal cell deformations and lies outside the affine
family of its saved fixed-box mesh. A rational metric and a positive-definite
dual slack are saved for each mesh. Reverification needs no optimizer.

The admissible class has curved faces, free positive cell volumes and
homogeneous box strains of each represented period. It does not implicitly
include all larger supercells or all FK structures. Fixed graph constraints
alone cannot strengthen the isolated-cell sphere bound. The A15 example
also shows finite energy as two cell volumes vanish; it does not exclude
other degenerating families near the observed minimum.

Neither local stationarity, a local Hessian, agreement of restarts, nor
Richardson extrapolation supplies the missing global lower bound. The
manuscript gives the first variation and a sufficient divergence-field
lower-certificate template, without claiming to have constructed the fields
needed to settle the comparison.

## Corrections to the inherited account

- Floating-point mesh energies were labelled “exact upper bounds”. The
  formulas are exact, but those decimal evaluations were not certified.
  The CSV column is now `F_maxlevel_mesh_energy`; its numeric rows are unchanged.
- Sixteen entries fail the Voronoi test at the supplied coordinates. This
  does not prove the nonexistence of another Voronoi realization.
- The catalogue has 105 such Voronoi-realised inputs; 108 includes P, M, BCC.
- Eight entries are represented only by unrelaxed starts after degenerate
  runs. Their high values do not exclude their topology classes.
- The five printed-number comparisons were initially limited to arXiv v1.
  The published SI was subsequently retrieved and visually checked; its
  values agree. The manual transcription and PDF hash are in
  `data/pnas_si_crosscheck.json`.
- Non-dyadic target boxes now have a rigorous affine transport argument.
- The original catalogue checker now fails when the required 84/37 inputs
  are absent. CGD rational parsing uses `Fraction` instead of `eval`.
- Three spring-embedding energies (part-II Y04–Y06) changed slightly on
  rerun; these floating-point diagnostics do not affect the energy certificates.

## Reproduction

Environment used: Windows host, WSL Fedora 44; Python 3.12.15 with NumPy,
SciPy, Numba, python-flint, SymPy and gemmi 0.7.5. From the repository root:

```sh
python problems/frank-kasper-dfm-screening/scripts/audit_certificates.py
python problems/frank-kasper-dfm-screening/scripts/certify_catalogue.py
python problems/frank-kasper-dfm-screening/scripts/test_audit_certificates.py
python problems/frank-kasper-dfm-screening/scripts/check_a15_family.py
python problems/frank-kasper-dfm-screening/scripts/check_uniform_gap.py
python problems/frank-kasper-dfm-screening/scripts/check_affine_orbits.py
python problems/frank-kasper-dfm-screening/scripts/test_affine_orbits.py
python problems/frank-kasper-dfm-screening/scripts/check_source_comparison.py
python problems/frank-kasper-dfm-screening/scripts/check_milestone1.py
python problems/frank-kasper-dfm-screening/scripts/check_restarts.py
```

The catalogue verification needs no external source files. Its regeneration
uses `--build` and the Zenodo archive. The energy identity checks use the
repository's crystallographic registry; energies are independently computed.
Logs are in `certificates/logs/`; machine-readable results are in
`audit_report.json`, `catalogue_exact_report.json`, `a15_family_report.json`,
`uniform_gap_report.json`, `affine_orbits_report.json`
and `data/symmetry_breaking_restarts.csv`.

Most historical optimization checkpoints and literature API snapshots were
absent from this checkout. The historical screen is preserved and labelled
as such; this audit does not claim to have recomputed all 124 high-resolution
relaxations or independently re-enumerated the catalogue.
