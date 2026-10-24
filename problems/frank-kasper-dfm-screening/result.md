# Diblock foam model on the Frank–Kasper catalogue: method, validation, screen

*Status: partial result (numerical screen with exact upper-bound certificates; not independently verified). Summary
table: README.md; all data: `data/`.*

## 1. Question and framing

Reddy, Buckley, Arora, Bates, Dorfman & Grason [1] introduced the **diblock foam model (DFM)**: for a periodic
partition X of space into n cells (per period) with a prescribed cell topology,

  𝓕(X) = (𝒜(X)² ℐ(X))^{1/3},  𝒜 = (n⁻¹ Σᵢ Aᵢ)/(4πR₀²),  ℐ = (n⁻¹ Σᵢ Iᵢ)/(4πR₀⁵/5),  (4π/3)R₀³ = n⁻¹ Σᵢ Vᵢ,

where Aᵢ is the area of cell i, Iᵢ = ∫_{Vᵢ}|x − xᵢ|² d³x its second moment about its centroid xᵢ ([1], Eqs. (1)–(2)
and Appendix A). Cells are curved, volumes free, topology fixed; [1] minimises 𝓕 with Surface Evolver at **fixed
periodicity** ([1], Appendix B: "for a given set of constraints, fixed periodicity, number of cells …"; the box
aspect ratios are given in the captions of Figs. 7–18). Among 11 FK phases (+ BCC) σ is the minimiser ([1], p. 7:
"Among these nearly degenerate, fully unconstrained DFM structures, the σ phase overtakes A15 (minimal for fixed,
equal volume) as the minimal energy phase (with next lowest energy for P)").

**No source poses an open problem here.** The question addressed is the natural extension: is σ still the DFM
minimiser over the catalogue of Frank–Kasper structures of Dutour Sikirić, Delgado-Friedrichs & Deza [2] (84
structures with ≤ 20 cells in a reduced fundamental domain; 9 of them are among the structures of [1]) and its
part-II extension [3] (37 structures with 21 cells)? Answering it numerically requires a curved-cell minimiser that
reproduces [1]; flat (power-diagram) cells cannot (flat optimum of A15 ≈ 1.0700 versus 1.069037 curved, [1] App. A2).

## 2. Method

**Discretisation.** The foam of a periodic triangulation (Delaunay complex of the cell centres, or the 4-cliques of
the catalogue's contact graph) is built with periodic shifts; every foam face (pentagon/hexagon) is fan-triangulated
from its centroid (level 0) and refined by 4-to-1 midpoint subdivision (level L has 4^L times more triangles; foam
edges are polylines shared by three faces, foam vertices are shared by four cells). For a closed triangulated cell
surface and any reference point r,

  Vᵢ = Σ det(a,b,c)/6,  Mᵢ = Σ Vₜ(a+b+c)/4,  Sᵢ = Σ Vₜ(|a|²+|b|²+|c|²+a·b+b·c+c·a)/10,  Iᵢ = Sᵢ − |Mᵢ|²/Vᵢ

(sums over the oriented boundary triangles, a, b, c relative to r) are **exact** for the polyhedral cell. Hence
every embedded triangulated partition is an admissible competitor and its exactly evaluated 𝓕 is an **upper bound**
for the infimum of 𝓕 over all (curved) partitions with that topology. Note the identity
𝓕³ = 20 A_tot² Σᵢ Iᵢ / (27 V_box³), A_tot = total interface area (π cancels).

**Minimisation** (`scripts/dfm/`): analytic gradient (numba), single-thread L-BFGS. Vertex motion is *restricted*
as in Surface Evolver's normal motion: foam vertices 3 dof, foam-edge vertices 2 dof (normal plane of the edge),
face vertices 1 dof (normal); between three L-BFGS cycles a tangential Laplacian smoothing (projected on the tangent
plane / edge tangent; Evolver's vertex averaging) keeps the mesh regular (minimum angles stay ≈ 45°). Full 3-dof
motion was tried first: it lets triangles degenerate and needs ~10⁴ iterations. Each level starts from the
subdivided optimum of the previous level. Optional **box relaxation**: a symmetric homogeneous strain of the box (6
dof; 𝓕 is scale invariant) with the exact virial gradient (checked by finite differences).

**Continuum extrapolation.** The sequences 𝓕_L converge as O(h²): successive-difference ratios are 3.88–3.95 for
every structure. We report R_L = 𝓕_L − (𝓕_{L−1} − 𝓕_L)/3 with an uncertainty max(|R_K − R_{K−1}|, |R_K − Aitken|);
for BCC, R₄ and R₅ (an extra level-5 run) agree to 7×10⁻⁸. This error estimate is heuristic, not rigorous.

**Exact certificates** (`scripts/certify_upper_bound.py`). A relaxed mesh is rounded to the dyadic grid 2⁻⁴⁴. In
exact integer/rational arithmetic we check: (i) every cell boundary is a closed oriented surface; (ii) every cone
tetrahedron (apex = rounded centroid of the cell, base = a boundary triangle) has positive volume; (iii) Σᵢ Vᵢ = det L.
By the degree argument (the cones form a closed oriented pseudomanifold mapping to the torus with all local degrees
+1 and total degree Σ Vᵢ / det L = 1), the cells then form a genuine partition with the stated topology. Then
𝓕³ ≤ 20 A_up² Σ Iᵢ / (27 V_box³) with exact rational Iᵢ and triangle areas bounded above by integer ceiling square
roots. This certifies inf 𝓕 ≤ 𝓕_up for that topology and (rounded) box; for non-dyadic boxes the rounding changes the
box by ~10⁻¹⁴ relative, which moves 𝓕 by O(10⁻¹³), far below every margin used. **No lower bounds** are certified.

## 3. Validation against [1] (Milestone 1) — `scripts/check_milestone1.py` (0 FAIL)

*(a) Inputs.* The flat Voronoi values ("Initial Voronoi" rows) are reproduced to the printed 6 decimals for BCC, A15,
C15, H, P, δ, pσ, M with our coordinates (crystallographic files of the archive of [2]) and the caption boxes,
except: H requires the Z12 sites at Wyckoff 4f (as the caption of [1] Fig. 9 says; 4e gives a degenerate Delaunay
complex); M is Pnma with box (1.89, 1, 3.3) (= Pnam 1.89:3.3:1); Z matches with c/a = 0.9935 (caption 0.993). For σ,
μ, C14 our coordinates differ from the paper's (their Voronoi rows are not reproduced); the relaxed DFM value does
not depend on the starting points, only on topology and box.

*(b) Kelvin problem* (equal volumes, minimal area; augmented Lagrangian, `scripts/dfm/constrained.py`): BCC (Kelvin
foam) 1.0972512 (paper 1.097251), A15 (Weaire–Phelan) 1.0935413 (paper 1.093541). A15 equal-volume DFM: ≈1.069742
(paper 1.069743).

*(c) Unconstrained DFM, fixed paper boxes, levels 0–4* (`data/dfm_fixedbox_paper_structures.csv`):

| structure | 𝓕₄ (exact, upper bound) | extrapolated | published [1] | difference |
|---|---|---|---|---|
| σ (a/c = 1.9) | 1.0685800 | 1.0685789(1) | 1.068583 | −4.1×10⁻⁶ |
| P | 1.0688820 | 1.0688810(1) | 1.068909 | −2.8×10⁻⁵ |
| M | 1.0689992 | 1.0689983(1) | 1.069038 | −4.0×10⁻⁵ |
| A15 | 1.0690387 | 1.0690374(1) | 1.069037 | +3.8×10⁻⁷ |
| pσ | 1.0690665 | 1.0690656(1) | 1.069068 | −2.4×10⁻⁶ |
| C14 | 1.0690711 | 1.0690704(1) | 1.069071 | −5.6×10⁻⁷ |
| H | 1.0692046 | 1.0692035(1) | 1.069219 | −1.6×10⁻⁵ |
| C15 | 1.0692441 | 1.0692434(1) | 1.069243 | +4.3×10⁻⁷ |
| Z | 1.0692511 | 1.0692501(1) | 1.069250 | +9.5×10⁻⁸ |
| μ | 1.0695339 | 1.0695330(1) | 1.069662 | −1.3×10⁻⁴ |
| δ | 1.0696766 | 1.0696755(1) | 1.069677 | −1.5×10⁻⁶ |
| BCC | 1.0711078 | 1.0711023(7) | 1.071111 | −8.7×10⁻⁶ |

Reproduced within 2×10⁻⁶ (relative 2×10⁻⁶): **A15, C15, Z, C14, δ** (the cubic A15 and C15, whose box is
unambiguous, within 4.3×10⁻⁷). Our value is **lower** for σ, pσ, H, P, M, μ and BCC. For BCC, P, M, H and μ an
exact certificate shows that the published number is not the infimum for the stated topology and box
(`certificates/claims.json`, `scripts/check_certificates.py`):

| claim | certified upper bound | published |
|---|---|---|
| BCC (cubic box, level 4) | inf 𝓕 ≤ 1.071107810 | 1.071111 |
| P (1.91:3.57:1, level 2) | inf 𝓕 ≤ 1.068897166 | 1.068909 |
| M (1.89:3.3:1, level 2) | inf 𝓕 ≤ 1.069011740 | 1.069038 |
| H (1:3.88:1, Z12 at 2a+4f, level 3) | inf 𝓕 ≤ 1.069208084 | 1.069219 |
| μ (c/a = 5.2, level 2) | inf 𝓕 ≤ 1.069546495 | 1.069662 |

For σ the box is ambiguous: the paper's own CVP row (A 1.094909, I 1.023207) is matched by a/c ≈ 1.93 (the ICSD
ratio 1.9297), not by the caption's 1.9; with a/c = 1.9297 our σ value is 1.0685970 (level 4: 1.0685982; higher than
published by 1.4×10⁻⁵), with 1.9 it is 1.0685789 (lower by 4.1×10⁻⁶). The paper's large low-symmetry cells (P, M, μ, H) also have CVP stretching values
above ours for the same box (P: ours I = 1.025510, paper 1.025525), and our optimal cell volumes are more
dispersed than the paper's histograms (μ: Z16 cells 1.227 V₀ versus ≈1.20). The most likely explanation is
incomplete relaxation in [1] for large cells (their stopping rule is a per-step change below 10⁻⁶), but we cannot
verify their runs. The ranking among the structures of [1] at fixed boxes is unchanged except that **M (1.068998) now
lies below A15 (1.069037)**.

A diagnostic (`scripts/diagnostics/se_quadrature_bias.py`) shows that the one-point facet-centroid rule written in
[1] Eq. (B2) would bias ℐ by ≈ −1.5×10⁻⁴ at our level 4; the agreement for A15/C15/Z/C14 implies that [1] actually
used an exact facet integral (Evolver's `facet_general_integral`), as their text says.

## 4. Screen of the catalogue

### 4.1 Inputs and topology

The 84 CGD files (`ListDATApaper`) and the 37 part-II files (`ListDATApaperStr`) of the Zenodo archive of [2,3]
list, for each FK structure, the space group, Systre cell, the cell centres ("atoms", with coordination 12/14/15/16)
and their neighbours ("edges" = shared faces). We expand the orbits with gemmi (merge tolerance 10⁻³ Å-units: the
files print special positions with a +10⁻⁶ offset, e.g. 0.500001), locate every edge endpoint on an expanded site
(no miss in 121 files), and build the periodic **contact graph**. In three part-II files (N22, Y02, Y06) every edge
line is duplicated (declared coordination 24/28/30/32); the de-duplicated graph is used. For all 121 structures the
4-cliques of the contact graph form a periodic triangulation of the 3-torus (each triangle in exactly two
tetrahedra, Euler characteristic 0, all tetrahedra non-degenerate and Σ|vol| = box volume with the given
coordinates, i.e. embedded), whose dual is a foam with only pentagonal and hexagonal faces and Z − 12 hexagons in
every cell: all 121 are genuine FK tilings (`scripts/check_catalogue.py`: 121 PASS).

*Voronoi test.* For 78/84 (part I) and 27/37 (part II) structures the Delaunay triangulation of the given points
equals the contact triangulation (the initial foam is then their Voronoi diagram). The others — N28, N55, N57, N58,
N59, N60 (as reported by the exploration) and part-II N03, N05, N17, N22, Y01–Y06 — are not Voronoi-realised by the
CGD coordinates; we start them from the barycentric dual (foam vertices at tetrahedron centroids, always embedded),
and, as a second start, from a spring re-embedding of the contact graph (`scripts/dfm/embed.py`: rest lengths
r_Zi + r_Zj with r12 = 0.519, r14 = 0.594, r15 = 0.630, r16 = 0.651 fitted on 11 real FK crystals, log barrier keeping
all tetrahedra positive). After re-embedding, the spring energy per contact is −0.001…0.006 for the
Voronoi-realised structures, 0.006–0.009 for N28 and part-II N03, N05, N17, N22, Y01 (mildly frustrated) and
0.020–0.099 for N55, N57, N58, N59, N60 and part-II Y02–Y06 (strongly frustrated: no embedding with nearly regular
tetrahedra was found). All checks: `scripts/check_catalogue.py` (121 PASS), table `data/catalogue_checks.csv`.

*Cross-check.* Of the nine structures of [1] present in the catalogue, eight (Phys01 C14, Phys02 C15, Phys13 pσ,
Phys14 μ, Phys18 Z, Phys25 H, Phys26 σ, Phys27 A15) were run from both inputs: the catalogue's Systre coordinates and
the crystallographic coordinates/boxes of [1] give the same box-relaxed level-2 value to ≤ 2.3×10⁻¹⁰ (δ = Phys20 was
run once). This validates the CGD pipeline and shows that the H of [1] (Z12 at 4f) is the catalogue's H_comp.

### 4.2 Box

[1] fixes the periodic box at the crystallographic aspect ratios of metallic prototypes. For new structures there is
no such box, and in a melt the lattice parameters relax; we therefore relax the **box shape** (6 strain dof; volume is
irrelevant since 𝓕 is scale invariant) for every structure, including those of [1]. The relaxation preserves the
crystal metric symmetry in all cases inspected (tetragonal, hexagonal and cubic metrics stay so). Box relaxation
lowers 𝓕 by 4×10⁻⁶ (σ, relaxed a/c = 1.909), 2×10⁻⁴ (Z, c/a 0.9935 → 1.027), 3.9×10⁻⁴ (H, b/a 3.88 → 3.73) and about
7×10⁻⁴ for μ (c/a = 5.2 in [1]).

### 4.3 Protocol

`scripts/screen_pipeline.py`: levels 0–2 for all structures; level 3 for those whose level-2 extrapolation lies
within 5×10⁻⁴ of σ's; level 4 for those within 3×10⁻⁴. Results: `data/dfm_screen_boxrelaxed.csv`
(`scripts/make_tables.py`), with all level values (each an exact upper bound for its topology at the relaxed box).

### 4.4 Results (box relaxed)

All 121 catalogue structures plus P, M (not in the catalogue) and BCC were relaxed; the complete table is
`data/dfm_screen_boxrelaxed.csv` and README.md. Leading structures (𝓕_L is an exact upper bound for that topology
in the relaxed box; extrapolated values with the heuristic error of §2; relative difference to σ):

| structure | group, cells | Z12/14/15/16 | L | 𝓕_L | 𝓕 extrapolated | Δ vs σ | 𝒜 | ℐ | V_rms/V₀ |
|---|---|---|---|---|---|---|---|---|---|
| **σ** (Phys26) | P4₂/mnm, 30 | 10/16/4/0 | 4 | 1.068575909 | 1.0685748 ± 1×10⁻⁷ | 0 | 1.090482 | 1.026074 | 0.097 |
| **E11** (new) | I4₁/amd, 44 | 16/20/8/0 | 4 | 1.068809869 | 1.0688088 ± 1×10⁻⁷ | +0.0219% | 1.089788 | 1.028056 | 0.115 |
| H (Phys25) | Cmmm, 30 | 10/16/4/0 | 4 | 1.068810390 | 1.0688092 ± 1×10⁻⁷ | +0.0219% | 1.090394 | 1.026917 | 0.102 |
| μ (Phys14) | R3̄m, 39 | 21/6/6/6 | 4 | 1.068839594 | 1.0688388 ± 1×10⁻⁷ | +0.0247% | 1.087001 | 1.033423 | 0.167 |
| zra-d K₇Cs₆ (Phys12) | P6₃/mmc, 26 | 14/4/4/4 | 4 | 1.068840461 | 1.0688396 ± 1×10⁻⁷ | +0.0248% | 1.086999 | 1.033429 | 0.167 |
| **N56** (new) | I4₁/amd, 76 | 24/44/8/0 | 4 | 1.068869452 | 1.0688683 ± 1×10⁻⁷ | +0.0275% | 1.090757 | 1.026403 | 0.093 |
| P | Pnma, 56 | 24/20/8/4 | 4 | 1.068875254 | 1.0688742 ± 1×10⁻⁷ | +0.0281% | 1.089004 | 1.029727 | 0.133 |
| N16 / N27 (new) | P6₃/mmc 40 / R3̄m 60 | 20/8/8/4 | 3 | 1.06892318 | 1.0689197 ± 6×10⁻⁷ | +0.0323% | 1.087445 | 1.032820 | 0.160 |
| M | Pnma, 52 | 28/8/8/8 | 3 | 1.06892324 | 1.0689199 ± 5×10⁻⁷ | +0.0323% | 1.087130 | 1.033420 | 0.166 |
| N53 (new) | P4₂/mmc, 36 | 14/14/8/0 | 3 | 1.06894210 | 1.0689380 ± 7×10⁻⁷ | +0.0340% | 1.089248 | 1.029458 | 0.128 |
| N03, N02 (new) | P6₃/mmc 38, P6̄m2 19 | 22/4/4/8 | 3 | 1.0689494, 1.0689501 | 1.0689463, 1.0689471 | +0.035% | 1.08655 | 1.03460 | 0.174 |
| N09, N04 (new) | R3̄m 57, P6₃/mmc 38 | 22/4/4/8 | 3 | 1.0689925, 1.0689933 | 1.0689894, 1.0689902 | +0.039% | 1.08649 | 1.03483 | 0.175 |
| pσ (Phys13) | Pbam, 26 | 14/4/4/4 | 3 | 1.06903229 | 1.0690289 ± 5×10⁻⁷ | +0.0425% | 1.087263 | 1.033482 | 0.166 |
| A15 (Phys27) | Pm3̄n, 8 | 2/6/0/0 | 4* | 1.06903865 | 1.0690374 ± 1×10⁻⁷ | +0.0434% | 1.091986 | 1.024581 | 0.054 |
| Z (Phys18) | P6/mmm, 7 | 3/2/2/0 | 3 | 1.06905073 | 1.0690469 ± 6×10⁻⁷ | +0.0442% | 1.088253 | 1.031656 | 0.147 |
| C14 (Phys01) | P6₃/mmc, 12 | 8/0/0/4 | 3 | 1.06907259 | 1.0690700 ± 4×10⁻⁷ | +0.0463% | 1.085652 | 1.036669 | 0.189 |
| C15 (Phys02) | Fd3̄m, 24 | 16/0/0/8 | 4* | 1.06924409 | 1.0692434 ± 4×10⁻⁸ | +0.0628% | 1.085412 | 1.037628 | 0.192 |
| δ (Phys20) | P2₁2₁2₁, 56 | 24/20/8/4 | 2 | 1.06968163 | 1.069665 ± 7×10⁻⁶ | +0.101% | 1.089371 | 1.031364 | 0.139 |
| BCC | Im3̄m, 2 | — | 4* | 1.07110781 | 1.0711023 ± 7×10⁻⁷ | +0.237% | 1.097266 | 1.020647 | 0 |

(* cubic: the relaxed box is the cubic box; values from the fixed-box runs to level 4.)

**Answer (numerical).** σ remains the DFM minimiser over the whole catalogue (84 + 37 structures) and over the 11
structures of [1], once every box is relaxed. The runner-up is a tie between the catalogue entry **E11** (never
considered before) and **H**, 2.34×10⁻⁴ (0.022 %) above σ, followed by μ and zra-d (K₇Cs₆ type, never evaluated in the
DFM) at +0.025 %, N56 (new) at +0.027 % and P at +0.028 %. The gap σ → E11 is about 2000 times the extrapolation
uncertainty and 200 times the whole level-4 discretisation correction of σ (𝓕₄ − 𝓕_∞ ≈ 1.1×10⁻⁶); the exact level-4
upper bound for σ, inf 𝓕 ≤ 1.068575909 (certified, `scripts/certify_upper_bound.py`), is below every extrapolated
competitor value by more than 2.3×10⁻⁴. No structure falls below σ and none lies within 0.02 % of it (E11 and H lie
just outside, at 0.0219 %).

**Robustness.** For σ, H, μ, zra-d, N16, N53 and N02 a second relaxation started from a different initial geometry
(spring re-embedding of the contact graph instead of the catalogue/crystal coordinates) reaches the same level-2 value
to ≤ 2.3×10⁻¹⁰ (`work/logs/robustness_reembed_L2.log`). The relaxed boxes keep the crystal metric (σ: a/c = 1.909;
H: b/a = 3.73, c/a = 1.012; E11: c/a = 5.42 for the body-centred conventional cell).

**Box relaxation reorders the structures of [1].** At fixed boxes [1] ranks σ < P < A15 ≈ M < pσ < C14 < H < C15 < Z <
μ < δ (our fixed-box values: σ < P < M < A15 < pσ < C14 < H < C15 < Z < μ < δ); with relaxed boxes the order is
σ < H < μ < P < M < pσ < A15 < Z < C14 < C15 < δ: H and μ, which were penalised by the metallic aspect ratios, become
the closest competitors of σ.

**Polytype families.** Several pairs are degenerate to 10⁻⁶ or better and have identical cell-volume spectra: μ and
zra-d (8.7×10⁻⁷ apart), N16 and N27 (7×10⁻¹²), N02 and N03, N04 and N09, part-II N12/N21, N10/N11, N20/N24. We
presume these are stacking variants built from the same layers (as C14/C15/C36 are, whose values differ more, by up
to 1.7×10⁻⁴), but we have not checked the layer decompositions.

**Frustrated (non-Voronoi) entries.** *(see §4.5)*

### 4.5 Non-Voronoi (frustrated) entries

Relaxed from the spring re-embedding (levels 0–1, 10-min cap; N28 to level 2) and, for N57/N58, also from the catalogue
coordinates (level 2): mildly frustrated entries stay embedded and relax to 𝓕 ≈ 1.080–1.089 (N28 1.0851 at level 2;
part-II N05, N17, N22, N03, Y01 1.080–1.089 at level 1), i.e. ≥ 1.0 % above σ. Strongly frustrated entries drive a
foam edge to zero length: N57 (1.1030) and N58 (1.0992) stay embedded; for N55, N59, N60 and part-II Y02–Y06 the
relaxation became degenerate (cone volumes ≤ 10⁻¹³, or negative with a collapsing box for Y04–Y06), so we report only
the valid unrelaxed starting partitions (𝓕 = 1.13–1.20 for N55, N59, N60, Y02, Y03; 1.57–1.62 for Y04–Y06;
`data/frustrated_start_values.csv`). All of these lie far above σ; their reported values are upper bounds reached,
not converged infima of their (unstable) topologies.

**Further check.** σ in the box of the paper's CVP row (a/c = 1.9297, fixed): 1.0685970; the E11 restart from the
re-embedded geometry reproduces level 2 (1.0688260).

### 4.6 Candidates for SCFT

Ranked by DFM value (all within 0.03 % of σ): **E11** (I4₁/amd, 44 cells per conventional cell, composition
Z12:Z14:Z15 = 4:5:2, the composition of the J complex Phys22 of [2]; new), **H** (studied as metastable by SCFT at
χN = 25 in [1]), **μ**, **zra-d / K₇Cs₆** (P6₃/mmc, 26 cells; nearly degenerate with μ, same composition and cell-volume
spectrum), **N56** (I4₁/amd, 76 cells, new) and **P**. Next (≤ 0.035 %): N16/N27, M, N53, N03/N02. SCFT is out of scope here.

## 5. What is and is not established (caveats)

* **Certified** (exact integer/rational arithmetic, `scripts/check_certificates.py`): six explicit embedded partitions
  with the stated topologies and boxes, giving inf 𝓕 ≤ 1.071107810 (BCC, cubic), ≤ 1.068897166 (P), ≤ 1.069011740 (M),
  ≤ 1.069208084 (H), ≤ 1.069546495 (μ), all below the values published in [1] for the same topology and box, and
  inf 𝓕(σ) ≤ 1.068575909 (σ, relaxed box). The embeddedness argument (positive cones + Σ Vᵢ = det L ⇒ partition) is a
  standard degree argument written out in §2, not taken from a reference. The boxes of the certificates are the
  dyadic roundings of the intended boxes (relative change ~10⁻¹⁴).
* **Numerical, not certified**: every extrapolated value, every ranking statement, and in particular "σ is the DFM
  minimiser over the catalogue". We have no lower bounds: a structure whose relaxation stopped in a non-global local
  minimum *within its topology* could have a lower true infimum. Mitigations: restarts from different initial geometries
  give identical values for the leading structures (§4.4); every level is a monotone-in-refinement exact upper bound with
  clean O(h²) convergence; the gap σ → next structure (2.3×10⁻⁴) is large compared with all observed numerical
  effects (extrapolation ≤ 10⁻⁶, start dependence ≤ 10⁻⁹).
* **Error estimates** of the extrapolations are heuristic (stability of Richardson values between levels; checked for BCC
  with a fifth level: R₄ and R₅ agree to 7×10⁻⁸).
* **Topology constraint.** As in [1], the cell topology is fixed. The DFM optimum of a topology may sit at the boundary
  of its class (vanishing foam edge): this happens for the strongly frustrated non-Voronoi entries (§4.5), for which the
  reported numbers are the values reached, not converged infima. For all 108 Voronoi-realised structures the shortest
  foam edge stays ≥ 0.29 R₀ (≥ 0.51 R₀ for all structures within 0.1 % of σ), and every level sequence decreases
  monotonically.
* **Box.** [1] fixes the box; we report both (fixed boxes for the validation, relaxed boxes for the screen). Box
  relaxation is restricted to homogeneous strains of the given periodic cell; no supercells were tried, so symmetry
  breaking with a larger period is not excluded.
* **Model scope.** The DFM is the strong-stretching, polyhedral-interface limit of [1]; its ranking need not coincide
  with SCFT at finite segregation (SCFT is out of scope; §4.6 lists candidates).
* **Inputs.** The catalogue files are used as published on Zenodo (the 2010/2012 articles themselves were not read;
  NEEDS.md); the identification of `ListDATApaperStr` with part II rests on the archive's own LaTeX tables.

## References

[1] Reddy et al., PNAS 115 (2018) 10233, arXiv:1806.01726. [2] Dutour Sikirić, Delgado-Friedrichs, Deza, Acta Cryst.
A 66 (2010) 602; data Zenodo 10.5281/zenodo.13691342. [3] Dutour Sikirić, Deza, Struct. Chem. 23 (2012) 1103.
Full list and novelty check: references.md.


## 6. State at stop and how to resume

Work was stopped on request (24 Oct 2026, ~12:45). Done: Milestone 1 (all 12 structures of [1], fixed boxes, levels
0–4); box-relaxed screen of all 124 structures at level 2, level 3 for the 20 within 5×10⁻⁴ of σ, level 4 for the 7 within
3×10⁻⁴ (σ, E11, H, μ, zra-d, N56, P) and for BCC; frustrated entries at levels 0–2; six exact certificates.
Not done: level 3+ for structures more than 0.047 % above σ (their level-2 extrapolations carry ±6×10⁻⁶); converged
values for the degenerate frustrated entries; CVP checks of the paper's δ, pσ, M, μ, H rows; an independent
Surface Evolver cross-check. To resume (checkpoints in `work/runs/<name>__<box>/state.pkl`, git-ignored):
`python3 scripts/run_dfm.py <names> --box relax --maxlev K` continues from the last level;
`python3 scripts/screen_pipeline.py data/screen_order.txt --win3 X --win4 Y` reruns the adaptive selection;
`python3 scripts/make_tables.py` rebuilds the tables. If `work/` is lost, the runs are deterministic and can be
regenerated (hours).
