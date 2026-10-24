# References and novelty check

## Primary sources

1. **A. Reddy, M. B. Buckley, A. Arora, F. S. Bates, K. D. Dorfman, G. M. Grason**, *Stable Frank–Kasper phases of
   self-assembled, soft matter spheres*, Proc. Natl. Acad. Sci. USA **115** (2018) 10233–10238,
   doi:10.1073/pnas.1809655115; arXiv:1806.01726 (v1, 5 June 2018).
   **Read: full text of arXiv v1** (`certificates/sources/reddy2018/`, git-ignored), including Appendix A (model,
   Eqs. (1)–(2), (A1)–(A17)), Appendix B (Surface Evolver method, Eqs. (B1)–(B3), "fixed periodicity", convergence
   criterion 1e-6) and the tables of Figs. 7–18 (transcribed to `data/published_reddy2018.csv`). The PNAS version
   and SI were not read (see NEEDS.md).
2. **M. Dutour Sikirić, O. Delgado-Friedrichs, M. Deza**, *Space fullerenes: a computer search for new Frank–Kasper
   structures*, Acta Cryst. A **66** (2010) 602–615, doi:10.1107/S0108767310022932.
   **Not read** (IUCr returns HTTP 403 to the container, ANU repository copy embargoed): only abstract/metadata.
   **Data used**: *Dataset and programs for the work on space-fullerenes*, Zenodo record 13691342
   (doi:10.5281/zenodo.13691342, 5 Sep 2024, file `SpaceFullerene.tar.gz`, 182,836,409 bytes,
   md5 582371cee665c62e01a2f47a73dd2c9d), folder `ListDATApaper` (84 CGD files: E01–E11, N01–N60, Phys01, 02, 03, 04,
   06, 12, 13, 14, 18, 20, 25, 26, 27) and the LaTeX tables `LatexFiles/*.tex` (names of the physical phases:
   Phys01 C14, 02 C15, 03 C36, 04 6-layer, 06 9-layer, 12 zra-d (K7Cs6), 13 pσ (Th6Cd7), 14 μ, 18 Z, 20 δ, 25 H,
   26 σ, 27 A15). Note: the repository is on Zenodo (linked from the author's publication page
   https://mathieudutsik.github.io/Publications/), not on GitHub as stated in the brief.
3. **M. Dutour Sikirić, M. Deza**, *Space fullerenes: computer search for new Frank–Kasper structures II*,
   Struct. Chem. **23** (2012) 1103–1114, doi:10.1007/s11224-012-0005-3. **Not read** (paywalled). The same Zenodo
   archive contains folder `ListDATApaperStr` with 37 CGD files (E01–E05, N01–N26, Y01–Y06) whose groups match the
   21-cell tables `LatexFiles21/*.tex` (e.g. N01 = entry 1, `idx=[21,1]`, Pmc2₁); we treat them as the part-II
   structures (37 new structures with 21 cells in the reduced fundamental domain, as stated in the part-II abstract).
   Contrary to the brief, the part-II data **are** in the archive.

## Methods and validation references (read: abstract or known content only)

4. K. A. Brakke, *The Surface Evolver*, Exp. Math. **1** (1992) 141–165 — the tool used in [1]; not used here (we use
   our own exact polyhedral evaluation, see result.md).
5. D. Weaire, R. Phelan, *A counter-example to Kelvin's conjecture on minimal surfaces*, Phil. Mag. Lett. **69**
   (1994) 107–110 — Weaire–Phelan foam (normalised area 1.093541 in [1]); reproduced (claim (b) of
   `scripts/check_milestone1.py`).
6. R. Kusner, J. M. Sullivan, *Comparing the Weaire–Phelan equal-volume foam to Kelvin's foam*, Forma **11** (1996)
   233–242 — flat weighted-Voronoi A15 formulae (A14)–(A16) quoted in [1].
7. N. Kashyap, D. L. Neuhoff, *On quantization with the Weaire–Phelan partition*, IEEE Trans. Inf. Theory **47**
   (2001) 2538 — second moment of A15 partitions, (A15) of [1].
8. Crystallographic inputs: coordinates from the physical-phase CGD files of [2]'s archive
   (`PhysicalPhases_V1/DONE/*`, ICSD-type data) and the cell ratios printed in [1] Figs. 7–18. The original
   crystallographic papers cited by [1] (Bergman–Shoemaker 1954 σ; Wilson et al. 1960 Z; Ye–Li–Kuo 1984 H;
   Shoemaker et al. 1957 P, 1963 δ, 1967 M; Fornasini et al. 1984 pσ; De Graef–McHenry 2012 A15/C14/C15/μ) were not
   read (NEEDS.md).

## Novelty check (24 Oct 2026)

Question: has anyone evaluated the diblock foam model (DFM) of [1] on the Frank–Kasper catalogue [2,3], or otherwise
screened FK structures beyond the 11 of [1] with the DFM?

- **Citers of [1]** (Semantic Scholar API, `paper/DOI:10.1073/pnas.1809655115/citations`, 141 returned; OpenAlex
  reports 158 but its shared daily budget was exhausted from this IP): all titles read
  (`work/novelty/s2_10_1073_pnas_1809655115.txt`, git-ignored). None evaluates the DFM on further structures. The
  closest are SCFT studies of the classical FK phases (e.g. He & Wang, Polymers 16 (2024) 372; Zhuang et al.,
  Macromolecules 2024, doi:10.1021/acs.macromol.4c02119; Xie & Shi 2021–2023), interfacial-geometry analyses
  (Collanton & Dorfman, Phys. Rev. Mater. 6 (2022) 015602; Greenvall, Dimitriyev & Grason, J. Polym. Sci. 2023) and
  reviews (Dorfman, Macromolecules 2021; Magruder & Dorfman 2024).
- **Citers of [2]** (Semantic Scholar: 60) and **of [3]** (16): clathrate hydrates, MOFs, foams (Cox, Graner, Mosseri,
  Sadoc, J. Phys.: Condens. Matter 29 (2017) 114001 — Kelvin-problem areas, equal volumes), alloys, graph indices; two
  block-copolymer papers (Magruder & Dorfman 2021 C36; Mueller et al. 2021) that do not use the DFM. No DFM screen.
- **Groups** (Grason, Dorfman, Bates, W. Li, A.-C. Shi): later papers found through the citer lists concern SCFT, chain
  trajectories, networks; none re-evaluates the DFM over a structure catalogue.
- **AI-driven discovery 2024–2026**: web search found AI/SCFT automated phase-diagram construction
  (Chinese J. Polym. Sci. 2025, doi:10.1007/s10118-025-3396-5) and inverse-design studies, all SCFT on standard
  candidate lists; no DFM catalogue screen.
- **arXiv** (API: `all:"diblock foam"`, `abs:"Frank-Kasper" AND abs:foam`, `abs:"Frank-Kasper" AND abs:"unequal"`):
  only [1] and Cox et al. 2016.
- **Zenodo** (`"Frank-Kasper"`, `diblock foam model`, `space fullerenes Frank-Kasper`, Evidence Press): only the data
  archive of [2] and unrelated records.
- **GitHub** user listings (`mmislan`, `ipitchford`) are blocked by the session's API restrictions (repository-scoped
  endpoints only); raw-file access to `ipitchford/ipitchford` returned 404. Zenodo mirrors of Evidence Press showed
  nothing related.
- **Web search** ("diblock foam model" Frank-Kasper; FK enumeration + block copolymer screening; ML discovery of FK
  phases): nothing beyond the above.

Conclusion: no prior DFM evaluation of the catalogue structures was found (bounded search; OpenAlex and GitHub-user
checks incomplete for the reasons stated).
