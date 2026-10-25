# References and novelty check

| Source | Read | Use |
|---|---|---|
| M. Banaji, B. Boros, J. Hofbauer, *Bifurcations in planar, quadratic mass-action networks with few reactions and low molecularity*, Nonlinear Dyn. 112 (2024) 21425–21448, doi:10.1007/s11071-024-10068-1, arXiv:2406.13451 | full text (arXiv v1 and published version, open access) | Open statement (§7, p. 22 of v1 / p. 21446 published); definitions (§2.1–2.2); parametrisation of positive equilibria (§2.4–2.6, Lemma 11); Lemma 24 (838 classes), Theorem 25, Remark 26, Lemma 27, Table 1 (orange bullets read on the rendered PDF, p. 12 of v1) |
| T. Wilhelm, *The smallest chemical reaction system with bistability*, BMC Syst. Biol. 3 (2009) 90, doi:10.1186/1752-0509-3-90 | abstract (through BBH and a web summary) | Wilhelm's network (9) of BBH, in one of the ten classes |
| G. Craciun, C. Pantea, *Identifiability of chemical reaction networks*, J. Math. Chem. 44 (2008) 244–259, doi:10.1007/s10910-007-9307-x | not read (cited by BBH for dynamical equivalence) | definition only, taken from BBH §2.2 |
| M. A. Krasnosel'skii, P. P. Zabreiko, *Geometrical Methods of Nonlinear Analysis*, Springer 1984, doi:10.1007/978-3-642-69409-7 | not consulted here (classical result) | index of an asymptotically stable equilibrium = 1; the manuscript reproves the needed half-index version |
| L. Perko, *Differential Equations and Dynamical Systems*, 3rd ed., Springer 2001, doi:10.1007/978-1-4613-0003-8 | standard textbook | stable manifold theorem (§2.7), index theory |
| R. T. Rockafellar, *Convex Analysis*, Princeton 1970 | standard | relative interior of a sum of cones (Cor. 6.6.2) |
| B. Boros, GitHub `balazsboros/reaction_networks` | not used | BBH's Mathematica notebooks (not used by our computations) |

## Novelty check (25 Oct 2026)

- **Citing works.** OpenAlex lists 6 works citing the Nonlinear Dyn. paper (W4401634851); titles and abstracts were
  checked, and none addresses multistability of bimolecular (2, 4, 2) networks:
  - Plesa, *Mapping dynamical systems into chemical reactions* (Nonlinear Dyn. 2025);
  - a paper on epidemiologic strains, siphons and autocatalytic cores (Mathematics 2025);
  - two papers on the equivalence of zero-one reaction networks (2025, 2026);
  - BBH, *The inheritance of local bifurcations in mass action networks* (J. Nonlinear Sci. 2025);
  - *First-order endotactic reaction networks* (J. Math. Biol. 2026).

  Semantic Scholar was rate-limited.
- **Later work of the authors.** arXiv listings for Banaji, Boros and Hofbauer up to Oct 2026 were checked.
  - Banaji, arXiv:2609.07265 (fold and cusp bifurcations via dual equations): the full text was grepped. It mentions
    (2, 4, 2) networks only for cusp bifurcations and three positive equilibria (BBH Remark 38), not for this question.
  - Boros et al., arXiv:2510.03621 (disguised toric locus): unrelated.
- **Related work.** Jiao–Tang–Zeng, *Multistability of small zero-one reaction networks*, arXiv:2406.11586: the full
  text was grepped. It treats zero-one networks (no 2X or 2Y complexes) and positive steady states, so it does not
  answer the question.
- **AI campaigns.**
  - `mmislan` (27 Palomar repositories): nothing on this question.
  - `ipitchford` / Evidence Press (193 repositories): the CRN-related titles were checked (`quadratic-mass-action-sharpness`,
    `reverse-reaction-attracting-oscillations`, `tcell-exactly-five`, …). Nothing on bimolecular multistability or BBH
    Table 1.
  - Zenodo API, searched for the key phrases, the DOI and the arXiv id: only unrelated records, the closest being
    "Minimal canard networks in planar mass-action systems" (2026).
  - `gh search code "2406.13451"`: only this repository and RSS mirrors.
- **Web search** on the key phrases (open problem, Wilhelm smallest bistable system): no prior answer found.
- The literature-search helper (Gemini via agy) was unavailable (quota exhausted), so no Gemini output was used.

Conclusion: no prior answer was found; the question appears open as of 25 Oct 2026.
