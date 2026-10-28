# Sources and novelty audit

Audit date: 28 October 2026 (Europe/Paris). Solver: GPT-6-Astra.

## Primary source

Joseph Klobusicky and Matthew Rakauskas, *Coagulation equations with particle emission*, arXiv:2604.14076v1, 15 April 2026. DOI: [10.48550/arXiv.2604.14076](https://doi.org/10.48550/arXiv.2604.14076).

Read: full PDF extracted locally, with detailed examination of Definition 1, Theorem 1, Theorem 5, equation (37), the monodisperse formulas in Section 6, and Sections 7--8. Relevant locations: pp. 5--6 (solution and gelation definitions), p. 8 (large-cluster equations), p. 9 (short-time theorem), pp. 17--19 (moments and polynomial solutions), p. 21 (conjectured threshold), p. 22 (conclusion).

The quote used in the manuscript is the complete sentence proposing the monodisperse threshold, beginning “Formula (96) suggests”, Section 7, p. 21. Its formula is `t_gel = k/(2k-ell)`. The same page conjectures the dimer threshold `2/3`. The displayed second-moment formula is already in the source; the contribution here is an explicit full distribution and a proof of existence, summability and uniqueness in the pre-gelation class all the way to the proposed threshold. The numerical value itself is not a new prediction.

## Public search

- [arXiv record](https://arxiv.org/abs/2604.14076): only v1 was displayed on the audit date.
- [Author's publication list](https://joeklobusicky.net/publications/): the paper is listed as submitted; no follow-up proving this result appeared in the inspected list.
- [OpenAlex work W7154644347](https://openalex.org/W7154644347): returned zero indexed citations. The broad title-search result count is **not** a count of citing works; the relevant work was selected by its exact title and DOI.
- Semantic Scholar API by arXiv identifier: unavailable, HTTP 429. No claim that its citing literature was checked.
- Zenodo API exact identifier search: zero records.
- GitHub repository searches for the exact arXiv identifier and `coagulation emission`: zero repositories returned.
- Public GitHub repository metadata for `mmislan`: 27 repositories, all pages exhausted. Public metadata for `ipitchford`: 204 repositories, all pages exhausted. No name or description matched coagulation, particle emission, Klobusicky, Rakauskas or gelation. This was a metadata screen, not inspection of every file in every repository.
- Web queries on the exact paper title, author names with gelation/solution, and particle emission with coagulation/explicit solution found the primary paper and summaries, but no later resolution. Third-party summaries were not used as mathematical evidence.

The API queries, URLs, results and failures are saved in `certificates/novelty-api.json`; the read-only script is `scripts/novelty_audit.py`. A lack of indexed hits does not establish absolute priority. The justified statement is **no prior resolution was found in this bounded audit**.

The method uses classical Lagrange inversion and Stirling's formula, with the exact coefficient extraction and asymptotic calculation displayed in the proof. No originality is claimed for those tools or for the previously proposed moment formulas.
