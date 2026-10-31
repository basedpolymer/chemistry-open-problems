# References and novelty check (all-beta stage of result 23)

## Primary source

- **[GMD26]** K. K. Gottwald, S. Majstorović Ergotić, T. Došlić, *Šoltés problem for the Kirchhoff index of a graph*,
  Appl. Math. Comput. 510 (2026) 129694, doi:10.1016/j.amc.2025.129694; arXiv:2609.13751v1.
  **Read:** journal text (`../certificates/sources/`), definitions, Problems 3, 4, 5 and 17, Theorem 15 and the concluding remarks.
  The quotations of Problems 4, 5 and 17 in `paper.tex` and `README.md` were re-copied from that text on 31 October 2026
  (journal pp. 3 and 11).

## Earlier stages in this repository

- **Result 23, earlier stage** (`../paper.tex`). **Read in full.** Used: the splice lemma (re-proved in `paper.tex`, Lemma 3)
  and the statement of what was open. Its near-miss section on rook graphs motivated the search. Theorem 1 of the present
  folder does not depend on any computation of that stage.
- **Result 21** (`../../kirchhoff-soltes-problem/`), the case β = 1. **Not used.**

## Tools and standard facts

- D. J. Klein, M. Randić, *Resistance distance*, J. Math. Chem. 12 (1993) 81–95, doi:10.1007/BF01164627; I. Gutman, B. Mohar,
  J. Chem. Inf. Comput. Sci. 36 (1996) 982–985, doi:10.1021/ci960007t. Standard definitions, not re-read.
- Ľ. Šoltés, *Transmission in graphs: a bound and vertex removing*, Math. Slovaca 41 (1991) 11–16. Not read (cited through [GMD26]).
- S. Majstorović Ergotić, K. K. Gottwald, Bol. Soc. Mat. Mex. 31 (2025) 88, doi:10.1007/s40590-025-00780-2. **Read at the earlier stage;**
  it has no result on proportions of good vertices.

## Novelty check, 31 October 2026

Done before the attack, by the solver, on the coverage window 25 October to 31 October 2026 (the earlier check, in
`../references.md`, covers 2025-01-01 to 2026-10-25). Raw outputs: `certificates/novelty_check.log`,
`certificates/novelty_check_zenodo.log`; scripts: `scripts/novelty_check.py`, `scripts/novelty_check_zenodo.py`.

| Source | Query | Outcome |
|---|---|---|
| OpenAlex W4413802645 | `filter=cites:` | cited_by_count 0; no citing work |
| Semantic Scholar | citations of the DOI and of arXiv:2609.13751 | 0 and 0 |
| Crossref | `is-referenced-by-count` of the DOI | 0 |
| OpenAlex search | "Kirchhoff Šoltés", "Kirchhoff Soltes", "Kirchhoff index good vertices", "Kirchhoff index vertex removal preserving", sorted by date | only [GMD26] and [MG25] relevant |
| arXiv API | "Kirchhoff"+"Soltes"; "Kirchhoff index"+ good vertices / vertex removal / vertex deletion / removal of a vertex; "resistance distance"+"Soltes"; "Soltes" | only 2609.13751 relevant |
| Zenodo API | Soltes (26 hits, all listed), Šoltés (24, all listed), Soltes Kirchhoff, Šoltés Kirchhoff, Soltes Wiener index, "Wiener index" "vertex removal", "resistance distance" "vertex deletion", "Kirchhoff index" + "good vertex", "deleting a vertex", "vertex deleted", "rook graph", "open problem", proportion; `+Kirchhoff +Soltes`, `+Kirchhoff +"good vertices"`, "Soltes problem" | nothing relevant. Hits with Kirchhoff and vertices in title or description (first 100 per query, most recent first, screened by keyword): an author named Soltes, finite-element software, Picón's interactome atlas (21977553), Ou's spectral functional (21023394). No AI-assisted or author-less working paper on this problem |
| GitHub ipitchford | 204 repositories, names and descriptions | nothing on Kirchhoff, Šoltés or resistance distance |
| GitHub mmislan | 27 repositories, names and descriptions | nothing relevant |
| `gh search repos` | "kirchhoff soltes", soltes, "Kirchhoff index" (by update), "soltes graph" | StijnCambie/Soltes (Wiener index) only; other Kirchhoff-index repositories concern edge augmentation, network vulnerability, simplicial complexes |
| `gh search code` | "Kirchhoff Soltes"; "Kirchhoff index" "good vertices" | this repository; an RSS mirror of the arXiv listing of [GMD26] |
| google-deepmind/formal-conjectures | code search Kirchhoff, Soltes | no hit |
| Web search (2 queries, 1 extended) | Kirchhoff Šoltés proportion of good vertices; β-Kirchhoff Šoltés solved construction | [GMD26], [MG25], the 2024 Croatian Mathematical Congress abstract, Wiener-index papers; nothing beyond proportion 2/3 |

Not completed: the Semantic Scholar *keyword* search returned HTTP 429 three times (rate limit); its citation lookups above did
answer. Google Scholar was not queried. This is a bounded search, not a proof that no such result exists.

**Conclusion.** No work found goes beyond proportion 2/3 for the Kirchhoff index. On the evidence above, Problem 4 for
β ≥ 2/3, Problem 5 for γ > 2/3, and Problem 17 were open on 31 October 2026.
