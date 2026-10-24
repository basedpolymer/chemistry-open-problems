# References and novelty check (kirchhoff-soltes-problem)

## Primary source (read in full)

- K. K. Gottwald, S. Majstorović Ergotić, T. Došlić, *Šoltés problem for the Kirchhoff index of a graph*,
  Appl. Math. Comput. 510 (2026) 129694, doi:10.1016/j.amc.2025.129694 (online 29 Aug 2025 per OpenAlex);
  arXiv:2609.13751v1 (12 Sep 2026). **Read: full text of arXiv v1** (PDF in
  `certificates/sources/arxiv-2609.13751/`, git-ignored). The journal version was not accessible (paywall);
  see NEEDS.md. Open statement: Problem 3 (p. 5), search ranges (p. 9), concluding remarks (p. 21), quoted in
  README.md.

## Novelty check (24 Oct 2026)

| Source | Query | Outcome |
|---|---|---|
| OpenAlex W4413802645 | work record | cited_by_count = 0 (citing-works endpoint rate-limited that day) |
| Semantic Scholar | citations of DOI and of arXiv id | none |
| arXiv API | abs:"Kirchhoff index", newest 40 (to 2026-10-05) | no paper on the Kirchhoff Šoltés problem other than 2609.13751 |
| Crossref | "Šoltés problem Kirchhoff index" | only the source and unrelated works |
| Web search (2 queries) | Kirchhoff Šoltés / good vertex | source; Majstorović Ergotić–Gottwald 2025 (below); TUM seminar announcement (30 Jul 2024) of the same work |
| Zenodo API | "Kirchhoff Soltes", "Kirchhoff index vertex removal", "Soltes problem" | nothing relevant |
| GitHub ipitchford (Evidence Press) | repository search kirchhoff / soltes / graph | nothing on this problem (9 unrelated graph repos) |
| GitHub mmislan (Palomar) | repository list | CRNT only |
| google-deepmind/formal-conjectures, openai/math | — | not reachable from this session (GitHub API restricted); web code search requires login |

Conclusion: the problem was open as of 24 Oct 2026 as far as could be checked.

## Related work on the Kirchhoff version

- S. Majstorović Ergotić, K. K. Gottwald, *Graphs that preserve Kirchhoff index after the removal of a
  vertex*, Bol. Soc. Mat. Mex. 31 (2025) 88, doi:10.1007/s40590-025-00780-2. **Read: abstract only**
  (Crossref). Relaxed version (unicyclic and cactus graphs with one or two good vertices); states that only
  C₅ is known for the full problem.

## Wiener-index Šoltés problem (background, methods)

- Ľ. Šoltés, *Transmission in graphs: a bound and vertex removing*, Math. Slovaca 41 (1991) 11–16. Not read
  (cited through the source).
- M. Knor, S. Majstorović, R. Škrekovski, *Graphs whose Wiener index does not change when a specific vertex
  is removed*, Discrete Appl. Math. 238 (2018) 126–132, doi:10.1016/j.dam.2017.12.012. Abstract (web search).
- M. Knor, S. Majstorović, R. Škrekovski, *Graphs preserving Wiener index upon vertex removal*, Appl. Math.
  Comput. 338 (2018) 25–32, doi:10.1016/j.amc.2018.05.047. Abstract (web snippet: dense graphs are not
  solutions — the Wiener analogue of our Theorem 8).
- J. Bok, N. Jedličková, J. Maxová, *A relaxed version of Šoltés's problem and cactus graphs*, Bull. Malays.
  Math. Sci. Soc. 44 (2021) 3733–3745, doi:10.1007/s40840-021-01144-5. Not read.
- M. Akhmejanova et al., *Wiener index and graphs, almost half of whose vertices satisfy Šoltés property*,
  Discrete Appl. Math. 325 (2023) 37–42, doi:10.1016/j.dam.2022.09.021. Not read.
- N. Bašić, M. Knor, R. Škrekovski, *On regular graphs with Šoltés vertices*, arXiv:2303.11996. Abstract
  (computer search on vertex-transitive graph collections found no 1-Šoltés graph).
- M. Knor, R. Škrekovski, A. Tepeh, *Selected topics on Wiener index*, arXiv:2303.11405. Abstract.
- S. Cambie, *Towards the essence of Šoltés' problem*, arXiv:2406.03451; S. Cambie, A. Tiwari, *Uniform
  Šoltés' hypergraphs and Šoltés' weighted graphs*, arXiv:2506.07511. Abstracts.

## Tools and facts used

- D. J. Klein, M. Randić, *Resistance distance*, J. Math. Chem. 12 (1993) 81–95, doi:10.1007/BF01164627
  (Kf = n tr L^+; not re-read, standard).
- I. Gutman, B. Mohar, *The quasi-Wiener and the Kirchhoff indices coincide*, J. Chem. Inf. Comput. Sci. 36
  (1996) 982–985, doi:10.1021/ci960007t (Kf = n Σ 1/μ_i; standard).
- M. Fiedler, *Algebraic connectivity of graphs*, Czechoslovak Math. J. 23 (1973) 298–305,
  doi:10.21136/CMJ.1973.101168 (a(G − v) ≥ a(G) − 1 and a(G) ≤ κ(G); both re-proved or standard; the first
  is re-proved in result.md).
- N. Biggs, *Potential theory on distance-regular graphs*, Combin. Probab. Comput. 2 (1993) 243–255,
  doi:10.1017/S096354830000064X (resistances in DRGs from the intersection array; our layer-current
  derivation is self-contained). Not read.
- A. E. Brouwer, W. H. Haemers, *Spectra of Graphs*, Springer 2012 (SRG eigenvalue facts: integrality unless
  conference graph; interlacing). Standard; not re-read.
- B. D. McKay, A. Piperno, *Practical graph isomorphism, II*, J. Symb. Comput. 60 (2014) 94–112,
  doi:10.1016/j.jsc.2013.09.003 — nauty 2.8.9 (`geng`, `multig`), downloaded from
  users.cecs.anu.edu.au/~bdm/nauty (sha256 in certificates/sources/nauty/nauty.sha256).
- D. Holt, G. Royle, *A census of small transitive groups and vertex-transitive graphs*, J. Symb. Comput. 101
  (2020) 51–60, doi:10.1016/j.jsc.2019.06.006; data: *Vertex-transitive graphs on fewer than 48 vertices*,
  Zenodo, doi:10.5281/zenodo.4010122 (all tar files downloaded; sha256 in
  certificates/sources/holt-royle/sha256.txt).
- OEIS A002218 (2-connected graphs), A006799 (vertex-transitive graphs): counts used as cross-checks.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms* — not used: all certificates are exact
  (modular or rational); floating point appears only in labelled explorations.
