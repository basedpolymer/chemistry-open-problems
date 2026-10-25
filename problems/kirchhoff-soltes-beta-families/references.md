# References and novelty check (kirchhoff-soltes-beta-families)

## Primary source

- **[GMD26]** K. K. Gottwald, S. Majstorović Ergotić, T. Došlić, *Šoltés problem for the Kirchhoff index of a graph*,
  Appl. Math. Comput. 510 (2026) 129694, doi:10.1016/j.amc.2025.129694 (received 29 Jul 2024, online 29 Aug 2025);
  arXiv:2609.13751v1 (12 Sep 2026).
  - **arXiv v1: read in full.**
  - **Journal version** (publisher PDF read): definitions, Problems 3, 4, 5, 17, Sections 3–4 and the concluding
    remarks checked. Their wording is identical to arXiv v1.
  - Problem 4: journal p. 3, arXiv p. 6. Theorem 15 and Corollary 16: journal pp. 8–10, arXiv pp. 16–20. Problem 17:
    journal p. 11, arXiv pp. 21–22.
  - The definition of a β-Kirchhoff Šoltés graph is "|K(G)| ≥ β|V(G)|", i.e. at least β.
  - Theorem 15 already gives proportions 2k/(3k + 1) → 2/3, hence Problem 4 for every β < 2/3.

## Follow-up by the same authors

- **[MG25]** S. Majstorović Ergotić, K. K. Gottwald, *Graphs that preserve Kirchhoff index after the removal of a vertex*,
  Bol. Soc. Mat. Mex. 31 (2025) 88, doi:10.1007/s40590-025-00780-2 (open access, CC-BY). **Read in full** (publisher
  PDF).
  - Unicyclic and cactus graphs with exactly one or at least two good vertices: Theorems 10, 12, 14, 16.
  - The "bouquet" Theorem 16 uses the splice condition Rt_{G_m}(v5) = 2(m − 1). This is our Lemma 1 with B = C₅, which
    has (δ, ρ) = (0, 2) at the neighbours of the root.
  - No result on proportions of good vertices.
  - Its open problems (its own Problems 17–20) concern forbidden transmission values, cacti, removing vertex sets, and
    mixed Wiener/Kirchhoff good vertices.

## Wiener-index Šoltés literature (context)

- **[KMS18]** M. Knor, S. Majstorović, R. Škrekovski, *Graphs preserving Wiener index upon vertex removal*, Appl. Math.
  Comput. 338 (2018) 25–32, doi:10.1016/j.amc.2018.05.047. Publisher PDF. **Read: first page and
  abstract** (dense graphs and relaxed versions for the Wiener index); not used in the proofs.
- **[BJM21]** J. Bok, N. Jedličková, J. Maxová, *A relaxed version of Šoltés's problem and cactus graphs*, Bull. Malays.
  Math. Sci. Soc. 44 (2021) 3733–3745, doi:10.1007/s40840-021-01144-5, arXiv:1911.10502. Publisher PDF. **Read:
  abstract** (cacti with exactly 2k or k good vertices for the Wiener index); not used in the proofs.
- Akhmejanova et al., *Wiener index and graphs, almost half of whose vertices satisfy Šoltés property*, Discrete Appl.
  Math. 325 (2023) 37–42, doi:10.1016/j.dam.2022.09.021, arXiv:2012.08786. **Read: arXiv abstract** (an infinite
  series of graphs whose proportion of Šoltés vertices tends to 1/2).
- A. A. Dobrynin, K. V. Vorob'ev, *Some results on the Wiener index related to the Šoltés problem of graphs*,
  arXiv:2311.15144. **Read: arXiv abstract** (graphs with a proportion > 1/2 of vertices with W(G) − W(G − v) = m for
  every m ≥ 0, and new families with proportion → 1/2 for m = 0).
- S. Cambie, *Towards the essence of Šoltés' problem*, arXiv:2406.03451. **Read: arXiv abstract** (the variant where
  the total distance decreases by a fixed constant). The sub-agent reported finite 2/3-Šoltés graphs in the full text;
  not checked, not used.
- N. Bašić, M. Knor, R. Škrekovski, *On regular graphs with Šoltés vertices*, arXiv:2303.11996. **Read: arXiv
  abstract** (cubic 2-connected graphs with ≥ 2^r Šoltés vertices; large 1/3-Šoltés graphs).
- Ľ. Šoltés, *Transmission in graphs: a bound and vertex removing*, Math. Slovaca 41 (1991) 11–16. Not read (cited
  through [GMD26]).
- All arXiv ids above were confirmed with the arXiv API (titles and authors), 25 Oct 2026.

## Tools and standard facts

- D. J. Klein, M. Randić, *Resistance distance*, J. Math. Chem. 12 (1993) 81–95, doi:10.1007/BF01164627; I. Gutman,
  B. Mohar, J. Chem. Inf. Comput. Sci. 36 (1996) 982–985, doi:10.1021/ci960007t: Kf = n tr L⁺ (standard, not re-read).
- E. H. Bareiss, *Sylvester's identity and multistep integer-preserving Gaussian elimination*, Math. Comp. 22 (1968)
  565–578, doi:10.1090/S0025-5718-1968-0226829-0 (DOI checked on Crossref). Fraction-free elimination in
  `gadgetbound.c`.
- B. D. McKay, A. Piperno, *Practical graph isomorphism, II*, J. Symb. Comput. 60 (2014) 94–112,
  doi:10.1016/j.jsc.2013.09.003. nauty 2.8.9 `geng`.
- OEIS A001349 (connected graphs): used as a count check.
- Sibling folder `../kirchhoff-soltes-problem/` (Problem 3, partial result). Its modular arithmetic routine (`kscheck.c`)
  is reused in `goodcount.c`.

## Novelty check (25 Oct 2026)

Done by a read-only Claude Opus 5.5 sub-agent (effort not recorded), on the coverage window 2025-01-01 to 2026-10-25.
The Gemini helper was not used: its quota was exhausted.

| Source | Query | Outcome |
|---|---|---|
| OpenAlex W4413802645 | work record, `filter=cites:` | cited_by_count 0; no citing work |
| OpenAlex search | "Kirchhoff Šoltés", "Kirchhoff Soltes" | only [GMD26] and [MG25] relevant |
| Semantic Scholar | citations of the DOI and of arXiv:2609.13751 | 0 |
| arXiv API | "Kirchhoff index" + Šoltés / good vertex / removal / unchanged; "resistance distance" + removal (2025–2026) | only 2609.13751 |
| Crossref | Kirchhoff Šoltés | [GMD26], [MG25]; Cambie's SSRN preprint on z-Šoltés digraphs (Wiener) |
| Web search (4 queries) | β-Kirchhoff Šoltés, Kirchhoff Šoltés graph, … | nothing new (TUM seminar announcement of 2024) |
| GitHub ipitchford (193 repos), mmislan (27 repos) | names/descriptions | nothing on Kirchhoff/Šoltés/resistance |
| `gh search repos` | "Kirchhoff index", soltes | StijnCambie/Soltes (Wiener-index notebooks) only |
| google-deepmind/formal-conjectures | tree + code search | only WienerIndex.lean definitions |
| Zenodo API | "+Kirchhoff +Soltes", "good vertices" Kirchhoff | 0 relevant |

Conclusion: Problem 4 for β ≥ 2/3, exact proportions, and Problem 17 were open as of 25 Oct 2026, as far as could be
checked. [MG25] does not overlap with our results.
