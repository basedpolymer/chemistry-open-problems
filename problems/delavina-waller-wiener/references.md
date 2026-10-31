# Sources and novelty audit

Audit date: 31 October 2026. Auditor: Claude Opus 5.5 (effort not recorded).

**Outcome.** The conjecture is claimed solved by Mingchang Liu (Zenodo, 19 August 2026). This is the only full proof
found. It is unrefereed. No refutation or erratum was found. By PROTOCOL §1.3 the work stops here.

## Sources read

| Source | DOI / id | What was read |
|---|---|---|
| E. DeLaViña, B. Waller, *Spanning trees with many leaves and average distance*, Electron. J. Combin. 15 (2008) #R33 | [10.37236/757](https://doi.org/10.37236/757) | Full PDF pp. 12–14; verbatim source of Conjecture 7 (p. 14) |
| D. Pandey, P. S. Ravi, *On a conjecture of DeLaViña and Waller*, arXiv:2605.24855v1, 24 May 2026 | [arXiv](https://arxiv.org/abs/2605.24855) | Full text, 18 pages |
| M. Knor, R. Škrekovski, A. Tepeh, *Selected topics on Wiener index*, Ars Math. Contemp. 24 (2024) #P4.07 | [10.26493/1855-3974.3077.63a](https://doi.org/10.26493/1855-3974.3077.63a), arXiv:2303.11405 | arXiv version, passage on p. 13 only (Conjecture 29 there; the brief gives 5.2 for the journal version, not checked) |
| Mingchang Liu, *The DeLaViña–Waller conjecture on the Wiener index*, working paper, Zenodo, 19 August 2026 | [10.5281/zenodo.22015517](https://doi.org/10.5281/zenodo.22015517) (concept DOI 10.5281/zenodo.22015516) | Full text, 37 pages; `verify_certificates.py` and `verify_direct_bellman_majorant.py` downloaded as text, header and imports read, **not executed** |
| VibeMathed entry on the conjecture | [page](https://vibemathed.com/problem/the-delavina-waller-conjecture-on-the-wiener-index) | Web page; lists Liu's claim as "Candidate (review pending)", verification "Unreviewed", submitted 20 August 2026 |
| GitHub `Kevin-OS7/math`, README | [page](https://github.com/Kevin-OS7/math) | Web page; reports an exhaustive check up to 11 vertices, no proof, and says "an unreviewed proof claim exists on Zenodo" |

Cited through the above and not read: J. Plesník, J. Graph Theory 8 (1984) 1–21, doi:10.1002/jgt.3190080102; H. Hua,
Bull. Aust. Math. Soc. 98 (2018) 14–17, doi:10.1017/S0004972718000163; Q. Sun, B. Ikica, R. Škrekovski, V. Vukašinović,
Appl. Math. Comput. 356 (2019) 438–448, doi:10.1016/j.amc.2019.03.025. The first two are not needed, since Lemma 3.1 of
Liu is self-contained. No paywalled paper is needed, so there is no `NEEDS.md`.

## Exact sentences

DeLaViña–Waller, p. 14:

> "Conjecture 7. Let G be a graph with diameter d > 2 and order 2d + 1. Then W(G) ≤ W(C(2d + 1))."

Pandey–Ravi, abstract and p. 1 of arXiv:2605.24855v1:

> "In this context, a conjecture proposed by DeLaViña and Waller has remained open for the last 18 years."

Knor–Škrekovski–Tepeh, arXiv:2303.11405, p. 13:

> "Conjecture 29. Let G be a graph with diameter d > 2 and order 2d + 1. Then W(G) ≤ W(C2d+1), where C2d+1 denotes the
> cycle of length 2d + 1."

Liu, abstract:

> "We prove the conjecture for every d≥3 and determine all equality cases"

## Public search

- Semantic Scholar: citers of arXiv:2605.24855: empty list. Citers of the EJC paper since 2024: four, namely Pandey–Ravi,
  arXiv:2507.17885 (trees), a conditional-diameter paper, and a 2026 proof of Graffiti.pc Conjecture 2, which is a
  different conjecture from the same paper.
- OpenAlex, first run only: this search is what found Liu's record (work W7203882583, 0 citers); arXiv:2605.24855 had 0
  citers. Knor–Škrekovski–Tepeh has 17 citers, of which only the first few were looked at. **Limitation:** the anonymous
  daily budget of OpenAlex, shared by this IP, then ran out, so `certificates/novelty-api.json` holds error entries for
  OpenAlex list queries. Rerun `scripts/novelty_audit.py` after midnight UTC to complete them. No API key was requested.

- arXiv API (queries in `certificates/novelty-api.json`, plus title, author and abstract queries repeated by hand on the
  audit date): no version of Liu's paper and no other claimed proof. The only abstract match on the conjecture is
  Pandey–Ravi. `Wiener index` titles since June 2026: only arXiv:2607.20383, on another subject.
- Zenodo API: the exact phrase returns the one Liu record, version 1 only. The record's concept-DOI number resolves only
  through doi.org, not through `/api/records/22015516` (HTTP 410). Identifier 2605.24855 returns no record. The hit "A
  Proof of Written on the Wall II Conjecture 100" (10.5281/zenodo.21914031) came only from a loose keyword match; it was
  judged from its title and not opened.
- GitHub: `ipitchford`, 204 repositories, and `mmislan`, 27 repositories, all pages exhausted; names and descriptions
  screened for Wiener, DeLaViña, average distance, Plesník and Graffiti. The only hits are TxGraffiti repositories on
  another conjecture. This screens metadata and does not inspect every file. `gh search repos` on the key phrases found no
  proof. `google-deepmind/formal-conjectures` has Wiener-index definitions and Written on the Wall tests, but a code search
  shows no DeLaViña–Waller statement.
- Web search on the key phrases and on the author's name: the two pages in the table above, and nothing later.

A lack of hits does not establish that nothing else exists. The justified statement is that **Liu's preprint is the only
resolution found in this bounded audit, and no objection to it was found**.

The scripts `audit_identities.py`, `audit_scalar.py` and `audit_graphs.py` were written from the paper's statements. The
reduction to a block with attached trees and the per-level Cauchy bound were found independently in this session before the
paper was opened, and they are the same as in Liu's Sections 3 and 6, so they add nothing to the record.
