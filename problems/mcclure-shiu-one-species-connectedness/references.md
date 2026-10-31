# Sources and novelty audit

Audit date: 31 October 2026. Auditor: Claude Opus 5.5 (effort not recorded).

**Outcome.** The question is claimed solved by Patrick White (preprint dated 7 September 2026; code public since
25 June 2026). This is the only proof found. It is unrefereed. No refutation or erratum was found. By PROTOCOL §1.3 the
work stops here. No paywalled paper is needed, so there is no `NEEDS.md`.

## Sources read

| Source | DOI / id | What was read |
|---|---|---|
| A. McClure, A. Shiu, *On the connectedness of multistationarity regions of small reaction networks*, J. Symbolic Comput. 125 (2024) 102323 | [10.1016/j.jsc.2024.102323](https://doi.org/10.1016/j.jsc.2024.102323), arXiv:2303.03960v2 | Full text of arXiv v2, 16 pages, rendered pages (the arrow glyphs were read on p. 14); journal version not read. Crossref confirms authors, volume and article number |
| P. White, *Connectedness of multistationarity regions of one-species reaction networks*, draft, 7 September 2026 | [PDF](https://read.byclaude.net/crnt-connectivity.pdf) | Full text, 14 pages, rendered pages, line by line |
| GitHub `pw/crnt-connectivity` | [repository](https://github.com/pw/crnt-connectivity) | `README.md`; commit list; `lean/CRNTConnectivity.lean` as text only (definitions, the two main theorem statements, word search for `sorry`, `axiom`, `admit`, `native_decide`); **not compiled**. The Python scripts were neither read nor run |
| B. Joshi, A. Shiu, *Which small reaction networks are multistationary?*, SIAM J. Appl. Dyn. Syst. 16 (2017) 802–833 | arXiv:1603.06441v2 | Example 2.4, Definitions 3.3–3.4, Theorem 3.6, Corollary 3.10 |

Cited through the above and not read: S. Otake, T. Shaska (quadrinomial discriminants, reference [21] of the source);
E. Feliu, M. L. Telek, Adv. Math. 408 (2022) 108582; M. L. Telek, E. Feliu, PLoS Comput. Biol. 19 (2023) e1010970.

## Exact sentences

McClure–Shiu, Section 6, p. 14 of arXiv v2:

> "Indeed, it is an open question whether, for one-species networks with four or five reactions (e.g., {0 ⇆ A, 2A ⇆ 3A}), the multistationarity regions are always connected."

The next sentence reads "(An answer to this question may involve discriminants of quadrinomials [21].)" Table 1 lists
"1 | 6 or more | Yes or No".

Definition 3.3(1) of the source: a system "is multistationary ... if there exists a stoichiometric compatibility class
(3) that contains two or more positive steady states". Definition 4.2(1): the multistationarity-allowing region "is the
set of vectors of positive rate constants κ for which (G, κ) is multistationary". By Remark 4.1 the two kinds of region
coincide for one species.

Joshi–Shiu, Theorem 3.6, part 2: G "is multistationary (i.e., cappos(G) ≥ 2) if and only if G has a 2-alternating
subnetwork or the arrow diagram of G has the form" with a two-way arrow at every source.

White, Theorem 1.1:

> "Let G be a one-species reaction network under mass-action kinetics. The multistationarity region M(G) is connected unless G has exactly three distinct source complexes and each of them is mixed. In that case M(G) has exactly two connected components."

## Novelty search

| Channel | Query | Result |
|---|---|---|
| **GitHub code search** | `"multistationarity region" "one-species"`, `multistationarity McClure`, `multistationarity "one-species" connected` | **Found `pw/crnt-connectivity`**, described as resolving "McClure-Shiu arXiv:2303.03960 §6". Other hits: `ehijano/rss_fetch` (arXiv feed dumps), `rick-stevens-ai/Pure-Thought` (a review listing a one-species benchmark network); neither addresses the question |
| Semantic Scholar | citations of arXiv:2303.03960 | Six listed, none answers it: arXiv:2608.18509, 2607.25456, 2512.07560, 2406.11586 (Jiao–Tang–Zeng, J. Math. Biol. 2025), 2403.16556 (Kaihnsa–Telek), 2307.04186 (Kaihnsa–Nguyen–Shiu). Full texts were searched for the citation and for "four or five", "six reactions" and one-species connectedness, not read through; the only hits are about three-species zero-one networks |
| OpenAlex | citers of W4392868789 | Count 2. **The list could not be fetched**: HTTP 429, daily budget exhausted |
| Crossref | the DOI; title search for White's paper | 1 citation recorded; no record of White's paper |
| arXiv API | eight phrase, title and author queries (see the script) | Nothing relevant beyond the source; White's paper is not listed. Closest: arXiv:2607.25456 and 2403.16556, on other networks |
| arXiv:2407.01760, 2310.05466 | word search of the full text, as for the citers | No mention of the question |
| Zenodo API | four queries, plus `"one-species reaction networks"` | Nothing on the question. `10.5281/zenodo.22858244` concerns a phosphorylation cascade |
| `mmislan` | names and descriptions of all 27 repositories, keyword filter | Three reaction-network hits; none on this question |
| `ipitchford` | names and descriptions of all 204 repositories, keyword filter | `gk-shared-phosphatase-connectedness` concerns a cascade; none on this question |
| `gh search repos` | four key phrases | Nothing. **Repository search misses White's repository; only code search finds it** |
| Web search | six queries, including the exact title of White's paper | Did not surface White's paper |

Not reached: Telek's PhD thesis (`noter.math.ku.dk/phd24mlt.pdf`, download failed); an AMS abstract
(`community.ams.org/amsmtgs/2237_abstracts/1117-92-271.pdf`, HTTP 403), which appeared in two searches; Google Scholar. A dissertation by the first author was searched for and not found. G. Voitiuk's 2022 West Virginia thesis came
up in a search; only its catalogue page was read, and that page does not mention this question.

The raw API output is in `certificates/novelty-api.json`. Local copies of the sources are in `certificates/sources/`, which
git ignores.
