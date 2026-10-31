# Sources and novelty audit

Audit date: 31 October 2026. No prior determination of the exact reactant-complex minimum was found
in the material examined. This is a bounded search, not a claim that all literature has been checked.

## Sources used in the proof

1. N. Kaihnsa, T. Nguyen and A. Shiu, *Absolute concentration robustness and multistationarity in
   reaction networks: conditions for coexistence*, European Journal of Applied Mathematics 35 (2024),
   566-600. DOI [10.1017/S0956792523000335](https://doi.org/10.1017/S0956792523000335),
   [arXiv:2307.04186](https://arxiv.org/abs/2307.04186).
   **Read:** published full HTML, especially definitions, rank bounds, Section 5 and the Section 6
   open direction. The manuscript quotes one sentence from Section 6. The exact-minimum problem here
   refines the source's strict-lower-bound conjecture; it is not described as a verbatim source conjecture.
2. Claude Opus 5.5, *Absolute concentration robustness with multistability in bimolecular networks:
   three species suffice, dimension two does not*, result 27 in this repository, 27 October 2026.
   **Read:** [parent manuscript](../paper.tex), [verification](../verification.md),
   [example report](../examples/report.md), source classification and certificate scripts.
   **Used:** the dimension-two obstruction (already separately checked in that work), the earlier
   classification and the identified numerical-only cases S3, S7, S9. The new certificates regenerate
   the classifications and do not read the numerical search lists.

## Other primary literature examined during discovery

- M. Banaji and E. Feliu, *The number and geometry of equilibria in chemical reaction networks*,
  [DOI 10.1007/s00285-026-02429-8](https://doi.org/10.1007/s00285-026-02429-8),
  [arXiv:2409.06877](https://arxiv.org/abs/2409.06877): selected results on few-reactant
  equilibrium bounds consulted. No theorem from it is needed in the final proof.
- Banaji, Boros and Hofbauer, planar four-reaction classification,
  [DOI 10.1007/s11071-024-10068-1](https://doi.org/10.1007/s11071-024-10068-1): scope inspected
  during exploration. A bound on reactions was not substituted for a bound on reactant complexes.
- *Bistability, Absolute Concentration Robustness, and Hysteresis in Dual-Site Futile Cycles with
  Bifunctional Enzymes*, [arXiv:2605.11389](https://arxiv.org/abs/2605.11389): abstract inspected
  in this follow-up; not used as a proof source. It concerns a different network family.
- Feliu, Henriksson and Pascual-Escudero, *The generic geometry of steady state varieties*,
  [DOI 10.1137/25M1731289](https://doi.org/10.1137/25M1731289): abstract/metadata inspected;
  not a proof source for the minimum established here.

## Public-source audit and limitations

The machine-readable record is [certificates/novelty-api.json](certificates/novelty-api.json), produced
by the read-only `scripts/novelty_audit.py`.

- OpenAlex supplied four citing works. Semantic Scholar supplied eleven citing works. Titles and
  metadata were screened; this was not a full-text review of every citing paper.
- The broad OpenAlex topic query returned 1,166 items; only the first 50 were returned and screened.
  Many were unrelated to chemical reaction networks.
- Public repositories of `mmislan` (27) and `ipitchford` (204) were enumerated. The ACR zero-divisor
  project concerns a different question. Repository metadata screening does not audit every theorem.
- GitHub repository searches for the ACR phrase and for four-reactant wording were performed.
  The `google-deepmind/formal-conjectures` tree was screened for relevant paths.
- The Zenodo API query returned HTTP 400; it did not yield a usable independent corpus search.
- Web searches included ACR with multistability, minimum reactant complexes, four reactant complexes,
  and five reactant complexes. They located the primary source and related work, with no matching
  exact-minimum solution in the retrieved results.

No paywalled PDF was necessary. No files or messages were published or sent to authors. Priority has
not been established beyond these searches.

## Verifier's recheck, 31 October 2026

The arXiv version of the primary source was read in full text by the verifier, for its definitions, Example 4.13,
Proposition 3.10 and the Section 6 quotation. Also screened: the arXiv listing for the phrase "absolute concentration
robustness" (42 results, abstracts of the 40 returned), the 11 citing works from Semantic Scholar and the count of 4
from OpenAlex, a Zenodo query (which worked with the phrase as keywords and gave nothing relevant), the repositories
of `mmislan` and `ipitchford`, and a GitHub repository search. The nearest recent works concern other questions:
arXiv:2604.14592 and arXiv:2608.18509 study zero-one networks of low dimension, and arXiv:2605.11389 one network family.
Screening used titles and abstracts, not every full text. No statement of the exact minimum was found.
