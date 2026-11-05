# Sources and bounded novelty audit

Audit: 5 November 2026, Europe/Paris. Solver: GPT-6-Astra (effort xhigh).
**Outcome: partial model audit; no claim to have resolved the intended open problem.**

## Primary sources read

1. Dimitrios Letsios, Radu Baltean-Lugojan, Francesco Ceccon, Miten Mistry,
   Johannes Wiebe, Ruth Misener. *Approximation Algorithms for Process Systems
   Engineering*. [arXiv:1909.12328v1](https://arxiv.org/abs/1909.12328),
   26 September 2019. Journal metadata: *Computers & Chemical Engineering* 132
   (2020), 106599, [doi:10.1016/j.compchemeng.2019.106599](https://doi.org/10.1016/j.compchemeng.2019.106599).
   Downloaded the 35-page arXiv PDF and extracted its text. Read §§6.1--6.4 and
   Appendix A.3 in detail, plus the introduction and concluding future directions.
   Rendered and visually inspected printed p. 20 to verify the <= symbols in (5j,k)
   and the complete constraint list. Open sentence is §6.4, printed p. 22, copied
   verbatim in PROGRESS.md. The journal metadata and publisher preview were read;
   the journal's full equation (5) was **not** obtained. No assumption that its
   display is identical to the preprint is made. NEEDS.md requests that DOI.

2. Dimitrios Letsios, Georgia Kouyialis, Ruth Misener. *Heuristics with Performance
   Guarantees for the Minimum Number of Matches Problem in Heat Recovery Network
   Design*. [arXiv:1709.04688v2](https://arxiv.org/abs/1709.04688), 11 April 2018;
   *Computers & Chemical Engineering* 113 (2018), 57--85.
   Downloaded PDF and extracted text; read problem definitions, complexity discussion,
   and Appendix A, especially Theorem 1 (printed p. 41). The reduction uses a bound
   on the **number of positive matches** to make item heat indivisible. Utility
   minimization does not have that bound. The new rational packing counterexample
   only refutes a direct transfer of that construction; it does not challenge the
   paper's hardness result for minimum matches.

3. Ruth Misener, *Approximation Algorithms for Process Systems Engineering*,
   [CMU EWO seminar slides](https://egon.cheme.cmu.edu/ewo/docs/RuthMisener_2020_03_30.pdf).
   Despite the filename, the displayed slides say Tuesday 30 March **2021**.
   Read the heat-recovery slides (approximately slides 12--31), which distinguish
   sequential utility LPs from minimum-matches MILPs and discuss the packing
   reduction. No resolution of the non-isothermal multistage utility question was
   found in those slides. The seminar does not repair equation (5).

4. Christos N. Markides et al., *Working fluid and system optimisation of organic
   Rankine cycles via computer-aided molecular design: A review*, *Progress in
   Energy and Combustion Science* 107 (2025), 101201,
   [doi:10.1016/j.pecs.2024.101201](https://doi.org/10.1016/j.pecs.2024.101201).
   OpenAlex indexes this as 2024; the downloaded PDF states 2025. Downloaded the
   [Ghent repository PDF](https://backoffice.biblio.ugent.be/download/01JJ24PKC9R6M9GNC1DYH05AB6/01JK643C615F8A3NV7PFB8B43A),
   searched the full extracted text, and read the passage citing Letsios et al.
   (reference 182) and surrounding optimization discussion. It cites the survey
   generally for approximation algorithms, not a multistage utility complexity
   result. This was targeted inspection, not a full mathematical audit of the review.

Downloaded third-party PDFs, extracted text, and the equation screenshot are in
`certificates/sources/`, which the repository ignores. They are not deliverables.

## Citation and author follow-up

[OpenAlex W2981102500](https://openalex.org/W2981102500) returned seven citing works
on this run, rather than the approximately nine mentioned in the brief. The exact
DOI record and `filter=cites:W2981102500` response are saved in `novelty-api.json`.
All seven records were screened; access/read depth is stated here:

| Citing work | DOI | Read / relevance |
|---|---|---|
| Markides et al., Rankine-cycle review (2024 online/2025 issue) | 10.1016/j.pecs.2024.101201 | Full text searched and citation context read, as above |
| Ceccon--Misener, Solving the pooling problem at scale with extensible solver GALINI (2022) | 10.1016/j.compchemeng.2022.107660 | Abstract and metadata; pooling software, not a utility-complexity theorem |
| Approaches and application of heat and water network integration in chemical process system engineering: A review (2022) | 10.1016/j.cep.2022.109263 | Metadata only; full text and abstract unavailable in the saved OpenAlex response |
| Optimization in systems engineering: A review of how data analytics and optimization algorithms are applied (2024) | 10.51594/csitrj.v5i4.1027 | Abstract in API; publisher PDF access failed (403) |
| Approximate and robust bounded job start scheduling for Royal Mail delivery offices (2021) | 10.1007/s10951-021-00678-7 | Abstract in API; scheduling, not heat integration |
| Towards scalable quantum annealing for pooling and blending problems: A methodological proof-of-concept (2025) | 10.1016/j.cherd.2025.08.031 | API abstract and publisher preview; discretized pooling, not a utility-complexity result |
| Marousi--Charitopoulos, Accelerated robust spatial Branch-and-Bound algorithm via multi-parametric programming: Application to pooling problems (2026) | 10.1016/j.compchemeng.2026.109746 | API abstract, publisher abstract/highlights and background excerpts; robust pooling algorithm |

An additional web hit was Bloor et al., *A survey and tutorial of reinforcement
learning methods in Process Systems Engineering*, C&CE 206 (2026), 109515,
[doi:10.1016/j.compchemeng.2025.109515](https://doi.org/10.1016/j.compchemeng.2025.109515).
Read the publisher abstract and tutorial overview, not the full paper. No resolution
was identified. Search hits attached to a 2018 scheduling paper are not counted as
citers merely because a secondary aggregator also lists the 2020 survey.

The [Misener publication list](https://wp.doc.ic.ac.uk/rmisener/publications/) was
read through its displayed 2024--2025 entries and the relevant 2018--2022 entries.
It links the 2020 paper back to arXiv; no corrected model or later resolution was
found. The displayed list ends in 2025 and is not evidence of complete 2026 coverage.
Letsios's indexed CV/search results were screened; attempts at the guessed current
homepage/publications URLs returned 404 and are not counted as a completed author
publication audit.

Semantic Scholar's DOI lookup returned 404; its arXiv lookup returned 429. A broad
OpenAlex title search also returned 429, but the direct DOI record and citation
query succeeded. **The full set of roughly nine citers was not independently
reconstructed, and the metadata-only review remains an audit limitation.**

## GitHub and AI-campaign sweep

Read-only commands, with exact queries/results saved in `certificates/novelty-api.json`:

- `gh search code` and `gh search repos`, limit 100, for
  `"multistage minimum utility"`, `"1909.12328"`, `"minimum utility cost"`, and
  `"heat exchanger" complexity`.
- The first, second, and fourth queries returned no code/repository results.
  The third returned one code hit and five repository hits. The code hit was an
  energy/mass-integration **course syllabus** in
  `chauhanpriyanshu001/Khudkibook-Web.github.io/data/syllabus_cache/2711606.json`;
  its raw file was downloaded and read. It lists minimum utility cost as a topic
  and contains no theorem. Repository names/descriptions were unrelated portfolio,
  road-network, project-list, inventory-software, and blockchain material.
- `gh api "users/ipitchford/repos?per_page=100" --paginate --slurp`: 251 repositories,
  all returned pages screened by name/description for heat, utility, pooling,
  Letsios, and the identifier. Broad heat/pool hits concerned internally heated
  convection, quantum thermal cumulants, and pooled pricing, not this question.
- `gh api "users/mmislan/repos?per_page=100"`: 27 repositories, no relevant name or
  description. Since fewer than 100 were returned, no further page was needed.

This is a public code search plus repository-metadata screen, not an inspection
of every file in every campaign repository. No cloning, committing, messaging,
posting, or uploading was done.

## Zenodo and web queries

Initial Zenodo API requests with size 100 returned 400. A first retry batch using
a browser User-Agent returned 403. These failures are preserved, not converted to
zero hits. Later requests with the default urllib User-Agent, percent-encoded
spaces, and size 10 succeeded; see `zenodo-identifier.json` and `zenodo-topics.json`:

| Query | Total returned by API | Inspected scope |
|---|---:|---|
| `1909.12328` and `"1909.12328"` | 0 | Complete exact-identifier result |
| `"minimum utility cost"` | 0 | Complete exact-phrase result |
| `multistage AND utility` | 79 | First ten titles/metadata only |
| `"heat exchanger" AND complexity` | 50 | First ten titles/metadata only |
| Unquoted multistage utility / minimum utility cost / heat exchanger complexity | 31311 / 161490 / 101780 | First ten only; broad OR-like retrieval, weak evidence |

No returned title/description inspected claimed a resolution. The large broad-query
counts are **not** exhaustive screens. The exact phrase/identifier searches are the
most relevant negative evidence.

Web searches included the exact title with 2024, 2025, 2026, PDF, citation and
correction terms; `"multistage minimum utility"` with complexity, polynomial,
hardness and solved; `"minimum utility" "complexity" heat exchanger`;
`"Letsios" "utility" "NP-hard"`; and site-restricted Zenodo searches. Results
pointed to the source, the older minimum-matches work, author metadata, generic
reviews, and simultaneous HENS optimization. Abstracts for non-isothermal stagewise
HENS and simultaneous synthesis were also screened; their general HENS hardness
claims do not isolate utility minimization. Secondary automated summaries were
used only as discovery leads, never as mathematical evidence.

**Conclusion of the bounded search:** no prior resolution found in inspected
sources. This is not proof of priority or proof the problem remains open. More
fundamentally, the displayed source model is defective, so interpreting an intended
physical variant needs information not supplied by the preprint.

## Files and reproduction

- `PROGRESS.md` (git-ignored): source quote, complete limited proofs, model audit,
  missing work, and next idea. No `paper.tex`, `paper.pdf`, or solved README was made.
- `NEEDS.md`: journal DOI and precise model question.
- `scripts/verify_model_audit.py`: deterministic exact arithmetic; all 12 reported
  claims pass. Finite examples are corroboration, not the general proof.
- `certificates/exact-checks.json`: full rational witnesses, LP vertices, PASS results.
- `scripts/novelty_audit.py`, `novelty_followup.py`, `novelty_zenodo_topics.py`:
  read-only, network-dependent audit scripts. Request failures print FAIL and exit
  nonzero; this is separate from mathematical verification.
- `scripts/inspect_audit.py`: concise display of the saved API audit.
- `certificates/novelty-api.json`, `novelty-followup.json`, `zenodo-identifier.json`,
  `zenodo-topics.json`: successful responses and access failures.

The source-to-equation comparison was done by the solver, not an independent
verifier. No independent verification is claimed.
