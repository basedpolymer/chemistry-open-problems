# References and novelty audit

Audit date: **2026-10-20**. Scope: the discrete system with `a(i,j)=ij`,
`b(i,j)=1`, specifically Jang--Tran Conjecture 2. No earlier solution was
identified in the sources and searches below. This is a bounded literature
audit, not a guarantee against unindexed or unpublished prior work.

Raw third-party PDFs, extracted text and API responses are in
`certificates/sources/`, which is git-ignored. No third-party full text is
intended for publication with this repository.

## Primary mathematical sources

1. **Jiwoong Jang and Hung Vinh Tran**, *Discrete coagulation-fragmentation
   equations with multiplicative coagulation kernel and constant fragmentation
   kernel*, Advances in Continuous and Discrete Models **2025**, article 86,
   published 16 April 2025.
   [Journal / DOI](https://doi.org/10.1186/s13662-025-03946-4);
   [arXiv:2409.17974v1](https://arxiv.org/abs/2409.17974v1).
   **Read:** full 25-page arXiv PDF, with page 7 also visually inspected;
   journal HTML for bibliographic metadata and statement. The journal PDF
   endpoint returned an HTML challenge rather than a PDF, so no journal-PDF
   reading is claimed.
   Relevant locations: system (1.1), kernels (1.2), finite-second-moment
   assumption (1.3), weak formulation (1.5), number identity (1.6),
   transform equations (1.7)--(1.8), recursion (1.9), Conjecture 2
   (arXiv printed p. 7, Section 1.4), uniqueness Proposition 1.2,
   small-mass existence Theorem 1.3, transform convergence Lemma 5.1,
   stationary Lemmas 5.5--5.7.
   The quotation in the manuscript/README reproduces the statement supplied
   in the mission brief, checked against the arXiv page.

2. **Truong-Son P. Van**, *A convexity-type invariant for the critical
   coagulation--fragmentation Hamilton--Jacobi equation*,
   [arXiv:2607.02293v2](https://arxiv.org/abs/2607.02293v2),
   submitted 2 July 2026, revised 5 July 2026 (PDF dateline 7 July 2026).
   **Read:** full downloaded PDF, with Sections 1--3 reread for the invariant;
   the parallel agent independently read the full paper.
   This is a **continuous-size** result. In Section 1.2, printed p. 4,
   the exact relevant clause is:
   > “so the discrete analogue of the window 1/2 ≤ m ≤ 1 remains open as of this writing.”
   Thus this source does not pre-empt the discrete conjecture.

3. **Philippe Flajolet and Robert Sedgewick**, *Analytic Combinatorics*,
   Cambridge University Press, 2009,
   [DOI](https://doi.org/10.1017/CBO9780511801655),
   [author-hosted book PDF](https://algo.inria.fr/flajolet/Publications/book.pdf),
   [author's book page](https://sedgewick.io/books/analytic-combinatorics/).
   **Read:** relevant results, not the whole book: Pringsheim's theorem IV.6;
   coefficient transfer VI.1, VI.3--VI.4; parameter-uniform singularity transfer,
   Lemma IX.2, Section IX.7. Used in the (discarded) stationary continuation proof;
   the final manuscript uses only Theorem IV.6, Definition VI.1 (Delta-domain,
   Section VI.3) and Theorems VI.1, VI.3 (see the last section).

## Citing-work databases

Queries were made directly to public APIs. See `scripts/novelty_audit.py`
and the raw responses in `certificates/sources/`.

| Service / exact query | Observed result |
|---|---|
| OpenAlex `works/https://doi.org/10.1186/s13662-025-03946-4` | Work `W4409503603`; one indexed citation. |
| OpenAlex `works?filter=cites:W4409503603&per-page=200` | One work: numerical schemes for a fully nonlinear wave-kinetic coagulation-fragmentation model, DOI `10.1098/rspa.2025.0197`. |
| OpenAlex search `Discrete coagulation fragmentation multiplicative constant Jang Tran` | Five records, including Van 2026, the numerical paper, and inverse/wave-kinetic problems. No resolution identified. |
| Semantic Scholar `paper/ARXIV:2409.17974` with citation fields | HTTP 429; not treated as a negative result. |
| Semantic Scholar `paper/DOI:10.1186/s13662-025-03946-4` with citation fields | Successful; maps DOI and arXiv to the same paper; four citing works listed below. |

Semantic Scholar citing records:

- Van, [2607.02293](https://arxiv.org/abs/2607.02293): full text read;
  continuous model, explicitly leaves discrete window open.
- Le, Tran and Nguyen, [2603.21185](https://arxiv.org/abs/2603.21185),
  DOI [10.1088/1361-6420/ae98ad](https://doi.org/10.1088/1361-6420/ae98ad):
  title/abstract screen; inverse initial-value problem with diffusion,
  not a resolution of this homogeneous discrete existence problem.
- [2509.20747](https://arxiv.org/abs/2509.20747),
  *State-Constrained Chemical Reactions: Discrete-to-Continuous
  Hamilton--Jacobi Equations and Large Deviations*:
  citation/title screen, no claim of solving this conjecture identified.
- [2412.05402](https://arxiv.org/abs/2412.05402),
  *Numerical schemes for a fully nonlinear coagulation-fragmentation model
  coming from wave kinetic theory*: abstract screen, different collision
  model and numerical objective.

These databases have different coverage and dates. Their citation counts
must not be interpreted as a complete bibliography.

## Authors and arXiv 2024--2026

- Read [Jang's research list](https://sites.google.com/view/jiwoongjang/research)
  and [Tran's research list](https://people.math.wisc.edu/~htran24/research.html).
  Both list the discrete paper; the later listed papers do not announce a
  solution. Their current lists include 2026 preprints on other PDE topics.
- arXiv API:
  `all:coagulation AND all:fragmentation AND submittedDate:[202401010000 TO 202610202359]`,
  sorted by submission date, limit 200: **70 records**, all titles screened.
  Van's July 2026 paper is the relevant later critical-window advance.
- Supplementary arXiv author queries `au:"Jang_Jiwoong"` and
  `au:"Tran_Hung"` returned name collisions (23 and 50 records respectively);
  these were not treated as clean or exhaustive author bibliographies.
  Several Van author-query spellings returned empty or unrelated records.
  The arXiv author search page also failed to load through the browser service.
  Van's identified 2026 paper and the broad 70-record subject query were
  checked; no claim of an exhaustive Van author search is made.

Additional potentially relevant records were examined:

| Work | What was read / why it does not pre-empt the claim |
|---|---|
| Bhattacharya--Throm, [2602.10059](https://arxiv.org/abs/2602.10059) | Abstract: perturbations of **constant** coagulation/fragmentation kernels, not the unbounded product kernel. |
| Kerr--Langer, [2504.21665v2](https://arxiv.org/abs/2504.21665v2), DOI [10.1007/s00020-026-02838-w](https://doi.org/10.1007/s00020-026-02838-w) | PDF abstract, introduction and hypotheses in Section 4 read. The coagulation growth assumptions relative to fragmentation do not include `k(n,1)=n` with fragmentation of order `n` and interpolation exponent strictly less than one. No critical mass-window resolution stated. |
| Bouton, [2609.22399](https://arxiv.org/abs/2609.22399) | PDF abstract and theorem/hypothesis screen: critical **additive** bound `a(i,j)<=C(i+j)` with spatial diffusion; excludes `ij`. |
| Matveev, [2606.10785v2](https://arxiv.org/abs/2606.10785v2) | Downloaded review; introduction, aggregation/fragmentation sections and keyword/reference screen. No Jang--Tran resolution located. Full line-by-line reading is not claimed. |
| Folena--Kruszewski, [2607.07959](https://arxiv.org/abs/2607.07959) | Abstract screen: reaction network of assemblies, no stated resolution of the discrete critical-mass conjecture. |

## GitHub and AI campaigns

Read-only commands:

```text
gh api "users/ipitchford/repos?per_page=100" --paginate --slurp
gh api "users/mmislan/repos?per_page=100" --paginate --slurp
gh search repos "coagulation fragmentation" --limit 100 --json fullName,description,url
gh search repos '"Jang" "Tran" coagulation' --limit 100 --json fullName,description,url
gh api "repos/google-deepmind/formal-conjectures/git/trees/main?recursive=1"
```

- **Evidence Press / ipitchford:** 171 public repositories. Names and
  descriptions screened for coagulation, fragmentation, Jang, polymer and
  gelation; no relevant candidate found.
- **Palomar / mmislan:** 24 public repositories. Two polymer-model hits,
  `palomar-p002-binary-polymer-critical-window` and
  `palomar-p012-sparse-raf-exclusion`, concern RAF emergence, not this system.
- Repository search `coagulation fragmentation`: five hits. The relevant-looking
  [linear-polymer-relaxation](https://github.com/yspennstate/linear-polymer-relaxation)
  README was read. It concerns a **finite stochastic chain with pairwise
  constant merging rate `2/n`**, not this product-kernel deterministic equation.
  Other hits were numerical/simulation/control repositories.
- Combined Jang/Tran/coagulation repository search: zero hits.
- `google-deepmind/formal-conjectures`: complete, untruncated tree with
  2,090 paths; filename screen found no relevant conjecture. This was
  **not** a full-text search of every Lean file.

No clone, commit, push, issue, pull request, message, or publication was made.

## Zenodo

Requests used `https://zenodo.org/api/records?q=...`.
Initial larger-page compound queries returned HTTP 400. Retrying with
`size=10` produced usable results:

| Exact `q` | Results |
|---|---|
| `"2409.17974"` | 0 |
| `"10.1186/s13662-025-03946-4"` | 0 |
| `coagulation AND fragmentation AND multiplicative` | 0 |
| `"coagulation-fragmentation"` | 10; all titles screened |
| `coagulation` | 2,659; only first 10, mostly blood/coagulation records |
| `coagulation fragmentation` | 18,675; broad default matching, only first 10 |
| `multiplicative` | 4,333; broad and noisy, only first 10 |

Among the ten exact-phrase hits,
[record 23329462](https://zenodo.org/records/23329462),
*Sharp relaxation, cutoff, and graph diffusion for a linear-polymer chain*,
was inspected through its full API description and the linked GitHub README.
It is the different constant-merging finite-chain model noted above.
The remaining hits concern simulations, control, other stochastic population
balances or other aggregation settings; no solution was identified.
Broad partial result pages are explicitly not claimed to be exhaustive.

## Web queries

The following queries were issued and their relevant hits inspected:

```text
"Jang" "Tran" "coagulation" "Conjecture 2"
"2409.17974"
"2607.02293"
"coagulation" "Jang" "stationary" 2026
"coagulation" "fragmentation" "multiplicative" "stationary" arxiv 2025 2026
"Jang" "Tran" "conjecture" "discrete" solved
site:zenodo.org "coagulation" "Jang"
site:github.com "Jang" "coagulation"
"Jiwoong Jang" publications coagulation 2026
"Truong-Son" "Van" publications 2026 discrete
"Jiwoong Jang" mathematician research
"coagulation" "fragmentation" "Conjecture 2" proof
site:zenodo.org "coagulation" "multiplicative"
```

The primary source, Van's continuous advance, unrelated PDE works and a
database entry still listing the discrete question as open were found.
Search-result snippets and third-party profiles were used only to locate
primary sources, not as mathematical evidence.

## Interpretation

The novelty check found **no earlier solution of the discrete Conjecture 2**.
It supports proceeding with the mission, subject to the explicit coverage
limitations above. A new candidate proof is not an independently verified
resolution, and none of the API retrieval checks is a mathematical certificate.

## Additional sources used in the final manuscript (Claude Opus 5.5)

- **Jang–Tran, arXiv:2409.17974v1**: read in full independently (text extracted from the arXiv PDF). Used:
  system (1.1)–(1.2), assumption (1.3), weak form (1.5), number identity (1.6), transform equation (1.7),
  recursion (1.9), Definition 1, Conjectures 1 and 2, Theorem 1.1, Proposition 1.2, Theorem 1.3, Lemma 3.6
  (curvature barrier, restricted to m < 1/2), Lemmas 5.6–5.7, Corollary 5.2 and Remark 1.
- **Van, arXiv:2607.02293v2**: read in full independently. Used: the half-slope invariant (Section 1.4,
  Lemmas 3.1–3.2, equation (3.3)), the improved barrier (Theorem 4.1), and the quotation of Section 1.2, p. 4.
- **H. V. Tran and T.-S. Van**, *Coagulation-fragmentation equations with multiplicative coagulation kernel and
  constant fragmentation kernel*, Comm. Pure Appl. Math. 75 (2022) 1292–1331, doi:10.1002/cpa.21979.
  **Not read**; cited for the continuous result m < 1/2, as described by Jang–Tran and Van.
- **M. Escobedo, S. Mischler, B. Perthame**, *Gelation in coagulation and fragmentation models*, Comm. Math. Phys.
  231 (2002) 157–188, doi:10.1007/s00220-002-0680-9; and **M. Escobedo, Ph. Laurençot, S. Mischler, B. Perthame**,
  *Gelation and mass conservation in coagulation-fragmentation models*, J. Differential Equations 195 (2003)
  143–174, doi:10.1016/S0022-0396(03)00134-7. **Not read**; cited only because Jang–Tran attribute their
  Conjecture 1 to these papers (bibliographic data from Van's and Jang–Tran's reference lists).
- **Flajolet–Sedgewick**: the final manuscript uses only Pringsheim's theorem (IV.6) and the fixed-parameter
  transfer theorems VI.1 and VI.3 (standard statements); the parameter-uniform Lemma IX.2 was needed only by the
  discarded continuation-in-m argument.

No new database or web search was made for the final manuscript; the audit above (2026-10-20) applies.

## Independent verification (Claude Opus 5.5, effort xhigh, 20 October 2026)

- **Jang–Tran, arXiv:2409.17974v1**: full extracted text re-read by the verifier; Conjecture 2 also checked on
  the rendered page image (printed p. 7). Conjecture 1, Definition 1, (1.3), (1.5)–(1.9), Theorem 1.1,
  Proposition 1.2, Theorem 1.3, Lemmas 3.6, 5.1, 5.6, 5.7, Corollary 5.2 and Remark 1 compared with the manuscript.
- **Jang–Tran, journal version** (doi:10.1186/s13662-025-03946-4): the open-access article HTML was retrieved by
  the verifier (the PDF endpoint was not used). Conjecture 1 ("For m > 1, gelation occurs") and Conjecture 2 have
  the same content as in arXiv v1, slightly reworded, both in Section 1; Definition 1 and Theorem 1.1 unchanged.
- **Van, arXiv:2607.02293v2**: abstract, Section 1.2 (quotation, printed p. 4) and Section 3 (Lemmas 3.1–3.2,
  equation (3.3)) re-read; the discrete invariant W reduces to Van's 2M − xM_x when τ is replaced by x.
- **Flajolet–Sedgewick**: Theorem IV.6 (Pringsheim), Definition VI.1 (Δ-domain, Section VI.3, p. 389),
  Theorems VI.1 and VI.3 checked in the book text.
- **Novelty re-check** (read-only, 20 October 2026): OpenAlex citing works of W4409503603 (1 record, the wave-kinetic
  numerical paper); Semantic Scholar citations of the DOI (4 records: Van 2607.02293, 2603.21185, 2509.20747,
  2412.05402, none a resolution); arXiv API queries `all:coagulation AND all:fragmentation` (40 most recent, to
  2609.29842, 24 Sep 2026), `all:multiplicative AND all:coagulation`, `all:Bernstein AND all:coagulation`,
  `all:gelation AND all:discrete`, `all:Smoluchowski AND all:fragmentation`, `"critical mass" AND coagulation`,
  author queries `au:Jang_Jiwoong` and `au:Tran_Hung_V` (most recent, to 2610.23095 and 2610.17142: no
  coagulation paper after 2409.17974); Zenodo `"coagulation-fragmentation"`, `"2409.17974"`,
  `Jang AND Tran AND coagulation`; GitHub `ipitchford` (171 repositories) and `mmislan` (24) names and descriptions
  screened for coagulation/fragmentation/gelation/Jang/Bernstein (no hit); `gh search repos "coagulation
  fragmentation"`; one web search. No resolution of the discrete conjecture was found. `scripts/novelty_audit.py`
  was also re-run in a scratch copy (same retrieval failures as recorded above: journal PDF challenge, Semantic
  Scholar HTTP 429 on the arXiv id, Zenodo HTTP 400 on two compound queries).
