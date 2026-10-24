# Open problems in chemistry: manuscripts and certificates

This repository collects solutions, produced with AI agents, to **problems explicitly left open in recent
chemistry papers**: questions stated as "we leave open", "we conjecture" or "it remains an open question" in
reaction-network theory, polymer physics and chemical engineering. It follows the model of OpenAI's
[mathematics collection](https://github.com/openai/math): each result comes with the verbatim statement it answers, a
manuscript, and scripts that re-check every computational claim.

## Models

All mathematics, code and manuscripts in this collection were produced by AI agents:

- **Claude Opus 5.5** (Anthropic): solvers, helpers, verifiers, editors and the orchestrating agent; reasoning effort
  *xhigh* unless a problem's README says otherwise (*high* for the scouting phase and a few early solvers).
- **GPT-6-Astra** (OpenAI) through the Codex CLI, reasoning effort *xhigh*: solver of results 4 and 10, contributor
  to results 11 and 12.

Every verifier is a Claude Opus 5.5 agent with a fresh context. Each problem's README and the title
page of each manuscript give the detailed credits (solver, helpers, verifier, editor) with the model and effort of
every agent.

> **Status of this collection.** All results are at most *independently verified by another AI agent*; none has been
> refereed or formally checked by humans yet. Statuses are tracked per problem (see the table below and
> [PROTOCOL.md](PROTOCOL.md)). Corrections are welcome and will be recorded in [history.md](history.md).

## Results

| # | Problem (source) | Field | Answer | Status |
|---|---|---|---|---|
| 1 | [Sharpness of the (n, n, n−1) equilibrium bound](problems/banaji-feliu-remark-4-23/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 4.23) | Reaction networks | Bound attained for every n; minimum source molecularity is 3 | **solved** (independently verified, with corrections) |
| 2 | [Limit cycle in a three-reaction example](problems/bbh-three-reaction-no-limit-cycle/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No periodic orbit for any rates; extends to a subfamily | **solved** (independently verified, with corrections) |
| 3 | [Algebraic limit cycle in a cubic weakly reversible system](problems/craciun-erban-cubic-wr-algebraic-cycle/) (Craciun–Erban, J. Math. Biol. 2025) | Reaction networks | Yes: 10-reaction network with integer rates, W^a(3) ≥ 1 | **solved** (independently verified, with corrections) |
| 4 | [Critical temperature of the DGH charged di-block polymer](problems/dgh-charged-polymer-diblock/) (Poisat, arXiv:2502.18964, question (3)) | Polymer physics | β_c is the root of Σ_{ℓ≥1} e^{−βℓ²/2} = 1 (value already conjectured by Derrida–Griffiths–Higgs 1992) | **solved** (independently verified, with corrections; value predicted in 1992) |
| 5 | [Four limit cycles with three reactions](problems/boros-hofbauer-four-limit-cycles/) (Boros–Hofbauer, J. Dyn. Diff. Equat. 2024, §5.1, conjecture) | Reaction networks | Conjecture true: stable focus surrounded by four limit cycles (two independent computer-assisted proofs) | **solved** (independently verified, with corrections) |
| 6 | [Jäschke–Skogestad temperature conjecture](problems/jaschke-skogestad-temperature-conjecture/) (Comput. Chem. Eng. 2014, heat-exchanger networks with stream splits) | Chemical engineering | Conjecture false for every branch with ≥ 3 exchangers (exact counterexample); corrected recursion proved for all sizes | **solved** (independently verified, with corrections) |
| 7 | [Quadratic (n, n+2, n) networks with n+1 equilibria](problems/banaji-feliu-remark-5-8/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 5.8) | Reaction networks | Bound attained for every n (explicit family). Already answered on 4 Oct 2026 by an anonymous unrefereed manuscript (Zenodo 10.5281/zenodo.23146115), on the Phillipson–Rojas (2013) polynomial chain | already solved (independent rediscovery, verified) |
| 8 | [Hopf criticality in Becker–Döring equations with atomization](problems/becker-doring-atomization-hopf/) (Pego–Velázquez, Nonlinearity 2020) | Nucleation kinetics | ℓ₁ < 0 certified for 16 ≤ M ≤ 10⁵ on branches 1–3, **subcritical at M = 15**; stable cycles proved for 16 ≤ M ≤ 160 (spectral condition certified there), conditional beyond; M → ∞ only formal | partial result (independently verified, with corrections) |
| 9 | [Upper bound for the connective constant of self-avoiding walks on Z³](problems/cubic-saw-upper-bound/) (Pönitz–Tittmann, Electron. J. Combin. 2000: 4.7387, unchanged since) | Polymer physics | μ(Z³) ≤ 4.722048 by an exact certificate; also new bounds for d = 4, 5, 6 | **record** (independently verified, with corrections) |
| 10 | [Critical parameter of the perturbed Gelfand equation on the disk](problems/perturbed-gelfand-2d-critical-epsilon/) (Korman–Li–Ouyang, J. Differential Equations 2017, conjecture) | Combustion theory (thermal explosion) | Conjecture true: unique ε₀ = 0.2421061655952…, exactly S-shaped below, monotone above (computer-assisted proof) | **solved** (independently verified, with corrections) |
| 11 | [Discrete critical coagulation–fragmentation](problems/jang-tran-discrete-coagulation-fragmentation/) (Jang–Tran, Adv. Contin. Discrete Models 2025, Conjecture 2) | Aggregation and polymerization kinetics | Conjecture true: for kernels ij and 1, a unique stationary state of every mass 0 < m < 1 (recursion (1.9) positive), and a unique global classical mass-conserving solution for every mass 0 < m ≤ 1, without the second-moment assumption; with Jang–Tran's Theorem 1.1 this gives the critical-mass statement of their Conjecture 1 in the mass-conserving sense | **solved** (independently verified, with corrections) |
| 12 | [The generalized Petlyuk arrangement minimizes vapor](problems/halvorsen-skogestad-petlyuk-minimum-vapor/) (Halvorsen–Skogestad, Ind. Eng. Chem. Res. 2003, conjecture) | Chemical engineering (distillation) | Conjecture true for every number of components N in the class 𝒜 (any finite network of constant-molar-overflow equilibrium stages, splitters, mixers, recycles, multi-feed columns): every arrangement needs total vaporization ≥ max_j H_j (stage-level lower bound), and the generalized Petlyuk arrangement attains it with finitely many stages (forced binary columns, nondegenerate limit system, exponentially small perturbation). In Underwood's infinite-stage model, the permissive formulation used in the optimization literature admits a non-physical counterexample for N = 5 | **solved** (independently verified, with corrections) |
| 13 | [Cases 9 and 10 without Andronov–Hopf bifurcation](problems/bbh-cases-9-10-limit-cycles/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No: a case-9 or case-10 network that admits no Andronov–Hopf bifurcation has no periodic orbit for any rate constants (reduction to three families, Liénard form, comparison lemma); so a network of these cases can have a limit cycle only if it admits an Andronov–Hopf bifurcation. Generalizes result 2 | solved (preliminary; independent verification pending) |

## Under investigation

Work started on 24 October 2026. Everything in this table is **preliminary and unverified**: it records what the
solver agents report in their folders, not checked results. Statuses will follow [PROTOCOL.md](PROTOCOL.md) once an
independent verifier has checked each claim; progress is logged in [history.md](history.md).

| Problem (source) | Field | Current state (unverified) |
|---|---|---|
| [Swelling exponent of self-avoiding walks on Z³](problems/saw-3d-swelling-exponent/) (E_n‖ω_n‖² = n^{2ν+o(1)}; open since Flory, see Madras–Slade 1993, Bauerschmidt–Duminil-Copin–Goodman–Slade 2012) | Polymer physics | First wave done (literature, lower bound, upper bound; computation lab still running). No rung is passed: the best rigorous bounds remain E_N‖ω_N‖² ≥ N^{4/9}/36 (Madras 2014) and o(n²) (Duminil-Copin–Hammond 2013). Two conditional theorems: Σ_T B_T(z_c)² < ∞ ⇒ E_n‖ω_n‖² = O(n²/(log n)²); a sublinear bound on the multiplicity of Madras's reflection ⇒ an exponent above 4/9. Exact enumeration to n = 19 refutes several candidate lemmas (e.g. positivity of all step correlations fails at n = 16) |
| [Becker–Döring with atomization for every M](problems/becker-doring-atomization-hopf/all-M/) (extension of result 8) | Nucleation kinetics | Spectral hypothesis H2 swept for all M ≥ 161 on branch 1 (201 interval boxes in ε = M^{−1/2}, all pass); a formula for ℓ₁ regular as M → ∞ gives ℓ₁ < 0 near ε = 0; combined sweep in progress |
| [Knotting probability of equilateral hexagons](problems/equilateral-hexagon-knot-probability/) (Hake, J. Knot Theory Ramif. 2019: ≤ (14 − 3π)/192 ≈ 0.0238; still open per Geisz–Peterson–Shonkwiler 2026) | Polymer physics (ring polymers) | Interval branch-and-bound running; provisional enclosure [6.9·10⁻⁷, 8.95·10⁻³] (Monte Carlo ≈ 1.35·10⁻⁴) |
| [Minimum number of matches in heat-exchanger networks](problems/hen-minimum-matches/) (Letsios–Kouyialis–Misener, Comput. Chem. Eng. 2018) | Chemical engineering | Block-partition lower bound plus exactly checked solutions: 37sp-yfyv = 36, 20sp1 = 19 and 23sp1 = 22 (published bounds 23 / 17) reported closed; lower bound 32 (was 31) for balanced15 and b15_rand0–2; exact SCIP + VIPR pipeline built |
| [Skewness of Taylor dispersion in regular-polygon ducts](problems/taylor-dispersion-polygon-skewness/) (Aminian–Camassa–McLaughlin, Stud. Appl. Math. 2018: the triangle "appears to be the only regular polygon" with negative long-time skewness, from finite-element numerics) | Chemical engineering (transport) | Long-time sign reduced to the sign of ∫ u g₁²; first interval certificates for 4 ≤ n ≤ 8 being produced; large n by expansion around the disk |
| [Frank–Kasper tilings in the diblock foam model](problems/frank-kasper-dfm-screening/) (Reddy et al., PNAS 2018; catalogue of Dutour Sikirić–Delgado-Friedrichs–Deza 2010) | Polymer physics (block copolymers) | Curved-cell DFM engine written; flat-cell values reproduce the paper; reproducing the published curved-cell values (first milestone) in progress. Not posed as open in the literature: the expected result is the first exhaustive screen of the catalogue |
| [Kirchhoff analogue of Šoltés's problem](problems/kirchhoff-soltes-problem/) (Gottwald–Došlić–Majstorović Ergotić, Appl. Math. Comput. 2026: graphs with Kf(G − v) = Kf(G) for every v; only C₅ known) | Mathematical chemistry | No new graph: only C₅ among 2-connected graphs with n ≤ 10 and vertex-transitive graphs with 10 ≤ n ≤ 31. A thread lemma reduces each cyclomatic number c to finitely many kernels; none for c = 2, 3, 4 |
| [Spectral conjecture on fullerene duals](problems/fullerene-dual-cospectral/) (Bille–Buchstaber–Spodarev, J. Math. Chem. 2021) | Mathematical chemistry | Conjecture 1(a) of arXiv v1 is false (exact cospectral non-isomorphic pairs from n = 42), but the authors already replaced it in the published version, whose Table 1 matches; the published Conjecture 1 (a)–(d) is being tested in exact arithmetic |

## Layout

Each problem lives in `problems/<slug>/`:

| File | Content |
|---|---|
| `README.md` | Verbatim open statement and source, result, status, how to reproduce |
| `paper.tex`, `references.bib` | The manuscript, written directly in LaTeX (TeXtured-based preamble `templates/textured-cr/chemop-preamble.tex`) |
| `paper.pdf` | The compiled manuscript (`python templates/build_tex.py problems/<slug>/paper.tex`) |
| `scripts/` | Deterministic scripts; each prints `PASS`/`FAIL` for every claim |
| `certificates/` | Exact data checked by the scripts, logs |
| `references.md` | Sources read and novelty check |
| `verification.md` | Report of the independent verifier (once available) |

## How the results were produced

An orchestrating agent (Claude Opus 5.5) selected problems from a scouting phase that mined recent papers for
explicitly open statements, then assigned each problem to a solver agent (Claude Opus 5.5 sub-agents, or GPT-6-Astra
through the Codex CLI). Every solver follows [PROTOCOL.md](PROTOCOL.md): novelty check, attack with exact or interval
arithmetic, deliverables, then verification by a separate agent with a fresh context that re-derives the result from
the original statement. The human owner of the repository chose the scope and approves every external
action; no author has been contacted yet.

## Reproducing

Scripts were run under Fedora 44 (WSL) with Python 3.12 and sympy 1.14, mpmath, python-flint 0.9, numpy/scipy,
pyscipopt/HiGHS. PDFs are built with latexmk, pdflatex and biber (MiKTeX or TeX Live):

```bash
python templates/build_tex.py problems/<slug>/paper.tex
```

## License

- **Code** (`scripts/`, `verification/`, `templates/*.py` and the LaTeX preamble code in `templates/`): [MIT](LICENSE).
- **Manuscripts, reports, documentation and certificate data** (`paper.tex`, `paper.pdf`, `*.md`, `certificates/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0). Please cite the repository and the problem folder.
- **Exceptions.** Verbatim quotations of the source papers remain the property of their authors and are reproduced
  for citation only. The PDF template is derived from [TeXtured](https://github.com/jdujava/TeXtured), which is
  CC0 (`templates/textured-cr/LICENSE-TeXtured`).
