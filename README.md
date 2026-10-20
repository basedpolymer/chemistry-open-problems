# Open problems in chemistry: manuscripts and certificates

This repository collects solutions, produced with AI agents, to **problems explicitly left open in recent
chemistry papers**: questions stated as "we leave open", "we conjecture" or "it remains an open question" in
reaction-network theory, polymer physics and chemical engineering. It follows the model of OpenAI's
[mathematics collection](https://github.com/openai/math): each result comes with the verbatim statement it answers, a
manuscript, and scripts that re-check every computational claim.

> **Status of this collection.** All results are at most *independently verified by another AI agent*; none has been
> refereed or formally checked by humans yet. Statuses are tracked per problem (see the table below and
> [PROTOCOL.md](PROTOCOL.md)). Corrections are welcome and will be recorded in [history.md](history.md).

## Results

| # | Problem (source) | Field | Answer | Credits (model, effort) | Status |
|---|---|---|---|---|---|
| 1 | [Sharpness of the (n, n, n−1) equilibrium bound](problems/banaji-feliu-remark-4-23/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 4.23) | Reaction networks | Bound attained for every n; minimum source molecularity is 3 | Claude Opus 5.5 (xhigh); n = 3 case: Claude Opus 5.5 (high); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 2 | [Limit cycle in a three-reaction example](problems/bbh-three-reaction-no-limit-cycle/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No periodic orbit for any rates; extends to a subfamily | Claude Opus 5.5 (high); certificate: Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 3 | [Algebraic limit cycle in a cubic weakly reversible system](problems/craciun-erban-cubic-wr-algebraic-cycle/) (Craciun–Erban, J. Math. Biol. 2025) | Reaction networks | Yes: 10-reaction network with integer rates, W^a(3) ≥ 1 | Claude Opus 5.5 (high, xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 4 | [Critical temperature of the DGH charged di-block polymer](problems/dgh-charged-polymer-diblock/) (Poisat, arXiv:2502.18964, question (3)) | Polymer physics | β_c is the root of Σ_{ℓ≥1} e^{−βℓ²/2} = 1 (value already conjectured by Derrida–Griffiths–Higgs 1992) | GPT-6-Astra via Codex CLI (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections; value predicted in 1992) |
| 5 | [Four limit cycles with three reactions](problems/boros-hofbauer-four-limit-cycles/) (Boros–Hofbauer, J. Dyn. Diff. Equat. 2024, §5.1, conjecture) | Reaction networks | Conjecture true: stable focus surrounded by four limit cycles (two independent computer-assisted proofs) | Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 6 | [Jäschke–Skogestad temperature conjecture](problems/jaschke-skogestad-temperature-conjecture/) (Comput. Chem. Eng. 2014, heat-exchanger networks with stream splits) | Chemical engineering | Conjecture false for every branch with ≥ 3 exchangers (exact counterexample); corrected recursion proved for all sizes | Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 7 | [Quadratic (n, n+2, n) networks with n+1 equilibria](problems/banaji-feliu-remark-5-8/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 5.8) | Reaction networks | Bound attained for every n (explicit family). Already answered on 4 Oct 2026 by an anonymous unrefereed manuscript (Zenodo 10.5281/zenodo.23146115), on the Phillipson–Rojas (2013) polynomial chain | Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | already solved (independent rediscovery, verified) |
| 8 | [Hopf criticality in Becker–Döring equations with atomization](problems/becker-doring-atomization-hopf/) (Pego–Velázquez, Nonlinearity 2020) | Nucleation kinetics | ℓ₁ < 0 certified for 16 ≤ M ≤ 10⁵ on branches 1–3, **subcritical at M = 15**; stable cycles proved for 16 ≤ M ≤ 160 (spectral condition certified there), conditional beyond; M → ∞ only formal | Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | partial result (independently verified, with corrections) |
| 9 | [Upper bound for the connective constant of self-avoiding walks on Z³](problems/cubic-saw-upper-bound/) (Pönitz–Tittmann, Electron. J. Combin. 2000: 4.7387, unchanged since) | Polymer physics | μ(Z³) ≤ 4.722048 by an exact certificate; also new bounds for d = 4, 5, 6 | Claude Opus 5.5 (xhigh); verified by Claude Opus 5.5 (xhigh) | **record** (independently verified, with corrections) |
| 10 | [Critical parameter of the perturbed Gelfand equation on the disk](problems/perturbed-gelfand-2d-critical-epsilon/) (Korman–Li–Ouyang, J. Differential Equations 2017, conjecture) | Combustion theory (thermal explosion) | Conjecture true: unique ε₀ = 0.2421061655952…, exactly S-shaped below, monotone above (computer-assisted proof) | GPT-6-Astra via Codex CLI (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 11 | [Discrete critical coagulation–fragmentation](problems/jang-tran-discrete-coagulation-fragmentation/) (Jang–Tran, Adv. Contin. Discrete Models 2025, Conjecture 2) | Aggregation and polymerization kinetics | Conjecture true: for kernels ij and 1, a unique stationary state of every mass 0 < m < 1 (recursion (1.9) positive), and a unique global classical mass-conserving solution for every mass 0 < m ≤ 1, without the second-moment assumption; with Jang–Tran's Theorem 1.1 this gives the critical-mass statement of their Conjecture 1 in the mass-conserving sense | Claude Opus 5.5 (xhigh); GPT-6-Astra via Codex CLI (xhigh); verified by Claude Opus 5.5 (xhigh) | **solved** (independently verified, with corrections) |
| 12 | [The generalized Petlyuk arrangement minimizes vapor](problems/halvorsen-skogestad-petlyuk-minimum-vapor/) (Halvorsen–Skogestad, Ind. Eng. Chem. Res. 2003, conjecture) | Chemical engineering (distillation) | Conjecture true for every number of components N in a stated class: networks of two-section columns (containing Agrawal's regular-column space), strict infinite-stage Underwood model, no reverse flows; the Petlyuk minimum equals the highest V_min peak (also proved). In the permissive (TG-literal) model: true for N ≤ 4, false for N = 5 (interval-certified, non-physical counterexample) | Claude Opus 5.5 (xhigh): theory, computation, general N (three agents); verified by Claude Opus 5.5 (xhigh) | **solved in the stated class and model** (independently verified, with corrections); open outside it |

Each problem's README and the title page of each manuscript give the detailed credits (solver, helpers, verifier,
editor) with the reasoning effort of every model.

Problems under investigation are listed in [history.md](history.md).

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
