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

| # | Problem (source) | Field | Answer | Status |
|---|---|---|---|---|
| 1 | [Sharpness of the (n, n, n−1) equilibrium bound](problems/banaji-feliu-remark-4-23/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 4.23) | Reaction networks | Bound attained for every n; minimum source molecularity is 3 | solved (preliminary) |
| 2 | [Limit cycle in a three-reaction example](problems/bbh-three-reaction-no-limit-cycle/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No periodic orbit for any rates; extends to a subfamily | solved (preliminary) |
| 3 | [Algebraic limit cycle in a cubic weakly reversible system](problems/craciun-erban-cubic-wr-algebraic-cycle/) (Craciun–Erban, J. Math. Biol. 2025) | Reaction networks | Yes: 10-reaction network with integer rates, W^a(3) ≥ 1 | solved (preliminary) |
| 4 | [Critical temperature of the DGH charged di-block polymer](problems/dgh-charged-polymer-diblock/) (Poisat, arXiv:2502.18964, question (3)) | Polymer physics | β_c is the root of Σ_{ℓ≥1} e^{−βℓ²/2} = 1 (value already conjectured by Derrida–Griffiths–Higgs 1992) | solved (preliminary) |

Problems under investigation are listed in [history.md](history.md).

## Layout

Each problem lives in `problems/<slug>/`:

| File | Content |
|---|---|
| `README.md` | Verbatim open statement and source, result, status, how to reproduce |
| `notes.md` | Full technical write-up (theorem and proof) |
| `paper.pdf` | The same write-up typeset with the TeXtured template (`templates/build_pdf.py`) |
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
pyscipopt/HiGHS. PDFs are built with pandoc 3 and pdflatex:

```bash
python templates/build_pdf.py problems/<slug>/notes.md
```

## License

Manuscripts and scripts: to be decided before the repository is made public. The PDF template is derived from
[TeXtured](https://github.com/jdujava/TeXtured) (CC0).
