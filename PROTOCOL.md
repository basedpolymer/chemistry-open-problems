# Protocol for solving and verifying an open problem

Every problem in this repository goes through the same pipeline. A result is only called **solved** once an
independent verifier, who did not see the solver's reasoning, has re-derived or re-checked it from the original
statement.

## 1. Novelty check (at most half a day)

1. Read the primary source (arXiv PDF or publisher PDF) and copy, **verbatim**, the sentence that makes the
   problem open ("we leave open", "we conjecture", "it remains an open question"), with page/section.
2. Check that it is still open: citing works 2024–2026 (Semantic Scholar / OpenAlex / Google Scholar / arXiv
   listing), the authors' later papers, and AI-driven campaigns:
   - `gh api "users/mmislan/repos?per_page=100"` (Palomar project: CRNT and autocatalytic-set problems);
   - google-deepmind/formal-conjectures (graph theory);
   - a web search on the problem's key phrases.
3. If the problem turns out to be already solved, that is a result: record who solved it, where, and stop.

## 2. Attack

- Prefer exact arithmetic (rationals, algebraic numbers, resultants, Gröbner bases), interval arithmetic
  (python-flint `arb`, mpmath `iv`, IntervalArithmetic.jl) and certified solvers over floating point.
- Floating-point exploration is fine for finding the answer; it is never the certificate.

## 3. Deliverables, in `problems/<slug>/`

| File | Content |
|---|---|
| `README.md` | Problem statement (verbatim quote + full reference), status (solved / partially solved / already solved / failed), one-paragraph result, **attribution** (see below), how to reproduce, list of files. |
| `paper.tex` + `references.bib` | The manuscript, written directly in LaTeX in English with the shared preamble (`\input{chemop-preamble}`, see `templates/paper-skeleton.tex`): precise theorem statement, complete proof, every computation needed to check it, limitations; readable by a referee without the scripts' source code. Compiled to `paper.pdf` by `python templates/build_tex.py problems/<slug>/paper.tex`. |
| `scripts/` | Self-contained, deterministic scripts. Each certificate script prints `PASS` or `FAIL` for each claim and exits non-zero on failure. State the runtime and the environment at the top of each script. |
| `certificates/` | Data that the scripts check (exact parameters, rational vectors, solver logs, interval enclosures). |
| `references.md` | Every source used, with DOI/arXiv id, and what was read (full text / abstract only). |

Scripts run under WSL Fedora (`outils/wsl "<command>"` from Git Bash at the project root runs a command in WSL
with the scientific Python environment activated). Python environment: numpy, scipy, sympy, mpmath,
python-flint, networkx, pyomo, pyscipopt, highspy, cvxpy, ortools, numba. Julia 1.13 with
HomotopyContinuation, IntervalArithmetic, IntervalRootFinding, Catalyst, BifurcationKit, JuMP. Also PARI/GP,
Singular, Macaulay2, msolve, nauty, plantri, buckygen, Lean 4 (no Mathlib yet), CUDA 12.8.

**Attribution.** Every piece of work is credited to the model that produced it, with its reasoning effort: solver,
helpers, verifier, editor (for example "Solver: Claude Opus 5.5 (effort xhigh)", "Solver: GPT-6-Astra via Codex
CLI (effort xhigh)"). The credits appear in `\PaperCredits` of `paper.tex` (printed under the title), at the top of
the problem's `README.md`, and in the results table of the root `README.md`. Unknown efforts are written "effort not
recorded", never guessed.

## 4. Independent verification

A separate agent with a fresh context receives only: the original statement and source, the claimed result, and
the deliverables. It re-derives the key steps itself, re-runs every certificate, looks for gaps, and writes
`verification.md` (verdict: CONFIRMED / CONFIRMED WITH CORRECTIONS / REFUTED / INCONCLUSIVE). Corrections are
applied by the verifier, never by the original solver.

## 5. Vocabulary and honesty rules

- **solved**: complete proof or certificate, confirmed by an independent verifier.
- **solved (preliminary)**: complete argument, not yet independently verified.
- **partial result**: e.g. a certified answer on a parameter range, or a special case.
- **record**: best known value, not proven optimal.
- Never present numerical evidence as a proof. State every assumption and every non-rigorous step explicitly.
- Always quote the exact sentence that makes the problem open.
- No external action (email to authors, arXiv, public GitHub, OEIS, databases) without the repository owner's explicit approval.
