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
- **GPT-6 through Codex** (OpenAI), effort not recorded: Frank–Kasper continuation, exact certificate audit,
  additional analysis and manuscript; no sub-agents used in that continuation.

Independent reviews recorded so far were performed by Claude Opus 5.5 agents with a fresh context. The Frank–Kasper
continuation includes a same-agent audit and still awaits independent review. Each problem's README and the title
page of each manuscript give the detailed credits (solver, helpers, verifier, editor) with the model and effort of
every agent.

> **Status of this collection.** Results have at most been checked by a second AI agent, and those not yet checked by
> an agent say so in the table; none has been refereed or formally checked by humans yet. Statuses are tracked per
> problem (see the table below and [PROTOCOL.md](PROTOCOL.md)). Corrections are welcome and will be recorded in
> [history.md](history.md).

## Results

| # | Problem (source) | Field | Answer | Status |
|---|---|---|---|---|
| 1 | [Sharpness of the (n, n, n−1) equilibrium bound](problems/banaji-feliu-remark-4-23/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 4.23) | Reaction networks | Bound attained for every n; minimum source molecularity is 3 | **solved** |
| 2 | [Limit cycle in a three-reaction example](problems/bbh-three-reaction-no-limit-cycle/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No periodic orbit for any rates; extends to a subfamily | **solved** |
| 3 | [Algebraic limit cycle in a cubic weakly reversible system](problems/craciun-erban-cubic-wr-algebraic-cycle/) (Craciun–Erban, J. Math. Biol. 2025) | Reaction networks | Yes: 10-reaction network with integer rates, W^a(3) ≥ 1 | **solved** |
| 4 | [Critical temperature of the DGH charged di-block polymer](problems/dgh-charged-polymer-diblock/) (Poisat, arXiv:2502.18964, question (3)) | Polymer physics | β_c is the root of Σ_{ℓ≥1} e^{−βℓ²/2} = 1 (value already conjectured by Derrida–Griffiths–Higgs 1992) | **solved** (value predicted in 1992) |
| 5 | [Four limit cycles with three reactions](problems/boros-hofbauer-four-limit-cycles/) (Boros–Hofbauer, J. Dyn. Diff. Equat. 2024, §5.1, conjecture) | Reaction networks | Conjecture true: stable focus surrounded by four limit cycles (two separate computer-assisted proofs) | **solved** |
| 6 | [Jäschke–Skogestad temperature conjecture](problems/jaschke-skogestad-temperature-conjecture/) (Comput. Chem. Eng. 2014, heat-exchanger networks with stream splits) | Chemical engineering | Conjecture false for every branch with ≥ 3 exchangers (exact counterexample); corrected recursion proved for all sizes | **solved** |
| 7 | [Quadratic (n, n+2, n) networks with n+1 equilibria](problems/banaji-feliu-remark-5-8/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 5.8) | Reaction networks | Bound attained for every n (explicit family). Already answered on 4 Oct 2026 by an anonymous unrefereed manuscript (Zenodo 10.5281/zenodo.23146115), on the Phillipson–Rojas (2013) polynomial chain | already solved (rediscovered) |
| 8 | [Hopf criticality in Becker–Döring equations with atomization](problems/becker-doring-atomization-hopf/) (Pego–Velázquez, Nonlinearity 2020) | Nucleation kinetics | First branch: the Hopf bifurcation is **supercritical, with orbitally stable periodic solutions, for every M ≥ 16**, and **subcritical at M = 15** (no Hopf bifurcation for M ≤ 14). Branches 2 (M ≥ 44) and 3 (M ≥ 87): supercritical on the centre manifold, but the periodic orbits are unstable in the full system. Computer-assisted: ball arithmetic for M ≤ 10⁵ (verified), and, for every M ≥ 161 at once, interval boxes in ε = M^{−1/2} down to ε = 0 ([extension](problems/becker-doring-atomization-hopf/all-M/)). Branches j ≥ 4 open | solved (preliminary; the extension to every M is not yet verified by an agent) |
| 9 | [Upper bound for the connective constant of self-avoiding walks on Z³](problems/cubic-saw-upper-bound/) (Pönitz–Tittmann, Electron. J. Combin. 2000: 4.7387, unchanged since) | Polymer physics | μ(Z³) ≤ 4.722048 by an exact certificate; also new bounds for d = 4, 5, 6 | **record** |
| 10 | [Critical parameter of the perturbed Gelfand equation on the disk](problems/perturbed-gelfand-2d-critical-epsilon/) (Korman–Li–Ouyang, J. Differential Equations 2017, conjecture) | Combustion theory (thermal explosion) | Conjecture true: unique ε₀ = 0.2421061655952…, exactly S-shaped below, monotone above (computer-assisted proof) | **solved** |
| 11 | [Discrete critical coagulation–fragmentation](problems/jang-tran-discrete-coagulation-fragmentation/) (Jang–Tran, Adv. Contin. Discrete Models 2025, Conjecture 2) | Aggregation and polymerization kinetics | Conjecture true: for kernels ij and 1, a unique stationary state of every mass 0 < m < 1 (recursion (1.9) positive), and a unique global classical mass-conserving solution for every mass 0 < m ≤ 1, without the second-moment assumption; with Jang–Tran's Theorem 1.1 this gives the critical-mass statement of their Conjecture 1 in the mass-conserving sense | **solved** |
| 12 | [The generalized Petlyuk arrangement minimizes vapor](problems/halvorsen-skogestad-petlyuk-minimum-vapor/) (Halvorsen–Skogestad, Ind. Eng. Chem. Res. 2003, conjecture) | Chemical engineering (distillation) | Conjecture true for every number of components N in the class 𝒜 (any finite network of constant-molar-overflow equilibrium stages, splitters, mixers, recycles, multi-feed columns): every arrangement needs total vaporization ≥ max_j H_j (stage-level lower bound), and the generalized Petlyuk arrangement attains it with finitely many stages (forced binary columns, nondegenerate limit system, exponentially small perturbation). In Underwood's infinite-stage model, the permissive formulation used in the optimization literature admits a non-physical counterexample for N = 5 | **solved** |
| 13 | [Cases 9 and 10 without Andronov–Hopf bifurcation](problems/bbh-cases-9-10-limit-cycles/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No: a case-9 or case-10 network that admits no Andronov–Hopf bifurcation has no periodic orbit for any rate constants (reduction to three families, Liénard form, comparison lemma); so a network of these cases can have a limit cycle only if it admits an Andronov–Hopf bifurcation. Generalizes result 2 | solved (preliminary; not yet verified by an agent) |
| 14 | [Sign of solute skewness in regular-polygon ducts](problems/taylor-dispersion-polygon-skewness/) (Aminian–Camassa–McLaughlin, Stud. Appl. Math. 2018: the triangle "appears to be the only regular polygon" with negative skewness, from finite-element numerics) | Chemical engineering (transport, Taylor dispersion) | Long-time sign of the skewness = sign of ∫ u g₁² for every Péclet number (proved); negative for the equilateral triangle (exact) and positive for every regular n-gon with n ≥ 4; short-time (geometric) skewness positive for every n ≥ 3 (computer-assisted: ball arithmetic for small n, uniform expansion around the disk for large n). "Positive at all times": only certified windows of short and long times (partial) | solved (preliminary) for the long- and short-time signs; not yet verified by an agent |

## Under investigation

Work started on 24 October 2026 and was paused the same day at the owner's request. Frank–Kasper was subsequently
resumed by Codex (GPT-6; effort not recorded), without sub-agents, at the owner's request. Each row gives the state reached;
**fresh independent review is pending for these entries**. Each folder's README and `result.md` state what is proved, what is
certified by computer, what is numerical only, and how to resume. Progress is logged in [history.md](history.md).

| Problem (source) | Field | State reached | Status |
|---|---|---|---|
| [Swelling exponent of self-avoiding walks on Z³](problems/saw-3d-swelling-exponent/) (E_n‖ω_n‖² = n^{2ν+o(1)}; open since Flory, see Madras–Slade 1993, Bauerschmidt–Duminil-Copin–Goodman–Slade 2012) | Polymer physics | No rung passed: the best rigorous bounds remain E_N‖ω_N‖² ≥ N^{4/9}/36 (Madras 2014) and o(n²) (Duminil-Copin–Hammond 2013). Two conditional theorems: Σ_T B_T(z_c)² < ∞ ⇒ E_n‖ω_n‖² = O(n²/(log n)²); a sublinear bound on the multiplicity of Madras's reflection ⇒ an exponent above 4/9. Both reduce to sparse renewal levels of critical bridges (B_T(z_c) → 0, open on Z^d). Exact enumeration to n = 19 refutes several candidate lemmas (e.g. E₁₆[s₁·s₁₆] < 0) | no result on the question; conditional theorems, not yet verified by an agent |
| [Knotting probability of equilateral hexagons](problems/equilateral-hexagon-knot-probability/) (Hake, J. Knot Theory Ramif. 2019: ≤ (14 − 3π)/192 ≈ 0.0238; open per Geisz–Peterson–Shonkwiler 2026) | Polymer physics (ring polymers) | Certified enclosure 3.0675·10⁻⁵ < p₆ < 8.4945·10⁻⁴ (interval branch-and-bound, exact dyadic sums; Monte Carlo 1.364·10⁻⁴); first positive lower bound found; the upper bound is 28 times smaller than Hake's. Hake's Lemmas 3.9 and 3.12 are false (certified knotted counterexamples), so his proof is incomplete, though his bound holds | record (preliminary; not yet verified by an agent) |
| [Minimum number of matches in heat-exchanger networks](problems/hen-minimum-matches/) (Letsios–Kouyialis–Misener, Comput. Chem. Eng. 2018) | Chemical engineering | Closed exactly (block-partition lower bound + exactly checked networks): 37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22 (published upper bound 23). Lower bound 32 (was 31) for balanced15 and balanced15_random0–2. unbalanced10 = 25 by exact SCIP (24 matches infeasible), VIPR certificate not checked. Remaining instances open | partial result and records (preliminary; not yet verified by an agent) |
| [Frank–Kasper tilings in the diblock foam model](problems/frank-kasper-dfm-screening/) (Reddy et al., PNAS 2018; catalogues 2010/2012; project-framed question) | Polymer physics (block copolymers) | Exact embeddedness certificates for all 121 catalogue inputs; six energy witnesses audited with exact box transport, including five improvements over the published SI values. Global enclosure 1 ≤ inf 𝓕 ≤ 1.068575909. σ has the lowest reported numerical value in the 124-entry screen; eight entries remain unrelaxed. A restricted flat A15 family has certified 𝓕 > 1.0699. [Manuscript](problems/frank-kasper-dfm-screening/paper.pdf), [audit](problems/frank-kasper-dfm-screening/verification.md), [progress journal](problems/frank-kasper-dfm-screening/JOURNAL.md). Original screen: Claude Opus 5.5; continuation: Codex GPT-6 (efforts not recorded). | partial result; global ordering unresolved; fresh independent review pending |
| [Kirchhoff analogue of Šoltés's problem](problems/kirchhoff-soltes-problem/) (Gottwald–Majstorović Ergotić–Došlić, Appl. Math. Comput. 2026, Problem 3: graphs with Kf(G − v) = Kf(G) for every v; only C₅ known) | Mathematical chemistry | No new graph. C₅ is the only solution among all 910 914 360 two-connected graphs with n ≤ 11 (published range n ≤ 9), among graphs of cyclomatic number ≤ 5, and among vertex-transitive graphs with n ≤ 39 and n = 41, 43, 45, 47 (published n ≤ 23). Proved: thread lemma (n ≤ 8(c − 1)), expander and dense bounds, complete multipartite graphs, strongly regular graphs and triangle-free distance-regular graphs (C₅ only) | partial result (preliminary; not yet verified by an agent) |
| [Spectral conjecture on fullerene duals](problems/fullerene-dual-cospectral/) (Bille–Buchstaber–Spodarev, J. Math. Chem. 2021) | Mathematical chemistry | Conjecture 1(a) of arXiv v1 is false (smallest pair C42:6 / C42:22) and was already replaced by the authors in the journal version. Exact scan of all isomers for even n ≤ 138. Published Conjecture 1(b) proved for every n ≡ 20 (mod 24), n ≥ 44 and every n ≡ 2 (mod 24), n ≥ 50 (cap rotation on the (6,0) nanotube), verified for 60 ≤ n ≤ 138; 1(c), 1(d) verified for n ≤ 138; literal 1(a) false (C34:2). The authors' 2025 Conjecture 3 (J. Chem. Inf. Model.) fails for α = β (C44:72 / C44:80) | already solved (arXiv-v1 conjecture); partial result on the published conjecture (preliminary; not yet verified by an agent) |

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
| `verification.md` | Report of the verifying agent (once available) |

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
