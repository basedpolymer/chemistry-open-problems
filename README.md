# Open problems in chemistry: manuscripts and certificates

This repository collects solutions, produced with AI agents, to **problems explicitly left open in recent
chemistry papers**: questions stated as "we leave open", "we conjecture" or "it remains an open question" in
reaction-network theory, polymer physics and chemical engineering. It follows the model of OpenAI's
[mathematics collection](https://github.com/openai/math): each result comes with the verbatim statement it answers, a
manuscript, and scripts that re-check every computational claim.

## Models

- **Claude Opus 5.5** (Anthropic), effort *xhigh*: orchestration, solvers, verifiers and editors.
- **GPT-6-Astra** (OpenAI), effort *xhigh*: solver of results 4 and 10, contributor to results 11 and 12 and to the
  Frank–Kasper problem.

> **Status of this collection.** Results have at most been checked by a second AI agent with a fresh context; none
> has been refereed by humans. Statuses follow [PROTOCOL.md](PROTOCOL.md); corrections are recorded in
> [history.md](history.md).

## Results

| # | Problem (source) | Field | Answer | Status |
|---|---|---|---|---|
| 1 | [Sharpness of the (n, n, n−1) equilibrium bound](problems/banaji-feliu-remark-4-23/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 4.23) | Reaction networks | Bound attained for every n; minimum source molecularity is 3 | [**solved**](problems/banaji-feliu-remark-4-23/paper.pdf), 18 Oct 2026 |
| 2 | [Limit cycle in a three-reaction example](problems/bbh-three-reaction-no-limit-cycle/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No periodic orbit for any rates; extends to a subfamily | [**solved**](problems/bbh-three-reaction-no-limit-cycle/paper.pdf), 18 Oct 2026 |
| 3 | [Algebraic limit cycle in a cubic weakly reversible system](problems/craciun-erban-cubic-wr-algebraic-cycle/) (Craciun–Erban, J. Math. Biol. 2025) | Reaction networks | Yes: 10-reaction network with integer rates, W^a(3) ≥ 1 | [**solved**](problems/craciun-erban-cubic-wr-algebraic-cycle/paper.pdf), 18 Oct 2026 |
| 4 | [Critical temperature of the DGH charged di-block polymer](problems/dgh-charged-polymer-diblock/) (Poisat, arXiv:2502.18964, question (3)) | Polymer physics | β_c is the root of Σ_{ℓ≥1} e^{−βℓ²/2} = 1 (value predicted by Derrida–Griffiths–Higgs 1992) | [**solved**](problems/dgh-charged-polymer-diblock/paper.pdf), 18 Oct 2026 |
| 5 | [Four limit cycles with three reactions](problems/boros-hofbauer-four-limit-cycles/) (Boros–Hofbauer, J. Dyn. Diff. Equat. 2024, §5.1, conjecture) | Reaction networks | Conjecture true: stable focus surrounded by four limit cycles | [**solved**](problems/boros-hofbauer-four-limit-cycles/paper.pdf), 18 Oct 2026 |
| 6 | [Jäschke–Skogestad temperature conjecture](problems/jaschke-skogestad-temperature-conjecture/) (Comput. Chem. Eng. 2014, heat-exchanger networks with stream splits) | Chemical engineering | Conjecture false for every branch with ≥ 3 exchangers; corrected recursion proved for all sizes | [**solved**](problems/jaschke-skogestad-temperature-conjecture/paper.pdf), 18 Oct 2026 |
| 7 | [Quadratic (n, n+2, n) networks with n+1 equilibria](problems/banaji-feliu-remark-5-8/) (Banaji–Feliu, J. Math. Biol. 2026, Rem. 5.8) | Reaction networks | Bound attained for every n; already answered on 4 Oct 2026 (Zenodo 10.5281/zenodo.23146115) | [already solved (rediscovered)](problems/banaji-feliu-remark-5-8/paper.pdf), 18 Oct 2026 |
| 8 | [Hopf criticality in Becker–Döring equations with atomization](problems/becker-doring-atomization-hopf/) (Pego–Velázquez, Nonlinearity 2020) | Nucleation kinetics | First branch supercritical with stable periodic orbits for every M ≥ 16, subcritical at M = 15; branches 2 and 3 supercritical on the centre manifold, orbits unstable; branches j ≥ 4 open | [solved (preliminary)](problems/becker-doring-atomization-hopf/paper.pdf): M ≤ 10⁵ verified, 18 Oct 2026; [every M](problems/becker-doring-atomization-hopf/all-M/result.md) not yet verified, 24 Oct 2026 |
| 9 | [Upper bound for the connective constant of self-avoiding walks on Z³](problems/cubic-saw-upper-bound/) (Pönitz–Tittmann, Electron. J. Combin. 2000: 4.7387) | Polymer physics | μ(Z³) ≤ 4.722048 by an exact certificate; new bounds for d = 4, 5, 6 | [**record**](problems/cubic-saw-upper-bound/paper.pdf), 18 Oct 2026 |
| 10 | [Critical parameter of the perturbed Gelfand equation on the disk](problems/perturbed-gelfand-2d-critical-epsilon/) (Korman–Li–Ouyang, J. Differential Equations 2017, conjecture) | Combustion theory | Conjecture true: unique ε₀ = 0.2421061655952…, exactly S-shaped below, monotone above | [**solved**](problems/perturbed-gelfand-2d-critical-epsilon/paper.pdf), 18 Oct 2026 |
| 11 | [Discrete critical coagulation–fragmentation](problems/jang-tran-discrete-coagulation-fragmentation/) (Jang–Tran, Adv. Contin. Discrete Models 2025, Conjecture 2) | Aggregation kinetics | Conjecture true for kernels ij and 1: unique stationary state for every mass 0 < m < 1, unique global mass-conserving solution for 0 < m ≤ 1, without second moment | [**solved**](problems/jang-tran-discrete-coagulation-fragmentation/paper.pdf), 20 Oct 2026 |
| 12 | [The generalized Petlyuk arrangement minimizes vapor](problems/halvorsen-skogestad-petlyuk-minimum-vapor/) (Halvorsen–Skogestad, Ind. Eng. Chem. Res. 2003, conjecture) | Chemical engineering (distillation) | Conjecture true for every number of components, for any network of constant-molar-overflow equilibrium stages: every arrangement needs at least max_j H_j, attained by the Petlyuk arrangement with finitely many stages | [**solved**](problems/halvorsen-skogestad-petlyuk-minimum-vapor/paper.pdf), 20 Oct 2026 |
| 13 | [Cases 9 and 10 without Andronov–Hopf bifurcation](problems/bbh-cases-9-10-limit-cycles/) (Banaji–Boros–Hofbauer, Stud. Appl. Math. 2024, §6) | Reaction networks | No: such a network has no periodic orbit for any rate constants (generalizes result 2) | [**solved**](problems/bbh-cases-9-10-limit-cycles/paper.pdf), 24 Oct 2026 |
| 14 | [Sign of solute skewness in regular-polygon ducts](problems/taylor-dispersion-polygon-skewness/) (Aminian–Camassa–McLaughlin, Stud. Appl. Math. 2018) | Chemical engineering (transport) | Long-time skewness negative for the equilateral triangle and positive for every regular n-gon with n ≥ 4; short-time skewness positive for every n ≥ 3; "at all times" only on certified windows | [solved (preliminary)](problems/taylor-dispersion-polygon-skewness/result.md) for long and short times, 24 Oct 2026; not yet verified |
| 15 | [Minimum number of matches in heat-recovery networks](problems/hen-minimum-matches/) (Letsios–Kouyialis–Misener, Comput. Chem. Eng. 2018, §9.2) | Chemical engineering (heat integration) | Four benchmark instances closed: 37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22, unbalanced10 = 25; lower bound 32 for the balanced15 family | [solved (preliminary)](problems/hen-minimum-matches/paper.pdf) for four instances, 24–25 Oct 2026; not yet verified |
| 16 | [Bimolecular (2, 4, 2) networks with two stable equilibria](problems/bbh-242-bimolecular-multistability/) (Banaji–Boros–Hofbauer, Nonlinear Dyn. 2024, §7) | Reaction networks | No further networks: the ten found by the authors are the only ones, boundary equilibria included; at most two equilibria are stable, the origin and a positive one | [solved (preliminary)](problems/bbh-242-bimolecular-multistability/paper.pdf), 25 Oct 2026; not yet verified |

## Under investigation

| Problem (source) | Field | State reached | Status |
|---|---|---|---|
| [Swelling exponent of self-avoiding walks on Z³](problems/saw-3d-swelling-exponent/) (open since Flory; Madras–Slade 1993) | Polymer physics | Known bounds not improved; two conditional theorems | [conditional results](problems/saw-3d-swelling-exponent/upper-bound/result.md), 24 Oct 2026; not yet verified |
| [Knotting probability of equilateral hexagons](problems/equilateral-hexagon-knot-probability/) (Hake, J. Knot Theory Ramif. 2019: p₆ ≤ 0.0238) | Polymer physics (ring polymers) | Certified 3.0675·10⁻⁵ < p₆ < 8.4945·10⁻⁴; two lemmas of Hake's proof are false | [record (preliminary)](problems/equilateral-hexagon-knot-probability/result.md), 24 Oct 2026; not yet verified |
| [Frank–Kasper tilings in the diblock foam model](problems/frank-kasper-dfm-screening/) (Reddy et al., PNAS 2018) | Polymer physics (block copolymers) | σ lowest numerically; certified enclosure 1.00000025 ≤ inf 𝓕 ≤ 1.068575909; global ordering unresolved | [partial result](problems/frank-kasper-dfm-screening/paper.pdf), 24–25 Oct 2026; not yet verified |
| [Kirchhoff analogue of Šoltés's problem](problems/kirchhoff-soltes-problem/) (Gottwald–Majstorović Ergotić–Došlić, Appl. Math. Comput. 2026, Problem 3) | Mathematical chemistry | C₅ still the only solution: exhaustive for n ≤ 11 and for vertex-transitive graphs with n ≤ 39; several classes excluded | [partial result (preliminary)](problems/kirchhoff-soltes-problem/paper.pdf), 24 Oct 2026; not yet verified |
| [Spectral conjecture on fullerene duals](problems/fullerene-dual-cospectral/) (Bille–Buchstaber–Spodarev, J. Math. Chem. 2021) | Mathematical chemistry | arXiv-v1 Conjecture 1(a) false (already corrected by the authors); published 1(b) proved for two residue classes mod 24; 1(c), 1(d) verified for n ≤ 138 | [partial result (preliminary)](problems/fullerene-dual-cospectral/result.md), 24 Oct 2026; not yet verified |
| Infinite families of β-Kirchhoff Šoltés graphs (Gottwald–Majstorović Ergotić–Došlić 2026, Problem 4) | Mathematical chemistry | The source gives families for every β < 2/3 (its Theorem 15); open for β ≥ 2/3. Started 25 Oct 2026 | open |
| Duality between siphons and preclusters (Iosif, [arXiv:2603.08767](https://arxiv.org/abs/2603.08767), 2026, conjecture) | Reaction networks | Started 25 Oct 2026 | open |
| Upper bounds for the connective constants of the diamond and face-centred cubic lattices (follow-up of result 9) | Polymer physics | Started 25 Oct 2026 | open |
| Becker–Döring with atomization, Hopf branches j ≥ 4 (Pego–Velázquez, Nonlinearity 2020, §3.3; follow-up of result 8) | Nucleation kinetics | Started 25 Oct 2026 | open |

## Candidate problems

Shortlisted, not yet started; statements to be checked against the sources before work starts.

| Problem (source) | Field | What is open |
|---|---|---|
| [Uniqueness of the Brownian coagulation similarity profile](problems/brownian-coagulation-profile-uniqueness/) (Smoluchowski equation, full Brownian kernel) | Aggregation kinetics | Uniqueness at fixed mass of the self-similar profile for K(x,y) = 2 + (x/y)^(1/3) + (y/x)^(1/3) |
| Smallest networks with absolute concentration robustness and multistationarity (Kaihnsa–Nguyen–Shiu, Eur. J. Appl. Math. 2024, [arXiv:2307.04186](https://arxiv.org/abs/2307.04186), Question 1.1) | Reaction networks | Which small networks, not necessarily reversible, have both properties |
| Knotting probability of equilateral heptagons (Geisz–Peterson–Shonkwiler, [arXiv:2609.24860](https://arxiv.org/abs/2609.24860), §5) | Polymer physics (ring polymers) | No rigorous bound on p₇ known |
| Energy penalty of the Kaibel column against the Petlyuk arrangement | Chemical engineering (distillation) | Whether the gap is strictly positive for every four-component feed, and a closed form (literature to be checked) |
| Conjectures 1(c) and 1(d) on fullerene duals for every n (Bille–Buchstaber–Spodarev 2021) | Mathematical chemistry | Verified for n ≤ 138 only |

## Layout

Each problem lives in `problems/<slug>/`:

| File | Content |
|---|---|
| `README.md` | Verbatim open statement and source, result, status, credits, how to reproduce |
| `paper.tex`, `references.bib`, `paper.pdf` | The manuscript (preamble `templates/textured-cr/chemop-preamble.tex`) |
| `scripts/` | Deterministic scripts; each prints `PASS`/`FAIL` for every claim |
| `certificates/` | Exact data checked by the scripts, logs |
| `references.md` | Sources read and novelty check |
| `verification.md` | Report of the verifying agent |

Every problem follows [PROTOCOL.md](PROTOCOL.md): novelty check, exact or interval arithmetic, deliverables, then
verification by a second agent with a fresh context.

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
