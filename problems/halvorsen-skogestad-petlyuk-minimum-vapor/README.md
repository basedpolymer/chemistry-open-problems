# Halvorsen–Skogestad conjecture: the generalized Petlyuk arrangement minimizes total vaporization

**Attribution.** *Version 3 (2026-10-20, third wave: finite-stage attainment).* Design and limit configuration:
Claude Opus 5.5 (effort xhigh) (Agent D). Limit system and its nondegeneracy: Claude Opus 5.5 (effort xhigh)
(Agent L). Analysis and integration (`attainment/theorem.md`): Claude Opus 5.5 (effort xhigh) (Agent A).
Computation and red team: Claude Opus 5.5 (effort xhigh) (Agent C). Proof architecture (coordinator): Claude Opus
5.5 (effort not recorded). Independent verification of version 3: **pending**.
*Version 2 (2026-10-20, second wave).* Foundations, class 𝒜, red team and integration: Claude Opus
5.5 (effort xhigh) (Agent F). Sections and stage-level lemmas, binary attainment: Claude Opus 5.5 (effort xhigh)
(Agent S). Network argument, device decomposition, tie structure: Claude Opus 5.5 (effort xhigh) (Agent N). Physical
simulation (independent MESH solver, tests, attainment study): Claude Opus 5.5 (effort xhigh) (Agent P). Literature:
Claude Opus 5.5 (effort xhigh). Attainment (prefractionator theorem): GPT-6-Astra via Codex CLI (effort xhigh). Independent verification of version 2: **pending**.
*Version 1 (independently verified).* Theory: Claude Opus 5.5 (effort xhigh) (Agent T). Computation: Claude Opus 5.5
(effort xhigh) (Agent C). General N: Claude Opus 5.5 (effort xhigh) (Agent G). Independent verification: Claude Opus
5.5 (effort xhigh).

**Status (version 3): conjecture solved (preliminary) for every N, in the class 𝒜; independent verification
pending.**
- **Lower bound (version 2, every N).** Every steady state of every arrangement of the class 𝒜 with purities ≥ 1 − ε
  has R ≥ H_j − c_j ε for every j; hence R*(𝒮) ≥ max_j H_j for every family 𝒮 ⊆ 𝒜 (an early independent check
  confirmed it: `total/verification-early/verification-lower-bound.md`).
- **Attainment (version 3, every N ≥ 3; N = 2 in version 2).** For every δ, ε > 0 there is n_0 such that, for all
  section lengths ≥ n_0, the generalized Petlyuk arrangement has an operating specification with R = max_j H_j + δ
  and a steady state with all mole fractions positive and all purities ≥ 1 − ε. Hence R*(𝒮_Petlyuk) = max_j H_j.
  Complete proof: `attainment/theorem.md` (design `attainment/design/design.md`, Theorem L
  `attainment/limit/limit.md`, analysis `attainment/analysis/analysis.md`); summary `attainment/FINAL-REPORT.md`.
- **Conjecture.** No arrangement of 𝒜 needs less total vaporization, in the infinite-stage pure-product limit, than
  the generalized Petlyuk arrangement, whose minimum is exactly HS's value max_j H_j (`attainment/theorem.md`,
  Corollary 2).
- The manuscript `paper.tex` is still version 2 (attainment for N ≥ 3 stated there as open); it will be rewritten
  as a single paper. Version 1's results in the Underwood model are unchanged (independently verified, with
  corrections).

## Problem

I. J. Halvorsen, S. Skogestad, *Minimum energy consumption in multicomponent distillation. 3. More than three
products and generalized Petlyuk arrangements*, Ind. Eng. Chem. Res. 42 (2003) 616–629, doi:10.1021/ie0108651.

Verbatim open statements:
- Abstract: "We also conjecture that the minimum energy requirement for the generalized extended Petlyuk arrangement
  is lower than the minimum-energy requirement for any distillation configuration when we consider conventional
  adiabatic sections and no internal heat exchange."
- §6.3 (p. 626): "We assume constant relative volatilities, constant molar flows, constant pressure, and no internal
  heat integration. Then the generalized adiabatic Petlyuk arrangement has the lowest need for vaporization compared
  to any other adiabatic distillation arrangement for separation of an arbitrary feed mixture into its pure
  components. We have not presented a complete proof, so the above conclusion is a conjecture."
- §5.3 (p. 624), on the Petlyuk arrangement itself: "We have not carried out a detailed proof for the general
  N-component M-product case for other than the far left and right peaks."

Setting: ideal N-component feed (constant relative volatilities α₁ > … > α_N, composition z, liquid fraction q),
constant molar overflow, constant pressure, infinitely many stages, sharp separation into pure products; objective =
total vaporization (heat input). With θ_j the roots of Σ α_p z_p/(α_p − θ) = 1 − q, the Petlyuk value is
V_Petlyuk = max_j H_j, H_j = Σ_{p>j} α_p z_p/(θ_j − α_p) (the highest peak of the V_min diagram, referred to the
reboiler).

## Result (version 3: finite-stage attainment)

Idea: start from HS's preferred-split point (paper.tex Theorem `achieve`), give every column with ≥ 3 components a
small slack in GPT-6-Astra's direction (which keeps its internal common Underwood roots), let the margin δ enter
through the reboiler only, and **force** every binary column: a small net flow η_p of one key is sent into the
adjacent side product. Forcing pins the binary junction compositions (without it the infinite-stage limit has a
free composition front, the neutral mode of `total/network/network.md` §10.4). In the infinite-stage limit every
section is then singly pinched, and the limit system, in junction compositions and column splits, has an invertible
Jacobian: block triangular, with a top-down block for the major variables (Underwood interpolation and a Cauchy-type
matrix per column) and one triangular M-matrix per trace component (Theorem L). The finite-stage equations differ
from the limit system by terms of size ρ^{n+1} in C¹, so a quantitative inverse-function step gives a unique
steady state near the limit for all section lengths ≥ n_0; its stage profiles have S(x) > 0 everywhere, and the
M-matrix structure of the component balances makes every mole fraction positive. Limits: n → ∞, then η → 0, then
δ → 0. Numerical evidence (not used in the proof): 61 of 62 cases pass for N = 3..8 with an independent solver at up
to about 1500 digits (`attainment/computation/REPORT.md`); the failing case lies outside the design's parameter
region on purpose and shows that this region is needed.

## Result (version 2)

**Class 𝒜** (`total/foundations/class.md`): every finite network of constant-molar-overflow equilibrium stages,
single-phase splitters and mixers, total condensers and total vaporizers, with one process feed (split, preheated,
precooled or pre-flashed at will), heat duties at any stage, recycles, and products of either phase; R = total heat
supplied. "Infinitely many stages and pure products" is the limit R*(𝒮) = lim_{ε→0} inf{R : steady state with purity
≥ 1 − ε} over a family 𝒮 of arrangements (all stage numbers of a configuration, or all of 𝒜). Every device of the
arrangements HS consider (columns, thermal couplings, side draws, dividing walls, Kaibel columns, several feed stages,
stream splitters, remixing, bypasses, intermediate exchangers) and Murphree trays decompose into these units
(`total/network/network.md` §1.2); 𝒜 even contains non-adiabatic sections with external utilities.

**Theorem 1 (lower bound, complete proof).** At every steady state of every arrangement of 𝒜, with any numbers of
stages and any purities, and for every j,

  R − Σ_P max(0, φ_P(θ_j)) = (heat supplied to light-side units) + (heat removed from heavy-side units)
                              + Σ_{streams between the two sides} |φ_e(θ_j)| ≥ 0,

where φ_e(θ) = v_e − Σ_i α_i n_{e,i}/(α_i − θ) is the Underwood flux of a stream (v_e vapor flow, n_e component flows)
and a unit is heavy-side if its outflows have φ ≥ 0. With purity ≥ 1 − ε: R ≥ H_j − c_j ε,
c_j = max(α_j/(α_j − θ_j), α_{j+1}/(θ_j − α_{j+1})). Hence every arrangement needs R*(𝒮) ≥ max_j H_j.

*Proof idea.* (1) CMO and component balances give, at every unit, Σ_out φ − Σ_in φ = Q_u (+ feed term, zero at θ_j).
(2) A liquid of composition x has φ = −l G_x(θ) and the vapor in equilibrium with it φ = −(θ v/S(x)) G_x(θ),
G_x(θ) = Σ_i α_i x_i/(α_i − θ): all streams leaving a stage have the same sign of flux (single-phase units trivially).
(3) Summing (1) over the heavy-side units gives the identity. No section physics, no Underwood roots of sections, no
assumption on topology, net-flow signs, distributions or numbers of stages. Found independently by three agents
(N, S, F); complete proofs in `class.md` §5, `total/sections/sections.md`, `total/network/network.md`; checked by four
independent solvers and in exact rational arithmetic.

**Consequences.** The excess over the bound is three nonnegative "pinch losses" (heating on the light side, cooling on
the heavy side, flux across the hyperplane G_x(θ_j) = 0 through the feed-flash liquid); approaching max_j H_j forces
Underwood's active-root pinches; all gaps of version 1 are closed for the lower bound (degenerate points, reverse
flows, several feed stages, splitters, recycles); the version-1 permissive counterexample cannot be approached by
physical operation (it has a light-to-heavy sign path C → CDE → BCDE → BCD → D, impossible at a steady state); with
heat exchange between process streams the external heat can drop below max_j H_j only by moving heat from light-side
to heavy-side units, which is why HS exclude it.

**Attainment.** By Theorem 1 the conjecture follows from R*(𝒮_Petlyuk) = max_j H_j (HS's formula for the Petlyuk
minimum, in the finite-stage-limit sense). In Underwood's infinite-stage model the FTC attains it for
every N (version 1, Theorem achieve). Finite stages: **N = 2 proved** (S, `total/sections/sections.md` Prop. 10);
steady states exist for any prescribed flows (N, `total/network/network.md` Theorem 10.1); an isolated prefractionator
with any number of components admits exact positive finite-stage profiles with traces ≤ C η⁻¹ ρⁿ (GPT-6-Astra,
`total/attainment-astra/attainment.md`, checked by S); the closed Petlyuk network at fixed flows R = max_j H_j + δ
reaches purities → 1 geometrically in three independent solvers (N, P, F; N = 3 impurities down to 1e-47, also
q = 0, 0.5, 1.3; N = 4 indicative). **Open:** a proof of finite-stage attainment for the closed network, N ≥ 3.

## Result (version 1, independently verified)

In Underwood's infinite-stage model (TG's formulation), for networks of two-section columns with interval streams
(containing Agrawal's regular-column space, 8/152/6128 configurations for N = 3/4/5), without reverse flows: in the
strict model R ≥ max_j H_j for every N; in the permissive model, true for N ≤ 4, for extreme peaks, for merge-free
mixed networks and for the FTC itself, false for N = 5 (configuration 161, j = 3, R/V_FTC ≈ 0.99946,
interval-certified); the FTC attains max_j H_j (Theorem achieve). Version 2 shows that the permissive point is not
physical. Details: `paper.tex` Sections 9–16 (and the archived `paper-v1.pdf`), `verification.md`.

## Reproduction

From the project root (Git Bash), with `outils/wsl` running WSL Fedora with the scientific Python environment; paths
below are relative to `chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/`.

```
# version 3 (attainment)
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/limit/scripts && python check_limit_jacobian.py"   # Theorem L checks
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/design/scripts && python check_design.py"         # (H1)-(H8) on random data
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/computation/scripts && python paper_tables.py"   # finite-stage tables
# version 2
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/network/scripts && python unit_identities_exact.py && python section_identity_exact.py"   # exact, seconds
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/foundations/scripts && python unitnet.py"            # ~30 s, ALL PASS
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/sections/scripts && python check_stage_cut.py 25 2"   # minutes
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/simulation/scripts && python test_stage_bound.py"   # see header
outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/foundations/scripts && python petlyuk_attain.py fixed 0.02 10,20,40"   # fixed-flow attainment, ~2 min
# version 1
outils/wsl "bash chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/theory/scripts/run_theory_checks.sh"          # ~2 min
outils/wsl "bash chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/computation/scripts/run_certificates.sh"      # ~40 s
py -3 templates/build_tex.py problems/halvorsen-skogestad-petlyuk-minimum-vapor/paper.tex   # from chemistry-open-problems/
```

## Files

| Path | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript, version 2 (Sections 1–8 new; Sections 9–16 = version 1) |
| `paper-v1.tex`, `paper-v1.pdf` | Version 1 of the manuscript as independently verified (archived copy) |
| `references.md` | Sources and novelty checks (version 1 and version 2) |
| `attainment/` | Version 3: `theorem.md` (integrated proof of finite-stage attainment), `FINAL-REPORT.md`; `design/`, `limit/`, `analysis/`, `computation/` (agents D, L, A, C; briefs `BRIEF-*.md`) |
| `total/foundations/` | Class 𝒜 (`class.md`), integration (`total-result.md`), checks (`scripts/unitnet.py`, `scripts/petlyuk_attain.py`) |
| `total/sections/` | Agent S: `sections.md` (proofs, section and junction lemmas, binary attainment), scripts |
| `total/network/` | Agent N: `network.md` (unit graph, device decomposition, cut identity, ties, sign potential, permissive sign path), exact scripts |
| `total/simulation/` | Agent P: independent MESH solver, bound tests, attainment study (`runs/`) |
| `total/literature/` | Targeted novelty search of the stage-level argument (`novelty.md`) |
| `total/attainment-astra/` | GPT-6-Astra: isolated-prefractionator attainment theorem (`attainment.md`, `FINAL-REPORT.md`), checks |
| `total/verification-early/` | Early independent check of the lower bound |
| `theory/`, `computation/`, `general/` | Version 1 agents (T, C, G) |
| `verification.md`, `verification/` | Independent verification of version 1 |
| `certificates/sources/` | Third-party PDFs (not distributed) |
