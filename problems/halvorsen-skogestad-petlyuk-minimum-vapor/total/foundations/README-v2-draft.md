# Halvorsen–Skogestad conjecture: the generalized Petlyuk arrangement minimizes total vaporization

**Attribution.** *Version 2 (2026-10-20, second wave).* Foundations, class 𝒜, red team and integration: Claude Opus
5.5 (effort xhigh) (Agent F). Sections and stage-level lemmas, binary attainment: Claude Opus 5.5 (effort xhigh)
(Agent S). Network argument, device decomposition, tie structure: Claude Opus 5.5 (effort xhigh) (Agent N). Physical
simulation (independent MESH solver, tests, attainment study): Claude Opus 5.5 (effort xhigh) (Agent P). Literature:
Claude Opus 5.5 (effort xhigh). Attainment (prefractionator theorem): GPT-6-Astra via Codex CLI (effort xhigh). Independent verification of version 2: **pending**.
*Version 1 (independently verified).* Theory: Claude Opus 5.5 (effort xhigh) (Agent T). Computation: Claude Opus 5.5
(effort xhigh) (Agent C). General N: Claude Opus 5.5 (effort xhigh) (Agent G). Independent verification: Claude Opus
5.5 (effort xhigh).

**Status (version 2; independent verification pending).** 
- **Lower bound: solved (preliminary) in the broadest class, for every N.** For every arrangement of the class 𝒜
  (any finite network of constant-molar-overflow equilibrium stages, single-phase splitters and mixers, total
  condensers and vaporizers; one feed split and conditioned at will; heat duties anywhere; recycles; products of
  either phase), at every steady state with any numbers of stages and any purities, R ≥ Σ_P max(0, φ_P(θ_j)) for
  every j; with purity ≥ 1 − ε, R ≥ H_j − c_j ε; hence, with infinitely many stages and pure products,
  R ≥ max_j H_j, the Petlyuk value. This closes every gap left by version 1 (degenerate points, reverse flows,
  several feed stages, splitters, recycles) and proves that version 1's permissive counterexample is not physical.
- **Attainment.** The generalized Petlyuk arrangement attains max_j H_j in Underwood's infinite-stage model
  (version 1, verified). In the finite-stage model: **proved for N = 2**; for N ≥ 3, steady states exist for any
  prescribed flows, the isolated prefractionator is proved to admit exact finite profiles with exponentially small
  traces (any number of components), and the full Petlyuk network reaches purities → 1 at R = max_j H_j + δ in
  simulations (impurities down to 1e-47 for N = 3); @@ATTAIN_N3_STATUS@@
- **Conjecture.** Solved (preliminary) for N = 2. For N ≥ 3: every arrangement of 𝒜 needs at least the Petlyuk
  minimum max_j H_j computed by Halvorsen and Skogestad (proved); that the Petlyuk arrangement realizes this value as
  a finite-stage limit is proved in Underwood's model and supported numerically, @@N3_PHRASE@@
- Version 1's results in the Underwood model are unchanged (independently verified, with corrections).

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
minimum, in the finite-stage-limit sense). @@ATTAINMENT@@

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
# version 2
outils/wsl "cd <problem>/total/network/scripts && python unit_identities_exact.py && python section_identity_exact.py"   # exact, seconds
outils/wsl "cd <problem>/total/foundations/scripts && python unitnet.py"            # ~30 s, ALL PASS
outils/wsl "cd <problem>/total/sections/scripts && python check_stage_cut.py 25 2"   # minutes
outils/wsl "cd <problem>/total/simulation/scripts && python test_stage_bound.py ..." # seed and trial count as in its header
# version 1
outils/wsl "bash <problem>/theory/scripts/run_theory_checks.sh"          # ~2 min
outils/wsl "bash <problem>/computation/scripts/run_certificates.sh"      # ~40 s
py -3 templates/build_tex.py problems/halvorsen-skogestad-petlyuk-minimum-vapor/paper.tex   # from chemistry-open-problems/
```

## Files

| Path | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript, version 2 (Sections 1–8 new; Sections 9–16 = version 1) |
| `paper-v1.tex`, `paper-v1.pdf` | Version 1 of the manuscript as independently verified (archived copy) |
| `references.md` | Sources and novelty checks (version 1 and version 2) |
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
