# Total result — the Halvorsen–Skogestad conjecture in the class 𝒜 (Agent F, integration)

**Final version, 2026-10-20 15:10** (integrated into the top-level README.md and paper.tex v2). Integrates Agents S (`../sections/sections.md`), N (`../network/network.md`),
P (`../simulation/`), F (`class.md`, `scripts/unitnet.py`), and, where stated, the literature agent
(`../literature/`) and GPT-6-Astra (`../attainment-astra/`). Independent verification of the new claims:
**pending**. Paths are relative to `total/foundations/` unless stated.

## 1. Statement of the problem and of the class

HS3 (Ind. Eng. Chem. Res. 42 (2003) 616–629), abstract and §6.3: the generalized Petlyuk arrangement has the lowest
need for vaporization of any adiabatic distillation arrangement without internal heat exchange, for constant relative
volatilities, constant molar flows, constant pressure, an arbitrary feed and pure products, with infinitely many
stages. HS gave no definition of "any arrangement" and no complete proof; they also proved the Petlyuk formula
max_j H_j only for the extreme peaks (§5.3).

Formalization (`class.md`, self-contained): 𝒜 = all finite networks of CMO equilibrium stages (E), single-phase
mixers/splitters (M), total condensers (C) and total vaporizers (R), one process feed (composition z, quality q) split
and conditioned at will, heat duties at any unit, recycles, products of either phase; objective R = Σ_u max(Q_u, 0)
(total heat supplied = total vaporization); infinitely many stages and pure products as the limit
R*(𝒮) = lim_{ε→0} inf{R : steady state of a member of 𝒮 with purity ≥ 1 − ε}. Every device of HS-type arrangements
decomposes into units of 𝒜 (network.md §1.2), including Murphree trays with a common efficiency; 𝒜 also contains
non-adiabatic sections with external utilities. Excluded (and excluded by HS): uncounted heat exchange between process
streams, several pressures, non-ideal VLE, decanters/membranes/component-dependent efficiencies, heat pumps, several
feeds.

## 2. Theorem 1 (lower bound; complete proof)

For θ > 0, θ ∉ {α_i}, the Underwood flux of a stream is φ_e(θ) = v_e − Σ_i α_i n_{e,i}/(α_i − θ). With θ_j the feed
roots, H_j = Σ_{p>j} α_p z_p/(θ_j − α_p), 𝒰₊ the units whose outflows have φ(θ_j) ≥ 0 and 𝒰₋ the others:

**Theorem 1.** At every steady state of every arrangement of 𝒜, for every j,

  R − Σ_P max(0, φ_P(θ_j)) = R(𝒰₋) + Q(𝒰₊) + Σ_{e: 𝒰₊→𝒰₋} φ_e(θ_j) − Σ_{e: 𝒰₋→𝒰₊} φ_e(θ_j) ≥ 0,

and, if every product stream has purity ≥ 1 − ε, R ≥ H_j − c_j ε with c_j = max(α_j/(α_j − θ_j),
α_{j+1}/(θ_j − α_{j+1})). Hence R*(𝒮) ≥ max_j H_j for every 𝒮 ⊆ 𝒜.

Proof: `class.md` §5 (Lemma 5.1 node balance from CMO and linearity; Lemma 5.2 sign uniformity: liquid
φ = −l G_x(θ), vapor in equilibrium φ = −(θ v/S(x)) G_x(θ), G_x(θ) = Σ α_i x_i/(α_i − θ), so all outflows of a unit
share a sign; Theorem 5.3 summation over 𝒰₊, the feed having zero flux at θ_j; Corollaries 5.4–5.5). Same proof, with
section and junction lemmas: sections.md §2–4; with the device decomposition, the cut identity for arbitrary unit
sets and the tie structure: network.md §1–5. Found independently by N (13:38), F (13:44–13:47) and S (13:45); algebra
re-derived by P.

Checks (evidence only): exact rationals and sympy (N: `unit_identities_exact.py`, `section_identity_exact.py`); four
independent MESH solvers (F `scripts/unitnet.py`: ALL PASS, 396 random networks + binary columns, bound tight to
1e-11 at 80+80 stages below V_min; S `check_stage_cut.py`; N `cut_bound_check.py`; P `test_stage_bound.py`).

## 3. Consequences

- **Excess = three pinch losses** (Theorem 1's identity): heat supplied on the light side, heat removed on the heavy
  side, flux across the hyperplane {x : G_x(θ_j) = 0} (which passes through the feed-flash liquid for 0 ≤ q ≤ 1).
- **Ties / optimal structure** (network.md §5): approaching max_j H_j requires vanishing losses: no reboiler/heater
  delivering a light product p ≤ j*, no condenser/cooler delivering a heavy product p > j*, θ_{j*} asymptotically a
  root of every section whose liquid changes side (Underwood's active root and pinch, derived for every
  arrangement); finite stages never tie.
- **Sections** (sections.md Lemma 6; network.md Cor. 6.2): the side changes at most once along a section.
- **Gaps of version 1 closed for the lower bound**: degenerate points, reverse flows, several feed stages,
  intermediate sections, splitters, remixing, recycles, bypasses, non-interval streams, products of either phase.
- **The permissive counterexample of version 1 is unphysical**: it has the light→heavy sign path
  C → CDE → BCDE → BCD → D at θ_3 (network.md §6.4), impossible at a steady state; by Theorem 1 no family of steady
  states with purity → 1 approaches it.
- **Internal heat exchange**: the bound is on the heat delivered to heavy-side units, whatever its origin; external
  heat can fall below max_j H_j only by moving heat from light-side to heavy-side units (Halvorsen's examples).

## 4. Attainment (the remaining step)

By Theorem 1, attainment R*(𝒮_Petlyuk) = max_j H_j (HS's Petlyuk formula in the finite-stage-limit sense) implies
the conjecture R*(𝒮_Petlyuk) ≤ R*(𝒮) ∀𝒮 (sufficient, not necessary).
- Infinite-stage Underwood model: proved (version 1, Theorem achieve, independently verified; HS's preferred-split
  construction).
- **N = 2, finite stages: proved** (S, sections.md Prop. 10: shooting on x_B, monotone pinch maps, IVT; F re-checked
  the quadratic f(s) with roots H_1 and s_−). So for N = 2 the conjecture holds in 𝒜 completely.
- **Existence for prescribed flows** (N, network.md Theorem 10.1: M-matrix structure per component, automatic
  normalization, a priori bounds, degree along α^t): R is fixed by the flows, so attainment = purity of some steady
  state at flows F_δ with R = max_j H_j + δ.
- Isolated prefractionator with any number of components (Astra, attainment.md Theorem 1): exact nonnegative finite
  profiles with exponentially small traces. Network coupling: open (zero modes of the limit system at excess-carrying,
  doubly pinched sections: N §10.4, S §8.2).
- Numerics: fixed flows F_δ (N §10.3): N = 3 impurities 1.4e-2 → 1.3e-14 (n = 20 → 160, δ = 0.002); N = 4 → 2e-6
  (n = 120, indicative); N = 5 decreasing to n = 40. Optimization at fixed purity (P): R − bound(ε) → 0 as n grows
  (3e-8 at n = 80, ε = 0.01, N = 3); F's separate solver/optimizer reproduces P's optimum at n = 20 (R = 1.372571)
  and, at fixed flows R = max H + 0.04, gives impurities 6.7e-5 at n = 40 for q = 1 and ≤ 1.4e-4 for q = 0, 0.5, 1.3.
- N = 3 closure (S, N, not completed): sharp limit at fixed flows fixes all net flows, Astra's feed-stage vector and
  pure B at the side draw; the two junction ratios are neutral modes; with the prefractionator's two trace net flows
  they form four unknowns fixed by exponentially small balances (network.md §10–11, sections.md §8.5–8.6).

## 5. Final statement

**Theorem (HS conjecture in 𝒜, lower bound; complete proof).** For every arrangement of 𝒜, every steady state with
purity ≥ 1 − ε satisfies R ≥ max_j (H_j − c_j ε); hence R*(𝒮) ≥ max_j H_j for every family 𝒮 ⊆ 𝒜. The generalized
Petlyuk arrangement attains max_j H_j in Underwood's infinite-stage model (every N) and in the finite-stage model for
N = 2. So no arrangement of 𝒜 needs less vaporization than the Petlyuk minimum computed by HS; the conjecture is
completely proved (finite-stage sense) for N = 2; for N ≥ 3 the finite-stage attainment by the Petlyuk is supported
by three independent simulations but not proved.

Status in PROTOCOL vocabulary: lower bound **solved (preliminary)**; conjecture **solved (preliminary)** for N = 2 and
for N ≥ 3 relative to the Petlyuk value as computed by HS; finite-stage attainment for N ≥ 3 **open**. Independent
verification: pending (early check of the lower bound: CONFIRMED, provisional, `../verification-early/`).

## 6. What remains open

- Finite-stage attainment of max_j H_j by the closed Petlyuk network for N ≥ 3 (reduction for N = 3 described above).
- Minimum *external* heat with internal heat exchange (HS's "unknown"); non-ideal mixtures (sign uniformity uses
  constant α); several feeds and non-sharp products (bounds hold, attainment not studied).

## 7. Suggested row for the root README (orchestrator)

| 12 | ... | Chemical engineering (distillation) | Version 2: every equilibrium-stage arrangement (class 𝒜: any network of CMO stages, splitters, mixers, condensers, vaporizers, recycles, duties anywhere) needs at least the Petlyuk value max_j H_j (elementary stage-level proof, any N); conjecture fully proved for N = 2; for N ≥ 3 finite-stage attainment by the Petlyuk supported numerically, not proved. Version 1 (Underwood model) verified | Claude Opus 5.5 (xhigh): foundations, sections, network, simulation, literature (wave 2) and theory, computation, general N (wave 1); GPT-6-Astra (xhigh): prefractionator attainment theorem; verified (v1) by Claude Opus 5.5 (xhigh) | **lower bound solved (preliminary)**, v2 verification pending |
