# Early independent verification of the stage-level lower bound (wave 2)

Verifier: Claude Opus 5.5 (effort xhigh), fresh
context, no sub-agent. Date 2026-10-20. Status: **final for this early round** (first written as provisional at the
coordinator's request, then completed with all numerical runs). Attainment (the upper half) was not verified here.

Files read: PROTOCOL.md, this folder's BRIEF-verifier.md, README.md, `total/foundations/class.md` (v2),
`total/sections/sections.md` (v1), `total/network/network.md` (v1), HS3 text (`certificates/sources/HS2003_part3.txt`,
§1, §6), Halvorsen's thesis §6.4–6.5 (text file). Not opened: any PROGRESS.md, FINAL-REPORT.md, total-result.md, other
BRIEF files, paper-v2-* drafts.

## Verdicts

| Claim | Verdict |
|---|---|
| 1. Stage-level bound R ≥ Σ_P max(0, φ_P(θ_j)) at every steady state of every arrangement of 𝒜 | **CONFIRMED** |
| 2. Pure products ⇒ R ≥ max_j H_j; R ≥ H_j − c_j ε; Murphree stages; permissive point not a physical limit | **CONFIRMED** (minor editorial corrections) |
| 3. Scope = HS's "any adiabatic arrangement" (and more); conjecture ⇔ attainment | **CONFIRMED WITH CORRECTIONS** (wording and stated scope boundaries; no counterexample, no gap in the proof) |

No modelling loophole was found within HS's assumptions (constant α, ideal VLE, CMO, constant pressure, no internal
heat exchange, one feed, continuous steady-state operation, equilibrium-stage thermodynamics).

## 1. Re-derivation (done by hand, then checked symbolically: `scripts/vexact.py`, all PASS)

Notation as in class.md: φ_e(θ) = v_e − Σ_i α_i n_{e,i}/(α_i − θ), G_x(θ) = Σ_i α_i x_i/(α_i − θ), S(x) = Σ α_i x_i.

(a) **Node balance.** φ is linear in (v, n). For any unit, Σ_out φ − Σ_in φ = (Σ_out v − Σ_in v) − Σ_i c_i(Σ_out n_i − Σ_in n_i),
c_i = α_i/(α_i − θ). The CMO energy balance makes the first bracket Q_u + (1 − q) f_u, the component balances make the second
f_u z_i; hence Σ_out φ − Σ_in φ − f_u φ_F(θ) = Q_u with φ_F(θ) = (1 − q) − Σ_i α_i z_i/(α_i − θ). This holds for **every real q**
(the "vapor" (1 − q) f is the CMO enthalpy of the feed; subcooled q > 1 and superheated q < 0 included), for split feeds
(each portion has the same z, q, so each has zero flux at θ_j), for pre-flashed feeds (an (E) unit with Q = 0 or a duty
receiving the feed), and for total condensers/vaporizers (Q_u defined by the same energy balance; note that if a feed
portion enters a (C)/(R) unit, "Q_u = −(total vapor in)" must include the feed's (1 − q) f — editorial). θ_j ∈ (α_{j+1}, α_j)
exists for every q because Σ α_i z_i/(α_i − θ) runs from −∞ to +∞ on each interval. Checked symbolically (vexact (B)) and
at 50 digits for q ∈ {−1.5, −0.3, 0, 0.4, 1, 1.7} (vexact (C1)).

(b) **Sign uniformity.** Liquid of flow l: φ = −l G_x. Vapor in equilibrium: α_i/(α_i − θ) = 1 + θ/(α_i − θ) and Σ y = 1 give
φ = −θ v Σ_i y_i/(α_i − θ) = −(θ v/S(x)) G_x(θ). So φ_vap/v ÷ φ_liq/l = θ/S(x) > 0 (θ_j > α_N > 0): both outflow phases of an
(E) unit, split into any number of streams, have the sign of −G_{x_u}(θ). (M), (C), (R): one composition and one phase out,
so one sign. Degenerate cases: G_{x_u}(θ_j) = 0 gives φ = 0 for all outflows (put the unit in 𝒰₊); a phase with zero flow
or a component absent changes nothing; an (E) unit releasing only one phase is trivially uniform; θ_j never equals an
α_i. Checked symbolically for N = 2..5 (vexact (A)).

(c) **Cut.** Summing (a) over 𝒰₊ (units whose outflows all have φ(θ_j) ≥ 0): internal streams cancel; streams and products
leaving 𝒰₊ contribute φ ≥ 0; streams entering 𝒰₊ come from 𝒰₋ units and contribute −φ ≥ 0; feed portions contribute
−φ_F(θ_j) = 0. Products leaving 𝒰₋ have φ ≤ 0, so Σ_{P from 𝒰₊} φ_P = Σ_P max(0, φ_P). Hence
Σ_{𝒰₊} Q_u = Σ_P max(0, φ_P) + (nonnegative crossing terms) and R ≥ Σ_{𝒰₊} max(Q_u, 0) ≥ Σ_{𝒰₊} Q_u. The exact defect
identity R − Σ_P max(0, φ_P) = R(𝒰₋) + Q(𝒰₊) + Σ_cross |φ_e| (sections.md Thm 3(c), network.md Prop. 5.1) is correct.
The sign of every term was re-checked; no sign or convention error.

**Claim 2.** Pure liquid product p: φ = α_p P/(θ_j − α_p) (> 0 iff p > j); pure vapor product: θ_j P/(θ_j − α_p) ≥ the liquid
value when p > j. Component balances over the whole network give Σ_P max(0, φ_P) ≥ H_j. Finite purity: I re-derived
c_j = max(α_j/(α_j − θ_j), α_{j+1}/(θ_j − α_{j+1})) (both write-ups' proofs are correct; network.md's version via
φ_vap ≥ φ_liq is slightly tighter). Checked at 50 digits on 4000 random impure product sets (liquid and vapor products,
several streams per component, overall balances exact): min slack −1e-50 (rounding) (vexact (C2)); for the worst binary
placement the ratio (H_j − B_j)/(c_j ε) is 0.87 ≤ 1, so the constant is of the right order. Murphree trays with a common
efficiency E ∈ [0, 1]: the decomposition (E unit on the fraction E of the vapor + (M) bypass + (M) mixer) reproduces the
tray balances exactly (re-derived). Permissive counterexample of result 12: it has pure products and R < max_j H_j, so by
Corollary 5 it is not the limit (in R) of any family of 𝒜 steady states with purity → 1. Correct, and trivially so.

## 2. Modelling red team (main risk)

The proof uses exactly: CMO energy balance at every unit; every unit releases either streams of one composition and
one thermal state, or saturated phases in mutual equilibrium; one feed composition and thermal state; R counts every
positive duty; steady state. I examined the candidate loopholes:

| Object | Status |
|---|---|
| Partial condenser/reboiler with equilibrium outlet; side reboilers/condensers; heaters/coolers on stages | in 𝒜 ((E) with duty) |
| "Non-equilibrium" partial condensation: condensate removed as formed (dephlegmator, Rayleigh) | closure of 𝒜 (chains of (E) + (M) mixers) |
| Fog-type partial condensation / partial vaporization giving phases of the **same** composition | in 𝒜: (M) splitter + (C) or (R). **Correction**: sections.md §3 Remark (a) calls such a unit "not sign-uniform, not a CMO equilibrium-stage object"; it is decomposable into sign-uniform units and is covered. |
| Entrainment, weeping, vapor/liquid bypasses, streams between columns at any stages, phase maldistribution at tees, mixing of two-phase streams | in 𝒜 (extra arcs; per-phase (M) splitters) |
| Cross-flow trays (non-mixed liquid; Murphree E_MV > 1 possible) | in 𝒜: liquid cells in series each contacted by a parallel portion of the vapor. **Correction**: network.md Lemma 1.1 last row and sections.md Remark (b) list "Murphree efficiency > 1" as not covered; E_MV > 1 arising from cross-flow is covered; only a well-mixed tray with E > 1 (unphysical) is not. |
| Rate-based contacting (packed sections, two-film model) with **component-independent** film coefficients | closure of 𝒜: in each height element equilibrate a fraction ε of the vapor with a fraction ε' of the liquid and remix; this gives exactly N = k_y(y − y_I) = k_x(x_I − x), y_I = K(x_I) x_I, k_y = εV/Δh, k_x = ε'L/Δh (equimolar, CMO-consistent). **Correction**: sections.md Remark (b) "rate-based models are not covered" is too pessimistic; only component-specific coefficients (Maxwell–Stefan with unequal diffusivities, diffusion distillation) are not. These are outside HS's equilibrium-stage model; whether the bound survives them is open (not a gap in HS's sense). |
| Subcooled/superheated internal streams (pumparound coolers, superheaters, subcooled reflux) | covered: either move the duty to the receiving unit (class.md) or, cleaner, use the enthalpy flux φ_e = h_e − Σ c_i n_{e,i}: a heater/cooler on a stream is a single-outflow unit, hence sign-uniform; R only grows. Suggest stating this generalization explicitly. |
| Feed subcooled/superheated, split, pre-flashed, partially heated portions | in 𝒜 (Section 1(a)) |
| Several feeds of different composition / thermal state | outside HS; bound becomes R ≥ Σ_P max(0, φ_P(θ)) − Σ_F max(0, φ_F(θ)) (re-derived; tested T4) |
| Internal heat exchange / heat pumps / several pressures / non-ideal VLE / decanters / membranes | outside HS (explicitly excluded). With internal exchange counted in R the bound holds. |
| **Non-steady (batch, cyclic "controlled-cycling") operation** | **not covered and not mentioned** in class.md/network.md: with holdups, time-averaged outflows of a tray need not be sign-uniform (vapor and liquid periods at different compositions). Outside HS's continuous steady-state setting; should be listed in class.md §7 as an exclusion. |

Consistency with HS's own numbers (Th §6.4–6.5, α = 4, 2, 1, z = 1/3, q = 1, V_Petlyuk = 1.366): the reversible
(non-adiabatic, in the closure of 𝒜) arrangement needs total heat 15/9 = 1.667 ≥ 1.366; Example 3 (feed preheated by the
sidestream to q = 0.762) has reboiler 1.181 plus 0.238 internally delivered = 1.419 ≥ 1.366; the external-heat values
1.222, 1.181, 1.000 and the second-law limit 0.7925 are below 1.366 only because internal heat is not counted, as class.md
says. No conflict.

## 3. Claim 3 (scope and "conjecture ⇔ attainment")

- 𝒜 is a faithful (and larger: non-adiabatic sections allowed) formalization of HS's arrangements, given the limit
  definition of "infinitely many stages and pure products" (class.md §4). The limit definition deliberately excludes
  Underwood-level objects that are not limits of finite-stage steady states (such as the permissive point); I agree with
  this reading (HS3 and the 2002 paper say infinite stages are used as the limit/lower bound of real columns).
- HS's conjecture has two parts: (A) the generalized adiabatic Petlyuk (HS3 Fig. 1: one reboiler, one condenser, all
  couplings direct, N pure products) needs no more vapor than any other arrangement; (B) its minimum is max_j H_j (HS3 §5,
  proved there only for the extreme peaks). Given the lower bound for all of 𝒜: attainment (R*(𝒮_Petlyuk) ≤ max_j H_j) ⇒ (A)
  and (B). Conversely, if attainment failed, (B) would be false while (A) could still hold. So "conjecture ⇔ attainment" is
  correct for the conjecture **including HS's formula**; for (A) alone, attainment is sufficient, not necessary.
  **Correction**: state this explicitly. Also "lowest" must be read as "lower or equal" (many arrangements tie: HS3 §7.1).
- What attainment must show: a sequence of steady states of the HS Fig. 1 topology (any stage numbers, any internal
  splits; HS allow some heat at product outlets — still in 𝒜) with product purities → 1 and R → max_j H_j; equivalently
  lim_{ε→0} inf_{stages, operation} R ≤ max_j H_j. No order-of-limits issue for the lower bound (it is uniform in the stage
  numbers); for attainment the natural route is V = max_j H_j + δ, stages → ∞ (purity → 1), then δ → 0. Result 12's
  Underwood-model attainment does not by itself give this (finite-stage convergence for N ≥ 3 is the missing step
  identified in sections.md §8.2; N = 2 is proved there, Prop. 10, which I did not re-check line by line).

## 4. Editorial corrections (team files not edited)

1. network.md Prop. 6.3 uses σ = + for **light-ish** units while network.md Lemma 3.2 / sections.md use 𝒰₊ / σ_u = +1 for
   **heavy-ish** units: rename to avoid sign confusion.
2. class.md §3 "R also equals the total cooling plus that constant": the "constant" (vapor in products) − (1 − q) is
   constant only for fixed product phases (e.g., all liquid products).
3. sections.md §3 Remark (a) and Remark (b), network.md Lemma 1.1 last row: see the table in Section 2 (same-composition
   two-phase units, cross-flow E_MV > 1, scalar-coefficient rate models are all covered).
4. class.md §7: add non-steady (batch/cyclic) operation and component-specific mass-transfer rates as exclusions.
5. Literature: some reviews state that HS "proved" FTC optimality (e.g., a 2013 Caballero–Grossmann review, found by web
   search); HS3 itself calls it a conjecture — worth a sentence in the paper.

## 5. Numerical tests (own solver, `scripts/vnet.py`, `scripts/vtests.py`)

Own CMO/constant-α network solver written from scratch (total flows from a linear system given split fractions and
duties; component flows from one linear system per component given one scalar S_u per two-phase (E) unit; Newton on
Σ_i y_{u,i} = 1 with analytic Jacobian, checked against finite differences: error 8e-8 on a 91-unknown column).
At every converged steady state and every θ: node-balance residual, sign uniformity of every unit, exact defect identity,
bound.

Logs: `scripts/vtests.log` (T1, T2, T6; seed 1), `scripts/vtests_random.log` (T3–T5; seed 3). Runtime ≈ 10 min + 1 min.

- **T1** binary columns (α, z, q) ∈ {(2, .5, 1), (1.5, .3, 0), (3, .7, .5), (1.4, .5, 1.4 subcooled), (2.5, .4, −0.3
  superheated)}, V_B = (1 + δ) H_1, δ ∈ {0.05, 0.01}, 10–120 stages per section (continuation in the stage number):
  65 converged steady states (5 long-column runs did not converge: solver issue, not a violation), **0 violations**,
  residuals ≤ 8e-13; R − (H_1 − c_1 ε) ≥ 6.1e-3 everywhere. At n = 120 the impurity is 1e-5 … 1e-15 and
  R − B_1 → δ H_1 (e.g. α = 2.5, q = −0.3, δ = 0.01: impurity 9e-13, R − B_1 = 6.1e-3 = δ H_1): the bound is tight up
  to the excess vapor, as Corollary 5.4 of network.md predicts.
- **T2** ternary Petlyuk (prefractionator at HS's preferred split, main column with liquid side draw), three feeds
  (α = 4,2,1 equimolar q = 1; α = 6,2.5,1, z = .2,.5,.3, q = 1; α = 4,2,1, z = .3,.3,.4, q = 0.7), V_B = (1 + δ) max_j H_j,
  8–48 stages per section: all steady states satisfy the bound for both j, residuals ≤ 1e-13. With δ = 0.01 and 48 stages
  per section, impurity ≈ 1–3e-4 and the slack at the highest peak is 0.011–0.014 ≈ δ V_Petlyuk: numerical evidence that
  the finite-stage Petlyuk approaches max_j H_j (attainment), not a proof.
- **T6** binary column with Murphree trays (E = 1, 0.7, 0.4; (E) unit + vapor bypass + vapor mixer), V_B = 1.02 H_1, up
  to 120 trays: bound holds, R − B_1 → 0.03 = 0.02 H_1 as purity → 1.
- **T3** 3000 random networks (N = 2..5, 2–12 units of types E/ML/MV/C/R, 1–3 destinations per outflow phase including
  products of both phases, 1–3 feed portions, random duties on E units, q ∈ [−0.6, 1.6]); 2354 contain recycles. All
  solved; sign uniformity never violated; node-balance and defect-identity residuals ≤ 1.1e-9; min relative slack
  −4e-13 (rounding; trivial ties such as all-light-ish networks).
- **T4** 1000 random networks with several feeds of different compositions and thermal states, 6 random θ each: the
  bound with the feed term, R ≥ Σ_P max(0, φ_P(θ)) − Σ_F max(0, φ_F(θ)), always holds (min relative slack −2e-13).
- **T5** adversarial Nelder–Mead minimisation of the slack: only a few topologies had a feasible random start (weak
  test); no negative slack found.
- Symbolic/high-precision checks: `scripts/vexact.py` ALL PASS (log `scripts/vexact.log`).
- Team scripts re-run: `network/scripts/section_identity_exact.py` PASS, `unit_identities_exact.py` PASS,
  `sections/scripts/check_stage_cut.py 25 2` PASS (columns, Petlyuk, random networks).

Reproduction (from the project root, Git Bash): `outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/verification-early/scripts && python vexact.py && python vtests.py 1 0.5"`.

## 6. Novelty (quick)

OpenAlex: 119 works cite HS3 (doi:10.1021/ie0108651); titles 2003–2026 scanned: none proves the conjecture or a
network-level stage bound (closest: Agrawal–Tawarmalani group MINLP/global-optimization papers, which compute rather than
prove; "Liquid-only transfer ... analytical minimum vapor expressions" 2025). Web searches found no prior proof of a
stage-level Underwood bound for arbitrary networks. The single-column statement (Underwood's V_min is a lower bound at
any number of stages) is classical in spirit; the network cut argument appears new as far as this quick search goes
(the dedicated literature agent should confirm).
