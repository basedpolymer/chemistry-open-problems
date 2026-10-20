# The global argument on arbitrary networks (Agent N)

Agent N (network), Claude Opus 5.5 (effort xhigh). Version 1, 2026-10-20. Paths are relative to the problem folder
`problems/halvorsen-skogestad-petlyuk-minimum-vapor/`. Reference definitions: `total/foundations/class.md` (class 𝒜,
finite-stage model, limit definition of "infinitely many stages and pure products"). The stage-level bound of
Section 4 was found independently, within a few minutes of each other, by Agents N (13:38), S (13:45) and F (13:47);
this file gives the network-side write-up: the unit graph and the decomposition of physical junctions (Section 1),
every node balance and the cut identity for arbitrary partitions (Section 2), the sign lemma (Section 3), the bound
(Section 4), the exact slack and the tie structure (Section 5), sections, the sign potential and the relation with
result 12, including what the permissive counterexample violates (Section 6), extensions (Section 7), checks
(Section 8), what is not covered (Section 9), and, for the attainment work, existence of steady states for
prescribed flows (Section 10).

**Status of the claims.** Sections 2–7 and Theorem 10.1 are complete proofs (hand proofs; the scripts of Section 8
only check them). They establish, for every arrangement of 𝒜, the lower bound R*(𝒮) ≥ max_j H_j. Together with the
attainment R*(𝒮_Petlyuk) ≤ max_j H_j, which is **not** proved in this file (proved in the infinite-stage Underwood
model in result 12, Theorem achieve; in the finite-stage model for N = 2 by Agent S; general N: Agent Astra), this is
the Halvorsen–Skogestad conjecture in class 𝒜. Sections 6.5, 10.3 and 10.4 are remarks or numerical observations,
marked as such. Section 11 is a SKETCHED closure plan for N = 3 (forced impurities pin the neutral modes), with
numerical support; it is not a proof.

## 0. Statements

Data: α_1 > … > α_N > 0, feed F = 1 of composition z (z_p > 0) and liquid fraction q ∈ ℝ, θ_j ∈ (α_{j+1}, α_j) the
roots of Σ_p α_p z_p/(α_p − θ) = 1 − q, H_j = Σ_{p>j} α_p z_p/(θ_j − α_p), V_FTC = max_j H_j.

For a stream e with component flows n_e (≥ 0) and vapor flow v_e (= Σ_i n_{e,i} for a vapor, 0 for a liquid), and
θ > 0, θ ∉ {α_i}, the **stream flux** is φ_e(θ) = v_e − Σ_i α_i n_{e,i}/(α_i − θ). For a liquid of composition x,
G_x(θ) := Σ_i α_i x_i/(α_i − θ).

- **Theorem 4.1 (stage-level cut bound).** At every steady state of every arrangement of 𝒜 and for every θ > 0,
  θ ∉ {α_i}: let 𝒰₊ be the set of units all of whose outflows have φ ≥ 0 at θ (every unit is of one of the two
  kinds, Lemma 3.2). Then
  R ≥ R(𝒰₊) ≥ Σ_{u∈𝒰₊} Q_u = Σ_P max(0, φ_P(θ)) − φ_F(θ)·f(𝒰₊) + Σ_{e crossing} |φ_e(θ)|,
  where P runs over the product streams, f(𝒰₊) ∈ [0, 1] is the part of the feed entering 𝒰₊, and the last sum runs
  over the internal streams between 𝒰₊ and its complement (each such term is nonnegative with the sign shown).
- **Corollary 4.2.** At θ = θ_j the feed term vanishes: R ≥ Σ_P max(0, φ_P(θ_j)) ≥ H_j − c_j ε when every product
  stream has purity ≥ 1 − ε, with c_j = max(α_j/(α_j − θ_j), α_{j+1}/(θ_j − α_{j+1})). For pure products R ≥ H_j.
- **Corollary 4.3.** R*(𝒮) ≥ max_j H_j = V_FTC for every family 𝒮 ⊆ 𝒜 (class.md §4).
- **Proposition 5.1 (exact slack).** R − Σ_P max(0, φ_P(θ_j)) = R(𝒰₋) + Q(𝒰₊) + Σ_{e crossing} |φ_e(θ_j)|.
  Corollary 5.2–5.4: necessary conditions for an arrangement to approach V_FTC (no heating of light-ish units, no
  cooling of heavy-ish units, θ_j asymptotically an Underwood root of every section whose liquid changes side).
- **Proposition 6.3 (sign potential).** Give each unit the sign σ = + if G_x(θ_j) > 0 for its outflows ("light-ish",
  φ < 0), σ = − otherwise. Along every arc of result 12's sign graph (flux flowing from one node to another, through
  a stream or a section) σ never goes from + to −. Hence no light product is joined to a heavy product by a sign
  path: result 12's reachability argument holds at every physical steady state, with no column lemma, and the
  permissive counterexample, which has the sign path C → CDE → BCDE → BCD → D at θ_3, is not a limit of physical
  steady states (Section 6.4).
- **Theorem 10.1 (existence, for attainment).** For every prescribed set of flows of any arrangement satisfying a
  reachability condition (G), a steady state exists (Brouwer degree, homotopy in α); hence R is fixed by the flows and
  attainment reduces to the purity of steady states at fixed flows (Corollary 10.2). Numerical evidence and a caution
  (a neutral "front" mode of the infinite-stage limit system) in 10.3–10.4.

## 1. The unit graph

### 1.1 Units and streams (class.md §3)

An arrangement is a finite directed multigraph. Vertices: units of four types — (E) equilibrium stage (any inflows,
duty Q_u ∈ ℝ, outflows = its saturated liquid of composition x_u and its saturated vapor y_u = K(x_u) x_u, each
divided into any number of streams; K_i(x) = α_i/S(x), S(x) = Σ_k α_k x_k), (M) single-phase mixer/splitter (Q = 0),
(C) total condenser, (R) total vaporizer — plus the feed source (the process feed divided into portions of its
composition and thermal state) and the product outlets. Arcs: streams, each with a phase. A steady state satisfies at
every unit the component balances, the CMO energy balance Σ_out v − Σ_in v = Q_u (a feed portion of size f brings
vapor (1 − q) f) and the composition rule of its type. R := Σ_u max(Q_u, 0). Cycles (recycles) are allowed; nothing is
assumed about the signs of net flows inside sections, about stability or uniqueness of the steady state, or about
the topology.

### 1.2 Decomposition of physical devices into units (Lemma 1.1)

**Lemma 1.1.** Every device in the list below is, exactly (same balances, same compositions of every stream leaving
it), a sub-network of units of types (E), (M), (C), (R), or a limit of such sub-networks; hence every arrangement built
from them belongs to 𝒜 (or is a limit of members of 𝒜, to which the bounds pass because they are continuous in the
stream data).

| Device | Units |
|---|---|
| tray, packing element (equilibrium stage) | (E), Q = 0 |
| feed stage, side-draw stage, stage receiving several streams of any origin (mixing of different streams, thermal couplings, merging of links) | one (E) receiving all the streams |
| liquid or vapor withdrawn to several destinations (side draws, stream splitters dividing a section product among several columns, dividing-wall liquid/vapor splits, bypasses) | the (E) outflow divided into several streams (no extra unit), or an (M) splitter on a transferred stream |
| liquid or vapor streams of one phase mixed before entering a unit | (M) mixer |
| two-phase transfer stream whose phases are in equilibrium (flash, partial vaporizer, partial condenser, kettle reboiler, thermosiphon) | (E) with duty; its two outflows travel together or separately |
| two-phase stream whose phases are not in equilibrium (e.g. a liquid and a vapor from different units carried in one pipe) | two single-phase streams; if they are brought into contact and separated, an (E) |
| phase maldistribution at a tee (the two phases of a stream divided in different proportions) | one (M) splitter per phase |
| total condenser with reflux drum; partial condenser (TG's saturated-vapor transfer) | (C) followed by an (M) splitter; partial condenser = (E) with Q < 0 |
| total reboiler / total vaporizer of a liquid; product withdrawn before it | (M) splitter then (R) |
| intermediate reboiler or condenser on a stage, side heater, pumparound cooler returning to the column (in CMO a subcooled return is a condensing duty at the receiving stage) | duty on an (E) unit |
| feed: split among several stages, preheated, precooled, partially vaporized, pre-flashed and its phases sent to different places | feed portions (external inflows of flux φ_F per unit), (E) units with duties for conditioning, an (E) with Q = 0 for a flash |
| products of either phase, several streams for one product, products at any purity | outflows of any unit sent to outlets |
| Murphree stage with vapor efficiency E ∈ [0, 1] (Agents S and F) | (E) fed with the fraction E of the incoming vapor + (M) vapor bypass of the fraction 1 − E + (M) mixer: L x_in + E V y_in = L x + E V K(x) x |
| differential (Rayleigh) vaporization or condensation, dephlegmator with heat removal along it | limit of chains of (E) units with duties whose vapor (liquid) outflows are collected by (M) mixers |
| several pressures, non-ideal VLE, liquid–liquid split (decanter), membranes, reactions, heat exchange between process streams not counted in R, Murphree efficiency > 1 | **not** decomposable; outside HS's assumptions (class.md §7) |

*Proof.* Each row is checked by writing the component balances and the CMO energy balance of the device and of the
sub-network: they coincide, and every stream leaving the device has the same phase, flow and composition in both. For
the Murphree row: a Murphree stage is defined by y_out = y_in + E (K(x)x − y_in) and the stage balances
L_in x_in + V y_in = L x + V y_out (stage index omitted); the sub-network has the (E) unit balance
L_in x_in + E V y_in = L x + E V K(x) x and the mixer output (1 − E) V y_in + E V K(x) x = V y_out; adding gives the
same balances, and the liquid leaving is x in both. For the differential devices: the Rayleigh integral is the limit of
a chain of m equilibrium flashes, each removing a fraction 1/m of the vapor, as m → ∞; the chain is a member of 𝒜 for
every m. ∎

The decisive property used below (Lemma 3.2) is that every unit releases streams of a **single composition per
phase, the two phases being in equilibrium if both occur**. A device that releases two streams of different
compositions not in phase equilibrium is the only kind of building block that could escape the bound; the last row
lists the ones known, all excluded by HS.

### 1.3 Junctions and sections

A **section** is a chain of (E) units with Q = 0, stages 1..n from bottom to top, with constant liquid L and vapor V
(CMO), between a bottom-end unit X and a top-end unit Y: X sends vapor V to stage 1, stage n sends vapor V to Y, Y
sends liquid L to stage n, stage 1 sends liquid L to X, and the stages exchange no other stream. Its **net upward
flow** is w = V y_k − L x_{k+1} (constant along the section by the stage balances, also at the X/1 and n/Y
interfaces, with y_0 = composition of the vapor leaving X and x_{n+1} = composition of the liquid leaving Y). A
**junction** is any unit (or set of units) that is not inside a section. Result 12's networks of two-section columns,
TG's and Agrawal's regular-column configurations, Fidkowski's (2006) configurations with any of his ten types of
interconnecting streams, Christiansen–Skogestad–Lien's extended Petlyuk arrangements, multi-feed and side-stream
columns (Nikolaides–Malone 1987), side strippers/enrichers (Carlberg–Westerberg 1989) and the Kaibel column are
arrangements of this form.

## 2. Fluxes, node balances and the cut identity

**Lemma 2.0 (phase relation).** For the same component flows n, φ_vapor(n) − φ_liquid(n) = Σ_i n_i ≥ 0.

**Lemma 2.1 (node balance).** At a steady state, for every unit u and every θ > 0, θ ∉ {α_i},
  Σ_{e out of u} φ_e(θ) − Σ_{e into u} φ_e(θ) = Q_u + f_u φ_F(θ),
where f_u is the size of the feed portion entering u (0 if none), φ_F(θ) := (1 − q) − Σ_i α_i z_i/(α_i − θ) is the
flux of a unit of feed, and product streams count as outflows. Written out per type:
- (E): Σ_out v − Σ_in v − (1 − q) f_u = Q_u (energy) and Σ_out n_i − Σ_in n_i = f_u z_i (components); multiply the
  i-th component balance by α_i/(α_i − θ), sum over i and subtract from the energy balance.
- (M): Q_u = 0, no feed, all streams of one phase: Σ_out φ − Σ_in φ = (Σ_out v − Σ_in v) − Σ_i α_i(Σ_out n_i −
  Σ_in n_i)/(α_i − θ) = 0.
- (C): vapor in, liquid out of the same component flows: Σ_out φ − Σ_in φ = −Σ_in v = Q_u (Lemma 2.0).
- (R): liquid in, vapor out: Σ_out φ − Σ_in φ = +Σ_in Σ_i n_i = Q_u.
- Feed source: the external inflow of flux φ_F(θ) per unit of feed (each portion has composition z and thermal state
  q, so the same flux per unit).
- Product outlet: an external outflow (its flux is that of the stream).
*Proof.* Linearity of φ in (v, n); the energy balance is the CMO balance of class.md §3. ∎

**Proposition 2.2 (cut identity for an arbitrary set of units).** For every set W of units,
  Σ_{u∈W} Q_u = Σ_{e: W→W^c} φ_e − Σ_{e: W^c→W} φ_e + Σ_{P from W} φ_P − φ_F · f(W),
with f(W) = Σ_{u∈W} f_u. *Proof.* Sum Lemma 2.1 over W; streams with both ends in W cancel. ∎

Exact rational check on random unit graphs (arbitrary subsets W, arbitrary θ): `scripts/unit_identities_exact.py`.

**Lemma 2.3 (sections).** For a section with flows (V, L, w), at every interface between two consecutive units of
the chain X, 1, …, n, Y (vapor V going up, liquid L going down),
  φ(vapor up)(θ) − φ(liquid down)(θ) = V − Σ_i α_i w_i/(α_i − θ) =: Ψ_s(θ),
independent of the interface (result 12's Λ_s). *Proof.* φ(vapor) − φ(liquid) = V − Σ_i α_i (V y_i − L x_i)/(α_i − θ)
and V y − L x = w at every interface (Section 1.3). ∎

**Corollary 2.4 (section form of the cut identity).** Let W be a set of junction units, and add to W the stages of
every section whose two ends lie in W (stages of the other sections stay outside). Then
  Σ_{u∈W} Q_u = Σ_{s: bottom end in W, top end outside} Ψ_s − Σ_{s: top end in W, bottom end outside} Ψ_s
               + Σ_{other streams W→W^c} φ − Σ_{other streams W^c→W} φ + Σ_{P from W} φ_P − φ_F f(W).
*Proof.* Proposition 2.2. A section with bottom end X ∈ W and top end outside crosses the cut at the interface X/1:
vapor leaves W (+φ_vapor), liquid enters W (−φ_liquid), together Ψ_s by Lemma 2.3. A section with top end Y ∈ W and
bottom end outside crosses at the interface n/Y: vapor enters W, liquid leaves W, together −Ψ_s. ∎
With W = the complement of result 12's light region (heavy pure liquid products give Σ_{P from W} φ_P = H_j) and
θ = θ_j this is result 12's Proposition cut, now for arbitrary networks (several feed stages, splitters, recycles,
bypasses, non-section transfer streams, products of either phase, feed splitting and conditioning).

## 3. Sign uniformity

**Lemma 3.1 (equilibrium outflows).** For a liquid of flow l and composition x and the vapor of flow v in
equilibrium with it (y = K(x) x), at every θ > 0:
  φ_liquid = −l G_x(θ),   φ_vapor = −(θ v / S(x)) G_x(θ).
*Proof.* φ_liquid = −Σ α_i l x_i/(α_i − θ). For the vapor, Σ_i y_i = 1 and Σ_i α_i y_i/(α_i − θ)
= (1/S) Σ_i α_i² x_i/(α_i − θ) = (1/S)(Σ_i α_i x_i + θ Σ_i α_i x_i/(α_i − θ)) = 1 + θ G_x/S, using
α_i² = α_i(α_i − θ) + θ α_i; so φ_vapor = v − v(1 + θ G_x/S) = −θ v G_x/S. ∎ (Symbolic check for N = 2..6:
`scripts/unit_identities_exact.py`.)

**Lemma 3.2 (every unit is sign-uniform).** For θ > 0, θ ∉ {α_i}, the outflows of every unit of 𝒜 either all have
φ ≥ 0 (the unit is called *heavy-ish at θ*, set 𝒰₊) or all have φ < 0 or zero flow (*light-ish*, set 𝒰₋).
*Proof.* (E): by Lemma 3.1 both outflow phases have φ of the sign of −G_{x_u}(θ) (θ > 0, S > 0), with value 0 only if
G_{x_u}(θ) = 0 or the flow is 0. (M), (C), (R): all outflows have the same phase and composition, so φ is a positive
multiple of one number. ∎

*Remark.* The sign is not a property of the composition alone: a liquid of composition x has the sign of −G_x(θ),
a vapor of the same composition (outlet of an (R) unit) the sign of 1 − G_x(θ) (Lemma 2.0). This is why total
condensers and vaporizers are units of their own, and why the only possible loophole is a device releasing two
streams of different compositions not in phase equilibrium (Section 1.2).

## 4. The bound

**Theorem 4.1.** At every steady state of every arrangement of 𝒜, for every θ > 0, θ ∉ {α_i}:
  Σ_{u∈𝒰₊} Q_u = Σ_P max(0, φ_P(θ)) − φ_F(θ) f(𝒰₊) + Σ_{e: 𝒰₊→𝒰₋} φ_e + Σ_{e: 𝒰₋→𝒰₊} (−φ_e),
and every term of the last two sums is ≥ 0. Consequently
  R ≥ Σ_{u∈𝒰₊} max(Q_u, 0) ≥ Σ_P max(0, φ_P(θ)) − max(0, φ_F(θ)).
*Proof.* Proposition 2.2 with W = 𝒰₊. Streams leaving 𝒰₊ have φ ≥ 0 (Lemma 3.2), streams leaving 𝒰₋ have φ ≤ 0, so
the crossing terms are ≥ 0. A product stream leaving a 𝒰₊ unit has φ_P ≥ 0, one leaving a 𝒰₋ unit has φ_P ≤ 0
(< 0 or zero flow), hence Σ_{P from 𝒰₊} φ_P = Σ_P max(0, φ_P). Finally 0 ≤ f(𝒰₊) ≤ 1 and R = Σ_u max(Q_u, 0) ≥
Σ_{𝒰₊} max(Q_u, 0) ≥ Σ_{𝒰₊} Q_u. ∎

Note the stronger form: the vapor generated in the heavy-ish units alone is at least Σ_P max(0, φ_P(θ)) −
max(0, φ_F(θ)).

**Corollary 4.2 (single feed, θ = θ_j, purity).** φ_F(θ_j) = 0 by definition of θ_j, so
  R ≥ Σ_P max(0, φ_P(θ_j)).
If every product stream P is assigned to a component p(P) and has purity ≥ 1 − ε (mole fraction of p(P)), then
  R ≥ H_j − c_j ε,   c_j = max( α_j/(α_j − θ_j), α_{j+1}/(θ_j − α_{j+1}) ),
for every j, every arrangement and every steady state. For pure products R ≥ H_j.
*Proof.* Keep only the products assigned to heavy components (p(P) > j) and use max(0, φ) ≥ φ and Lemma 2.0
(φ_vapor ≥ φ_liquid):
  Σ_P max(0, φ_P) ≥ Σ_{p(P)>j} φ_P^{liq} = Σ_{i>j} α_i n^h_i/(θ_j − α_i) − Σ_{i≤j} α_i n^h_i/(α_i − θ_j),
where n^h_i is the flow of component i in the heavy-assigned products. Steady state: Σ_P n_{P,i} = z_i, so for i > j,
n^h_i = z_i − (flow of i in light-assigned products). Hence
  Σ_P max(0, φ_P) ≥ H_j − Σ_{i>j} (α_i/(θ_j − α_i)) m_i − Σ_{i≤j} (α_i/(α_i − θ_j)) m_i,
with m_i ≥ 0 the flow of i in products assigned to the other side; each m_i is impurity, the total impurity of all
products is ≤ ε Σ_P (flow of P) = ε, and the coefficients are at most c_j (α_i/(θ_j − α_i) increases with α_i < θ_j,
α_i/(α_i − θ_j) decreases with α_i > θ_j). ∎ (Agent S states the same constant, sections.md (S5').)

**Corollary 4.3 (infinitely many stages).** For every family 𝒮 ⊆ 𝒜, R*(𝒮; ε) ≥ max_j H_j − ε max_j c_j, hence
R*(𝒮) = lim_{ε→0} R*(𝒮; ε) ≥ max_j H_j = V_FTC. ∎

**What the theorem does not use.** No section physics, no Underwood roots of sections, no number of stages, no sign
of net component flows (reverse flows), no assumption on which components are present where (degenerate
distributions), no restriction on the topology (several feed stages per column, intermediate sections, stream
splitters, remixing of different streams at one stage, recycles, vapor or liquid bypasses, non-interval streams,
splits in which a product keeps every component), no restriction on heat duties (heaters and coolers at any unit,
including non-adiabatic sections with external utilities), products of either phase and any purity. All gaps (1)–(3)
of result 12 (degenerate points, reverse net flows, arrangements outside the class of two-section column networks)
are therefore closed for the lower bound.

## 5. Exact slack and the tie structure

**Proposition 5.1 (exact slack).** At θ = θ_j (single feed),
  R − Σ_P max(0, φ_P(θ_j)) = R(𝒰₋) + Q(𝒰₊) + Σ_{e: 𝒰₊→𝒰₋} φ_e + Σ_{e: 𝒰₋→𝒰₊} (−φ_e),
where R(𝒰₋) = Σ_{u∈𝒰₋} max(Q_u, 0) (heating of light-ish units) and Q(𝒰₊) = Σ_{u∈𝒰₊} max(−Q_u, 0) (cooling of
heavy-ish units); all four terms are ≥ 0. *Proof.* R = Σ_{𝒰₋} max(Q, 0) + Σ_{𝒰₊} Q + Σ_{𝒰₊} max(−Q, 0) and Theorem
4.1. ∎ (Agent S, sections.md (S4), has the same identity.)

**Corollary 5.2 (necessary conditions for approaching the bound).** Let s_n be steady states (of any members of 𝒜)
with purity ≥ 1 − ε_n, ε_n → 0, and R(s_n) − H_j → 0. Then, at θ_j:
(a) the heating of light-ish units and the cooling of heavy-ish units tend to 0;
(b) every stream between a light-ish and a heavy-ish unit (in either direction) has flux φ_e(θ_j) → 0; for a
    stream leaving an (E) unit this means flow × G_x(θ_j) → 0 (Lemma 3.1): its flow, or the G-value of its source
    composition, tends to 0;
(c) for every section in which the side changes (some interface with a heavy-ish unit below and a light-ish unit
    above, or the reverse), Ψ_s(θ_j) → 0: θ_j is asymptotically an Underwood root of that section, and at the
    interface where the side changes both G-values tend to 0 (the change happens inside a zone where the liquid
    composition is close to the hyperplane G_x(θ_j) = 0).
*Proof.* By Corollary 4.2 the left side of Proposition 5.1 is ≤ R − H_j + c_j ε_n → 0, so each nonnegative term on
the right tends to 0. (b) is the crossing term; for (c), at an interface with heavy-ish stage X below and light-ish
stage Y above, the two crossing terms are φ(vapor up) = −(θ_j V/S_X) G_X ≥ 0 and −φ(liquid down) = L G_Y > 0, whose
sum is Ψ_s(θ_j) (Lemma 2.3); both tend to 0 (similarly for the reverse case with sum −Ψ_s). ∎

**Corollary 5.3 (where heat exchangers may sit in an optimal arrangement).** Suppose the highest peak H_{j*} is
strict. If an arrangement topology 𝒮 attains R*(𝒮) = V_FTC along steady states s_n as in Corollary 5.2, then along
s_n the vapor generated at units that are light-ish at θ_{j*} and the vapor condensed at units heavy-ish at θ_{j*}
tend to 0. In particular a reboiler (or heater) whose outflows are, in the limit, a pure light product p ≤ j*, or a
condenser (cooler) whose outflow is a pure heavy product p > j*, must have vanishing duty: components 1..j* are
withdrawn through condensers or side draws, components j*+1..N through reboilers or side draws, at the optimum.
(Pure light products are light-ish: G = α_p/(α_p − θ) > 0; pure heavy products heavy-ish, both as liquids and, by
Lemma 2.0 and the vapor formula, as vapors.) This is the stage-level form of result 12's Corollary tie; it holds in
the whole class 𝒜.

**Corollary 5.4 (finite stages never tie).** In a steady state with finitely many stages and some section whose
liquid changes side with L > 0 at the changing interface, R > Σ_P max(0, φ_P(θ_j)) strictly. (The liquid going down
from the light-ish unit has −φ = L G > 0.) The bound is reached only in the limit, consistently with F's
observation (Agent F's own solver) that a binary column with 80+80 stages below minimum reflux has R equal to the
finite-purity bound to 1e-11.

## 6. Sections, the sign potential, and result 12

**Lemma 6.1 (Underwood's identity at any θ).** In a section with V, L > 0, at the interface between consecutive
units a (below, sending vapor up) and b (above, sending liquid down), with a an (E) unit,
  Ψ_s(θ) = L G_{x_b}(θ) − (θ V / S(x_a)) G_{x_a}(θ).
At the roots of Ψ_s this is Halvorsen's eq. (2.49) (thesis §2.5.2). *Proof.* Lemma 2.3 and Lemma 3.1. ∎ (Exact
check on random stage-to-stage profiles with net flows of both signs: `scripts/section_identity_exact.py`.)

**Corollary 6.2 (one change of side per section).** Along a section (bottom-end unit, stages, top-end unit; every
unit sign-uniform, V, L > 0), going up, the side (light-ish / heavy-ish) changes at most once: from heavy-ish to
light-ish if Ψ_s(θ) > 0, from light-ish to heavy-ish if Ψ_s(θ) < 0, never if Ψ_s(θ) = 0. *Proof.* At an interface
with heavy-ish a below and light-ish b above, Ψ_s = φ(vapor from a) − φ(liquid from b) > 0 (first term ≥ 0, second < 0
since L > 0); with light-ish a below and heavy-ish b above, Ψ_s < 0. Ψ_s is the same at every interface (Lemma 2.3).
This uses only sign uniformity of the end units, so it holds also when X or Y is an (M), (C) or (R) unit. ∎

**Proposition 6.3 (sign potential).** Fix θ = θ_j and a steady state. Put σ(u) = + for light-ish units, σ(u) = − for
heavy-ish units; product outlets inherit the sign of their stream (pure light products +, pure heavy products −, also
at purity ≥ 1 − ε for ε small). Define the sign graph: for every stream e: a → b an arc a → b if φ_e > 0 and b → a if
φ_e < 0; for every section s with bottom end X and top end Y an arc X → Y if Ψ_s > 0 and Y → X if Ψ_s < 0 (result 12's
graph when the nodes are the junctions). Then no arc goes from a + node to a − node.
*Proof.* Stream e: a → b with φ_e > 0: the source a has σ = − (Lemma 3.2), so the arc starts at a − node; with φ_e < 0
the arc b → a ends at the source a, which has σ = +. Sections: by Corollary 6.2, Ψ_s > 0 excludes (σ(X), σ(Y)) =
(+, −) and Ψ_s < 0 excludes (σ(Y), σ(X)) = (+, −). ∎

**Corollary 6.4 (result 12's reachability argument at every physical point).** In every steady state of every
arrangement of 𝒜, the set of nodes reachable in the sign graph from the light products contains only + nodes, so it
contains no heavy product, and choosing U = this set (or U = 𝒰₋) in the cut identity makes every term nonnegative. In
result 12 this was established, in the strict Underwood model only, from the column Lemmas sign and interior; at a
physical point it is automatic. (The node sign is the invariant of the general reachability argument: the partition is
U = {light-ish units}, the graph is the unit graph, no order on roots and no induction on junctions is needed.)

### 6.4 What the permissive counterexample violates

At the verifier's certified permissive point (N = 5, basic configuration 161, q = 0; `verification/
permissive_counterexample_N5.json`), the Underwood fluxes at θ_3 ≈ 2.53005 are (`scripts/permissive_sign_path.py`):

| column | Λ_T (top, to distillate) | Λ_B (bottom, from residue) |
|---|---|---|
| ABCDE | 1.0e-6 | 1.0e-6 |
| ABCD | 2.7e-6 | 1.7e-6 |
| BCDE | +0.005946 | +0.005947 |
| ABC | 0.772 | 0.772 |
| BCD | +0.002452 | **−0.003492** |
| CDE | **−0.002711** | +0.003236 |
| AB, BC, DE | > 0 | > 0 |

The sign graph has exactly one light → heavy path: **C → CDE → BCDE → BCD → D** (flux flows down the top section of
CDE from the light product C, up the bottom and the top sections of BCDE, which carries the excess vapor e ≈ 0.006 in
both sections, and down the bottom section of BCD into the heavy product D). By Proposition 6.3 the signs along it
would have to be + (C), + (CDE), + (BCDE), + (BCD), − (D) with a final + → − arc, which is impossible. So no family
of finite-stage steady states of this configuration can converge (flows and fluxes) to this point; equivalently, any
physical approximation satisfies R ≥ H_3 − c_3 ε (Corollary 4.2), whereas the point has R = H_3 − 2.6·10⁻⁴.

**What exactly the permissive model lacks.** The permissive Underwood model describes each column by its fluxes at
its own common roots only; it has no junction compositions. Physics adds, at every junction, a composition whose
G-sign is shared by all streams leaving it (Lemma 3.2), and the section rule of Corollary 6.2. These are global
constraints: the violation above is not visible in any single column (BCDE with both fluxes positive is locally
allowed: it only forces σ(CDE-node) ≤ σ(BCDE) ≤ σ(BCD-node)); it appears along a chain of four sections. Result 12's
strict model enforces a column-local sufficient condition (Lemma interior: no column of type (d) passes flux upward
through both sections), which happens to cut every such chain in two-section column networks; the stage model needs
no such condition.

### 6.5 A minimal repair of the permissive model (remark for global optimization)

Proposition 6.3 gives a necessary condition for an Underwood-level point (section vapor flows and net flows, link
duties) to be a limit of physical steady states, which can be imposed in a global-optimization model without
removing the physically legitimate excess vapor: **sign consistency** — for every j there are junction signs
σ_j(X) ∈ {+, −}, + at the light products and − at the heavy products, such that every section with Λ_s(θ_j) > 0 has
(σ_j(bottom end), σ_j(top end)) ≠ (+, −), every section with Λ_s(θ_j) < 0 has (σ_j(bottom), σ_j(top)) ≠ (−, +), and
every non-section transfer stream (e.g. a link through a condenser or a reboiler, modelled as a stream leaving its
source junction) has φ(θ_j) ≥ 0 if σ_j(source) = − and φ(θ_j) ≤ 0 if σ_j(source) = +. (Disjunctive constraints: one
binary variable per junction and per j.)
*Sufficiency.* Every sign-consistent point satisfies R ≥ H_j for every j: Corollary 2.4 with W = {σ_j = −} has
every crossing term ≥ 0, as in Theorem 4.1.
*Necessity (when each junction of the Underwood description is realized by one equilibrium stage, as in result 12's
class: a feed or side-draw stage).* If finite-stage steady states s_n converge to the point (section flows and net
flows, hence the Λ_s(θ_j) and the link fluxes), take σ_j(X) = the sign of the junction stage X in s_n for n large
along a subsequence on which these signs are constant (finitely many junctions). By Proposition 6.3 every section
and link with nonzero limit flux has, for n large, a flux of the same strict sign in s_n, hence satisfies the
constraint; links with zero limit flux impose no constraint. Pure light (heavy) products are + (−).
The permissive counterexample violates sign consistency (Section 6.4); result 12's strict model satisfies it at its
feasible points (its reachability set gives the signs). So "permissive model + sign consistency" is a relaxation of
physics (it keeps the excess vapor of TG wherever it is physically admissible at the level of signs) in which the FTC
is optimal, unlike the permissive model. This is a modelling remark; nothing else in this file depends on it.

## 7. Extensions

**7.1 Several process feeds.** If the arrangement has feeds F_1, …, F_m (compositions z^(k), thermal states q_k,
amounts f_k), Theorem 4.1 gives, for every θ > 0 not in {α_i},
  R ≥ Σ_P max(0, φ_P(θ)) − Σ_k f_k max(0, φ_{F_k}(θ)),
and hence R ≥ sup_θ of the right side. For one feed and pure products the supremum is max_j H_j (on I_j the right side
is Σ_{i≤j} α_i z_i/(α_i − θ) − (1 − q), increasing, for θ ≤ θ_j, and H(θ) = Σ_{i>j} α_i z_i/(θ − α_i), decreasing, for
θ ≥ θ_j; for θ < α_N it is ≤ 0, and for θ > α_1 it is ≤ max(0, q − 1) < max_j H_j, since H_j = P_j^T + (q − 1)),
so the one-feed bound is not improved by other θ. Whether the multi-feed bound is attained is not studied.

**7.2 Non-sharp products.** For arbitrary product specifications the bound R ≥ max_j Σ_P max(0, φ_P(θ_j)) holds as it
stands. For a two-product column it is Underwood's classical bound V_B ≥ Σ_i α_i b_i/(θ_j − α_i) at the feed roots,
here proved for any number of stages and any arrangement.

**7.3 Vapor products.** A product withdrawn as vapor contributes φ^vap = φ^liq + flow (Lemma 2.0): heavy vapor
products raise the bound, light ones do not enter (they come from light-ish units).

## 8. Checks (scripts in `total/network/scripts/`, run from that folder under WSL)

| Script | What it checks | Arithmetic | Result |
|---|---|---|---|
| `section_identity_exact.py` | Lemma 6.1 and Corollary 6.2 on 400 random stage-to-stage profiles (2–5 components, net flows of both signs), 13490 interface × θ tests | exact rationals | PASS (~1 s) |
| `unit_identities_exact.py` | Proposition 2.2 on random unit graphs with recycles (1798 graph × θ × subset tests); Lemma 3.1 symbolically for N = 2..6 | exact rationals, sympy | PASS (~6 s) |
| `cut_bound_check.py 40 1` | MESH steady states (own solver `stagenet.py`, `secnet.py`): random section/junction networks (splitters, mixing, recycles, split feeds, vapor products, heaters/coolers), ternary Petlyuk grid (8–40 stages per section): identity residual ≤ 3e-12, crossing signs, R ≥ H'_j | floating point | PASS (~3 min) |
| `permissive_sign_path.py` | the light → heavy sign path of the permissive point (Section 6.4) | mpmath, 40 digits | path found |
| `ftc_fixed_flows.py N δ n-list [seed]` | finite-stage FTC at fixed flows F_δ (Section 10.3): impurities vs n | floating point | evidence |
| `ftc_neutral_mode.py N δ n-list [seed]` | smallest singular value of the steady-state Jacobian vs n (Section 10.4) | floating point | evidence |
| `ftc3_structure.py δ n-list`, `ftc3_forced.py δ η_A η_C n-list` | N = 3 junction ratios, net flows, traces; forced-impurity closure (Section 11) | floating point | evidence |

Independent numerical checks of the same bound by Agents S (`total/sections/scripts/`), F (`total/foundations/
scripts/unitnet.py`) and P (`total/simulation/scripts/test_stage_bound.py`) all PASS.

## 9. Not covered here

- **Attainment** R*(𝒮_Petlyuk) ≤ max_j H_j in the finite-stage limit sense (HS's own formula for the Petlyuk minimum).
  Given Corollary 4.3 it implies the conjecture (it is sufficient, not logically necessary: the conjecture compares
  R*(𝒮_Petlyuk) with every R*(𝒮)). Proved in the infinite-stage Underwood model (result 12, Theorem achieve); for N = 2 in
  the finite-stage model by Agent S (sections.md, Prop. 10); general N: Agent Astra (`total/attainment-astra/`),
  numerics: Agent P. Section 10 below gives the network-level ingredients: existence of steady states for prescribed
  flows (so that R is fixed by the flows and attainment reduces to purity) and two numerical observations.
- Devices excluded by HS (Section 1.2, last row; class.md §7).
- The bound concerns total vaporization R (= total heat input under CMO with the heat of vaporization as unit). With
  internal heat exchange counted in R the bound still holds; HS exclude internal heat exchange because the *external*
  heat can then be lower.

## 10. Steady states exist for prescribed flows; attainment reduces to purity

This section serves the attainment work (S, Astra). It is about existence, not about the bound.

**Setting.** Every unit type of Section 1.1 is an (E) unit with prescribed outflows: an (M) liquid mixer/splitter or
a (C) unit is an (E) unit with zero vapor outflow (its liquid has the mixed composition); an (M) vapor unit or an (R)
unit is an (E) unit with zero liquid outflow (its vapor has the mixed composition; x is then the "virtual" liquid in
equilibrium with it, x ∝ y/α). An **operating specification** of an arrangement prescribes the flow of every stream
(internal streams and product streams) and the feed portions f_k, such that every unit satisfies the total mass
balance Σ_in + f_k = Σ_out; the duties are then Q_k = V_out,k − V_in,k − (1 − q) f_k (CMO), and R = Σ_k max(Q_k, 0)
**depends on the flows only**. Unknowns: one composition x_k ∈ Δ (open simplex) per unit. Steady-state equations:
for every unit k and component i,
  F_{k,i}(x) := Σ_{k'} (ℓ_{k'k} x_{k',i} + v_{k'k} K_i(x_{k'}) x_{k',i}) + f_k z_i − (L_k + V_k K_i(x_k)) x_{k,i} = 0,
with ℓ_{k'k}, v_{k'k} the liquid and vapor flows from k' to k, L_k, V_k the total liquid and vapor outflows of k
(products included) and K_i(x) = α_i / Σ_j α_j x_j.
Hypothesis (G): in the directed graph of streams with positive flow, every unit is reachable from a unit receiving
feed, and from every unit a unit delivering a product can be reached.

**Theorem 10.1 (existence).** Under (G), for every α_1 > … > α_N > 0 the equations F = 0 have a solution with all
mole fractions positive. Every solution satisfies c ≤ x_{k,i} ≤ 1 for a constant c > 0 depending only on the flows,
z and κ := α_1/α_N.

*Proof.* (i) *Linear structure.* For fixed values K_{k,i} > 0 the equations for component i read −A_i(K) X_{·,i} =
f z_i, where −A_i(K) has diagonal entries L_k + V_k K_{k,i} > 0 and off-diagonal entries −(ℓ_{k'k} + v_{k'k} K_{k',i})
≤ 0 (row k, column k'). The excess of column k' (diagonal minus the sum of the absolute off-diagonal entries) is the
flow sent by k' to product outlets (liquid, plus vapor times K_{k',i}), ≥ 0, and by (G) every column is joined, along
the graph of nonzero off-diagonal entries, to a column with positive excess. So −A_i(K) is a weakly chained
diagonally dominant Z-matrix with positive diagonal, hence a nonsingular M-matrix: (−A_i(K))^{-1} ≥ 0, and its (k, m)
entry is > 0 exactly when k can be reached from m. Hence any solution with positive x has
X_{·,i} = (−A_i(K))^{-1} f z_i, positive at every unit (reachability from a feed unit).
(ii) *Normalization is automatic.* Work with unnormalized X ∈ (0, ∞)^{n×N}, K_i(X_k) := α_i / S(X_k),
S(X_k) := Σ_j α_j X_{k,j} / Σ_j X_{k,j}. Summing the equations of unit k over i and using Σ_i K_i(X_k) X_{k,i} = σ_k :=
Σ_i X_{k,i} gives Σ_{k'} (ℓ_{k'k} + v_{k'k}) σ_{k'} + f_k − (L_k + V_k) σ_k = 0: the same M-matrix system with K ≡ 1,
whose unique solution is σ ≡ 1 by the total mass balances. So every zero X of F is a steady state (Σ_i X_{k,i} = 1).
(iii) *A priori bounds.* For every composition, S ∈ [α_N, α_1], so K_{k,i} ∈ [κ^{-1}, κ]. The positive entries of
(−A_i(K))^{-1} are continuous in K on this compact box, hence bounded below by some c_0 > 0; therefore every zero
satisfies X_{k,i} ≥ c := c_0 min_{feed units} f_k · min_i z_i, and X_{k,i} ≤ σ_k = 1.
(iv) *Degree.* Consider the homotopy α(t) = (α_1^t, …, α_N^t), t ∈ [0, 1]; κ(t) ≤ κ, so (iii) holds uniformly in t.
Let Ω = {X : c/2 < X_{k,i} < 2}. F(·; t) is continuous on Ω̄ × [0, 1] and, by (ii)–(iii), has no zero on ∂Ω. At
t = 0, K ≡ 1, the system is linear with the nonsingular block-diagonal matrix diag(A, …, A) (A = A_i(1)), so its unique
zero X_{k,i} = z_i is nondegenerate and deg(F(·; 0), Ω, 0) = (sign det A)^N ≠ 0. By homotopy invariance of the
Brouwer degree, deg(F(·; 1), Ω, 0) ≠ 0, so F(·; 1) has a zero in Ω. ∎

**Corollary 10.2 (attainment reduces to purity).** For an arrangement topology 𝒮 (any stage numbers),
R*(𝒮) ≤ ρ if and only if for every η > 0 there are stage numbers and an operating specification with R ≤ ρ + η
for which *some* steady state has all purities ≥ 1 − η. Steady states exist for every specification (Theorem 10.1),
so only the purity of the steady states with given flows is at stake: for the FTC it suffices to exhibit flows F_δ
with R(F_δ) = max_j H_j + δ and to show that, for these fixed flows, steady states with purities → 1 exist as the
stage numbers grow (step (C) of the plan; S's Lemma 11 gives the limit relations of long sections; Astra works on the
general case).

**10.3 Numerical evidence for step (C)** (`scripts/ftc_fixed_flows.py`; floating point). Flows F_δ: result 12's
attainment point (internal columns at their preferred split, binary columns with excess e_r = max H − H_r) plus a
vapor/liquid circulation of size δ from the reboiler to the condenser spread over all monotone paths of the FTC grid
(every section gets a margin; all net flows and junction balances unchanged; R = max_j H_j + δ exactly). Then the
flows are FIXED and the stage number n per section grows:

| case | δ | n | product impurities (1 − x_p, products 1..N) |
|---|---|---|---|
| N = 3, α = (4, 2, 1), z = 1/3, q = 1 | 0.02 | 10 / 20 / 40 / 80 | (1.6e-2, 5.0e-2, 3.5e-2) / (1.6e-3, 9.4e-3, 7.8e-3) / (1.4e-6, 1.8e-4, 1.8e-4) / (6e-15, 1.9e-8, 8e-15) |
| same | 0.002 | 20 / 40 / 80 / 160 | (2.6e-3, 1.4e-2, 1.1e-2) / (2.2e-5, 1.1e-3, 1.1e-3) / (1e-14, 4.5e-7, 1e-14) / (0, 1.3e-14, 0) |
| N = 4, seed 1 (α ≈ (10.57, 6.48, 2.22, 1), q = 0.54) | 0.01 | 10 / 40 / 80 / 120 | (7.8e-2, 2.7e-2, 3.1e-2, 3.2e-2) / (3.0e-3, 9.3e-4, 7.3e-5, 2.3e-8) / (4e-10, 5.9e-5, 9e-10, 0) / (0, 2.1e-6, 1e-14, 0) |
| N = 5, seed 2 | 0.01 | 10 / 20 / 40 | (4.8e-2, 0.16, 2.7e-2, 4.4e-2, 5.8e-2) / (1.8e-2, 5.7e-2, 3.8e-3, 1.1e-2, 1.9e-2) / (3.5e-3, 1.1e-2, 2.5e-5, 1.2e-3, 2.6e-3) |

For N = 3 the Newton residuals are ≤ 1e-13 (n ≤ 40) and ≤ 2e-8 (n = 80); for N = 4 the residuals are 1e-5–5e-7 at
n ≥ 20 and for N = 5 my solver does not converge beyond n = 40 (see 10.4): those rows are indicative only. Agent P's
optimization runs are the reference numerics.

**10.4 A neutral mode of the limit system (caution for the transversality route).** At these fixed flows the smallest
singular value of the steady-state Jacobian decreases geometrically with n (`scripts/ftc_neutral_mode.py`, N = 3,
δ = 0.02: 7.0e-3, 1.2e-3, 3.2e-4, 9.0e-5 at n = 5, 10, 15, 20), and the singular vector is an exchange A ↔ B
concentrated in the middle of the top section of the binary column [A, B], which carries excess vapor (e_1 = H_2 −
H_1 = 0.29); for N = 5 (seed 2, n = 15) the singular value is 4.5e-6, the vector a B ↔ C exchange in the bottom
section of the binary column [B, C] (excess 0.50). In the N = 3 case at n = 20 the profile of that top section
lingers near the junction composition (x_A = 0.337 at the junction, then 0.348, 0.364, 0.387, …), crosses over in the
middle (x_A from 0.46 to 0.92 between stages 5 and 13) and approaches pure A at the condenser: the section is
pinched at both ends and carries a composition front in between; 78% of the weight of the singular vector lies in this
section. *Interpretation (not proved):* the position of such a front is fixed only by exponentially small terms
(how far the junction composition is from the section's pinch), so in the infinite-stage limit it is a free
parameter, and the limit (pinch-matching) system has a zero mode for every excess-carrying section that is pinched at
both ends. An implicit-function argument must then add these positions as variables fixed by the exponentially small
terms, or avoid them (monotonicity and intermediate-value arguments, as in S's binary proof), or use flows in which no
excess-carrying section is pinched at both ends.

**10.5 Fixed-flow N = 3 data to n = 256** (`ftc_fixed_flows_N3.log`; α = (4, 2, 1), z = 1/3, q = 1, R = max H + δ
fixed; residuals ≤ 2e-12 except 5e-7 / 4e-6 at n = 64): impurities (distillate, side product, bottoms)

| δ | n = 32 | n = 64 | n = 128 | n = 256 |
|---|---|---|---|---|
| 0.05 | (1.0e-5, 3.2e-4, 3.1e-4) | (4.7e-11, 1.3e-7, 1.3e-7) | (4.6e-24, 2.0e-14, 2.8e-25) | (1.6e-47, 5.7e-28, 3.8e-52) |
| 0.02 | (3.0e-5, 1.0e-3, 9.9e-4) | (5.0e-12, 7.6e-7, 2.0e-11) | (1.2e-23, 3.2e-13, 5.9e-25) | (7.0e-47, 5.8e-26, 1.2e-51) |
| 0.005 | (1.3e-4, 2.5e-3, 2.4e-3) | (7.3e-12, 5.0e-6, 3.2e-11) | (2.0e-23, 3.7e-12, 9.4e-25) | (1.8e-46, 2.0e-24, ~0) |

## 11. Closing the ternary Petlyuk network: plan (SKETCHED, not proved)

Status: only Theorem 10.1 (existence for prescribed flows) is proved. What follows is a closure plan for N = 3 with
numerical support; it is **not** a proof of attainment.

**11.1 Why the naive limit is degenerate (observation, consistent with Astra's rank count and §10.4).** At fixed
flows with sharp products, the infinite-stage limit has all net flows fixed by the flows (prefractionator sharp A/C
split with d_B = D_pre − z_A), F = Astra's X, J2 = pure B, but **J1 is any point of the AB face and J3 any point of the
BC face**: a two-parameter family. Reason (S's Lemma 11): the top section of the binary column [A, B] is pinched only
at the condenser end, its bottom section only at the side-draw end, and the prefractionator top accepts any AB return
(Astra's Theorem 2). Numerically ξ = A/(A+B) at J1 drifts with n (0.364, 0.337, 0.328, 0.319, 0.313 at n = 16..64,
δ = 0.05): the A/B front of the binary column sits at a fixed FRACTION of its total length (where the distillate's B
leak and the side product's A leak balance), not at J1.

**11.2 Pinning the fronts by prescribed forced impurities (the key idea).** Take product rates D = z_A − η_A,
B = z_C − η_C, S = z_B + η_A + η_C (side product collects an O(η) impurity by design). Then A flows down the bottom
section M2 of [A, B] at net rate ≈ η_A, so M2's largest Underwood root φ_A' lies just above α_A and its eigenvector
e_A' ∝ (η_A/(φ_A' − α_A), d_B/(φ_A' − α_B), ε_C/(φ_A' − α_C)) is a positive AB pinch with O(1) A content (as η_A → 0,
η_A/(φ_A' − α_A) → (V_2 − α_B d_B/(α_A − α_B))/α_A > 0, positive exactly when A rises in M2, L_2/V_2 < α_A/α_B).
Going up M2 the profile converges to e_A', so **J1 → e_A'(M2) is pinned**; symmetrically J3 → e_C'(M3). With
η_A, η_C > 0 fixed, every net flow has its natural sign: PT ≥ 0, PB ≤ 0, M1 ≥ 0, M2 ≤ 0, M3 ≥ 0, M4 ≤ 0
(ε_C = C leak up the prefractionator top and ε_A = A leak down its bottom are > 0 and exponentially small), so every
section can be propagated in a **positivity-preserving direction** (M_T^{-1} ≥ 0 going down a section with w ≥ 0,
M_B ≥ 0 going up a section with w ≤ 0): down M1 from the distillate to J1, up M2 from J2 to J1, down the
prefractionator top from J1 to F, up its bottom from J3 to F, down M3 from J2 to J3, up M4 from the reboiler to J3.
Unknowns: distillate impurities (2), bottoms impurities (2), ε_A, ε_C; J2 = (z − D x_D − B x_Reb)/S; equations:
matching at J1, J3, F (2 each). All six sections are then in S's Lemma 12 charts with nonzero pivots and bounded
coefficients (each singly pinched), and at leading order the six equations are block-triangular: J1's A:B ratio
fixes the scaled distillate B impurity, J1's C content (M1 side ∝ distillate C, M2 side ∝ ε_C, prefractionator side
∝ (α_C/u)^{n_T}) fixes the scaled distillate C impurity and ε_C (a linear 2×2 system with positive solution),
symmetrically at J3, and the F matching is Astra's Lemma 2 intersection X (positive weights). **Missing for a proof:**
the precise exponential scalings when section lengths differ, uniform C¹ bounds (Lemma 12) for the scaled system,
the invertibility of its limit Jacobian (expected from the triangular structure), then Poincaré–Miranda/IFT. Since
the side-product purity is 1 − O(η) and R = max H + δ is fixed by the flows (Theorem 10.1, Corollary 10.2), letting
n → ∞, then η → 0, then δ → 0 would give R*(𝒮_FTC) ≤ max_j H_j for N = 3.

**11.3 Numerical support for 11.2** (`scripts/ftc3_forced.py 0.05 0.003 0.003 …`; α = (4, 2, 1), z = 1/3, q = 1,
δ = 0.05, η_A = η_C = 0.003):

| n | impurities (D, S, B) | ξ at J1 | ζ at J3 | ε_C | ε_A | smallest singular value of the Jacobian |
|---|---|---|---|---|---|---|
| 16 | (2.1e-4, 2.1e-2, 3.3e-3) | 0.68629 | 0.59985 | 1.5e-3 | 2.3e-3 | 2.6e-4 |
| 24 | (6.6e-6, 1.77e-2, 2.5e-5) | 0.69613 | 0.56209 | 2.7e-4 | 1.7e-4 | 2.4e-4 |
| 32 | (2.2e-7, 1.768e-2, 4.7e-7) | 0.696212 | 0.554915 | 3.9e-5 | 8.1e-6 | 2.2e-4 |
| 48 | (2.3e-10, 1.768e-2, 2.0e-10) | 0.696156 | 0.553945 | 7.2e-7 | 1.8e-8 | 1.8e-4 |
| 64 | (2.5e-13, 1.768e-2, 8.3e-14) | 0.696154 | 0.553931 | 1.3e-8 | 3.9e-11 | — |
| 96 | (< 1e-15, 1.768e-2, < 1e-15) | 0.696154 | 0.553931 | 4.0e-12 | 1.9e-16 | — |

ξ and ζ now **converge** (pinned), the smallest singular value **stays bounded away from zero** (2e-4, versus
geometric decay without forced impurities, §10.4), distillate and bottoms purify geometrically, the side-product
impurity tends to (η_A + η_C)/S = 0.0177, and w_M2 ≈ (−3.0e-3, −0.111, −ε_C) ≤ 0, w_M3 ≈ (ε_A, 0.222, 3.0e-3) ≥ 0, as
predicted. (n = 128 did not converge numerically.) This supports the closure plan; it is evidence, not a proof.
