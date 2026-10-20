# Sections and junctions at every physical operating point (Agent S)

Agent S (sections and junctions), Claude Opus 5.5 (effort xhigh). Version 2, 2026-10-20 14:25. Paths relative to
the problem folder `problems/halvorsen-skogestad-petlyuk-minimum-vapor/`. Reference definition of the class 𝒜 and
of the physical model: `total/foundations/class.md` (Agent F); Section 1 below restates everything this file uses,
so that it can be read alone. Sources: HS3 = Halvorsen & Skogestad 2003 Part 3 (doi:10.1021/ie0108651);
Th = Halvorsen's thesis (NTNU 2001:43).

**Status of this file.** Complete proofs (no numerical step) of: the stream-flux node balance (Lemma 1), sign
uniformity of every unit (Lemma 2), the stage-level cut identity and bound R ≥ Σ_P max(0, φ_P(θ_j)) at every steady
state of every arrangement of 𝒜 (Theorem 3), the exact defect identity (Theorem 3(c)), the pure-product and
explicit finite-purity forms R ≥ H_j and R ≥ H_j − c_j ε (Corollaries 4, 5), the evaluation at any θ and several
feeds (Corollary 6), the section-level and column-level lemmas valid at every finite-stage operating point with
arbitrary signs of the net flows (Lemmas 6–9), which replace result 12's sign lemmas S3/S4/S5 and settle the
degenerate-point and reverse-flow gaps of the verifier, **attainment for two components** (Proposition 10: with
Theorem 3 the conjecture is completely proved in 𝒜 for N = 2), and the limit relation of a long section
(Lemma 11) and its uniform exponential closeness to the Underwood limit relation (Lemma 12), two tools for
attainment. What is **not** proved here: attainment for N ≥ 3 (finite-stage generalized
Petlyuk steady states with purity → 1 and R → max_j H_j); see Section 8. The same stage-level bound was found
independently by Agents N (13:38) and F (13:44); this file gives the complete proof with every unit type and the
explicit constant.

---

## 1. Model (as in class.md §2–3)

Components 1..N, α_1 > … > α_N > 0 constant; vapor–liquid equilibrium y_i = α_i x_i / S(x), S(x) = Σ_k α_k x_k.
Molar heat of vaporization 1 (constant molar overflow, CMO). Feed F = 1, composition z (z_p > 0, Σ z_p = 1),
liquid fraction q ∈ ℝ. Feed roots θ_j ∈ (α_{j+1}, α_j) of Σ_p α_p z_p/(α_p − θ) = 1 − q; peaks
H_j = Σ_{p>j} α_p z_p/(θ_j − α_p).

An **arrangement** is a finite directed multigraph of *units* joined by *streams*; each stream e has a phase
(saturated liquid or vapor), component flows n_{e,i} ≥ 0, total flow |e| = Σ_i n_{e,i}, and vapor flow v_e = |e| for
a vapor, 0 for a liquid. Unit types (class.md §3):

- (E) equilibrium stage: any inflows; external duty Q_u ∈ ℝ; its content leaves as liquid of composition x_u and
  vapor of composition y_u, y_{u,i} = α_i x_{u,i}/S(x_u), each phase divided into any number of outflow streams
  (a phase may have zero flow). Trays, feed stages, flashes, partial condensers and reboilers, stages with
  intermediate exchangers.
- (M) single-phase mixer/splitter: inflows of one phase, outflows of that phase with the mixed composition; Q_u = 0.
- (C) total condenser: vapor inflows, liquid outflows of the mixed composition; Q_u = −(vapor in).
- (R) total vaporizer: liquid inflows, vapor outflows of the mixed composition; Q_u = +(liquid in).
- Feed portions (composition z, liquid fraction q, sizes f_k ≥ 0, Σ f_k = 1) enter units; product streams leave
  units, each assigned to a component p, in either phase.

**Steady state:** at every unit u, (i) component balances Σ_{out} n_{e,i} = Σ_{in} n_{e,i} (feed portions and
products included); (ii) CMO energy balance Σ_{out} v_e − Σ_{in} v_e = Q_u, a feed portion f contributing the vapor
inflow (1 − q) f; (iii) the composition rule of the unit type. Objective R = Σ_u max(Q_u, 0) (total heat supplied =
total vaporization); Q := Σ_u max(−Q_u, 0).

Nothing is assumed about the sign of net flows in sections, the distribution of components, the topology
(recycles allowed), the numbers of stages, stability or uniqueness of the steady state.

---

## 2. The stream flux and the node balance

For θ ∉ {α_1, …, α_N} define the **flux** of a stream e and of a feed portion of size f:

  φ_e(θ) = v_e − Σ_i α_i n_{e,i}/(α_i − θ),   φ_F(θ) = f [(1 − q) − Σ_i α_i z_i/(α_i − θ)].

(In PROGRESS.md I first wrote T_e = −φ_e.) For a liquid stream of flow l and composition x:
φ = −l G_x(θ), G_x(θ) := Σ_i α_i x_i/(α_i − θ). For a vapor stream of flow v and composition y:
φ = v − v Σ_i α_i y_i/(α_i − θ) = −θ v Σ_i y_i/(α_i − θ), using α_i/(α_i − θ) = 1 + θ/(α_i − θ) and Σ y_i = 1.
By definition of θ_j, **φ_F(θ_j) = 0** for every feed portion.

**Lemma 1 (node balance).** At every unit u of a steady state and every θ ∉ {α_i}:

  Σ_{e out of u} φ_e(θ) − Σ_{e into u} φ_e(θ) = Q_u,

where "into u" includes the feed portions entering u (with flux φ_F) and "out of u" includes the products leaving u.

*Proof.* φ is linear in (v_e, n_e): φ_e = v_e − Σ_i c_i n_{e,i} with c_i = α_i/(α_i − θ). For a feed portion, the
same formula holds with vapor (1 − q) f and component flows f z_i. Hence Σ_out φ − Σ_in φ =
(Σ_out v − Σ_in v) − Σ_i c_i (Σ_out n_i − Σ_in n_i) = Q_u − 0 by (ii) and (i). ∎

Summing Lemma 1 over all units at θ = θ_j (internal streams cancel): Σ_P φ_P(θ_j) = Σ_u Q_u = R − Q. For pure
liquid products this is the overall energy balance R − Q = Σ_p α_p z_p/(θ_j − α_p) = H_j − P^T_j = −(1 − q).

**Section flux.** An adiabatic CMO section is a chain of (E) units 1..n (bottom to top, n ≥ 0) with Q = 0, constant
vapor V flowing up and liquid L flowing down, and no other streams; the net upward flows w_i = V y_i − L x_i are the
same at every cross-section (component balances). At any cross-section, with vapor y going up and liquid x coming
down,

  Ψ_s(θ) := V − Σ_i α_i w_i/(α_i − θ) = φ(vapor up) − φ(liquid down)   (linearity),

which is result 12's Λ_s. If the vapor going up leaves an (E) stage of liquid composition x_b and the liquid coming
down has composition x_a:

  Ψ_s(θ) = L G_{x_a}(θ) − (θ V / S(x_b)) G_{x_b}(θ)          (Underwood's identity; Th §2.5.2, eq. 2.49).   (2.1)

No sign of w is used.

---

## 3. Sign uniformity of every unit

**Lemma 2 (sign uniformity).** Let θ > 0, θ ∉ {α_i}. At any steady state, for every unit u there is
σ_u ∈ {−1, 0, +1} such that every stream leaving u (products included) has flux of sign σ_u or zero flux.
Precisely:
- (E): every liquid outflow of flow l has φ = −l G_{x_u}(θ); every vapor outflow of flow v has
  φ = −(θ v/S(x_u)) G_{x_u}(θ). Hence σ_u = −sign G_{x_u}(θ).
- (M), (C), (R): all outflows have the same phase and the same composition c, hence fluxes proportional to their
  flows with a common factor (−G_c(θ) for liquids, −θ Σ_i c_i/(α_i − θ) for vapors).

*Proof.* (E): for the vapor, Σ_i y_i/(α_i − θ) = Σ_i α_i x_i/(S(x)(α_i − θ)) = G_x(θ)/S(x), and θ > 0,
S(x) > 0. (M), (C), (R): immediate. ∎

θ > 0 is the only place where α_N > 0 matters; every θ_j lies in (α_N, α_1) ⊂ (0, ∞).

*Remarks.* (a) A unit releasing a liquid and a vapor of compositions not in phase equilibrium (e.g. a hypothetical
"partial vaporizer" returning vapor and liquid of the same composition) is not sign-uniform; such a unit is not a
CMO equilibrium-stage object. Every physical heat exchanger is (E) (phases leave in equilibrium), (C) or (R) (single
outflow phase). (b) **Murphree stages are in 𝒜.** A tray with Murphree vapor efficiency E ∈ [0, 1],
y_n = y_{n−1} + E (K(x_n) x_n − y_{n−1}), is exactly the combination: (M) splitter sending the fraction E of the
vapor from below to an (E) unit and 1 − E to a bypass; the (E) unit (liquid from above, vapor E V from below,
outflows L and E V); an (M) mixer of its vapor with the bypass. Indeed the mixed vapor is
E K(x_n) x_n + (1 − E) y_{n−1} and the balance L x_{n+1} + V y_{n−1} = L x_n + V y_n holds. Murphree liquid
efficiency in [0, 1]: the same with a liquid bypass. (Efficiencies above 1 and rate-based models are not covered.)

---

## 4. The stage-level cut theorem

Fix j and put θ = θ_j. Let 𝒰₊ be the set of units with σ_u ≥ 0 (all outflows have φ ≥ 0, "heavy-ish") and
𝒰₋ the units with σ_u = −1 ("light-ish"). By Lemma 2 this is a partition, and:

- every stream from a unit of 𝒰₊ (to a unit of 𝒰₋ or to the outside) has φ ≥ 0;
- every stream from a unit of 𝒰₋ has φ ≤ 0.

**Theorem 3 (stage-level cut).** At every steady state of every arrangement of 𝒜, for every j:

(a) R(𝒰₊) − Q(𝒰₊) = Σ_{P from 𝒰₊} φ_P(θ_j) + Σ_{e: 𝒰₊→𝒰₋} φ_e(θ_j) − Σ_{e: 𝒰₋→𝒰₊} φ_e(θ_j),
    where R(𝒲) = Σ_{u∈𝒲} max(Q_u, 0) and Q(𝒲) = Σ_{u∈𝒲} max(−Q_u, 0);

(b) R ≥ Σ_{P} max(0, φ_P(θ_j));

(c) (exact defect) R − Σ_P max(0, φ_P(θ_j)) = R(𝒰₋) + Q(𝒰₊) + Σ_{e crossing between 𝒰₊ and 𝒰₋} |φ_e(θ_j)|.

*Proof.* (a) Sum Lemma 1 over u ∈ 𝒰₊. Streams with both ends in 𝒰₊ cancel; streams from 𝒰₊ to 𝒰₋ and products
leaving 𝒰₊ appear with +φ; streams from 𝒰₋ into 𝒰₊ with −φ; feed portions entering 𝒰₊ with −φ_F(θ_j) = 0; the
right-hand side is Σ_{u∈𝒰₊} Q_u = R(𝒰₊) − Q(𝒰₊).
(c) Products leaving 𝒰₊ have φ_P ≥ 0 and products leaving 𝒰₋ have φ_P ≤ 0 (Lemma 2), so
Σ_{P from 𝒰₊} φ_P = Σ_P max(0, φ_P). Streams 𝒰₊→𝒰₋ have φ_e = |φ_e|, streams 𝒰₋→𝒰₊ have −φ_e = |φ_e|. Insert in
(a) and add R(𝒰₋) to both sides (R = R(𝒰₊) + R(𝒰₋)).
(b) Every term on the right of (c) is ≥ 0. ∎

**Corollary 4 (pure products).** If every product stream is pure (it contains only its assigned component), then
R ≥ H_j for every j, i.e. R ≥ max_j H_j = V_Petlyuk.

*Proof.* A pure liquid product of component p with flow P_p has φ = α_p P_p/(θ_j − α_p): positive iff p > j. A pure
vapor product has φ = θ_j P_p/(θ_j − α_p), positive iff p > j and then ≥ α_p P_p/(θ_j − α_p) (θ_j > α_p). The total
flow of component p in its products is z_p (component balances over the whole arrangement). Hence
Σ_P max(0, φ_P) ≥ Σ_{p>j} α_p z_p/(θ_j − α_p) = H_j. ∎

**Corollary 5 (finite purity, explicit constant).** If every product stream has purity ≥ 1 − ε (its assigned
component has mole fraction ≥ 1 − ε), then for every j

  R ≥ H_j − c_j ε,   c_j = max( α_j/(α_j − θ_j), α_{j+1}/(θ_j − α_{j+1}) ).

*Proof.* Let a(P) be the component assigned to P, P_i its component flows, |P| = Σ_i P_i. Put κ = α_j/(α_j − θ_j)
and λ = α_{j+1}/(θ_j − α_{j+1}). For i ≤ j, α_i/(α_i − θ_j) ≤ κ (decreasing in α_i > θ_j); for p > j,
α_p/(θ_j − α_p) ≤ λ (increasing in α_p < θ_j). Let P be a product with a(P) = p > j.
Liquid: φ_P = Σ_i α_i P_i/(θ_j − α_i) ≥ α_p P_p/(θ_j − α_p) − κ Σ_{i≤j} P_i (the heavy impurities i > j, i ≠ p,
contribute ≥ 0). Vapor: φ_P = Σ_i θ_j P_i/(θ_j − α_i) ≥ α_p P_p/(θ_j − α_p) − κ Σ_{i≤j} P_i (for i > j:
θ_j/(θ_j − α_i) ≥ α_i/(θ_j − α_i) > 0; for i ≤ j: θ_j/(α_i − θ_j) < κ). Summing over heavy-assigned products and
dropping the light-assigned ones (max(0, ·) ≥ 0):

  Σ_P max(0, φ_P) ≥ Σ_{p>j} (α_p/(θ_j − α_p)) m_p − κ I_L,

where m_p = Σ_{a(P)=p} P_p and I_L = Σ_{a(P)>j} Σ_{i≤j} P_i (light impurities of heavy-assigned products). By the
overall component balance m_p = z_p − r_p, r_p = Σ_{a(P)≠p} P_p (component p carried as an impurity). So the right
side is ≥ H_j − κ I_L − λ Σ_{p>j} r_p. The quantities I_L (light components in heavy products) and Σ_{p>j} r_p (heavy
components in products assigned to another component) count disjoint parts of the total impurity
I = Σ_P (|P| − P_{a(P)}) ≤ ε Σ_P |P| = ε (the products add up to the feed, Σ_P |P| = 1). Hence
R ≥ H_j − max(κ, λ)(I_L + Σ_{p>j} r_p) ≥ H_j − c_j ε. ∎

**Corollary 6 (any θ, several feeds; θ_j is the best choice).** Lemma 2 and Theorem 3(a) hold for every θ > 0,
θ ∉ {α_i}, if the feed portions entering 𝒰₊(θ) are kept: R ≥ Σ_P max(0, φ_P(θ)) − Σ_F max(0, φ_F(θ)), the last sum
over all feed portions (and over several process feeds of different compositions and thermal states, if any). For
one feed and pure products, on (α_{j+1}, α_j) the right-hand side equals h(θ) := Σ_{p>j} α_p z_p/(θ − α_p)
(decreasing) for θ ≥ θ_j, where φ_F ≤ 0, and h(θ) − φ_F(θ) = Σ_{p≤j} α_p z_p/(α_p − θ) − (1 − q) (increasing) for
θ ≤ θ_j, where φ_F ≥ 0 (φ_F is decreasing and vanishes at θ_j). So its maximum over the interval is H_j, attained at
θ = θ_j: the feed roots are the optimal evaluation points of the stage-level bound.

*Proof.* In (a), feed portions entering 𝒰₊ contribute −φ_F(θ) ≥ −max(0, φ_F(θ)); the rest is unchanged. ∎

**Consequence for the conjecture (lower half).** With class.md §4, for every family 𝒮 ⊆ 𝒜 and every ε > 0,
R*(𝒮; ε) ≥ max_j (H_j − c_j ε), hence R*(𝒮) ≥ max_j H_j = V_Petlyuk. In particular the Petlyuk (FTC) arrangement
itself needs at least max_j H_j: the part of HS's formula that HS3 §5.3 left unproved ("We have not carried out a
detailed proof for the general N-component M-product case for other than the far left and right peaks") holds in
the physical model, for every N, every q and every peak. Together with attainment (Section 8), this is the
conjecture in 𝒜.

*What the proof uses, and nothing else:* (1) CMO at every unit (Lemma 1); (2) phase equilibrium with constant α at
(E) units and single outflow composition and phase at (M), (C), (R) units (Lemma 2); (3) θ_j > 0, i.e. α_N > 0;
(4) φ_F(θ_j) = 0 for every feed portion (one feed composition and thermal state, θ_j its root); (5) R counts every
positive duty (no heat delivered internally without being counted). It does not use: infinite stages, adiabatic
sections (duties may sit at any stage), signs of net flows, distributions, topology, purities (except in
Corollaries 4–5), stability.

---

## 5. Sections: physics at every operating point (brief task 1)

**Lemma 6 (one sign change per section; arbitrary net flows).** Let s be an adiabatic CMO section of n ≥ 0 (E)
stages, numbered 1..n from the bottom, with L, V ≥ 0 and arbitrary net flows w. Fix θ > 0, θ ∉ {α_i}. Let γ_0 be the
sign of φ of the vapor entering stage 1 from below, γ_k = σ_{(stage k)} for 1 ≤ k ≤ n, and γ_{n+1} the sign of φ
of the liquid entering stage n from above (for n = 0: the vapor entering and the liquid entering the section).
Then:
- if Ψ_s(θ) > 0: as soon as some γ_k ≤ 0, all γ_{k'} with k' > k are −1. The sequence is (≥0, …, ≥0, −1, …, −1):
  at most one change, from heavy-ish below to light-ish above ("normal" orientation);
- if Ψ_s(θ) < 0: as soon as some γ_k ≥ 0, all γ_{k'} with k' > k are +1 ("inverted" orientation);
- if Ψ_s(θ) = 0: φ(liquid coming down) = φ(vapor going up) at every cross-section, and the γ_k are all equal
  where the flows are positive.

*Proof.* At the cross-section above stage k (0 ≤ k ≤ n): φ(liquid down) = φ(vapor up) − Ψ_s, the vapor coming from
stage k (or the entering vapor for k = 0) and the liquid from stage k + 1 (or the entering liquid for k = n). If
Ψ_s > 0 and the vapor from stage k has φ ≤ 0, the liquid from stage k + 1 has φ < 0, so (L > 0 and) G_{x_{k+1}} > 0,
so every outflow of stage k + 1 has φ ≤ 0 (Lemma 2), in particular its vapor; induct. If γ_k ≤ 0 for a stage k ≥ 1,
its vapor has φ ≤ 0; for k = 0 this is the hypothesis. The other cases are symmetric. ∎

This is the stage-level content of "the pinch zones of an infinite section": on the side of a section where the
stages are light-ish at θ the section can only continue light-ish upward when Ψ_s(θ) > 0, and a section with
Ψ_s(θ) < 0 transports heavy-ish material up. Agent N states the same fact as (N1) and checks (2.1) in exact
arithmetic (`total/network/scripts/section_identity_exact.py`).

*Pinch compositions.* If x* is a pinch of the section (x_{k+1} = x_k = x*), then x*_i = (S*/V) w_i/(α_i − φ_k) for
a root φ_k of Ψ_s (Th eq. 3.39; at a pinch w_i = x*_i (V α_i/S* − L), so either w_i = 0 = x*_i or
x*_i ∝ w_i/(α_i − φ_k) with φ_k = L S*/V). Then
G_{x*}(θ) = (S*/V) Σ_i α_i w_i/((α_i − φ_k)(α_i − θ)) = −(S*/V) Ψ_s(θ)/(θ − φ_k): the sign of a pinch zone at θ is
−sign Ψ_s(θ) · sign(θ − φ_k). (A pinch containing a component h of zero net flow, x*_h > 0 with w_h = 0, has
L S*/V = α_h, and then G_{x*}(θ) = −(S*/V) Ψ_s(θ)/(θ − α_h) + α_h x*_h/(α_h − θ), Ψ_s being computed without the
h-term, which vanishes.)

**(a) Absent components, (b) degenerate distributions, (c) reverse (negative) net flows, (d) intermediate
sections between two feeds or side streams, (e) finite stages.** Lemmas 1, 2, 6 and Theorem 3 contain no hypothesis
on (a)–(d): a component may be absent, present at one end only, flow against the section's direction, and a
section may join any two units (feed stages, side-draw stages, mixers, stages of other columns). (e) The model is the
finite-stage one; infinitely many stages and pure products are the limit of class.md §4, and every statement above is
an identity or an inequality at a finite-stage steady state, so it passes to any limit of steady states (the fluxes
are continuous in flows and compositions for θ ∉ {α_i}).

---

## 6. Junctions: the physical forms of result 12's column lemmas (brief task 2)

Result 12 used, for a two-section column at θ = θ_j with A = Ψ_T(θ_j), B = Ψ_B(θ_j), a = A − B (the column's feed
flux): S3/S4 ("A, B never both negative") and S5 (strict model only: "A, B never both positive when both products
are mixed"). The verifier showed that S5 fails at degenerate points of the Underwood model. The physical statements,
valid at every finite-stage operating point with any net flows and any distribution, are:

**Lemma 7 (two-section column).** Let a column consist of a top section, an (E) junction stage J (feed stage,
with any inflows, outflows and duty) and a bottom section, sections as in Lemma 6. Fix θ > 0. Let s_T be the sign of
φ of the liquid entering the top of the top section and s_B the sign of φ of the vapor entering the bottom of the
bottom section. Then:
- (S3/S4, physical) if s_T ≤ 0 ≤ s_B (a light-ish return at the top and a heavy-ish return at the bottom), then A
  and B are not both negative;
- (S5, physical) if s_T ≥ 0 ≥ s_B, then A and B are not both positive.

*Proof.* Suppose A < 0 and B < 0. Top section (bottom end J, top end the returned liquid), Ψ < 0: by Lemma 6 a
value γ_J ≥ 0 would force s_T = +1; so s_T ≤ 0 gives γ_J = −1. Bottom section (bottom end the returned vapor, top end
the liquid leaving J, of sign γ_J), Ψ < 0: s_B ≥ 0 forces γ_J = +1. Contradiction. The second statement is the same
with the orientations of Lemma 6 for Ψ > 0. ∎

At a column whose distillate is light (only components ≤ j) and residue heavy, the returns satisfy s_T ≤ 0 ≤ s_B
whenever the returned streams are themselves light and heavy, which is the situation of S3/S4. S5 is the case of a
column whose two returned streams have the "wrong" signs; in result 12's strict model it was imposed through the
labels (both products mixed), and at degenerate points the labels do not determine the signs of the returned
streams — which is why S5 fails there. The global argument (Theorem 3) needs neither: it uses the signs of the units
themselves.

**Lemma 8 (general junction).** Let J be an (E) unit at which any number of sections and streams meet. For every
section s with its bottom end at J: if Ψ_s(θ) > 0 and σ_J ≤ 0 then all its stages and the liquid it receives at its
top are light-ish (σ = −1); if Ψ_s(θ) < 0 and σ_J ≥ 0, all are heavy-ish (+1). For every section s with its top end
at J: if Ψ_s(θ) > 0 and σ_J ≥ 0 then all its stages and the vapor it receives at its bottom are heavy-ish; if
Ψ_s(θ) < 0 and σ_J ≤ 0, all are light-ish. (Lemma 6 applied with γ_0 = σ_J, resp. γ_{n+1} = σ_J.)

**Lemma 9 (where the section fluxes enter the cut).** In Theorem 3(c), a section whose sign sequence
(γ_0, …, γ_{n+1}) of Lemma 6 is not contained in one class (the classes being {≥ 0} for 𝒰₊ and {−1} for 𝒰₋; by
Lemma 6 the sequence changes class exactly once) contributes, at that change, the two crossing streams of one
cross-section, whose
contribution is |φ(vapor up)| + |φ(liquid down)| = |Ψ_s(θ_j)|, with Ψ_s(θ_j) ≥ 0 if the change is from 𝒰₊ below to
𝒰₋ above, ≤ 0 otherwise. So the defect of Theorem 3(c) contains Σ |Ψ_s(θ_j)| over the sections that change sign.
At an exact tie (R = Σ_P max(0, φ_P)) no section changes sign (if Ψ_s(θ_j) = 0 the signs are constant along s by
Lemma 6); along a family of steady states whose defect tends to 0, Ψ_s(θ_j) → 0 for every section that changes sign:
in the limit θ_j is a root of the section flux (an active, pinched Underwood root) of every section whose stages
change sign at θ_j. This is the stage-level meaning of HS's "active roots" and of result 12's Corollary tie.

---

## 7. Reverse flows and degenerate points (brief task 3; verifier's gaps (1), (2))

- **Reverse flows never lower R below the bound**: Theorem 3 holds with arbitrary signs of the net flows (no
  monotonicity of Ψ_s in θ is used; the only inequality is Lemma 2). Halvorsen's reverse-flow regimes (Th §4.3.3:
  B flowing down through the top of a directly coupled prefractionator, D = V_T − L_T < 0 allowed) are steady states
  of 𝒜 and satisfy R ≥ max_j H_j at pure products.
- **Degenerate points** (a nominally distributing component absent from one product; excess vapor in an
  infinite-stage column with fewer than two distributing components): physical points are steady states (or limits
  of steady states) of 𝒜, and Theorem 3 holds at all of them. No relabelling or perturbation is needed.
- **The permissive counterexample** (verification/, N = 5, R/V_FTC ≈ 0.99946) is not the limit of any family of
  physical steady states with purity → 1: by Corollary 5 any such family has lim inf R ≥ max_j H_j = V_FTC. This
  settles the verifier's question ("an argument, not a proof that no exotic adiabatic arrangement could emulate the
  point"): no arrangement of 𝒜, however exotic, emulates it.

---

## 8. Attainment

The upper half of the conjecture, R*(𝒮_FTC) ≤ max_j H_j (the generalized Petlyuk arrangement with enough stages
and purities → 1 needs no more than the highest peak), is HS's construction (every internal column at its
preferred split), proved in the infinite-stage Underwood model by result 12 (Theorem achieve). A proof in the
finite-stage model requires the convergence of finite-stage steady states to the Underwood point. Below: a complete
proof for N = 2 (Proposition 10), and the route and the missing step for N ≥ 3 (8.2). Agent P tests attainment
numerically for N = 3, 4, 5.

### 8.1 Two components: complete proof

Take α_1 = α > 1, α_2 = 1 (only ratios matter), feed light fraction z ∈ (0, 1), any q; θ = θ_1 ∈ (1, α),
H_1 = (1 − z)/(θ − 1). Write x, y for light mole fractions and e(x) = αx/(1 + (α − 1)x) (increasing, strictly
concave, e(0) = 0, e(1) = 1). The column: a partial reboiler (an (E) unit with duty R = V_B, liquid product B of
composition x_B), n_B stages, a feed stage, n_T stages, a total condenser (a (C) unit, reflux and distillate of
composition x_D). Flows: D = z, B = 1 − z, V_B = s, V_T = s + 1 − q, L_T = V_T − D, L_B = s + 1 − z.

**Proposition 10.** If s > H_1, then for every ε > 0 there is n_0 such that for all n_T, n_B ≥ n_0 the column has a
steady state with x_B ≤ ε and 1 − x_D ≤ ε, and R = s. Hence R*(𝒮_column) = H_1 and, with Corollary 5, for N = 2
the conjecture holds in 𝒜 with complete proofs: every arrangement of 𝒜 needs R*(𝒮) ≥ H_1, and one column attains it.

*Proof.* (i) Feasibility of the flows: H_1 + 1 − q = P^T_1 = αz/(α − θ) > z, so L_T > 0 for s ≥ H_1; and
H_1 > (1 − z)/(α − 1) since θ < α.
(ii) Pinch points. The "ideal" bottom map β_0(x) = (s/L_B) e(x) has the fixed points 0 and
x̂_B = (αs/L_B − 1)/(α − 1), which lies in (0, 1) because s(α − 1) > 1 − z by (i); β_0(x) > x on (0, x̂_B) (e(x)/x
is decreasing). The "ideal" top map going down, τ_0(x) = e^{-1}((L_T x + D)/V_T), has the fixed points 1 and
x̂_T = (V_T/L_T − 1)/(α − 1); between x̂_T and 1 the line (L_T x + D)/V_T lies below the concave curve e, so
τ_0(x) < x on (x̂_T, 1). Now x̂_T < x̂_B ⇔ V_T L_B < α s L_T ⇔ f(s) := α s (s + 1 − q − z) − (s + 1 − q)(s + 1 − z) > 0.
Eliminating q with the feed equation 1 − q = αz/(α − θ) + (1 − z)/(1 − θ), the quadratic f (leading coefficient
α − 1 > 0) has the roots H_1 and s_− = 1/(α − 1) − zθ/(α − θ) (sympy check, `scripts/check_binary_attainment.py`),
and H_1 − s_− > z α (θ − 1)/((α − θ)(α − 1)) > 0. Hence s > H_1 gives x̂_T < x̂_B (and x̂_T < 1, so x̂_T ∈ (0, 1)).
(iii) Shooting. For η ∈ [0, 1] put x_B = η, x_D(η) = 1 − (1 − z)η/z (light balance with D = z). Bottom profile
upward: x_0 = η (reboiler liquid), x_{k+1} = β_η(x_k) := (s e(x_k) + Bη)/L_B (the cross-section balance
s y_k − L_B x_{k+1} = −B x_B); after n_B + 1 steps this gives the feed-stage liquid x_f^bot(η). Top profile downward:
the top stage vapor equals x_D (total condenser: s + 1 − q = L_T + D), x = e^{-1}(y) on each stage and
y_{next} = (L_T x + D x_D)/V_T (the cross-section balance V_T y − L_T x = D x_D); after n_T + 1 stages this gives the
feed-stage liquid x_f^top(η). All iterates stay in [0, 1] (β_η([0,1]) ⊂ [0,1] since L_B = s + B; the top step is a
convex combination followed by e^{-1}). A pair of profiles with x_f^bot(η) = x_f^top(η) is a steady state: every
stage balance holds by construction, the feed-stage balance follows from the two cross-section balances and the
overall light balance, and the reboiler and condenser balances hold. Both x_f^bot and x_f^top are continuous in η.
At η = 0: x_f^bot = 0 and x_f^top = 1. Fix η_1 ∈ (0, 1] with η_1 max(1, (1 − z)/z) ≤ ε. Since β_{η_1} ≥ β_0 and
both are increasing, x_f^bot(η_1) ≥ β_0^{n_B+1}(η_1) → x̂_B (monotone convergence on (0, x̂_B), or x_f^bot ≥ x̂_B if
η_1 ≥ x̂_B); since the top step with x_D(η_1) ≤ 1 is ≤ τ_0 and τ_0 is increasing, x_f^top(η_1) ≤ max(x̂_T,
τ_0^{n_T}(e^{-1}(x_D(η_1)))) → x̂_T. So for n_T, n_B ≥ n_0(η_1): x_f^bot(η_1) > (x̂_T + x̂_B)/2 > x_f^top(η_1). By the
intermediate value theorem there is η* ∈ (0, η_1) with x_f^bot(η*) = x_f^top(η*): a steady state with
x_B = η* ≤ ε and 1 − x_D = (1 − z)η*/z ≤ ε. ∎

### 8.2 N ≥ 3: the route and the missing step

(N's projective form, network.md.) In a CMO section, x_{k+1} = M(w) x_k / (L S(x_k)) with M(w) = V diag(α) − w α^T;
M is analytic in w, its eigenvalues are V φ_r (the roots of Ψ_s) and its eigenvectors e_r(w) ∝ (w_i/(α_i − φ_r))_i
are the pinch compositions; when an impurity flow w_m → 0 the eigenvalue tends to V α_m and the eigenvector to the
pinch containing m at zero net flow, without singularity. For a section from X to Y with n stages,
x_Y ∝ M^{n+1} x_X. With the pinch index r* of the section and t_r = (φ_{r*}/φ_r)^{n+1} (upper modes) or
(φ_r/φ_{r*})^{n+1} (lower modes), the section relation is, in the chart where the e_{r*}-coefficient is 1, a smooth
function of (coefficients, w, t), and at t = 0 it is the **limit relation x_X ∈ E^−_{r*}(w), x_Y ∈ E^+_{r*}(w)**
(spans of the lower and upper eigenvectors). If w_m = 0 the upper span is the face without m: the limit condition at
a section's top end is Halvorsen's rule (Th §4.3.1: the returned stream contains no component that the section
removes). So the finite-stage FTC system is a smooth, exponentially small perturbation of the Underwood limit system,
and the implicit function theorem gives finite-stage steady states converging to every limit solution at which the
limit system is transversal; Halvorsen's insensitivity makes the limit system block-triangular over the FTC grid, so
transversality reduces to single columns. **Missing:** a proof that a single column's limit system is transversal at
points strictly inside a distribution region, for all data, and the double limit (margins → 0 at the HS knot,
n → ∞). Not done in this wave. *Correction (14:27, see the warning at the end of 8.4):* at sharp splits the limit
system written in the chart of Lemma 12 is singular for doubly pinched sections (a union of branches), so
"transversality" must be formulated in logarithmic coordinates of the exponentially small coefficients; the route
stands, but this step is more delicate than stated above.

### 8.3 The limit relation of a long section (for N's compactness route)

Agent N's route (network.md): (A) a finite-stage steady state exists for every consistent set of flows (degree
argument), so R is fixed by the flows and attainment reduces to purity; (B) along n → ∞ the junction compositions
and net flows have limit points, at which every section satisfies the limit relation below; (C) every limit state
with suitable flows is sharp. Here is (B), with the degenerate cases.

**Lemma 11 (limit relation).** Let a CMO section have fixed flows L, V > 0 and n_k → ∞ stages, and let the
steady states have net flows w^(k) → w, bottom-end "virtual liquid" x_X^(k) → x_X (the liquid in equilibrium with
the vapor entering the bottom: x ∝ y/α) and top-end entering liquid x_Y^(k) → x_Y. Let M(w) = V diag(α) − w α^T
and suppose that the limit matrix M(w) has N eigenvalues of pairwise distinct moduli, |λ_1| > … > |λ_N|, with left
eigenvectors ℓ_1, …, ℓ_N. Then there is r* ∈ {1..N} with

  ℓ_r · x_X = 0 for every r < r*,   ℓ_s · x_Y = 0 for every s > r*.

The left eigenvectors are explicit: for an eigenvalue λ = Vφ with φ a root of Ψ_s (Ψ_s(φ) = 0, φ ∉ {α_i}),
ℓ ∝ (α_i/(α_i − φ))_i, so ℓ · x = G_x(φ) (Underwood's function); for a component h with w_h = 0, λ = V α_h is an
eigenvalue with ℓ = e_h, so ℓ · x = x_h. (Every eigenvalue of M is of one of these two kinds: ℓ^T M = λ ℓ^T reads
ℓ_i (V α_i − λ) = (ℓ · w) α_i.)

*Proof.* At every cross-section x_above = M(w) x_below / (L S(x_below)) (component balance with
y_below = diag(α) x_below / S(x_below)); hence x_Y^(k) = c_k M(w^(k))^{n_k+1} x_X^(k) with c_k > 0. Simple
eigenvalues and their left eigenvectors depend continuously on the matrix, so for k large M(w^(k)) has eigenvalues
λ_r^(k) → λ_r and left eigenvectors ℓ_r^(k) → ℓ_r, and ℓ_r^(k) · x_Y^(k) = c_k (λ_r^(k))^{n_k+1} ℓ_r^(k) · x_X^(k).
For r < s, eliminating c_k: (ℓ_s^(k) · x_Y^(k))(ℓ_r^(k) · x_X^(k)) = (ℓ_r^(k) · x_Y^(k))(ℓ_s^(k) · x_X^(k))
(λ_s^(k)/λ_r^(k))^{n_k+1}. The right side tends to 0 (bounded factors, |λ_s/λ_r| < 1 in the limit), so
(ℓ_s · x_Y)(ℓ_r · x_X) = 0 for all r < s. Since x_X ≠ 0 and the ℓ_r form a basis, r* := min{r : ℓ_r · x_X ≠ 0}
exists, and then ℓ_s · x_Y = 0 for all s > r*. ∎

*Reading.* r* is the pinch of the long section. At the bottom end the Underwood functions vanish at every root above
the pinch root (the "active-root" conditions on a feed-stage composition); at the top end they vanish at every root
below it, and, for a component h with zero net flow and α_h below the pinch root, x_{Y,h} = 0: **the entering
liquid contains no component that the section removes** (Halvorsen's rule, Th §4.3.1, here a theorem about limits of
finite-stage steady states). Example: the top section of an FTC column [i, m] at the preferred split has w_m = 0,
eigenvalues V θ_i > … > V θ_{m−1} > V α_m and pinch r* at θ_{m−1}: the limit relation says G_{x_f}(θ_r) = 0 for
r = i..m−2 at the feed stage, and x_{Y,m} = 0 at the junction above.

*Hypothesis.* Distinct moduli hold whenever the roots of Ψ_s are real and simple and differ from the α_h of the
zero-flow components, e.g. for net flows of one sign in the section (Underwood's interlacing). With reverse flows
the roots can coincide or become complex (Th §4.3.4); then the lemma holds for the modes of distinct moduli and the
others need a Jordan/complex-pair version (not written).

*Caution for (C).* Petlyuk-type columns with specified flows may have several steady states; a statement "every
limit state with the flows F_δ is sharp" would exclude such multiplicity in the infinite-stage limit. If it fails,
the alternative is local: a nonzero local degree (e.g. transversality) of the limit system at the sharp HS solution,
which gives finite-stage steady states converging to it (8.2).

### 8.4 Uniform exponential closeness of a long section to its limit relation (tool for the local route)

**Lemma 12.** Fix L, V > 0 and a compact set 𝒲 ⊂ ℝ^N of net-flow vectors such that for every w ∈ 𝒲 the matrix
M(w) = V diag(α) − w α^T has N real nonzero eigenvalues of distinct moduli |λ_1(w)| > … > |λ_N(w)|. Let e_r(w),
ℓ_r(w) be right and left eigenvectors with ℓ_r · e_s = δ_rs, and ρ := max_{w∈𝒲} max_r |λ_{r+1}(w)/λ_r(w)| < 1.
Fix r* ∈ {1..N}. For a = (a_r)_{r>r*}, b = (b_r)_{r<r*} and t = (t_r)_{r≠r*} put

  u_X(a, b, w, t) = e_{r*} + Σ_{r>r*} a_r e_r + Σ_{r<r*} t_r b_r e_r,   u_Y(a, b, w, t) = e_{r*} + Σ_{r<r*} b_r e_r + Σ_{r>r*} t_r a_r e_r,

and τ_n(w)_r = (λ_{r*}/λ_r)^{n+1} for r < r*, (λ_r/λ_{r*})^{n+1} for r > r*. Then:
(i) (exact representation) M(w)^{n+1} u_X(a, b, w, τ_n(w)) = λ_{r*}(w)^{n+1} u_Y(a, b, w, τ_n(w)) for all a, b; and
every pair of compositions (x_X, x_Y) related by an n-stage section with net flows w (x_Y ∝ M(w)^{n+1} x_X, Lemma 11)
with ℓ_{r*} · x_X ≠ 0 is (x_X, x_Y) = (u_X/(1·u_X), u_Y/(1·u_Y)) for exactly one (a, b) (1 = the all-ones vector);
(ii) (uniform closeness) |τ_n(w)| ≤ ρ^{n+1} and |∂_w τ_n(w)| ≤ C (n + 1) ρ^n on 𝒲, with C = max over 𝒲 and r, s of
|∂_w log|λ_r/λ_s||. Hence, on every compact set of (a, b, w) on which 1·u_X and 1·u_Y do not vanish at t = 0, the
normalized maps (a, b, w) ↦ (u_X/(1·u_X), u_Y/(1·u_Y)) at t = τ_n(w) converge in C^1, at the rate O(n ρ^n), to the
same maps at t = 0, which parametrize the **limit relation** x_X ∈ span{e_r : r ≥ r*}, x_Y ∈ span{e_r : r ≤ r*}.

*Proof.* (i) M^{n+1} e_r = λ_r^{n+1} e_r gives M^{n+1} u_X = λ_{r*}^{n+1} e_{r*} + Σ_{r>r*} a_r λ_r^{n+1} e_r +
Σ_{r<r*} t_r b_r λ_r^{n+1} e_r = λ_{r*}^{n+1} u_Y when t = τ_n(w). Conversely, given x_X with ℓ_{r*} · x_X ≠ 0, put
a_r = (ℓ_r · x_X)/(ℓ_{r*} · x_X) for r > r* and b_r = (ℓ_r · x_X)/((ℓ_{r*} · x_X) τ_n(w)_r) for r < r*; then
x_X = (ℓ_{r*} · x_X) u_X, and x_Y ∝ M^{n+1} x_X ∝ u_Y. Uniqueness: the e_r form a basis. (ii) Simple eigenvalues
and their eigenvectors are real-analytic in w near 𝒲 (implicit function theorem on the characteristic polynomial),
so C < ∞; |τ_n| ≤ ρ^{n+1} by definition of ρ, and ∂_w (λ_r/λ_s)^{n+1} = (n + 1)(λ_r/λ_s)^{n+1} ∂_w log|λ_r/λ_s|. The
maps are polynomial in t and smooth in (a, b, w), so the C^1 distance at t = τ_n(w) to t = 0 is bounded by
(sup of the first and second derivatives in t on the compact set) × (|τ_n| + |∂_w τ_n|). ∎

*Use.* A finite-stage network steady state is a zero of junction balances in the junction compositions and net flows;
replacing every section relation by its parametrization (i) gives a system G_n(v) = 0 in the variables v =
(junction compositions, w, a, b), with G_n → G_∞ in C^1 on compacts (ii), G_∞ being Underwood's infinite-stage
system written as pinch-matching conditions. If G_∞ has a zero v_0 with invertible derivative, then G_n has zeros
v_n → v_0 (C^1-perturbation of the implicit function theorem). Two points remain to be checked in any application:
the transversality at v_0 (the crux), and **physical admissibility** of v_n — compositions that vanish at v_0 (a
component absent from a limit product or junction) must be shown to be ≥ 0 at v_n (the stage equations do not
preserve positivity by themselves, M having negative entries), e.g. by identifying the sign of their leading
exponentially small term.

*Warning: doubly pinched sections.* Lemma 12 needs bounded coefficients (a, b) in the chart, i.e. a single pinch
mode r* carrying an O(1) coefficient at both ends. At sharp splits near minimum vapor a section typically connects
two pinch zones, one at each end (in the binary column of Proposition 10 with n_T = n_B, when the bottom profile
leaves the pure-heavy pinch faster than the top profile leaves the pure-light one, the bottom section lingers near
x = 0 at the reboiler and near x̂_B at the feed stage, and the top section carries the transition). Then, for two
modes r_1 < r_2 (|λ_{r_1}| > |λ_{r_2}|), ℓ_{r_1} · x_X → 0 and ℓ_{r_2} · x_Y → 0, and the exact identity of the proof
of Lemma 11 reads (ℓ_{r_1} · x_X)(ℓ_{r_2} · x_Y) = (λ_{r_2}/λ_{r_1})^{n+1} (ℓ_{r_2} · x_X)(ℓ_{r_1} · x_Y): a product of
two small coefficients equals an exponentially small number times O(1). Its limit is the *union* of the two branches
ℓ_{r_1} · x_X = 0 and ℓ_{r_2} · x_Y = 0: the naive limit system is not a manifold there and is not transversal. A proof must use logarithmic coordinates for the exponentially small coefficients (as the binary proof
uses monotonicity instead). In naive variables the binary limit system is even discontinuous: the feed-stage
mismatch of Proposition 10 tends to −1 at η = 0 and to x̂_B − x̂_T > 0 for every fixed η > 0.

---

### 8.5 Critical check of GPT-6-Astra's prefractionator theorem (total/attainment-astra/attainment.md)

Requested by the orchestrator (14:52), since the N = 3 work builds on it. I re-derived every step of Sections 2–4
(Theorem 1) for m = 3 by hand, and the algebra of the redistribution for general m:
- (2): Ψ_T^0 has poles a_1..a_{m−1}, roots θ_1..θ_{m−1} and limit V_T^0, hence the product form; same for Ψ_B^0.
- (3)–(5): K(∞) = 1, residues a_i h_i at a_i, K(0) = β Π θ_j / Π a_i = 1 so Σ h_i = 0; for m = 3, β = a_2, K ≡ 1,
  h = 0 (the ternary construction keeps every net flow of the preferred split and adds η to the four phase flows).
  K(θ_1) > 0 (all factors positive) and K(θ_{m−1}) > 0 (m − 2 negative factors above and below the line). d, b stay
  positive for η < min d^0_i/|h_i| (the h_i may have either sign; only smallness of η is used).
- (8): at a simple root, Ψ^0 + ηK with K > 0 moves a decreasing (top) root up and an increasing (bottom) root down;
  middle roots stay common because K vanishes there.
- (11)–(13): e_u, e_v are the face pinches (normalized by Ψ_T(u) = 0, Ψ_B(v) = 0); (e_c)_m = Ψ_T(c)/L_T > 0 because
  c < u and Ψ_T decreases on (−∞, a_{m−1}); (e_a)_1 = Ψ_B(a)/L_B > 0 likewise (my own computation; Astra writes the
  complement form). The identity Ψ(t) = −(V/S(e))(t − φ) G_e(t) for a pinch e (with φ = L S(e)/V) gives G_e = 0 at
  every root of Ψ other than φ, including the zero-net-flow pinch.
- Lemma 2: G_{e_c}(s) = −S Ψ_T(s)/(V(s − c)) < 0 since c < s < u; G_{e_u}(s) > 0 since e_u lives on components lighter
  than s; symmetrically at r for the bottom. Hence strictly positive cone weights and X = A e_c + B e_u = C e_v + D e_a
  by the uniqueness in Lemma 1 (interlacing of (r, θ_2, …, θ_{m−2}, s) holds by (8)).
- (16)–(18): forward orbit of M_T, backward orbit of M_B, convex combinations with a full-support vector, so strictly
  positive; the balances (17) and the feed balance follow; x_{j,m} ≤ (A/B)(c/u)^j, x_{−j,1} ≤ (D/C)(v/a)^j; ρ can be
  taken η-independent (c/u < c/θ_{m−1}, v/a < θ_1/a); B, C ~ η (B'(0), C'(0) > 0), whence the factor 1/η.
- Theorem 2 (m = 3): the top eigenvector at r has components of opposite signs on A, B (d_A/(a − r) > 0 >
  d_B/(b − r)), so f_U e_r ≠ 0; the scaled rows converge to G(r), G(s); M_T^{−1} ≥ 0 (Sherman–Morrison, L_T > 0) and
  M_B ≥ 0 propagate positivity from the two endpoints. Correct.
Independent numerical re-implementation from the formulas only (`scripts/check_astra_prefrac.py`, 60 digits; m = 3,
4, 5; q = −0.5, 0, 0.5, 1, 1.4; η = 1e-2, 1e-4, 1e-6): every claim (a)–(f) of the script holds, balance residuals
≤ 7e-60, B/η and C/η converge to positive limits. PASS (`runs/check_astra_prefrac.log`).
**Verdict: Theorem 1 and Theorem 2 are correct as stated.** Their scope is exactly what Astra says: the return
compositions are outputs of the construction (Theorem 2 prescribes two major-component ratios), so closing the FTC
network remains open; the trace bounds carry the factor η^{−1}, so the double limit needs n ≫ log(1/η)/log(1/ρ).

### 8.6 N = 3: the main column of the Petlyuk arrangement (building blocks for the closure; joint work with N)

**Status (15:07, final for this wave).** PROVED: the count and the exact net-flow bookkeeping; Proposition 10'
(binary two-section column with general ends: existence, uniqueness and continuity of the matching parameter); the
equivalence of its pinch conditions with the HS vapor margins computed with the prefractionator's actual roots
(sympy, `scripts/check_finishing_pinch.py`); Lemma 14 (exact trace relations, positivity rules) and the explicit
solution of each trace block under the stated affine form of the prefractionator's trace relation. SKETCHED / NOT
PROVED: that the exact ternary main column (traces inside S(x) of the major profiles) inherits the binary
matching — the block-triangular "leading-order" structure below is a description, not a proof — and the closure
with the prefractionator (X depending on the junction ratios, Astra's construction with ε ≠ 0), i.e. attainment for
N = 3 remains open here. A route for the closure: block Gauss–Seidel map on (σ_A, σ_C, trace profiles, ε, X) with
Brouwer's theorem — the A-recursion of the AB column with a GIVEN C-profile is still a scalar map increasing in x_A
and monotone in σ_A, so the AB block has a unique continuous root (same for BC), and the trace blocks are uniquely
solvable (Lemma 14); what is missing is an invariant box (a priori bounds) and the handling of the signs of the trace
net flows w_{U,C}, w_{W,A} at the ends of the shooting ranges (positivity).

Notation (ternary A, B, C, volatilities a > b > c): prefractionator P (top net flows d_P, bottom b_P = z − d_P, top
vapor V_P, bottom vapor V_PB, liquids L_P = V_P − Σd_P, L_PB = V_PB + Σb_P); main column, top to bottom: total
condenser, top section T (V_T, L_T), junction J1 (an (E) stage receiving P's top vapor and returning liquid to P),
upper-middle section U (V_U = V_T − V_P, L_U = L_T − L_P), side-draw stage S (liquid product S), lower-middle section W
(V_W = V_U, L_W = L_U − S), junction J2 (receives P's bottom liquid, sends vapor to P), bottom section Bo
(V_Bo = V_W + V_PB, L_Bo = L_W + L_PB), partial reboiler (product B). Product flows D = z_A, S = z_B, B = z_C.

**Count.** Given the flows, a steady state is determined by: the product compositions x_D, x_S, x_B (6 numbers, minus
the 2 independent overall balances D x_D + S x_S + B x_B = z), P's feed-stage composition X (2) and P's two trace net
flows ε_C = d_{P,C} (C up P's top) and ε_A = b_{P,A} (A down P's bottom): 8 unknowns. All net flows follow:
w_T = D x_D, w_U = w_T − d_P, w_W = w_U + S x_S, w_Bo = w_W − b_P = −B x_B. Equations: at J1 the composition reached
by T (from x_D downward), by U (from x_S upward) and by P's top section (P's upper return) coincide (4); the same at
J2 for W (from x_S downward), Bo (from x_B upward) and P's lower return (4). Freezing ε_C = ε_A = 0 (exactly sharp
prefractionator net flows, as in Astra's Theorems 1–2) leaves 6 unknowns for 8 equations: the two "doubled rows".
Write x_D = (1 − ε_D − κ_D, ε_D, κ_D), x_S = (σ_A, 1 − σ_A − σ_C, σ_C), x_B = (κ_B, ε_B, 1 − κ_B − ε_B); the balances
give z_A(ε_D + κ_D) = z_B σ_A + z_C κ_B and z_C(κ_B + ε_B) = z_A κ_D + z_B σ_C, so (σ_A, σ_C, κ_D, κ_B) are free.

**Leading-order structure (block triangular).** (i) The AB finishing column (T + J1 + U, components A, B; C is a
trace) is the binary problem of Proposition 10' below with shooting parameter σ_A; (ii) the BC finishing column
(W + J2 + Bo; A a trace) the same with σ_C; (iii) P with prescribed A:B ratio at J1 and B:C ratio at J2 is Astra's
Theorem 2 (unknown X); (iv) the C trace of the AB column (unknowns κ_D, x_{J1,C}, ε_C) and the A trace of the BC
column (κ_B, x_{J2,A}, ε_A) satisfy the exact affine relations of Lemma 14 below. The perturbative coupling of these
blocks (traces in S(x) of the major profiles, X depending on the junction ratios) is the closure (Agent N).

**Proposition 10' (binary two-section column with general ends).** Binary, light fraction x,
e(x) = αx/(1 + (α − 1)x), α > 1. Upper section Σ_1 (n_1 stages, V_1 > L_1 > 0) and lower section Σ_2 (n_2 stages,
L_2 > V_2 > 0) meet at an (E) junction J. A parameter σ ∈ [0, σ̄] fixes the ends: the vapor leaving the top stage of
Σ_1 has light fraction y_1(σ) = 1 − c_1 σ, the net light flow up Σ_1 is w_1(σ) = (V_1 − L_1) − c_2 σ, the bottom
end of Σ_2 is an (E) stage with liquid light fraction c_3 σ, and the net light flow up Σ_2 is −c_4 σ, with constants
c_i > 0 and c_2 = c_4 (so that the net light feed received at J, w_1 − w_2 = V_1 − L_1, does not depend on σ). Profiles: downward in Σ_1, x = e^{-1}(y) on each stage and y_next = (L_1 x + w_1(σ))/V_1; upward in Σ_2,
x_{k+1} = (V_2 e(x_k) + c_4 σ)/L_2. Let x̂_1 = (V_1/L_1 − 1)/(α − 1) and x̂_2 = (α V_2/L_2 − 1)/(α − 1). If
0 < x̂_1 < x̂_2 < 1, then for every ε > 0 there is n_0 such that for all n_1, n_2 ≥ n_0 some σ ∈ (0, ε] makes the two
profiles reach the same liquid at J (a binary steady state with both end impurities O(σ)).
*Proof.* As in Proposition 10 (iii): both maps are increasing in x; the Σ_1 step is ≤ the σ = 0 step, the Σ_2 step
≥ the σ = 0 step; at σ = 0 the profiles stay at 1 and 0 (mismatch −1); for fixed σ_1 > 0 the σ = 0 maps converge
monotonically to x̂_1 (from above) and x̂_2 (from below), so for n_1, n_2 large the mismatch at σ_1 is
≥ (x̂_2 − x̂_1)/2 > 0; continuity in σ and the intermediate value theorem. (σ̄ small keeps every iterate in [0, 1].) ∎
*Uniqueness and continuity (also for Proposition 10).* The Σ_2 value at J is strictly increasing in σ (strictly
increasing start c_3 σ, maps increasing in σ and strictly increasing in x) and the Σ_1 value strictly decreasing, so
the mismatch is strictly increasing: the matching σ* is unique and depends continuously on the flows and on any
parameter on which the maps depend continuously (the root of a continuous strictly monotone family). In particular
the binary column of Proposition 10 with the flows given there has exactly one steady state (every steady state is
the shooting profile of its own x_B), consistent with the classical uniqueness for molar flow specifications.
For the Petlyuk main column: AB column, Σ_1 = T, Σ_2 = U, c_1 = S/D, c_2 = S, c_3 = 1, c_4 = S (S's A content
σ_A is the impurity carried down U); BC column with B as the light component, Σ_1 = W (top end S:
y_1 = 1 − σ_C L_U/V_U, w_{W,B} = (V_W − L_W) − S σ_C), Σ_2 = Bo (reboiler, x_{B,B} = S σ_C / B, c_4 = S). The pinch
conditions are exactly the vapor margins over the Underwood minima computed with P's actual section roots
(`scripts/check_finishing_pinch.py`, sympy, PASS): with P's top root r ∈ (b, a) and bottom root s ∈ (c, b),
x̂_T < x̂_U ⇔ V_T > a z_A/(a − r), and x̂_W < x̂_Bo ⇔ V_Bo > c z_C/(s − c) (in each case the other root of the
quadratic is smaller). With P at the preferred split plus slack η (Astra), r ↓ θ_1 and s ↑ θ_2 as η → 0, so the two
conditions become V_T > P^T_1 and V_Bo > H_2, i.e. R = V_Bo > max(H_1, H_2) in the limit: the margins of HS.

**Lemma 14 (exact trace relations of the main column).** For any steady state and any component h, in every
section the h-balance at consecutive cross-sections is the scalar affine recursion
x_{h,k+1} = μ_k x_{h,k} − w_h/L with μ_k = V α_h/(L S(x_k)) > 0 (going up), equivalently
x_{h,k} = (x_{h,k+1} + w_h/L)/μ_k (going down). Hence, with Π = Π μ_k > 0 and Σ = Σ_k Π_{j>k} μ_j / L > 0 over the
section (all computed on the actual profile):
- (T, C trace) going down from the condenser, x_{J1,C} = Γ_T κ_D with Γ_T > 0 (explicit in the S(x_k) of T): the
  C-content of J1 is a positive multiple of the C-content of the distillate, and every C-content in T is positive;
- (U, C trace) x_{J1,C} = Π_U σ_C − w_{U,C} Σ_U with w_{U,C} = D κ_D − ε_C;
- eliminating: κ_D (Γ_T + D Σ_U) = Π_U σ_C + ε_C Σ_U, so κ_D ≥ 0 iff Π_U σ_C + ε_C Σ_U ≥ 0, and then every
  C-content in U is positive (going down from a positive end value with w_{U,C} > 0, or going up with w_{U,C} ≤ 0);
- symmetric statements for the A trace of the BC column: going up from the reboiler, x_{J2,A} = Γ_Bo κ_B with Γ_Bo > 0;
  in W, x_{J2,A} relates to σ_A, w_{W,A} = ε_A − B κ_B by the same affine rule (from S downward:
  x_{J2,A} = Π'_W σ_A + w_{W,A} Σ'_W with the downward multipliers 1/μ_k).
Proof: the h-component of the cross-section balance w = V y(x_k) − L x_{k+1} with y_h = α_h x_h/S(x_k); the condenser
gives y_C(top stage) = κ_D, so x_C(top stage) = S κ_D/c > 0 and each downward step adds w_{T,C}/L ≥ 0 before
dividing by μ_k > 0; the reboiler is the mirror image. ∎
These relations, together with P's trace relation, close the two doubled rows. If, with P's major profile frozen,
the C-content of P's upper return is affine in ε_C, x^P_{J1,C} = P_0 − ε_C Σ_P with P_0 ≥ 0 and Σ_P > 0 (as P's own
C recursion gives near its top end), then the C block {x_{J1,C} = Γ_T κ_D, x_{J1,C} = Π_U σ_C − (D κ_D − ε_C) Σ_U,
x_{J1,C} = P_0 − ε_C Σ_P} has the unique solution
x_{J1,C} = (Π_U σ_C + P_0 Σ_U/Σ_P)/(1 + D Σ_U/Γ_T + Σ_U/Σ_P) > 0, κ_D = x_{J1,C}/Γ_T > 0,
ε_C = (P_0 − x_{J1,C})/Σ_P; the A block is the mirror image. The sign of ε_C (resp. ε_A) is not fixed a priori; a
reverse trace flow (ε_C < 0: C net downward in P's top section) is allowed in 𝒜 and keeps P's C-contents positive
going up. Turning this leading-order block structure into a proof (the traces enter S(x) of the major profiles; X
depends on the junction ratios; Astra's P with ε ≠ 0) is the closure step (Agent N's 4-dimensional reduction).

## 9. Numerical checks (evidence only; every statement above is proved)

`scripts/stagenet.py` (own MESH/CMO solver: Gauss–Seidel with exact stage flash, damped Newton with analytic
Jacobian) and `scripts/check_stage_cut.py` (Python 3, numpy, scipy; WSL Fedora sci environment; a few minutes):
random columns (N = 2..5, 0–14 stages per section, random D, V_B, q), three-product Petlyuk arrangements (0–9
stages per section, random internal flows), random networks of equilibrium stages and mixers with recycles, split
feeds, vapor products and duties anywhere. At every converged steady state and every θ_j: global identity
Σ_P φ_P(θ_j) = R − Q (residual ≤ 1e-11), sign uniformity of every (E) unit, cut identity (a) for 𝒰₊ (residual ≤
1e-12), bound (b) never violated (`runs/check_stage_cut_25_2.log`: 24 columns, 18 Petlyuk arrangements, 75 random
networks). Command (from the project root, Git Bash):
`outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/total/sections/scripts && python check_stage_cut.py 25 2"`.

`scripts/check_binary_attainment.py` (Python 3, sympy, mpmath, numpy, scipy; seconds): the two exact identities used
in Proposition 10 (sympy: f(H_1) = 0, f(s_−) = 0, the formula for H_1 − s_−), and an illustration of the shooting
construction (120-digit arithmetic, four data sets with q = 1, 0, 0.5, 1.3, s = 1.02 H_1, n = 20, 40, 80, 160 stages
per section; every profile re-checked as a steady state by `stagenet.py`, residual ≤ 3e-15; impurities decrease
geometrically in n, e.g. 5e-3, 2e-5, 1e-10, 6e-21 for α = 2, z = 1/2, q = 1). PASS (`runs/check_binary_attainment.log`).

`scripts/petlyuk_fixed_flows.py` (exploration only, not used in any proof): ternary Petlyuk with all flows fixed at
the HS point plus a 2 % margin, α = (4, 2, 1), z = 1/3, q = 1. Converged steady states (residual < 1e-14) for
n = 4, 8, 16 stages per section: product impurities 1.1e-1, 4.8e-2, 1.3e-2 (largest), decreasing; for n ≥ 32 the
solver of `stagenet.py` did not converge (rows with residual > 1e-10 in `runs/petlyuk_fixed_flows_421.log` are not
steady states). Systematic finite-stage attainment runs are Agent P's (`total/simulation/`).
