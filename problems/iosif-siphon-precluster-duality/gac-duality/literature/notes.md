# Literature map: dual reaction networks, reverse-duals, and the Global Attractor Conjecture

Agent L — Claude Opus 5.5 (effort xhigh). Started 27 October 2026. Working notes for agents P (petri/), G (gac/),
E (computation/) and for the paper. Written incrementally: the most urgent entries (Blanchini–Giordano,
Teruel–Silva, then Craciun v3, Hack/Desel–Esparza/Jančar et al.) come first.

Conventions. N has species X₁..Xₙ, reactions R₁..R_r, source matrix Ye, product matrix Yp, Γ = Yp − Ye (n × r),
s = rank Γ. Petri dictionary (exploration brief, Fact 4): places = species, transitions = reactions,
Pre = Ye, Post = Yp, incidence C = Γ. Iosif's dual N* = weighted reverse-dual N^rd: places T, transitions P,
Pre^rd = Preᵀ, Post^rd = Postᵀ, incidence Cᵀ = Γᵀ.

Status labels used below: **[read]** = statement checked in the full text (page given); **[abstract]** = only the
abstract or a summary was available; **[derived]** = an observation of mine, proved in two lines here but not
checked by computation — teammates should re-derive before relying on it.

### Key points (read these first)

1. **GAC status.** The GAC is claimed in full by Craciun (arXiv:1501.02860 v3, 23 Sep 2026; induction on the number of
   species) and follows from OpenAI's boundedness-and-persistence theorem for all WR mass-action systems (25 Sep 2026,
   github.com/openai/math, with a Lean statement and proof catalogued). Treat it as settled; aim at transfers of
   structure, not at "new GAC cases" (§3, §13).
2. **Blanchini–Giordano** transfer only *structural local stability* (Thm 1, via SG vs GS spectra) and *polyhedral
   Lyapunov norms* (Thm 3, x-system of N ↔ w-system of N*), for general monotone non-autocatalytic kinetics; no
   persistence, siphons, complex balance or GAC; they already state that the dual of a reversible network is not even
   WR and that mass action is not preserved (§1).
3. **Petri side.** Hack 1972 (reverse-dual, duality theorem for FC), Desel–Esparza Thm 6.14 (rank theorem: FC net
   well-formed ⇔ connected + positive S- and T-invariants + rank = #clusters − 1), Thm 6.17 (live and bounded ⇔
   well-formed + every proper siphon marked), Thm 6.21 (duality); Teruel–Silva Thm 28 (same rank theorem for weighted
   **equal-conflict** nets, *pure* nets = no catalysts), Cor. 32 (duality only for **balanced-joins** EC nets); no
   weighted siphon theory exists (§2, §5–7). Equal-conflict sets = reactions with the same source complex = seed of
   Iosif's precluster graph (§2, §17).
4. **Siphons ↔ T-traps**: definitional; named "minimal T-traps" by Jančar et al. 2026 §4 for ordinary FC nets (§6).
5. **Known small facts to cite**: WR + free choice ⇒ disjoint complexes and δ = 0 (Mairesse–Nguyen 2009, §11);
   toric steady states from a 1-dimensional kernel (Pérez Millán et al. 2012, §16); conservation pools as extreme
   pathways of −Sᵀ (Famili–Palsson 2003, §10).
6. **"Small codimension"**: three readings in §12 (n − s; CDSS's codimension δ of the toric locus; codimension of
   boundary faces, where Anderson–Shiu's facet theorem is the classical result).
7. **Novelty** (§18): nobody has connected the CRN dual with Petri reverse-duality, siphons, persistence or the GAC;
   no work cites or answers Iosif's Conjecture 22 / Question 24.

---

## 1. F. Blanchini, G. Giordano, *Dual chemical reaction networks and implications for Lyapunov-based structural stability*, IEEE Control Systems Letters 6 (2022) 488–493, DOI 10.1109/LCSYS.2021.3081369

**What was read.** Full text (accepted version as distributed by IEEE Xplore, 6 pages), on disk in
`certificates/sources/Blanchini-Giordano-2022-LCSS.{pdf,txt}`; formulas checked on the rendered PDF pages 2–5.

**Setting (p. 1–2) [read].** General (not mass-action) networks ẏ = S g̃(y) + g₀, S ∈ ℤ^{n×m} (n species,
m reactions). They assume an equilibrium ȳ exists and shift it to 0: ẋ = S g(x) (2). Standing assumptions:
- A1: g(0) = 0;
- A2: g is C¹;
- A3: g_j depends on x_i **iff S_ij < 0**, and g_j = 0 when x_i = 0 (this "rules out autocatalytic reactions, such
  as A → 2A or A + B → 2B");
- A4: if g_j depends on x_i then ∂g_j/∂x_i > 0.

**Definition 1 (dual CRN) [read].** The dual CRN* of (2) is ż = Sᵀ h(z) (3), "where function h satisfies the same
assumptions as g, but for the new stoichiometric matrix S* = Sᵀ". "If a CRN has n species and m reactions, then its
dual CRN has m species and n reactions." So the dual is a *class* of systems determined by the **net** stoichiometric
matrix only; reactant sets are read off from the negative entries of Sᵀ.

*Relation to Iosif's dual [derived].* Iosif's N* has source matrix Yeᵀ and product matrix Ypᵀ (exploration brief,
Fact 4); Blanchini–Giordano's dual is determined by Γᵀ alone. (i) They always have the same stoichiometric matrix Γᵀ.
(ii) They have the same *reactant supports* (the kinetic support of A3) iff Ye_ij > 0 ⇔ Γ_ij < 0 for all i, j, i.e.
no species is a non-consumed reactant (catalyst, Γ_ij = 0, or autocatalyst, Γ_ij > 0) of any reaction; partial
consumption such as 2A → A + B is allowed, and this condition is transpose-invariant (agent G makes the same
point). (iii) They coincide as networks iff N is pure (Ye_ij · Yp_ij = 0; then Ye = max(−Γ, 0), Yp = max(Γ, 0)): for
2A → A + B Iosif's dual has the arc pattern 2k → k + …, the net-stoichiometry reading k → …. A3 excludes mass-action
catalysis and autocatalysis in any case. Example 1 of the paper (A+B →p C, C+D →q E+A, E →r B+D; dual P →a Q, P →b R,
Q →c P, Q →d R, R →e Q) and Example 5 (2A+B ⇄ C ⇄ D; dual 2P → 2Q, P → Q, Q+R → P+S, S → R) agree with Iosif's
definition (checked by hand: e.g. species C of Example 5 is a reactant of q and r and a product of p and s, so the
dual reaction c is k_q + k_r → k_p + k_s).

**Qualitative classes (eq. (8)–(10), p. 3) [read].** Σ = {M ∈ ℝ^{n×m}: sign(M) = sign(max{−S, 0})}, Σᵀ = {M: Mᵀ ∈ Σ}.
Jacobian of the primal: S ∂g/∂x = SG with G ∈ Σᵀ; of the dual: Sᵀ ∂h/∂z = SᵀH with H ∈ Σ.

**Definition 2 [read].** (2) is SSC (stable in the stoichiometric compatibility class) if every trajectory with
x(0) ∈ im S converges to 0; LSSC if the linearisation at 0 is stable in im S; S-SSC / S-LSSC = "structurally": the
property holds for **all** g satisfying A1–A4.

**Theorem 1 (p. 3) [read].** "System (2) is S-LSSC if and only if its dual system (3) is S-LSSC."
Proof idea: S-LSSC ⇔ SG Hurwitz on im S for all G ∈ Σᵀ ⇔ the r = rank S nonzero eigenvalues of SG have negative
real part (n − r zero eigenvalues); choose H = Gᵀ ∈ Σ; SG and GS have the same nonzero eigenvalues, and
SᵀH = SᵀGᵀ = (GS)ᵀ. So the statement is really pointwise: **for each admissible Jacobian factor G of N, the
dual with H = Gᵀ has the same nonzero spectrum.** Nothing beyond linear algebra of SG vs GS is used.

**Theorem 2 (p. 4) [read].** Jacobian decompositions SG = BΔC (concentrations) and GS = EΔF (rates), Δ diagonal
positive of size q = number of nonzero partial derivatives. For the dual: B*ΔC* = FᵀΔEᵀ and E*ΔF* = CᵀΔBᵀ
(Table I): "the dual CRN* in concentrations has the same structure as the transpose of the primal CRN in rates", and
vice versa.

**Theorem 3 (p. 5) [read].** Four representations (Table II): x-system ẋ = BD(t)Cx (N, concentrations),
r-system ṙ = ED(t)Fr (N, rates), z-system ż = FᵀD(t)Eᵀz (N*, concentrations), w-system ẇ = CᵀD(t)Bᵀw (N*, rates),
D(t) diagonal positive (differential inclusions absorbing the nonlinear system, eq. (14)–(17)). "The CRN in
concentrations (x-system) admits a PLF [polyhedral Lyapunov function] if and only if the CRN* in rates (w-system)
admits the dual PLF; the same holds for the z-system and the r-system." Dual norm V*(z) = max_{V(x)≤1} xᵀz; vertex
representation V_X(x) = inf{‖p‖₁ : Xp = x} is dual to plane representation U_Z(z) = ‖Zz‖_∞ with Z = Xᵀ. Proof:
matrix measure μ_V(A) = μ_{V*}(Aᵀ).
Remark 2: a Lyapunov norm alone does not give asymptotic stability; one also needs 0 < ε ≤ D_i(t) ≤ μ (18) and a
non-singular inclusion (citing their Automatica 2017 paper). Remark 4: the result extends to any norm (Goebel et
al. 2006, conjugate convex Lyapunov functions for dual linear differential inclusions).

**Propositions 1–2 (p. 5) [read].** Unitary networks (S_ij ∈ {−1, 0, 1}). Prop. 1 (from their Automatica 2014
paper): at most two nonzero entries per **column** of S ⇒ structurally stable with ‖x‖₁ as PLF. Prop. 2 (new, by
duality): at most two nonzero entries per **row** of S ⇒ structurally stable with ‖w‖_∞ as PLF in rates.
Examples 2, 4, 6 use the dual to certify stability of the primal (Example 6: 2A+B ⇄ C ⇄ D, non-unitary, whose dual
becomes unitary after lumping 2a(p)+b(p) = g(p)).

**What they say about mass action, weak reversibility, GAC (p. 6) [read].** Verbatim, short: "the dual of a
reversible network (such that all reactions are reversible) is not reversible, and not even weakly reversible", and
for mass-action rates "these functional expressions are not preserved for the dual" (Example 5: primal rates
p = a²b, q = c, r = c, s = d; dual rates a = p², b = p, c = qr, d = s, i.e. the dual is again mass action in its own
variables but is a different system, not the primal's rates). Closing sentence: future work on "how the relation
between primal and dual chemical reaction networks reflects on the structure of other types of graph representations,
such as the Species-Reaction Graph", and on "exploiting dual networks to assess properties other than stability".
**No mention** (grep of the full text and reading) of complex balance, detailed balance, persistence, siphons,
toric systems, the Global Attractor Conjecture, or Petri nets. Deficiency appears only in the introduction
(zero- and one-deficiency theorems as background). Earlier credit: Famili and Palsson (Biophys. J. 85 (2003)
16–26) proposed the "compound map" associated with −Sᵀ; Blanchini–Giordano say their stability results "cannot be
obtained for the compound map".

**Consequences for this project.**
1. The equivalence is for *structural* (all monotone kinetics) *local* stability of a *given* equilibrium, and for
   polyhedral Lyapunov functions of differential inclusions. It says nothing about existence or uniqueness of
   equilibria in the dual, boundary behaviour, or persistence. The dual system is *assumed* to satisfy A1–A4
   (in particular to have an equilibrium at 0 after shifting).
2. Mass-action complex-balanced systems: the Lyapunov function there is the pseudo-Helmholtz free energy, not
   polyhedral, and structural stability over all monotone kinetics is a much stronger property than what mass
   action needs. Nothing in the paper transfers persistence.
3. [derived, for agent G; to be checked] Under mass action at a positive point x with flux vector v = κ ∘ x^{Ye},
   G = diag(v) Yeᵀ diag(x)⁻¹. Then Gᵀ = diag(x)⁻¹ Ye diag(v), which has the mass-action form diag(w) (Yeᵀ)ᵀ diag(z)⁻¹
   of a Jacobian factor of N* (source matrix Yeᵀ) with w = 1/x (fluxes of N*) and z = 1/v (concentrations of N*).
   So the *set* of mass-action Jacobian factors at arbitrary positive points transposes onto the corresponding set
   for N* (pure networks), via (x, v) ↦ (z, w) = (1/v, 1/x). But an equilibrium of N (Γv = 0) is mapped to a point
   with Γᵀw = Γᵀ(1/x) — zero only if 1/x is a conservation law of N. So Theorem 1's pointwise mechanism does not map
   equilibria of N to equilibria of N* under mass action.
4. Fact 2 of the exploration brief (weak reversibility is not preserved) is **already stated** by
   Blanchini–Giordano (p. 6, Example 5, reversible ⇒ dual not even WR). Cite them for it.

---

## 2. E. Teruel, M. Silva, *Structure theory of equal conflict systems*, Theor. Comput. Sci. 153 (1996) 271–300, DOI 10.1016/0304-3975(95)00124-7

**What was read.** Full text (`certificates/sources/Teruel-Silva-1996-TCS.{pdf,txt}`), all sections; Definitions 5–6
and Theorems 28, 30, Corollary 32 checked on the rendered pages 278, 292, 293.

**Conventions (p. 274) [read].** P/T net N = (P, T, W), Pre, Post, C = Post − Pre. "Without loss of generality, nets
are assumed to be pure" (Pre[p,t]·Post[p,t] = 0): **in CRN terms, no species on both sides of a reaction.** Reverse
N^r (incidence −C), dual N^d (−Cᵀ), reverse-dual N^rd (incidence Cᵀ). Consistent: ∃X ≥ 𝟙 with CX = 0;
conservative: ∃Y ≥ 𝟙 with YC = 0. Well-structured = str. bounded and str. repetitive ⇔ consistent and conservative
(Thm 3.3, Memmi–Roucairol, Shields). Well-formed = str. bounded and str. live (some marking makes it live).
Minimal semiflow = support-minimal, gcd 1.

**Relations on transitions (Def. 5, p. 277–278) [read].**
- choice relation: t = t′ or •t ∩ •t′ ≠ ∅ (not transitive) — *two reactions share a reactant species*;
- coupled conflict relation: transitive closure; classes = coupled conflict sets, set 𝒞;
- equal conflict relation: t = t′ or Pre[P,t] = Pre[P,t′] ≠ 0; classes = equal conflict sets, set ℰ — *reactions
  with the same nonzero source complex*. **This is exactly the initial edge set of Iosif's precluster graph
  (Definition 19 of the manuscript: reactions with equal source)**, except that reactions with source complex 0 are
  singletons here.
- reverse-dual on places: join relation (common output transition = *two species are reactants of a common
  reaction*), coupled synchronisation sets 𝒮. The input places of a coupled conflict set form a coupled
  synchronisation set; in nets without sink places or source transitions (e.g. strongly connected), |𝒞| = |𝒮|.
  Footnote 3: a coupled conflict set plus its coupled synchronisation set is what Desel calls a *cluster*.

**Subclasses (Def. 6, p. 278) [read].** CF (choice-free): |p•| ≤ 1 (every species is a reactant of at most one
reaction); JF (join-free): |•t| ≤ 1 (every source complex is 0 or a multiple k·X_i); **EC (equal conflict):
•t ∩ •t′ ≠ ∅ ⇒ Pre[P,t] = Pre[P,t′]** (two reactions sharing a reactant species have the same source complex);
TEFC: •t ∩ •t′ ≠ ∅ ⇒ •t = •t′ (same reactant support). CF, JF, EC generalise marked graphs, state machines and
extended free choice to weighted nets. In EC nets ℰ = 𝒞 ("all conflicts are equal").
[derived] EC for a CRN ⇔ distinct source complexes have disjoint supports.

**Components and allocations (§4) [read].** P-component (Def. 16): strongly connected, join-free, conservative
P-subnet (subnet generated by a set of places with all their input and output transitions); T-component by
reverse-duality (N′ is a T-component of N iff (N′)^rd is a P-component of N^rd). Thm 17: the places of a
P-component are the support of a minimal P-semiflow, and P-components are maximal strongly connected JF subnets
(converses false in general, Fig. 1). T-allocation (Def. 18): for each coupled conflict set c, a maximal subset of c
with no two transitions in choice relation; T-allocatable (Def. 19): for every T-allocation the T-subnet generated
by the allocated transitions has a T-semiflow. Thm 20: N conservative and **rank(C) < |𝒞| ⇒ T-allocatable**.
Thm 21: in a conservative net the support of such a T-semiflow generates a T-component. Thm 24, Cor. 25:
T-allocatable strongly connected TEFC ⇒ consistent, and T-decomposable if conservative. Thm 26: T-allocatable
conservative strongly connected TEFC ⇒ every minimal T-semiflow generates a T-component, and maximal strongly
connected CF subnets are T-components.

**Theorem 28 — rank theorem for EC nets (p. 292) [read].** Let N be an EC net. Equivalent:
1. N is well-formed; 2. N is well-structured (consistent and conservative) and **rank(C) = |𝒞| − 1**;
3. N is strongly connected and P-allocatable; 4. N is strongly connected and T-allocatable.
(Necessity via Thm A.2, valid for *all* P/T nets: well-formed ⇒ well-structured and rank(C) < |ℰ|.)

**Theorem 29–30 (p. 292–293) [read].** A well-behaved (live and bounded) EC system has a well-formed net.
Thm 30: an EC system (N, M₀) is live and bounded iff N is well-structured, rank(C) = |𝒞| − 1, and for every
P-component N′ the system (N′, M₀[P′]) is live. For ordinary EC (= extended free choice) the last condition reduces
to "no P-component is unmarked" (no semipositive P-semiflow Y with Y·M₀ = 0). **There is no siphon condition in
the weighted case**: the authors say (p. 280, p. 296) that no satisfactory generalisation of siphons and traps to
weighted nets was available, and "in particular — and unfortunately — Commoner's theorem" is not generalised.

**Corollary 31 (decomposition) and Corollary 32 (duality theorem, p. 293) [read].** Well-formed EC ⇒ T- and
P-decomposable. "EC nets are not closed under the reverse-dual transformation, since joins need not be balanced,
while choices do"; the subclass of **balanced-joins EC nets** (every transition with several input places has equal
weights on its input arcs) is self-reverse-dual, and: **"N is a well-formed balanced-joins EC net iff N^rd is."**
Proof from Thm 28: well-structuredness is self-dual, rank Cᵀ = rank C, |𝒞^rd| = |𝒮|, and |𝒮| = |𝒞| in strongly
connected nets.

**Translation to CRNs [derived — for agent P, re-derive before use].**
- (H1) conservative and (H2) consistent of Conjecture 22 are exactly Teruel–Silva's "conservative" and "consistent";
  together = well-structured. Under the duality they swap (Γᵀ).
- Balanced joins (CRN): every source complex has all its species with the same coefficient (A+B, 2A+2B; not 2A+B).
- EC ∧ balanced joins ⇒ two species occurring in a common source complex have identical rows of Ye, i.e. they are
  **twins** in the manuscript's sense (paper.tex, l. 318). So for networks without twins, EC with balanced joins
  forces every source complex to be 0 or k·X_i (join-free). With (H3) and weak reversibility there are no twins
  (paper.tex, l. 339–344), so **a WR, EC, balanced-joins network satisfying (H3) is join-free**, and its dual is
  choice-free.
- Rank condition in CRN terms: |𝒞| = (number of distinct nonzero source complexes) + (number of reactions with
  source 0), for an EC network. For a **WR EC** network without the zero complex, distinct complexes have disjoint
  nonempty supports, so the complex matrix has independent columns and the deficiency is 0; then s = m − ℓ, and the
  rank condition s = |𝒞| − 1 = m − 1 holds iff **ℓ = 1 (one linkage class)**. Such a network is automatically
  conservative (choose w on each complex's support with w·y = 1) and consistent (WR). So by Thm 28: *a WR EC CRN
  without zero complex is well-formed iff it has a single linkage class* — the hypothesis of Anderson's 2011 GAC case.
  (Teruel–Silva assume pure nets; a WR EC network is automatically pure, since distinct complexes have disjoint
  supports. The "disjoint complexes, δ = 0" part is **known** for ordinary nets: Mairesse–Nguyen 2009, §11.)
  Without balanced joins: in a WR EC network a species occurs in one complex only, so X_i, X_j are twins iff they lie
  in the same complex with equal coefficients; (H3) + WR therefore forces pairwise distinct coefficients inside each
  complex.
- Teruel–Silva's results are about Petri-net **liveness/boundedness of the discrete token game**, not about ODE
  persistence; the bridge to persistence is Angeli–De Leenheer–Sontag (siphons), which Teruel–Silva explicitly could
  not generalise to weighted nets.

---
## 3. G. Craciun, *Toric differential inclusions and a proof of the Global Attractor Conjecture*, arXiv:1501.02860 — v2 (9 Jan 2016, 41 p.) versus v3 (23 Sep 2026, 93 p.)

**What was read.** Both PDFs from arXiv (open access; local copies in `literature/sources/1501.02860v{2,3}.pdf`,
git-ignored). v3: abstract, §1–4 in full, §8 in full, table of contents and theorem list of §5–7, references and
acknowledgements. v2: abstract, §1, §4.5 and §5 in full, section list. arXiv comment of v3 [read]: "Version 3 is a
major revision that completes and clarifies the extension from the 3D and 4D cases to the general n-dimensional
case." Submission history: v1 13 Jan 2015, v2 9 Jan 2016, v3 23 Sep 2026. Not refereed as far as the arXiv page shows
(no journal reference).

**Claimed theorems (v3, §1, p. 8) [read].**
- Theorem A (embedding, proved in Craciun, SIAM J. Appl. Alg. Geom. 3 (2019) 87–106): a variable-k weakly reversible
  system (rate "constants" k_{y→y′}(t) ∈ [ε, 1/ε]) embeds in a toric differential inclusion dx/dt ∈ F_{F,δ}(log x),
  F a complete pointed polyhedral fan in ℝⁿ.
- Theorem B: toric differential inclusions have an *exhaustive* family of zero-separating hypersurfaces (ZSH).
- Theorem C: toric dynamical systems are persistent.
- Theorem D: toric dynamical systems have a globally attracting point in each positive linear invariant subspace.
- Remark at the end of §8: any *bounded* trajectory of a variable-k weakly reversible (or variable-k endotactic,
  via Craciun–Deshpande SIADS 2020) system is persistent — a partial result towards the (Extended) Persistence
  Conjecture. v2 already stated this in footnote 52.

**Structure of the proof (v3) [read].** Persistence + LaSalle with the Horn–Jackson Lyapunov function (classical:
Siegel–MacLean, Sontag 2001). Persistence is obtained by compact forward-invariant regions K_{x₀} ∋ x₀ bounded by a
polyhedral ZSH (near the boundary of the orthant) and a level set of the Horn–Jackson function (far away) (§4,
Fig. 2). Plan (§4, p. 17): Step 0 embed the toric system T_{n+1} ⊂ ℝ^{n+1}_{>0} into a toric differential inclusion
given by a fan F_{n+1}; Step 1 build the "x₁–xₙ-subdivision" Sₙ of ℝⁿ generated by F_{n+1}; Step 2 build a ZSH
blueprint Bₙ for Sₙ; Step 3 build an exhaustive family of polyhedral ZSHs from Bₙ; Step 4 conclude. §5 (2D), §6 (3D,
first a special case, then general), §7 (subdivisions in ℝⁿ, binary words for faces, ZSH pre-blueprints and
"faithful" ZSH blueprints, with the 1D, 2D and nD cases in §7.4; Theorem 7.1 and Lemma 7.3 proved by induction),
§8 (nD case), §9 appendix (finer/coarser fans, necessary and sufficient conditions for ZSH, asymptotic level sets of
the Horn–Jackson function).

**What changed from v2 [read].** In v2 the general case is a sketch: §4.1–4.4 construct zero-separating surfaces in
dimensions 1–4 by explicit "patches" and "dotted surfaces"; §4.5 (one page, "can be carried on analogously in higher
dimensions") subdivides the simplex {X₁+…+Xₙ = −1} into n! ordered pieces and builds the surface tunnel by tunnel;
§5 (one page) excludes neighbourhoods of coordinate axes, then planes, … ("proceed as in Section 7 of [CNP 2013]").
v2 footnote 48 even notes that one construction "fails because the system of restrictions given by these alignments
is overdetermined in ℝ⁴". v3 replaces this with a systematic combinatorial machinery (blueprints, x₁–xₙ-subdivisions,
binary words, pre-blueprints → faithful blueprints, an explicit induction in §7.4.3 and §8). The ZSH notion is
changed ("different from (but still very similar with)" the identically named notion of v2, footnote 17). The
abstract changes from "prevents it from approaching the origin" (v2) to "each bounded positive solution … is contained
in an invariant region that prevents it from approaching the boundary of the positive orthant" (v3).

**Where dimension enters [read].** The induction is on the **ambient dimension** (number of species n, log
coordinates X = log x), not on the stoichiometric dimension s: along each boundary piece {x₁ ≤ … ≤ xᵢ = xᵢ₊₁ ≤ … ≤
x_{n+1}} of the fundamental region "there is a n-dimensional subproblem of the same type, and we can solve it as an
induction hypothesis" (§8, p. 80); in v2 §5 the restriction of the system near a coordinate axis is "an (n−1)-
dimensional k-variable toric dynamical system". The stoichiometric class enters only through the Horn–Jackson level
sets. The small cases cited are CNP 2013 (n = 3) and Pantea 2012 (dim of the stoichiometric subspace ≤ 3).

**Relevance.** (i) If v3 is correct the GAC is a theorem and "transferring known cases through duality" has no
target in the GAC itself; the duality question becomes one about *which proofs or certificates* transfer (e.g.
siphon-based persistence certificates, Lyapunov functions), or about non-complex-balanced extensions (Persistence
Conjecture for WR or endotactic networks, which remains open for unbounded trajectories). (ii) No mention of
duality, transposed stoichiometric matrices, siphons as such (Shiu–Sturmfels is only cited), Petri nets, or
codimension (grep of the full v3 text: "codimension" occurs only as "facets are codimension-one faces"). (iii) Status
for citation: an arXiv preprint, not peer-reviewed, 93 pages; the community has not (to my knowledge, see §11
citing works) confirmed it. Cite as "claimed".

---
## 4. C. Wiuf, *A proof of the Global Attractor Conjecture in a special case*, arXiv:2609.24553 v1 (21 Sep 2026, 17 p.)

**What was read.** Full text (arXiv PDF, `literature/sources/2609.24553.pdf`, git-ignored). The acknowledgement
states that two AI agents checked examples and computations and drafted a version that the author "completely
rewrote"; "This paper will only appear on arXiv."

**Setting (§2–3) [read].** N weakly reversible mass-action, m species. For a set Σ of species write x = (z, s) with
z = x_Σ. The **embedded network frozen at s** N_Σ(s): complexes = Σ-projections of the complexes of N, reaction
r: α_r → α′_r with rate constant k_r s^{β_r} (β_r = Σᶜ-part of the source); parallel reactions kept. Σ is a nonempty
siphon (zero sets of boundary ω-limit points are siphons, by Angeli–De Leenheer–Sontag 2007). Σ-blind linkage class:
one (hence every, by the siphon property and WR) complex has zero Σ-projection; the others are Σ-active.
Σ-minimal complex: componentwise minimal among Σ-active embedded complexes; Σ-minimal embedded linkage class: one
containing a Σ-minimal complex.
- **Condition (⋄)**: S_ℓ ∩ ℝ^{|Σ|}_{≥0} ≠ {0} for every Σ-minimal embedded linkage class L_ℓ (S_ℓ its stoichiometric
  space). By Stiemke: no Σ-minimal class has a strictly positive conservation law of its own. Remark 3.3: (⋄) for L
  ⇔ for every ω ∈ ℝ^{|Σ|}_{>0} some reaction whose source minimises a·ω over L has (α′ − α)·ω > 0 — i.e. the
  "strongly endotactic" condition of Gopalkrishnan–Miller–Shiu restricted to each Σ-minimal class and to directions
  w ∈ −ℝ^{|Σ|}_{>0}.
- S_bl(Σ) = span of reaction vectors of Σ-blind reactions (only Σᶜ coordinates); S₀(Σ) = {v : (0, v) ∈ S}, the
  stoichiometric vectors internal to the face {z = 0}. Always S_bl(Σ) ⊆ S₀(Σ) (17).
- Realisable Σ: some point of the closure of the compatibility class of x₀ has zero set exactly Σ; S(x₀) = nonempty
  realisable siphons. A positive conservation law supported on I forbids realisable Σ ⊇ I.

**Results [read].**
- Theorem 1.2 (local, §5): Σ a nonempty siphon of a WR network, s(t) → b > 0, the frozen embedded network has a
  positive complex-balanced equilibrium c, and (⋄) holds ⇒ liminf z(t) > 0. Tool: explicit Chetaev function
  U(z) = Σ_j z_j (1 − log(z_j/c_j)) (4), with a dissipation bound per linkage class (Lemmas 4.1–4.2, Prop. 4.3).
- Lemma 6.1: a boundary ω-limit point (0, b) with zero set Σ is an equilibrium, the Σ-blind subsystem is complex
  balanced at b, and log(b/s*) ∈ S_bl(Σ)^⊥ (s* = Σᶜ-part of the positive complex-balanced equilibrium).
- Prop. 6.2 / Cor. 6.3: if S_bl(Σ) = S₀(Σ) (18) then the frozen embedded network is complex balanced at every boundary
  ω-limit point with zero set Σ; (18) holds in particular if S₀(Σ) = {0} (a compatibility class meets the face in at
  most one point). Prop. 6.4: a pointwise criterion in the rate vector (Examples 6.5–6.6: Z+S₁ → 2Z+S₂ → 3Z+S₃ →
  Z+S₁, embedded Z → 2Z → 3Z → Z of deficiency 1).
- Corollary 7.1: (18) for all Σ ∈ S(x₀) ⇒ either x(t) → x* or x(t) converges to a single boundary equilibrium.
- **Theorem 7.2 (main).** N WR and complex balanced, x₀ > 0. If for every Σ ∈ S(x₀): (a) S₀(Σ) = S_bl(Σ), and
  (b) (⋄) holds for every Σ-minimal linkage class of the embedded network, then x(t) → x*.
- Examples: 8.1 (A+B ⇄ 2A, 3A+C ⇄ 3A+D, C ⇄ D; deficiency 1; shows (⋄) only on Σ-minimal classes is weaker);
  8.2 (A+C ⇄ B+C, 2A+C ⇄ 2A+D; deficiency 0; satisfies Thm 7.2, **not strongly endotactic**); 8.3 (A+B → A+C → 2A →
  A+B; single linkage class, hence strongly endotactic and covered by Anderson/GMS; **hypothesis (a) fails**: Σ = {A},
  S_bl = {0}, S₀ = span{(−1, 1)}). So Theorem 7.2 and strong endotacticity are incomparable (stated p. 3).

**Known cases listed by Wiuf (p. 1) [read]:** single linkage class (Anderson 2011), strongly endotactic (GMS 2014),
dim S ≤ 3 (Pantea 2012), first-order endotactic (Xu, J. Math. Biol. 92 (2026) 74, DOI 10.1007/s00285-026-02378-2);
persistence criteria ADS 2007, Johnston–Siegel 2011, Anderson–Shiu 2010. **Craciun v3 is not cited** (posted two days
later). No mention of duality, transposed matrices, Petri reverse-duals or codimension.

**Relevance for G [derived; check].** Wiuf's hypotheses are siphon-indexed, so under the dictionary siphon of N ↔
"T-trap" of N* (exploration brief, Fact 5) they become conditions on reaction sets of N*: Σ-blind reactions of N are
the species of N* touched by no dual reaction in Σ; S₀(Σ) and S_bl(Σ) are subspaces of im Γ, i.e. statements about
im Γ restricted to the rows Σᶜ, which in N* is the row space of Γᵀ restricted to columns Σᶜ … The objects do not
obviously become natural for N*. Note also that hypothesis (a) concerns *realisable* siphons, a notion that depends on
the conservation laws of N (= flux modes of N*).

---
## 5. M. H. T. Hack, *Analysis of production schemata by Petri nets*, MS thesis, MIT, Feb. 1972; Project MAC TR-94 (MIT-LCS-TR-094); errata: Computation Structures Note 17 (June 1974)

**What was read.** Scan from MIT DSpace (official, open: https://dspace.mit.edu/handle/1721.1/149406; local copy
`literature/sources/hack-1972-MAC-TR-94.pdf`, git-ignored, OCR text `hack.txt`): §5.2 and §5.5 in full; the rest
skimmed. **Not paywalled**, so it does not go in NEEDS.md.

**Content [read].** Hack's free choice (FC) is the strict one: every arc from a place to a transition is either the
unique input arc of the transition or the unique output arc of the place. "Deadlock" = siphon in modern terms.
§5.2 defines the reverse, the dual (place–transition interchange) and the **reverse-dual** of a net (the composite),
notes that the reverse-dual of the reverse-dual is the primal, and Lemma 6: the reverse-dual of an FC net is FC (also:
of a state machine is a marked graph; strong connectivity is preserved). Closed consistent subnets (defined by
places) are dual to open ones (defined by transitions); a dictionary "place ↔ transition, input ↔ output, covering by
strongly connected state machines ↔ by strongly connected marked graphs, SM-allocation ↔ MG-allocation". Hack remarks
that the reverse of a trap is a deadlock (siphon) but that he has no interpretation yet "for the dual or the
reverse-dual notion of a trap" (p. 71).

**Well-Formedness Theorem (Thm 9, §5.5, p. 82–83) [read].** For an FC net, equivalent: (a) well-formed (every minimal
deadlock is a strongly connected state machine and there is a covering by such); (b) there is a live and safe
marking; **(c) the reverse-dual is well-formed**; (d) every SM-reduction is a non-empty strongly connected state
machine and they cover the net; (e) the same for MG-reductions. This is the original **duality theorem**.

**Relevance.** Iosif's dual of a network with 0/1 stoichiometry and no catalysts is exactly Hack's reverse-dual
(exploration Fact 4). Hack's own words show that the reverse-dual image of traps/siphons was left without
interpretation in 1972; the modern answer is "T-traps"/semi-T-components (entry 6).

---

## 6. P. Jančar, E. Best, R. Devillers, M. Ošťádal, *Well-formed free-choice Petri nets revisited*, arXiv:2602.02352 v2 (30 Jul 2026, 30 p.)

**What was read.** Full PDF of v2 (`literature/sources/2602.02352.pdf`): §1–3 and §4 in full; appendices skimmed.

**Setting [read].** Plain (ordinary, unweighted) nets only (Remark 2.2); free choice in the extended sense
(•t₁ ∩ •t₂ ≠ ∅ ⇒ •t₁ = •t₂, Def. 3.1); clusters C(x). Thm 2.10 (well-formed ⇒ every bottom SCC is a top SCC, i.e.
strongly connected up to disjoint union); its proof "remains valid for nets with weighted arcs as well".

**Definitions [read].** Semi-T-component Y (Def. 3.5): T_Y ≠ ∅, strongly connected, |s• ∩ T_Y| = 1 for s ∈ S_Y,
T_Y• ⊆ S_Y. T-component (Def. 3.9): additionally |•s ∩ T_Y| = 1 and T_Y• = S_Y = •T_Y. Proper semi-T-component: of
Type I (a place with two input arcs in Y) or Type II (an inbound arc from outside). Reverse-dual rd(N) = (T, S, F⁻¹)
(Def. 3.14); rd(rd N) = N, strong connectivity and free choice preserved, clusters map to clusters (Prop. 3.15).
Semi-S-component X (Def. 3.16) is the reverse-dual notion: S_X ≠ ∅, strongly connected, |•t ∩ S_X| = 1 for t ∈ T_X,
•S_X ⊆ T_X. Prop. 3.17: X is a (semi-, proper semi-)S-component of N iff the corresponding subnet of rd(N) is a
(semi-, proper semi-)T-component.
- Prop. 3.7 / 3.19: semi-T-components = bottom SCCs of the subnets N_α induced by (transition-)allocations α
  (one transition per cluster); semi-S-components = top SCCs of N_β for place-allocations β.
- Prop. 3.8 / 3.20: every strongly connected FC net is covered by its semi-T-components and by its semi-S-components.
- Prop. 3.11: in a well-formed FC net every semi-T-component is a T-component.

**Main theorem [read].** **Theorem 3.23** (the duality theorem is its (a)⇔(d)): for an FC net N, equivalent:
(a) N well-formed; (b) N covered by T-components and no proper semi-T-components; (c) N covered by S-components and
no proper semi-S-components; (d) rd(N) well-formed. The introduction says the duality result was "originally
established by Hack". Prop. 3.25: at a DL-marking (every transition dead or live, one dead) of an FC net some
semi-S-component is unmarked. **Prop. 3.27**: an FC net with no proper semi-S-component of Type II is structurally
live (any marking with a token in each semi-S-component is live). *Numbering caveat:* the exploration brief cites
"Thm 3.27 (duality)"; in the v2 PDF the duality is Theorem 3.23(a)⇔(d) and 3.27 is the proposition just quoted.

**Concluding remarks (§4) [read].** For an FC net, the place sets of semi-S-components (•S_X ⊆ S_X•) "correspond to
minimal siphons"; dually, semi-T-components "correspond to what might be termed minimal T-traps, that is,
transition subsets T_Y ⊆ T satisfying T_Y• ⊆ •T_Y", and neither T-traps nor semi-T-components "appear to have
played such a prominent role in the history of T-coverability results". Appendix A: a polynomial algorithm for
well-formedness via semi-T-components; Appendix B: Commoner's theorem.

**Relevance.** The siphon ↔ T-trap dictionary of the exploration brief (Fact 5) is, for ordinary nets, stated here
(as a correspondence of minimal objects, under free choice). In CRN language: a T-trap of N* is a set Σ of species
of N such that every reaction of N producing a species of Σ consumes one — i.e. a siphon of N. The authors treat only
plain nets; nothing on weights, ODEs, persistence or chemistry.

---

## 7. J. Desel, J. Esparza, *Free Choice Petri Nets*, Cambridge Tracts in Theoretical Computer Science 40, CUP 1995 (ISBN 0-521-46519-2)

**What was read.** Author-hosted PDF (https://www7.in.tum.de/~esparza/fcbook-middle.pdf, 256 p.; local copy
`literature/sources/desel-esparza-fcbook.pdf`, git-ignored): Ch. 6 in full (p. 111–134), Prop. 5.4, Ch. 10
introduction and bibliographic notes.

**Conventions [read].** Ordinary nets; "free-choice" = what older texts call *extended* free choice (footnote 5,
Ch. 1): •t ∩ •u ≠ ∅ ⇒ •t = •u. Clusters as in entry 6.

**Theorems [read].**
- Prop. 6.1: a connected FC net is structurally live iff every proper siphon contains a proper trap; well-formed iff
  structurally live and has a positive S-invariant.
- **Theorem 6.14 (Rank Theorem)**: an FC net N is well-formed iff (a) connected with ≥ 1 place and ≥ 1 transition,
  (b) a positive S-invariant, (c) a positive T-invariant, (d) Rank(N) = |C_N| − 1 (number of clusters minus one).
  Prop. 6.15: if (a)–(c) hold and N is not well-formed then Rank(N) ≥ |C_N|.
- **Theorem 6.17**: an FC *system* (N, M₀) is live and bounded iff N is well-formed and **M₀ marks every proper
  siphon**. *Correction to the exploration brief:* the siphon condition belongs to Thm 6.17 (systems), not to the
  Rank Theorem (nets).
- Prop. 6.20: the dual net N^d = (T, S, F⁻¹) ("in some texts N^d is called the reverse-dual") has transposed incidence
  matrix, the same clusters, S-invariants ↔ T-invariants, S-components ↔ T-components; FC is preserved.
- **Theorem 6.21 (Duality Theorem)**: N is a well-formed FC net iff N^d is. Proof: a corollary of the Rank Theorem
  (all four conditions are self-dual). Bibliographic notes: proved by Hack [42]; other proofs by Desel and by
  Esparza–Silva; "Silva observed that the Duality Theorem is a corollary of the Rank Theorem". Exercise 6.6: exhibit a
  well-formed non-FC net whose dual is not well-formed (so duality fails outside FC).
- Prop. 5.4: a minimal siphon of a well-formed FC net is a trap and generates an S-component.
- Ch. 10 (Generalizations): asymmetric-choice nets (one direction of Commoner's theorem), a necessary condition for
  well-formedness of arbitrary nets (Thm 10.9, Colom–Campos–Silva, corrected by Teruel–Silva), a sufficient condition
  (Thm 10.16, Desel). Notes: the Coverability Theorems, Rank Theorem etc. "have been extended to equal-conflict
  systems by Teruel and Silva", equal conflict being the generalisation of free choice "when weighted arcs are
  allowed, first defined by Starke".

**Translation for CRNs [derived].** For an ordinary CRN (all stoichiometric coefficients 0/1) without catalysts:
condition (b) = (H1), (c) = (H2) of Conjecture 22, Rank(N) = s, and clusters ↔ classes of reactions with a common
source complex (together with the species of that complex) when N is free choice. So the Rank Equation reads
**s = (number of clusters) − 1**. Cluster conventions differ for places with empty postset: Desel–Esparza define
[x] as the least set containing x and closed under s ↦ s• and t ↦ •t, so a species consumed by no reaction (and a
reaction with source 0) is a singleton cluster; Jančar et al. (p. 9) lump all such places into one cluster. In a
well-formed net (strongly connected) there are no such nodes and both conventions give: clusters = distinct source
complexes. Desel–Esparza Prop. 6.13 / 6.15 [read, p. 128]: well-formed ⇒ Rank ≤ |C_N| − 1; (a)–(c) and not
well-formed ⇒ Rank ≥ |C_N|. For weighted CRNs see Teruel–Silva (entry 2), where EC replaces FC and the rank equation
is the same with coupled conflict sets.

---
## 8. D. Angeli, P. De Leenheer, E. D. Sontag, *A Petri net approach to the study of persistence in chemical reaction networks*, Math. Biosci. 210 (2007) 598–618, DOI 10.1016/j.mbs.2007.07.003; arXiv:q-bio/0608019

**What was read.** arXiv v2 (10 Aug 2006, `literature/sources/q-bio_0608019.pdf`): §3–5 in full, examples skimmed.
(The journal version is paywalled; the arXiv version suffices for the statements below.)

**Statements [read].** Species–reaction Petri net with weights W(S_j, R_i) = α_ij (reactant coefficient),
W(R_i, S_j) = β_ij (product coefficient) — the weighted net of exploration Fact 4. General kinetics (monotone rates
vanishing when a reactant is absent), not only mass action.
- Def. 3.1: P-semiflow = row vector c ⪰ 0 (nonnegative, nonzero) with cΓ = 0; conservative = there is a P-semiflow
  c ≫ 0. Def. 3.2: T-semiflow v ⪰ 0 with Γv = 0; consistent = there is a T-semiflow v ≫ 0. Weakly reversible ⇒
  consistent.
- Def. 3.3: siphon Σ ≠ ∅: every reaction producing a species of Σ consumes a species of Σ. Prop. 5.5: Σ is a siphon
  iff the closed face {x_Σ = 0} is forward invariant.
- **Theorem 1**: conservative and persistent ⇒ consistent (time-average argument on the ω-limit set).
- **Theorem 2**: (1) conservative and (2) every siphon contains the support of a P-semiflow ⇒ **persistent**.
  Prop. 5.4: the zero set of a boundary ω-limit point is a siphon.

**Relevance [derived].** Theorem 1 is a "dual-flavoured" statement already: (H1) + persistence ⇒ (H2). Under N ↦ N*,
(H1)(N) ⇔ (H2)(N*), and siphons of N are the T-traps of N* (exploration Fact 5), P-semiflows of N are T-semiflows of
N*. So Theorem 2 becomes, on N*: "N* consistent (H2) and every T-trap of N* contains the support of a T-semiflow
(flux mode)". No paper found that states this dual form (see §12).

---

## 9. A. Shiu, B. Sturmfels, *Siphons in chemical reaction networks*, Bull. Math. Biol. 72 (2010) 1448–1463, DOI 10.1007/s11538-010-9502-y; arXiv:0904.4529

**What was read.** arXiv v2 (5 Jan 2010, `literature/sources/0904.4529.pdf`): §1–3 in full.

**Statements [read].** Prop. 2.1 (= ADS Prop. 2): Σ siphon ⇔ face V₀(P_Σ) forward invariant. Lemma 2.2: the zero set
of a boundary steady state is a siphon; Lemma 2.3: for strongly connected G, a boundary point is a steady state iff
its zero set is a siphon. **Theorem 3.1**: minimal siphons = inclusion-minimal sets {i : c_i ∈ P}, P a minimal prime
of the ideal T_G (replaced by J_G for weakly reversible G, by a monomial ideal M_G if G is strongly connected).
Relevant (= "critical" in ADS, = not "stoichiometrically infeasible") siphons: the face F_Σ of some invariant
polyhedron P_{c(0)} = (c(0) + L_stoi) ∩ ℝ^s_{≥0} is nonempty. Theorem 3.2 computes relevant minimal siphons by
saturations; Theorem 3.5 ("extends Theorem 2 in Angeli et al."): no relevant siphon ⇔ Sat(T_G, B) = ⟨1⟩ ⇔ all siphons
satisfy property (⋆), and then no invariant polyhedron has a boundary steady state. Algorithm 3.6 uses **Alexander
duality** of squarefree monomial ideals (the only "duality" in the paper; unrelated to Γ ↦ Γᵀ). §5 mentions that
persistence of toric dynamical systems "is the content of the Global Attractor Conjecture".

**Relevance.** Gives the standard algebraic computation of (relevant) minimal siphons, useful for agent E. Nothing
on transposed networks.

---
## 10. Other "dual" constructions built on Γᵀ in the reaction-network and metabolic literature

Searched (27 Oct 2026): OpenAlex title/abstract search on "dual reaction network", "transposed stoichiometric
matrix", "dual network" + "stoichiometric matrix", "reverse-dual" + reaction, "global attractor conjecture" + duality;
web search on the same phrases; citing works of Blanchini–Giordano (OpenAlex: 4 citing works). Each item below was
opened and the relevant passage read.

**(a) I. Famili, B. Ø. Palsson, *The convex basis of the left null space of the stoichiometric matrix leads to the
definition of metabolically meaningful pools*, Biophys. J. 85 (2003) 16–26, DOI 10.1016/S0006-3495(03)74450-6
(open: PMC1303061).** [read: abstract, the passages on the compound map, figure captions.] "By analogy to the
traditional reaction map formed by S, a compound map can be derived from −Sᵀ": reactions become nodes and
metabolites become directed edges; the convex basis (extreme rays) of the left null space of S — the extreme
conservation laws, "metabolic pools" — are drawn on the compound map "in an analogous fashion to extreme pathways",
and classified into three types by analogy with flux analysis. **This is the earliest instance found of the idea
"conservation laws of N = flux modes of the transposed network", i.e. of MIPS(N) being read on N*.** No dynamics,
no siphons. (Blanchini–Giordano credit it and note the sign: −Sᵀ reverses every arrow relative to their Sᵀ.)

**(b) N. Jamshidi, B. Ø. Palsson, *Flux-concentration duality in dynamic nonequilibrium biological networks*,
Biophys. J. 97 (2009) L11–L13, DOI 10.1016/j.bpj.2009.06.049 (open: PMC2749767).** [read: abstract and main text.]
With the gradient matrix G = ∂v/∂x, "dual Jacobian matrices": J_x = S·G for concentrations and J_v = G·S for fluxes,
with the same nonzero spectrum. This is exactly the mechanism of Blanchini–Giordano's Theorem 1 (their "rates"
representation), stated 12 years earlier for a single network, without forming a dual network.

**(c) R. M. T. Fleming, N. Vlassis, I. Thiele, M. A. Saunders, *Conditions for duality between fluxes and
concentrations in biochemical networks*, J. Theor. Biol. 409 (2016) 1–10, DOI 10.1016/j.jtbi.2016.06.033;
arXiv:1512.02690.** [read: abstract and introduction.] Duality = a one-to-one map between unidirectional flux
vectors and concentration vectors (log-linear, via the source complex matrix); they give a necessary and sufficient
stoichiometric condition and a combinatorial sufficient condition ("for every two disjoint sets of molecular
species, there is at least one reaction complex that involves species from only one of the two sets"). Not a dual
network. Cites Ballerstein et al. (d) and Jamshidi–Palsson (b).

**(d) K. Ballerstein, A. von Kamp, S. Klamt, U.-U. Haus, *Minimal cut sets in a metabolic network are elementary
modes in a dual network*, Bioinformatics 28 (2012) 381–387, DOI 10.1093/bioinformatics/btr674 (open access).**
[abstract and the description in (c).] A "dual network" built from Sᵀ and the target constraints (Farkas
duality) whose elementary modes are the minimal cut sets of the primal. A different duality (cuts ↔ modes), also
built on Sᵀ; not Iosif's.

**(e) A. Duvall, *Computationally checking if a reaction network is monotone or non-expansive*, arXiv:2409.12487 v2
(27 Jan 2025).** [read: §11 "Duality" and the example.] Defines, for an *irreversible* network with stoichiometric
matrix Γ, "the dual network is the irreversible network with stoichiometric matrix Γᵀ" (reversible reactions split
first). Lemma 23: A restricted to im Aᵀ → im A is an isomorphism. **Lemma 24: if Γ is monotone with respect to a
cone K ⊆ im Γ, the dual network is monotone with respect to ΓᵀK\***; the same argument gives non-expansivity for
the dual norm ball. Explicitly "related to theorem 3" of Blanchini–Giordano. Motivation: computational (a polytope
with many vertices has a dual with few). Nothing on siphons, persistence, complex balance or the GAC.

**(f) F. Blanchini, C. A. Devia, G. Giordano, *Structural polyhedral stability of a biochemical network is
equivalent to finiteness of the associated generalised Petri net*, arXiv:2109.01709 (2021).** [abstract and the
passages mentioning duality.] A Petri net with negative tokens associated with the BDC decomposition; uses the
Blanchini–Giordano duality only to run their procedure on the dual system ẇ = CᵀDBᵀw. No reverse-dual of the
species–reaction net, no siphons.

**(g) Other works citing Blanchini–Giordano (OpenAlex):** M. A. Al-Radhawi, *Graphical characterizations of robust
stability in biological interaction networks*, Math. Control Signals Syst. (2023), DOI 10.1007/s00498-023-00350-9,
and *Graphical construction of stability certificates for biomolecular interaction networks*, CDC 2022, DOI
10.1109/CDC51059.2022.9993105; G. Nilsson, G. Giordano, *Stability analysis of biological networks through a
dynamical flow network modelling approach*, ECC 2025, DOI 10.23919/ECC65951.2025.11187128. [titles and abstracts
only, through OpenAlex.] All are about robust/structural stability certificates; none is about siphons, persistence
or the GAC.

**Conclusion of §10.** "Dual network = network with stoichiometric matrix Γᵀ" exists in three independent lines:
metabolic pools (Famili–Palsson 2003), structural stability / monotonicity (Blanchini–Giordano 2021–22, Duvall 2024),
and Iosif 2026. All transfer *linear-algebraic or local* properties (conservation pools, Jacobian spectra, polyhedral
norms, invariant cones). None transfers persistence, siphon conditions or boundary behaviour.

---

## 11. J. Mairesse, H.-T. Nguyen, *Deficiency zero Petri nets and product form*, Petri Nets 2009, LNCS 5606, 103–122, DOI 10.1007/978-3-642-02424-5_8; Fundam. Inform. 105 (2010) 237–261, DOI 10.3233/FI-2010-366; arXiv:0905.3158

**What was read.** arXiv v3 (31 May 2010): §1–2 skimmed, §3.3 in full. Companion: S. Haddad, J. Mairesse, H.-T.
Nguyen, *Synthesis and analysis of product-form Petri nets*, arXiv:1104.0291 (Fundam. Inform. 2013) [grep only].

**Statements [read].** Markovian Petri nets with race policy; Feinberg's deficiency-zero theorem gives a structural
sufficient condition for a product-form stationary distribution (Thm 2.3; also Anderson–Craciun–Kurtz). For
(ordinary, extended) free-choice nets: **Lemma 3.2**: in a weakly reversible FC net the nonempty complexes are
pairwise disjoint subsets of places, and the partition of places they induce equals the partition induced by the
clusters. Def. 3.5 / Lemma 3.3–3.4: collapsing each complex to one place gives a weakly reversible generalised state
machine with the same marking graph and the same deficiency. **Corollary 3.2: weakly reversible free-choice Petri
nets have deficiency 0.** Haddad–Mairesse–Nguyen (p. 4) recall: "for free-choice nets, WR implies deficiency zero.
The converse does not" hold.

**Relevance.** This settles, for ordinary nets, my [derived] remark in §2 (WR + EC ⇒ disjoint complex supports ⇒
deficiency 0): it is **known** (Mairesse–Nguyen 2009). The weighted (EC) version is the same two-line argument.
Consequence for P and G: in the free-choice/EC world the weakly reversible networks are deficiency-zero networks
whose complexes behave as "big places" — a very small class. [derived, unchecked] In the mass-action ODE of such a
network each species i lies in exactly one complex y (or in none), so dx_i/dt = y_i·(inflow − outflow of y) and
x_i/y_i − x_j/y_j is constant for i, j ∈ supp y; whether this gives an easy GAC proof for the class (beyond the
single-linkage-class case of Anderson 2011, which by §2 is exactly the well-formed case) is for agent G to decide.
The Mairesse–Nguyen paper is, as far as found, the only published bridge between free-choice theory and CRNT
deficiency; it does not use the reverse-dual.

---
## 12. G. Craciun, A. Dickenstein, A. Shiu, B. Sturmfels, *Toric dynamical systems*, J. Symb. Comput. 44 (2009) 1551–1565, DOI 10.1016/j.jsc.2008.08.006; arXiv:0708.3431 — the source's reference [4] for "small dimensions"

**What was read.** arXiv v2: §1, §2 (Theorems 7, 9), §3 (statement of the GAC), §5 (partial results).

**Statements [read].**
- Theorem 7: (1) is toric (complex balanced) iff the rate vector lies in the moduli variety V(M_G); then the positive
  steady states are the positive points of the toric variety V(T_G).
- **Theorem 9**: M_G is the toric ideal of the Cayley matrix Cay_G(Y) (complex matrices Y₁..Y_l of the linkage
  classes stacked over the 0/1 linkage-indicator rows), in the coordinates K₁..K_n given by the maximal minors of the
  Laplacian (tree constants, one per complex); **"The codimension of this toric ideal equals the deficiency δ"**, and
  "the deficiency δ is precisely the codimension of this moduli space".
- §3: the conjecture "is currently open, even for deficiency zero systems".
- §5: Prop. 19 (bounded polyhedron + trajectory at positive distance from the boundary ⇒ convergence to the Birch
  point); Prop. 20 (vertices of the invariant polyhedron are never ω-limit points); **Theorem 23: for detailed
  balancing systems with two-dimensional stoichiometric subspace and bounded invariant polygon, the Birch point is
  a global attractor** ("dimension two", "bounded polyhedron" and "detailed balancing" all assumed).

**Relevance — two readings of "small codimension" [derived].** The source writes: "It is well known that the Global
Attractor Conjecture holds in small dimensions [4], and we wonder whether there is a similar argument, through
duality, in small codimensions" (v2, p. 7–8). [4] is CDSS, whose "small dimension" result is Theorem 23 (s = 2).
(i) *Codimension of the stoichiometric subspace*, n − s: this is the reading of the exploration brief; under N ↦ N*
it becomes dim ker Γ(N*) = r(N*) − s, the cycle-space dimension of N*. (ii) *Codimension in the CDSS sense*: in CDSS
the only "codimension" is the codimension δ of the toric (complex-balanced) locus in parameter space (Thm 9).
Iosif's own framing of the dual (§3 of the source: the duality is "between chemical species and rate constants", and
the dual of x^{Ye} →ₖ x^{Yp} is k^{Yeᵀ} →ₓ k^{Ypᵀ}) puts the **rate constants of N in the role of the species of N***,
so the parameter space of N, where CDSS's codimension-δ toric locus lives (in tree-constant coordinates), is the state
space of N*. This second reading is a lead for agent G, not a claim; I found no paper connecting the CDSS moduli
variety with the transposed network.
(iii) *Codimension of boundary faces of the invariant polyhedron.* CDSS Prop. 20 excludes the **vertices**
(0-dimensional faces) as ω-limit points in every dimension, and Theorem 23 handles s = 2; Anderson–Shiu (SIAM J.
Appl. Math. 70 (2010) 1840–1858, arXiv:0903.0901, abstract read) prove that **facets ("codimension-one faces") are
repelling** for WR systems, hence the GAC when all boundary equilibria lie in relative interiors of facets, and for
s = 2. Under polarity of polytopes vertices ↔ facets, so "small dimension → small codimension" has a classical
realisation at the level of faces, predating the source; Pantea's s ≤ 3 closes the gap (edges). Again a lead, not a
claim; it does not involve N ↦ N*.

---

## 13. Status of the Global Attractor Conjecture and of the Persistence Conjecture, October 2026

Read for this entry: the four documents below (PDF first pages, main theorems, the remarks on GAC), and the
READMEs/catalogue entries quoted. None is refereed.

| Date | Source | Claim | Verification status found |
|---|---|---|---|
| 21 Sep 2026 | Wiuf, arXiv:2609.24553 (§4 above) | GAC for complex-balanced networks satisfying siphon conditions (a), (b) | preprint, "will only appear on arXiv" |
| 23 Sep 2026 | Craciun, arXiv:1501.02860 v3 (§3 above) | full GAC (Theorem D), via exhaustive families of zero-separating hypersurfaces for toric differential inclusions | preprint, 93 p. |
| 25 Sep 2026 | OpenAI, *Boundedness and persistence of weakly reversible mass-action systems*, github.com/openai/math, preprints/Boundedness-and-persistence-of-weakly-reversible-mass-action-systems-September-25-2026/paper.pdf | **Theorem 1.1: for every finite WR mass-action system with positive constant rates and every x₀ > 0, the solution is global and ε ≤ xᵢ(t) ≤ 1/ε for t ≥ 0 (ε depends on network, rates, x₀)** — boundedness and persistence conjectures; Remark 3.5 derives the GAC (complex balance ⇒ WR, then Horn–Jackson entropy on a compact forward-invariant convex polytope K); Prop. 3.3: each x₀ lies in a compact convex forward-invariant polytope; Cor. 3.4: short proof of Boros's existence theorem | the repository's catalogue entry 149 links a **Lean formalization**: `lean/ComparatorChallenges/MassAction.lean`, theorem `OAI.Problem326.global_bounded_persistent_solution` (statement read: WR, positive κ, positive x₀ ⇒ ∃ε ∈ (0,1), a global forward solution exists and every global forward solution stays in [ε, 1/ε]). Produced by "an unnamed internal OpenAI model"; the repo README warns that "some of the unformalized results could have issues". |
| 5 Oct 2026 | OpenAI, *Uniform permanence in weakly reversible mass-action systems*, same repository | permanence conjecture: every positive compatibility class (even unbounded) has one compact convex forward-invariant absorbing set | same catalogue family 149 |
| 10 Oct 2026 | Palomar (mmislan), github.com/mmislan/palomar-endotactic-permanence | uniform permanence for endotactic networks with measurable rates bounded above and away from 0 — "resolves the Endotactic Permanence Conjecture" (Conj. 1.3 of Johnston–Pantea–Donnell, arXiv:1412.4662), adapting the OpenAI affine barrier construction | Lean (Comparator + kernel checks recorded in the repo) |
| 11 Oct 2026 | Palomar, github.com/mmislan/palomar-quantitative-permanence | quantitative permanence and recovery bounds with vanishing reverse rates | Lean (18 declarations) |
| 27 Oct 2026 | BitMiracle AI Lab, github.com/BitMiracle-AI/detailed-balance-gac | GAC for detailed-balanced mass-action systems | Lean 4 |

The OpenAI paper (§1.2) states that Craciun's v3 "gives persistence for bounded positive trajectories of weakly
reversible systems with uniformly positive and bounded variable rates", that "the general weakly reversible statement
assumes a trajectory upper bound", and that August–Barahona's earlier general claim (their Theorem 6) has a gap (an
explicit counterexample to one implication, 2A ⇄ B). It recalls that Pantea proved persistence of bounded WR
trajectories when dim S = 2 and the GAC for dim S ≤ 3.

**Consequence for the project.** As of 27 October 2026 the GAC is claimed by at least two independent routes, one
of them with a Lean formalization of the stronger WR boundedness-and-persistence statement. Agent G should treat the
GAC as (very probably) settled, and Iosif's question as one about **which structure transfers through the duality**
(certificates, siphon conditions, Lyapunov functions, special classes), not about proving new cases of the GAC. None
of these documents mentions the transposed/dual network or preclusters.

**Classical special cases (for G's table) [abstracts read, arXiv API, 27 Oct 2026].**
- Anderson, SIAM J. Appl. Math. 71 (2011) 1487–1508, arXiv:1101.0761: GAC for a single linkage class.
- Anderson–Shiu, SIAM J. Appl. Math. 70 (2010) 1840–1858: boundary points in the relative interior of facets are not
  ω-limit points (WR, bounded); GAC when the boundary equilibria lie on facets / s = 2.
- Craciun–Nazarov–Pantea, SIAM J. Appl. Math. 73 (2013) 305–329, arXiv:1010.3050: two-species WR systems persistent
  and permanent; endotactic ⇒ permanent for n = 2 (also with time-dependent rates); GAC for n = 3 species.
- Pantea, SIAM J. Math. Anal. 44 (2012) 1636–1673, arXiv:1103.0603: lower-endotactic persistence; WR with dim S = 2:
  bounded trajectories persistent; **GAC for dim S ≤ 3**.
- Gopalkrishnan–Miller–Shiu, SIAM J. Appl. Dyn. Syst. 13 (2014) 758–797, arXiv:1305.5303: GAC for strongly endotactic
  complex-balanced systems (projection argument reducing to lower dimension).
- Deshpande–Gopalkrishnan, Bull. Math. Biol. 76 (2014), arXiv:1309.3957, DOI 10.1007/s11538-014-0024-x: every
  minimal critical siphon is drainable or self-replicable; networks without drainable siphons are persistent;
  **non-autocatalytic WR networks are persistent**.
- Xu, J. Math. Biol. 92 (2026) 74, arXiv:2409.01598: characterisation of first-order endotactic reaction graphs; every
  first-order endotactic mass-action system has a WR deficiency-zero realisation (hence GAC for them, per Wiuf p. 1).
- Johnston–Pantea–Donnell, arXiv:1412.4662: MILP tests for (strong) endotacticity; Conjecture 1.3 (endotactic
  permanence).
- CDSS 2009 (entry 12): GAC for detailed balance, s = 2, bounded polygon.

---
## 14. Persistence via siphons without conservation laws: Angeli–De Leenheer–Sontag 2009/2011

D. Angeli, P. De Leenheer, E. D. Sontag, *On persistence and cascade decompositions of chemical reaction networks*,
arXiv:0905.1332 (extended version of *Persistence results for chemical reaction networks with time-dependent
kinetics and no global conservation laws*, SIAM J. Appl. Math. 71 (2011) 128–146, and CDC 2009). [read: §3–5 statements.]

Def. 3.2: a nonempty set Σ of species is **stoichiometrically constrained** if some P-semiflow (semipositive
conservation law) has support in Σ. Def. 3.4: a set Λ of reactions is consistent if some T-semiflow has support ⊇ Λ.
**Theorem 1**: every conditionally persistent CRN is consistent (time-varying rates allowed). **Theorem 2**: if every
siphon is stoichiometrically constrained, the CRN is bounded-persistent (no conservation hypothesis). Theorem 4:
every extinction set is a siphon. Theorem 3: every conditionally persistent set of reactions is consistent.

**Dual reading [derived].** Under N ↦ N* (siphon of N = T-trap of N*, P-semiflow of N = T-semiflow of N*),
Theorem 2 reads: "if every T-trap of N* contains the support of a T-semiflow (flux mode) of N*, then N is
bounded-persistent". Theorem 1 reads: "conditionally persistent N ⇒ N* is conservative". The second is a statement
about N*, not about the dynamics of N*; no source found states either reading.

---

## 15. AI-campaign and repository scans (PROTOCOL.md §1), 27 October 2026

- `gh api users/mmislan/repos` (27 repos, Palomar): CRN-related repos are on RAFs, autocatalytic cores, D-unstable
  cores, affinity bounds, shared-resource realisations, chemical heredity, and two on persistence:
  **palomar-endotactic-permanence** (10 Oct 2026; uniform permanence of endotactic networks with measurable rates,
  claims Conjecture 1.3 of Johnston–Pantea–Donnell; Lean) and **palomar-quantitative-permanence** (11 Oct 2026; Lean).
  Both adapt OpenAI's affine-barrier construction. Nothing on duality, transposed networks, preclusters, MIPS, Petri
  reverse-duals.
- `gh api users/ipitchford/repos --paginate` (201 repos, Evidence Press): CRN-related: quadratic-mass-action-sharpness,
  reverse-reaction-attracting-oscillations, mprk-nonnegative-two-cycles, bilateral-deficiency (graph theory, not CRNT).
  Grep of names and descriptions for dual, Iosif, siphon, precluster, attractor, Petri, transpose, toric, persist,
  MIPS: **no hit on the topic**.
- github.com/openai/math (751 preprint folders, "unnamed internal OpenAI model"): two CRN preprints, *Boundedness and
  persistence of weakly reversible mass-action systems* (25 Sep 2026) and *Uniform permanence in weakly reversible
  mass-action systems* (5 Oct 2026), catalogue family 149 with a Lean statement (§13). Grep of the catalogue for
  Petri, siphon, stoichiometric, reaction network, free choice: only these.
- github.com/BitMiracle-AI/detailed-balance-gac (updated 27 Oct 2026): Lean 4 proof of the GAC for detailed-balanced
  systems.
- `gh search repos` on "global attractor conjecture", "dual reaction network", "reverse-dual", "preclusters",
  "siphon trap duality", "transpose stoichiometric matrix", "Iosif duality", "mass-action duality", "weakly
  reversible persistence" (name/description and README): only the repos above, Lean CRNT libraries
  (marpaia/crnt-lean, mdnestor/ChemicalReactionNetworkTheory), a siphon-computation script
  (curielma/Siphons-of-a-chemical-reaction-network) and this project's own public repository
  (basedpolymer/chemistry-open-problems). No repository addresses Iosif's questions.
- google-deepmind/formal-conjectures: `gh search code` for "global attractor", "reaction network", "siphon": none.
- Zenodo: the API blocks curl from this machine (HTTP 403); through the web fetcher, queries "global attractor
  conjecture", "reaction network" ∧ (dual ∨ duality ∨ siphon ∨ transpose), and (preclusters ∨ "reverse-dual" ∨ "dual
  reaction network" ∨ "transposed stoichiometric"): nothing relevant (a tropical-geometry note, a bifurcation
  dataset).
- arXiv API (sorted by date): "global attractor conjecture" (latest: Wiuf 2609.24553; Craciun v3 is listed under its
  2015 id), "reverse-dual" (only Jančar et al. among Petri/CRN papers), "dual reaction network" and "transposed
  stoichiometric" (none), siphon ∧ duality ∧ reaction (only Iosif), "free-choice" ∧ "reaction network" (only
  Mairesse–Nguyen), "equal conflict" (none), preclusters (only Iosif among CRN papers).
- Citing works: Iosif arXiv:2603.08767 / HAL hal-03058670 — OpenAlex: 0. Wiuf 2609.24553 — 0. Blanchini–Giordano —
  OpenAlex 4 (listed in §10(g)); a web search adds the survey Proverbio–Katz–Giordano, arXiv:2512.01462 (no
  duality content). Craciun 1501.02860 — OpenAlex records no 2026 citing work; the OpenAI preprint (§13) discusses v3.
  Semantic Scholar's API was rate-limited (HTTP 429) throughout; Google Scholar not reachable. **Citation searches are
  therefore incomplete for very recent (Sept–Oct 2026) works.**

---
## 16. Prior art for statements in agent P's `candidates.md` (as of 17:28) and agent G's progress log

Read-only check of teammates' files; each pointer below was verified in the cited text.

- **P-K7 (WR + EC ⇒ disjoint complexes, δ = 0)**: for ordinary nets this is Mairesse–Nguyen 2009, Lemma 3.2 and
  Corollary 3.2 (§11 above); the weighted case is the same argument. Cite them; the siphon part of K7 (minimal siphons =
  supports of extreme conservation laws, no critical siphon) was not found there.
- **P-K3 (WF-FC structure)**: (a) is Hack 1972 Thm 9 / Desel–Esparza Thm 6.21; (b) "minimal siphons = S-components"
  is Desel–Esparza Prop. 5.4(2) plus the S-coverability theorem (Thm 5.6 of the book) and Exercise 6.5 (supports of
  minimal S-invariants generate S-components); (c) by duality (Thm 5.17 of the book for T-invariants: "Let J be a
  minimal T-invariant of a well-formed FC net; the subnet generated is a T-component", and Jančar et al. Prop. 3.11);
  (f) no critical siphon: every minimal siphon of a well-formed FC net is the support of an S-component hence of a
  P-semiflow (Prop. 5.4); (h) is Thm 6.17 (M₀ = 𝟙 marks every proper siphon).
- **P-Theorem T (17:27: EC + connected + (H1) + (H2) + s = m − 1 ⇒ for every k the positive steady states are
  {x : x^{y_c} = λ u_c(k)}, toric)**: closely related to M. Pérez Millán, A. Dickenstein, A. Shiu, C. Conradi,
  *Chemical reaction systems with toric steady states*, Bull. Math. Biol. 74 (2012) 1027–1065, DOI
  10.1007/s11538-011-9685-x, arXiv:1102.1590 [read: §3]. Their Condition 3.1 asks that ker Σ (Σ = complex-to-species
  rate matrix, columns indexed by complexes) have a basis with disjoint supports; Theorem 3.3: then the steady-state
  ideal is generated by the binomials b_{j1} x^{y_{j2}} − b_{j2} x^{y_{j1}}; Theorem 3.8: positive toric steady states
  exist iff Conditions 3.4 (same-sign coefficients) and 3.6 (solvability of a log-linear system) hold; Corollary
  3.12: the same partition for all rate constants. **[derived]** If the source complexes are the only complexes with
  nonzero columns and dim ker Σ restricted to them is 1 (which the rank hypothesis s = m − 1 gives when Σ has rank s
  — to be checked), Condition 3.1 holds trivially (d = 1, plus singletons for non-source complexes); in an EC network
  the source complexes have disjoint supports, so the log-linear system of Condition 3.6 is always solvable; positivity
  of the kernel vector (Condition 3.4) is where (H2) and the EC structure must enter. So Theorem T is probably a
  special case of PMDSC Theorem 3.8 / Cor. 3.12 together with the EC structure; P should cite PMDSC and make the
  reduction explicit. The *self-duality* of Theorem T under balanced joins, and its converse, were not found in the
  literature.
- **G (16:26): "persistence is not self-dual" (3A ⇄ 2A+B globally attracting, dual 3U+2V → 2U+3V, V → U with a
  bounded extinction trajectory)**: no prior statement found; Blanchini–Giordano only state that WR and mass-action
  form are not preserved. G's "ADLS/T-trap cycle coverage" corresponds to the dual reading of ADS Theorem 2 in §8
  and §14 above (not found in print).

---
## 17. Origin of preclusters (for the paper)

Preclusters are Iosif's: thesis 2019 (`certificates/sources/thesis.txt`, Definition 4.3 and Algorithm 4.4; "by
definition any precluster is contained in a cluster"), HAL hal-03058670 v1 (11 Dec 2020) and arXiv:2603.08767
(Definition 19). The span rule comes from the **clusters J_i of C. Conradi, D. Flockerzi, *Multistationarity in mass
action networks with applications to ERK activation*, J. Math. Biol. 65 (2012) 107–156, DOI
10.1007/s00285-011-0453-1, Fact 2** [read]: a *doubling vector* u = e_ℓ − e_k (two reactions with the same educt
complex) "induces the cluster" of all row indices j with n_j ∈ span{n_ℓ, n_k}, rows n of the matrix E of extreme
flux rays. The seed of both constructions — reactions with equal source complex — is Teruel–Silva's **equal
conflict relation** (§2). No source found links preclusters or Conradi–Flockerzi clusters to Petri-net clusters.

Priority of the dual: Famili–Palsson 2003 (compound map −Sᵀ, §10a) < Iosif HAL v1 (11 Dec 2020) < Blanchini–Giordano
(IEEE L-CSS, accepted May 2021; Iosif v2 p. 8 says the two were developed independently) < Duvall 2024. The Petri
reverse-dual is Hack 1972.

---

## 18. Novelty verdict (27 October 2026)

**Searched.** Primary sources §1–17 (full texts where marked [read]); OpenAlex citing works and title/abstract
searches; arXiv API (date-sorted) on "global attractor conjecture", "reverse-dual", "dual reaction network",
"transposed stoichiometric", siphon ∧ duality, "free-choice" ∧ "reaction network", "equal conflict", preclusters;
web searches on the same phrases; the AI campaigns of PROTOCOL.md §1 (mmislan, ipitchford, Zenodo, `gh search`,
formal-conjectures) and the OpenAI math repository; Semantic Scholar and Google Scholar were not reachable, so very
recent citing works may be missed.

**Already known — must be cited, not claimed.**
1. *The dual itself.* A network with stoichiometric matrix Γᵀ: Famili–Palsson 2003 (compound map, conservation pools
   as its extreme pathways), Blanchini–Giordano 2021/22, Duvall 2024 (irreversible networks), Iosif 2020/2026. With
   arcs and weights: the Petri **reverse-dual** of the species–reaction net (Hack 1972; Desel–Esparza Prop. 6.20;
   Teruel–Silva 1996 for weighted pure nets, incidence Cᵀ; Jančar et al. 2026). That Iosif's dual equals the weighted
   reverse-dual of the Angeli–De Leenheer–Sontag species–reaction Petri net is definitional; **no paper states the
   identification** (no CRN paper on the dual mentions Petri reverse-duals; no Petri paper mentions CRN duals).
2. *Linear and local transfers.* Spectra of SG and GS (Jamshidi–Palsson 2009; Blanchini–Giordano Thm 1), structural
   local stability, polyhedral Lyapunov norms (Blanchini–Giordano Thm 3), invariant cones and non-expansivity (Duvall
   Lemma 24), conservation laws ↔ flux modes (Famili–Palsson; Iosif Thm 21).
3. *Non-preservation of (weak) reversibility and of mass-action form*: Blanchini–Giordano p. 6, Example 5.
4. *Siphon ↔ T-trap*: Hack (1972) noted that the reverse of a trap is a siphon and had no interpretation of the
   reverse-dual of a trap; Jančar et al. (2026, §4) name "minimal T-traps" (T_Y• ⊆ •T_Y) as the reverse-dual
   counterpart of minimal siphons for ordinary free-choice nets. The general statement "siphons of N = T-traps of
   N^rd" is immediate from the definitions.
5. *Free-choice / equal-conflict structure theory*: duality theorem (Hack 1972; Desel–Esparza Thm 6.21; Jančar et al.
   Thm 3.23; Teruel–Silva Cor. 32 for balanced-joins EC nets, the only weighted version), rank theorem (Desel–Esparza
   Thm 6.14; Teruel–Silva Thm 28: EC + consistent + conservative + rank C = |𝒞| − 1 ⇔ well-formed), minimal siphons
   of well-formed FC nets are traps and S-components (Desel–Esparza Prop. 5.4). Teruel–Silva state that siphons and
   Commoner's theorem have **no** satisfactory weighted generalisation.
6. *WR free-choice ⇒ deficiency zero, complexes = clusters*: Mairesse–Nguyen 2009 (Lemma 3.2, Cor. 3.2).
7. *Toric steady states from a one-dimensional kernel of the complex-to-species rate matrix*: Pérez Millán–
   Dickenstein–Shiu–Conradi 2012 (Condition 3.1, Theorems 3.3, 3.8, Cor. 3.12).
8. *GAC status*: claimed in full by Craciun (arXiv v3, 23 Sep 2026) and, via a Lean-formalized proof of boundedness
   and persistence for all WR mass-action systems, by OpenAI (25 Sep 2026, §13); special cases by Wiuf (siphon
   conditions, 21 Sep 2026) and the classical list (§13). "Small dimension" in the source = CDSS 2009 Thm 23 (s = 2);
   facets (codimension-one faces) are repelling by Anderson–Shiu 2010.
9. *Partitions by occurrence pattern in minimal invariants*: maximal common transition sets (Sackmann–Heiner–Koch
   2006, §19) = support reading of MIPS(N*); the span rule of preclusters comes from Conradi–Flockerzi 2012 (§17).

**Not found anywhere — open as far as this search can tell.**
- Any answer, partial or negative, to Iosif's Conjecture 22 or Question 24 other than result 17 of this repository;
  any paper citing Iosif's duality.
- Any use of the transposed/reverse-dual network to study **siphons, persistence, boundary behaviour, complex
  balance or the GAC** (all existing transfers are linear/local; Blanchini–Giordano list "properties other than
  stability" as future work).
- The dictionary preclusters ↔ equal-conflict sets/clusters, and any theorem relating P(N), P(N*) to Petri clusters,
  S-/T-components or minimal siphons/T-traps (agent P's K1–K3, K7-siphon part, Theorem T's self-duality).
- The dual readings of the Angeli–De Leenheer–Sontag criteria (§8, §14) and of Wiuf's hypotheses (§4).
- A precise meaning of "a similar argument, through duality, in small codimensions": three readings are recorded in
  §12 (codimension n − s ↔ cycle space of N*; CDSS codimension δ of the toric locus in parameter space, which
  Iosif's species↔rate-constants framing makes the state space of N*; codimension of boundary faces, where
  Anderson–Shiu's facet theorem is the classical "small codimension" result). None is developed in the literature in
  connection with N ↦ N*.
- Counterexamples showing that persistence/GAC-type properties are not self-dual (agent G's 3A ⇄ 2A+B example) —
  not found in print.

**Bottom line.** The ingredients (transposed network, reverse-dual, free-choice/EC duality and rank theorems,
siphon–persistence criteria, the GAC special cases) are all classical, but they live in three literatures (CRNT,
control/metabolic networks, Petri nets) that have **never been connected for this duality**. A paper from this
project can claim as new: the identification of Iosif's dual with the weighted reverse-dual and the resulting
dictionary (presented as an observation with credits to Hack, Desel–Esparza, Teruel–Silva, Jančar et al.), every
theorem relating preclusters to clusters/components/siphons, every transfer or non-transfer of persistence and GAC
hypotheses through the dual, and the precise readings of "small codimension". It cannot claim: the dual, WR
non-preservation, the spectral/Lyapunov transfers, WR FC ⇒ δ = 0, the free-choice duality and rank theorems, or the
toric-steady-state criterion, and it should not present any result as "a new case of the GAC" without stating that
the GAC is claimed in full (Craciun v3; OpenAI 2026 with Lean).

---

## 19. Prior art for MIPS-type partitions: maximal common transition sets (MCT-sets)

A. Sackmann, M. Heiner, I. Koch, *Application of Petri net based analysis techniques to signal transduction
pathways*, BMC Bioinformatics 7 (2006) 482, DOI 10.1186/1471-2105-7-482 (open: PMC1686943) [read: the MCT-set
section]; follow-up E. Grafahrend-Belau, F. Schreiber, M. Heiner, A. Sackmann, *Modularization of biochemical networks
based on classification of Petri net t-invariants*, BMC Bioinformatics 9 (2008) 90, DOI 10.1186/1471-2105-9-90
(PMC2277402) [read: passages on MCT-sets and t-clusters].

Definition (Sackmann et al., eq. (6)–(7)): transitions t_i, t_j are in the same **maximal common transition set** iff
they occur in exactly the same minimal t-invariants (identical zero patterns over all minimal t-invariants); MCT-sets
partition T, "transitions not contained in any t-invariant form their own MCT-set", and they need not be connected.

**Relevance [derived].** Iosif's MIPS(N) are classes of species with positively proportional columns in the matrix of
extreme conservation laws (manuscript Prop. 3); MIPS(N*) are therefore classes of reactions of N with positively
proportional columns in the matrix of extreme flux modes (minimal T-invariants) of N. The **support reading** of
MIPS(N*) is exactly the MCT-set partition of N (and the support reading of MIPS(N) is its place analogue); the ray
reading refines it. Agent P's K3(d) ("i ~ j iff i and j lie in exactly the same minimal siphons", for WF-FC) is the
place-side MCT construction. The paper should cite Sackmann–Heiner–Koch for the support-pattern partition.

---

## 20. Reference list (what was read)

| Ref | Identifier | Read |
|---|---|---|
| Angeli, De Leenheer, Sontag, Math. Biosci. 210 (2007) 598–618 | DOI 10.1016/j.mbs.2007.07.003; arXiv:q-bio/0608019 | arXiv v2 §3–5 |
| Angeli, De Leenheer, Sontag, SIAM J. Appl. Math. 71 (2011) 128–146 (ext. version arXiv:0905.1332) | arXiv:0905.1332 | §3–5 statements |
| Anderson, SIAM J. Appl. Math. 71 (2011) 1487–1508 | DOI 10.1137/11082631X; arXiv:1101.0761 | abstract |
| Anderson, Shiu, SIAM J. Appl. Math. 70 (2010) 1840–1858 | DOI 10.1137/090764098; arXiv:0903.0901 | abstract |
| Ballerstein, von Kamp, Klamt, Haus, Bioinformatics 28 (2012) 381–387 | DOI 10.1093/bioinformatics/btr674 | abstract + citation context |
| Blanchini, Giordano, IEEE L-CSS 6 (2022) 488–493 | DOI 10.1109/LCSYS.2021.3081369 | full |
| Blanchini, Devia, Giordano (2021) | arXiv:2109.01709 | abstract + duality passages |
| Conradi, Flockerzi, J. Math. Biol. 65 (2012) 107–156 | DOI 10.1007/s00285-011-0453-1 | Fact 2 |
| Craciun, arXiv:1501.02860 v2 (2016) and v3 (2026) | arXiv:1501.02860 | v3 §1–4, §8, TOC; v2 §1, §4.5, §5 |
| Craciun, Dickenstein, Shiu, Sturmfels, J. Symb. Comput. 44 (2009) 1551–1565 | DOI 10.1016/j.jsc.2008.08.006; arXiv:0708.3431 | §1–3, §5 |
| Craciun, Nazarov, Pantea, SIAM J. Appl. Math. 73 (2013) 305–329 | arXiv:1010.3050 | abstract |
| Desel, Esparza, *Free Choice Petri Nets*, CUP 1995 | ISBN 0-521-46519-2 | Ch. 5 (5.4, 5.6, 5.8, 5.17), Ch. 6, Ch. 10 notes |
| Deshpande, Gopalkrishnan, Bull. Math. Biol. 76 (2014) | DOI 10.1007/s11538-014-0024-x; arXiv:1309.3957 | abstract |
| Duvall (2024/25) | arXiv:2409.12487 | §11 |
| Famili, Palsson, Biophys. J. 85 (2003) 16–26 | DOI 10.1016/S0006-3495(03)74450-6 | compound-map passages (PMC) |
| Fleming, Vlassis, Thiele, Saunders, J. Theor. Biol. 409 (2016) 1–10 | DOI 10.1016/j.jtbi.2016.06.033; arXiv:1512.02690 | abstract + intro |
| Gopalkrishnan, Miller, Shiu, SIAM J. Appl. Dyn. Syst. 13 (2014) 758–797 | DOI 10.1137/130928170; arXiv:1305.5303 | abstract |
| Grafahrend-Belau, Schreiber, Heiner, Sackmann, BMC Bioinformatics 9 (2008) 90 | DOI 10.1186/1471-2105-9-90 | MCT passages |
| Hack, MIT MAC TR-94 (1972) | hdl 1721.1/149406 | §5.2, §5.5 |
| Haddad, Mairesse, Nguyen (2012) | arXiv:1104.0291 | grep |
| Iosif, *Duality in mass-action networks* | arXiv:2603.08767 v2; HAL hal-03058670 | §3 (on disk) |
| Jamshidi, Palsson, Biophys. J. 97 (2009) L11–L13 | DOI 10.1016/j.bpj.2009.06.049 | full (PMC) |
| Jančar, Best, Devillers, Ošťádal (2026) | arXiv:2602.02352 v2 | §1–4 |
| Johnston, Pantea, Donnell (2014) | arXiv:1412.4662 | abstract |
| Mairesse, Nguyen, LNCS 5606 (2009); Fundam. Inform. 105 (2010) | DOI 10.1007/978-3-642-02424-5_8; arXiv:0905.3158 | §3.3 |
| OpenAI, *Boundedness and persistence of WR mass-action systems* (25 Sep 2026) | github.com/openai/math | §1, Remark 3.5, Lean statement |
| OpenAI, *Uniform permanence in WR mass-action systems* (5 Oct 2026) | github.com/openai/math | abstract, Thm 1.1 |
| Pantea, SIAM J. Math. Anal. 44 (2012) 1636–1673 | DOI 10.1137/110840509; arXiv:1103.0603 | abstract |
| Pérez Millán, Dickenstein, Shiu, Conradi, Bull. Math. Biol. 74 (2012) 1027–1065 | DOI 10.1007/s11538-011-9685-x; arXiv:1102.1590 | §3 |
| Sackmann, Heiner, Koch, BMC Bioinformatics 7 (2006) 482 | DOI 10.1186/1471-2105-7-482 | MCT section |
| Shiu, Sturmfels, Bull. Math. Biol. 72 (2010) 1448–1463 | DOI 10.1007/s11538-010-9502-y; arXiv:0904.4529 | §1–3 |
| Teruel, Silva, Theor. Comput. Sci. 153 (1996) 271–300 | DOI 10.1016/0304-3975(95)00124-7 | full |
| Wiuf (2026) | arXiv:2609.24553 | full |
| Xu, J. Math. Biol. 92 (2026) 74 | DOI 10.1007/s00285-026-02378-2; arXiv:2409.01598 | abstract |
| Palomar (mmislan) repos; BitMiracle-AI/detailed-balance-gac | GitHub | READMEs |
