# Well-formed equal-conflict networks and their duals: siphons, uniqueness, stability, absolute consistency

Agent Q3 — **Claude Opus 5.5 (effort xhigh)**, 27 October 2026, second round on the questions left open by result 17
(`problems/iosif-siphon-precluster-duality/`). Working notes, not a public document; nothing here quotes the
author's e-mail. Builds on agent P (`../../petri/notes.md`, cited as P), agent E (`../../computation/notes.md`)
and agent L (`../../literature/notes.md`). Status tags: **proved** (complete proof here), **checked** (exact
computation on finite families, not a proof), **numerical** (floating point, not a proof), **refuted**, **open**.

## Summary

| id | Statement (N WF-EC = EC, connected, (H1), (H2), s = m − 1; no purity, no balanced joins unless said) | Status |
|---|---|---|
| C | **Theorem C without purity.** Every WF-EC network has no critical siphon; minimal siphons = supports of extreme conservation laws. Elementary proof (no Teruel–Silva) | proved; checked exactly on all 8 757 WF-EC networks of P's families (7 538 impure, 52 unbalanced) |
| CP | **Persistence on both sides.** N and N* are persistent for every rate vector (ADLS Thm 2 + C + P's C′) | proved |
| PC | **The self-dual class: proportional conflict** (sources sharing a species are proportional; Ye block diagonal with rank-one blocks). Well-formed PC networks (connected, (H1), (H2), s = #clusters − 1) and their duals (again well-formed PC) are absolutely consistent, have exactly one steady state per class, Hurwitz class-Jacobian everywhere, no critical siphon, persistence. Contains EC (a constant) and BH = EC* (b = 1) | proved; checked exactly on random genuine PC networks |
| C^gen | **Siphon theorem without equal conflict.** If N is connected, (H1), (H2) and s = #clusters − 1 (a self-dual condition), then N and N* have no critical siphon and are persistent for every rate vector; also s ≥ #clusters − 1 always (K5gen) | proved; checked exactly (2 981 networks with the rank condition, 1 229 of them non-EC; 40 908 for K5gen) |
| U | **Uniqueness and stability (Birch type).** For every k > 0, every positive compatibility class of N contains exactly one positive steady state; at **every** positive point x (not only steady states) the Jacobian restricted to S is Hurwitz. Sign form: σ(S) ∩ σ(D^⊥) = {0}, D = span{y_c − y_c′} | proved; checked (exact Routh–Hurwitz, maximal minors, direct sign vectors) |
| GS | **Global stability.** For every k > 0 the unique positive steady state of each class attracts every positive initial condition of the class. Lyapunov functions: max_c x^{y_c}/u_c(k) is nonincreasing and min_c x^{y_c}/u_c(k) nondecreasing along every trajectory | proved; identity checked exactly (claim UN), simulations §5 |
| T*′ | **Absolute consistency transfers to the dual without balanced joins**: N WF-EC ⇒ N* has a positive steady state for every κ | proved |
| T*-cex | The converse fails: N0 = {2A+B → 3A, 2A+B → A+2B} (EC, (H1)–(H3) on both sides, not WF) is not absolutely consistent, N0* is. Smallest possible. Also N1 = {2A+B → C, C → 3A, C → 3B} | proved, checked exactly |
| T*-dim | EC, connected, (H1), (H2): if s ≥ m + 1, N* is not absolutely consistent; if s = m it is not when joins are balanced; with s = m and unbalanced joins both cases occur (N0, N1 yes; NT = {X₁+X₂ → 2X₂, X₁+X₂ → 2X₃+X₄, 2X₃+X₄ → 2X₁} no) | proved, checked exactly |
| U* | **Uniqueness and stability for the dual.** For every κ > 0, every positive class of N* contains exactly one positive steady state; at every positive point of N* the Jacobian restricted to S* = im Γᵀ is Hurwitz | proved; checked exactly |
| GS*≤2 | N* is globally stable in each class when s ≤ 2; also N via the same argument (no need of GS) | proved (Bendixson–Dulac + Poincaré–Bendixson) |
| GS* | Global stability of N* for s ≥ 3 with unbalanced joins | open; numerically true on 100 random networks with s = 3, 4 (§5); balanced joins: N* is WF-EC, so GS applies |
| tor* | For unbalanced joins the steady-state set of N* is not toric (N2* example) | proved |
| K8s | K8 survives with rows and columns of Γ pairwise linearly independent on both networks, (H1)–(H3) on both sides, weakly reversible, pure (M3, M4; species side M3*, M4*) | checked exactly |
| BH | **Beyond EC through duality**: block-homogeneous networks (sources a·1_B on a partition into blocks; TEFC, not EC in general) with connected, (H1), (H2), s = #blocks − 1 have all the properties of U*, T*′, C′, CP | proved (they are exactly the duals of WF-EC networks) |
| B | Without the cluster rank equation every transfer fails: U, U*, 3A ⇄ 2A+B, E's s = 3 drain example, and new exact examples NM (TEFC, multistationary, critical siphon), NB (BH one above the rank equation, multistationary), NT (EC, dual not absolutely consistent); in E's families every non-EC network satisfying it is EC or BH up to pure catalysts (Conjecture CR for the dynamical part) | proved (examples), checked (scan) |

**In one sentence:** for well-formed equal-conflict networks the duality N ↦ N* transfers *everything that the
GAC is about* — absence of critical siphons, persistence, existence of positive steady states for all rate
constants, uniqueness per class, exponential stability — with no purity and no balanced-joins hypothesis; the
only asymmetries are toricity of the steady-state set (N always; N* with balanced joins, not in general otherwise) and the
global-attractor property, proved for N in every dimension and for N* when s ≤ 2. By duality the same properties
hold for the block-homogeneous networks (BH), a class outside EC; both classes sit in the self-dual class of
well-formed **proportional-conflict** networks (§4b), where existence, uniqueness, stability and persistence hold on
both sides and the duality acts as a symmetry. The siphon part (no critical siphon, hence
persistence, on both sides) holds much more generally, for every connected network with (H1), (H2) and the
self-dual cluster rank equation s = #clusters − 1 (Theorem C^gen); without that equation every transfer fails (§7).

**For the author's question** ("a meaningful way to dualise objects and use known cases of the GAC"): the
reverse-dual maps the well-formed equal-conflict networks — a Petri-net class that contains the weakly reversible EC
networks with one linkage class, where the GAC was classical — onto the block-homogeneous networks (Corollary BH),
which are in general neither EC, nor weakly reversible, nor complex balanced, and whose steady states are not toric.
Every dynamical property that the GAC asserts transfers along this map (existence for all rate constants,
uniqueness per class, local exponential stability, persistence), and global attraction transfers in dimension
s ≤ 2; for s ≥ 3 it is the one open point (GS*). The transfer is proved by dualising the *objects* the author
named: siphons ↔ T-traps (Theorems C, C′), the cluster/precluster seeds ↔ twin classes (the matrices G and G* are
built on them), and the conservation laws ↔ flux modes (Lemma IM and its use in T*′, U*).

## 0. Setting and standing lemmas

Notation of the manuscript and of P: species X₁..Xₙ, reactions t: y_t → y′_t, Γ = Yp − Ye, s = rank Γ, S = im Γ,
L = (im Γ)^⊥ (conservation laws). N is **EC** if two reactions sharing a reactant species have the same source;
equivalently the distinct source complexes y_1, …, y_m have pairwise disjoint supports S_1, …, S_m. T_c = set of
reactions with source y_c. Under (H2) every species is consumed (P §0), so every species lies in exactly one S_c;
write c(p) for it and a_p = (y_{c(p)})_p ≥ 1. For k ∈ ℝʳ_{>0} let A_k be the n × m matrix with columns
a_c(k) = Σ_{t∈T_c} k_t (y′_t − y_c). Mass action: ẋ = f(x) = A_k φ(x), φ_c(x) = x^{y_c}.
**WF-EC**: EC, connected, (H1), (H2), s = m − 1. No purity (species may occur on both sides of a reaction) and no
balanced joins are assumed anywhere unless said.

Two facts about the sign pattern of A_k are used throughout (EC disjointness):

(Z) for p ∈ S_c and d ≠ c: A_k[p, d] = Σ_{t∈T_d} k_t (y′_t)_p ≥ 0, and A_k[p, d] > 0 iff some reaction of T_d
produces p ("d → c" in the cluster graph);
(Σ) A_k 1 = Γk (sum of the columns), so for v > 0 in ker Γ, A_v 1 = 0.

**Cluster graph.** c → d iff some reaction of T_c produces a species of S_d. By P's Lemma SC (connected + (H1) +
(H2) ⇒ the Petri net is strongly connected) the cluster graph is strongly connected: a directed path
t₁ → p₁ → t₂ → … of the Petri net maps to c(t₁) → c(p₁) = c(t₂) → ….

**Lemma CF (P).** N EC, every species consumed, strongly connected. If u ≥ 0, u ≠ 0 and A_k u = 0 then u > 0 and
ker A_k = ℝu. (P's proof is valid without purity: in row p the term of the own cluster may have any sign, only (Z)
is used.)

**Theorem T (P).** N WF-EC. For every k > 0, ker A_k = ℝu(k) with u(k) > 0; the positive steady states are
E_k = {x > 0 : x^{y_c} = λu_c(k) ∀c, λ > 0}. (P's proof needs no purity.)

**Lemma IM.** N WF-EC, k > 0. Then im A_k = S. *Proof.* im A_k ⊆ S (each column is Γ applied to k·1_{T_c});
rank A_k = m − 1 by Theorem T; s = m − 1. ∎ Consequently the left kernel of A_k is exactly L, for every k.

**Lemma M (Perron–Frobenius for Metzler matrices).** Let G be an irreducible Metzler matrix (off-diagonal entries
≥ 0) with Gw = 0 for some w > 0 (or ξᵀG = 0 for some ξ > 0). Then 0 is a simple eigenvalue of G, every other
eigenvalue has negative real part, and G has nonnegative left and right null vectors (both positive).
*Proof.* For c > max|G_ii|, M = G + cI is nonnegative, irreducible, with positive diagonal, hence primitive;
Mw = cw with w > 0 forces ρ(M) = c (a nonnegative irreducible matrix has a positive eigenvector only for its
spectral radius). Primitivity: ρ(M) is simple and every other eigenvalue μ of M has |μ| < c, so λ = μ − c has
Re λ < 0. The left Perron vector of M is a positive left null vector of G. ∎
(If G is Metzler with Gw = 0, w > 0, but possibly reducible, the left kernel still contains a nonzero ξ ≥ 0: apply
Perron–Frobenius to Mᵀ, which has spectral radius c.)

## 1. Theorem C without purity (task 1)

**Theorem C.** Let N be WF-EC (no purity, no balanced joins). Every siphon of N contains the support of a nonzero
z ≥ 0 with zᵀΓ = 0. Hence the minimal siphons of N are exactly the supports of the extreme conservation laws, and
N has no critical siphon.

*Proof.* Let R be a minimal siphon contained in the given siphon. All species of S_c have the same consuming
reactions T_c, so by Theorem A of P (transversality, no hypothesis) R meets each S_c in at most one species; let
C_R = {c : R ∩ S_c ≠ ∅} and p_c the species of R ∩ S_c, so c ↦ p_c is a bijection C_R → R. Fix v > 0 in ker Γ
(H2) and put A = A_v, B = A[R, C_R] (rows p_c, columns c ∈ C_R; square).

1. *B is Metzler*: B[p_c, d] = A[p_c, d] ≥ 0 for d ≠ c by (Z).
2. *A vanishes on R outside C_R*: for p ∈ R and d ∉ C_R, A[p, d] = Σ_{t∈T_d} v_t (y′_t)_p (p ∉ S_d); if a
   reaction t ∈ T_d produced p ∈ R, the siphon property would give a reactant of t in R, i.e. R ∩ S_d ≠ ∅,
   d ∈ C_R. So A[p, d] = 0.
3. *B has zero row sums*: by 2 and (Σ), B1 = (A1)_R = (Γv)_R = 0.

By Lemma M (reducible form) there is z ≥ 0, z ≠ 0, on R with zᵀB = 0. Extend z by 0. Then zᵀA = 0: on columns
of C_R this is zᵀB = 0, on the other columns it is 2. So z lies in the left kernel of A_v, which is L by Lemma IM:
zᵀΓ = 0. The consequence (minimal siphons = supports of extreme laws) follows exactly as in P's proof of Theorem C
(the support of a nonnegative conservation law is a siphon). ∎

Remarks. (i) This proof uses neither purity nor Teruel–Silva (P's proof went through Thm 28, Thm 15 and Cor. 13,
stated for pure nets). The only inputs are transversality, the sign pattern (Z), Lemma CF and the rank equation.
(ii) Where the rank equation enters: z ⊥ im A_v always (for any EC network with (H2)); z ⊥ im Γ needs
im A_v = im Γ. On non-WF EC networks the certificate fails on every network tested (logs: "C [nonWF] 0/…").
(iii) A second proof, from P's Theorem C for pure nets: split every reaction t with a self-loop into
t₁: y_c → Q_t, t₂: Q_t → y′_t (new species Q_t). This keeps EC, connectedness, (H1) (z_Q = z·y_c), (H2), purity
becomes true, and rank and number of clusters both increase by the number of split reactions, so the rank equation
is kept. If R is a siphon of N, R ∪ {Q_t : t ∈ •R split} is a siphon of the split network; a conservation law z′
supported there restricts to a nonzero conservation law of N supported in R (z′_Q = z′·y_c = z′·y′_t). The direct
proof above is preferred (self-contained).

**Corollary CP (persistence on both sides).** If N is WF-EC then N and N* are persistent for every choice of rate
constants (mass action, and every kinetics satisfying the assumptions of Angeli–De Leenheer–Sontag 2007, §2).
*Proof.* N is conservative and has no critical siphon (Theorem C); N* is conservative (N consistent) and has no
critical siphon (P's Theorem C′, no purity, no balanced joins). ADLS Theorem 2. ∎ This removes the hypothesis
"pure or balanced joins" from P's Corollary P. For mass action, §3 gives a second, siphon-free proof for N.

*Checked:* claim C of `scripts/check_wfec_families.py` (exact: Metzler, zero row sums, nonnegative left kernel
vector, conservation law with support exactly R, for every minimal siphon) on all WF-EC networks of agent P's six
EC families: (2,3,5,0) 71 (62 impure), (3,2,5,0) 958 (793 impure), (3,3,6,1) 7, (4,2,5,1) 34, (4,2,4,0) 692 (436
impure), (3,3,4,0) 6 995 (6 230 impure, 52 unbalanced): 8 757 networks, all PASS (`logs/wfec_*.log`). E had
already found no failure on 2 116 impure sides (K14c); E2 confirms independently (`../computation/logs/q3_tests_*.log`).

### 1.1 The siphon theorem without equal conflict: the cluster rank equation

For an arbitrary network let C_1, …, C_q be its Petri clusters (components of the bipartite graph species–reactions
with an edge for each reactant; Teruel–Silva's coupled conflict sets), with species sets S_c and reaction sets T_c
(under (H1), (H2) every cluster has both). Let Ā_k be the n × q matrix with columns Γ(k ⊙ 1_{T_c}). The sign pattern
(Z) holds for Ā_k for **every** network (a species of S_c is not a reactant of any reaction outside T_c), and
Ā_k 1 = Γk.

**Lemma CF′.** N with every species consumed and strongly connected. If u ≥ 0, u ≠ 0, Ā_k u = 0, then u > 0 and
ker Ā_k = ℝu. *Proof.* Row p ∈ S_c reads (Σ_{t∈T_c} k_t (y_t)_p) u_c = Σ_d (Σ_{t∈T_d} k_t (y′_t)_p) u_d with a
positive coefficient on the left (p is consumed by some reaction of T_c); so u_d > 0 and d → c imply u_c > 0, the
support of u is closed in the strongly connected cluster graph, and the kernel argument of Lemma CF applies. ∎

**Proposition K5gen.** N connected, (H1), (H2) ⇒ s ≥ q − 1. (Lemma SC, Lemma CF′ with k = v > 0 in ker Γ and u = 1:
rank Ā_v = q − 1 ≤ s.) Checked on all 18 220 connected networks with (H1), (H2) of E's general families.

**Theorem C^gen (no equal conflict).** Let N be connected with (H1), (H2) and s = q − 1 (*cluster rank
equation*). Then every siphon of N contains the support of a nonzero z ≥ 0 with zᵀΓ = 0, and so does every siphon of
N*. Hence N and N* have no critical siphon, their minimal siphons are the supports of their extreme conservation laws
(for N*: the minimal T-traps of N are its internal cycles), and both are persistent for every rate vector (ADLS).

*Proof.* Let R be any siphon, C_R the clusters whose species meet R, v > 0 in ker Γ, B = Ā_v[R, C_R] (rows R,
columns C_R; not necessarily square). (1) For p ∈ R ∩ S_c and d ≠ c, B[p, d] ≥ 0 (Z). (2) Ā_v[p, d] = 0 for p ∈ R,
d ∉ C_R (a reaction producing p consumes a species of R, so its cluster is in C_R). (3) B1 = (Γv)_R = 0. If no
z ≥ 0, z ≠ 0 had zᵀB = 0, Gordan's alternative would give w with Bw > 0 componentwise; take c* maximising w over C_R
and p ∈ R ∩ S_{c*}: (Bw)_p = Σ_d B[p, d](w_d − w_{c*}) ≤ 0 by (1) and (3), a contradiction. So such z exists; by (2)
zᵀĀ_v = 0 (z extended by 0), and im Ā_v = S because rank Ā_v = q − 1 = s (Lemma CF′). So zᵀΓ = 0. For N*: its
clusters are those of N with species and reactions exchanged, its rank is s, it is connected, and (H1), (H2) swap
(P, Prop. 1.2); so N* satisfies the same hypotheses. ∎

Remarks. (i) This contains Theorem C (EC: clusters = source classes, q = m) and P's Theorem C′ (C^gen for N*), and
gives both without Theorem T, Teruel–Silva, minimality or transversality. (ii) **The class "connected + (H1) + (H2) +
s = #clusters − 1" is self-dual** (EC is not; EC with balanced joins is): it is the natural "well-formed" class for
the siphon part of the duality. It contains WF-EC, the block-homogeneous networks of §7.1 (their duals) and the
well-formed PC networks of §4b. (Adding a catalyst can merge clusters and break the equation; the catalyst remark
of §7.1 treats that case separately.) (iii) Every counterexample of §7.2 violates the cluster rank equation; in particular the
cluster rank equation, not the source-count equation s = m_src − 1, is what persistence of both sides needs.
(iv) Checked (`scripts/check_cgen.py`, exact): on E's general families, 1 250 EC and 755 non-EC networks with the
cluster rank equation, no critical siphon on either side, and the constructive certificate (B with zero row sums and
the sign pattern, a nonnegative left kernel vector found by exact Fourier–Motzkin, for every minimal siphon of N and
of N*) holds (`logs/check_cgen_all.log`); on seven weakly reversible families ((3,2,6), (4,3,4), (5,3,4), (5,2,5),
(6,2,5) homogeneous, (2,2,5) and (3,2,5) with zero complex): 502 EC and 474 non-EC networks with the cluster rank
equation, all PASS, and K5gen on 22 688 networks (`logs/check_cgen_wr.log`; the run on all WR files exceeded 30 min).
(v) What C^gen does not give outside EC: absolute consistency, uniqueness, toricity (the mass-action field no
longer factors through one monomial per cluster). Conjecture CR (§7.3) asks for them.

## 2. Uniqueness and local stability for N (task 3)

**Theorem U.** Let N be WF-EC with mass-action kinetics and k > 0.
(a) Every positive stoichiometric class (x⁰ + S) ∩ ℝⁿ_{>0} contains exactly one positive steady state.
(b) At every positive point x (steady state or not), det(tI − J(x)) = t^{n−s} q(t) with q Hurwitz; J(x) maps
ℝⁿ into S and J(x)|_S has all its eigenvalues in the open left half-plane. In particular every positive steady
state is exponentially stable relative to its class, and no Hopf bifurcation can occur, for any k.
(c) σ(S) ∩ σ(D^⊥) = {0}, where D = span{y_c − y_d} (dim m − 1 = s) and σ denotes sign vectors.

*Proof of (c) and of "at most one" in (a).* Let x, x′ be positive steady states for k with x − x′ ∈ S. By
Theorem T, φ(x) = λu, φ(x′) = λ′u, so d = log x − log x′ satisfies y_c·d = δ := log(λ/λ′) for all c: d ∈ D^⊥ up
to the common value δ (indeed d ∈ D^⊥ means exactly y_c·d independent of c). Write x − x′ = A_v w (Lemma IM,
v > 0 in ker Γ). Let c* maximise w. For p ∈ S_{c*}, by (Σ),
(x − x′)_p = Σ_d A_v[p, d](w_d − w_{c*}) ≤ 0
since A_v[p, d] ≥ 0 for d ≠ c* (Z) and w_d ≤ w_{c*}. So x_p ≤ x′_p on S_{c*}, hence δ = Σ_{p∈S_{c*}} a_p d_p ≤ 0.
At an argmin cluster the same computation gives δ ≥ 0. So δ = 0, and for every argmax cluster c the sum
Σ_{p∈S_c} a_p d_p = 0 of nonpositive terms forces x = x′ on S_c, so 0 = (x − x′)_p = Σ_d A_v[p, d](w_d − w_c) with
all terms ≤ 0: w_d = w_c for every d → c. The argmax set is closed under predecessors in the strongly connected
cluster graph, so w is constant and x − x′ = A_v w = 0.
The same argument, read with s ∈ S and d ∈ D^⊥ having the same sign vector (x, x′ replaced by s, d), proves (c):
s = A_v w; at argmax c*, s ≤ 0 on S_{c*}, so d ≤ 0 there and y_{c*}·d ≤ 0; at argmin, ≥ 0; the common value
y_c·d is 0, so d = 0 = s on S_{c*} and the propagation gives w constant, s = 0. ∎
(For mass action, (c) is equivalent to "at most one point of E_k per class for every k" by the classical
sign-vector argument — Müller–Regensburger 2012, Prop. 3.9 — since u(k) ranges over all of ℝᵐ_{>0} as k varies:
k_t = v_t/u_{c(t)} realises any u.)

*Proof of (b).* J = A_k Dφ(x) with Dφ(x)[c, p] = a_p φ_c(x)/x_p for p ∈ S_c, 0 otherwise. Put
G(x) = Dφ(x) A_k (m × m). For c ≠ d, G[c, d] = Σ_{p∈S_c} (a_p φ_c/x_p) A_k[p, d] ≥ 0 by (Z), and > 0 iff d → c;
so G is Metzler and irreducible (cluster graph strongly connected). And **G(x)u(k) = Dφ(x)A_k u(k) = 0 for every
x**. By Lemma M, 0 is a simple eigenvalue of G and the other m − 1 eigenvalues have negative real part. By
det(tIₙ − A_k Dφ) = t^{n−m} det(tI_m − Dφ A_k), J has 0 with algebraic multiplicity n − m + 1 = n − s and its
other eigenvalues are those of G. J maps into im A_k = S and induces 0 on ℝⁿ/S (dimension n − s), so the
eigenvalues of J|_S are the m − 1 = s nonzero eigenvalues of G. ∎

*Existence in (a).* The class P = (x⁰ + S) ∩ ℝⁿ_{≥0} is compact (H1), convex and forward invariant, so the flow has
an equilibrium x̄ ∈ P (Brouwer applied to the time-τ maps, τ = 2^{−j}, and a limit). If x̄ were on the boundary, its
zero set Z would be a siphon (for i ∈ Z, ẋ_i = Σ_t rate_t(x̄)(y′_t)_i = 0, so every reaction producing i has a
reactant in Z), which contains the support of some z ≥ 0, z ≠ 0, zᵀΓ = 0 (Theorem C); but z·x̄ = z·x⁰ > 0.
So x̄ > 0. (Also from §3.) ∎

**Generalisation (cluster kinetics).** Everything above uses mass action only through "all reactions of T_c have
rates k_t φ_c(x) with φ_c depending on the species of S_c". For any C¹ φ_c ≥ 0 depending on S_c only, strictly
increasing in each of them on the interior and vanishing iff one of them vanishes, (a) and (b) hold (in the
uniqueness proof replace "δ" by the common ratio φ_c(x)/φ_c(x′) = λ/λ′ and use strict monotonicity). This covers
equal-conflict Petri nets with routing rates under any common "server" function per conflict set.

*Checked* (exact, `scripts/check_wfec_families.py`, `scripts/check_dual_families.py`): (c) by the maximal-minor
criterion of Müller–Feliu–Regensburger–Conradi–Shiu–Dickenstein (2016) on all WF-EC networks of the families and by
direct enumeration of sign vectors (exact Fourier–Motzkin) for n ≤ 4: PASS everywhere (claims U1, U1d); (b) by
exact characteristic polynomial and Routh array at x = 1 and two random rational steady states (claim U2) and at
two random rational *non-steady* points with random k (claim UN, which also checks G Metzler and G u(k) = 0): PASS
on all 8 757 WF-EC networks of P's six families (`logs/wfec_*.log`, `logs/dual_*.log`).
On non-WF EC networks U1d and U2 fail on every network tested (degenerate steady states), as the proof predicts
(J has 0 with multiplicity n − m + 1 > n − s when s ≥ m).

**So WF-EC networks are never multistationary.** E2's search for multistationarity in WF-EC (brief E2, task 3) is a
red-team test of Theorem U.

## 3. Global stability of WF-EC networks (consequence of the structure of §2)

**Theorem GS.** Let N be WF-EC with mass-action kinetics and k > 0; u = u(k) > 0 as in Theorem T. For a positive
solution x(t) put ψ_c(t) = x(t)^{y_c}, M(t) = max_c ψ_c/u_c, μ(t) = min_c ψ_c/u_c. Then M is nonincreasing, μ is
nondecreasing, and x(t) converges to the unique positive steady state of its class. In particular N is persistent
(every coordinate is bounded below by a positive constant along each positive trajectory), with no siphon argument.

*Proof.* ψ̇ = Dφ(x)ẋ = Dφ(x)A_kφ(x) = G(x)ψ with G(x) Metzler and G(x)u = 0 (§2). If c attains M(t),
ψ̇_c = Σ_d G[c, d]ψ_d ≤ Σ_{d≠c} G[c, d] M u_d + G[c, c] M u_c = M (G u)_c = 0, so the upper Dini derivative of M is
≤ 0; similarly μ is nondecreasing. Hence μ(0) ≤ ψ_c(t)/u_c for all c, t; with x bounded above (H1) and every
species in some S_c, every x_p(t) is bounded below by a positive constant: the trajectory stays in a compact set
K of the open orthant, inside its class. Let W = log M − log μ ≥ 0 (W = 0 iff φ(x) ∈ ℝ_{>0}u iff x is a steady
state, since A_kφ(x) = 0 iff φ(x) ∈ ker A_k = ℝu). W(x(t)) decreases to a limit W*, and on the ω-limit set Ω ⊆ K
(compact, invariant) W ≡ W*. Suppose W* > 0 and take a trajectory x̄(t) in Ω: along it log M − log μ is constant
with M nonincreasing and μ nondecreasing, so M ≡ M̄ is constant. The functions h_c = ψ_c/u_c are real-analytic in t
(mass action), so {t : h_c(t) = M̄} is discrete unless h_c ≡ M̄; these sets cover [0, ∞), so I = {c : h_c ≡ M̄} is
nonempty. For c ∈ I, 0 = ψ̇_c = −Σ_{d≠c} G[c, d](M̄u_d − ψ_d) with nonnegative terms, and G[c, d] > 0 for d → c
in the interior: so every predecessor d of c is in I. Strong connectivity gives I = all clusters, ψ ≡ M̄u, W = 0:
contradiction. So W* = 0 and Ω consists of steady states of the class, i.e. Ω = {x*} by Theorem U. ∎

Remarks. (i) The Lyapunov functions M and μ are piecewise linear in the *cluster rates* ψ_c and depend on k
through u(k); they are of the type "piecewise-linear in rates" of Blanchini–Giordano (2014) and Al-Radhawi–Angeli
(2016), but rate-dependent, so the robust (kinetics-independent) results of those papers do not apply as such; this
is consistent with the fact that with independent kinetics per reaction a conflict set no longer splits its flux in
fixed ratios. (ii) WF-EC networks are in general neither weakly reversible nor complex balanced (W*, F, B_k of P),
so Theorem GS is a global-attractor statement outside the complex-balanced setting. It re-proves the GAC for
weakly reversible EC networks without the zero complex: their complexes are pairwise disjoint (P's Prop. F), so
each linkage class lives on its own species and is a WF-EC network (connected, conservative, consistent, δ = 0,
s = m − 1), to which GS applies. (iii) The same proof works for cluster
kinetics with analytic φ_c. (iv) Novelty: not found in the sources read by L (§18 of their notes) — see §8.

## 4. The dual: absolute consistency, uniqueness, stability (tasks 2, 3)

N* has species k_t (t = 1..r) and one reaction ρ_p per species p of N: ρ_p consumes a_p of every k_t, t ∈ T_{c(p)},
and produces Σ_t (y′_t)_p k_t. With π_c = Π_{t∈T_c} k_t, the mass-action rate of ρ_p is ω_p = κ_p π_{c(p)}^{a_p},
and k̇ = Γᵀω. A positive k is a steady state iff ω ∈ L. N* is EC iff N has balanced joins (P, Prop. 1.2).

### 4.1 Absolute consistency without balanced joins

**Theorem T*′.** Let N be WF-EC. For every κ ∈ ℝⁿ_{>0}, N* has a positive steady state. More precisely, for every
M > 0 there is π > 0 with Σ_c π_c = M and ω(π) ∈ L; then every k > 0 with Π_{t∈T_c} k_t = π_c is a steady state.

*Proof.* Fix v > 0 in ker Γ. By Lemma IM, ω ∈ L iff H_c(π) := ω(π)·a_c(v) = 0 for all c. (i) Σ_c H_c(π) =
ω·Γv = 0 identically. (ii) For d ≠ c, ∂H_c/∂π_d = Σ_{p∈S_d} a_c(v)_p a_p κ_p π_d^{a_p−1} ≥ 0 by (Z): H is
cooperative. (iii) If π_c = 0 then ω_p = 0 on S_c (a_p ≥ 1), so H_c(π) = Σ_{p∉S_c} ω_p a_c(v)_p ≥ 0. (iv)
H_c(π) ≥ −V_c Σ_{p∈S_c} a_p κ_p π_c^{a_p} ≥ −Cπ_c on the simplex Δ_M = {π ≥ 0, Σπ_c = M} (V_c = Σ_{t∈T_c} v_t).
For 0 < ε ≤ 1/C, π ↦ π + εH(π) maps Δ_M continuously into itself; Brouwer gives π* ∈ Δ_M with H(π*) = 0. If
Z = {c : π*_c = 0} were nonempty, for c ∈ Z the identity 0 = H_c(π*) = Σ_{p∉S_c} ω_p(π*) a_c(v)_p (all terms ≥ 0)
gives π*_{c(p)} = 0 for every p produced by T_c: Z is closed under successors, so Z = all clusters, contradicting
Σπ* = M > 0. ∎

**Corollary.** For N EC, connected, (H1), (H2): N absolutely consistent ⇒ N* absolutely consistent (N absolutely
consistent ⇒ s = m − 1 by P's Theorem T(2) ⇒ Theorem T*′). With balanced joins this is P's T*, an equivalence.

**Proposition T*-cex (the converse fails; smallest example).** N0 = {2A+B → 3A, 2A+B → A+2B} is EC (one source),
unbalanced, satisfies (H1)–(H3), and so does N0* = {2K₁+2K₂ → 3K₁+K₂, K₁+K₂ → 2K₂}; s = m = 1. N0 is not
absolutely consistent: A_k = (k₁ − k₂)(1, −1)ᵀ has a nonzero kernel only if k₁ = k₂ (ẋ_A = (k₁ − k₂)x_A²x_B). N0* is:
K̇₁ = κ_A(K₁K₂)² − κ_B K₁K₂ vanishes on K₁K₂ = κ_B/κ_A, e.g. at (1, κ_B/κ_A). Minimality: an unbalanced source needs
two species and molecularity ≥ 3; (H2) needs r ≥ 2. A second example with two clusters: N1 = {2A+B → C, C → 3A,
C → 3B} (s = m = 2; positive steady states only for k₂ = 2k₃), N1* = {2K₁ → 3K₂, K₁ → 3K₃, K₂+K₃ → K₁} with the
positive steady state (κ_B/κ_A, 3κ_AK₁²/κ_C, 1) for every κ. (`scripts/check_examples.py` X0, X1: exact.)

**Proposition T*-dim.** N EC, connected, (H1), (H2). (i) If s ≥ m + 1, N* is not absolutely consistent. (ii) If
s = m and joins are balanced, N* is not absolutely consistent. (iii) If s = m and joins are unbalanced, both cases
occur: N0, N1 have absolutely consistent duals; NT = {X₁+X₂ → 2X₂, X₁+X₂ → 2X₃+X₄, 2X₃+X₄ → 2X₁} (EC, unbalanced
source 2X₃+X₄, (H1)–(H3) on NT and NT*, s = m = 2) does not: every conservation law of NT has z₁ = z₂, while the
reactions ρ₁, ρ₂ of NT* have the same source K₁+K₂, so ω₁/ω₂ = κ₁/κ₂ at every point and NT* has no positive steady
state when κ₁ ≠ κ₂. With one cluster and molecularity 4: NTc = {2A+B+C → 2A+2B, 2A+B+C → 2A+2C}. Minimality: with
molecularity ≤ 3, n = 3 forces either one cluster (then the extra equation of S^⊥ involves two exponents — needs
molecularity 4) or m = 2, where L is a line and ω(π) ∝ z is always solvable; r = 2 forces m = 1. So NT (n = 4,
r = 3) is smallest with molecularity ≤ 3. (`scripts/check_examples.py` X3, X3'.) Agent E2 found N* absolutely
consistent for all 1 989 unbalanced EC networks with s = m of the family (3,3,4,0): the failure needs the extra
equation to live on species of one exponent, as in NT.
*Proof.* The κ admitting a positive steady state are κ_p = w_p π_{c(p)}^{−a_p}, w ∈ L ∩ ℝⁿ_{>0}, π ∈ ℝᵐ_{>0}: a
semialgebraic image of a set of dimension (n − s) + m. (i) This is < n, so the image has measure zero. (ii) With
balanced joins (w, π) ↦ (tw, t^{1/a}π) leaves κ fixed, so the image has dimension ≤ n − 1. (iii) Examples. ∎
So for EC networks absolute consistency of N* is decided by s − m and the joins, except in the unbalanced case
s = m, where it depends on the network (a characterisation there is open).

### 4.2 Uniqueness and stability of the dual

**Theorem U*.** Let N be WF-EC (balanced joins or not) and κ > 0. (a) At every positive point k of N*,
det(tI − J*(k)) = t^{r−s} q(t) with q Hurwitz: J*(k) restricted to S* = im Γᵀ has all its eigenvalues in the open
left half-plane. (b) Every positive class of N* contains exactly one positive steady state, exponentially stable
relative to the class.

*Proof.* (a) ∂ω_p/∂k_t = a_p ω_p/k_t for t ∈ T_{c(p)}, 0 otherwise, so J* = Γᵀ diag(aω) E_S E_Tᵀ diag(1/k) = P*Q*,
with E_S (n × m) and E_T (r × m) the species- and reaction-to-cluster incidence matrices. G* = Q*P* =
E_Tᵀ diag(1/k) Γᵀ diag(aω) E_S has entries G*[c, d] = Σ_{t∈T_c} k_t^{−1} Σ_{p∈S_d} Γ[p, t] a_p ω_p, ≥ 0 for d ≠ c
(Γ[p, t] = (y′_t)_p) and > 0 iff c → d: irreducible Metzler. For ξ ∈ ℝᵐ,
ξᵀG* = (A_{1/k} ξ)ᵀ diag(aω) E_S, where A_{1/k} is the matrix A of N with rate constants 1/k_t. By Theorem T
(applied to N with these rates) ξ = u(1/k) > 0 gives ξᵀG* = 0. Lemma M: 0 is simple, the other eigenvalues have
negative real part. As in §2, J* has 0 with multiplicity r − m + 1 = r − s and maps into S* (dimension s = m − 1),
so J*|_{S*} carries the m − 1 nonzero eigenvalues of G*. (b) N* is conservative (N consistent), so each class P* is
a compact convex polytope of dimension s; it has no steady state on its boundary (the zero set would be a siphon of
N*, which contains the support of a conservation law of N* by P's Theorem C′); every positive steady state k* is
nondegenerate, with index sign det(J*(k*)|_{S*}) = (−1)^s by (a) (s eigenvalues with negative real parts). The
homotopy λF + (1 − λ)(k̂ − k) (k̂ an interior point of P*, F the vector field) has no zero on ∂P* for λ ∈ [0, 1]
(F(k) lies in the tangent cone of P* at boundary points and is nonzero there, k̂ − k points strictly inward), so the
sum of the indices of the zeros of F in P* equals the degree of k ↦ k̂ − k, which is (−1)^s. Every zero has index
(−1)^s, so there is exactly one. (Equivalently: Conradi–Feliu–Mincheva–Wiuf, PLoS Comput. Biol. 13 (2017) e1005751,
Thm 1.) Existence also follows (degree ≠ 0). Stability from (a). ∎

*Checked* (exact, `scripts/check_dual_families.py`, claim Ustar): G* Metzler, u(1/k) a positive left null vector,
det(tI − J*) = t^{r−s}q(t) with q Hurwitz (exact Routh array), at k = 1 and two random rational points with random
rational κ, on all 8 757 WF-EC networks of P's families (52 unbalanced) and on the unbalanced networks of
`dual_search.py` (n = 3): all PASS (`logs/dual_*.log`).

Remark. Theorem U* for N is Theorem U(b) without the toric argument: the same degree argument proves U(a) from U(b)
and Theorem C. The sign-vector proof of U(a) is kept because it is elementary and gives the equivalent condition (c).

### 4.3 Toricity and global stability of the dual

**Proposition tor*.** With unbalanced joins the steady-state set of N* need not be toric: N2 = {2A+B → C, C → 2A+B}
is WF-EC and unbalanced; N2* = {2K₁ → 2K₂, K₁ → K₂, K₂ → K₁} has steady states K₂ = (2κ₁K₁² + κ₂K₁)/κ₃, whose
logarithm is not an affine line. With balanced joins N* is WF-EC and Theorem T makes it toric (P's T*).

**Proposition GS*≤2.** If N is WF-EC and s ≤ 2, N* is globally asymptotically stable in each positive class.
*Proof.* s = 1: a scalar field on a segment with a unique, attracting zero (U*), persistent (CP). s = 2: the class
is a convex polygon; the divergence of the field restricted to it is tr(J*|_{S*}) < 0 at every interior point
(U*(a)), so by Bendixson–Dulac there is no periodic orbit and no cycle of orbits in the (simply connected) interior;
ω-limit sets lie in the interior (persistence, CP); by Poincaré–Bendixson and the uniqueness of the equilibrium,
which is a hyperbolic sink, every ω-limit set is {k*}. ∎

**Open (GS* for s ≥ 3, unbalanced).** Numerically every simulated trajectory converges (§5). The Hurwitz property
at every point (U*(a)) does not by itself imply global stability for s ≥ 3 (Markus–Yamabe fails in dimension ≥ 3);
the structure that gave GS for N (a constant positive right null vector of G) is replaced in N* by a k-dependent
left null vector u(1/k) of G*, and the Hilbert-metric argument does not transfer directly.

## 4b. The self-dual class behind these results: proportional conflict

**Definition.** N is *proportional-conflict* (PC) if two reactions sharing a reactant species have positively
proportional source complexes. Equivalently (sharing is transitive inside a Petri cluster, proportionality is
transitive): on every Petri cluster c (species S_c, reactions T_c) the block Ye[S_c, T_c] is a rank-one matrix
b_c a_cᵀ with b_c ∈ ℕ^{S_c}_{>0} primitive and a_c ∈ ℕ^{T_c}_{>0} (integrality of a_t follows from primitivity of
b_c), and Ye vanishes outside the blocks. **PC is self-dual**: Ye(N*) = Yeᵀ has the transposed rank-one blocks on the
same clusters. EC = PC with a_c constant on each cluster; BH (§7.1) = PC with b_c = 1 (the duals of EC); balanced EC
= both. Write θ_c = x^{b_c}; mass action reads w_t = k_t θ_c^{a_t} (t ∈ T_c), ẋ = Γw.

**Lemma R (routing; the cluster version of P's Theorem T(1)).** N connected, (H1), (H2), s = q − 1 (q = number of
clusters). For every r ∈ ℝʳ_{>0}, the n × q matrix Ā_r with columns Γ(r ⊙ 1_{T_c}) has ker Ā_r = ℝu(r) with
u(r) > 0, and im Ā_r = S. *Proof.* P's closed–open argument for Theorem T(1), with Lemma CF′ in place of Lemma CF:
rank Ā_r ≤ s < q; the set of r with a nonzero nonnegative kernel vector is closed, contains r = v (u = 1), and is open
by Lemma CF′ and semicontinuity of the rank. ∎

**Theorem PC.** Let N be PC, connected, (H1), (H2), s = q − 1 ("well-formed PC"). Then N* is well-formed PC, and for
both N and N*, for every rate vector: (a) no critical siphon, persistence (Theorem C^gen); (b) absolutely consistent;
(c) at every positive point the Jacobian restricted to the stoichiometric subspace is Hurwitz; (d) every positive
class contains exactly one positive steady state, exponentially stable.
*Proof.* Self-duality and (a): above and §1.1. (b) Let z > 0 be a conservation law. By Lemma R applied to N* (with
routing z), im Γᵀ is spanned by the vectors Γᵀ(z ⊙ 1_{S_c}), so w ∈ ker Γ iff H_c(θ) := Σ_{p∈S_c} z_p (Γw(θ))_p = 0
for all c. Then Σ_c H_c ≡ zᵀΓw = 0; H is cooperative (for d ≠ c, ∂H_c/∂θ_d = Σ_{p∈S_c} z_p Σ_{t∈T_d}
(y′_t)_p a_t k_t θ_d^{a_t−1} ≥ 0); H_c ≥ 0 when θ_c = 0, and H_c ≥ −Cθ_c on a simplex. Brouwer on
{θ ≥ 0, Σθ_c = M} gives a zero θ*, positive by strong connectivity exactly as in T*′; then x with x^{b_c} = θ*_c
exists because the b_c have disjoint supports. (c) J = PQ with P (n × q) the columns Γ((a∘w) ⊙ 1_{T_c}) and Q (q × n)
the rows (b_c/x)ᵀ; G = QP is Metzler (G[c, d] = Σ_{p∈S_c} (b_p/x_p) P[p, d], P[p, d] ≥ 0 for p ∉ S_d), irreducible,
and Pη = Ā_{a∘w}η = 0 for η = u(a∘w) > 0 (Lemma R): by Lemma M and the multiplicity count of §2 (J has 0 with
multiplicity n − q + 1 = n − s and maps into S), J|_S is Hurwitz. (d) Degree argument of U*(b), with (a) excluding
boundary steady states and (H1) giving compact classes. N*: apply (a)–(d) to N*. ∎

Remarks. (i) Theorems T(1) (existence part), U(a)(b), T*′, U*, Corollary BH are the special cases a constant (EC) and
b = 1 (BH). The toric description (T) and the global attractor property (GS) use a constant on clusters, i.e. EC;
for general PC networks, GS is open for s ≥ 3 (for s ≤ 2 it follows from (c) as in Prop. GS*≤2). (ii) So the
duality acts as a symmetry of the class of well-formed PC networks, exchanging the "shape" vectors b_c and the
"multiplicity" vectors a_c; this is the precise sense in which the Petri-net reverse-dual transfers persistence,
absolute consistency, uniqueness and stability. (iii) Checked (`scripts/check_pc.py`, exact): on the example
NP = {A+2B → C, 2A+4B → C+D, C → A+2B, D → A+2B} (PC, neither EC nor BH) and on random well-formed PC networks that
are neither EC nor BH (NP plus 4 random networks with n = 4, `logs/check_pc_n4.log`, and 3 with n = 5,
`logs/check_pc_n5.log`; such networks are rare among random draws): N* well-formed PC, G Metzler with positive null
vector u(a∘w), J|_S Hurwitz at two random points, on N and N*, and no critical siphon on either side. A small sample:
this part of the evidence is thin, the proof is what carries the claim.
(iv) In E's general families every non-EC network with the cluster rank equation is EC or BH up to pure catalysts
(§7.3); genuine PC networks need molecularity ≥ 6 (two proportional sources with a non-0/1 shape).

## 5. Numerical exploration of the dual (task 3, red team of U*)

`scripts/dual_search.py` (floating point, not a proof): all unbalanced WF-EC networks with source 2X₁+X₂ and the
other species in sources with coefficients ≤ 2, molecularity ≤ 3; for each, sampled κ and sampled classes of N*,
25 Newton-type starts; distinct positive solutions counted, Jacobian eigenvalues on S* computed.
- n = 3, r ≤ 3: 22 networks, 1 solution in every sampled class, all stable.
- n = 3, r ≤ 5: 23 networks (only one extra network with r ∈ {4, 5} passes WF-EC), 207 (network, κ, class)
  triples: exactly 1 solution each, all stable (`logs/dual_search_n3_r5.log`).
- n = 4, r ≤ 4 (`logs/dual_search_n4_r4.log`): the run hit its 30-minute limit before printing its summary; it printed
  no case with two solutions or an unstable solution, and 50 "no solution found" cases. Integrating N* from the
  class point in each of these 50 cases (`scripts/nosol_recheck.py`, `logs/nosol_recheck.log`) converges to a
  positive steady state of the class (residual ≤ 2e-15): failures of the multistart solver, not of T*′.
These agree with Theorem U*. Independent exact confirmation by agent E2 (`../computation/logs/ms_dual_3_3_4.log`,
msolve): for the 52 unbalanced WF-EC networks of the family (3,3,4,0), 208 sampled (κ, class) pairs of N*, exactly one
positive solution each. E2 also found no multistationary WF-EC network N in its families (`ms_wfec_*.log`), and
re-checked Q3-C, U(b), U(c) and T*′ exactly (`q3_tests_*.log`), as Theorems C, U, T*′ say.

**Simulations (`scripts/gs_sim.py`, floating point).** 60 random unbalanced WF-EC networks with n = 5 (all with
s = 3) and 40 with n = 6 (all with s = 4), random k, κ; for each, three initial points in one class of N and three
in one class of N*, LSODA to t = 2000 (rtol 1e-10). N: the largest relative increase of M(t) along the computed
trajectories is 1e-11 and the largest relative decrease of μ(t) 6e-12 (rounding level), as Theorem GS says; the
three end points coincide to 1e-14 and are steady states. **N* (the open case GS*, s = 3 and 4, unbalanced):** the
three end points coincide to ≤ 8e-13 and |f| ≤ 1e-15 at them, in all 100 networks. No periodic or non-convergent
behaviour was seen (`logs/gs_sim_n5.log`, `logs/gs_sim_n6.log`).

## 6. K8, strong form (task 4) — answered: the obstruction survives

**Proposition K8s (checked exactly).** There are pairs of weakly reversible networks, each satisfying (H1)–(H3)
together with its dual and having the rows of Γ pairwise linearly independent and the columns of Γ pairwise linearly
independent, with the same doubling classes and the same minimal T-traps (up to a common relabelling of the
reactions) but different preclusters P(N) — for **every** relabelling that matches the two keys. Pure example on one
species set (n = 4, r = 6, s = 3, deficiency 1):
- M3 = X₄ → 2X₃, X₄ → X₂+X₃, 2X₃ → 2X₂, X₂+X₃ → 2X₁, 2X₂ → X₄, 2X₁ → 2X₃;
  doubling classes {R₁,R₂} + singletons; minimal T-traps {R₁,R₃,R₅}, {R₂,R₃,R₅}, {R₄,R₆}; P(M3) = {R₁,R₂,R₅}, {R₃},
  {R₄}, {R₆};
- M4 = X₄ → X₂+X₃, 2X₃ → X₄, 2X₃ → 2X₂, X₂+X₃ → 2X₁, 2X₂ → X₄, 2X₁ → 2X₃;
  doubling classes {R₂,R₃} + singletons; minimal T-traps {R₁,R₂,R₅}, {R₁,R₃,R₅}, {R₄,R₆}; P(M4) = {R₂,R₃} +
  singletons.
(8 relabellings match the keys; none matches the preclusters.) Smallest example found (n = 3 and 4, r = 5, impure):
M1 = 2X₃ → 2X₂, 2X₃ → X₁+X₂, 2X₂ → X₁+X₃, X₁+X₃ → 2X₃, X₁+X₂ → 2X₂ and M2 = 2X₄ → X₃, X₃ → X₂+X₄, X₃ → X₁+X₄,
X₂+X₄ → 2X₄, X₁+X₄ → 2X₄ (minimal T-traps are two autocatalytic self-loop reactions; P(M1) has the doubling pair as
only non-singleton block, P(M2) is one block).
**Species side** (the form asked by P): by the dictionary (P, Prop. 1.2) N₁ = M3* and N₂ = M4* are two networks on
the same 6 species with (H1)–(H3), pairwise independent rows and columns, the same twin classes and the same minimal
siphons, and P(N₁*) = P(M3) ≠ P(M4) = P(N₂*) under every matching relabelling. So the profile rule of P's Theorem
D(e) has no analogue outside the well-formed free-choice class, even under the strongest non-degeneracy hypotheses,
and even for weakly reversible networks (reaction side).
Search: `scripts/k8_strong.py` over all 431 148 distinct networks of E's record files (31 892 strong with (H1)–(H3);
reaction side: 130 860 colliding labelled keys; species side: none in the records, whose species-side twins are
rare — answered by duality instead); with `PURE=1` (1 285 strong pure networks): 7 875 colliding keys on the
reaction side. Exact verification: `scripts/check_k8s.py` (all PASS, `logs/check_k8s.log`).

## 7. Beyond equal conflict (task 5)

### 7.1 How far the theory reaches: WF-EC and its dual class

What the proofs use: the sign pattern (Z) (a species is a reactant of one source complex only) in every argument
(Metzler structure of A_k, G, G*; cooperativity of H), the rank equation through Lemma IM; never purity nor balanced
joins. The duality itself gives a second class, outside EC, where everything except toricity and (for s ≥ 3) the
global attractor property holds.

**Corollary BH (block-homogeneous networks; a class beyond EC reached through the duality).** Call K
*block-homogeneous* if its species are partitioned into blocks B_1, …, B_m and every source complex is a·1_B (a
positive integer multiple of the indicator vector of one block), every block being the support of some source. If K is
block-homogeneous, connected, (H1), (H2) and s = m − 1, then: K is absolutely consistent; every positive
compatibility class contains exactly one positive steady state, exponentially stable relative to the class (the
Jacobian restricted to S is Hurwitz at every positive point); K has no critical siphon and is persistent; K is
globally stable in each class when s ≤ 2. K is EC iff the multiples a are constant on each block (then K is WF-EC
with balanced joins); otherwise K is TEFC (sources sharing a species have the same support) but not EC.
*Proof.* K* is EC: two reactions ρ_q, ρ_q′ of K* share a reactant t (a reaction of K) iff q, q′ lie in the block B(t)
consumed by t; then for every reaction t′ of K, q ∈ B(t′) ⟺ q′ ∈ B(t′), so rows q and q′ of Ye(K) coincide, i.e. ρ_q
and ρ_q′ have the same source. The clusters of K* are the blocks, connectedness and the rank are the same, (H1) and
(H2) swap (P, Prop. 1.2): K* is WF-EC with m clusters. As K = (K*)*, Theorems T*′, U*, P's C′ (with ADLS) and
Prop. GS*≤2 applied to N = K* give the claims. ∎
Examples: N0* (not WF, s = m, so not covered) and N1*, N2*: N2* = {2K₁ → 2K₂, K₁ → K₂, K₂ → K₁} is
block-homogeneous with blocks {K₁}, {K₂}; N1* = {2K₁ → 3K₂, K₁ → 3K₃, K₂+K₃ → K₁} is block-homogeneous with blocks
{K₁}, {K₂, K₃} but s = 2 = m (not covered; indeed N1 is not WF). So **the transfers go beyond EC exactly in the
direction of TEFC networks whose conflict sets consume a block uniformly** — Teruel–Silva's TEFC class restricted to
"homogeneous" weights, the reverse-dual image of EC.
(Open, conjecture: the common generalisation "TEFC, every source of the conflict set of block B a multiple a·b_B
of one primitive vector b_B with support B" contains both classes (b_B arbitrary and one a per block: EC; b_B = 1_B:
BH). Its Jacobian has the same P·Q factorisation with a Metzler m × m matrix G; a positive null vector of G would
follow from a Theorem-T-type argument for the linear problem "(ξ_B b_q/x_q)_q ∈ L", which I sketched but did not
write; absolute consistency and the siphon part were not attempted.)

**Remark (pure catalysts).** If some species are pure catalysts (zero rows of Γ) and deleting them gives a WF-EC
network N̂, then in each class of N the catalyst concentrations are constant and the dynamics is that of N̂ with
rate constants multiplied by catalyst monomials; Theorems C (a siphon containing a catalyst E contains supp e_E),
T, U, GS transfer to N, and N* differs from N̂* by trivial reactions y → y. Catalysts that are also consumed or
produced elsewhere break the cluster factorisation and are not covered.

### 7.2 Counterexamples bounding the class

Two "rank equations" make sense outside EC: s = m_src − 1 (m_src = number of distinct source complexes, the count
used for EC) and s = m_cl − 1 (m_cl = number of Petri clusters = Teruel–Silva's coupled conflict sets; for EC the
two counts agree, for block-homogeneous networks m_cl = number of blocks). Keeping connectedness and (H1)–(H3) on
both N and N* (true for every network of the table; for NM and NB: NM* = k₃ → k₁, 2k₃ → k₁+k₂, 2k₁+3k₂ → 2k₂+3k₃ and
NB* = k₃ → k₂, k₃ → k₁, k₁+2k₂ → 3k₃, each with two critical siphons, as N), each conclusion fails outside EC ∪ BH:

| conclusion | network where it fails (connected, (H1)–(H3)) | m_src, m_cl, s | source |
|---|---|---|---|
| absolute consistency of N | U* = K₃ → K₂, K₃ → K₁, K₂ → K₁, K₁+K₂ → 2K₃: no positive steady state for κ₂ ≥ κ₁ | 3, 2, 2 | P, T-cex |
| transfer of absolute consistency N → N* | U = X₄ → X₂+X₃, X₃+X₄ → X₁, X₁+X₂ → 2X₄ | 3, 2, 2 | manuscript, P |
| no critical siphon in N*; persistence of N* | 3A ⇄ 2A+B: N globally stable, N* = {3U+2V → 2U+3V, V → U} has the drain species V (extinction for all rates) | 2, 1, 1 | G, E (Prop. E.4) |
| same, intermediate dimension | X₄ → 3X₃ → X₁+2X₂ → 3X₂ → X₄: N* loses k₄ for every rate vector | 4, 3, 3 | E |
| uniqueness per class; no critical siphon in N | NM = 2X₃ → X₁+X₂, 3X₃ → X₂+2X₃, X₁+2X₂ → 3X₃ (TEFC, not EC, not BH): two positive steady states in the class x₁+x₂+x₃ = 5/2 for k = (1, 2, 1); critical siphon {X₁, X₃} | 3, 2, 2 | here, X4 |
| uniqueness per class, BH one above the rank equation | NB = X₃ → X₂, 2X₃ → X₁, X₁+X₂ → 3X₃ (BH, blocks {X₃}, {X₁,X₂}): two positive steady states in the class 2x₁+x₂+x₃ = 6 for k = 1 | 3, 2, 2 | here, X5 |
| absolute consistency of N* (EC, one above the rank equation) | NT = X₁+X₂ → 2X₂, X₁+X₂ → 2X₃+X₄, 2X₃+X₄ → 2X₁ | 2, 2, 2 | here, X3 |

Every one of these networks satisfies s = m_src − 1 or s ≥ m_cl but **violates the cluster rank equation
s = m_cl − 1**. (NM: the class contains exactly the points x₃ = 1/2, x₁x₂² = 1/4, x₁ + x₂ = 2, i.e. the two roots of
x₂³ − 2x₂² + 1/4 in (0, 2); NB: x₃ = 1, x₁x₂ = 1, 2x₁ + x₂ = 5. Exact sign certificates in `scripts/check_examples.py`.)

### 7.3 Scan of the general families (`scripts/beyond_ec.py`, `scripts/cluster_rank.py`)

On agent E's general families (all networks with (H1), (H2) in (2,2,5), (2,2,6), (2,3,4), (3,2,4), (3,3,3), (4,1,12),
(4,2,3)), restricted to connected networks with (H1)–(H3):
- **with s = m_src − 1, outside EC**: all 687 TEFC networks and 9 413 of the 9 442 other networks have a critical
  siphon (and 687 / 9 412 a critical siphon in N*); a multistart Newton search (numerical) finds two positive steady
  states in one class for 85 of 300 TEFC and 97 of 300 other networks tried (`logs/beyond_ec.log`; NM above is an
  exact instance). So the source-count rank equation is not the right hypothesis outside EC.
- **with s = m_cl − 1, outside EC**: 737 networks (503 BH, 64 TEFC, 170 other). None has a critical siphon in N or
  in N*, and the multistart search finds no second steady state in any of them (all 300 + 64 + 170 tried). Reason
  (`logs/cluster_rank.log`): after deleting pure catalysts (species with a zero row of Γ) every one of them is EC or
  BH (BH 13 + 490, TEFC → BH 58, TEFC → EC 6, other → BH 129, other → EC 41; genuinely new: 0), so they are covered by
  Corollary BH and the catalyst remark.

**Conjecture CR (for the next round).** Let N be connected with (H1), (H2) and s = m_cl − 1. Then N (and N*) is
absolutely consistent and has exactly one, exponentially stable, positive steady state in each class, for every rate
vector; perhaps because such N is always, up to pure catalysts, EC or "block-power" (every source of a cluster a
multiple of one primitive vector b_B). Status: the siphon part is **proved** (Theorem C^gen, §1.1); the rest is true
on all 737 + 1 198 networks of the families above (all EC or BH up to catalysts); the block-power case with b_B not
0/1 needs molecularity ≥ 6 and is not in the families.

## 8. Credits and novelty (to be checked by the verifier and by agent L)

- Theorem T (toric steady states) — P, with priors Pérez Millán–Dickenstein–Shiu–Conradi 2012 and
  Campos–Chiola–Silva 1991 (P §7).
- Lemma M — classical Perron–Frobenius. Sign-vector criterion — Müller–Regensburger 2012 (Prop. 3.9); maximal-minor
  form — Müller, Feliu, Regensburger, Conradi, Shiu, Dickenstein, Found. Comput. Math. 16 (2016). Degree argument —
  standard (Conradi–Feliu–Mincheva–Wiuf 2017, Thm 1).
- Donnell–Banaji, *Local and global stability of equilibria for a class of chemical reaction networks*, SIAM J.
  Appl. Dyn. Syst. 12 (2013) 899–920, arXiv:1211.2153 (read here, `certificates/sources/`, git-ignored): local stability of every
  positive equilibrium and global stability for networks with a factorisation Γ = ΛΘ (each row of Λ with one nonzero
  entry, columns of Θ with one negative and one positive entry, ker Θᵀ one-dimensional with a positive vector),
  general kinetics, via monotone systems with increasing integrals. Their class forces the rows of Γ of species of
  one group to be proportional; WF-EC networks do not satisfy this in general (e.g. 2A+B → C+D type products into
  two clusters), and our results are rate-structure (cluster kinetics) specific. Related in spirit (monotone
  structure, Lyapunov functions from Perron–Frobenius), different class and proof.
- Piecewise-linear-in-rates Lyapunov functions: Blanchini–Giordano (Automatica 2014), Al-Radhawi–Angeli (IEEE TAC
  2016) — robust versions; ours (Theorem GS) is rate-dependent.
- Node-balanced steady states (Feliu–Cappelletti–Wiuf, Math. Biosci. 2018, arXiv:1711.07074; abstract read):
  uniqueness and stability for steady states balanced at the nodes of a "reaction graph" whose nodes carry
  complexes; a WF-EC network with a product complex that is not a source complex has a node without outgoing edge,
  so its steady states are not node balanced: not covered.
- New here, as far as checked: Theorem C without purity (elementary proof), Theorems U, GS, T*′, U*, Corollary BH,
  the counterexamples N0, N1, NT, Prop. T*-dim, K8s, the boundary §7. The literature search on "equal conflict" /
  "disjoint reactant complexes" + mass action + uniqueness / global stability was short (two web searches, the two
  papers above); a verifier should repeat it.

## 9. What is left

- GS* / GS for proportional-conflict networks with non-constant multiplicities and s ≥ 3 (global stability of N*
  for unbalanced WF-EC N); numerically true (100 random networks with s = 3, 4).
- Conjecture CR (§7.3): absolute consistency and uniqueness for every connected network with (H1), (H2) and the
  cluster rank equation (the siphon part is Theorem C^gen); possibly such networks are always PC up to catalysts.
- Characterisation of absolute consistency of N* when N is EC with s = m and unbalanced joins (both cases occur).
- Novelty: a proper literature search for Theorems U, GS, C^gen and PC (two short web searches only; Donnell–Banaji
  2013 and node balancing do not cover them).
- Larger exact families for PC networks that are neither EC nor BH (they need molecularity ≥ 6).

## Files

- `scripts/eclib.py` — library (imports crn.py and P's petri.py): EC structure, A_k, Theorem C certificate,
  sign condition (minors and direct), exact Jacobian, characteristic polynomial, Routh–Hurwitz.
- `scripts/check_examples.py` — X0, X1, X2, X3, X3′, C0, U0 (exact; `logs/check_examples.log`).
- `scripts/check_wfec_families.py` + `scripts/run_families.sh` — claims C, U1, U1d, U2 on P's six EC families
  (`logs/wfec_*.log`).
- `scripts/check_dual_families.py` — claims Ustar (U*(a) and the left null vector u(1/k) of G*) and UN (U(b) and
  G u(k) = 0 at non-steady points) (`logs/dual_*.log`).
- `scripts/dual_search.py` + `scripts/run_dual.sh` — numerical search for multistationarity/instability of N*
  (`logs/dual_search_*.log`).
- `scripts/k8_strong.py` — K8 strong form search in E's records (`logs/k8_strong.log`, `logs/k8_strong_pure.log`);
  `scripts/check_k8s.py` — exact check of the witnesses (`logs/check_k8s.log`); `scripts/k8s_search.py` — random
  species-side search (none found, `logs/k8s_search_n4.log`; answered by duality instead).
- `scripts/gs_sim.py` — simulations for GS and GS* (`logs/gs_sim_n5.log`, `logs/gs_sim_n6.log`).
- `scripts/beyond_ec.py` — scan of the general families outside EC (`logs/beyond_ec.log`); `scripts/cluster_rank.py`
  — classification of the non-EC networks with the cluster rank equation (`logs/cluster_rank.log`).
- `scripts/check_pc.py` — Theorem PC on random genuine PC networks (`logs/check_pc_n4.log`, `logs/check_pc_n5.log`).
- `scripts/check_cgen.py` — Theorem C^gen and K5gen on E's records, with the constructive certificate on N and N*
  (`logs/check_cgen_all.log`, `logs/check_cgen_wr.log`).
- `scripts/run_round2b.sh` — batch: beyond_ec, k8s_search, check_dual_families on the families.
- `candidates.md` (for E2), `PROGRESS.md`. Sources: `certificates/sources/` (Donnell–Banaji arXiv PDF and text; git-ignored).
