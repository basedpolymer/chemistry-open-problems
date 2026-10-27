# Verification V1: Petri-net dictionary (agent P), corrected preclusters (agent Q2), generalised readings (agents E, E2)

Verifier V1, **Claude Opus 5.5 (effort xhigh)**, fresh context, 27 October 2026. Briefs: `../../BRIEF-V-common.md`,
`../../BRIEF-common.md`, `../../BRIEF-V1.md`; protocol `PROTOCOL.md` §4. Prior-art input from agent L2
(`../../literature2/novelty.md`) was checked against the sources where possible (§4).

Scope: `petri/notes.md` (all sections except §7, which is V2's; results of §7 used by my claims are flagged),
`round2/preclusters/notes.md` (all numbered results and the table of witnesses), `computation/notes.md` §3 (agent E)
and Proposition E2.6 of `round2/computation/notes.md`.

Method. I re-derived every proof in scope by hand. Independently of the authors' code (`crn.py`, `petri.py`,
`pclust.py`, `ecore.py` are **not** imported by my checks) I wrote an exact library (`scripts/v1lib.py`, Fractions,
extreme rays by minimal supports, siphons and T-traps by brute force, Definition 19 and the thesis Algorithm 4.4
implemented separately, κ_V, π, P_cs) and checked every example exactly, plus exhaustive families and seeded random
families. Symbolic certificates use sympy 1.14 (WSL). The authors' example scripts were re-run only for comparison
(they agree). All sources cited by the notes were opened on disk (PDF pages rendered where the text extraction was
garbled).

## 0. Verdicts at a glance

| Claim | Verdict |
|---|---|
| P 1.1 (Iosif's dual = weighted reverse-dual) | CONFIRMED |
| P 1.2 (dictionary) | CONFIRMED |
| P 2.1, 2.2, Cor. 2.3, 2.4 | CONFIRMED (2.3: implicit hypotheses made explicit) |
| P Remark 2.5 (B₂, B₃) | CONFIRMED |
| P Theorem A (transversality) | CONFIRMED (FC case is Desel–Esparza Thm 4.30(a)) |
| P Theorem B (exact criterion for Q3) | CONFIRMED WITH CORRECTION (definition of twins; degenerate case) |
| P Theorem C (pure WF-EC: no critical siphon) | CONFIRMED (ordinary case is classical; weighted case new) |
| P Lemma SC | CONFIRMED |
| P Theorem C′ (duals of WF-EC: no critical siphon) | CONFIRMED (uses Theorem T(1) of §7; see §1.9) |
| P Theorem D (WF-FC structure, (a)–(h)) | CONFIRMED (one citation imprecise) |
| P K8 (obstruction pairs) | CONFIRMED (checked exactly; add caveat on M₂) |
| P Theorem E (i), (ii) | CONFIRMED (failure lists reproduced exactly) |
| P table "counterexamples in Petri terms" | CONFIRMED for EC/balanced/ordinary/WF/crit columns ("crit" counts all critical siphons); abs. cons. column = V2 |
| P Proposition H | CONFIRMED |
| P Proposition G and the refutations by D, N₄.₁₂ | CONFIRMED |
| P Proposition F | CONFIRMED (with "no zero complex", as restricted) |
| P Corollary P | CONFIRMED WITH CORRECTION (one sentence of the proof) |
| Q2 Lemma 1.2 | CONFIRMED (part (b) rests on a reading of Conradi–Flockerzi 2012, stated explicitly) |
| Q2 Props 2.1, 2.2 and the table of strict inclusions | CONFIRMED (exact, with symbolic certificates) |
| Q2 Theorem 3.1 (invariance of CF) | CONFIRMED |
| Q2 Props 3.2, 3.3 | CONFIRMED |
| Q2 Theorems 4.1, 4.2, Prop. 4.3 | CONFIRMED (4.1 uses P's Lemma CF of §7, re-derived here) |
| Q2 Prop. 5.1, Cors 5.2, 5.3 | CONFIRMED |
| Q2 Theorem 6.1, Cors 6.2, 6.3, Prop. 6.4 | CONFIRMED; "always" means "for every consistent network" (reaction side) / "every conservative network" (species side) |
| Q2 Lemmas 7.1, 7.2 | CONFIRMED (7.2 needs nonzero sources, e.g. (H1)) |
| Q2 table §7 (counts) | three rows re-computed with an exact CF: identical |
| E Props E.1, E.2, E.3 | CONFIRMED |
| E refutations E-a3…a8, E-b1…b4, E-c1…c4, E-d1…d4 | CONFIRMED (smallest counterexamples re-checked exactly) |
| E §3(a) summary sentence on "MIPS refine preclusters" | CORRECTION: true on the species side only |
| E2 Proposition E2.6 | CONFIRMED |

No claim in scope is refuted. Two corrections affect statements (Theorem B, E's §3(a) sentence), one affects a proof
sentence (Corollary P); the others are attributions, hypotheses to make explicit and one reading of a source.

## 1. Agent P (`petri/notes.md`)

### 1.1 Proposition 1.1 — CONFIRMED
Re-derived: arc Xᵢ → Rⱼ with weight (Ye)ᵢⱼ becomes kⱼ → ρᵢ with the same weight, and Rⱼ → Xᵢ with weight (Yp)ᵢⱼ
becomes ρᵢ → kⱼ; this is (T, S, F⁻¹) with weights, i.e. Desel–Esparza Def. 6.19 (p. 130, "dual net"; footnote: "In
some texts N^d is called the reverse-dual") and Teruel–Silva §2 ("reverse-dual N^rd with incidence matrix Cᵀ").
Checked by building the arc sets and reversing them (`check_P.py`, claim "1.1 reverse-dual") on 152 231 networks
(three exhaustive families, zero complex allowed in one).

### 1.2 Proposition 1.2 (dictionary) — CONFIRMED
Every row with a proof was re-derived. Siphons of N = T-traps of N*, traps = T-siphons, clusters swapped, FC ↔ FC,
"EC + balanced joins" ↔ itself (Ye[S_c, T_c] = a_c·J is transposition-invariant), connectivity, (H1) ↔ (H2),
(H3) ↔ distinct columns, rank, supports of extreme conservation laws = internal cycles of N*, twin pairs = doubling
pairs of N*, and (under (H1)) ray-MIPS(N) = full-coupling classes of N*: all checked mechanically on the same
152 231 networks, no exception. The "EC alone is not self-dual" example 2X₁+X₂ ⇄ X₃ was checked. Desel–Esparza
Prop. 6.20(2), (6) and Thm 6.21 say what is cited (p. 130–131 of the book PDF).

### 1.3 Propositions 2.1, 2.2, Corollary 2.3, Proposition 2.4 — CONFIRMED
Proofs re-derived. Prop. 2.2 checked on both sides of all family networks satisfying its hypotheses. Corollary 2.3
inherits the hypotheses of 2.2: it needs every reaction to have a nonzero source and every species to be consumed
(both hold under (H1)+(H2)); state them. Prop. 2.4 is definitional.

### 1.4 Remark 2.5 (B₂, B₃) — CONFIRMED
Exact: P(B₂) is one block; P(B₃) = {X→Y_a, Y_a→Z (a = 1,2,3)}, {Z→X}; also P(B₄) has two blocks. Definition 19 (one
step), the thesis Algorithm 4.4 (iterated) and the closure Ḡ agree on all of them, and the class-span closure is one
block for B₂, B₃, B₄. B₂, B₃, B₄ are WF-FC with (H1)–(H3). Independently of these examples, Def. 19 = thesis
Alg. 4.4 = Ḡ on all 22 829 consistent networks of three exhaustive families and on 699 random consistent networks with
a doubling pair (`check_misc.py`), consistent with Lemma 5 of the manuscript. The remark's conclusion (Definition 19
depends on how a choice is split) is correct and is made rigorous by Q2's Prop. 3.2 (§2.4 below).

### 1.5 §3 (classical theorems) — citations checked
Desel–Esparza (`petri/certificates/sources/fcbook/`): Def. 4.4 (clusters), Thm 4.30 (minimal siphons: at most one
place per cluster), Def. 5.1, Prop. 5.2(1)(3), Prop. 5.4, Prop. 5.7, Thm 5.17, Thm 6.14 (Rank Theorem, nets only),
Prop. 6.15, Thm 6.17 (live and bounded iff well-formed and every proper siphon marked), Def. 6.19, Prop. 6.20,
Thm 6.21, Exercise 6.6: all say what is cited. **One imprecision**: "Prop. 6.1(1) / Commoner: a structurally live FC
net is live from the marking with one token per place" — Prop. 6.1(1) says "structurally live iff every proper
siphon contains a proper trap"; the stated consequence follows with Commoner's theorem, and Theorem D(h) is anyway
proved from Thm 6.17, which is the right citation.
Teruel–Silva 1996 (`certificates/sources/`): Def. 6 item 3 (EC), Thm 7 (CF nets), Thm 12, Cor. 13 (live EC systems:
potentially reachable markings are live), Thm 15 (liveness monotonicity, bounded EC systems), Thm 28 (rank theorem
for EC nets: well-formed ⇔ well-structured and rank C = |𝓔| − 1), Cor. 32 (duality for balanced-joins EC nets), and
the remarks on siphons (p. 272: "no satisfactory generalisation of siphons and traps is available"; p. 296: siphon
material and Commoner's theorem not generalised). Teruel–Silva assume pure nets ("without loss of generality, nets are
assumed to be pure", §2): P's purity hypothesis in Theorem C is therefore the right one.
Hack 1972 (`literature/sources/hack.txt`): Well-Formedness Theorem (Thm 9), (a) ⇔ (c) "the reverse-dual is
Well-Formed"; remark "no interpretation yet for the dual or the reverse-dual notion of a trap". Correct.
Jančar–Best–Devillers–Ošťádal arXiv:2602.02352v2: Thm 3.23 (a) ⇔ (d) is the duality; §4 names T_Y• ⊆ •T_Y "minimal
T-traps". Correct (and the brief's "Thm 3.27" is indeed Prop. 3.27, a different statement).

### 1.6 Theorem A (transversality) — CONFIRMED
Proof re-derived (R \ {i} is a siphon; U \ {t} is a T-trap). Checked on both sides of all 152 231 family networks
and of all EC family networks below. Attribution: for free-choice nets this is Desel–Esparza Thm 4.30(a) with the
same three-line proof (P says so, "the first half of Thm 4.30"); the general weighted statement (places with equal
postsets) is the same argument. Present it as a mild generalisation with that credit.

### 1.7 Theorem B — CONFIRMED WITH CORRECTION
The proof is right provided a "twin pair" is a pair of species with **identical rows of Ye, zero rows included** (the
manuscript's definition and `crn.py`'s): the doubling pairs of N* include the pairs of species of N that are never
consumed (both reactions of N* have the zero source). P's §0 defines twin classes with *nonzero* rows; with that
definition the "no hypothesis" statement fails in a degenerate case: X₁ ⇄ X₂ with two further species X₃, X₄ that
occur in no reaction has P(N*) = {X₁}, {X₂}, {X₃, X₄}, minimal siphons {X₃}, {X₄}, {X₁, X₂}; the left side is false,
the right side (no nonzero-row twins, every species in a minimal siphon) is true (`check_misc.py`, claim 5).
**Corrected statement:** "Every precluster of N* lies in a minimal siphon of N iff N has no two species with equal
rows of Ye (zero rows included) and every species lies in a minimal siphon." The discrepancy needs two species
absent from every reaction, which (H3) excludes; under (H3), or if every species occurs, both readings agree. With the
manuscript's definition, Theorem B holds on both sides of all 152 231 family networks (no hypothesis) and on all EC
family networks.

### 1.8 Theorem C — CONFIRMED (pure WF-EC networks)
Every step re-derived: Gordan's alternative on Γ[R,·]; the integer vector σ ≥ 0 with Γ[R,·]σ < 0 (add a multiple of
a positive flux); the marking M₀′ ≥ M₀ with M = M₀′ + λΓσ ≥ 0, M[R] = 0; Teruel–Silva Thm 15 (needs boundedness: (H1))
gives liveness of M₀′, Cor. 13 gives liveness of M (M ∈ PR(N, M₀′)); an empty siphon stays empty and R• ≠ ∅ by
Lemma SC; contradiction. The consequence "minimal siphons = supports of extreme conservation laws" is correct (it
holds in every network without critical siphons). Checks: 144 random pure WF-EC networks (88 ordinary, 56 weighted,
seeded) and all 198 pure WF networks of the four re-enumerated EC families of §1.15; no exception. Impure WF-EC networks
(156 random, 62 of them ordinary): no critical siphon either, consistent with P's open question and E's 2 116.
**Attribution (L2, checked):** the ordinary case is classical — Desel–Esparza Prop. 5.4 + Prop. 5.2(3) for FC nets,
and for every ordinary net (impure included) with positive S- and T-invariants and Rank = #clusters − 1 it follows in
two lines from Desel–Esparza Thms 10.16–10.17 (Desel 1994): steps (ii)–(vi) of the proof of 10.17 show that any
marking marking every S-component is live, so a siphon containing no S-component would give a live marking with an
empty siphon. What is new in Theorem C is the weighted EC case (Teruel–Silva give well-formedness but no siphon
theory).

### 1.9 Lemma SC and Theorem C′ — CONFIRMED
Lemma SC: the terminal-strong-component argument is correct (inner sums ≥ 0, total 0, so no arc enters K).
Theorem C′: I re-derived the limit argument: the weighted balance of R_U = ∪_{c∈C_U} S_c gives total outside
production into R_U = ε Σ(|T_c| − 1)u_c z·y_c ≤ Cε; outside clusters producing into R_U have u = O(ε); shortest
paths to R_U through outside clusters propagate u = O(ε) (constants independent of ε); the limit flux is ≥ 0, nonzero
(Σ_{c∈C_U} v*_{t_c} = 1) and supported in U. The proof uses Theorem T(1) of §7 (a positive kernel vector of A_k for
every k), which is V2's scope; I read its proof (Lemma CF: support closed under successors, strong connectivity;
openness/closedness of the set of good k, rank A_k ≤ s = m − 1) and found no gap. Checks: 300 random WF-EC networks
(150 ordinary, 150 weighted, unbalanced joins included) and the 1 070 WF networks of §1.15: N* never has a critical siphon and the minimal T-traps are
the internal cycles. The consequences stated after C′ (bijective siphon/cycle duality for pure WF-EC; failure for D,
whose only internal cycle of D* is {X₁, X₂, X₃}) were checked. The remark "if N has balanced joins, C′ for N* gives C
for N without purity" is correct (N* is then WF-EC: checked on the 256 balanced-joins networks among the 300).

### 1.10 Theorem D (WF-FC networks) — CONFIRMED
(a)–(h) re-derived; each checked on R, T, W, W*, F, B₂, B₃, B₄, K₃, Y5 and on 150 random WF-FC networks
(`check_P.py`), including the profile rule (e) on both sides and (h) by an exhaustive token game from (1,…,1)
(liveness = every bottom strongly connected component of the reachability graph fires every transition). The only
citation remark is the one on Prop. 6.1 in §1.5.

### 1.11 K8 — CONFIRMED (checked exactly), with a caveat to add
N₁/N₂: both (H1)–(H3), same twin classes {X₁,X₄},{X₂},{X₃}, same species parts of clusters, same minimal siphons
{X₁,X₃},{X₂},{X₃,X₄}; P(N₁*) = {X₁,X₃,X₄},{X₂}, P(N₂*) one block. X₂ is a zero row of Γ(N₁) (P says so). E's pair
(E₁/E₂, `computation/notes.md` §6) is also confirmed; E₁ also has a zero row (X₁). M₁/M₂: both (H1)–(H3), no zero
row, same doubling classes {R₃,R₄} and minimal T-traps {R₄},{R₁,R₃},{R₂,R₃}; P(M₁) one block, P(M₂) = {R₁},{R₂},{R₃,R₄}.
**Caveat to add:** M₂ has two reactions with the same reaction vector (2X₃ → X₂+X₃ and X₁+X₃ → X₁+X₂), so M₂* violates
(H3); the reaction-side pair satisfies the hypotheses on the N side only (P's `check_k8.py` prints this, the notes do
not say it).

### 1.12 Theorem E — CONFIRMED
(i) R, T, W, W*, F are WF-FC with (H1)–(H3) and WF-FC duals with (H3). I recomputed every formulation of Table 1
(ray MIPS) on the five networks; the failure lists of the notes are reproduced exactly (A1, A4, B1, B4, Q1, Q6 fail on
all five; A2 on W*, F; A3 on R, T, W; B2 on W, F; B3 on R, T, W*; C on R, T, W*, F; Q2/Q4/Q5 on R, T, W (species side)
and R, T, W* (reaction side); Q3 on W*, F (species) and W, F (reaction); R23 on R, T, W*).
(ii) Proof re-derived (the CF case: every nonzero semiflow is positive, so ker Γ is a line); checked on the ten WF-FC
examples and through Theorem D(g) on 150 random WF-FC networks.

### 1.13 Table "manuscript's counterexamples in Petri terms" — CONFIRMED (columns in scope)
EC/balanced/ordinary/WF columns reproduced for R, T, W, W*, F, D, U, U*, N₃, N₃*. "crit" is the number of **all**
critical siphons (not only minimal ones): U 3 (2 minimal), U* 1, N₃ 1, N₃* 3 (2 minimal), D and D* 2; say so in the
caption. The "abs. cons." column depends on §7 (V2).

### 1.14 Propositions H, G — CONFIRMED
H: proof re-derived (s in a twin pair or adjoined through one; w_{s′} = μw_s lies in the same span). Checked on both
sides of every family network with (H1) (species side) and (H2) (reaction side): 152 231 networks.
G: (iii) ⇒ (i), (ii) re-derived; checked on all (H1) networks of the families. D refutes (i)–(iii) (M(D) one block,
minimal siphons {X₁,X₂},{X₁,X₃}); N₄.₁₂ refutes (iv) (law X₁+X₂, only minimal siphon {X₁}). Both exact.

### 1.15 Proposition F — CONFIRMED (no zero complex)
Proof re-derived. Checked on WR EC examples and on all WR members of the re-enumerated EC families; the zero-complex
restriction is needed (0 → X₂ → X₁ → 0 has no siphon at all). Credit: Mairesse–Nguyen arXiv:0905.3158, Lemma 3.2
(WR free-choice: nonempty complexes disjoint, = cluster places) and Cor. 3.2 (deficiency 0), checked in the source.
The remark "a WR EC network is automatically pure and conservative, well-formed iff one linkage class" is correct.

Re-enumeration of P's EC families with my own enumerator (`check_P_families.py`, definition copied from P's
docstring, independent code): (n, m, R, ordinary) = (3, 3, 6, 1): EC 8, H3 7, WF 7, WF-FC 7, WR 6; (2, 3, 5, 0): EC 75,
H3 75, WF 71, WF-FC 1, WR 6; (4, 2, 5, 1): EC 66, H3 55, WF 34, WF-FC 34, WR 15; (3, 2, 5, 0): EC 1 143, H3 1 141,
WF 958, WF-FC 6, WR 28 — all identical to P's logs. Theorems A and B (both sides, 1 292 networks), C (198 pure WF
networks), C′ (1 070 WF networks), D (48 WF-FC networks) and Prop. F (55 WR networks) pass on every member.

### 1.16 Corollary P — CONFIRMED WITH CORRECTION
The conclusion is correct (N* is conservative because N is consistent; no critical siphon by C′; Angeli–De
Leenheer–Sontag Thm 2, p. 607, applies to mass action, which their §2 includes explicitly). One proof sentence is
wrong as a justification: "Under (H1) for N every reaction of N* has a nonzero source". Reaction i of N* has a nonzero
source iff species i of N is consumed, which (H1) does not give (X → Y is conservative and Y is not consumed); it is
(H2) for N (part of WF-EC) that gives it — or simply: N* is conservative, so a zero-source reaction of N* would be
0 → 0. Replace "(H1)" by "(H2)".

### 1.17 Other statements checked in §§8–9
Agent G's example quoted in §9 (2A ⇄ 2B, 2A ⇄ A+B satisfies the ADLS criterion, its dual has the critical minimal
siphon {k₄}; neither EC): confirmed exactly. "The GAC holds for every WR EC network by classical results": correct
(δ = 0, conservative, no critical siphon ⇒ persistent by ADLS ⇒ global stability of complex-balanced equilibria).

## 2. Agent Q2 (`round2/preclusters/notes.md`)

### 2.1 Definitions used
P is Definition 19 of arXiv:2603.08767 v2 (one step: adjunction only through doubling pairs). The thesis
(Iosif 2019, Def. 4.3, Alg. 4.4) iterates through every edge of the current graph; the manuscript's Lemma 5 proves
both give the same partition, and my implementation of each agrees on all networks tested (§1.4). P_cs and π are
Q2's own constructions (whole-class spans) and are defined explicitly; CF = ∧_{V∈C°} κ_V. All are defined only under
(H2) (C° ≠ ∅); on the species side they need (H1) for N.

### 2.2 Lemma 1.2 — CONFIRMED, with the reading of Conradi–Flockerzi made explicit
(a) and (c) re-derived: D of the thesis (Def. 4.6, D = {λ ∈ Λ(E) : Eλ/k constant on doubling sets}) corresponds to
ψ ∈ ker A_k ∩ ℝᵐ_{>0} via Eλ = k∘Qψ, and the ratio identities over a nonempty relatively open subset of ker A_k are
the proportionality classes; κ(k) = κ_{k∘Qψ}. (b): Conradi–Flockerzi 2012 (PDF p. 118–121) define the clusters by
"n_jν/n_jλ ≡ n_{s_k}ν/n_{s_k}λ on Λ²(E)", where Λ²(E) = Λ(E) × Λ(E) (eq. above (3.4)), inside Algorithm 1 which starts
from "κ = ln(Eν/Eλ) ∈ ℝʳ ∩ im(𝒴ᵀ)", and Fact 2(a) derives the clusters from the doubling relations only (Example 1
uses only the doubling vectors). The domain of the identity is therefore not stated precisely: read over all of
Λ²(E) it would give the full-coupling classes; read over the pairs satisfying the doubling relations (Q2's reading,
which is the thesis' Def. 4.6) it gives CF; read over the pairs with κ ∈ im(𝒴ᵀ) (the full Fredholm constraint,
log-ratios of two steady states for the same rates) it can only be coarser, and coincides with CF when the source
complexes are linearly independent. **Correction for the paper:** state Lemma 1.2(b) as "with the clusters of
Conradi–Flockerzi read as in Fact 2(a) (identities on the pairs satisfying the doubling relations), equivalently the
thesis' Def. 4.6".

### 2.3 Propositions 2.1, 2.2 and the witnesses — CONFIRMED
Proofs re-derived (the class-wise ratio argument for π ⪯ κ_V is correct; flats of the vector configuration (n_t) are
complements of supports of flux vectors). Checked: P ⪯ P_cs ⪯ π ⪯ κ_V (two random V) and π = P_cs ∨ M* on every
consistent network of three exhaustive families (22 829 networks). Witnesses, exact (`check_Q2.py`,
`check_Q2_sym.py`):
- B₃: P two blocks, P_cs one block. R: P_cs singletons, π = M* one block.
- π ≺ CF: B₂dup (π = {R1–R4},{R5,R6}), K₃ (π = the 3 source classes; K₃ is WF-FC, WR, δ = 0), and the two-species
  witness {2X₂→X₁+X₂, 2X₂→2X₁, X₁+X₂→2X₂, X₁+X₂→2X₁}. CF = one block is certified for every V ∈ C° by an
  (m−1)-minor of A_{E(λ)} with sign-definite coefficients: K₃: Σ of all monomials of degree 2 except l₁², l₂², l₃²
  (positive); B₂dup: −(l₁+l₂)(l₃+l₄)(l₁+l₂+l₃+l₄); two-species witness: l₁ + 2l₂. It is among the smallest in the
  two-species families (no witness with r ≤ 3; `smallest_pi_cf.py`).
- XY4: A_V = (l₂ − l₃)·[[1, −1], [−1, 1]] in the ray parametrisation, so κ_V is one block off the hypersurface
  l₂ = l₃ and A_V = 0 on it (e.g. V = (1,1,2,2)); π = two blocks ⪯ CF ⪯ κ_{(1,1,2,2)} = two blocks, so CF = π = two
  blocks ≺ κ_V (generic). Z₄: CF = π = κ_V(generic) = {R1},{R2},{R3,R4}, κ_{(1,1,2,2)} = {R1,R2},{R3,R4}.
- XY4, Z₄, the two-species witness, B₂dup and K₃ satisfy (H1)–(H3).

### 2.4 Theorem 3.1, Propositions 3.2, 3.3 — CONFIRMED
Theorem 3.1 re-derived: for (O2) A′_{V′} = A_V; for (O3) the new row W forces ω′_{c′} = ω′_A and the two new columns add
up to the old column c, so ker A′_{V′} ≅ ker A_V with c′ and A merged; V′ ↦ V is onto C°(N) (O2) and bijective (O3).
The statement is pointwise in V, so I checked κ′_{V′} = φ⁻¹(κ_V) exactly at random V′ for every duplication of every
reaction and every splitting of every nonempty subset of every source class: 180 616 checks on all consistent networks
of the family n = 2, molecularity ≤ 2, r ≤ 4 (zero complex allowed) and 4 674 on 400 networks of n = 3, molecularity
≤ 2, r ≤ 3 (`check_Q2.log`), no exception.
Prop. 3.2: N₄.₁₂ vs N₄.₁₂ with R1 duplicated (P: two blocks → one), B₃ vs B₃ split (two → one), B₂ vs B₂dup (P_cs,
π: one → two), B₂dup vs its split (two → one): all exact; CF is one block for all of them (certificates). Prop. 3.3:
proof re-derived; the mechanism (after duplicating t, s fully coupled to t is adjoined to {t, t′}) checked on 9 447
fully coupled pairs. The remark "loops change P but not π or CF" was checked on R with a loop at X₁.
Note: Prop. 3.3 shows that every duplication-invariant extension of P contains P ∨ M*; it does not show that CF is the
finest invariant extension. The words "canonical" and "forced" in the summary should be read accordingly.

### 2.5 Theorems 4.1, 4.2, Proposition 4.3 — CONFIRMED
4.1: uses P's Lemma SC (§1.9) and Lemma CF (§7 of P; re-derived here: a nonnegative kernel vector of A_k has support
closed under successors, hence is positive, and spans the kernel). Checked at random V on every connected EC network
with (H1), (H2) in the families (438) and on the EC examples. 4.2: Γ = YI, A_V = Y L_V, ker A_V = ker L_V (δ = 0),
Laplacian kernel of a WR graph; checked on 1 762 random WR networks with δ = 0 (ker A_V = span{1_L}, κ_V = linkage
classes). 4.3: re-derived; it is what makes Theorem 6.1 work.

### 2.6 Proposition 5.1, Corollaries 5.2, 5.3 — CONFIRMED
Re-derived. 5.1 is the statement that a kernel with a positive vector has dimension ≤ number of proportionality
classes, with equality iff it has a basis with disjoint supports (PMDSC Condition 3.1 restricted to the source
complexes: PMDSC's Σ has zero columns for non-source complexes, so this restriction is exact). Checked dim ker A_V ≤
#κ_V on all consistent family networks; and "isolation with respect to π ⇒ dim ker A_V = #π = #κ_V" (the proof of
5.3 with π, which is finer than CF). Credit: "isolation ⇒ toric" is Iosif's thesis Thm 4.17 (rate-dependent form),
which Q2 cites; the rate-free form with P̂ and the equivalence (i) ⇔ (ii) of 5.2 are new as statements.

### 2.7 Theorem 6.1, Corollaries 6.2, 6.3, Proposition 6.4 — CONFIRMED
Re-derived (M* ⪯ π by 2.2; equal numbers of blocks plus refinement give equality; 4.3 gives the converse).
**On "MIPS refine P̂ always"** (attention point of the brief): it is proved for every *consistent* network on the
reaction side (M(N*) ⪯ P̂(N), Theorem 6.1, needs (H2)) and for every *conservative* network on the species side
(M(N) ⪯ P̂(N*), Cor. 6.2, needs (H1)); no further hypothesis ((H3), weak reversibility, deficiency) is used. Since P̂
is only defined under these hypotheses, "always" is accurate once this is said; the summary table and Remark 6.5(1)
should say "for every consistent (resp. conservative) network". The refinement is built in (π = P_cs ∨ M*), as Q2
acknowledges.
Checks: chain M* ⪯ π ⪯ κ_V, the equivalences of 6.1 for π and κ_V, the species-side statements of 6.2 (with π(N*),
κ_V(N*)), the "product complexes are sources + (H3)" clause, 6.4(b), all on the three families; 6.3 on 2 650 random WR
networks (no loops). Table 1 for P̂ (§6 of the notes): the surviving counterexamples listed (U, F, W* species side;
W, N₃, F reaction side; XY4, Z₄ for R23; T now satisfies R23) were checked; N₃ is the smallest reaction-side
counterexample (r ≤ 2 forces a doubling pair to be fully coupled under (H2)).

### 2.8 Lemmas 7.1, 7.2 — CONFIRMED
7.1 re-derived (f₁ > 0 on C° = {Eλ : λ > 0}). 7.2 re-derived (convex-combination/maximum argument, including the
case where β(p) also produces p); it needs the maximal block to consume a species, i.e. nonzero sources ((H1)), as Q2
notes in passing. Both are only used as sound merge rules.

### 2.9 Table of §7 (counts) — three rows re-computed exactly, identical
For two species every network with (H1) has s ≤ 1, Γ = dγᵀ and A_V = d f(V)ᵀ, so CF is decided **exactly** (c ≈ c′
iff f vanishes on every extreme ray outside {c, c′} and f_c is sign-definite and nonzero on the rays), with no bounds.
Families as E's `enum_all.py` (re-enumerated independently): rows all (2,2,5), (2,2,6), (2,3,4) — networks 100, 116,
2 784; P ≠ P_cs 0, 0, 27; P_cs ≠ π 10, 10, 63; π ≠ CF (N/N*) 2/0, 2/0, 38/0; κ_gen ≠ CF 5 in (2,3,4); isolation N/N*
20/100, 20/116, 224/2 784; species-side Q2 100, 116, 2 784 — **all identical** to Q2's table (`check_Q2_counts.py`).
The other rows (n ≥ 3, WR families) were not re-run; E2 re-ran several with an independent κ_V.

## 3. Agent E (`computation/notes.md` §3) and E2 (Prop. E2.6)

- **Prop. E.1** — CONFIRMED (Cor. 10(1) of the manuscript plus "no doubling pair ⇒ P discrete"); checked on 2 014
  random WR networks.
- **Prop. E.2** — CONFIRMED (Cor. 9 and Prop. 3(3) of the manuscript; the rank test s − rank Γ_{[n]∖{i,j}} ≥ 1 is
  right); checked on 566 random WR networks with (H1), (H3).
- **Prop. E.3 (Galois connection)** — CONFIRMED. Proof re-derived; brute force over all partitions P, Q and all
  union-closed families F ∋ ∅ of a ground set of ≤ 3 elements (5 806 triples): the adjunction holds, and c_P, i_Q are
  inverse isomorphisms iff P = Q.
- **Refuted readings** (E-a3…a8, E-b1…b4, E-c1, c1′, c2, c3, c4, E-d1…d4): every smallest counterexample printed in
  `tests_final.log` was re-checked exactly with my code, with the stated hypotheses ((H1), (H2); (H3) holds except for
  E-c4's network, whose rows X₁, X₂ of Γ are equal — E only claims (H1), (H2) there). The Y5 illustration (P(Y5),
  M(Y5*), join with two blocks; WR, δ = 0, (H1)–(H3)) is exact. E-a5's network X₃ → X₂, X₁+X₂ → X₁+X₃ has a pure
  catalyst, so its dual has a loop k₂ → k₂; allowed by the manuscript's conventions, worth a word in a paper.
- **Correction (§3(a) summary sentence).** "'MIPS refine preclusters' holds exactly when the MIPS are trivial" is
  proved (E.2) and true **on the species side only**. On the reaction side it is false: W is WR, δ = 0, (H1)–(H3),
  M(W*) = {R1,R3},{R2,R4} is not discrete, yet M(W*) ⪯ P(W) (one block); 38 of 376 random WR networks with (H1), (H3)
  behave like this. The reaction-side criterion is the manuscript's Theorem 7(2).
- **Prop. E2.6** — CONFIRMED. The proof is two lines and correct; it needs only that each linkage class is consistent
  on its own (true for WR). Checked: κ_V at V = Σ_L V_L refines the linkage classes on every family network where each
  linkage class is consistent (10 713 networks) and on 2 650 random WR networks.

## 4. Attribution and citations

- Checked against the sources: every Desel–Esparza, Teruel–Silva, Hack, Jančar et al., Angeli–De Leenheer–Sontag,
  Mairesse–Nguyen and Conradi–Flockerzi citation in scope (§1.5, §1.15, §1.16, §2.2). Only imprecisions: Desel–Esparza
  Prop. 6.1(1) (§1.5) and the domain of the Conradi–Flockerzi identity (§2.2).
- To credit in the paper (L2's findings, verified in the book where on disk): Theorem A for FC nets = Desel–Esparza
  Thm 4.30(a), same proof; Theorem C for ordinary nets = Desel–Esparza Prop. 5.4 (FC) and the two-line consequence of
  Thms 10.16–10.17 (Desel 1994) for all ordinary nets with the rank condition (I read Thm 10.17's proof: steps
  (ii)–(vi) use only "every S-component is marked"); new: weighted EC. Theorem D is assembled from Desel–Esparza
  (P says so). Prop. F's δ = 0 part = Mairesse–Nguyen (credited). "Iosif's dual = weighted reverse-dual" is
  definitional; its ingredients are in Hack, Desel–Esparza Prop. 6.20, Teruel–Silva §2, Jančar et al. §4 (P credits
  them). Isolation ⇒ toric = Iosif's thesis Thm 4.17 (credited by Q2). I found no statement of Theorem B, Prop. H,
  Lemma 1.2 or Theorem 3.1 in the sources on disk (consistent with L and L2; I did not repeat the web searches).

## 5. Corrected statements to use in the paper

1. **Theorem B:** "Every precluster of N* lies in a minimal siphon of N iff N has no two distinct species with equal
   rows of Ye (zero rows included) and every species lies in a minimal siphon. Dually …". (Or keep P's twin definition
   and add "every species occurs in a reaction".)
2. **Corollary 2.3:** add "every reaction has a nonzero source and every species is consumed (e.g. (H1) and (H2))".
3. **Corollary P, proof:** "Under (H2) for N (part of WF-EC) every reaction of N* has a nonzero source" (not (H1)).
4. **K8:** add "M₂* does not satisfy (H3) (two reactions of M₂ have the same reaction vector); N₁ and E's E₁ have a
   pure catalyst (zero row of Γ)".
5. **Theorem C:** add the credit "ordinary case: Desel–Esparza Prop. 5.4 (FC), Thms 10.16–10.17 (all ordinary nets
   with the rank condition)"; the weighted statement is the new part. **Theorem A:** credit Desel–Esparza Thm 4.30(a).
6. **Table of §6 of P:** "crit = number of critical siphons, minimal or not".
7. **Q2 Lemma 1.2(b):** "with the Conradi–Flockerzi clusters read as identities on the pairs (ν, λ) satisfying the
   doubling relations (Fact 2(a); = thesis Def. 4.6)".
8. **Q2 Theorem 6.1 / Cor. 6.2 / Remark 6.5(1):** "MIPS(N*) ⪯ P̂(N) for every consistent network; MIPS(N) ⪯ P̂(N*) for
   every conservative network" instead of "always / for every network"; and "CF is an invariant repair containing
   P ∨ M*; minimality among invariant repairs is not proved".
9. **E §3(a):** "'MIPS refine preclusters' holds exactly when the MIPS are trivial — on the species side (Prop. E.2);
   on the reaction side the criterion is Theorem 7(2) of the manuscript (W is a WR δ = 0 network with non-discrete
   M(W*) ⪯ P(W))".
10. **Desel–Esparza citation for D(h):** Thm 6.17 (not Prop. 6.1(1)).

## 6. What I could not check, and limits

- §7 of P (Theorems K5, T, T*, the U example) is V2's; I only read the parts of T(1) and Lemma CF that Theorem C′ and
  Q2's Theorem 4.1 use.
- Theorem C without purity for weighted, unbalanced EC networks remains open (no counterexample in 156 random impure
  networks here; the ordinary impure case is covered by Desel–Esparza Ch. 10).
- Counts of P's large families (11 748 EC networks, 8 757 WF-EC) and Q2's rows with n ≥ 3 were not re-run in full; I
  re-enumerated four of P's six EC families (1 292 of the 11 748 networks; identical class counts), and three of
  Q2's rows with an exact CF.
- Whether CF is the finest duplication- and splitting-invariant extension of Definition 19, and an exact finite
  algorithm for CF when s ≥ 2, are open (as Q2 says).
- No web search was repeated for novelty; attribution relies on the sources on disk and on L2's table, whose two
  main points (Thm 4.30(a), Thms 10.16–10.17) I verified in the book.

## 7. Scripts and logs (all exact; environment printed at the top of each log)

| script | what | log | runtime |
|---|---|---|---|
| `scripts/v1lib.py` | independent library (no import of the authors' code) | – | – |
| `scripts/examples.py` | named networks re-typed from the manuscript and the notes | – | – |
| `scripts/check_P.py` | P: 1.1, 1.2, 2.2, Rem. 2.5, A, B, C, C′, D, E, G, H, F, K8, table; three exhaustive families (152 231 networks, both sides), 300 random WF-EC | `logs/check_P.log` (93 PASS, 0 FAIL) | 565 s (Windows Python 3.14) |
| `scripts/check_P_families.py` | re-enumeration of four of P's EC families + A, B, C, C′, D, F | `logs/check_P_families_*.log` | 4 s, 5 s, 184 s, 42 s |
| `scripts/check_Q2.py` | Q2: hierarchy witnesses, 3.1, 3.2, 3.3, 4.1, 4.2, 5.1, 5.3, 6.1–6.4, E2.6 | `logs/check_Q2.log` (72 PASS) | 269 s |
| `scripts/check_Q2_sym.py` | sympy one-block certificates (B₂, B₃, B₂dup, splits, K₃, two-species witness), XY4, Z₄ | `logs/check_Q2_sym.log` | 0.6 s (WSL, sympy 1.14) |
| `scripts/check_Q2_counts.py` | exact CF for two species; three rows of Q2's table | `logs/check_Q2_counts.log` | 5 s |
| `scripts/smallest_pi_cf.py` | smallest π ≠ CF in the two-species families | `logs/smallest_pi_cf.log` | 5 s |
| `scripts/check_E.py` | E.1, E.2, E.3, smallest counterexamples of E's readings | `logs/check_E.log` (17 PASS) | 4 s |
| `scripts/check_misc.py` | Def. 19 = thesis Alg. 4.4 = Ḡ; loops; G's example; E §3(a) sentence; Theorem B twins | `logs/check_misc.log` | 47 s |

Authors' scripts re-run for comparison only: `petri/scripts/check_dictionary_examples.py` (all PASS, 0.2 s),
`petri/scripts/check_k8.py` (all PASS), `round2/preclusters/scripts/check_examples.py` (259 PASS, 0 FAIL, WSL).
