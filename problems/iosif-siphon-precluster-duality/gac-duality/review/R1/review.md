# Review R1 (mathematics) of the rewritten manuscript of result 17

Reviewer R1: **Claude Opus 5.5 (effort xhigh)**, fresh context, 27–28 October 2026. Brief: `../../BRIEF-R.md`, scope
"R1: mathematics". Read: `PROTOCOL.md` §3–§5, `../../BRIEF-common.md`, `paper.tex` (1 711 lines, all of it) and
`paper.pdf` (38 pages, text extracted), `../../writer/CHANGES.md`, the three reports `../../verification/V1..V3/verification.md`,
and the notes of agents P, E, G, Q1, Q2, Q3, E2, D2, A2 wherever a statement had to be compared with its source.
Paths below are relative to the problem folder `problems/iosif-siphon-precluster-duality/` unless they start with
`scripts/` or `logs/` of this review folder. Nothing outside this folder was edited.

Line numbers are those of `paper.tex` as reviewed (1 711 lines). Item numbers are the printed numbers (one counter
shared by theorem-like items, tables, figures and numbered equations; computed by `scripts/numbering.py`,
`logs/numbering.log`, and checked against the numbers quoted in `CHANGES.md`).

## 0. Summary

| part | verdict |
|---|---|
| Certificates | `scripts/certify_duality.py`: 129/129 PASS on Windows (`py -3`, Python 3.14) and on WSL (Python 3.12); `scripts/certify_examples.py`: 142/142 PASS (WSL). `scripts/figures.py` (copy writing into this folder) runs; its log is identical to `certificates/figures.log`; the figure is captioned as an illustration only. |
| Own exact re-checks | `scripts/r1_exact.py` (sympy, own code, no import of the authors' code): 61/61 PASS (`logs/r1_exact.log`). NC3, NC, NM, NB, 3A ⇄ 2A+B (class T = 2), Theorem O, Theorem A (identity and margin), Theorem AF, the examples of §5–§6, and two counterexamples found in this review. |
| CHANGES.md, A-items (new proofs) | all five **OK** (re-derived by hand and in exact arithmetic) |
| CHANGES.md, B-items (merged) | every merged statement has the hypotheses and conclusions of its verified sources, none stronger (item 11 is weaker, as CHANGES says). Presentation: B16 cosmetic **CORRECTION**; B18, B20, B21 **GAP** (no proof or no sketch label); B14, B17, B22 pointer missing |
| CHANGES.md, C-items (dropped) | all four **OK**; the writer's correction of P's remark on B_k is right |
| Corrections of V1/V2/V3 | applied, except V1 correction 2 (hypotheses of the cluster/source correspondence: only half applied, **CORRECTION** R-4) |
| Numbered items (76) and abstract | no statement refuted. **ERROR** in 3 places (R-1, R-5, R-6: one claim of the abstract, one displayed inequality in the proof of Theorem 68, one percentage); **CORRECTION** in 16 places (R-2–R-4, R-7–R-16, R-19, R-21b, R-21c); **GAP** (missing proof, pointer or "sketch" label) in 15 places (R-17, R-18 ×2, R-20, R-21, R-22 ×10); everything else **OK** |

No theorem of the paper is false. The three ERRORs are local: a factor `c` lost in a displayed inequality (the
conclusion stands), an abstract sentence that drops the hypothesis "equal-conflict" (counterexample below), and a
percentage copied from agent Q1's notes that contradicts Q1's own table.

## 1. Certificates and scripts

| run | result | log |
|---|---|---|
| `py -3 problems/iosif-siphon-precluster-duality/scripts/certify_duality.py` (repository root, Windows) | 129 PASS, 0 FAIL, exit 0, 1.0 s | `logs/certify_duality_windows.log` |
| same under WSL (`python -B`) | 129 PASS, 0 FAIL, exit 0, 0.8 s | `logs/certify_duality_wsl.log` |
| `scripts/certify_examples.py` (WSL) | 142/142 PASS, exit 0 | `logs/certify_examples_wsl.log` |
| `scripts/figures_copy.py` (= `scripts/figures.py` with the output folder redirected here) | runs in 17 s; status 0, 43 719 steps, drift 2.71e-11, min b 3.278e-20, min c 7.001e-13, period 2.5126e6: identical to `certificates/figures.log` | `logs/figures_copy.log` |
| `scripts/r1_exact.py` (own, WSL, sympy 1.14) | 61 PASS, 0 FAIL, 5 s | `logs/r1_exact.log` |

The caption of Figure 55 (l. 1232) says "Numerical illustration of Theorem 54, not part of the proof": **OK**. (The
docstring of `scripts/figures.py` says "about 1-3 minutes"; it takes 16–17 s, as §9 says. Not a paper issue.)

**Coverage of the certificate.** Every example and finite claim of the paper that I listed is in
`certify_duality.py`, except: (i) the converse example after Proposition 49 (agent E: X₄ ⇄ X₂+X₃, 2X₄ ⇄ X₁+X₃,
"criterion holds for N* but not N", l. 1135), which I certified separately (`r1_exact.py`, claim "Sec. 5.1"; N has the
critical siphons {X₃,X₄}, {X₂,X₃,X₄}, N* has none, both satisfy (H1)–(H3)); (ii) the two "smallest" claims (l. 1115:
"smallest conservative example is 2X₂ ⇄ X₁"; l. 1350: X "is the smallest weakly reversible network with (H1)–(H3) and a
NONE dual"), which are enumeration results of agent E and of V3. See R-21.

## 2. Own exact re-checks (`scripts/r1_exact.py`)

Written from scratch (own parser, Γ, dual, extreme rays by minimal supports, siphons by brute force, clusters, exact
restricted Jacobians). Positive steady states are counted on the **full** polynomial system (field and class
equations, no hand reduction): Rabinowitsch variable to discard boundary solutions, lex Gröbner basis in shape
position, exact isolating intervals of the univariate polynomial, exact signs of the other coordinates by interval
refinement. Results:

- **NC3 (Prop. 44).** Rank 2; extreme laws B+D, A+C+D; flux rays (0,1,1,0), (1,2,0,1); (H1)–(H3) on NC3 and NC3*;
  connected; clusters {D},{C},{A,B}, q = 3; not PC; minimal siphons {B,D}, {A,C,D}; no critical siphon on either
  side. f_B = −f_D and f_A + f_C + f_D = 0 for symbolic rates. At k = (1/6,1/3,6,6): D = 36A²B, C = 18AB + 36A²B;
  class → 72A³ − 204A² + 122A − 19 = (2A − 1)(36A² − 84A + 19). The full system has univariate polynomial
  (2A−1)(36A²−84A+19)/72 and exactly three real solutions, all positive. det J|_S = −7/9 at (1/2,1/3,6,3); at
  A = (7 ∓ √30)/6, trace J|_S = −61/3 ± (4/3)√30 < 0 and det J|_S = 10/3 ∓ (4/9)√30 > 0, decided exactly in ℚ(√30).
- **NC.** (H1)–(H3) both sides, laws A+2D and B+C, x = (3,1/2,3,1/4) steady at k = (2/3,5/18,24,64), det J|_S = −88/9.
  The class B+C = A+2D = 7/2 contains **exactly three** positive steady states (univariate
  (A−3)(800A⁶ − … + 297381)/800): V2's "exactly two more" is confirmed; the paper's weaker "at least one further" is true.
- **NM.** (H1)–(H3) both sides, PC, q = 2, s = 2, critical siphons {X₁,X₃}, {X₂,X₃}; x₃ = 1/2, x₁x₂² = 1/4; signs of
  t³ − 2t² + 1/4 at −1, 0, 1, 2 are −, +, −, +; the full system has exactly two positive solutions.
- **NB.** x₃ = 1, x₁ = (5 ± √17)/4, x₂ = 1/x₁ in 2x₁+x₂+x₃ = 6, and no other positive steady state.
- **3A ⇄ 2A+B (Example 50).** Dual, field v̇ = v(κ₁u³v − κ₂); unit rates, u+v = 2: h − 1 = −(v−1)(v³ − 5v² + 7v − 1),
  discriminant −44, g(1/8) = −103/512 < 0 < g(1/2) = 11/8, g(1) ≠ 0, exactly two positive equilibria; h′ = (2−v)²(2−4v),
  h(1/8) = 3375/4096, decay rate 721/4096.
- **Theorem 54 (O).** Dual and field; ḃ − 2ċ identity; the parametrisation of equilibria solves all three equations;
  the two derivative numerators; fluxes (17,4,7,1) at (1,1,1/100), T = 201/100; J = [[−59,3658],[−25,79]], trace 20,
  det 86 789, discriminant −4·86 689; H positive definite (det 1 177 283 151 421); JᵀH + HJ = 3 471 560 I.
- **Theorem 68 (A).** Identity (69) exact. **(log V)˙ = c·[γa³ − (3δ + 2βγa³/α)c − 3βδc²/α]** (factor c); the bracket is
  ≥ γT³/8 − Eε ≥ γT³/16 on {b + c < ε} at 400 exact random points and rates; at unit rates, T = 1, ε = 1/128, the point
  (a,b,c) = (1 − ε/2 − 10⁻⁶, ε/2, 10⁻⁶) has (log V)˙ ≈ 9.9·10⁻⁷ < γT³/16 = 1/16: the displayed inequality of the paper is
  false as written (R-5).
- **Theorem 70 (AF)**, u = 3..7, R = 1..7: law, (H1)–(H3), rank, dual, the (log V)˙ identity and the decompositions used
  for R = 1, 2.
- Examples X, S, B, Q (rank one), C (drain set), E's 4-cycle (drain species K₃ in path order), the (L) example, Example
  67 (U) identities, N1* steady state.
- **Two counterexamples found in this review**: (a) N = A+B → C, B+D → C, C → A+D, C → 2B (connected, (H1)–(H3), not EC)
  has the minimal siphon {A,C,D}, which meets the Petri cluster {A,B,D} twice (R-1); (b) 0 ⇄ X₁, 0 ⇄ X₂ is EC with (H2)
  and has m = 3 sources but q = 4 clusters (R-4).

## 3. CHANGES.md, item by item

| item | paper | verdict | note |
|---|---|---|---|
| A1 NC3 | Prop. 44 | **OK** | Reduction, cubic, factorisation, roots, positivity and stability re-derived by hand and in ℚ(√30) (§2). The count is now exhaustive on the full system, not only on the reduction. |
| A2 NC | text after Prop. 44 | **OK** | "At least one further" follows from the degree argument; the argument should be cited (R-14). Exact count confirms V2's "exactly two more". |
| A3 NM | Table 45, l. 1057 | **OK** | Re-derived (x₃ = 1/2 from r₁ = r₂, x₁x₂² = 1/4, x₁ = 2 − x₂). |
| A4 NB | Table 45 | **OK** | Closed form checked in ℚ(√17). |
| A5 class T = 2 of D | Example 50 | **OK** | §2. |
| B6 | Thm 9 | **OK** | Same hypotheses and conclusions as archived Thm 7, Cors 8–10 (compared with `archive/paper-precise-forms.tex`). |
| B7 | Thm 11, Table 10 | **OK** | Each entry of Table 10 re-read against the partitions printed by the certificate (S3 claims). |
| B8 | Prop. 15 | **OK** | |
| B9 | Thm 16 | **OK** | |
| B10 | Thm 17 | **OK** | (1) needs (H2) (assumed in the theorem) and (H1); (3) is E2.6 as confirmed by V1. |
| B11 | Prop. 18 | **OK** | Weaker than Q2 Cor. 5.2 (converse dropped), not stronger. |
| B12 | Thm 19 | **OK** | V1's wording applied; Q2 Prop. 6.4 is the sentence at l. 691–693, faithfully. |
| B13 | Thm 24 | **OK** | Sketch label and pointer present. |
| B14 | Cor. 32 | **OK** (GAP-level: one-line proof missing, R-22) | Follows from Thm 31 with q = m and Lemma 2. |
| B15 | Thm 33 | **OK** | |
| B16 | Prop. 41 | **CORRECTION** (cosmetic, R-13) | N1* steady state written with K₁ inside its own expression. |
| B17 | Thm 47 | **OK** | Hypotheses equivalent (Yp = Γ + Ye). Manifold part is a condensed sketch; add pointer (R-22). |
| B18 | Prop. 48 | **GAP** (R-20) | Statement correct (V3 §1), but no proof and no pointer in the paper; a five-line proof is supplied. |
| B19 | Example 53 | **OK** | Each of (1)–(4) re-derived (dual reactions written out by hand; drain species K₃ resp. K₄ in path order). |
| B20 | Cor. 61 | **OK**, but **GAP** (R-21): no proof in the paper; a proof is supplied. |
| B21 | Prop. 66 | **GAP** (R-17) | Restricted to μ = 1 as V3 requires, but not labelled as a sketch and without pointer; §9 overstates the completeness of (CM), (L), (S). |
| B22 | Props 72, 74 | **OK** (Prop. 74: pointer missing, R-22) | |
| B23 | Thm 58 | **OK** | V3's wording of (4) applied. |
| C24 | parenthesis l. 907–908 | **OK** | B_k: k+2 complexes, one linkage class, s = k+1, δ = 0, weakly reversible, so complex balanced for every rate vector (Horn–Jackson/Feinberg); P's remark was wrong; certified ("S4 B2, B3, B4 …"). |
| C25 | Cor. P dropped | **OK** | Superseded by Cor. 32 (persistence of N and N* for WF-EC without purity). |
| C26 | Conjecture CR not stated | **OK** | As V2 prescribes; open problem 4 keeps the absolute-consistency part. |
| C27 | E's 2X₂ ⇄ X₁ endotacticity claim dropped | **OK** | Dropping is safe. (For the record: the dual K₂ → K₁, 2K₁ → 2K₂ is indeed not endotactic: for w = (1,2) both sources have w-height 2 and 2K₁ → 2K₂ has w·(−2,2) = 2 > 0.) |
| D (notation) | | **OK**, one inconsistency | "twin classes of N = doubling classes of N*" (l. 324, Table 6) conflicts with "doubling classes = equal **nonzero** source" (l. 322–323) when two species have zero rows of Ye (R-3). |

## 4. Corrections of the verifiers

- **V1, §5 (ten statements).** 1 (twins with zero rows, Thm 21): applied. 2 (hypotheses of the cluster/source
  correspondence): **only half applied** — l. 340–342 derive "clusters correspond to sources and q = m" from (H2)
  alone; P's Prop. 2.2 also needs nonzero sources, which (H1) gives; counterexample 0 ⇄ X₁, 0 ⇄ X₂ (R-4). 3 (Cor. P
  proof): N/A (Cor. P dropped). 4 (K8 caveats on M₂*, N₁, E₁): N/A (Prop. 25 uses M₃/M₄, which satisfy (H1)–(H3) on
  both sides). 5 (credits of Thms A, C): applied (l. 655, l. 807–812). 6 (table "crit"): N/A. 7 (reading of
  Conradi–Flockerzi): applied (Lemma 14, l. 541–543). 8 ("for every consistent / conservative network"; "not shown to be
  finest"): applied (Thm 19, l. 576–577). 9 (E §3(a) sentence): N/A (the sentence is not in the paper). 10 (Desel–Esparza
  Thm 6.17 for D(h)): N/A because D(h) is not reproduced, but §9 l. 1626 lists it as applied (R-16).
- **V2, §9.** 1: unchanged statements, OK. 2: cluster-kinetics remark (Remark 36): applied verbatim. 3: Table 45 row 4
  ("species dual to X₁+2X₂ → 3X₂") and NM's two critical siphons, "NM is PC": applied. 4: CR replaced by Prop. 44 and
  open problem 4: applied. 5: D for span{y_c − y_c′}: applied. 6: WF-EC algebraic, Petri well-formedness for pure nets:
  applied (Def. 5, l. 342–345). 7: credits: applied (R2's scope).
- **V3, §8.** 1: (H1) in Cor. 61: applied. 2: Thm 58(4) wording: applied. 3: V_e(e) = 0 in (A2) of Thm 65: applied; the
  isolation step is mentioned but misstated (R-10). 4: Prop. E restricted to μ = 1: applied in Prop. 66 (E); the
  attribution sentence l. 1449–1450 should say who proved what (R-15). 5: Thm 63 in V3's general form: applied, but the
  convexity of λ was lost (R-8). 6: Conjecture F refuted by Thm 68: applied. 7: existence part of Thm 47 credited to E2:
  applied.

## 5. Verdict on every numbered item

| no. | label | verdict | comment |
|---|---|---|---|
| 1 | tab:results | OK (math) | Statuses consistent with the proofs; R2 checks the rest. |
| 2 | lem:transpose | OK | |
| 3 | prop:revdual | OK | |
| 4 | figure (3A ⇄ 2A+B) | OK | Matrices, arcs and weights checked against the dual written out by hand. |
| 5 | def:classes | OK | Equivalent forms of EC and PC re-derived (PC: transitivity inside a cluster, b_c primitive ⇒ a_t integral). |
| 6 | tab:dictionary | CORRECTION (R-3) | Row "twin classes / doubling classes". |
| 7 | lem:SC | OK | Proof re-read: both sums vanish, so no arc enters the terminal component, at places or at transitions. |
| 8 | prop:mips | OK | |
| 9 | thm:criteria | OK | |
| 10 | tab:formulations | OK | Each entry checked against the certified partitions. |
| 11 | thm:archive | OK | |
| 12 | prop:split | OK | Ray computation for B_k, the split of B₃ and the duplication in N₄.₁₂ re-derived by hand. |
| 13 | def:Phat | OK | Meet = coarsest common refinement, consistent with P̂ ⪯ κ_V. |
| 14 | lem:CF | OK | |
| 15 | prop:hierarchy | OK | |
| 16 | thm:invariance | OK | |
| 17 | thm:Phat-structure | OK | |
| 18 | prop:toric | OK | Equality case and isolation argument re-derived (V = Σ_B V\|_B gives 1_B ∈ ker A_V). |
| 19 | thm:duality-redone | OK | The step diag(V)⁻¹ker Γ = Q ker A_V when doubling pairs are fully coupled re-derived. |
| 20 | thm:transversal | OK | |
| 21 | thm:B | OK | Dual half uses doubling pairs including the zero source, consistent with l. 322. |
| 22 | prop:H | OK | |
| 23 | prop:G | GAP (R-22) | Correct; no proof or pointer (a three-line proof is supplied). |
| 24 | thm:WFFC | OK | Sketch, labelled, pointer present. Profile rule re-derived (0/1 columns with disjoint supports). |
| 25 | prop:K8 | OK | Status "checked"; certified. |
| 26 | eq:Z | OK | |
| 27 | lem:CFprime | OK | |
| 28 | lem:M | OK | Classical (Perron–Frobenius); V2 re-derived. |
| 29 | lem:R | OK | Closed–open argument re-read. |
| 30 | prop:K5gen | OK | |
| 31 | thm:Cgen | OK | Gordan step, vanishing outside C_R and B𝟏 = 0 re-derived. |
| 32 | cor:C | GAP (R-22) | One-line proof missing. |
| 33 | thm:T | OK | (2): dimension (r − s) + m − 1 < r iff s ≥ m, re-derived. |
| 34 | cor:Tbal | GAP (R-22) | One-line proof missing. |
| 35 | thm:U | OK | Uniqueness argument (max/min of w, predecessors) and the characteristic-polynomial count re-derived. |
| 36 | Remark (cluster kinetics) | OK | V2's correction applied. |
| 37 | thm:GS | OK | |
| 38 | thm:Tstar | OK | Cooperativity, invariance of the simplex, zero set closed under successors re-derived. |
| 39 | thm:Ustar | OK | G* formula re-derived (with Γ_pt it is the full matrix; for d ≠ c it equals V2's (y′_t)_p). |
| 40 | prop:GSstar | OK | |
| 41 | prop:Tstar-cex | CORRECTION (R-13, cosmetic) | All examples re-derived by hand; N1* steady state should not contain K₁. |
| 42 | thm:PC | OK | (c): null vector u(a∘w) of Ā_{a∘w} re-derived. |
| 43 | cor:BH | GAP (R-22) | Proof missing ("exactly the duals of WF-EC networks" needs two lines). |
| 44 | prop:NC3 | OK | §2. |
| 45 | tab:failures | OK | (m,q,s) of every row and each example re-derived; NC paragraph: CORRECTION R-14. |
| 46 | figure (classes) | CORRECTION (R-19) | "Global stability of N* for s ≥ 3 … is open" holds only for unbalanced joins (with balanced joins N* is WF-EC and Thm 37 applies). |
| 47 | thm:WRdual | OK (pointer, R-22) | Birch step re-derived; manifold part is a condensed sketch. |
| 48 | prop:nonWR | GAP (R-20) | No proof. |
| 49 | prop:ADLS | OK | Dual of H re-derived; {U₄} is a critical siphon. Converse example certified here only (§1). |
| 50 | ex:M | OK | §2. |
| 51 | prop:drain | CORRECTION (R-11) | "producing" must mean net production for the intended statement. |
| 52 | thm:drainsets | OK | Estimate x_Z^β ≤ U^{θ−1}Σu_i re-derived. |
| 53 | ex:drain | OK | |
| 54 | thm:O | OK | §2; Butler–McGehee step behind "lies in the ω-limit set of no other point" is correct. |
| 55 | fig:orbitO | OK | Illustration only; reproduced. |
| 56 | figure (face types) | OK | |
| 57 | lem:endpoint | GAP (R-22) | Pointer missing. |
| 58 | thm:S1 | OK | |
| 59 | cor:S1N | OK (pointer, R-22) | Dictionary (g = b, c = a, load, inward ⇔ net product) re-derived. |
| 60 | thm:S1WR | OK | Its hypothesis implies (H2) for N (both orientations occur), so Thm 58 applies. |
| 61 | cor:Cpure | GAP (R-21) | Correct; no proof in the paper. |
| 62 | ex:S1 | CORRECTION (R-7); the "smallest" claim needs its source (R-21b) | Example Q: "unless κ₂ = κ₃" claims more than Q1 proved; and the corollary fails in rank two also for weakly reversible pure networks (Ex. 53(4)). |
| 63 | thm:W | CORRECTION (R-8) | λ in the simplex dropped from D2's statement; the proof needs Σλ_q = 1. |
| 64 | lem:A | GAP (R-22) | Pointer missing. |
| 65 | thm:P2 | CORRECTION (R-9, R-10) | Standing hypothesis of §6.3 missing; isolation sentence of the sketch misstated. |
| 66 | prop:none-cert | GAP (R-17) | |
| 67 | ex:U | CORRECTION (R-12) | Reaction order not stated; with the path order the labels a = K₂, b = K₃ are swapped. |
| 68 | thm:A | **ERROR** (R-5) | Displayed inequality misses the factor c; conclusion unaffected. |
| 69 | eq:A1 | OK | |
| 70 | thm:AF | GAP (R-18) | Sketch labelled, pointer missing (A2's notes have the complete proof). |
| 71 | prop:L | GAP (R-22) | Pointer missing. |
| 72 | prop:dimension | OK | |
| 73 | prop:prime | GAP (R-18) | Sketch labelled, pointer missing. |
| 74 | prop:dictionary-boundary | GAP (R-22) | Pointer missing. |
| 75 | ex:BG | OK | det = 375 − 648 = −273; G is admissible for the first network's sign pattern. |
| 76 | tab:concordance | (R2) | |

Abstract and introduction (mathematical content only): l. 47–49 **ERROR** (R-1); l. 62 CORRECTION (R-2); the rest of the
abstract and §1.2 agree with the theorems. §6.1 "Counts and credits" l. 1362 **ERROR** (R-6). §9 l. 1597–1599 and
l. 1635–1637 CORRECTION (R-17, R-21c). §8 open problems: OK (NC3 is absolutely consistent: for any k choose A, B > 0,
then D = k₄A²B/k₁ and C = (k₃AB + 2k₄A²B)/k₂).

## 6. Corrections by severity, with exact replacements

Every replacement below is a literal substring of the current `paper.tex` (old) and its replacement (new). Where a
proof is added, "insert after line N" means after the line `\end{…}` that closes the item.

### 6.1 ERROR

**R-1 (abstract, l. 47–48).** Without equal conflict, the seeds of preclusters are the doubling classes, not the
clusters, and a minimal siphon can meet a cluster twice: N = A+B → C, B+D → C, C → A+D, C → 2B (connected, (H1)–(H3),
not EC) has the minimal siphon {A,C,D} and the cluster {A,B,D} (`r1_exact.py`). Theorem 20 itself says "In an EC
network".
- old: `Through it, siphons of \(\Net\) are the T-traps of \(\Nd\), the seeds of preclusters are the
Petri clusters, and minimal siphons cross every cluster at most once, so preclusters and siphons are related by`
- new: `Through it, siphons of \(\Net\) are the T-traps of \(\Nd\); for equal-conflict networks the seeds of
preclusters are the Petri clusters and minimal siphons cross every cluster at most once, so preclusters and siphons are related by`

**R-5 (proof of Theorem 68, §6.3, l. 1466).** (log V)˙ = c·[bracket]; the bracket, not (log V)˙, is ≥ γT³/16 (A2's
notes §2.3 say "the bracket"). As written the inequality is false near c = 0 (§2). The last-entry argument only needs
(log V)˙ > 0, so the theorem stands.
- old: `so \((\log V)^\cdot\ge\gamma T^3/8-E\varepsilon\ge\gamma T^3/16>0\) on \(\{b+c<\varepsilon\}\) with`
- new: `so \((\log V)^\cdot=c\,[\gamma a^3-(3\delta+2\beta\gamma a^3/\alpha)c-3\beta\delta c^2/\alpha]\ge c\,(\gamma T^3/8-E\varepsilon)\ge c\,\gamma T^3/16>0\) on the interior points of \(\{b+c<\varepsilon\}\), with`

**R-6 (§6.1 "Counts and credits", l. 1362).** Q1's table (and V3's recount) give ALL fractions 4/12, 115/460,
381/3 544, 662/2 742, 4/12, 101/426, 245/829 and 1/1, i.e. 11–33 % for the families with more than one network; "17–35 %"
is in Q1's prose and contradicts Q1's table (381/3 544 = 10.75 %).
- old: `\(\Nd\) is persistent for every rate vector for 17--35\,\% and for none for the rest,`
- new: `\(\Nd\) is persistent for every rate vector for 11--33\,\% of the networks of each family with more than one network (from 381 of 3\,544 to 4 of 12; the family \((4,3,4)\) has a single network, which is ALL) and for none for the rest,`

### 6.2 CORRECTION (hypotheses, statements, proofs)

**R-2 (abstract, l. 62).** The periodic orbit exists for one rate vector, not for all.
- old: `class, and a periodic orbit).`
- new: `class, and, for some rate vector, a periodic orbit).`

**R-3 (§2.2, l. 324; Table 6, l. 362).** Doubling classes are defined with nonzero sources (l. 322–323), twins with zero
rows included; two species never consumed are twins but give two reactions of N* with the zero source, which form a
doubling *pair* but no doubling *class*.
- old (l. 324): `the twin classes of \(\Net\) are the doubling classes of \(\Nd\).`
- new: `the twin pairs of \(\Net\) are the doubling pairs of \(\Nd\) (two zero rows give two reactions of \(\Nd\) with the zero source), and the twin classes of nonzero rows are its doubling classes.`
- old (l. 362): `twin classes / doubling classes & doubling classes / twin classes & definitions \\`
- new: `twin pairs / doubling pairs & doubling pairs / twin pairs & definitions \\`

**R-4 (§2.2, l. 342; V1 correction 2).** "q = m" needs every reaction to have a nonzero source (P's Prop. 2.2). (H2)
alone does not give it: 0 ⇄ X₁, 0 ⇄ X₂ is EC and consistent with m = 3, q = 4. Under (H1) a zero source can only occur
in 0 → 0. All later uses (WF-EC, Thm 33, Cor. 32) have (H1), so no theorem changes.
- old: `vector), so for EC networks clusters correspond to sources and \(q=m\).`
- new: `vector), and under (H1) a reaction with the zero source is \(0\to0\) (\(z\tp y'=0\) for a positive law \(z\)); so for EC networks with (H1), (H2) and no reaction \(0\to0\) clusters correspond to sources and \(q=m\) (with (H2) alone, \(0\rightleftarrows X_1,\ 0\rightleftarrows X_2\) has \(m=3\), \(q=4\)).`

**R-7 (Example 62, l. 1356–1358).** Q1 proved non-persistence only for κ₂ ≠ κ₃ ("the hypersurface κ₂ = κ₃ was not
analysed"); "unless" suggests persistence on κ₂ = κ₃. And the corollary fails in rank two for weakly reversible pure
networks as well (Example 53(4), E2's refutation of Q1-6), so "without weak reversibility" must not read as a
restriction.
- old: `ends). In rank two the corollary for pure networks fails without weak reversibility:`
- new: `ends). In rank two the corollary for pure networks fails, already without weak reversibility:`
- old: `non-persistent unless \(\kappa_2=\kappa_3\).`
- new: `non-persistent for every rate vector with \(\kappa_2\neq\kappa_3\) (the hypersurface \(\kappa_2=\kappa_3\) was not analysed); \zcref{ex:drain}(4) shows that it also fails for weakly reversible pure networks.`

**R-8 (Theorem 63, l. 1377).** D2's statement has "λ in the simplex"; the proof uses Σ_q λ_q = 1 (μ − μ = 0).
- old: `\((y_l)_o=\sum_q\lambda_q(y_l)_q\) for every \(l\notin Z\).`
- new: `\((y_l)_o=\sum_q\lambda_q(y_l)_q\) for every \(l\notin Z\), with \(\lambda_q\ge0\) and \(\sum_q\lambda_q=1\).`

**R-9 (§6.3, l. 1387).** The section never states its standing hypothesis; D2's setting is (H1)–(H3) for N, which makes
the classes of N* compact polygons (needed by Lemma 64, Thm 65, Prop. 66).
- old: `Let \(s=2\). Factor`
- new: `Let \(s=2\) and let \(\Net\) satisfy (H1)--(H3), so that the positive classes of \(\Nd\) are compact convex polygons. Factor`

**R-10 (sketch of Theorem 65, l. 1406–1407).** The ω-limit set of an orbit along which V_e increases has V_e = sup, not
V_e = 0; V3's isolation argument goes through ∂P and {e}.
- old: `(A2) or by the normal coordinate of an IN edge, and is isolated: an interior total orbit near \(e\) would have its
\(\alpha\)- and \(\omega\)-limits in \(\{V_e=0\}\cap\partial P\), which contains no invariant set besides \(e\) (V3's addition).`
- new: `(A2) or by the normal coordinate of an IN edge, and is isolated: along an interior total orbit in a compact invariant set near \(e\), \(V_e\) is strictly increasing, so its \(\alpha\)- and \(\omega\)-limit sets lie in \(\partial P\); near \(e\) the only invariant subset of \(\partial P\) is \(\{e\}\) (boundary orbits are monotone along the edges), so both limit sets are \(\{e\}\), which contradicts the strict increase of \(V_e\) between them since \(V_e(e)=0\) (V3's addition).`

**R-11 (Proposition 51, l. 1159 and l. 1163).** In the paper "producing" means "appearing in the product" (siphons,
l. 232). With that reading (b) excludes a draining reaction such as 2Xᵢ → Xᵢ + Y (its source coefficient is d₀). The
proof and E's Prop. E.4 use net production.
- old (l. 1159): `consumes \(X_i\), and (b) every reaction producing \(X_i\) has a source with \(X_i\)-coefficient \(>d_0\).`
- new: `consumes \(X_i\), and (b) every reaction with a net production of \(X_i\) (\(\Gamma_{ij}>0\)) has a source with \(X_i\)-coefficient \(>d_0\).`
- old (l. 1163): `production terms of \(x_i\) are \(\le Cx_i^{d_0+1}\)`
- new: `the terms of \(\dot x_i\) with \(\Gamma_{ij}>0\) are \(\le Cx_i^{d_0+1}\)`

**R-12 (Example 67, l. 1442).** The printed dual corresponds to the order 3X₃ → 2X₁+X₂, X₂+2X₃ → 3X₃, 2X₁+X₂ → X₂+2X₃
(as in the certificate); in the path order K₂ and K₃ are exchanged and "a = K₂, b = K₃" is wrong.
- old: `For \(3X_3\to2X_1+X_2\to X_2+2X_3\to3X_3\), the dual`
- new: `For \(3X_3\to2X_1+X_2\to X_2+2X_3\to3X_3\), with the reactions numbered \(3X_3\to2X_1+X_2\), \(X_2+2X_3\to3X_3\), \(2X_1+X_2\to X_2+2X_3\), the dual`

**R-13 (Proposition 41(1), l. 965).** The middle coordinate is written with K₁; at K₁ = κ_B/κ_A it is 3κ_B²/(κ_Aκ_C)
(checked symbolically).
- old: `with the steady state \((\kappa_B/\kappa_A,3\kappa_AK_1^2/\kappa_C,1)\).`
- new: `with the steady state \((\kappa_B/\kappa_A,3\kappa_B^2/(\kappa_A\kappa_C),1)\).`

**R-14 (paragraph after Prop. 44, l. 1031–1032).** "The index sum is +1" is the degree argument of Thm 39(b), which
needs a compact class without boundary steady states (here from Thm 31).
- old: `since the index sum is \(+1\), the
class contains at least one further positive steady state (V2: exactly two more, by exact root isolation).`
- new: `the class is compact without boundary steady states (\zcref{thm:Cgen}), so by the degree argument of \zcref{thm:Ustar}(b) the indices of its steady states add up to \((-1)^s=+1\), and the
class contains at least one further positive steady state (V2: exactly two more, by exact root isolation).`

**R-15 (§6.3, l. 1449–1450).** D2's written proof of the two-coloured case had a gap that V3 repaired for μ = 1.
- old: `(D2 proved it for
two-coloured and zero edge flows with \(\mu=1\)).`
- new: `(D2 proved it for zero edge flows with \(\mu=1\), and for two-coloured edge flows with \(\mu=1\) with V3's repair).`

**R-16 (§9.2, l. 1622 and l. 1626).** A missing backslash prints "zcrefthm:B" in the PDF (checked in `paper.pdf`); and
the Thm 6.17 citation concerns P's D(h), which the paper does not reproduce.
- old: `transversality, zcref{thm:B}, the WF-FC structure`
- new: `transversality, \zcref{thm:B}, the WF-FC structure`
- old: `the citation of Desel--Esparza Thm~6.17;`
- new: `the citation of Desel--Esparza Thm~6.17 (for P's liveness statement, not reproduced here);`

**R-19 (caption of Figure 46, l. 1084).** With balanced joins N* is itself WF-EC, so Thm 37 gives its global stability
in every rank; open problem 2 correctly says "with unbalanced joins".
- old: `Global stability of \(\Nd\) for \(s\ge3\), and of general PC`
- new: `Global stability of \(\Nd\) for \(s\ge3\) with unbalanced joins, and of general PC`

**R-21b ("smallest" claims, l. 1115 and l. 1350).** Both are enumeration results, not certified by
`certify_duality.py`.
- old (l. 1115): `the smallest conservative example is \(2X_2\rightleftarrows X_1\),`
- new: `the smallest conservative example in agent E's exhaustive families is \(2X_2\rightleftarrows X_1\),`
- old (l. 1350): `it is the smallest weakly
reversible network with (H1)--(H3) and a NONE dual.`
- new: `it is the smallest weakly
reversible network with (H1)--(H3) and a NONE dual (V3: one species is impossible, \(X_1\rightleftarrows X_2\) is ALL, and V3's enumeration of the family \((2,2,4)\) finds \(X\) first).`

**R-21c (§9.1, l. 1598).** The converse example after Proposition 49 is not in `certify_duality.py`. Preferred fix
when applying: add one claim to `certify_duality.py` (the check is in `scripts/r1_exact.py` of this review, "Sec. 5.1")
and change "129" to "130" in the abstract (l. 65), §9.1 (l. 1597) and the log. Minimal fix if the script is not touched:
- old: `It certifies every example and finite claim of this paper:`
- new: `It certifies every example and finite claim of this paper except the converse example after \zcref{prop:ADLS} (checked separately in exact arithmetic by the reviewer R1) and the enumeration counts cited as checked by computation:`

### 6.3 GAP (proofs, sketch labels and pointers)

**R-17 (Prop. 66, insert after l. 1422; §9.2 l. 1635–1637).** The brief requires the extinction certificates to be
labelled as sketches with a pointer; V3 checked (CM), (L), (S) only at sketch level, while §9.2 says that the sketched
items "have complete proofs in the working notes".
- insert after l. 1422:
```
\begin{proof}[Proof (sketch)]
Proofs in agent D2's notes (\texttt{gac-duality/round3/dichotomy/notes.md}, Propositions E, H, CM, L, S, B). (E) for
\(\mu=1\), as repaired by V3: on the facet there is a point \(p^*\) with zero set \(Z\), zero edge flow and normal
coefficient \(L(p^*)<0\); it is an equilibrium of its class, and the stable manifold of the eigenvalue \(L(p^*)\) is
transverse to the invariant edge and enters the interior. (H): Cauchy--Binet for a sum of rank-one \(2\times2\) matrices;
a Hurwitz vertex attracts an open set of its class. (B): the two blocks decouple. (CM), (L), (S) are standard
centre-manifold, Lyapunov and sector arguments, written as sketches in D2's notes and checked by V3 at sketch level and on
the worked examples below.
\end{proof}
```
- old (l. 1635–1637): `Items marked ``sketch'' (\zcref{thm:WFFC,thm:P2,thm:AF,prop:prime} and the certificates of \zcref{prop:none-cert}) have
complete proofs in the working notes, re-derived by the verifiers; V3 checked the tropical certificates (CM), (L), (S) at
sketch level and on their worked examples only.`
- new: `Items marked ``sketch'' (\zcref{thm:WFFC,thm:P2,thm:AF,prop:prime,prop:none-cert}) have complete proofs in the
working notes, re-derived by the verifiers, except the tropical certificates (CM), (L), (S) of \zcref{prop:none-cert},
whose proofs in D2's notes are themselves sketches; V3 checked them at sketch level and on their worked examples only.`

**R-18 (pointers for Thm 70 and Prop. 73).**
- old (l. 1480): `\(T=\min(1,\delta/(2(u-2)\gamma))\). Certified for \(u=3,\dots,7\), \(R=1,\dots,7\).`
- new: `\(T=\min(1,\delta/(2(u-2)\gamma))\) (using \(T^u\le T\) for \(T\le1\)). Certified for \(u=3,\dots,7\), \(R=1,\dots,7\). Full proof in agent A2's notes (Theorem AF), re-derived by V3 (§6).`
- old (l. 1532): `direction not orthogonal to \(S\) is nonconstant on every class, which gives strong endotacticity.`
- new: `direction not orthogonal to \(S\) is nonconstant on every class, which gives strong endotacticity. Full proof in agent G's notes (Proposition~4.1), re-derived by V3 (§1).`

**R-20 (Prop. 48, insert after l. 1120).**
```
\begin{proof}
If \(\Net\) is weakly reversible, every product complex is a source complex, so equal rows of \(\Ye\) would give equal
rows of \(\Yp\) and of \(\Gamma\), excluded by (H3); the reversed network gives the same for \(\Yp\). So \(\Nd\) has \(n\)
distinct sources and \(n\) distinct products. If \(\Nd\) is weakly reversible, each product is a source, which defines
an injective map \(\pi\) with \((\Yp)_{i\cdot}=(\Ye)_{\pi(i)\cdot}\); conversely the complex graph \(i\mapsto\pi(i)\) of a
permutation is a disjoint union of cycles. It has \(n\) complexes and \(c(\pi)\) linkage classes, so the deficiency of
\(\Nd\) is \(n-c(\pi)-s\), and \(\sum_{i\in C}\Gamma_{i\cdot}=\sum_{i\in C}((\Ye)_{\pi(i)\cdot}-(\Ye)_{i\cdot})=0\) for every
cycle \(C\) of \(\pi\). The examples are certified (agent G, Prop.~2.2; V3 §1).
\end{proof}
```

**R-21 (Cor. 61, insert after l. 1345).**
```
\begin{proof}
By (H1) the vector \(a\) has entries of both signs. Let \(Z\) be a face of \(\Nd\) (\zcref{cor:S1N}). A species with
\(a_i\neq0\) is a net reactant of every reaction of \(Z\) or a net product of every reaction of \(Z\). A net reactant
occurs in every source of \(Z\), so its load is at least \(|Z|\ge1\); a net product occurs in no source of \(Z\) by
hypothesis, so its load is \(0\), and net products exist by (H1). So every species of minimal load with \(a_i\neq0\) is a
net product, every face is IN, and \zcref{cor:S1N} applies. If \(n=2\), every reaction consumes one species and produces
the other, the faces of one orientation all consume the same species, loads add over \(Z\), so the singleton faces
decide, and \(\{j\}\) is IN iff \(y_j\) contains more of the consumed than of the produced species.
\end{proof}
```

**R-22 (one-line proofs or pointers; insert after the closing line given).**
- Prop. 23 (after l. 689): `\begin{proof} Conservation-coupled species (\(w_i=\mu w_j\)) lie in the same supports of extreme laws, and by (H1) every species lies in one. If every siphon contains the support of a nonnegative law, a minimal siphon equals such a support, which is then the support of an extreme law, and conversely. The examples are certified. \end{proof}`
- Cor. 32 (after l. 805): `\begin{proof} \zcref{thm:Cgen} with \(q=m\) (\zcref{sec:setting}); the internal cycles of \(\Nd\) are the supports of the extreme laws of \(\Net\) (\zcref{lem:transpose}). \end{proof}`
- Cor. 34 (after l. 851): `\begin{proof} \(\Nd\) is EC with balanced joins, connected, of rank \(s\), with the same clusters and (H1), (H2) exchanged (\zcref{tab:dictionary}); apply \zcref{thm:T} to \(\Net\) and to \(\Nd\). \end{proof}`
- Cor. 43 (after l. 1003): `\begin{proof} If \(\Net\) is BH, two species of one block have equal rows of \(\Ye\), so \(\Nd\) is EC with clusters the blocks, connected, of rank \(s\), with (H1), (H2) exchanged and \(s=m(\Nd)-1\): \(\Nd\) is WF-EC, and conversely the dual of a WF-EC network is BH. Apply \zcref{thm:Tstar,thm:Ustar,cor:C,prop:GSstar} to the WF-EC network \(\Nd\) (or \zcref{thm:PC}). \end{proof}`
- Thm 47 (end of proof, l. 1102): after `critical-point equation is the steady-state equation.` add ` Details: agent A2's notes (Theorems E, E+), re-derived by V3 (§6).`
- Lemma 57 (after l. 1297): `\begin{proof} Agent Q1's notes, Lemmas 1.1--1.2 (along a class \(F\) is a polynomial in \(\varepsilon\)); re-derived by V3 (§3). \end{proof}`
- Cor. 59 (after l. 1326): `\begin{proof} \zcref{thm:S1}(1) for \(\mathcal M=\Nd\), which is conservative by (H2), with the dictionary above. \end{proof}`
- Lemma 64 (after l. 1396): `\begin{proof} Agent D2's notes, Lemma A (the siphon edge is invariant in both time directions for the polynomial field); re-derived by V3 (§5). \end{proof}`
- Prop. 71 (after l. 1486): `\begin{proof} Agent A2's notes, Proposition L (opposite active vectors on a siphon edge would give a nonnegative conservation law of \(\Net\) vanishing on a complex, hence on the whole linkage class); re-derived by V3 (§6). \end{proof}`
- Prop. 74 (after l. 1550): `\begin{proof} Agent G's notes, Propositions 6.1--6.3 ((2): the tangent space of the face is the kernel of the projection of \(S\) onto \(\mathbb{R}^\Sigma\)); re-derived by V3 (§1). \end{proof}`

## 7. Points checked and found correct (selection, beyond §2)

- Lemma 7, Lemma 27, Lemma 29, Prop. 30, Thm 31, Thm 33, Thm 35 (all three parts), Thm 37, Thm 38, Thm 39, Thm 42: each proof
  in the paper re-read line by line; no gap beyond those listed (they are condensed versions of Q3's and P's proofs
  confirmed by V2).
- Thm 9, Thm 19, Prop. 18, Lemma 14: re-derived (equality of diag(V)⁻¹ker Γ and Q ker A_V under full coupling of doubling
  pairs; the isolation ⇒ equality argument).
- Thm 54 and Thm 68: the boundary analysis (strictly inward field on a = 0 and b = 0, nonzero field at the two vertices,
  backward escape on c = 0) and the last-entry argument re-derived; only the factor c is wrong (R-5).
- Thm 58–60: endpoint expansion, the sign conventions (upper face Z ⊆ N₋, inward ⇔ σc_j < 0), the dual dictionary
  (g = b, c = a, load, inward ⇔ net product) and the S1-WR argument (affine φ_o − φ_v with one zero per linkage class).
- Every dual network written in §§2–7 was recomputed by hand from Ye and Yp (Fig. 4, Thm 54, Thm 70, Ex. 50, 53, 62, 67,
  Props 41, 49, the rank-two examples, Table 45): all correct (Ex. 67 up to the unstated order, R-12).
- The Table 1 statuses agree with the proofs (P. with sketch for 24, 65, 66, 70, 73; Ch. for 25 and the counts).

## 8. What I could not check

- The general tropical certificates (CM), (L), (S) of Prop. 66 (only their worked examples, as V3); the rank-two counts
  (4 873 / 13 817 / 171) were not re-run.
- The large enumeration counts of §§3–6 (82 377 networks, 18 861, the rank-one families, 63 717/63 722) were not
  re-run; only their internal arithmetic was checked (e.g. 2 988 + 1 885 = 4 873, 1 048 + 12 769 = 13 817). The two
  drain-species counts (141 over the 82 377 networks, "all with s ≤ 2", and 199 over E2's 18 861 networks with s = 2)
  refer to different sets of files; the sentence at l. 1196–1202 does not say so (consistency question for R2).
- Thm 47 (manifold part), Thm 16 (invariance under splitting) and Thm 24 were checked by reading the proofs in the notes
  and the verifiers' reports, not re-derived in full.
- Citations, novelty wording, house-rule wording (including the lines 66, 163, 1363, 1364, 1438, 1620 that use a word
  the house rules restrict) and the build are R2's scope.

## 9. Files

| file | content |
|---|---|
| `scripts/r1_exact.py` | own exact checks, 61 claims (sympy 1.14, WSL; 5 s) → `logs/r1_exact.log` |
| `scripts/figures_copy.py` | `scripts/figures.py` with the output folder redirected to this folder → `logs/figures_copy.log`, `figures/theoremO-orbit.pdf` |
| `scripts/numbering.py` | printed numbers of all numbered items → `logs/numbering.log` |
| `scripts/check_replacements.py` | every "old" text of §6 occurs exactly once in `paper.tex` (28/28) → `logs/check_replacements.log` |
| `logs/certify_duality_windows.log`, `logs/certify_duality_wsl.log`, `logs/certify_examples_wsl.log` | re-runs of the certificates |
| `PROGRESS.md` | progress notes |
