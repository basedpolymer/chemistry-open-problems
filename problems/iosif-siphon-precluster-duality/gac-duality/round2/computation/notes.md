# Agent E2: computation and red team, second round of the gac-duality work

**Agent E2: Claude Opus 5.5, effort xhigh. 27 October 2026.** Working notes for result 17
(`problems/iosif-siphon-precluster-duality/`), second round (`../../BRIEF-round2-common.md`, `../../BRIEF-E2.md`).
Nothing here quotes the author's e-mail. Status tags: **proved** (complete proof here), **checked** (exact computation
on the finite families listed, not a proof), **refuted** (explicit counterexample, checked exactly), **open**.
Paths are relative to `problems/iosif-siphon-precluster-duality/gac-duality/round2/computation/`; record files live
outside the repository in WSL `~/gacE-data/` (first round) and `~/gacE2-data/` (this round; SHA-256 in
`logs/data_sha256_m3s4.txt`).

## 0. Conventions

As in `../../computation/notes.md` §0: N has species X_1..X_n, reactions y_j → y'_j (j = 1..r), Ye, Yp, Γ = Yp − Ye,
s = rank Γ; the dual N* has species k_1..k_r and one reaction per species of N, k^{Ye[i,:]} → k^{Yp[i,:]}, so
Γ(N*) = Γᵀ. (H1) conservative, (H2) consistent, (H3) rows of Γ distinct; H = (H1)–(H3) for N; WR = weakly
reversible. Under (H1), (H2) for N, the dual is consistent and conservative, so its compatibility classes are compact
polytopes of dimension s. "Persistent" (Angeli–De Leenheer–Sontag): no trajectory starting in the open orthant has
an ω-limit point on the boundary. ALL / SOME / NONE: N* persistent for every / some but not every / no rate vector κ.
WF-EC (agent P): equal conflict (distinct sources have disjoint supports), connected, (H1), (H2), s = m − 1 with m
the number of source complexes.

Every script prints PASS/FAIL per claim, its environment and runtime, and exits non-zero on a FAIL; it imports the
first-round libraries (`../../computation/scripts/ecore.py`, `endo.py`, `k8.py`, `../../petri/scripts/enumerate_ec.py`,
`petri.py`, `../../../scripts/crn.py`) and, for the comparisons, the round-2 libraries of Q1 and Q2, all unchanged.
Shared helpers: `scripts/e2lib.py` (exact Fourier–Motzkin feasibility with a witness point, record loading).

## 1. Summary

- **Log cleaning (task 5).** 29 absolute paths `/home/<user>/gacE-data/…` in 23 first-round logs replaced by
  `~/gacE-data/…`; nothing else changed (§9).
- **s = 1, persistence of N* decided exactly (task 2).** Proposition E2.1 gives a complete decision (ALL / SOME /
  HYPERPLANE / NONE) from the exponents of Ye; on all 4 754 first-round networks with s = 1 and (H1), (H2) it leaves
  nothing undecided: 2 029 ALL, 0 SOME, 5 HYPERPLANE (persistent only on a hyperplane of rates), 2 720 NONE. Weakly
  reversible with (H1)–(H3): **150 ALL, 411 NONE, 0 SOME** (agrees with Q1's independent Theorem S1-WR: for WR N
  with s = 1 the dual is persistent for all rates or for none). A generic SOME network exists outside WR (§2).
- **s = 4 with molecularity 3 (task 4).** Seeded samples n = 5 (20 000) and n = 6 (10 000), WR, conservative, s = 4:
  every claim marked proved by P, G and E passes; **238 duals with a drain species at s = 4** (first ones at s = 4),
  smallest N = X₅ → 3X₄ → X₂+X₄ → X₁+2X₃ → 3X₃ → X₅ (one 5-cycle, deficiency 0, (H1)–(H3), globally attracting by
  Anderson 2011), whose dual loses k₅ for every choice of rate constants; no dual is endotactic (0 of 28 883).
- **K12 settled for WR networks (task 4).** Proposition E2.5 (proved, a Birch argument): if N is conservative and
  every product complex lies in the linear span of the source complexes (e.g. N weakly reversible), then N* is
  absolutely consistent. So **no weakly reversible network has a dual that fails absolute consistency**; exact
  log-cone criterion (Prop. E2.4) on 464 331 networks confirms it (94 464 WR with (H1)–(H3): all absolutely
  consistent) and agrees with Boros 2019, Theorem T and Corollary T* on every decided case.
- **WF-EC networks are never multistationary in the families (task 3; Q3-U).** Exact sign-vector decision (Prop.
  E2.3) on all 11 185 WF-EC networks of seven exhaustive families (n ≤ 4): 0 multistationary, including the 7 254
  where the toric exponent space differs from S (Birch does not apply); Q3's maximal-minor criterion agrees on all;
  Jacobians at sampled steady states are Hurwitz on the class. For the dual of the 52 unbalanced ones, msolve finds
  exactly one positive steady state in all 208 sampled (κ, class) pairs and exact two-point LPs find no pair.
- **s = 2 (task 2).** Of 18 861 WR networks with (H1)–(H3): N* certified persistent for all κ 2 988, certified
  non-persistent for all κ 1 048 (drain species 199, drain sets 142, and 707 by a new vertex-scaling certificate,
  Prop. E2.7), undecided 14 825 (§5.1).
- **Q1**: Q1-1 … Q1-5 pass (S1-WR also on 3 000 random networks beyond the records); **Q1-6 refuted** (72 pure WR
  networks whose dual has a drain set; smallest X₄ ⇄ X₂+X₃, 3X₃ ⇄ X₁+2X₄, s = 2).
- **Q2**: Q2-1, Q2-5, Q2-6, Q2-8 pass with an independent κ_V; Q2's K5 check "π⁺ = linkage classes" fails on 3 WR
  deficiency-0 networks (Q2-6 itself holds); Q2-12 answered by a two-line proof (Prop. E2.6: CF never merges two
  linkage classes of a WR network); Q2-13 fails on 8 875 of 20 189 WR networks (smallest 2X₃ ⇄ X₁+2X₂).
- **Q3**: every claim marked proved passes: Q3-C, Q3-U(b)(c), Q3-T*′ (9 011 WF-EC networks), Q3-Ustar and Q3-GS
  (7 982, exact Hurwitz/Metzler tests at random points), Q3-Cgen and K5gen (every record of both rounds; Cgen on
  9 318 sides, 3 307 not EC), Corollary BH (88 directly generated networks, 22 not EC) and Q3-PC (NP and 60 random
  PC networks neither EC nor BH, both sides, msolve counts = 1); Q3's examples NM, NB, NT, M3/M4 re-checked
  exactly; conjecture Q3-CR holds on 4 651 networks; T*m: all 2 003 unbalanced EC networks with s = m (n ≤ 3) have
  an absolutely consistent dual. K8 strong form: no species-side witness among 48 296 strong networks, 384
  reaction-side witness keys, smallest pair checked exactly (5 reactions, smaller than Q3's M3/M4).
- **G's Propositions 6.2 and 6.3** hold on every tested face (555 430) and siphon (103 222) (§8).

## 2. s = 1: exact decision of persistence of the dual (`scripts/s1_persistence.py`, `logs/s1_persistence.log`)

**Setting.** s = 1, (H1), (H2) for N. Then Γ = a bᵀ with a ∈ ℤⁿ primitive and b ∈ ℚʳ, and N* has reactions
i = 1..n with source Ye[i,:] and reaction vector a_i b. On the class k⁰ + ℝb (a compact segment) the dynamics is
t′ = g(t), g = Σ_i κ_i a_i k^{Ye[i,:]}. For σ = +1 (the end where the species k_j with b_j < 0 vanish) and σ = −1
(b_j > 0), call reaction i *bad* at σ if sign a_i = σ (it pushes towards that end) and *good* if sign a_i = −σ. For a
nonempty set J of species of N* on side σ put d_i(J) = Σ_{j∈J} Ye[i,j] (the order of term i at the end with zero set
J), e_i(J) = (Ye[i,j])_{j∉J} (its exponent in the other coordinates), β_i(J) = Π_{j∈J} |b_j|^{Ye[i,j]}.

**Lemma E2.0 (proved).** If a_i b_j < 0 then Ye[i,j] > 0. *Proof.* The j-th entry of the reaction vector a_i b of
dual reaction i is negative, i.e. the reaction consumes k_j, so k_j is in its source. ∎ (Checked: claim L, 4 754/4 754.)

**Proposition E2.1 (proved).** With this notation:
1. N* is persistent for every κ (ALL) iff for every σ and every nonempty J on side σ, min_good d(J) < min_bad d(J).
2. N* is non-persistent for every κ (NONE) if for some σ, J either min_bad d(J) < min_good d(J), or the two minima
   are equal and some tied bad exponent e_i(J) lies outside the convex hull of the tied good exponents.
3. Otherwise N* is non-persistent for some κ; it is persistent for some κ in particular when the strict linear
   inequalities obtained by bounding every tied bad monomial by a convex combination of tied good monomials (weighted
   AM–GM) are feasible in κ (SOME). If all nontrivial dual reactions have the same source y, then g = (Σ κ_i a_i) k^y
   and N* is persistent iff Σ_i κ_i a_i = 0 (HYPERPLANE).

*Proof.* On a segment, a trajectory starting inside converges monotonically to a zero of g or to an end, and it
converges to the σ-end iff σg > 0 on the part of the segment between it and that end. So N* is persistent for κ iff,
in every class, σg is not positive on a one-sided neighbourhood of the σ-end (g is a polynomial, so it either
vanishes identically or has a sign near the end). Every boundary point x⁺ with zero set exactly J (J on side σ) is
the σ-end of its own segment, and near it k_j = ε|b_j| (j ∈ J), k_j = x⁺_j + O(ε) (j ∉ J); hence
g = Σ_{d} ε^{d} c_d(x⁺) with lowest coefficient Σ_{i: d_i(J) = d_min} κ_i a_i β_i (x⁺)^{e_i(J)}, and the x⁺_{J^c} range
over all positive vectors. (1) If the good minimum is strictly smaller, the lowest coefficient is a sum of good terms
only, so σg < 0 near every σ-end, for all κ and all classes. Conversely, if for some J the bad minimum is ≤ the good
one, the bad terms dominate at x⁺_{J^c} = 1 when the bad κ_i are large (strict case: always), so some κ fails.
(2) Strict case: bad terms alone at lowest order, for every κ. Vertex case: the convex hull of the tied exponents has
a vertex v that is a bad exponent and not a good one (if all vertices were good exponents the hull would equal the
hull of the good ones); choosing x⁺ = exp(tw) with w maximised over the hull only at v and t → ∞, the term at v
dominates and is bad, for every κ. (3) If every tied bad exponent is a convex combination Σ λ_k e_k of tied good
exponents, x^{e_i} ≤ Σ λ_k x^{e_k} (weighted AM–GM), so if for each good exponent e_k the bad coefficients
Σ_bad κ_i |a_i| β_i λ_{i,k} are smaller than the good coefficient at e_k, the lowest coefficient is negative for every
x⁺; these are homogeneous strict linear inequalities in κ, solved exactly (Fourier–Motzkin with back-substitution).
The single-source case is immediate. ∎

This is the same statement as Q1's Theorem S1 (round2/persistence/notes.md: IN = condition (1) at a face, OUT/ESC =
the two cases of (2), DOM = the tie of (3)); the implementation here is independent of Q1's and the two agree on all
networks (§5, test A).

**Results (exact, all first-round record files, 66 s).**

| class (N) | ALL | SOME | HYPERPLANE | NONE | total |
|---|---|---|---|---|---|
| (H1), (H2) | 2 029 | 0 | 5 | 2 720 | 4 754 |
| (H1)–(H3) | 1 920 | 0 | 5 | 2 674 | 4 599 |
| WR with (H1)–(H3) | 150 | 0 | 0 | 411 | 561 |
| WR, deficiency 0 | 17 | 0 | 0 | 4 | 21 |
| (H1)–(H3) for N and N* | 870 | 0 | 3 | 678 | 1 551 |

NONE certificates: 1 119 by a strictly lower bad order, 1 601 by a Newton-polytope vertex. Consistency claims (PASS):
ADLS(N*) ⇒ ALL (924 networks), drain species of N* ⇒ NONE (75), WF-EC N ⇒ ALL (141, P's Corollary P), and a sampled
exact evaluation of g at distance 10⁻⁴⁰ from 14 261 ends agrees with the leading-order prediction.

- Smallest WR network with a NONE dual: N = 2X₂ ⇄ X₁+X₂ (one linkage class, deficiency 0, (H1)–(H3)),
  N* = k₂ → k₁, 2k₁+k₂ → k₁+2k₂: on the class k₁+k₂ = T, k₂′ = k₂(κ₂k₁² − κ₁), so every class with κ₂T² < κ₁ loses
  k₂, for every κ (a vertex certificate: the bad exponent 0 lies outside the good hull {2}). This is a class-dependent
  extinction that is not a drain species of E.4.
- HYPERPLANE (5 networks, none WR), smallest N = X₁+X₂ → 2X₂, X₁+X₂ → 2X₁: N* = k₁+k₂ → 2k₂, k₁+k₂ → 2k₁,
  k₂′ = (κ₁ − κ₂)k₁k₂: persistent iff κ₁ = κ₂.
- **Generic SOME example (constructed, outside the families; claim S, PASS).** N = X₁+X₃+X₄ → 2X₁+X₂,
  X₁+2X₂+X₃+X₄ → X₂+2X₃+2X₄, X₁+X₂+2X₄ → X₃+3X₄ ((H1), (H2), s = 1, not WR), with
  N* = k₁+k₂+k₃ → 2k₁, 2k₂+k₃ → k₁+k₂, k₁+k₂ → 2k₂+k₃, k₁+k₂+2k₃ → 2k₂+3k₃: the only tie is at the face {k₂ = 0},
  where the lowest coefficient is k₁(κ_A k₃ − κ_C − κ_D k₃²); N* is persistent for κ = (1, 1, 3/2, 3/2) and not for
  κ = (3, 3, 1, 1) (both certified exactly). Q1 found independent SOME examples (their S and B).

**Translation for the source.** In Iosif's dual, the old rate constants are the variables and the old concentrations
the parameters. For a one-dimensional network N (all reaction vectors parallel to a), group the reactions of N by the
direction in which they move (b_j > 0 or < 0). The dual system is persistent for all parameter values iff for every
nonempty set J of reactions moving in one direction, some species *produced* by the reactions of that direction has a
smaller total coefficient in the sources of J than every species *consumed* by them (for n = 2 and J a single
reaction this is Q1's C-2: each source contains more of its consumed than of its produced species). When this fails
strictly, or at a Newton vertex, the dual loses the reactions of some such J (their rate variables tend to 0) for
every parameter choice, in suitable classes; for weakly reversible N there is no intermediate case (Q1's S1-WR).
(At the end where the k_j with b_j < 0 vanish, the good dual reactions are the species i with a_i < 0, i.e. those
produced by the reactions j with b_j < 0, since Γ_ij = a_i b_j > 0.)

## 3. s = 4 with molecularity 3 (`scripts/run_m3s4.sh`, `logs/sample_*_3_*.log`, `logs/tests_m3s4*.log`)

Seeded samples with the first-round sampler (`sample_wr.py`, unchanged): n = 5, molecularity ≤ 3, 5 ≤ r ≤ 8,
conservative, s = 4, seed 3: 20 000 WR networks (all with (H1)–(H3), 17 894 of deficiency 0, 0 with WR dual; 159 s);
n = 6, seed 4: 10 000 (8 883 with (H1)–(H3); 90 s). The full first-round suite (`tests.py`, 15 s and 207 s):

- every claim marked proved holds (E-a1, E-a2, P-K1, K2, K2n, K7, K7r, H, C, C′, K10, K11, K14a–b, G, G-2.1, G-2.2,
  X-endoN, X-sendoN, X-drainN, X-L5);
- transfers: N* WR 0 / 28 883, deficiency 0 preserved 0 / 26 519, N* endotactic 0 / 28 883;
- **X-drainD: 238 of 28 883 duals have a drain species** (s = 4; none had been seen at s ≥ 3 with molecularity
  ≤ 2, and 103 at s = 3 with molecularity 3). Smallest: N = X₅ → 3X₄, 3X₄ → X₂+X₄, X₂+X₄ → X₁+2X₃, X₁+2X₃ → 3X₃,
  3X₃ → X₅; N* = k₅ → k₄, k₄ → k₂, 3k₃+2k₅ → 2k₄+3k₅, 3k₂+k₄ → 3k₁+k₂, k₁ → k₃, in which k₅ is a drain species
  (consumed alone by k₅ → k₄, produced only by a reaction with 2k₅ in its source). N is a single cycle with
  deficiency 0, so complex balanced and globally attracting for every rate choice (Anderson 2011).
- K6 refuted again at s = 4 (smallest N = X₆ → X₄+X₅ → X₂+2X₅ → X₁+2X₅ → 2X₃ → X₆).

## 4. Absolute consistency of the dual: K12 (`scripts/k12_abscons.py`, `logs/k12_abscons.log`)

**Proposition E2.4 (log-cone criterion; proved).** Let X have (H2), V = ker Γ ∩ ℝʳ_{>0} with extreme rays
v¹..v^p, K = ker Ye ⊆ ℝʳ, d = dim K, c¹..c^d a basis of K. Then X has a positive steady state for the rate vector
k iff some v ∈ V has log v − log k ∈ im Yeᵀ = K^⊥; so X is absolutely consistent iff h: V → ℝ^d,
h(v) = (c^l · log v)_l is onto. (a) d = 0: absolutely consistent. (b) d = 1: writing v = Σ λ_k v^k (λ > 0, every
point of V has such a representation), h(v) = Φ(log λ) + O(1) uniformly, with Φ(α) = Σ_j c_j max_{k: v^k_j > 0} α_k;
hence h is onto iff Φ takes both signs, iff Σ_j c_j ≠ 0 or there are sets A, B of rays with Φ(1_A) > 0 > Φ(1_B)
(Φ(1_A) = Σ of c_j over the reactions covered by a ray of A).
*Proof.* Steady state ⇔ v := k ∘ x^{y_·} ∈ V, and log v = log k + Yeᵀ log x. For (b): log(Σ_k λ_k v^k_j) differs
from max_{k: v^k_j>0} log λ_k by at most a constant depending on the rays; Φ is positively homogeneous and, with
Σc = 0, invariant under α ↦ α + t1; on each cone of a fixed order of the α_k it is linear and the cone is generated
by ±1 and the indicators of initial segments, so sup Φ > 0 iff some Φ(1_A) > 0 (every set is an initial segment of
some order), and similarly for inf; h is onto iff it is unbounded in both directions (V is connected). ∎

**Proposition E2.5 (Birch; proved).** If X has (H2) and K = ker Ye ⊆ ker Γ, X is absolutely consistent.
*Proof.* Take v⁰ ∈ V. By Birch's theorem, for every k the set exp(log k + K^⊥) meets the affine space v⁰ + K, inside
the open orthant, in exactly one point v; and v⁰ + K ⊆ ker Γ. So v ∈ V and log v − log k ∈ K^⊥. ∎
**Corollary E2.5′ (answers K12(b)).** If N is conservative and every product complex of N lies in the linear span of
its source complexes (in particular if N is weakly reversible, or if every product complex is a source complex),
then N* is absolutely consistent. *Proof.* For X = N*, Ye(X) = Ye(N)ᵀ and Γ(X) = Γᵀ; cᵀYe = 0 implies cᵀYp = 0
because the columns of Yp are combinations of those of Ye, so cᵀΓ = 0: apply E2.5 (N* has (H2) because N has (H1)). ∎
So the question "is there a WR network N whose dual is not absolutely consistent?" has a negative answer, by a
one-line argument; positive steady states of N* need not be unique, nor complex balanced.

**Results (464 331 networks, both sides, 150 s).** Verdicts on side N: AC0 40 288, ACB 633, AC1 51 793, NOT1 1 496,
OPEN (d ≥ 2) 370 121; on side N*: AC0 96 639, ACB 30 725, AC1 936, NOT1 644, OPEN 33. Claims (PASS): every WR
network with a decided verdict is absolutely consistent (79 309; Boros, SIAM J. Math. Anal. 51 (2019)); EC, connected,
(H1), (H2): absolutely consistent iff s = m − 1 (5 707 decided; P's Theorem T); balanced EC: same verdict for N and N*
(1 722; Corollary T*); U is AC0 and U* is NOT1 (bounded above, c = (−1, 1, 0, 0)). **K12(a)**: the 177 EC networks
with unbalanced joins, connected, (H1), (H2), s = m − 1 all have an absolutely consistent dual (42 AC1, 135 ACB),
as Q3's Theorem T*′ states. **K12(b)**: all 94 464 WR networks with (H1)–(H3): N* absolutely consistent (79 678 AC0,
14 786 ACB), as Corollary E2.5′ proves.

## 5. Agent Q1's candidates (`round2/persistence/candidates.md`; `scripts/q1_tests.py`, `logs/q1_tests.log`, 198 s)

| id | statement (short) | result |
|---|---|---|
| Q1-1 / Q1-9 | Theorem S1 (face criterion, rank one) and its counts | **PASS**: my independent decision (§2) and Q1's `classify` agree on all 4 754 first-round networks with s = 1 and (H1), (H2): 2 029 ALL/ALL, 2 720 NONE/NONE, 5 HYPERPLANE/TWO-SIDED; WR with (H1)–(H3): 150 ALL, 411 NONE, as Q1 reports |
| Q1-2 | S1-WR: N weakly reversible, s = 1 ⇒ N* is ALL or NONE | **PASS** beyond the records: 3 000 seeded random WR networks with s = 1 and (H1)–(H3) (n = 2–4, up to 3 linkage classes, 2–5 complexes per class, molecularity ≤ 5): 923 ALL, 2 077 NONE, 0 SOME (both implementations) |
| Q1-3 | n = 2, s = 1: N* ALL iff every source has more of its consumed than of its produced species | **PASS** on 3 212 networks |
| Q1-4 | s = 1: no net product in a source ⇒ ALL | **PASS** on 892 networks |
| Q1-5 | drain set ⇒ extinction for every κ | consequence checked (**PASS**: in 2 856 WR duals with a drain set, the drain set never contains the support of a conservation law of N*); numerics: 30/30 drain-set duals decay from x_Z = 10⁻³ (LSODA, random κ) |
| Q1-6 | N WR and pure ⇒ N* persistent for every κ (conjecture) | **REFUTED**: 72 counterexamples among 8 195 pure WR networks with (H1)–(H3) (s = 1: 14, s = 2: 1 471, s = 3: 3 317, s = 4: 3 393). Smallest: N = X₄ ⇄ X₂+X₃, 3X₃ ⇄ X₁+2X₄ (pure, WR, two linkage classes, deficiency 0, (H1)–(H3), s = 2), N* = k₄ → k₂, k₃ → k₁, 3k₂+k₃ → k₁+3k₄, k₁+2k₄ → 2k₂+k₃: Z = {k₃, k₄} is a drain set (w = (3, 1), d = (1, 1): k₃ → k₁ and k₄ → k₂ decrease V = 3k₃+k₄, 3k₂+k₃ → k₁+3k₄ leaves it unchanged, k₁+2k₄ → 2k₂+k₃ has Σ y/d = 2 > 1), so N* loses k₃, k₄ for every κ |
| Q1-7, Q1-8 | dichotomy at s = 2; OUT facet ⇒ NONE | see §5.1 |

### 5.1 s = 2: certificates (`scripts/s2_persistence.py`, `logs/s2_persistence.log`)

All 18 861 distinct WR networks with (H1)–(H3) and s = 2 of the first-round files. Certificates, all exact and
rate-independent: ALL = N* has no critical siphon (ADLS Theorem 2); NONE = drain species (E.4), drain set (Q1's
Theorem D), or the vertex-scaling certificate of Proposition E2.7 below.

**Proposition E2.7 (vertex scaling; proved).** Let X be conservative with s = 2 and Z a set of species such that
(i) rank Γ[Z,:] = 2, (ii) every reaction source meets Z, (iii) W = {d_Z : d ∈ S} contains a strictly negative
vector. For the reactions k whose source is e_j + w_k (j ∈ Z, w_k supported off Z) put M_k = (Δ_k|_Z) e_jᵀ
restricted to W (rank one), and group them by |w_k|. If in the group of lowest (or highest) |w_k| every trace
tr M_k is ≤ 0, one is < 0, every mixed determinant coefficient m(M_k, M_l) = a₁₁b₂₂ + a₂₂b₁₁ − a₁₂b₂₁ − a₂₁b₁₂ is
≥ 0 and one is > 0, then X is non-persistent for every κ.
*Proof.* By (iii), for a suitable positive u the class of u has a vertex v with zero set exactly Z (by (i) the face
x_Z = 0 of a class is at most a point); by (ii) every reaction vanishes at v, so v is an equilibrium for every κ.
The only rates with a nonzero derivative at v are those with source e_j + w_k, and the linearisation, read in the
coordinates x_Z (injective on S by (i)), is L = Σ_k κ_k v^{w_k} M_k on W. In the class of Tu the vertex is Tv, so
L(T) = Σ_k κ_k T^{|w_k|} v^{w_k} M_k; as T → 0 (or ∞) the extreme group dominates: L(T) = T^{δ}(L₀ + o(1)). Here
tr L₀ = Σ κ_k v^{w_k} tr M_k < 0 and det L₀ = Σ_{k<l} κ_kκ_l v^{w_k+w_l} m(M_k, M_l) > 0 (det M_k = 0), for every
κ, so v is linearly asymptotically stable within the class for T small (large); its basin contains interior points,
whose ω-limit is v. ∎
Example: the smallest undecided network of the first pass, N = 2X₃ → X₂ → X₁+X₃ → 2X₃ (one cycle, deficiency 0),
N* = k₃ → k₂, k₂ → k₁, 2k₁+k₃ → k₂+2k₃: Z = {k₂, k₃}, the vertex (T, 0, 0) has eigenvalues −κ₂ and κ₃T² − κ₁, so
every class with κ₃T² < κ₁ loses k₂ and k₃, for every κ (traces −1, −1 and mixed coefficient 1 at T → 0).

**Results (18 861 networks, 245 s).** ALL (no critical siphon in N*): **2 988**. NONE, rate-independent: **1 048**
(drain species 199, drain set 142, vertex scaling 707). Undecided: 14 825; facet types present among them (IN /
OUT / ESC as in Prop. E2.1 with J = the facet's zero set): none 3 235, ESC only 3 370, ESC+IN 3 190, ESC+IN+OUT 609,
ESC+OUT 405, IN only 2 783, IN+OUT 585, OUT only 648. Smallest undecided: N = 2X₃ → X₂ → 2X₁+X₃ → 2X₃,
N* = 2k₃ → 2k₂, k₂ → k₁, 2k₁+k₃ → k₂+2k₃ (near k₃ = 0, k₃′ = k₃(κ₃k₁² − 2κ₁k₃) > 0: probably ALL, no certificate).

Sampled rates on the 200 smallest undecided networks (4 rational κ each, classes scaled by 10^{−6..6}): an exact
asymptotically stable boundary vertex (hence non-persistence for that κ) for all 4 κ in 15 networks, for some in
2, for none in 183. The LSODA runs (one start per coordinate, 3 scales) are **not** reliable evidence for Q1-7: on
the smallest network flagged "mixed" by numerics in a first pass, the exact analysis (Prop. E2.7) proves NONE —
extinction happens only in classes of small total mass, which a fixed-scale simulation misses. No network was found
that is certified non-persistent for one κ and certified persistent for another, so **Q1-7 (dichotomy at s = 2) is
neither refuted nor supported by a certificate**; one undecided network has an ESC facet with no extinction seen
(N = X₃+2X₄ → 2X₁+X₂ → X₂+X₃+X₄ → X₃+2X₄), a candidate to examine for Q1-8.

## 6. Agent Q2's candidates (`round2/preclusters/candidates.md`; `scripts/q2_tests.py`, `logs/q2_tests.log`)

Independent implementation of κ_V (classes c ~ c′ iff u_c = u_c′ for every u ∈ ker A_V; Q2 uses a rank test on the
row space), V = random positive integer combinations of the flux rays; Q2's π⁺ (lower bound of CF) imported. Run on
the families Q2 did not run: WR (2,3,6), (3,3,5), (4,2,6), (5,2,6) homogeneous, (5,3,4), the samples n = 5, 6
(molecularity 2) and the new molecularity-3 samples (s = 4), plus all-networks (3,3,3), (2,3,4); at most 4 000
networks per file (seeded): 42 784 networks, 73 485 sides with (H2).

| id | statement (short) | result (sides with (H2); 259 s) |
|---|---|---|
| Q2-1 | P ⪯ Pcs ⪯ π ⪯ CF ⪯ κ_V; tested: soundness π⁺ ⪯ κ_V at 2 generic V | **PASS** 73 485 / 73 485 |
| Q2-5 | EC + connected + (H1) + (H2) ⇒ κ_V one block | **PASS** 1 578 / 1 578 |
| Q2-6 | WR with deficiency 0 ⇒ κ_V = CF = linkage classes | **PASS** 25 867 / 25 867 (κ_V at 2 generic V) |
| Q2's K5 check | WR, deficiency 0 ⇒ π⁺ = linkage classes (i.e. CF determined by Q2's bounds) | **fails on 3 sides** (not a failure of Q2-6): smallest N = 2X₂ → 3X₂, 2X₂ → 3X₁, 3X₂ → 2X₂, 3X₂ → 3X₁, 3X₁ → 2X₂, 3X₁ → 3X₂ (complete digraph on three complexes, deficiency 0, not conservative) and X₂ → 2X₂, X₂ → X₁+2X₂, 2X₂ → X₂, 2X₂ → X₁+2X₂, X₁+2X₂ → X₂, X₁+2X₂ → 2X₂: κ_V is one block, π⁺ = Pcs = π are the three source classes. Family wr_2_3_6 (not in Q2's run) |
| Q2-8 | dim ker A_V ≤ #κ_V | **PASS** 73 485 / 73 485 |
| Q2-3 | M* ⪯ π, and π = M* ⇔ κ_V = M* (generic V) ⇔ every doubling pair fully coupled (Q2's `full_coupling`, `flat_closure`) | **PASS** 73 485 / 73 485 |
| Q2-10 | (H1): every block of CF(N*) lies in a minimal siphon of N ⇔ no twin pair and every ray-MIPS lies in a minimal siphon | **PASS** on the 30 701 sides where CF(N*) is determined (π⁺ = meet of 3 κ_V) |
| Q2-12 | for WR, CF merges two linkage classes only if δ ≥ 1; find the smallest | **never happens** (Proposition E2.6 below); 36 107 WR networks: CF ⪯ linkage classes certified for all |
| Q2-13 | every minimal siphon of N lies in one block of CF(N*) | **fails often**: holds (certified by π⁺) 20 764, fails (certified by the κ_V meet) 9 937 of 30 701 with (H1); with (H1)–(H3) 17 967 / 8 875; WR with (H1)–(H3) 11 314 / 8 875. Smallest failure (WR, (H1)–(H3)): N = 2X₃ ⇄ X₁+2X₂ (minimal siphons {X₁,X₃}, {X₂,X₃}; N* = k₂ → k₁, 2k₂ → 2k₁, 2k₁ → 2k₂ has three distinct sources, and the meet of three generic κ_V(N*), which CF(N*) refines, puts no minimal siphon of N inside one block) |

**Proposition E2.6 (answers Q2-12; proved).** If every linkage class L of N carries a positive flux of its own
(V_L > 0 with Γ_L V_L = 0; e.g. N weakly reversible), then CF(N) refines the linkage-class partition, whatever the
deficiency. *Proof.* V = (V_L)_L lies in C°. For the indicator u = 1_L of the source complexes of L,
A_V u = Σ_{t∈L} V_t (y′_t − y_t) = Γ_L V_L = 0, so u ∈ ker A_V; for c ∈ L and c′ ∉ L, (e_c − e_c′)·u = 1 ≠ 0, so
e_c − e_c′ is not in the row space (ker A_V)^⊥. Hence κ_V refines the linkage classes and CF, the meet of all κ_V,
does too. ∎ So "CF merges two linkage classes of a WR network" never happens (Q2-12 is vacuous: the conjectured
condition δ ≥ 1 is not needed and not sufficient); in the families the meet of three generic κ_V already refines the
linkage classes for all 36 107 WR networks tested.

## 7. Agent Q3's candidates (`round2/ec/candidates.md`; `scripts/q3_tests.py`, `scripts/ms_wfec.py`,
`scripts/ms_dual.py`, `scripts/k8_strong.py`)

**Proposition E2.3 (multistationarity criterion for WF-EC; proved, classical).** Let N be WF-EC with sources
y_1..y_m, L = span{y_c − y_c′}, S = im Γ. Some rate vector k and some compatibility class carry two distinct positive
steady states iff there are nonzero d ∈ S and w ∈ L^⊥ with sign d = sign w.
*Proof.* By P's Theorem T, the positive steady states for k are {x > 0 : x^{y_c} = λ u_c(k)}, i.e. log x ∈
log x* + L^⊥, and every u > 0 is u(k) for k_t = v_t/u_{c(t)} (v > 0 in ker Γ). Two steady states x₁ ≠ x₂ in one
class give d = x₁ − x₂ ∈ S and w = log x₁ − log x₂ ∈ L^⊥ with the same signs. Conversely, given d, w, put
x₂,i = d_i/(e^{w_i} − 1) (x₂,i = 1 if w_i = 0), x₁ = x₂ e^w, u_c = x₁^{y_c} and k as above; x₂^{y_c} = e^{−w·y_c}u_c
with w·y_c independent of c. ∎ This is the sign criterion of Pérez Millán–Dickenstein–Shiu–Conradi (2012) and
Müller–Feliu–Regensburger–Conradi–Shiu–Dickenstein (2016) specialised by Theorem T; Q3 proves the stronger statement
that the condition never holds (Q3-U). With w integral, ρ = 2^w gives an exact rational witness (claim M1).

**Exhaustive WF-EC families (P's generator, weighted; `logs/ms_wfec_*.log`, `logs/q3_tests_*.log`).**

| family (n, molecularity, r) | WF-EC | unbalanced | L = S (Birch) | L ≠ S | multistationary |
|---|---|---|---|---|---|
| (2, 3, ≤5) | 71 | 0 | 37 | 34 | 0 |
| (2, 4, ≤5) | 355 | 0 | 147 | 208 | 0 |
| (3, 2, ≤5) | 958 | 0 | 543 | 415 | 0 |
| (3, 2, ≤6) | 2 080 | 0 | 1 291 | 789 | 0 |
| (3, 3, ≤4) | 6 995 | 52 | 1 692 | 5 303 | 0 |
| (4, 2, ≤4) | 692 | 0 | 206 | 486 | 0 |
| (4, 2, ≤5) ordinary | 34 | 0 | 15 | 19 | 0 |
| **total** | **11 185** | 52 | 3 931 | 7 254 | **0** |

The families (5, 2, ≤5) ordinary and (3, 3, ≤5) were started and stopped by hand (machine overloaded; they would
have exceeded the 25-minute cap); no result is claimed for them.

Q3 claims on the same families (`q3_tests.py`; (2,3,5), (2,4,4), (3,2,5), (3,3,4), (4,2,4): 9 011 WF-EC networks):

| id | test | result |
|---|---|---|
| Q3-C | no critical siphon (minimal siphons = supports of extreme conservation laws), impure and unbalanced included | **PASS** 9 011 / 9 011 |
| Q3-U (c) | no common sign vector of S and L^⊥, decided by exhaustive sign vectors (exact LPs) and by Q3's maximal-minor criterion; the two agree | **PASS** 9 011 / 9 011, agreement 9 011 / 9 011 |
| Q3-U (b) | at 2 sampled steady states (x rational, v > 0 in ker Γ, k = v/x^y): det(tI − J) = t^{n−s} q(t), q Hurwitz (exact Hurwitz determinants) | **PASS** 9 011 / 9 011 |
| Q3-T*′ | N* absolutely consistent (exact verdict of Prop. E2.4) | **PASS** 9 011 / 9 011 (also the 177 unbalanced WF-EC networks of the first-round records, §4) |
| Q3-T*-cex / T*m | EC, connected, (H1), (H2), s = m: N* absolutely consistent? | balanced: none of 933 is certified absolutely consistent (907 NOT1, 26 OPEN with d ≥ 2); unbalanced: **all 2 003** are (AC1) — consistent with Q3's Prop. T*-dim; no counterexample to "unbalanced ⇒ N* absolutely consistent" in these families |
| Q3-Udual | multistationarity of N* for unbalanced WF-EC N (`ms_dual.py`) | (3,3,≤4): 52 networks; 300 exact "two-point" trials each (choose x₂, x₁ = x₂ + d in one class of N*, decide by an exact LP whether some κ makes both steady states): no certificate; msolve on 4 random (κ, class) each: **exactly one positive steady state in all 208 samples**. The run on (3,3,≤5) hit its 25-minute cap without output (machine overloaded) |
| Q3-T*m | EC, unbalanced, s = m: N* absolutely consistent? (`tstar_m.py`) | the exhaustive (4,3,≤4) run hit its cap twice (25 min under load, then 20 min on 10 processes; one partition of the species carries most of the work, as in P's first round) without output; only the n ≤ 3 evidence above (2 003 / 2 003) |
| Q3-Ustar | WF-EC N, κ random: at 3 random positive points k of N* (not steady states), det(tI − J*) = t^{r−s}q(t), q Hurwitz (`q3_tests2.py`) | **PASS** on 7 982 WF-EC networks ((3,3,≤4) 6 995 with 52 unbalanced, (4,2,≤4) 692, (2,4,≤4) 295) × 6 points |
| Q3-GS | G = Dφ(x)A_k Metzler with G u(k) = 0 at random rational (x, k) | **PASS** on the same 7 982 networks × 3 points |
| Q3-CR (conjecture) | connected, (H1), (H2), s = #clusters − 1 ⇒ no critical siphon in N, N* and N absolutely consistent | **holds** on all 4 651 first-round networks of that kind (2 033 not EC): 0 critical siphons; verdicts AC0 1 611, AC1 1 517, ACB 74, OPEN 1 449 (d ≥ 2, undecided), NOT1 0. Part (ii) (no multistationarity) not tested on this class except through PC below |
| Q3-PC | PC (sources sharing a species positively proportional), connected, (H1), (H2), s = #clusters − 1: N and N* without critical siphon, absolutely consistent, class Jacobian Hurwitz at every positive point, one steady state per class (`pc_tests.py`) | **PASS** on Q3's NP and 60 seeded random PC networks that are neither EC nor block-homogeneous, both sides: no critical siphon 61 + 61, absolutely consistent (decided) 61 + 61, Hurwitz at 6 random points 61 + 61, msolve exactly one positive steady state in 2 samples 61 + 61 |
| BH | block-homogeneous K (sources a·1_B), connected, (H1), (H2), s = #blocks − 1: K* WF-EC, no critical siphon in K, K*, absolutely consistent, Jacobian on S Hurwitz everywhere, one steady state per class (`bh_tests.py`, generated directly, seeded) | **PASS** on all 88 networks found (22 not EC): K* WF-EC 88, no critical siphon 88 + 88, absolute consistency 85 decided AC (3 OPEN), Hurwitz at 6 random points 88, msolve exactly one positive steady state in 2 samples 88 |
| Q3-B examples | NM, NB multistationary; NT* without positive steady state for κ₁ ≠ κ₂ (`q3_examples.py`) | **PASS**: msolve isolates exactly 2 positive steady states for NM (k = (1,2,1), class x₁+x₂+x₃ = 5/2) and for NB (k = 1, class 2x₁+x₂+x₃ = 6), and 0 for NT* (κ = (1,2,1,1), class of (1,2,3)); NT*'s log-cone verdict is OPEN (d = 2), so its non-absolute-consistency rests on Q3's argument plus this sample |
| Q3-Cgen, K5gen | connected, (H1), (H2): s ≥ #clusters − 1; with equality, no critical siphon in N nor N* (`final_retest.py`, every record of both rounds, 494 298 networks) | **PASS**: K5gen on 315 466 connected sides; Cgen on 9 318 sides with s = #clusters − 1, of which 3 307 are not EC |
| Q3-K8s | K8 strong form (both networks (H1)–(H3), rows and columns of Γ pairwise linearly independent) | species side: **no witness** among 48 296 strong networks (n ≤ 6, every labelling); reaction side: **384 witness keys** among 45 419 (r ≤ 6). Smallest, checked exactly (`k8_pair_check.py`): N₁ = 2X₃ → 2X₂, 2X₃ → X₁+X₂, 2X₂ → X₁+X₃, X₁+X₃ → 2X₃, X₁+X₂ → 2X₂ (s = 2, δ = 1) and N₂ = 2X₄ → X₃, X₃ → X₂+X₄, X₃ → X₁+X₄, X₂+X₄ → 2X₄, X₁+X₄ → 2X₄ (s = 3, δ = 0), both WR, (H1)–(H3), strong; under the 4 relabellings that match the doubling classes ({R₁,R₂} resp. {R₂,R₃}) and the minimal T-traps ({R₄}, {R₅}), P(N₁) has 4 blocks and P(N₂) one. Q3's pair M3/M4 re-checked independently (`q3_examples.py`): both (H1)–(H3), WR, pure, strong; 8 relabellings match the doubling classes and minimal T-traps, P differs under all 8 (P(M3) = {R₁,R₂,R₅},{R₃},{R₄},{R₆}; P(M4) = {R₁},{R₂,R₃},{R₄},{R₅},{R₆}) — **PASS** |

## 8. Agent G's Propositions 6.2 and 6.3 (`scripts/g6_props.py`, `logs/g6_props.log`, 207 s)

Seeded subsample of at most 1 500 conservative networks per first-round file (20 230 networks), two classes each
(x⁰ = 1 and a random rational x⁰).

- **6.2 (PASS, 555 430 realizable faces).** Independent computation: exact vertex enumeration of the class polytope,
  F_Σ = convex hull of the vertices with x_Σ = 0, realizable iff their barycentre has zero set exactly Σ, dimension
  = affine rank. In every case codim F_Σ = rank Γ[Σ,:] = rank of the dual subnetwork with the reactions Σ; every
  realizable zero set that is a siphon of N is a T-trap of N*.
- **6.3 (PASS, 103 222 siphons of 15 446 WR networks).** (i) products of J₀ avoid Σ (true for every siphon, WR not
  needed); (ii) Sbl ⊆ S₀; (a) dim S₀ = dim Sbl iff rank M[J₀,:] + rank M[:,Σ] = s, with S₀ computed independently as
  Γ(ker Γ[Σ,:]); (b) on every embedded linkage class (projection of the complexes to Σ), "im Γ[Σ,J_L] contains a
  nonzero nonnegative vector" and "ker M[J_L,Σ] contains a positive vector" are decided by two independent exact LPs
  and are always complementary (Stiemke). What is not tested: that (6.3)–(6.4) are exactly Wiuf's hypotheses
  (Σ-minimality, realizability in his sense), which needs his preprint's definitions; G's translation is linear
  algebra and is correct as stated.

## 9. Cleaning of the first-round logs (task 5)

In `../../computation/logs/`, the 23 files listed in `logs/logclean_counts.txt` (enum_all_*, enum_fc_*, enum_wr_*,
sample_*, facts.log) contained 29 occurrences of the absolute WSL path of the record directory
(`/home/<user>/gacE-data/`). Each was replaced by `~/gacE-data/` (the directory is outside the repository, so a
repository-relative path does not exist; `~` carries no user name). No other byte was changed; a final grep finds no
user name in `computation/logs/`. Not changed (not a log, other agent's file): `../../computation/notes.md` lines 59
and 329 still contain the same absolute path.

## 10. Scripts and logs

All under WSL Fedora 44, Python 3.12; run from `round2/computation/` (`python scripts/<name>.py …`, or from
`scripts/` for the ones that import siblings); exact rational arithmetic for every claim; runtimes in the logs.

| script | checks | log |
|---|---|---|
| `e2lib.py` | library: exact Fourier–Motzkin feasibility with witness (`fm_point`, `lp_point`, `in_convex_hull`), record loading, WF-EC test | – |
| `run_m3s4.sh` | s = 4, molecularity 3 samples (first-round `sample_wr.py`) and the first-round `tests.py` on them | `sample_5_3_5_8_s3.log`, `sample_6_3_5_8_s4.log`, `tests_m3s4.log`/`.json`, `tests_m3s4X.log`/`.json`, `data_sha256_m3s4.txt` |
| `s1_persistence.py` | Prop. E2.1 on all s = 1 networks; claims F, L, A, B, W, K, T, S | `s1_persistence.log` |
| `s2_persistence.py` | s = 2 certificates (ADLS, drain species, drain sets, Prop. E2.7), facet types, exact stable-vertex certificates for sampled κ, LSODA evidence | `s2_persistence.log` |
| `k12_abscons.py` | Props. E2.4, E2.5; claims B, T, S, U; K12(a), K12(b) | `k12_abscons.log` |
| `ms_wfec.py`, `run_ms.sh` | Prop. E2.3 on exhaustive WF-EC families, exact witnesses | `ms_wfec_*.log`, `run_ms_status.log` |
| `q3_tests.py`, `run_q3.sh` | Q3-C, Q3-U(b)(c), Q3-T*′, T*-cex report | `q3_tests_*.log`, `run_q3_status.log` |
| `q3_tests2.py`, `q3_examples.py` | Q3-Ustar, Q3-GS, Q3-CR, Q3-B examples (msolve), Q3-K8s pair | `q3_tests2.log`, `q3_examples.log` |
| `bh_tests.py`, `pc_tests.py` | Q3's Corollary BH and Q3-PC on directly generated networks (exact; msolve counts) | `bh_tests.log`, `pc_tests.log` |
| `final_retest.py` | final re-test of Q3-Cgen and K5gen on every record of both rounds | `final_retest.log` |
| `ms_dual.py` | Q3-Udual: exact two-point LP certificates and msolve counts for N* | `ms_dual_*.log` |
| `tstar_m.py` | Q3-T*m on the EC family (4, 3, ≤4) | `tstar_m_4_3_4.log` |
| `k8_strong.py`, `k8_pair_check.py` | K8 in the strong class (both sides), exact re-check of a pair | `k8_strong.log`, `k8_pair_check.log` |
| `q1_tests.py` | Q1-1 … Q1-6 (A–G) | `q1_tests.log` |
| `q2_tests.py` | Q2-1, Q2-5, Q2-6, Q2-8 (independent κ_V), K5 form, Q2-12, Q2-13 | `q2_tests.log` |
| `g6_props.py` | G's Props. 6.2, 6.3 | `g6_props.log` |

Data: `~/gacE2-data/sample_5_3_5_8_s3.pkl.gz` (SHA-256 0968f651…6e1108) and `sample_6_3_5_8_s4.pkl.gz`
(6ff4ec1e…bdf61db), full hashes in `logs/data_sha256_m3s4.txt`. `logs/logclean_counts.txt`: task 5.
Note: two scipy warning lines that printed a home path were deleted from `logs/s2_persistence.log` of the first
run; the script now suppresses warnings.

## 11. What is left

- s = 2: 14 825 of the 18 861 WR duals remain without a rate-independent certificate. Needed: a persistence
  certificate beyond ADLS (e.g. Q1's facet criterion IN at every facet plus a vertex condition, or a tropical
  version of Prop. E2.7 for repelling vertices), and Prop. E2.7 for edges (boundary equilibria inside facets) and
  for non-hyperbolic vertices. Q1-7 (dichotomy at s = 2) and Q1-8 (OUT facet ⇒ NONE) are open; the exact sampled
  vertex certificates found no counterexample.
- Q3-Udual (multistationarity of N* for unbalanced WF-EC N): no instance in the families tested; Q3 has since
  proved uniqueness (Q3-Ustar); its exact test at random points is in `q3_tests2.py`.
- Q3-T*m: all 2 003 unbalanced EC networks with s = m in the families have an absolutely consistent dual; a
  counterexample, if any, needs n ≥ 4 with molecularity ≥ 3 (`tstar_m.py` (4, 3, ≤4)).
- K12 for networks whose products are not in the span of the sources, with d = dim ker Ye ≥ 2 (33 duals and
  370 121 primal sides OPEN): the tropical criterion of Prop. E2.4 is exact only for d ≤ 1; for d ≥ 2, surjectivity of
  Φ is necessary, sufficiency was not proved.
- The exhaustive families stop at n ≤ 4 (WF-EC weighted) and r ≤ 6; s = 4 with molecularity 3 is covered by seeded
  samples only.
- `../../computation/notes.md` lines 59 and 329 still print the absolute record path (other agent's file, not a log).
