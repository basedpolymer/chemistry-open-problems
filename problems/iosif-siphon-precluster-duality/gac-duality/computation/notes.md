# Agent E: exhaustive computations and red team for the gac-duality round

**Agent E: Claude Opus 5.5, effort xhigh. 27 October 2026.** Working notes for result 17
(`problems/iosif-siphon-precluster-duality/`), questions Q-A and Q-D of `../BRIEF-exploration.md`. Nothing here
quotes the author's e-mail. Status tags: **proved** (complete proof below or in the manuscript), **checked**
(exact computation on the finite families listed, not a proof), **refuted** (explicit counterexample, checked
exactly), **open**.

## Summary

- **Families** (all exact, up to species permutation): 420 986 weakly reversible networks in 11 exhaustive
  families (n ≤ 6, molecularity ≤ 2 or 3, r ≤ 4–6; families overlap slightly), two seeded samples (n = 5, 6;
  40 000 conservative WR networks with s ∈ {2,3,4}), 18 527 general networks with (H1), (H2) (counted per family), 84 ordinary EC
  networks. Distinct WR networks with (H1)–(H3): 94 645, of which s = 2: 18 955, s = 3: 36 068, s = 4: 39 061.
  Generators cross-checked against result 17 (all counts agree).
- **Coordinator's facts** 1 and 5 hold on 769 316 networks; facts 2 and 6 hold, but the counts 32/9/24/8 and
  9/23 double-count the 3-cycle T: correct 31/8/23/7 and 9/22.
- **Generalised refinement**: in the GAC setting the exact exceptional classes are "some complex has two outgoing
  reactions" (reaction side, E.1) and "im Γ contains e_i − μe_j" (species side, E.2); a Galois connection between
  the lattices always exists and is trivial (E.3); every "union" and "singletons removed" reading is refuted already
  on two species; only P's K13 (non-singleton preclusters of N* are unions of MIPS) survives.
- **Agents P and G**: every claim marked proved holds on every applicable network (K1–K3, K10, K11, K13, K14, C, C′,
  G, K2 hypothesis-free, Lemma 5 variants, G's Props 2.1, 2.2, 4.1, 7.1, Examples §3, §5, 8.1, the 705-network
  counts), except **K7**, whose last clause fails on 26 WR networks with the zero complex (now restricted by P;
  the restricted form holds). K4 and K5 hold; K9(iv) holds on all 5 429 EC sides; K6 and every ADLS-transfer
  statement are refuted; K8 witnesses found.
- **Transfers** (82 377 WR networks with H: all families but (4,3,4), (5,3,4), and both samples): N* is WR for 98, has
  deficiency 0 for 43 of 46 709 (deficiency-0 ones), is endotactic for 122 of 84 103 and strongly endotactic for 43
  of 70 185 (these two counted per run, about 1 700 networks counted twice); the
  ADLS criterion is not invariant in either direction, even when both N and N* are WR with (H1)–(H3) (8 of 87).
- **Persistence is not self-dual, rate-independently** (Prop. E.4, "drain species"): G's D is non-persistent for
  every rate choice; 141 + 162 WR networks with (H1)–(H3) have a dual with a drain species, including s = 3
  networks with one linkage class and deficiency 0 (smallest X₄ → 3X₃ → X₁+2X₂ → 3X₂ → X₄), globally attracting by
  Anderson's theorem while the dual (also (H1)–(H3)) loses a species for every choice of rate constants.

## 0. Conventions

Network N: species 1..n, reactions y_j -> y'_j (j = 1..r), Ye, Yp, Γ = Yp − Ye, s = rank Γ, q = n − s. Dual N*:
species k_1..k_r, reaction i: k^{(Ye)_i·} -> k^{(Yp)_i·}; Γ(N*) = Γᵀ, N** = N. (H1) conservative, (H2) consistent,
(H3) rows of Γ pairwise distinct; "H" = (H1)–(H3) for N, "H*" = (H1)–(H3) for N*. P = preclusters (Definition 19,
one-step graph; all variants agree, manuscript Lemma 5), M = ray MIPS (Prop. 3), S = minimal siphons, SF =
supports of extreme semi-positive conservation laws ("semiflow supports"), IC = internal cycles = SF(N*).
A minimal siphon is *critical* if it contains no semiflow support; "ADLS" = no critical siphon = every siphon
contains a semiflow support (Angeli–De Leenheer–Sontag). Petri notions follow agent P (`../petri/candidates.md`):
ordinary, EC (two reactions sharing a reactant have the same source), FC = ordinary + EC, balanced joins, JF, CF,
clusters (components of the bipartite graph species–reactions with an edge iff (Ye)_ij > 0), rank equation
s = #clusters − 1, WF-FC = connected + FC + (H1) + (H2) + rank equation. Classes in the tables: "all" (every
network of the files, hypothesis of the candidate permitting), H, strong (H + rows and columns of Γ pairwise
linearly independent), WR (H + weakly reversible), WR0 (WR + deficiency 0). For candidates tested "on both
sides", side N is the network itself and side D is its dual; the class columns always refer to N.

All computations are exact (`fractions.Fraction`, the linear algebra of `../../scripts/crn.py`, imported and not
edited). Networks are counted up to permutation of the species (canonical form `ecore.canonical`, a complete
invariant: minimum over the species orderings that sort a permutation-invariant signature).

## 1. Families enumerated

Scripts in `scripts/`, logs in `logs/`. Records (all objects of N and N*) are pickled to the WSL directory
`~/gacE-data/ (WSL)` (outside the repository; 22 files, 26 MB in total; SHA-256 in `logs/data_sha256.txt`).

**Weakly reversible families** (`enum_wr.py`; generator: disjoint unions of strongly connected linkage classes on
distinct complexes, every species occurring; (n, m, R) = species, maximal molecularity, maximal number of
reactions). Cross-check: on the five families of result 17 the generator reproduces exactly the numbers of WR
networks with H found there by an independent subset enumeration (16, 112, 123, 84, 4; `facts.py` F0).

| family | WR networks | s = 1 | s = 2 | s = 3 | s = 4 | with H | H, s=2 | H, s=3 | H, s=4 | H & def 0 | H & WR dual | runtime |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (2,3,6) | 91 585 | 600 | 90 985 | – | – | 460 | 0 | – | – | 12 | 14 | 54 s |
| (3,2,6) | 29 120 | 16 | 3 389 | 25 715 | – | 2 585 | 2 573 | 0 | – | 129 | 10 | 25 s |
| (3,3,5) | 139 339 | 140 | 20 939 | 118 260 | – | 10 560 | 10 467 | 0 | – | 1 694 | 24 | 77 s |
| (4,2,6) | 95 087 | 6 | 1 288 | 37 887 | 55 906 | 19 984 | 684 | 19 300 | 0 | 7 084 | 24 | 146 s |
| (5,2,5) | 4 272 | 0 | 101 | 2 705 | 1 466 | 1 564 | 22 | 1 047 | 495 | 1 535 | 15 | 15 s |
| (5,2,6) homogeneous | 25 461 | 1 | 195 | 8 465 | 16 800 | 21 629 | 63 | 4 766 | 16 800 | 15 949 | 17 | 229 s |
| (6,2,5) homogeneous | 740 | 0 | 17 | 455 | 268 | 251 | 1 | 123 | 127 | 251 | 7 | 34 s |
| (2,2,5) with zero complex | 996 | 23 | 973 | – | – | 16 | 0 | – | – | 4 | 3 | 1 s |
| (3,2,5) with zero complex | 4 574 | 8 | 927 | 3 639 | – | 477 | 472 | 0 | – | 120 | 6 | 4 s |
| (4,3,4) | 14 480 | 48 | 6 191 | 8 241 | – | 6 139 | 3 373 | 2 765 | – | 6 032 | 49 | 7 s |
| (5,3,4) | 15 332 | 13 | 5 965 | 9 354 | – | 6 129 | 2 245 | 3 884 | – | 6 124 | 46 | 32 s |

"Homogeneous": only complexes of molecularity exactly 2 (all such networks are conservative). With H, s ≤ n − 1,
so intermediate dimensions s = 3 and s = 4 under (H1) need n ≥ 4 and n ≥ 5.

**Seeded random sample** beyond exhaustive reach (`sample_wr.py`, not uniform, deterministic): n = 5, molecularity
≤ 2, 4 ≤ r ≤ 8, conservative, s ∈ {2,3,4}, seed 1: 20 000 distinct WR networks (s = 2: 150, s = 3: 7 673, s = 4:
12 177), of which 17 776 satisfy H (42 / 5 557 / 12 177), 13 319 of these with deficiency 0, 14 with a WR dual.
n = 6, molecularity ≤ 2, 5 ≤ r ≤ 8, seed 2: 20 000 networks (s = 2: 28, s = 3: 4 137, s = 4: 15 835), 13 918
with H, 1 with WR dual among those with s = 4.

**All networks with (H1), (H2)** (`enum_all.py`, subsets of distinct reactions, as `enumerate_small.py` of result 17):
(2,2,5) 100, (2,2,6) 116, (2,3,4) 2 784, (3,2,4) 3 883, (3,3,3) 11 409, (4,1,12) 84, (4,2,3) 151 networks. On the
five families of result 17 the counts (H, H+WR, equal rows) agree exactly with its logs (PASS C2).
**Ordinary EC networks** (`enum_fc.py`): n = 3, r ≤ 6: 8 networks with (H1), (H2).

## 2. The coordinator's facts, re-derived (`facts.py`, log `logs/facts.log`, 316 s)

- **Fact 1** (dimension) — PASS on all 769 316 networks of Part A (every reaction set, no hypothesis, families
  (2,2,≤5), (2,3,≤4), (3,1,≤6), (3,2,≤4), (4,1,≤6), (4,2,≤3)) and on every WR record (`enum_wr.py` C3): N* has r
  species, n reactions, Γ(N*) = Γᵀ, the same rank, N** = N, dim ker Γ(N*) = n − s, (H1)(N) ⇔ (H2)(N*),
  (H2)(N) ⇔ (H1)(N*), (H3)(N*) ⇔ columns of Γ distinct. Proof: transposition (manuscript Lemma 1).
- **Fact 5** (siphon dictionary) — PASS on the same 769 316 networks, computed independently on N and on N*:
  siphons of N = T-traps of N*, traps of N = T-siphons of N*, semiflow supports of N = internal cycles of N*.
  Proof: species j of N* is produced (consumed) by reaction i of N* iff (Yp)_ij > 0 ((Ye)_ij > 0), which are the
  two sides of the siphon implication; the cones coincide since Γ(N*) = Γᵀ.
- **Fact 2** (WR not preserved) — the coordinator's numbers 32 / 9 / 24 / 8 are reproduced only if the three small
  families are added with multiplicity: the 3-cycle T belongs to two of them. **Corrected counts:** 31 distinct
  conservative WR networks, 8 with a WR dual; 23 of deficiency 0, 7 with a WR dual (29 / 7 / 21 / 6 with (H3)).
  Smallest conservative WR network with non-WR dual: 2X₂ ⇄ X₁ (dual k₂ → k₁, 2k₁ → 2k₂). On the larger families
  (H): (2,3,6) 14 of 460 have a WR dual; (3,2,6) 10 of 2 585; (3,3,5) 24 of 10 560; (4,2,6) 24 of 19 984;
  (5,2,5) 15 of 1 564; (5,2,6) homogeneous 17 of 21 629; at deficiency 0: 4/12, 4/129, 9/1 694, 9/7 084,
  15/1 535, 9/15 949. A WR dual is rare and becomes rarer as s grows (s = 3: 5 of 19 300 for n = 4; s = 4: 3 of
  16 800 for homogeneous n = 5).
- **Fact 6** (critical siphons) — corrected: of the 31 distinct conservative WR networks, 9 have a critical minimal
  siphon and 22 have minimal siphons = semiflow supports (coordinator: 9 / 23, again counting T twice). The two
  conditions are complementary in general (**proved**: no critical siphon ⇒ every minimal siphon contains, hence
  equals, a semiflow support; and every extreme semiflow support F is a siphon whose sub-siphons contain a
  semiflow support, which by minimality of extreme supports is F). With H: (2,3,6) 38 of 460 have a critical
  siphon; (3,2,6) 1 447 of 2 585; (3,3,5) 9 322 of 10 560; (4,2,6) 17 634 of 19 984 (s = 3: 17 487 of 19 300);
  (5,2,5) 1 224 of 1 564 (s = 4: 487 of 495); (5,2,6) homogeneous 20 633 of 21 629 (s = 4: 16 778 of 16 800).
  **In intermediate dimensions almost every WR network with H has a critical siphon**, so the ADLS criterion
  almost never applies there; the siphon route to persistence needs finer criteria (Anderson–Shiu faces, Wiuf).

## 3. Generalised readings of "refine" (task 3 of the brief)

Counts and smallest counterexamples for every item are in the table of §9 (ids E-…).

**(a) Refinement outside an exceptional class.** For each of the four refinements the smallest exceptional class
that makes it exact is the class where the criteria of Theorem 7 / Corollary 8 of the manuscript fail (they are
equivalences, re-confirmed on the old families by result 17). In the GAC setting the exceptional classes have a
graph description:

- **Proposition E.1 (proved; E-a1).** Let N be weakly reversible. Then P(N) ⪯ M(N*) if and only if no complex is
  the source of two reactions, i.e. every linkage class is a directed cycle; P(N) is then discrete.
  *Proof.* If no complex has two outgoing reactions, the graph G of Definition 19 has no edge, so nothing is
  adjoined and P(N) is discrete, which refines everything. Conversely, manuscript Corollary 10(1). ∎
- **Proposition E.2 (proved; E-a2).** Let N be WR with (H1)–(H3). Then P(N*) is discrete, P(N*) ⪯ M(N) always, and
  M(N) ⪯ P(N*) iff M(N) is discrete iff im Γ contains no nonzero vector supported on two species. *Proof.*
  Manuscript Corollary 9; a vector of im Γ vanishing outside {i, j} exists iff s − rank Γ_{[n]∖{i,j}} ≥ 1, which is
  the test used by the script. ∎

So, in the GAC setting, "preclusters refine MIPS" holds exactly outside the branching networks (reaction side)
and trivially (species side), and "MIPS refine preclusters" holds exactly when the MIPS are trivial. Neither has
content. Outside WR the exceptional classes are not small: see E-a5…E-a8 in §9 (e.g. P(N*) ⪯ M(N) fails exactly
when some twin pair is not conservation-coupled).

**(b) Unions.** "Every MIPS is a union of minimal siphons" (E-b1), "every minimal siphon is a union of MIPS"
(E-b2, = K9(i)), "every precluster of the dual is a union of minimal siphons" (E-b3), "every minimal siphon is a
union of preclusters of the dual" (E-b4): all refuted on both sides, already on two species (§9).

**(c) Galois connection.** *Proposition E.3 (proved).* Let X be finite, P a partition of X with Boolean lattice
L(P) of unions of blocks, and F ⊆ 2^X any family closed under unions and containing ∅ (for instance L(Q) for a
second partition Q, or the siphons together with ∅, or the traps together with ∅). Put c_P(U) = union of the
P-blocks meeting U and i_F(V) = union of the members of F contained in V (the largest member of F inside V).
Then for U ∈ F and V ∈ L(P): c_P(U) ⊆ V ⇔ U ⊆ V ⇔ U ⊆ i_F(V). So c_P ⊣ i_F is always a Galois connection
(monotone adjunction) between F and L(P). For F = L(Q) it is a pair of inverse isomorphisms iff P = Q.
*Proof.* The first equivalence holds because V is a union of P-blocks, the second because i_F(V) is the largest
member of F inside V. If the maps are inverse isomorphisms of Boolean lattices, c_P maps atoms to atoms, so every
Q-block lies in a single P-block (Q refines P) and |Q| = |P|; a refinement with the same number of blocks is an
equality. Conversely c_P and i_P are the identity on L(P). ∎
Consequently "there is a Galois connection between preclusters, MIPS and siphons" is true for every network and
carries no information; the information is in the closed elements, e.g. the join P ∨ M (unions of blocks of
both). Tested: E-c1 (WR: P(N) ∨ M(N*) = linkage classes) refuted; E-c1' (the same at deficiency 0) refuted;
E-c2 (the join refines the linkage classes) refuted outside WR; E-c3 (species side: P(N*) ∨ M(N) is one block for
connected N); E-c4 (adding the minimal siphons to the join changes nothing). Explicit illustration
(`e_examples.py` E1): Y5 = A→B, B→C, B→D, B→E, C→A, D→A, E→A (one linkage class, WR, deficiency 0, (H1)–(H3)) has
P(Y5) = {A→B}, {the six others} and M(Y5*) = {A→B}, {B→C, C→A}, {B→D, D→A}, {B→E, E→A}: the join has two blocks.

**(d) Removing singletons.** E-d1 (non-singleton blocks of P(dual) = non-singleton MIPS), E-d2 (non-singleton
blocks of P(dual) refine M), E-d3 (non-singleton MIPS lie in blocks of P(dual)), E-d4 (non-singleton minimal
siphons = non-singleton blocks of P(dual)): all refuted. What survives is agent P's K13 = P-H (**proved** by P,
checked here on both sides): every non-singleton precluster of N* is a union of MIPS of N (and dually under (H2)).
Its proof is one line from Definition 19: if t is in a doubling pair or adjoined through one, every s with
n_s ∈ ℝ_{>0} n_t lies in the same span and is adjoined too.

## 4. Red team of agent G (`scripts/g_examples.py`, log `logs/g_examples.log`; tests G-… and X-… in §9)

All explicit claims of `../gac/notes.md` that were checked hold exactly:

- **§3, extinction example (coordinator's point 1). Confirmed, and strengthened.** M = {3A → 2A+B, 2A+B → 3A},
  D = M* = {3U+2V → 2U+3V, V → U}: both Γ = [[−1,1],[1,−1]], both satisfy (H1)–(H3); M reversible, one linkage
  class, deficiency 0; D has 4 complexes, 2 linkage classes, deficiency 1, is not WR, and w = (0,1) witnesses that
  it is not endotactic. With unit rates v' = u³v² − v (exact), h(v) = (2−v)³v has h' = (2−v)²(2−4v) (polynomial
  identity), h(1/8) = 3375/4096, h(1) = 1, so v' ≤ −(721/4096)v on (0, 1/8] and v → 0 from (15/8, 1/8) (G3a–b PASS).
  The class u + v = 2 has exactly two positive equilibria (v = 1 and one root of h = 1 in (1/8, 1/2)), so the
  trajectory starts below the saddle-like interior equilibrium; there is no contradiction with the GAC since D is
  not WR. **Strengthening (E, proved):** D is non-persistent for *every* choice of rate constants and in *every*
  class: v' = v(k₁u³v − k₂) ≤ −(k₂/2)v whenever 0 < v ≤ v* = min(T, k₂/(2k₁T³)) (u ≤ T), so v → 0 from every
  point with v < v*. This is an instance of the general rate-independent certificate of §5 (V is a "drain
  species" of D). Since M is reversible with one linkage class it is globally attracting for all rate constants
  (Anderson 2011; here also by the explicit 1-D argument), so **persistence is not invariant under N ↦ N*, for any
  rate constants on the non-WR side**, within (H1)–(H3) on both sides.
- **§5, ADLS not self-dual.** H = {2A ⇄ 2B, 2A ⇄ A+B}: WR, (H1)–(H3), only minimal siphon {A,B} (ADLS holds); H* has
  (H1)–(H3) and the critical siphon {U4} (row (1,−1)) (G5 PASS).
- **§10, finite family.** The 705 labelled networks (575 unions of reversible pairs, 130 cycles on the six
  degree-two complexes of three species): 702 with (H3); ADLS on both sides / only N / only N* / neither =
  18 / 177 / 0 / 507; 13 WR duals; q = β* + δ* on all (G10 PASS). On the larger families "only N*" does occur
  (E-e3, §9: 924 of the 82 377 WR networks with H): smallest WR network with (H1)–(H3), N* also with (H1)–(H3),
  where N* satisfies ADLS but N does not: X₄ ⇄ X₂+X₃, 2X₄ ⇄ X₁+X₃ (n = 4, s = 2, deficiency 0).
- **Prop. 2.1** (q = β* + δ*; ℓ* + δ* ≤ q when N* is WR; q = 1 and N* WR ⇒ one directed cycle) and **Prop. 2.2**
  (WR + (H3): N* WR iff the sets of source rows and product rows coincide; m* = n, ℓ* = c(π), δ* = q − c(π),
  π-orbits give conservation laws; q = 1 ⇒ one n-cycle): no counterexample (G-2.1, G-2.2 in §9).
- **Prop. 4.1** (N, N* WR, (H3), no zero reaction, q = 1, n prime ⇒ N strongly endotactic): the argument was
  re-read and is correct (π is a single n-cycle by Prop. 2.2, each reaction maps y to Ty with T the cyclic shift,
  linkage classes are T-orbits, Φ_n irreducible makes every nonzero d ∈ V a cyclic vector, so each linkage class
  spans S = V, and a w not orthogonal to S is non-constant on every linkage class). Checked: the 5-species
  example (10 reactions, s = 4, δ = 4, dual one 5-cycle with δ* = 0, critical siphon {X1,X2,X3}, strongly
  endotactic by the exact test) and the 4-species example showing primality cannot be dropped (not strongly
  endotactic, witness w = (3,2,1,2)) (G4.1a–f PASS). On the enumerated families (X-G4.1, `logs/tests_X2.log`):
  all 27 networks satisfying the hypotheses (n = 2, 3, 5) are strongly endotactic; without primality (n = 4)
  4 of 5 are, the exception being G's own example A→B→C→D→A, A+C ⇄ B+D.
- **Example 8.1**: Γ, rank 2, deficiency 0, WR dual with deficiency 0, J_S = [[−15,27],[24,−25]], det −273,
  and the alternative realisation 10A+2C → 9A+3B, 2A+B → 3C, 10B+C → 3A+8B with the same Γ (G8.1a–c PASS).
- **Prop. 7.1** (WR + (H1) + rows of Ye with sum ≤ 1 ⇒ disjoint cycles, deficiency 0, ADLS; with (H3) every
  complex is one species): no counterexample (G-7.1).

## 5. Rate-independent non-persistence and other transfers (coordinator's priority)

**Proposition E.4 (drain species; proved).** Let X be conservative and i a species such that (a) some reaction
k₀ has source d₀e_i (only species i, d₀ ≥ 1) and consumes X_i, and (b) every reaction producing X_i has a source
with X_i-coefficient > d₀. Then for every choice of rate constants and every compatibility class, every positive
solution with x_i(0) small enough has x_i(t) → 0. In particular X is not persistent, whatever the rates.
*Proof.* On a class all coordinates are bounded by some T (H1). Production terms are κ_k x^{y_k}Δ_{k,i} with
y_{k,i} ≥ d₀ + 1, hence ≤ C x_i^{d₀+1}; consumption terms are ≤ 0 and include −κ_{k₀}|Δ_{k₀,i}| x_i^{d₀}. So
x_i' ≤ x_i^{d₀}(C x_i − c₀) with c₀ > 0. Below c₀/(2C), x_i decreases, converges to some L ≥ 0, and L > 0 would
give x_i' ≤ −(c₀/2)L^{d₀} forever. ∎
`endo.drain_species` implements (a)–(b). For a WR network the condition never holds on N itself (X-drainN, a
consistency check: otherwise the GAC would fail in its proved cases). Counted on N* for WR networks with
(H1)–(H3): X-drainD (all), X-drainD1 (one linkage class and s ≤ 2, where global attraction of N for
complex-balanced rates is classical: Anderson 2011, Pantea 2012). Each witness is a network N that is globally
attracting for complex-balanced rates while N* (conservative, consistent) loses a species for every rate choice.

**What the duality transfers, counted on the 82 377 WR networks with (H1)–(H3)** (all families and samples;
s = 2 / 3 / 4 in parentheses; ids in §9):

| property of N (WR + H) | holds for N* | ids |
|---|---|---|
| (H1)–(H3) | 72 342 (9 104 / 25 371 / 37 783 of 13 337 / 29 419 / 39 061); fails only through (H3), i.e. two reactions with the same reaction vector | E-e1 |
| weak reversibility | 98 (40 / 20 / 5) | E-e4, E-e5 applicability |
| deficiency 0 (46 709 networks) | 43 (17 / 14 / 5 of 1 802 / 14 090 / 30 797) | E-e6 |
| ADLS (no critical siphon) | fails for N* in 5 252 networks where it holds for N; holds for N* but not for N in 924 | E-e2, E-e3 |
| ADLS, N and N* both WR + H (87 networks) | differs in 8 | E-e4b |
| endotactic (X-endoD) | 122 of 84 103 (XA 92 / 34 927, XB 22 / 35 258, n = 6 sample 8 / 13 918) |  X-endoD |
| strongly endotactic (X-sendoD) | 43 of 70 185 (XA 39, XB 4) | X-sendoD |
| no drain species, i.e. no rate-independent extinction (X-drainD) | fails for 141 (all in the molecularity-3 families, s ≤ 2; 106 of them with one linkage class and s ≤ 2) | X-drainD, X-drainD1 |

So within the GAC setting the dual almost never stays in the GAC setting: weak reversibility survives in about
0.1 % of the networks (0.01 % at s = 4), deficiency zero in 0.09 %, and the siphon criterion is not invariant
in either direction, even when both sides are weakly reversible.

**Endotacticity** (exact Fourier–Motzkin decision of the definitions of Craciun–Nazarov–Pantea and
Gopalkrishnan–Miller–Shiu, `scripts/endo.py`). Validation: "WR ⇒ endotactic" holds on all 60 159 WR networks with
H and r ≤ 6 (X-endoN), "WR with one linkage class ⇒ strongly endotactic" on all 40 804 (X-sendoN); the code also
returns the expected answers on M, D, T, 2A → A, A → 2A and two reversible pairs. Result: the dual of a WR network
with (H1)–(H3) is almost never endotactic (122 of 84 103) and almost never strongly endotactic (43 of 70 185);
smallest case N = 2X₂ ⇄ X₁, N* = k₂ → k₁, 2k₁ → 2k₂. So none of the endotactic persistence theorems
(Craciun–Nazarov–Pantea, Gopalkrishnan–Miller–Shiu, the permanence results) can be applied to N* in general.

**Rate-independent extinction in the dual (X-drainD).** 141 WR networks with (H1)–(H3), all of molecularity 3
and s ≤ 2 (none in the molecularity-≤ 2 families, none with s ≥ 3), have a dual with a drain species; 106 of
them have one linkage class, and 45 of the 141 have deficiency 0. Smallest: G's M (3X₂ ⇄ X₁+2X₂); smallest with
three species (`e_examples.py` E3): N = 3X₃ → X₂ → X₁+2X₃ → 3X₃ (one cycle, deficiency 0, s = 2, (H1)–(H3), so
complex balanced and globally attracting for every choice of rate constants by Anderson 2011), whose dual
N* = k₃ → k₂, k₂ → k₁, 3k₁+2k₃ → 2k₂+3k₃ satisfies (H1)–(H3) and loses k₃ for every choice of rate constants.

**In intermediate dimension (molecularity-3 families (4,3,4), (5,3,4), `logs/tests_m3.log`, `logs/tests_m3X.log`;
12 268 WR networks with H).** 162 duals have a drain species, 103 of them with **s = 3** (all of these with one
linkage class and deficiency 0, and N* with (H1)–(H3)). Smallest s = 3 witness (`e_examples.py` E3):
N = X₄ → 3X₃ → X₁+2X₂ → 3X₂ → X₄ (one 4-cycle, n = 4, s = 3, deficiency 0, (H1)–(H3)), globally attracting for
all rate constants (Anderson 2011; also Pantea 2012 since s = 3), with dual N* = k₄ → k₂, 3k₃+2k₄ → 2k₂+3k₄,
3k₂ → 3k₁, k₁ → k₃ ((H1)–(H3)), in which k₄ is a drain species: N* is non-persistent for every choice of rate
constants. On these families too the dual is rarely endotactic (102 of 12 268) or strongly endotactic (16), and
every claim marked proved holds (K2n, K14a–b, and all of §9 rerun: PASS); E-e4b holds on the 95 networks of these
families with both sides WR + H.

## 6. Agent P's candidates (`../petri/candidates.md`, K1–K14; ids P-… in §9)

Proved claims re-tested (no counterexample unless stated): **K1** transversality (both sides, no hypothesis),
**K2** exact criterion for "every precluster lies in a minimal siphon" (also in P's latest hypothesis-free form,
P-K2n, on 869 968 sides, with twins read as crn.py does, i.e. including equal zero rows of Ye), **K3(a)–(g)** structure of WF-FC networks
(both sides), **K10** absolute consistency of EC networks (k = 1 and three seeded random rational k), **K11 / T\***
(for balanced EC networks: N has positive steady states for all tested k ⇔ N* for all tested κ ⇔ rank equation;
the dual side is decided exactly, because for EC with balanced joins the steady-state condition of N* is linear in
one positive unknown per source complex), **K13 = H** (non-singleton preclusters of the dual are unions of MIPS),
**C, C′, K14** (well-formed EC: no critical siphon on either side; minimal siphons = semiflow supports; minimal
T-traps = internal cycles), **G** ((iii) ⇒ (i), (ii)).

- **K7: one red-team finding (the only FAIL of the whole run).** The last clause ("with (H3) and ordinary, every
  complex is a single species") fails for weakly reversible EC networks containing the **zero complex**: 26
  networks, all in the zero-complex families, smallest 0 → X₂ → X₁ → 0 (WR, EC, ordinary, (H3), deficiency 0, not
  conservative). All other clauses hold on all 1 581 WR EC networks tested (incl. 26 + 174 with the zero complex),
  and the whole of K7 holds on the 556 with (H1)–(H3). Repair: add (H1) (or "no zero complex") to the last clause.
  P restricted K7 accordingly (18:15); the restricted version P-K7r holds on all 1 683 WR EC networks without the
  zero complex in all record files, molecularity-3 families included (704 with H; `logs/tests_K7.log`).
- **K3** holds on all 492 WF-FC sides (224 with H; 186 WR). **K3-need** (FC + (H1) + (H2) + connected, rank
  equation violated): (b)–(e) fail on all 76 such sides; smallest ordinary one: k₃ → k₂, k₃ → k₁, k₁+k₂ → k₂+k₃,
  k₁+k₂ → k₁+k₃ (the dual of X₃+X₄ → X₂+X₄, X₃+X₄ → X₁+X₃, X₁+X₂ → X₃+X₄). The rank equation is needed.
- **K4 (question)** well-formed EC ⇒ semiflow supports = minimal siphons: holds on all 4 632 sides with balanced
  joins and all 950 without.
- **K5 (conjecture, now proved by P)** s ≥ m − 1 for EC, connected, (H1), (H2): holds on all 6 208 sides.
- **K10 / K11** hold on all 6 208 / 4 982 sides; **C, C′, K14a** on all 5 582 well-formed EC sides, **K14b**
  (pure) on all 3 466; **K14c** (impure well-formed EC, minimal siphons = semiflow supports) holds on all 2 116:
  purity was not needed on the tested families.
- **K6 (question)** "no critical siphon" is not self-dual: refuted (smallest: X₂+X₃ → X₁, X₁+X₃ → X₂+2X₃, with
  (H1), (H2)); among the 82 377 WR networks with H the two sides differ for 6 176 (s = 2: 2 325 of 13 337; s = 3:
  2 740 of 29 419; s = 4: 627 of 39 061). Both directions fail (E-e2: ADLS(N) but not ADLS(N*) for 5 252, smallest
  2X₂ → X₁+X₂ → 2X₁ → 2X₂; E-e3: the converse for 924, smallest X₄ ⇄ X₂+X₃, 2X₄ ⇄ X₁+X₃), and even when **both**
  N and N* are WR with (H1)–(H3) the criterion differs for 8 of 87 networks (E-e4b; smallest
  N = 2X₂ ⇄ 2X₁, X₁+2X₂ ⇄ 2X₁+X₂, with N* = k₂+2k₃+2k₄ ⇄ 2k₁+2k₂+k₄ having critical siphons {k₂}, {k₄}).
- **K8 (obstruction)** found (`scripts/k8.py`, log `logs/k8.log`). Species side: N₁ = X₂+X₃ → 2X₃, X₂+X₃ → 2X₂,
  X₁+X₂+X₃ → X₁+2X₃ and N₂ = X₂+X₃ → 2X₃, X₁+X₂+X₃ → 2X₂+X₃, X₁+X₂+X₃ → 2X₁+X₂ both satisfy (H1)–(H3), have
  the same twin classes {X₁}, {X₂,X₃} and the same minimal siphons {X₁}, {X₂}, {X₃}, but P(N₁*) = {X₁},{X₂,X₃}
  and P(N₂*) = {X₁,X₂,X₃}. Reaction side: 169 colliding keys among the (H1)–(H3) networks with r ≤ 6. No witness
  pair consists of two WF-FC networks (consistent with K3(e)).
- **K9** (all four statements under (H1)–(H3) on the side's network): fail in general, smallest counterexample
  the dual side of X₂ → X₁, X₁+X₂ → 2X₂ (and for (i)–(iii) on EC: the dual side of X₁+X₂ → 2X₂, X₁+X₂ → 2X₁; on
  FC: the dual side of X₃+X₄ → X₂+X₄, X₃+X₄ → X₁+X₃, X₁+X₂ → X₃+X₄). Holds / applicable (sides):
  (i) 38 443 / 187 195, (ii) 56 854 / 187 195, (iii) 21 083 / 187 195, (iv) 84 767 / 187 195 under H; WR:
  30 264, 49 268, 15 051, 68 190 of 154 719; on EC sides 4 839, 4 895, 4 812, **5 429 of 5 429**; on FC sides
  372, 402, 370, **444 of 444**; on WF-FC sides all four hold (358 / 358). So **K9(iv) holds on every EC
  network tested** (P's observation confirmed on 5 429 sides), while (i)–(iii) need well-formedness.
- **K12 (question)** not settled. Remark (proved): for a WR network with (H1) the source complexes span a space of
  dimension ≥ s + 1 (if all complexes were in S, the positive conservation law would vanish on a nonzero
  nonnegative complex), so P's cheap necessary condition rank Ye − δ ≥ s holds automatically for N*; it cannot
  produce a WR example.

## 7. Scripts and logs

All scripts: Python 3.12 under WSL Fedora 44, standard library only, exact rationals, print PASS/FAIL per claim
and their runtime, exit non-zero on a FAIL of a claim marked proved. Run from `computation/` as
`python scripts/<name>.py …` (record files in `~/gacE-data/ (WSL)`, hashes in `logs/data_sha256.txt`).

| script | checks | log |
|---|---|---|
| `ecore.py` | library: cached network on top of `crn.py`, canonical form, WR generator, siphons/traps, Petri clusters, free-choice flags, the record of N and N* | – |
| `enum_wr.py` | WR families (C1 WR, C2 distinct orbits, C3 Lemma 1 on every network), tables by s | `enum_wr_*.log` |
| `enum_all.py` | all networks with (H1),(H2); C2 = agreement with result 17's logs | `enum_all_*.log` |
| `enum_fc.py` | ordinary EC networks with (H1),(H2) | `enum_fc_*.log` |
| `sample_wr.py` | seeded samples n = 5, 6 (conservative WR, s ∈ {2,3,4}) | `sample_*.log` |
| `facts.py` | coordinator's facts 1, 2, 5, 6; F0 WR generator vs result 17 | `facts.log` |
| `tests.py` | all candidates (E-…, P-…, G-…, X-…) on all record files; counts per class, smallest counterexample | `tests_final.log`, `tests_K14.log`, `tests_K2n.log`, `tests_e4b.log`, `tests_X2.log`, `tests_XA.log`, `tests_XB.log`, `tests_X6.log`, `tests_m3.log`, `tests_m3X.log` (+ `.json`); `tests_main.log` is an earlier full run superseded by `tests_final.log` |
| `endo.py` | exact endotactic / strongly endotactic / drain-species decisions (Fourier–Motzkin) | via `tests.py` |
| `k8.py` | agent P's K8 search | `k8.log` |
| `g_examples.py` | agent G's explicit examples and 705-network counts | `g_examples.log` |
| `e_examples.py` | E's explicit examples (Y5; D non-persistent for all rates; drain witnesses) | `e_examples.log` |

## 8. What is left

- K12 (absolute consistency of N* for EC networks with unbalanced joins, and for WR networks) is not decided;
  for unbalanced EC the steady-state system of N* is polynomial in one unknown per source complex and would need
  a real-root certificate.
- The exhaustive WR families stop at r ≤ 6 (n ≤ 4), r ≤ 5 (n = 5, 6) and homogeneous r ≤ 6 (n = 5); beyond, only
  seeded non-uniform samples. Linkage classes with more than 7 complexes and 8 reactions are not generated.
- Strong endotacticity of N* was not computed on the n = 6 sample (endotacticity and drain species were).
- Exhaustive molecularity-3 families stop at r ≤ 4 for n = 4, 5; s = 4 with molecularity 3 is not covered.
- Persistence itself is not decided by any script: the counts concern sufficient criteria (ADLS, drain species).
- G's Proposition 6.2 (face codimension = rank of a dual sub-network) and 6.3 (Wiuf translation) were read but not
  tested computationally.

## 9. Tables (from `logs/tests_final.log`, `logs/tests_K14.log`, 434 984 distinct networks, 523 s + 307 s)

Entries are "holds / applicable". Columns: all networks of the 20 record files (subject to the candidate's own
hypothesis); H; strong; WR (= H + WR); WR0; WR with s = 2, 3, 4. Candidates tested on both sides count each side
separately. 'proved' = claimed proved by its author (or by E): PASS = no counterexample; 'test' = conjecture,
question or reading: REFUTED = counterexample found (smallest one in the log). Statements: `tests.py` (each
candidate's docstring is printed in the log with its smallest counterexample).

```
candidate      status              all               H          strong              WR             WR0          WR,s=2          WR,s=3          WR,s=4  verdict
E-a1           proved   417174/417174    82377/82377     19496/19496     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    PASS
E-a2           proved    82377/82377     82377/82377     19496/19496     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    PASS
E-a3           test       6876/100084     6876/100084     1842/30344      2751/82377         0/46709      2077/13337       403/29419        96/39061    REFUTED
E-a4           test        798/100084      798/100084      336/30344       365/82377       362/46709        20/13337       305/29419        37/39061    REFUTED
E-a5           test      99932/100084    99932/100084    30341/30344     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    REFUTED
E-a6           test        950/100084      950/100084      339/30344       365/82377       362/46709        20/13337       305/29419        37/39061    REFUTED
E-a7           test      43446/100084    43446/100084    17774/30344     28303/82377     11755/46709      5784/13337      7602/29419     14722/39061    REFUTED
E-a8           test      42276/100084    42276/100084    13406/30344     35587/82377     20195/46709      9011/13337     11427/29419     14609/39061    REFUTED
E-b1           test      61855/223668    60201/200168    21888/60688     44960/164754    26925/93418      7465/26674     15192/58838     21776/78122    REFUTED
E-b2           test      54054/223668    44796/200168     7595/60688     33954/164754    13033/93418     11274/26674     16141/58838      5471/78122    REFUTED
E-b3           test      17899/223668    11817/200168     5687/60688      9660/164754     9591/93418       392/26674      1821/58838      7445/78122    REFUTED
E-b4           test     146181/223668   142895/200168    47182/60688    111379/164754    58494/93418     19468/26674     37340/58838     53810/78122    REFUTED
E-c1           test     180011/417174    38557/82377     11010/19496     38557/82377     38557/46709      1449/13337     10082/29419     27006/39061    REFUTED
E-c1'          test     180011/196870    38557/46709     11010/11958     38557/46709     38557/46709      1449/1802      10082/14090     27006/30797    REFUTED
E-c2           test      89263/100084    89263/100084    21173/30344     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    REFUTED
E-c3           test      75501/99354     75501/99354     28033/30328     59065/81647     28581/46152     12568/13278     19362/28983     26675/38826    REFUTED
E-c4           test     141263/223668   126432/200168    48985/60688     95159/164754    59740/93418     15828/26674     30324/58838     48392/78122    REFUTED
E-d1           test       7829/223668     7674/200168     2178/60688      3116/164754      362/93418      2097/26674       708/58838       133/78122    REFUTED
E-d2           test     146580/223668   143378/200168    48115/60688    110680/164754    58464/93418     19121/26674     37021/58838     53783/78122    REFUTED
E-d3           test      58105/223668    43226/200168    13745/60688     35952/164754    20557/93418      9031/26674     11732/58838     14646/78122    REFUTED
E-d4           test      12455/223668    12370/200168     6377/60688      3409/164754     1497/93418      1952/26674       656/58838       736/78122    REFUTED
E-e1           test      72342/82377     72342/82377     19496/19496     72342/82377     46709/46709      9104/13337     25371/29419     37783/39061    REFUTED
E-e2           test      77125/82377     77125/82377     19102/19496     77125/82377     46709/46709     11016/13337     27244/29419     38789/39061    REFUTED
E-e3           test      81453/82377     81453/82377     19453/19496     81453/82377     45793/46709     13333/13337     28854/29419     38706/39061    REFUTED
E-e4           test         85/98           85/98           12/12           85/98           45/45           37/40           19/20            5/5        REFUTED
E-e5           test         87/98           87/98           12/12           87/98           45/45           39/40           20/20            5/5        REFUTED
E-e6           test         43/46709        43/46709        12/11958        43/46709        43/46709        17/1802         14/14090         5/30797    REFUTED
E-e7           test       4388/7537       4388/7537        417/811        4388/7537       4127/4127        349/1730       2800/4300       1214/1434     REFUTED
P-K1           proved   869968/869968   200168/200168    60688/60688    164754/164754    93418/93418     26674/26674     58838/58838     78122/78122    PASS
P-K2           proved   223668/223668   200168/200168    60688/60688    164754/164754    93418/93418     26674/26674     58838/58838     78122/78122    PASS
P-K3a          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3b          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3c          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3d          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3e          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3f          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3g          proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3g'         proved      492/492         224/224          12/12          186/186         186/186          10/10          166/166           8/8        PASS
P-K3-need      test          0/76            0/72            0/26            0/0             0/0             0/0             0/0             0/0        REFUTED
P-K4i          test       4632/4632       3124/3124       1714/1714       1008/1008       1008/1008        174/174         708/708         114/114      holds on all tested
P-K4ii         test        950/950         770/770          40/40          236/236         164/164          67/67          124/124          38/38       holds on all tested
P-K5           test       6208/6208       4512/4512       2183/2183       1244/1244       1172/1172        241/241         832/832         152/152      holds on all tested
P-K6           test     101986/111834    91709/100084    29765/30344     76201/82377     45793/46709     11012/13337     26679/29419     38434/39061    REFUTED
P-K7           proved     1555/1581        556/556          98/98          556/556         556/556          93/93          388/388          66/66       FAIL
P-K9i          test      43334/197776    38443/187195     7595/60688     30264/154719    13033/93418      9092/22441     15187/54790      5393/76844    REFUTED
P-K9ii         test      62609/197776    56854/187195    10592/60688     49268/154719    28608/93418      8068/22441     22581/54790     18031/76844    REFUTED
P-K9iii        test      23655/197776    21083/187195     2999/60688     15051/154719     9170/93418      3076/22441      8378/54790      3057/76844    REFUTED
P-K9iv         test      89604/197776    84767/187195    27410/60688     68190/154719    45038/93418      8994/22441     24555/54790     34099/76844    REFUTED
P-H            proved   546818/546818   200168/200168    60688/60688    164754/164754    93418/93418     26674/26674     58838/58838     78122/78122    PASS
P-C            proved     5582/5582       3894/3894       1754/1754       1244/1244       1172/1172        241/241         832/832         152/152      PASS
P-C'           proved     5582/5582       3894/3894       1754/1754       1244/1244       1172/1172        241/241         832/832         152/152      PASS
P-G            proved   197776/197776   187195/187195    60688/60688    154719/154719    93418/93418     22441/22441     54790/54790     76844/76844    PASS
P-K10          proved     6208/6208       4512/4512       2183/2183       1244/1244       1172/1172        241/241         832/832         152/152      PASS
P-K11          proved     4982/4982       3466/3466       1898/1898       1008/1008       1008/1008        174/174         708/708         114/114      PASS
P-K9i-EC       test       4839/5429       3826/4412       1771/2183       1348/1348       1273/1273        256/256         901/901         172/172      REFUTED
P-K9i-FC       test        372/444         232/302          12/38          192/192         192/192          12/12          168/168          10/10       REFUTED
P-K9i-WFFC     test        358/358         224/224          12/12          186/186         186/186          10/10          166/166           8/8        holds on all tested
P-K9ii-EC      test       4895/5429       3880/4412       1777/2183       1348/1348       1273/1273        256/256         901/901         172/172      REFUTED
P-K9ii-FC      test        402/444         260/302          18/38          192/192         192/192          12/12          168/168          10/10       REFUTED
P-K9ii-WFFC    test        358/358         224/224          12/12          186/186         186/186          10/10          166/166           8/8        holds on all tested
P-K9iii-EC     test       4812/5429       3799/4412       1754/2183       1348/1348       1273/1273        256/256         901/901         172/172      REFUTED
P-K9iii-FC     test        370/444         230/302          12/38          192/192         192/192          12/12          168/168          10/10       REFUTED
P-K9iii-WFFC   test        358/358         224/224          12/12          186/186         186/186          10/10          166/166           8/8        holds on all tested
P-K9iv-EC      test       5429/5429       4412/4412       2183/2183       1348/1348       1273/1273        256/256         901/901         172/172      holds on all tested
P-K9iv-FC      test        444/444         302/302          38/38          192/192         192/192          12/12          168/168          10/10       holds on all tested
P-K9iv-WFFC    test        358/358         224/224          12/12          186/186         186/186          10/10          166/166           8/8        holds on all tested
G-2.1          proved   434984/434984   100084/100084    30344/30344     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    PASS
G-2.2          proved   398807/398807    82377/82377     19496/19496     82377/82377     46709/46709     13337/13337     29419/29419     39061/39061    PASS
G-7.1          proved       18/18            6/6             3/3             6/6             6/6             2/2             2/2             1/1        PASS
G-req1         test      77125/82377     77125/82377     19102/19496     77125/82377     46709/46709     11016/13337     27244/29419     38789/39061    REFUTED
G-req1'        test      69193/72342     69193/72342     19102/19496     69193/72342     46709/46709      7723/9104      23871/25371     37563/37783    REFUTED
P-K14a         proved     5582/5582       3894/3894       1754/1754       1244/1244       1172/1172        241/241         832/832         152/152      PASS
P-K14b         proved     3466/3466       1802/1802        601/601        1244/1244       1172/1172        241/241         832/832         152/152      PASS
P-K14c         test       2116/2116       2092/2092       1153/1153          0/0             0/0             0/0             0/0             0/0        holds on all tested
```

Further runs (same conventions): `logs/tests_K2n.log`, `logs/tests_e4b.log`, `logs/tests_X2.log`,
`logs/tests_X6.log` (n = 6 sample), `logs/tests_XA.log`, `logs/tests_XB.log` (transfers on the WR families):

```
P-K2n          proved   869968/869968   200168/200168    60688/60688    164754/164754    93418/93418     26674/26674     58838/58838     78122/78122    PASS
E-e4b          test         79/87           79/87           12/12           79/87           45/45           36/39           19/20            5/5        REFUTED
X-L5           proved   331653/331653   126998/126998    51648/51648     91584/91584     48606/48606      6862/6862      49542/49542     34220/34220    PASS
X-G4.1         proved       27/27           27/27           10/10           27/27           14/14            9/9             0/0             4/4        PASS
X-G4.1np       test          4/5             4/5             2/2             4/5             2/2             0/0             4/5             0/0        REFUTED
n = 6 sample:
X-endoD        test          8/13918         8/13918         0/2083          8/13918         8/12352         1/1             6/1531          1/12386    REFUTED
X-drainN       proved    20000/20000     13918/13918      2083/2083      13918/13918     12352/12352         1/1          1531/1531      12386/12386    PASS
X-drainD       test      13918/13918     13918/13918      2083/2083      13918/13918     12352/12352         1/1          1531/1531      12386/12386    holds on all tested
```

