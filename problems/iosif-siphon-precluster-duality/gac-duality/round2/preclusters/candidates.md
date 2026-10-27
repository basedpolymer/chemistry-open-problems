# Q2 candidates for E2 (corrected preclusters)

Agent Q2, Claude Opus 5.5 (effort xhigh), 27 October 2026. Definitions in `notes.md` §1; library
`scripts/pclust.py` (imports `crn.py` and `petri.py`; exact). All networks below satisfy (H2) on the side considered
(reaction side: N; species side: N*, i.e. (H1) for N). Notation: 𝒮 = source classes (reactions with the same source
complex; the doubling classes), P = Definition 19, Pcs = agent P's class-span closure, π = flat closure
(`flat_closure`), CF = Conradi–Flockerzi clusters = meet over V ∈ C° of κ_V (`kappa`), M* = full-coupling classes =
MIPS(N*), M = ray MIPS of N. Order: A ⪯ B means A refines B.

How to compute CF exactly: `pi_plus` is a sound lower bound (flat closure + two-block species rule + max-principle
rule), `cf_upper` an upper bound (meet of κ_V at random V and at "split" V); CF is determined when they coincide
(`cf_bounds`). For a one-block upper bound, `one_block_certificate` (sympy) certifies CF = one block when a
(q−1)-minor of the block matrix X(V) is a sign-definite polynomial in the ray weights.

| id | statement | hypotheses | status |
|---|---|---|---|
| Q2-1 | P ⪯ Pcs ⪯ π ⪯ CF ⪯ κ_V for every V ∈ C° | (H2) | proved (notes §2); checked `check_families.py` K1 |
| Q2-2 | π = Pcs ∨ M* (non-singleton Pcs blocks, plus full-coupling classes of the Pcs singletons) | (H2) | proved; K2 |
| Q2-3 | M* ⪯ π; and π = M* ⇔ CF = M* ⇔ κ_V = M* (any V) ⇔ every doubling pair is fully coupled | (H2) | proved (Thm 6.1); K3 |
| Q2-4 | species side: M(N) ⪯ π(N*) ⪯ CF(N*); equality ⇔ every twin pair is conservation-coupled; holds for every N with (H1), (H3) whose product complexes are source complexes (e.g. WR) | (H1) | proved (Cor 6.2); K3 on the dual side |
| Q2-5 | EC + connected + (H1) + (H2) ⇒ κ_V is one block for every V ∈ C°, so CF is one block (π may have several: K3) | | proved (Thm 4.1, from P's Lemma CF); K4 |
| Q2-6 | WR with deficiency 0 ⇒ κ_V = CF = linkage classes (as a partition of the reactions) for every V ∈ C° | | proved (Thm 4.2); K5 |
| Q2-7 | CF is invariant under relabelling, duplication of a reaction and splitting a choice through a new intermediate species (`duplicate`, `split_choice`); P, Pcs, π are not (B3+split for P; N4.12+dup for P, Pcs; B2dup+split of {R1,R2} for Pcs, π) | (H2) | proved (Thm 3.1, Prop 3.2); checked on examples (X7, X9) |
| Q2-8 | For V ∈ C°: dim ker A_V ≤ #κ_V (classes of source complexes), with equality iff ker A_V has a basis with pairwise disjoint supports (Pérez Millán et al. Condition 3.1) iff the positive steady states for k = V are {x : x^{y_c − y_c′} = 1 for c ~ c′} | (H2) | proved (Prop 5.1); K6 (one direction) |
| Q2-9 | Isolation (every internal cycle lies in a CF block) ⇒ κ_V = CF and dim ker A_V = #CF for every V ∈ C° (toric steady states for every k, rate-independent partition) | (H2) | proved (Cor 5.3); K6 |
| Q2-10 | Q3 for corrected preclusters: every block of CF(N*) lies in a minimal siphon of N ⇔ N has no twin pair and every ray-MIPS of N lies in a minimal siphon (dually for CF(N) and minimal T-traps) | (H1) | proved (Prop 6.4); K7 |
| Q2-11 (question) | How often is CF determined by `cf_bounds`, and how often π ≠ CF, κ_gen ≠ CF? Smallest π ≠ CF with (H1)–(H3) and rows/columns pairwise independent? | | first answers (notes §7): CF determined on 63 717 of 63 722 sides (rest: Thm 4.2); π ≠ CF on 100 sides (90 N, 10 N*), smallest n = 2: {2X₂→X₁+X₂, 2X₂→2X₁, X₁+X₂→2X₂, X₁+X₂→2X₁}; κ_gen ≠ CF on 5, smallest XY4. Larger families (samples n = 5, 6; (5,3,4)) not run |
| Q2-12 (conjecture) | For WR networks, CF ⪯ linkage-class partition fails only when δ ≥ 1 (find the smallest WR network where CF merges two linkage classes) | WR | open; please search |
| Q2-13 (question) | Species-side Q2 for CF: every minimal siphon of N lies in one block of CF(N*). Counts and smallest failure with (H1)–(H3), WR. | | reported in logs |

Smallest examples (all exact, `logs/check_examples.log`): B3 (P ≠ Pcs), R = X₁ ⇄ X₂ (Pcs ≠ π), K3 = complete digraph on
3 species (π = 3 blocks, CF = 1, WF-FC), B2dup (π ≠ CF), XY4 = {X+Y→2X, X+Y→2Y, X+2Y→2X+Y, X+2Y→3Y}
(κ_V = 1 block at generic V, CF = π = 2 blocks).

Note for E2: use random weights in a large range for generic V (weights 1..30 hit special V, where κ_V can be
coarser than CF: Z₄ = {2X₂→2X₁, 2X₁→2X₂, X₁+X₂→2X₁, X₁+X₂→2X₂} at V = (1,1,2,2)). Families not run by Q2:
wr_2_3_6, wr_3_3_5, wr_4_2_6, wr_5_2_6_homog, wr_5_3_4, the two samples.
