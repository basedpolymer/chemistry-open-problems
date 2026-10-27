# Changes made in writing the unified manuscript (writer W)

Writer: Claude Opus 5.5 (effort xhigh), 27–28 October 2026. This file lists every place where the manuscript
(`paper.tex`) merges, generalises, re-proves or changes a statement **beyond the corrections of the verifiers V1–V3**
(those are applied as written in their reports and listed in §9.2 of the paper). Numbers are the paper's
(theorem-like items, tables and figures share one counter); labels in brackets are the LaTeX labels.

## A. New proofs (elementary re-proofs of claims that the notes or verifiers settled by computer algebra)

1. **Proposition 44 [prop:NC3], NC3 has exactly three positive steady states in its class.** V2 proved the count
   with msolve (exact real-root isolation) and gave the stability of the two non-saddle points numerically
   (eigenvalues). The paper proves it by hand: f_B = −f_D and f_A + f_C + f_D = 0 reduce steady states to
   D = 36A²B, C = 18AB + 36A²B (at k = (1/6, 1/3, 6, 6)); the class gives B = (10/3)/(1 + 36A²) and
   72A³ − 204A² + 122A − 19 = (2A − 1)(36A² − 84A + 19) = 0, whose three roots A = 1/2, (7 ± √30)/6 are positive and
   give positive B, C, D. Stability on S of the two irrational points is decided by exact trace and determinant in
   ℚ(√30). Certified: claims "S4 NC3: ..." (four claims). The statement itself is V2's.
2. **NC (second CR counterexample, after Proposition 44).** The paper states only what the degree argument gives
   rigorously from the certified saddle: "at least one further positive steady state" (index sum +1, saddle index −1),
   and credits V2's exact count ("exactly two more") to V2's script. This is weaker than V2's statement, not stronger.
3. **NM (Table 45 and the paragraph after it): exactly two positive steady states in x₁+x₂+x₃ = 5/2.** Proved by the
   reduction x₃ = 1/2, x₁x₂² = 1/4, x₁ + x₂ = 2 and the sign pattern of t³ − 2t² + 1/4 at −1, 0, 1, 2 (three real roots,
   one negative). Q3/V2 used exact sign certificates and msolve; same statement.
4. **NB:** the two steady states x₁ = (5 ± √17)/4 are written in closed form and checked in ℚ(√17); same statement as
   Q3/V2.
5. **Example 50 [ex:M], class u + v = 2 of D = (3A ⇄ 2A+B)\* with unit rates has exactly two positive equilibria.**
   V3 checked it with sympy real-root isolation; the paper's certificate uses h(v) − 1 = −(v − 1)(v³ − 5v² + 7v − 1),
   the negative discriminant (−44) of the cubic and its sign change on (1/8, 1/2). Same statement.

## B. Merged statements (several results of the notes stated as one item; content unchanged unless said)

6. **Theorem 9 [thm:criteria]** merges the archived Theorem 7 and Corollaries 8–10 (criteria on both sides, the
   no-twin corollary for networks whose products are sources, the weakly reversible reaction-side corollary). The proof
   is condensed; the full proofs are in the archived report.
7. **Theorem 11 [thm:archive] and Table 10 [tab:formulations]** condense the archived Theorem 1 and Table 1: the table
   lists the counterexamples for the ray reading of MIPS only; the archive treats three readings.
8. **Proposition 15 [prop:hierarchy]** merges Q2 Props 2.1 and 2.2 with the table of witnesses.
9. **Theorem 16 [thm:invariance]** merges Q2 Theorem 3.1 and Propositions 3.2, 3.3 (invariance, non-invariance of
   P, P_cs, π, and "full coupling is forced").
10. **Theorem 17 [thm:Phat-structure]** merges Q2 Theorems 4.1, 4.2 and E2's Proposition E2.6 into parts (1)–(3).
11. **Proposition 18 [prop:toric]** merges Q2 Proposition 5.1 and Corollaries 5.2, 5.3. *Weakened in the statement:*
    of Corollary 5.2 the paper states the inequality rank A_V ≥ m − #P̂ and the implication "equality for every V ⇒
    toric steady states for every rate vector with the rate-independent partition P̂", but not the converse
    implication (ii) ⇒ (i) of Q2's equivalence.
12. **Theorem 19 [thm:duality-redone]** merges Q2 Theorem 6.1 and Corollaries 6.2, 6.3, with V1's wording ("for every
    consistent / conservative network"). Q2's Proposition 6.4 (Question 24 for P̂) appears as a sentence after
    Proposition 23 [prop:G].
13. **Theorem 24 [thm:WFFC]** condenses P's Theorem D (a)–(g) and the first sentence of P's Theorem E(ii);
    omitted: D(h) (liveness from the all-ones marking) and the JF/CF characterisations of E(ii). Proof given as a
    sketch with references to Desel–Esparza and to P's notes.
14. **Corollary 32 [cor:C]** (well-formed EC: no critical siphon on either side, persistence) is presented as a
    corollary of Q3's cluster-rank Theorem 31 [thm:Cgen] instead of P's Theorems C (pure, via Teruel–Silva) and C′
    (limit of rate constants) and Q3's direct proof of C; Q3's Remark (i) already notes that C^gen contains C and C′.
    The sentence "the minimal siphons of N are the internal cycles of N\* and conversely" is P's bijection stated
    without purity; it follows from Theorem 31 and Lemma 2 (transposition).
15. **Theorem 33 [thm:T], proof of (1)** refers to Q3's Lemma R (the cluster version of P's closed–open argument)
    with q = m, instead of repeating P's argument; same statement.
16. **Proposition 41 [prop:Tstar-cex]** merges Q3's Propositions T\*-cex, T\*-dim and tor\*.
17. **Theorem 47 [thm:WRdual]** merges E2's Corollary E2.5′ (existence, Birch) and A2's Theorems E, E+ (manifold),
    as V3 asks. The hypothesis is written as A2's "conservative and im Γ ⊆ im Ye"; E2's "every product complex lies
    in the span of the sources" is equivalent (Yp = Γ + Ye).
18. **Proposition 48 [prop:nonWR]** merges G's Propositions 2.1–2.2 (dual WR criterion, deficiency of the dual) with
    E's counts; the criterion is written with a permutation π of the species, (Yp)_{i·} = (Ye)_{π(i)·}, which is G's
    "b_i = a_{π(i)}".
19. **Example 53 [ex:drain]** gathers examples of E (4-cycle), E2 (5-cycle, refutation of Q1-6) and Q1 (example C).
    For the 5-cycle the drain species is stated as "the species dual to X₁+2X₃ → 3X₃" (computed by the certificate in
    path order), not with E2's labels.
20. **Corollary 61 [cor:Cpure]** merges Q1's C-pure (with V3's (H1)) and C-2 (with (H1)–(H3), as in Q1).
21. **Proposition 66 [prop:none-cert]** lists D2's NONE certificates E (restricted to μ = 1, V3), H, CM, L, S, B in
    one proposition without proofs (proof sketches in D2's notes; V3 checked CM, L, S at sketch level).
22. **Proposition 72 [prop:dimension]** merges G's Proposition 2.1(1)(3), the monomolecular-cycle example of G §4 and
    G's remark on rank-preserving operations. **Proposition 74 [prop:dictionary-boundary]** merges G's Propositions
    6.1–6.3.
23. **Theorem 58 [thm:S1]** is stated with V3's wording of (4); Q1's proof of (4) ("example B") is replaced by the
    semialgebraicity argument V3 gives (Tarski–Seidenberg) plus example B.

## C. Statements dropped or not repeated

24. **P §7, remark (i): "W\*, F and B_k are absolutely consistent and toric but not complex balanced."** The part
    about B_k is wrong: B_k is weakly reversible with deficiency zero (k + 2 complexes, one linkage class, s = k + 1),
    hence complex balanced for every rate vector. The paper says so in a parenthesis after Theorem 37 [thm:GS], and the
    certificate checks it (claim "S4 B2, B3, B4 are weakly reversible with deficiency 0"). No verifier flagged this
    remark; it does not affect any theorem.
25. P's Corollary P (persistence for pure or balanced WF-EC networks) is not stated: it is superseded by Corollary 32
    (V2 §7).
26. Q3's Conjecture CR is not stated as a conjecture; its siphon part is Theorem 31, its dynamical part is refuted
    (Proposition 44), and its absolute-consistency part is open problem 4 (§8), as V2 prescribes.
27. Agent E's claim that the smallest WR network with a non-endotactic dual is 2X₂ ⇄ X₁ is not repeated (I could not
    re-derive non-endotacticity of its dual by hand in the time available); only E's counts are cited, as "checked by
    computation".

## D. Notation unified across the notes (no change of content)

- D = span{y_c − y_c′} (Q3's letter; E2 wrote L), L = S^⊥ for conservation laws (V2 item 5).
- WF-EC defined algebraically (EC, connected, (H1), (H2), s = m − 1); Petri well-formedness mentioned only for pure
  nets (V2 item 6).
- q = number of Petri clusters (Q3, V2); the codimension is written n − s (G used q for it).
- Species of N\* are K_1..K_r with rate constants κ_1..κ_n (the notes use k_j); Theorems A, AF, O keep A2's a, b, c and
  α, β, γ, δ.
- Twins include equal zero rows of Ye (V1 correction 1).
