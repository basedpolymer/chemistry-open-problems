# Novelty check for the gac-duality claims (second literature round)

Agent L2 — **Claude Opus 5.5 (effort medium)**, 27 October 2026. Brief: `../BRIEF-L2.md`. Builds on agent L
(`../literature/notes.md`, cited "L §n"); does not repeat L's entries unless the verdict changes. Nothing here quotes
the author's e-mail.

Verdicts: **known** (the statement is in print; cite it), **partly known** (what is in print / what is new),
**not found** (searches listed). "Read" levels: *full* = full text, *sections* = the statements and proofs quoted,
*abstract* = abstract or catalogue record only.

## 0. Main changes with respect to agent L's map

1. **C^gen (no critical siphon under the cluster rank equation) is known for ordinary nets.** Desel–Esparza,
   *Free Choice Petri Nets* (1995), Theorem 10.16 (due to J. Desel, *Regular marked Petri nets*, WG'93, LNCS 790
   (1994) 264–275, DOI 10.1007/3-540-57899-4_58) and Theorem 10.17: an **arbitrary** (not necessarily free-choice)
   ordinary net that is weakly connected, has positive S- and T-invariants and Rank(N) = |C_N| − 1 is well-formed,
   and a marking is live and bounded iff it marks every proper siphon; the proof (steps (i)–(iv)) shows that a
   marking is live as soon as it marks every S-component. Consequence (two lines, not stated in the book): put no
   token on a siphon R and one token elsewhere; if R contained no S-component the marking would be live with an
   empty siphon, impossible. So **every siphon contains the place set of an S-component (a P-semiflow support): no
   critical siphon**, for every ordinary net with the rank equation; the hypotheses are self-dual (Prop. 6.20), so
   the same holds for the reverse-dual. Exercise 10.5: such nets are covered by S- and T-components. What remains
   new in C^gen: **weighted** nets (arbitrary stoichiometric coefficients) and the reading as ODE persistence
   (through Angeli–De Leenheer–Sontag). [read: book Ch. 10, p. 214–225, local copy of the author-hosted PDF in
   `../literature/sources/`.]
2. **The sufficient half of the rank-one face criterion S1 is in print.** X. Zhang, C. Gao, D. Dochain, *On
   persistence of some WI-endotactic chemical reaction networks*, IEEE Trans. Automat. Control (2023),
   DOI 10.1109/TAC.2023.3275880 (accepted version open on the UCLouvain repository, read: Def. 2.5, Def. 3.1,
   Theorem 3.10 and its proof, §4). Definition 3.1: N is *W_I-endotactic* if it is w_I-endotactic for the 0/1
   indicator w_I of every semilocking set (= siphon) I; **Theorem 3.10: every 1-dimensional W_I-endotactic
   mass-action system is persistent** (energy-like Lyapunov function, all rate constants). Their w-endotacticity is
   two-sided (w-maximal and w-minimal sources), so it is stronger than "every face IN", but it is the same idea:
   endotacticity tested only along the indicator directions of siphons, in rank one. S1 adds the exact converse (OUT
   or ESC ⇒ non-persistent for every rate vector), the DOM case (rate-dependent), and the dichotomy for WR duals.
   Cite Zhang–Gao–Dochain next to Pantea 2012 and Craciun–Nazarov–Pantea 2013.
3. **Continuous Petri nets (Zaragoza school) already contain the "routing fixed by structure and rates" half of
   Theorem T**, under infinite-server (min) semantics, not mass action: L. Recalde, J. Júlvez, M. Silva, *Steady
   state performance evaluation for some continuous Petri nets*, 15th IFAC World Congress (2002) [read: abstract,
   §1–3, text partly garbled in the PDF]: for a class containing marked graphs and equal-conflict nets "the visit
   ratio is fixed by the net structure and the rates", and the steady-state flow is a T-semiflow so determined; the
   paper remarks that with product ("marking products") semantics limit cycles or chaos may occur (citing Silva–
   Recalde 2000). Journal follow-up: Júlvez, Recalde, Silva, Automatica 41 (2005), DOI
   10.1016/j.automatica.2004.11.007 (not read; NEEDS). No mass-action (product-semantics) uniqueness, stability or
   global-attractor theorem for EC nets was found in this line.
4. **New citing works since agent L** (§2): Feliu, arXiv:2610.11679 (persistence of bounded trajectories in
   *positively transverse* networks; cites Craciun v3 and Wiuf) and the BitMiracle/ClawsGO paper arXiv:2610.30808
   (GAC for detailed balance, Lean). Neither concerns duals, Petri nets or equal conflict.

## 1. Verdict per claim

| # | Claim (where) | Verdict | Prior art (identifier) | Read |
|---|---|---|---|---|
| 1 | Iosif's dual = weighted Petri reverse-dual (petri §1, Prop. 1.1) | **partly known** | Construction: Hack 1972 (reverse-dual), Desel–Esparza Prop. 6.20, Teruel–Silva 1996 §2 (weighted, pure nets), Jančar et al. arXiv:2602.02352 Def. 3.14; Γᵀ duals: Famili–Palsson 2003, Blanchini–Giordano 2022, Duvall arXiv:2409.12487 (L §1, §5–7, §10). **The identification itself is stated nowhere** (no CRN dual paper mentions Petri reverse-duals; no Petri paper mentions CRN duals). It is definitional: present as an observation with credits. | L: full; L2: re-searched (arXiv "reverse-dual" 11 hits, only Jančar et al. relevant; arXiv siphon ∧ dual ∧ "reaction network": only Iosif) |
| 2 | Siphons of N = T-traps of N* (petri Prop. 1.2) | **partly known** | Ordinary FC nets: Jančar et al. 2026 §4 ("minimal T-traps", T_Y• ⊆ •T_Y) as the reverse-dual counterpart of minimal siphons; Hack 1972 p. 71 had no interpretation of the reverse-dual of a trap. General weighted statement: immediate, not found in print. Agent G's Prop. 6.1 (ADLS in T-trap language) is a re-expression of Angeli–De Leenheer–Sontag 2007 Thm 2, as G says. | L: sections; L2: web search "T-trap"/"transition trap" weighted: nothing else |
| 3 | Theorem A (transversality: a minimal siphon contains no two species with the same consumers; dual for T-traps) | **partly known** | Free-choice case = Desel–Esparza **Thm 4.30(a)** ("every cluster contains at most one place of R"), same 3-line proof. For general weighted nets (places with equal postsets) not found stated, but it is the same argument. Present as a mild generalisation, credit Thm 4.30. | sections (book p. 4.30) |
| 4 | Theorem B (every precluster of N* lies in a minimal siphon iff no twin pair and every species lies in one) | **not found** | Depends on Iosif's preclusters (Iosif 2019 thesis, arXiv:2603.08767); no other source uses preclusters (arXiv "preclusters": only Iosif). New. | — |
| 5 | Lemma 1.2: Conradi–Flockerzi clusters = proportional-coordinate classes of ker A_k (A_k = Γ·diag(k)·Q), CF = meet over k | **not found** (as a statement) | Ingredients: Conradi–Flockerzi 2012 (Fact 2, Def. 3, Lemma 4, Thm 2; DOI 10.1007/s00285-011-0453-1); Pérez Millán–Dickenstein–Shiu–Conradi 2012 Condition 3.1 (kernel basis with disjoint supports; arXiv:1102.1590); Conradi–Kahle, *Detecting binomiality*, Adv. Appl. Math. 71 (2015), arXiv:1502.04893 (linear-algebra binomiality test on the same complex rate matrix); Pérez Millán–Dickenstein, MESSI systems, arXiv:1612.08763. Checked the 41 works citing CF 2012 (OpenAlex list, titles) and Conradi–Iosif–Kahle arXiv:1810.08152 (grep: no "cluster"): nobody re-reads CF clusters as kernel proportionality classes. | CF 2012 sections (on disk); 1810.08152 grep; others abstract |
| 6 | Isolation property ⇔ "internal cycles refine P̂" ⇒ toric steady states with a fixed (rate-independent) partition (preclusters §5, Cor. 5.2–5.3) | **partly known** | Isolation property: Conradi–Flockerzi 2012 Def. 3 (used there for multistationarity, Lemma 4, Thm 2). "Isolation ⇒ toric positive steady states": **Iosif's thesis (2019) Theorem 4.17** (rate-dependent form). Disjoint-support kernel ⇒ binomial/toric: PMDSC 2012 Thm 3.3, 3.8, Cor. 3.12. New: the rate-free form with the partition P̂, the equivalence (i)⇔(ii) of 5.2 and the matroid-separator reading of isolation. | thesis grep + L §16–17 |
| 7 | Theorem C (pure WF-EC: no critical siphon; minimal siphons = supports of extreme conservation laws) and C without purity (round2/ec §1) | **partly known** | Ordinary FC: Desel–Esparza Prop. 5.4 + Thm 5.6 (minimal siphons of well-formed FC nets are S-components). Ordinary arbitrary nets with the rank equation: Desel–Esparza Thm 10.16–10.17 (see §0.1). Weighted EC: Teruel–Silva 1996 give well-formedness (Thm 28) but state that siphons have no satisfactory weighted theory (p. 280, 296); **weighted version not found**. | book Ch. 5, 10: sections; TS 1996 full (L) |
| 8 | C^gen (connected, (H1), (H2), s = #clusters − 1 ⇒ N and N* have no critical siphon, persistent for all rates; K5gen: s ≥ #clusters − 1) | **partly known** | Siphon part, ordinary nets: implied by Desel 1994 / Desel–Esparza Thm 10.16–10.17 (§0.1). K5gen (lower bound s ≥ #clusters − 1 under connected + (H1) + (H2)): for ordinary FC nets it follows from Desel–Esparza Prop. 6.13 + 6.15 (well-formed ⇒ Rank = |C| − 1; (a)–(c) and not well-formed ⇒ Rank ≥ |C|); the general-net results of Ch. 10 go the other way (Thm 10.9, via Prop. 10.8: well-formed ⇒ Rank ≤ |P_N| − 1, P_N = set of distinct presets; Colom–Campos–Silva, corrected by Teruel–Silva 1994). So **K5gen not found** outside FC (and P's K5 for EC not found either). Weighted siphon statement and persistence reading: new. | sections |
| 9 | Theorem T (WF-EC: ker A_k = ℝu(k), u(k) > 0 for all k; toric steady-state set; converse) and T* (self-duality for balanced joins) | **partly known** (credits already in P §7) | PMDSC 2012 Thm 3.8/Cor. 3.12 (toric part once the kernel is a positive line); Campos–Chiola–Silva, IEEE TAC 36 (1991) (flux part for ordinary FC with rational routing); **add** Recalde–Júlvez–Silva IFAC 2002 and Júlvez–Recalde–Silva Automatica 2005 (continuous EC nets, infinite-server semantics, visit ratios fixed by structure and rates). Teruel–Silva 1996 Thm 7 for CF nets (unique minimal T-semiflow). New: mass-action reading for weighted EC with real rates, the converse (2), T* and T*′. | IFAC 2002: abstract + §1–3 |
| 10 | Theorem U (WF-EC: exactly one positive steady state per class for every k; Jacobian Hurwitz on S at every positive point; σ(S) ∩ σ(D^⊥) = {0}) | **not found** (as a theorem); tools known | Uniqueness from the sign condition: Müller–Regensburger, SIAM J. Appl. Math. 72 (2012), DOI 10.1137/110847056, Prop. 3.9; Müller–Feliu–Regensburger–Conradi–Shiu–Dickenstein, Found. Comput. Math. 16 (2016), arXiv:1311.5493 (maximal minors). Spectral trick SG ↔ GS: Jamshidi–Palsson 2009, Blanchini–Giordano 2022 Thm 1. Metzler matrix with positive null vector: Perron–Frobenius. Related classes that do **not** contain WF-EC: Donnell–Banaji arXiv:1211.2153 (factorisation Γ = ΛΘ), node-balanced steady states (Feliu–Cappelletti–Wiuf arXiv:1711.07074), translated networks (Johnston arXiv:1305.5845, DOI 10.1007/s11538-014-9947-5: toric steady states via a translated WR δ = 0 GMA system — closest in spirit; a WF-EC network is generally not translatable to WR with sources as kinetic complexes, not checked in general), single-target networks (Boros–Craciun–Yu arXiv:2006.01192). Searched: "reactant complexes" + disjoint + uniqueness/stability (web, two queries), arXiv "equal conflict" (0 hits), "free-choice" ∧ reaction (only Mairesse–Nguyen). | abstracts; Donnell–Banaji abstract |
| 11 | Theorem GS (WF-EC: global attractor in every class; Lyapunov functions max_c x^{y_c}/u_c, min_c …) | **not found**; technique known | The functions are the max/min of ψ_c/u_c with ψ̇ = G(x)ψ, G Metzler, G u = 0: this is the classical consensus/compartmental argument (Maeda–Kodama–Ohta, IEEE Trans. Circuits Syst. 25 (1978), nonlinear compartmental systems, max-ratio Lyapunov functions; Moreau, IEEE TAC 50 (2005), time-dependent consensus; DOIs in NEEDS, not opened) and piecewise-linear-in-rates Lyapunov functions (Al-Radhawi–Angeli, IEEE TAC 61 (2016), DOI 10.1109/TAC.2015.2427691, robust version). Join-free first-order case (sources = single species, coefficient 1) is a linear compartmental system: classical. WR EC networks: GAC known (deficiency 0, linkage classes on disjoint species: Anderson 2011 single linkage class). New: the EC (weighted, non-WR, non-complex-balanced) global attractor statement. | abstracts |
| 12 | T*′, U*, BH (dual of WF-EC: absolute consistency, uniqueness, Hurwitz; block-homogeneous networks) | **not found** | Weighted Petri literature on "homogeneous" nets (Hujsa–Delosme–Munier-Kordon, Fundam. Inform. 2016; Hujsa–Devillers 2018): *homogeneous* = all output arcs of a place have equal weight — EC is homogeneous, but BH (sources a·1_B) and PC are different conditions; only discrete liveness/reversibility results. | abstracts |
| 13 | PC (proportional conflict: self-dual class; existence, uniqueness, Hurwitz, persistence on both sides) | **not found** | No "proportional conflict" class in the Petri literature found (web search; Teruel–Silva use EC, TEFC, balanced joins; Hujsa uses homogeneous). New. | — |
| 14 | ACR, injectivity, sign conditions, deficiency (brief's list) | — | Checked for overlap: injectivity criteria (Müller et al. 2016; Banaji–Craciun; Feliu–Wiuf) give "at most one steady state per class for all k" once the sign condition is verified — used, not replaced; ACR (Shinar–Feinberg) needs deficiency one with two non-terminal complexes differing in one species — WF-EC steady states x^{y_c} = λu_c(k) show ACR only when some source is k·X_i with i in no other source and λ fixed, not a theorem in the notes; deficiency theory: WR EC ⇒ δ = 0 is Mairesse–Nguyen 2009 (L §11). | — |
| 15 | S1 rank-one face criterion (persistent for all / some / no rates) and S1-N, C-pure, C-2 (round2/persistence §1–2) | **partly known** | Sufficient half: **Zhang–Gao–Dochain 2023, Thm 3.10** (§0.2). Lower-endotactic sweep tests: Pantea 2012 (arXiv:1103.0603) and Craciun–Nazarov–Pantea 2013 (arXiv:1010.3050); semilocking sets: Johnston–Siegel, SIAM J. Appl. Math. 71 (2011), DOI 10.1137/100807557 (title/abstract only). Stochastic rank one: Wiuf–Xu arXiv:2012.07954; rank-one multistationarity: Pantea–Voitiuk arXiv:2208.06310 (different question). New: the exact trichotomy with OUT/ESC ⇒ non-persistent for every rate vector, the DOM case, translations S1-N, C-pure, C-2 to N. | Zhang et al.: sections; others: abstract (via Q1 and arXiv) |
| 16 | S1-WR dichotomy at s = 1, Theorem W (no DOM siphon facet of N* for WR N, any rank), Prop. WR-Σ | **not found** | Nothing on persistence of duals of WR networks. WR-Σ (Yp = Ye·Σ for WR) is an elementary observation. | searches below |
| 17 | Drain species (E.4), drain sets (Theorem D), smallest examples X, C | **not found** (closest related listed) | Deshpande–Gopalkrishnan arXiv:1309.3957 (drainable / self-replicable siphons; reachability notion, opposite direction); Gopalkrishnan, *Catalysis in reaction networks*, Bull. Math. Biol. 73 (2011) (critical siphons ↔ catalysis; DOI in NEEDS); Agarwal–Craciun–Deshpande–Jin arXiv:2601.08801 (weak/strong extinction via linear Lyapunov functions; Thm 3.12: inconsistent ⇔ a linear Lyapunov function for all rates; Thm 1.2: δ = 0 non-WR with bounded class ⇒ weak extinction for all rates) — **not applicable to N*** because N* is consistent when N is conservative, and a consistent δ = 0 network is WR; Johnston–Avram arXiv:2506.01606 (boundary reproduction number: rate-dependent local stability of boundary steady states via siphons). The drain-set certificate (local dominance near a face, all rates) was not found. | 2601.08801: sections; others abstract |
| 18 | Rank-two ALL/NONE certificates (A, P2, E, H, CM, L, S, B), Conjecture F, vertex ties U (round3/dichotomy) | **partly known** (tools) / **not found** (statements) | Tools are standard planar techniques (Poincaré–Bendixson, Butler–McGehee/acyclicity: Butler–Waltman 1986, Hale–Waltman 1989, on disk; centre manifolds; Pantea 2012 for persistence of bounded 2-dimensional WR trajectories; Anderson–Shiu 2010 facets repel for WR). Applied to non-WR duals with the stated certificates: not found. | L §13; titles |
| 19 | Non-self-duality of persistence (3A ⇄ 2A+B, gac §3) and of the ADLS criterion (gac §5); counts E-e2/e3 | **partly known** | Blanchini–Giordano 2022 p. 6: dual of a reversible network is not WR, mass-action form not preserved (L §1). Non-self-duality of persistence, of global attraction and of the ADLS siphon criterion: not found (Blanchini–Giordano list "properties other than stability" as future work). | L: full |
| 20 | Prop. 4.1 of gac (q = 1, N and N* WR, n prime ⇒ N strongly endotactic ⇒ permanence; GAC for complex-balanced rates) | **not found**; ingredients classical | Irreducibility of the sum-zero space under the cyclic shift for prime n (Φ_n irreducible, Eisenstein): textbook. Conclusion is Gopalkrishnan–Miller–Shiu 2014 (DOI 10.1137/130928170). arXiv abs: "strongly endotactic" ∧ cyclic: 0 hits. New as a recognition criterion; it gives no new GAC case (GAC claimed in full, L §13). | — |
| 21 | P2-type statements "duality keeps s, swaps n and r; P(N*) discrete for WR" (exploration facts) | known/elementary | Γ(N*) = Γᵀ (Iosif; Blanchini–Giordano Def. 1). | L |

## 2. Citing works (refreshed 27 October 2026, evening)

Sources: Semantic Scholar Graph API (worked intermittently; HTTP 429 retried), OpenAlex (free budget), arXiv API.

- **Iosif, arXiv:2603.08767**: Semantic Scholar citationCount 0 (paper found, `citations: []`); OpenAlex: not
  indexed by arXiv id. **No citing work.**
- **Blanchini–Giordano 2022** (DOI 10.1109/LCSYS.2021.3081369): Semantic Scholar 7 citing works: Blanchini–Devia–
  Giordano arXiv:2109.01709, Al-Radhawi CDC 2022 and arXiv:2302.05090, Duvall arXiv:2409.12487, Nilsson–Giordano
  ECC 2025, Iosif arXiv:2603.08767, and **new: A. Stevens, N. Vassena, *Mathematical modeling and analysis of the
  Notch–Delta pathway*, arXiv:2604.05888 (DCDS 2026, DOI 10.3934/dcds.2026.105533)** [abstract: structural
  (parameter-rich) analysis of symmetry-induced bifurcations; no dual-network content in the abstract].
- **Craciun, arXiv:1501.02860** (all versions): 40+ citing works on Semantic Scholar; those after v3 (23 Sep 2026):
  **E. Feliu, *Persistence of bounded trajectories in positively transverse reaction networks*, arXiv:2610.11679
  (8 Oct 2026)** [abstract: a polytope condition on a subset of complexes ⇒ no bounded positive trajectory has a
  boundary ω-limit point] and **arXiv:2610.30808 (22 Oct 2026), *Global convergence of detailed balanced mass action
  systems*** (BitMiracle / "ClawsGO Science", Claude Opus 5.5 + GPT-6-Astra, Lean; builds on OpenAI's barrier) —
  the paper behind the BitMiracle repository already listed by L. Earlier 2026 citers (disguised toric locus,
  extinction, relaxation dynamics) are unrelated to duality.
- **Wiuf, arXiv:2609.24553**: one citing work, Feliu arXiv:2610.11679 (above).
- **OpenAI, *Boundedness and persistence of WR mass-action systems* (25 Sep 2026)**: not indexed by Semantic Scholar
  or OpenAlex (no DOI/arXiv id); cited by arXiv:2610.30808 (abstract) and by mmislan's Palomar repositories (L §15).
- **Conradi–Flockerzi 2012** (for Lemma 1.2): 41 citing works (OpenAlex titles read); none on clusters as kernel
  classes, Petri nets or duals.

## 3. AI-campaign scans (PROTOCOL.md §1), redone on the new key phrases

- `gh api users/mmislan/repos`: 27+ Palomar repositories; new since L: palomar-lienard-limit-cycles (22 Oct),
  template-free-heredity, joint-thermodynamic-compatibility (21 Oct), compartment-division-threshold (19 Oct),
  shared-resource-realizations (11 Oct; README read: two-species modules with a shared resource, realizations,
  permanence via endotactic permanence — not equal conflict, not duals). Nothing on duals, siphon/T-trap, equal
  conflict, drain sets or persistence of dual networks.
- `gh api users/ipitchford/repos --paginate`: 201 repositories; grep of names and descriptions for dual, Iosif,
  siphon, precluster, attractor, Petri, transpose, toric, persist, conflict, drain, reaction, chemical, mass-action,
  endotactic: CRN hits only quadratic-mass-action-sharpness (4 Oct) and reverse-reaction-attracting-oscillations
  (18 Oct); unrelated.
- `gh search repos` on "reverse-dual", "equal conflict", "equal-conflict", "drain species", "persistence dual
  network", "dual reaction network", "transpose stoichiometric", "Iosif duality", "preclusters", "siphon
  persistence", "free-choice reaction network": no relevant repository. GitHub search API with `in:readme` on
  "global attractor conjecture", "siphon reaction network", "petri net reverse dual", "equal conflict petri": only
  this project's public repository is relevant; one knowledge-base repository (harjassand/astra-research-knowledge,
  50 000+ files) surfaced; its tree has no Petri-siphon/CRN-duality card (one card on join dynamics vs fractional
  Petri firing, unrelated).
- github.com/openai/math: 751 preprint folders; CRN-related still only the two WR papers of 25 Sep and 5 Oct 2026.
- Zenodo (API through the web fetcher): `"reaction network" AND (dual OR siphon OR "equal conflict" OR
  persistence)`: 8 records, none on duals, reverse-duals, equal conflict or siphons (one on persistence of a
  kinetic-proofreading network, 30 Sep 2026).
- arXiv API (date-sorted): "reverse-dual" (11, only Jančar et al. relevant), "equal conflict" (0), siphon ∧ dual ∧
  "reaction network" (only Iosif), drain ∧ "reaction network" ∧ persistence (0), "reaction network" ∧ transpose (0),
  "dual network" ∧ "mass-action" (0), "free-choice" ∧ reaction (only Mairesse–Nguyen), "global attractor
  conjecture" (latest Wiuf 2609.24553), "reaction network" ∧ Petri ∧ persistence (4: ADLS 2006, Nabli–Fages–
  Martinez 2013, Szederkényi–Ács–Lipták 2022 (kinetic compartmental models, siphons, persistence; different
  class), Johnston–Avram 2025).

## 4. What could not be searched

- Google Scholar (not reachable); Semantic Scholar only intermittently (citation lists for Iosif, Wiuf and Craciun
  obtained once each; a later refresh may add very recent works).
- Paywalled sources not read (NEEDS.md): Desel 1994 (only through Desel–Esparza's account, which suffices for the
  statement), Júlvez–Recalde–Silva 2005, Teruel–Silva 1994, the classical Lyapunov references for GS
  (Maeda–Kodama–Ohta 1978, Moreau 2005), Gopalkrishnan 2011, Hujsa et al. 2016.
- The continuous-Petri-net literature (Silva, Recalde, Júlvez, Mahulea) is large; only the IFAC 2002 paper was
  opened. A targeted check of "product semantics" results for EC nets would settle whether Theorems U/GS have a
  continuous-Petri-net predecessor (none found in two searches).
- WF-EC versus network translation (Johnston 2014): whether every WF-EC network is translatable to a WR δ = 0
  generalized mass-action system (which would give Theorem T and part of U from Johnston's Theorem 5 and
  Müller–Regensburger) was not decided here; worth one check by the verifier of the EC block.

## 5. Recommended credits for the paper (short list)

Desel–Esparza Thm 4.30, Prop. 5.4, Thm 6.14/6.21, **Thm 10.16–10.17 (Desel 1994)**; Teruel–Silva 1996 Thm 28,
Cor. 32; Hack 1972; Jančar et al. 2026; Mairesse–Nguyen 2009; PMDSC 2012; Campos–Chiola–Silva 1991;
**Recalde–Júlvez–Silva 2002**; Conradi–Flockerzi 2012 (isolation) and **Iosif thesis Thm 4.17**; Müller–Regensburger
2012, Müller et al. 2016; Jamshidi–Palsson 2009, Blanchini–Giordano 2022; Al-Radhawi–Angeli 2016 and the
compartmental/consensus max-ratio Lyapunov functions; Angeli–De Leenheer–Sontag 2007; Pantea 2012,
Craciun–Nazarov–Pantea 2013, **Zhang–Gao–Dochain 2023**; Deshpande–Gopalkrishnan 2014; Agarwal–Craciun–Deshpande–Jin
2026 (to say why it does not apply to duals); Gopalkrishnan–Miller–Shiu 2014.
