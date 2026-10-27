# References and novelty check

Sources of the manuscript `paper.tex` (rewritten 28 October 2026) and of the archived report
`archive/paper-precise-forms.pdf` (25 October 2026, revised 27 October 2026). For each source: identifier, and what was
read (full text / sections / abstract only), and by whom when it was read in the working notes of `gac-duality/`
(agents L and L2 for the literature, P, Q1–Q3, E, E2, G, D2, A2 for the research, V1–V3 for the verification). Source
PDFs are kept, git-ignored, in `certificates/sources/` and in the `sources/` folders of `gac-duality/`.

## Primary source (full text read)

- A. Iosif, *Duality in mass-action networks*, arXiv:2603.08767. Read in full: v1 (9 March 2026) and v2 (15 July 2026).
  The mathematical content is identical; v2 adds Remark 7 and minor rewording. HAL hal-03058670 v3 (17 July 2026) is
  identical to arXiv v2. Quoted: Conjecture 22, Remark 23, Question 24 (§3, p. 7 of v2), the paragraph after Question 24
  on the GAC "in small codimensions" (p. 7–8), the abstract and p. 2, and Definition 19. Definitions copied: 13–14
  (siphons), 17 and Remark 18 (MIPS), 19 (preclusters), 20 (internal cycles), the dual network (§3).
- A. Iosif, *Duality in mass-action networks*, HAL hal-03058670v1 (11 December 2020). Read in full (4 pages). It states
  the duality as a theorem without proof and without the positive-steady-state hypothesis.
- A. Iosif, *Algebraic methods for the study of multistationarity in mass-action networks*, PhD thesis,
  Otto-von-Guericke-Universität Magdeburg (2019), DOI 10.25673/13949 (PDF from the Internet Archive copy of the
  repository file). Read: §1.6 (internal cycles), Chapter 4 (doubling sets, Definition 4.3 and Algorithm 4.4 for
  preclusters, Definition 4.6 for clusters, Theorem 4.17, Examples 4.1–4.14).
- A. Iosif, *Multistationarity in systems with toric states*, seminar slides, 19 November 2024,
  https://gasiull.webs.ull.es/beamer-2024-11-19.pdf. Read: slides 28–29.
- A. Iosif, L. Ananiadi, poster P1.H, Jóvenes Investigadores RSME 2025 (Bilbao). Abstract only.

## Archived report of this project

- *Preclusters, maximal invariant polyhedral supports and siphons: exact criteria and counterexamples for the precise
  forms of Iosif's duality conjecture*, `archive/paper-precise-forms.pdf` (source `archive/paper-precise-forms.tex`,
  bibliography `archive/references-precise-forms.bib`), certified by `scripts/certify_examples.py` (142 claims) and the
  enumeration scripts; verified by GPT-6-Astra (effort xhigh), `verification.md`. Section 3.1 of the new paper condenses
  it.

## Reaction-network background

- C. Conradi, D. Flockerzi, J. Math. Biol. 65 (2012) 107–156, DOI 10.1007/s00285-011-0453-1. Publisher PDF on disk.
  Read: §2, §3.1–3.2 (Assumptions 1–2, Fact 2, Algorithm 1, Remark 5, Example 1, Lemma 1, Remark 6), Definition 3
  (isolation). V1 re-read Fact 2 and Algorithm 1 for the reading of the clusters used in Lemma 14 of the paper (the
  domain of the ratio identity is not stated precisely; the paper states its reading).
- M. Pérez Millán, A. Dickenstein, A. Shiu, C. Conradi, *Chemical reaction systems with toric steady states*, Bull. Math.
  Biol. 74 (2012) 1027–1065, DOI 10.1007/s11538-011-9685-x, arXiv:1102.1590. Read: §3 (Condition 3.1, Theorems 3.3, 3.8,
  Corollary 3.12) and §5 (agents P, Q2; V2 checked Theorem 3.8 and Corollary 3.12).
- S. Schuster, C. Hilgetag, J. Biol. Syst. 2 (1994) 165–182, DOI 10.1142/S0218339094000131 (elementary flux modes).
  Not re-read; cited for terminology.
- A. Shiu, B. Sturmfels, *Siphons in chemical reaction networks*, Bull. Math. Biol. 72 (2010) 1448–1463, DOI
  10.1007/s11538-010-9502-y. §1–3 read by agent L.
- D. Angeli, P. De Leenheer, E. D. Sontag, *A Petri net approach to the study of persistence in chemical reaction
  networks*, Math. Biosci. 210 (2007) 598–618, DOI 10.1016/j.mbs.2007.07.003. Journal version on disk; Theorem 2
  (p. 605–607) and §5 read (agents P, G, Q3; V2, V3).
- A. P. Burgard et al., Genome Res. 14 (2004) 301–312, DOI 10.1101/gr.1926504; T. Pfeiffer et al., Bioinformatics 15
  (1999) 251–257, DOI 10.1093/bioinformatics/15.3.251 (full coupling, enzyme subsets). Not re-read; terminology.
- M. Feinberg, *Foundations of Chemical Reaction Network Theory*, Springer 2019, DOI 10.1007/978-3-030-03858-8;
  F. Horn, R. Jackson, Arch. Rational Mech. Anal. 47 (1972) 81–116, DOI 10.1007/BF00251225 (deficiency, complex balance,
  Birch's theorem). Standard facts, not re-read.
- F. Blanchini, G. Giordano, *Dual chemical reaction networks and implications for Lyapunov-based structural
  stability*, IEEE Control Syst. Lett. 6 (2022) 488–493, DOI 10.1109/LCSYS.2021.3081369. Publisher PDF on disk, read in
  full (agents L, G; V3): Assumptions 1–4 (p. 2), Theorems 1–3, Remark 2, p. 493 (the sixth page; the dual of a reversible network is not
  weakly reversible; mass action is not preserved).
- I. Famili, B. Ø. Palsson, *The convex basis of the left null space of the stoichiometric matrix leads to the
  definition of metabolically meaningful pools*, Biophys. J. 85 (2003) 16–26, DOI 10.1016/S0006-3495(03)74450-6
  (open access, PMC1303061). Compound-map passages read by agent L (`gac-duality/literature/notes.md` §10(a)): the
  compound map of −Sᵀ, with the extreme conservation laws drawn as extreme pathways.
- N. Jamshidi, B. Ø. Palsson, *Flux-concentration duality in dynamic nonequilibrium biological networks*, Biophys. J.
  97 (2009) L11–L13, DOI 10.1016/j.bpj.2009.06.049 (open access, PMC2749767). Abstract and main text read by agent L
  (§10(b)): the Jacobians S·G and G·S.
- A. Duvall, *Computationally checking if a reaction network is monotone or non-expansive*, arXiv:2409.12487 v2 (27
  January 2025). §11 ("Duality") read by agent L (§10(e)): the dual of an irreversible network has stoichiometric
  matrix Γᵀ; monotonicity and non-expansivity transfer.
- M. D. Johnston, *Translated chemical reaction networks*, Bull. Math. Biol. 76 (2014) 1081–1116, DOI
  10.1007/s11538-014-9947-5, arXiv:1305.5845. §4.1 (Definitions 4.1, 4.2, Theorem 4.1) read by V2 from a scratch copy;
  bibliographic data matched through Crossref on 28 October 2026.
- C. Conradi, E. Feliu, M. Mincheva, C. Wiuf, *Identifying parameter regions for multistationarity*, PLoS Comput. Biol.
  13 (2017) e1005751, DOI 10.1371/journal.pcbi.1005751 (degree argument, Theorem 1). Cited after agent Q3; not re-read;
  DOI matched through Crossref on 28 October 2026.
- S. Müller, G. Regensburger, SIAM J. Appl. Math. 72 (2012) 1926–1947, DOI 10.1137/110847056, arXiv:1209.6488
  (sign-vector criterion: Proposition 3.1 of arXiv v1, read by reviewer R2 on 28 October 2026 and checked again by
  reviewer R1 from the arXiv PDF; the numbering of the journal version was not checked, so the paper cites the arXiv
  version). Cited after agents Q3, L2; title and pages matched through Crossref.
- M. A. Al-Radhawi, D. Angeli, IEEE Trans. Automat. Control 61 (2016) 76–89, DOI 10.1109/TAC.2015.2427691
  (piecewise-linear-in-rates Lyapunov functions). Abstract (L2); title matched through Crossref on 28 October 2026.
- Y. Maeda, S. Kodama, Y. Ohta, IEEE Trans. Circuits Syst. 25 (1978) 372–378; L. Moreau, IEEE Trans. Automat. Control 50
  (2005) 169–182 (max/min consensus and compartmental Lyapunov functions). Not opened; cited after L2 for the technique
  of Theorem 37. DOIs 10.1109/TCS.1978.1084490 (full title: *Asymptotic behavior of nonlinear compartmental systems:
  nonoscillation and stability*) and 10.1109/TAC.2004.841888, matched through Crossref by reviewer R2 on 28 October 2026.
- X. Zhang, C. Gao, D. Dochain, *On persistence of some W_I-endotactic chemical reaction networks*, IEEE Trans. Automat.
  Control 68(12) (2023) 7420–7435, DOI 10.1109/TAC.2023.3275880. Accepted version read (Definitions 2.5, 2.11, 3.1, Theorem 3.10 and its
  proof; L2 and V3).

## Petri nets

- M. H. T. Hack, *Analysis of production schemata by Petri nets*, MS thesis, MIT, 1972; Project MAC TR-94, handle
  1721.1/149406. §5.2, §5.5 (Well-Formedness Theorem 9, the reverse dual; no interpretation of the reverse dual of a
  trap) read by agent L, checked by V1.
- J. Desel, J. Esparza, *Free Choice Petri Nets*, Cambridge Tracts in Theoretical Computer Science 40, CUP 1995, ISBN
  0-521-46519-2, DOI 10.1017/CBO9780511526558 (Crossref) (author-hosted PDF). Read: Def. 4.4, Defs 5.1, 5.11, Thm 4.30, Ch. 5 (Prop. 5.2, 5.4, 5.7, Thm 5.17), Ch. 6 (Thms 6.14,
  6.17, Def. 6.19, Prop. 6.20, Thm 6.21, Ex. 6.6), Ch. 10 (Thms 10.16–10.17 and the proof of 10.17) (agents P, L, L2;
  V1, V2).
- J. Desel, *Regular marked Petri nets*, WG'93, LNCS 790 (1994) 264–275, DOI 10.1007/3-540-57899-4_58 (Crossref). Not
  read; cited through Desel–Esparza's account of Theorem 10.16.
- E. Teruel, M. Silva, *Structure theory of equal conflict systems*, Theor. Comput. Sci. 153 (1996) 271–300, DOI
  10.1016/0304-3975(95)00124-7. Publisher PDF on disk, read in full
  (agent L; Def. 6, Thms 7, 12, 15, 28, Cor. 13, 32 checked by V1, V2).
- P. Jančar, E. Best, R. Devillers, M. Ošťádal, *Well-formed free-choice Petri nets revisited*, arXiv:2602.02352 v2 (30
  July 2026). §1–4 read (agent L; Theorem 3.23, §4 "minimal T-traps", checked by V1).
- J. Mairesse, H.-T. Nguyen, *Deficiency zero Petri nets and product form*, Petri Nets 2009, LNCS 5606, 103–122, DOI
  10.1007/978-3-642-02424-5_8, arXiv:0905.3158. §3.3 (Lemma 3.2, Corollary 3.2) read (agent L; V1).
- J. Campos, G. Chiola, M. Silva, *Properties and performance bounds for closed free choice synchronized monoclass
  queueing networks*, IEEE Trans. Automat. Control 36 (1991) 1368–1382, DOI 10.1109/9.106153 (Crossref). Author's copy
  (https://webdiis.unizar.es/~jcampos/web/pdf/ccs-tac-91.pdf) read: Lemma 5.1, Theorems 5.3–5.4 (agent P; V2).
- L. Recalde, J. Júlvez, M. Silva, *Steady state performance evaluation for some continuous Petri nets*, 15th IFAC World
  Congress (2002), IFAC Proc. Vol. 35(1) 455–460, DOI 10.3182/20020721-6-ES-1901.00078 (Crossref). Abstract and §1–3 read (L2).
- J. Júlvez, L. Recalde, M. Silva, *Steady-state performance evaluation of continuous mono-T-semiflow Petri nets*,
  Automatica 41 (2005) 605–616, DOI 10.1016/j.automatica.2004.11.007. Publisher PDF on disk; eq. (18), Definition 7 and
  the semantics read by the coordinator and V2.

## Persistence and dynamics

- G. Butler, P. Waltman, *Persistence in dynamical systems*, J. Differential Equations 63 (1986) 255–263, DOI
  10.1016/0022-0396(86)90049-5; J. K. Hale, P. Waltman, *Persistence in infinite-dimensional systems*, SIAM J. Math.
  Anal. 20 (1989) 388–395, DOI 10.1137/0520025. Publisher PDFs on disk; Butler–Waltman Theorem 4.1, Lemma 2.1,
  Theorem 3.1 and Hale–Waltman Lemma 4.3, Theorem 4.1 read (agent D2; V3).
- G. Craciun, F. Nazarov, C. Pantea, *Persistence and permanence of mass-action and power-law dynamical systems*, SIAM
  J. Appl. Math. 73 (2013) 305–329, DOI 10.1137/100812355 (Crossref), arXiv:1010.3050. Introduction and Theorems 6.4, 7.2 read (agent G).
- C. Pantea, *On the persistence and global stability of mass-action systems*, SIAM J. Math. Anal. 44 (2012) 1636–1673,
  DOI 10.1137/110840509, arXiv:1103.0603. Introduction, Definition 3.5, Theorem 5.1, §6 read (agents G, Q1).
- D. F. Anderson, *A proof of the global attractor conjecture in the single linkage class case*, SIAM J. Appl. Math. 71
  (2011) 1487–1508, DOI 10.1137/11082631X. Introduction, Theorem 4.10, Corollary 4.11 read (agent G; V3).
- D. F. Anderson, A. Shiu, *The dynamics of weakly reversible population processes near facets*, SIAM J. Appl. Math. 70
  (2010) 1840–1858, DOI 10.1137/090764098. Theorems 3.2, 3.4, 4.6, Corollary 4.7 read (agent G; V3).
- M. Gopalkrishnan, E. Miller, A. Shiu, *A geometric approach to the global attractor conjecture*, SIAM J. Appl. Dyn.
  Syst. 13 (2014) 758–797, DOI 10.1137/130928170. Definitions, Theorems 1.1–1.2, Def. 3.14 read (agent G; V3).
- G. Craciun, A. Dickenstein, A. Shiu, B. Sturmfels, *Toric dynamical systems*, J. Symb. Comput. 44 (2009) 1551–1565,
  DOI 10.1016/j.jsc.2008.08.006. §1–3, §5 read (agent L): the source's reference [4] for "small dimensions" (Theorem
  23); Theorem 9 (codimension δ of the toric locus).

## Global Attractor Conjecture, status in October 2026 (not checked)

- G. Craciun, *Toric differential inclusions and a proof of the global attractor conjecture*, arXiv:1501.02860, v3
  (23 September 2026), a major revision. Read by agents L and G: introduction, theorem statements, proof strategy,
  §7–8 construction statements. The proof was not checked.
- OpenAI, *Boundedness and persistence of weakly reversible mass-action systems*, preprint of 25 September 2026 in the
  `openai/math` repository (folder `preprints/`), with a Lean statement of its main theorem. Theorem 1.1, Remark 3.5 and
  the Lean statement read by agents L and G; not checked.
- C. Wiuf, *A proof of the Global Attractor Conjecture in a special case*, arXiv:2609.24553 (21 September 2026).
  Read in full by agent L; §2–3, 6–7, Theorem 7.2 by agent G; the reading of condition (⋄) checked on the PDF by V3.
  Not checked as a proof.

## Novelty checks

### 25 October 2026 (archived report)

- OpenAlex: no citing work of the source; Semantic Scholar rate-limited; Crossref: no published version.
- The author's later papers: arXiv:2609.12636 (unrelated).
- AI campaigns (`gh api users/ipitchford/repos`, `users/mmislan/repos`), Zenodo, `gh search repos`/`code` on
  "preclusters", "invariant polyhedral support": nothing relevant. arXiv search: only the source; two unrelated 2026
  preprints use "duality" in another sense (Banaji arXiv:2609.07265; Beldiman–Halanay arXiv:2610.27188).

### 27 October 2026 (agents L, effort xhigh, and L2, effort medium; `gac-duality/literature/notes.md`,
`gac-duality/literature2/novelty.md`)

- Primary sources above, full texts where stated; OpenAlex and Semantic Scholar citing works of the source (none),
  Blanchini–Giordano (7, none on duals beyond the source), Craciun's preprint, Wiuf's preprint, Conradi–Flockerzi (41
  titles; none on clusters as kernel classes, Petri nets or duals).
- arXiv API: "reverse-dual" (11 hits, only Jančar et al. relevant), "equal conflict" (0), siphon ∧ dual ∧ "reaction
  network" (only the source), drain ∧ "reaction network" ∧ persistence (0), "reaction network" ∧ transpose (0), "dual
  network" ∧ "mass-action" (0), "free-choice" ∧ reaction (only Mairesse–Nguyen), "global attractor conjecture",
  "preclusters" (only the source).
- AI campaigns of PROTOCOL.md §1 (Palomar/mmislan, Evidence Press/ipitchford, the OpenAI mathematics repository, Zenodo,
  GitHub search on "reverse-dual", "equal conflict", "drain species", "dual reaction network", "transpose
  stoichiometric", "Iosif duality", "preclusters", "siphon persistence"): nothing on duals, reverse duals, equal-conflict
  dynamics or drain sets.
- Google Scholar was not reachable. Prior art found and credited in the paper: Famili–Palsson 2003, Jamshidi–Palsson 2009 and Duvall 2024 (networks with stoichiometric matrix ±Γᵀ); Teruel–Silva p. 293 and Cor. 32; Desel–Esparza Thm 4.30(a), Prop. 5.4,
  Thms 10.16–10.17 (Desel 1994); Teruel–Silva; Hack; Jančar et al.; Mairesse–Nguyen; Pérez Millán et al.;
  Campos–Chiola–Silva; Recalde–Júlvez–Silva and Júlvez–Recalde–Silva; Iosif's thesis Theorem 4.17;
  Zhang–Gao–Dochain Theorem 3.10; Blanchini–Giordano p. 493; Butler–Waltman, Hale–Waltman; the max/min consensus
  technique. Other works examined and found not to contain the results: Donnell–Banaji arXiv:1211.2153, node-balanced
  steady states (Feliu–Cappelletti–Wiuf arXiv:1711.07074), Agarwal–Craciun–Deshpande–Jin arXiv:2601.08801,
  Deshpande–Gopalkrishnan arXiv:1309.3957, Feliu arXiv:2610.11679, arXiv:2610.30808.

Conclusion: no work cites or answers Conjecture 22 or Question 24 apart from this project; "not found" in the paper's
table of results refers to these searches.
