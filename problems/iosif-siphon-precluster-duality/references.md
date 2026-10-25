# References and novelty check

## Primary source (full text read)

- A. Iosif, *Duality in mass-action networks*, arXiv:2603.08767. Read in full: v1 (9 March 2026) and v2 (15 July 2026).
  The mathematical content is identical; v2 adds Remark 7 and minor rewording. HAL hal-03058670 v3 (17 July 2026) is
  identical to arXiv v2. Quoted: Conjecture 22, Remark 23, Question 24 (§3, p. 7 of v2), the abstract and p. 2.
  Definitions copied: 13–14 (siphons), 17 and Remark 18 (MIPS), 19 (preclusters), 20 (internal cycles), the dual
  network (§3), Example 25 and the final diagram.
- A. Iosif, *Duality in mass-action networks*, HAL hal-03058670v1 (11 December 2020). Read in full (4 pages). It
  states "Theorem. For conservative mass-action networks that do not have two species with exactly the same rates the
  set of preclusters is dual to the set of maximal invariant polyhedral supports", without proof and without the
  positive-steady-state hypothesis. It also asks: "Is there a duality relation between clusters and siphons?"
- A. Iosif, *Algebraic methods for the study of multistationarity in mass-action networks*, PhD thesis,
  Otto-von-Guericke-Universität Magdeburg (2019), DOI 10.25673/13949. The repository page is behind a bot check. The
  PDF was retrieved from the Internet Archive copy of
  `opendata.uni-halle.de/bitstream/1981185920/14078/1/Iosif_Alexandru_Dissertation_2019.pdf` (snapshot of 22 February
  2020). Read: §1.6 (internal cycles, Proposition 1.19, Example 1.22) and Chapter 4 (doubling sets, Definition 4.3 and
  Algorithm 4.4 for preclusters, Definition 4.6 for clusters, Examples 4.1–4.14).
- A. Iosif, *Multistationarity in systems with toric states*, seminar slides, 19 November 2024,
  https://gasiull.webs.ull.es/beamer-2024-11-19.pdf. Read: slides 28–29 ("Work in progress 2: a duality theory for
  mass-action networks (together with Lamprini Ananiadi)"). They state "Conjecture: Siphons and clusters are dual
  objects" and "Theorem (Ananiadi, I.) (Evidence for the conjecture): Maximal invariant polyhedral supports (MIPS)
  ... are dual to preclusters". Preclusters are described there as a "partition of the arrow set collecting relations
  between the coordinates of the cone ker(Ap − Ae) ∩ R^r_{≥0} (INFORMAL)".
- A. Iosif, L. Ananiadi, poster P1.H, Jóvenes Investigadores RSME 2025 (Bilbao). Only the abstract was read, via the
  conference web pages (GitHub repository `tonellicueto/jovenesRSME2025`, file `posteres_en.html`): "we prove that
  preclusters are dual to maximal invariant polyhedral supports".

## Background (cited for standard facts)

- C. Conradi, D. Flockerzi, *Multistationarity in mass action networks with applications to ERK activation*, J. Math.
  Biol. 65 (2012) 107–156, DOI 10.1007/s00285-011-0453-1. Publisher PDF read.
  Read: §2 (network N1, the two-site phosphorylation network), §3.1–3.2 (Assumption 1, the set Λ(E), Assumption 2,
  Fact 2, Algorithm 1, Remark 5, Example 1, Lemma 1, Remark 6) and Definition 3 (Isolation and Bridging Property).
  Comparison with the thesis definitions used here:
  - Doubling vectors e_l − e_k ∈ ker Ye (Assumption 2) are the thesis's doubling sets.
  - Fact 2(a) is exactly the span rule of the paper's Definition 19, applied once from a doubling pair.
  - Assumption 1(ii) (no zero row of E) corresponds to our use of consistency.
  - Clusters (Fact 2(b), Algorithm 1, Remark 6(5)) group the rows of E with equal ratio n_jλ/n_jμ on the pairs
    (λ, μ) ∈ Λ(E)² that satisfy the doubling relations (3.8c). They do not involve the rate constants, whereas the
    thesis fixes k (its Remark 4.9).
  - The original has no notion of precluster and no analogue of Algorithm 4.4; the thesis calls its algorithm "very
    similar" to Algorithm 1, which groups rows by equal ratios rather than building the span graph.
  - In Example 1, the clusters of N1 are {1,2,3,10,11,12} and {4,…,9}. These coincide with the preclusters computed
    by `motifs.py` for the same network, which uses the same reaction order.
  - For N3, the original definition also gives the single cluster {R1,R2,R3} (manuscript, Remark 13).

  Nothing in our definitions or results changes.
- S. Schuster, C. Hilgetag, J. Biol. Syst. 2 (1994) 165–182, DOI 10.1142/S0218339094000131 (elementary flux modes =
  internal cycles). Not re-read.
- A. Shiu, B. Sturmfels, *Siphons in chemical reaction networks*, Bull. Math. Biol. 72 (2010) 1448–1463, DOI
  10.1007/s11538-010-9502-y. Not re-read; cited with the definition of siphons, which was taken from the source.
- D. Angeli, P. De Leenheer, E. D. Sontag, Math. Biosci. 210 (2007) 598–618, DOI 10.1016/j.mbs.2007.07.003 (siphons and
  conservation laws). Not re-read; cited for the standard fact that the support of a semi-positive conservation law is
  a siphon.
- A. P. Burgard et al., *Flux coupling analysis*, Genome Res. 14 (2004) 301–312, DOI 10.1101/gr.1926504; T. Pfeiffer
  et al., *METATOOL*, Bioinformatics 15 (1999) 251–257, DOI 10.1093/bioinformatics/15.3.251 (full coupling, enzyme
  subsets). Not re-read; cited for terminology.
- M. Feinberg, *Foundations of Chemical Reaction Network Theory*, Springer 2019, DOI 10.1007/978-3-030-03858-8;
  F. Horn, R. Jackson, Arch. Rational Mech. Anal. 47 (1972) 81–116, DOI 10.1007/BF00251225 (deficiency, complex
  balancing). Standard facts, not re-read.
- F. Blanchini, G. Giordano, IEEE Control Syst. Lett. 6 (2022) 488–493, DOI 10.1109/LCSYS.2021.3081369 (dual networks,
  mentioned by the source). Not read.
- J. A. De Loera, E. D. Kim, S. Onn, F. Santos, J. Combin. Theory Ser. A 116 (2009) 1306–1325, DOI
  10.1016/j.jcta.2009.03.010 (cited by the source for the finiteness of the combinatorial types of invariant polyhedra). Not read; cited
  only as the source cites it.

## Novelty check (25 October 2026)

- **Citing works.** OpenAlex lists W4287859579 ("Duality in mass-action networks") with 0 citing works. The Semantic
  Scholar API was rate-limited (HTTP 429). The Crossref search found no published version.
- **The author's later papers.** arXiv lists one later preprint by A. Iosif, arXiv:2609.12636 (September 2026, on
  Sundman transformations of Painlevé equations), which is unrelated. There is no later duality paper.
- **AI campaigns.** `gh api users/ipitchford/repos` (Evidence Press, 193 repositories): the CRNT-related candidates
  (quadratic-mass-action-sharpness, reverse-reaction-attracting-oscillations, …) do not touch preclusters, MIPS or
  siphon duality. `gh api users/mmislan/repos` (Palomar, 27 repositories): no matching topic. Zenodo API on "preclusters",
  "maximal invariant polyhedral supports", "siphons preclusters" and "Iosif duality mass-action": 0 records.
  `gh search repos` and `gh search code` on "preclusters" and "invariant polyhedral support": only arXiv RSS mirrors of
  the abstract, the RSME 2025 poster abstract, this repository's candidate list, and unrelated physics code.
- **arXiv search** (`all:preclusters`, `all:"invariant polyhedral"`, `siphons AND duality`): only the source is
  relevant. Two other 2026 preprints mention duality in mass-action networks in a different sense. Banaji,
  arXiv:2609.07265, uses Gale duality for fold and cusp bifurcations. Beldiman and Halanay, arXiv:2610.27188, use a
  "vacancy duality" for canard points. Neither addresses preclusters or MIPS.
- **Web search** on the key phrases: nothing beyond the source.

Conclusion: the conjecture was open; no previous resolution was found.
