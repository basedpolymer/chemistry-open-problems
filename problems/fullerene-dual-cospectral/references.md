# References and novelty check

All third-party files are kept in `certificates/sources/` (git-ignored, not redistributed).

## Primary source (two versions with DIFFERENT conjectures)

1. **A. Bille, V. Buchstaber, E. Spodarev**, *Spectral clustering of combinatorial fullerene isomers based on their
   facet graph structure*, arXiv:1910.00846**v1** (2 Oct 2019; the only arXiv version). **Read in full**
   (`sources/bille2021/bille_arxiv.pdf`). Section 3, Conjecture 1, p. 9 (statement (a) on the dual graph T_n with
   the exception n = 44; Table 1 p. 9, MATLAB double precision, symbolic toolbox for suspected pairs).
2. **Same authors, same title**, J. Math. Chem. **59** (2021) 264–288, doi:10.1007/s10910-020-01193-4 (published
   online 18 Nov 2020; open access, CC BY 4.0). **Read in full** (`sources/bille2021_springer/bille2021_vor.pdf`,
   obtained through the publisher's open-access PDF link). Section 3, Table 1 (p. 274), Conjecture 1 (a)–(d)
   (p. 275); the computations are redone "in Python's int-precision". The conjecture differs from arXiv v1:
   (a) is about T_n^6, (b) asserts the *existence* of T_n-cospectral isomers for all n ≥ 60.
3. A. Bille, poster "Spectral clustering of combinatorial fullerene isomers ..." (CONVIE 2022, Seville),
   `https://congreso.us.es/convie/wp-content/uploads/2022/07/bille.pdf` — **read in full**: reproduces the VOR
   Table 1 and the VOR Conjecture (four bullets). Abstract at SGSIA 2021 (Prague),
   `https://sgsia21.karlin.mff.cuni.cz/abstracts/Bille.pdf` — read (abstract only, no conjecture).

## Later work of the same authors (checked for updates of the conjecture)

4. **A. Bille, V. Buchstaber, E. Spodarev**, *Some open mathematical problems on fullerenes*, arXiv:2410.19322v1
   (25 Oct 2024), published J. Chem. Inf. Model. **65** (2025) 2911–2923, doi:10.1021/acs.jcim.4c01997 (CC BY 4.0).
   **arXiv v1 read in full** (`sources/bille2024/`). Section 5.1, Conjecture 3 and Open Problem 4 (p. 13):
   "For α, β ∈ (0, 1) and feasible n1, n2, let Tn1, Tn2 be two dual fullerenes. It holds that
   Tn1 ≃ Tn2 ⇔ chα,β(Tn1) = chα,β(Tn2)", with ch_{α,β}(T) = tr exp(αA + βD). It does not mention the 2021
   Conjecture 1 or cospectral pairs. **The JCIM version could not be downloaded (ACS returns 403): see NEEDS.md.**
5. A. Bille, V. Buchstaber, S. Coste, S. Kuriki, E. Spodarev, *Random eigenvalues of graphenes and the
   triangulation of plane*, arXiv:2306.01462, doi:10.1088/1751-8121/ad9efb — title/abstract only (not about
   cospectrality of finite fullerene duals).
6. A. Bille, V. Buchstaber, P. Ievlev, S. Novikov, E. Spodarev, *Random eigenvalues of nanotubes*,
   arXiv:2408.14313, doi:10.1088/1751-8121/adb54a — title/abstract only.

## Fowler's announcement

7. **P. W. Fowler**, *Cospectrality and the fullerenes*, abstract, Canadian Mathematical Society Winter Meeting
   2005, graph theory session, `https://w7.cms.math.ca/Events/winter05/abs/pdf/gt-pf.pdf` — **read in full**
   (one paragraph). Verbatim: "As yet, no cospectral fullerenes are known, but many pairs of fullerenes with
   cospectral duals have now been found (P. W. Fowler and M. Duncan, to be published). The dual of a fullerene
   graph is also a polyhedral graph, having all faces triangular, 12 vertices of degree 5 and all others of degree
   6; ... A construction for (infinite) families of cospectral fullerene duals will be presented."
   Fowler's notion (adjacency spectrum of the triangulation dual to the fullerene) is exactly the T_n spectrum of
   Bille et al. **No publication by Fowler and Duncan was found** (web searches, Crossref title/author queries,
   Semantic Scholar); the announcement contains no explicit pair.

## Citing works (Semantic Scholar, 24 Oct 2026; OpenAlex budget was exhausted)

Six citing works of ref. 2: refs. 4, 5, 6 above; X. Liu et al. (Sci. Total Environ. 2024, adsorbents, irrelevant);
D. Stepenshchikov, S. Aksenov (J. Struct. Chem. 2022, symmetry groups, irrelevant); A. Egorov, A. Vesnin
(Comput. Math. Biophys. 2020, hyperbolic volumes, irrelevant). None reports cospectral fullerene duals.

## Other sources consulted

8. P. Schwerdtfeger, L. Wirz, J. Avery, *The topology of fullerenes*, WIREs Comput. Mol. Sci. 5 (2015) 96–145,
   doi:10.1002/wcms.1207 (Europe PMC full text, searched for "cospectr"/"isospectr": no hit).
9. G. Brinkmann, J. Goedgebeur, B. D. McKay, *The generation of fullerenes*, J. Chem. Inf. Model. 52 (2012)
   2910–2918, doi:10.1021/ci3003107, arXiv:1207.7010 — buckygen (software used, version 1.1,
   `https://caagt.ugent.be/buckygen/`, zip SHA-256 c151b330…e634).
10. B. D. McKay, A. Piperno, *Practical graph isomorphism, II*, J. Symb. Comput. 60 (2014) 94–112,
    doi:10.1016/j.jsc.2013.09.003 — nauty 2.8.8 (labelg, geng), taken from the pynauty 2.8.8.1 source distribution.
11. OEIS A007894 (number of fullerenes with 2n vertices), b-file by J. Goedgebeur, used to validate the isomer
    counts.
12. P. W. Fowler, D. E. Manolopoulos, *An Atlas of Fullerenes*, Clarendon Press, Oxford 1995 — spiral algorithm
    and isomer numbering (not consulted directly; the numbering convention is the one used by Bille et al., and
    our implementation reproduces C60:1812 = I_h and the C44:37/38 pair of Bille et al.).
13. C. D. Godsil, B. D. McKay, *Constructing cospectral graphs*, Aequationes Math. 25 (1982) 257–268 — GM
    switching (tested on the small pairs: no GM switching with |C| ≤ 6 explains them).
14. The cospectrality of C6 ∪ K1 and the subdivided claw K_{1,3}^{(2)} (spectrum {±2, ±1, ±1, 0}) is a classical
    small example of a disconnected graph cospectral with a connected one (see e.g. D. Cvetković, P. Rowlinson,
    S. Simić, *An Introduction to the Theory of Graph Spectra*, CUP 2010, ch. 1 — not re-read here); we re-derive it
    exactly and do not rely on the citation.

## AI campaigns and repositories

- google-deepmind/formal-conjectures (shallow clone, 24 Oct 2026): no file mentions fullerenes or cospectrality.
- openai/math (shallow clone, 24 Oct 2026): no relevant item (the only "fullerene" hit is a bibliography entry on the
  Barnette–Goodey conjecture).
- GitHub user listings for `mmislan` and `ipitchford`: not reachable from this session (API path refused); Zenodo
  searches ("fullerene cospectral", "fullerene dual spectrum", "Bille Buchstaber Spodarev", "fullerene conjecture
  spectrum"): no relevant record.
- Web searches ("cospectral duals" fullerenes; isospectral fullerenes; Fowler Duncan cospectral; fullerene
  cospectral mates): only refs. 1–4 and 7.
