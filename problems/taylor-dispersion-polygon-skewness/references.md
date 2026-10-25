# References and novelty check

Retrieved 24 October 2026. Third-party PDFs and texts are stored locally in `certificates/sources/` (git-ignored).

## Primary source

1. M. Aminian, R. Camassa, R. M. McLaughlin, *Mass distribution and skewness for passive scalar transport in pipes
   with polygonal and smooth cross sections*, **Stud. Appl. Math. 141(3) (2018) 399–417**,
   doi:[10.1111/sapm.12230](https://doi.org/10.1111/sapm.12230) (received 23 Feb 2018, revised 5 Aug 2018,
   accepted 9 Aug 2018, online 24 Sep 2018). **Read in full** (publisher PDF,
   `certificates/sources/supplied-aminian2018/`). Quoted: abstract (p. 399), introduction (p. 401), eqs. (4),
   (12)–(18) (pp. 402–406), Section 4 and Table 1 (pp. 406–407), Conclusion (p. 416).
2. Same authors, arXiv:[1708.05097](https://arxiv.org/abs/1708.05097), v1 (16 Aug 2017), v2 (18 Aug 2017),
   v3 (5 Feb 2018), physics.flu-dyn / math.AP. **Read in full** (all three versions; `certificates/sources/aminian2017/`).
   Differences v1 → v2: wording "long time solute skewness" → "long time solute negative skewness", a summary
   paragraph removed, figures moved. v2 → v3: abstract adds "both non-zero short time skewness and" (racetracks);
   the regular-polygon statement and Table 1 unchanged. arXiv → published: see result.md, Section 1.

## Earlier papers of the same group (definitions of the functionals)

3. M. Aminian, F. Bernardi, R. Camassa, R. M. McLaughlin, *Squaring the circle: geometric skewness and symmetry
   breaking for passive scalar transport in ducts and pipes*, Phys. Rev. Lett. 115 (2015) 154503,
   doi:[10.1103/PhysRevLett.115.154503](https://doi.org/10.1103/PhysRevLett.115.154503). **Not read** (publisher
   page and accepted manuscript returned HTTP 403); its definition of the geometric skewness `<u^3>/<u^2>^{3/2}`
   is restated in source 1, eq. (12), which is what we use.
4. M. Aminian, F. Bernardi, R. Camassa, D. M. Harris, R. M. McLaughlin, *How boundaries shape chemical delivery in
   microfluidics*, Science 354 (2016) 1252–1256, doi:[10.1126/science.aag0532](https://doi.org/10.1126/science.aag0532).
   **Not read** (supplementary material, which contains the long-time derivation, returned HTTP 403). The long-time
   formula is restated in source 1 (eqs. (15)–(17)); we re-derived it independently (result.md, Proposition 1).
5. R. Aris, *On the dispersion of a solute in a fluid flowing through a tube*, Proc. R. Soc. Lond. A 235 (1956)
   67–77, doi:[10.1098/rspa.1956.0065](https://doi.org/10.1098/rspa.1956.0065). Not read (classical moment
   hierarchy, restated in source 1 and used as in the source).

## Mathematical tools cited

6. M. Bebendorf, *A note on the Poincaré inequality for convex domains*, Z. Anal. Anwend. 22 (2003) 751–756,
   doi:[10.4171/ZAA/1170](https://doi.org/10.4171/ZAA/1170). **Read** (open access PDF): Theorem 3.2,
   `||u||_{L2} <= (d/pi) ||grad u||_{L2}` for mean-zero `u in H^1` on a convex domain of diameter d.
7. L. E. Payne, H. F. Weinberger, *An optimal Poincaré inequality for convex domains*, Arch. Rational Mech. Anal. 5
   (1960) 286–292, doi:[10.1007/BF00252910](https://doi.org/10.1007/BF00252910). Not read (cited through 6; the
   original proof is valid in dimension 2 according to 6).
8. P. Grisvard, *Elliptic Problems in Nonsmooth Domains*, Pitman 1985; SIAM Classics 2011,
   doi:[10.1137/1.9781611972030](https://doi.org/10.1137/1.9781611972030), Theorem 3.2.1.2 (H^2 regularity of the
   Dirichlet problem on bounded convex domains). **Not re-read**; statement recalled and consistent with its use in
   arXiv:2408.05461 (Theorem 6.1 there). Listed in NEEDS.md.

## Novelty check (citing works 2018–2026, later papers of the authors, AI campaigns)

* Semantic Scholar citation lists (OpenAlex rate-limited on 24 Oct 2026): 8 works cite source 1, 22 cite source 3,
  58 cite source 4 (lists in `certificates/sources/openalex/s2_*.json`). Titles/abstracts screened. None proves the
  sign of the long-time or the geometric skewness for regular polygons. Relevant ones:
  * J. Feng, X. Chu, *Transient and asymptotic Taylor–Aris dispersion of Brownian rods in arbitrary
    regular-polygonal ducts*, arXiv:[2605.22982](https://arxiv.org/abs/2605.22982) (May 2026). **Full text
    searched**: concerns the effective (Taylor) diffusivity of rods; the word "skew" appears only for initial
    profiles. No result on the skewness sign.
  * R. Camassa, Z. Kılıç, R. M. McLaughlin, Physica D (2019), arXiv:1802.03340 (random passive scalar) — title/abstract.
  * L. Ding, R. Hunt, R. M. McLaughlin, H. Woodie, Res. Math. Sci. (2021), arXiv:2008.10717 (oscillating wall) — title/abstract.
  * L. Ding, R. M. McLaughlin, Phys. Rev. Fluids 7 (2022) 074502, arXiv:2109.05617 — title/abstract.
  * L. Ding, J. Fluid Mech. (2025), arXiv:2504.14733 (periodically modulated channels) — title only.
  * A. R. Drumm, F. Bernardi, arXiv:2508.18262 (membrane filters) — title only.
  * F. Bernardi et al., APS DFD 2022 talk Q26.1 (isosceles triangles and trefoils; abstract via web search
    summary only): numerical/asymptotic, isosceles triangles change sign.
* arXiv API searches (author Camassa/McLaughlin/Aminian/Ding + "skewness"; "geometric skewness"; "skewness" +
  "polygon" + "dispersion"): no further relevant item.
* AI campaigns: GitHub `mmislan` (27 Palomar repositories, all reaction networks/Lean) and `ipitchford` (no
  repository matching "skew" or "dispersion"); Zenodo searches ("taylor dispersion skewness", "geometric skewness",
  "passive scalar skewness polygon", "Aminian Camassa McLaughlin"): nothing relevant. Web searches: nothing.
* Conclusion: to our knowledge the sign claims for regular polygons with n >= 4 were not proved before.
