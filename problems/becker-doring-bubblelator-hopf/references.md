# References and novelty check

## Primary source

- B. Niethammer, R. L. Pego, A. Schlichting, J. J. L. Velázquez, *Oscillations in a Becker–Döring model with
  injection and depletion*, SIAM J. Appl. Math. 82(4) (2022) 1194–1219, doi:10.1137/20M1398664,
  arXiv:2102.06751 (v2, 27 September 2021).
  - **Read in full**: the arXiv v2 LaTeX source and PDF, including the appendices, and the accepted manuscript on
    NSF PAR (par.nsf.gov/servlets/purl/10385863), whose open sentences are identical. The typeset SIAM version is
    paywalled and was not seen.
  - Used: limit model (1.12)–(1.14); §1.4 (open sentence, p. 4); §3 (open sentence, p. 13; change of variables
    (3.1)–(3.4); delay equation (3.7); steady state (3.8); characteristic equation (3.9)–(3.11); the Lemma of §3.1;
    criteria (3.13)–(3.15); direction formulas (3.17)–(3.23); nonresonance (3.28)); §4.1 (β = 0); §4.2 (Figures 3–4);
    §4.3 (odd β, Table 1); Appendix C (Tables 2–3).
  - Data files of the arXiv source, used only for comparison by `scripts/05_compare_source.py`:
    `BDoscillate-theta-kappa.csv`, `BDoscillate-resonance-transversality.csv`, `BDoscillate-odd_beta3.csv`
    (downloaded into `certificates/sources/`, which is git-ignored).

## Background used in the proofs

- M. G. Crandall, P. H. Rabinowitz, *The Hopf bifurcation theorem in infinite dimensions*, Arch. Rational Mech. Anal.
  67 (1977) 53–72, doi:10.1007/BF00280827. Cited for the Lyapunov–Schmidt approach to Hopf bifurcation; our
  proof is self-contained (not re-read).
- F. Johansson, *Numerical integration in arbitrary-precision ball arithmetic*, ICMS 2018, LNCS 10931, 255–263,
  doi:10.1007/978-3-319-96418-8_30, arXiv:1802.07942. The rigorous integrator `acb_calc_integrate` and the
  requirement on branch cuts (the "analytic" flag); documentation of python-flint 0.9 read.
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. 66
  (2017) 1281–1292, doi:10.1109/TC.2017.2690633 (cited, not re-read).

## Background cited for context (not used in the proofs)

- O. J. Staffans, *Hopf bifurcation for functional and functional differential equations with infinite delay*,
  J. Differential Equations 70 (1987) 114–151, doi:10.1016/0022-0396(87)90171-9. A general Hopf theorem with
  fading memory; found by web search, abstract only. It shows that the source's sentence of §1.4 ("only for finite
  delay horizons") is too strong as a statement about general theory; the source's sentence of §3 concerns "a model
  of this type", which is what we prove.
- Y. Hino, S. Murakami, T. Naito, *Functional Differential Equations with Infinite Delay*, LNM 1473, Springer 1991,
  doi:10.1007/BFb0084432 (cited, not read).
- O. Diekmann, M. Gyllenberg, *Equations with infinite delay: blending the abstract and the concrete*, J. Differential
  Equations 252 (2012) 819–851, doi:10.1016/j.jde.2011.09.038 (cited, not read).
- O. Diekmann, S. A. van Gils, S. M. Verduyn Lunel, H.-O. Walther, *Delay Equations*, Springer 1995,
  doi:10.1007/978-1-4612-4206-2 (cited by the source; not read).

## Novelty check (5 November 2026)

- **Citing works.** Semantic Scholar (7 citing works) and OpenAlex (8), titles and abstracts read:
  Kalinov et al. (J. Comput. Phys. 2022); Fortin (J. Phys. A 2022, Smoluchowski system with source); Fortin–Choi
  (J. Phys. A 2023, shattering models); Ferreira–Franco–Lukkarinen–Nota–Velázquez (J. Phys. A 2023); a cloud-chamber
  paper (Oxford Open Climate Change 2024); Matveev et al. (Comput. Math. Math. Phys. 2022); Ferreira–Lukkarinen–Nota–
  Velázquez (SIAM J. Math. Anal. 2024); Doumic–Fellner–Mezache–Velázquez (Nonlinearity 2025); S. Loin,
  arXiv:2508.07884 (injection and irreversible fragmentation; well-posedness and convergence, no Hopf theorem; read
  for the candidate list); P. Krapivsky, arXiv:2605.22383 (aggregation–fragmentation steady states);
  Deuschel–Köppl–Steenbeck–Worschech, arXiv:2609.07622 (Ornstein–Uhlenbeck processes). None treats the Hopf
  bifurcation of the bubblelator limit model.
- **Authors' later papers.** arXiv listings of R. L. Pego (to 2025-11) and A. Schlichting (to 2026-08): no
  follow-up on this model.
- **arXiv** full-text-index searches "bubblelator" (only the source) and "Becker–Döring" with "oscillations"
  (the source and the bi-monomeric paper).
- **Zenodo** API: "Becker-Döring Hopf", "bubblelator", "Becker-Döring oscillation": nothing relevant.
- **GitHub**: `gh search repos` "Becker-Doring Hopf" and "bubblelator", `gh search code` "bubblelator" (only this
  repository's candidate list); `ipitchford` and `mmislan` repositories: no matching names.
- **Web search** on the key phrases: only the source and Staffans's general theorem.

Conclusion: open as of 5 November 2026; no rigorous Hopf theorem for this model was found.
