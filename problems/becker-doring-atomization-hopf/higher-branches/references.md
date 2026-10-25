# References (higher branches, extension of result 8)

## Primary source

- R. L. Pego, J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
  Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605.
  **Read**: the publisher PDF (local copy, not distributed), §§3, 5 and the statements of §§7–10. Used: model
  (1)–(4) and matrix $B$ (27); Lemma 3.1 (simple zero eigenvalue); Theorems 3.2, 3.3 and the open sentence of §3.3
  (verbatim, publisher version, p. 1821: "We have not managed to determine analytically whether the bifurcating
  solutions are stable (the supercritical case) or not."); the reduction to $F(\varphi)$ ((63)–(76)), Prop. 5.3
  (eigenvalues ↔ zeros of $F$, an equivalence), Lemma 5.4 (simplicity); the limit function $Q$ (87) and its imaginary
  roots (89)–(91); Lemma 7.3 (transversality of the limit roots).  The numbering of the publisher version matches the
  arXiv v2 numbering used by result 8 and all-M for every item above.

## This repository

- Result 8, `../README.md`, `../paper.tex`, `../verification.md`, `../scripts/bdlib.py` (closed form of $\ell_1$ at integer
  $M$, Krawczyk test at integer $M$; verified there): `bdlib.py` copied unchanged.  Its formal limit coefficient $L_j$
  (Proposition B) is proved here to be the limit of $\ell_1(q_1=1)/M$ (Section 3 of the README).
- Result 8's verifier, `../verification/v03_sparse_float.py`, `v04_sparse_arb.py`: imported (not modified) by
  `scripts/xcheck_sparse.py` as a separate numerical cross-check.
- All-M extension, `../all-M/result.md` and `../all-M/scripts/` (not yet verified): `ad.py`, `allm_core.py`,
  `l1reg.py`, `h2_lib.py`, `h2j_lib.py`, `l1_sweep.py`, `verify_sweep.py`, `h2j_sweep.py` copied unchanged.  Used: the
  $\varepsilon$-regular formulation, Lemma 2.1 (parametric Krawczyk), Lemmas 3.1, 3.3, 3.3′, 4.1, 4.2′, 4.4,
  Theorem 4.5, Section 5 (regular formula for $\ell_1$) and Lemma 5.2.

## Methods (standard, not re-read)

- Yu. A. Kuznetsov, *Elements of Applied Bifurcation Theory*, 3rd ed., Springer 2004, doi:10.1007/978-1-4757-3978-7:
  projection formula for $\ell_1$, Hopf theorem with centre-manifold reduction (used as in result 8).
- R. Krawczyk, Computing 4 (1969) 187–201; R. E. Moore, SIAM J. Numer. Anal. 14 (1977) 611–615: Krawczyk operator.
- F. Johansson, *Arb*, IEEE Trans. Comput. 66 (2017) 1281–1292 (python-flint 0.9).
- Argument principle and Rouché-type counting in ball arithmetic: as in all-M Lemma 4.1.

## Novelty check (2026-10-25, short)

- OpenAlex, works citing doi:10.1088/1361-6544/ab6815: 17 records (2019–2026), the same list as result 8's
  `references.md` (latest: Characterizing the detailed balance property by means of measurements, ARMA 2026). None
  concerns the higher Hopf branches or their criticality.
- GitHub: `gh api users/ipitchford/repos --paginate` (193 repositories) and `users/mmislan/repos` (27): no repository on
  Becker–Döring equations or Hopf bifurcations.  `gh search repos "Becker-Doring"`: five simulation codes (2016–2026),
  none on this question.
- Zenodo API ("Becker-Döring" atomization; Becker Doring Hopf): no relevant record.
- arXiv API (all: Becker AND Doring): nothing relevant after result 8's check.
