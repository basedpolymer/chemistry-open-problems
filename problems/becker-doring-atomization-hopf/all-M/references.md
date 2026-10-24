# References (all-M extension of result 8)

## Primary source

- R. L. Pego, J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
  Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605.
  **Read in full** (arXiv v2, LaTeX source, downloaded to `certificates/sources/pv2020/`, git-ignored).
  Used: model (1)–(4); matrix $B$ (27); Lemma 3.1 (simple zero eigenvalue); Theorems 3.2–3.3 and the open sentence
  (§3.3, source lines 658–660: "We have not managed to determine analytically whether the bifurcating solutions are
  stable (the supercritical case) or not."); the reduction to $F(\varphi)$ (§5: (39)–(41), $f,g$, $\delta$,
  (76) $F=f_1g_2-g_1f_2$, the sorted representation (77) and the closed forms of $P_1,P_2$ ((e:P1-1), (e:P2-1))
  and $R_1,R_2$ ((e:r1), (e:r2)), all re-checked symbolically here, `check_formulas.py` F1); Prop. 5.1
  (eigenvalues ↔ zeros of $F$); Lemma 5.4 (simplicity; first step corrected by result 8's verifier); the limit
  function $Q(z;\kappa)$ (87) and its imaginary roots (89)–(91); Lemma 6.3 (Re λ ≥ 0 ⇒ Re φ ≥ 1); the structure
  (not the constants) of Props. 6.1, 6.4, 7.1, 7.2 and of the proof of Theorem 3.2 (§§6–9), whose constants
  $M_{0,k}$ are not explicit — replaced here by the explicit Lemma 3.3 of `result.md`. The typeset version
  (paywalled) was not seen; result 8's novelty sub-agent reported that the open sentence is identical in arXiv v1, v2
  and the accepted manuscript.

## Result 8 of this repository (the base of this extension)

- `problems/becker-doring-atomization-hopf/` (README.md, paper.tex, verification.md, scripts/, certificates/,
  verification/): read in full. Used: Theorem A (certified Hopf boxes and sign of $\ell_1$ for $15\le M\le10^5$,
  $j=1,2,3$), the verifier's H2 certificate for $M\le160$ (`acb_mat.eig`), the verifier's exact algebra for
  $3\le M\le57$, the formal limit coefficient $L_j$ (Proposition B). The certificate files
  `certificates/sweep_j*.csv.gz` are read by `scripts/check_vs_result8.py`.

## Novelty check (2026-10-24, bounded; result 8's full check of 2026-10-18 is not repeated)

- Semantic Scholar API, citations of arXiv:1905.02605: 16 citing works, the same list as in result 8's
  `references.md` (latest: Deuschel–Köppl–Steenbeck–Worschech, arXiv:2609.07622, unrelated). Nothing new.
- OpenAlex: rate limit exhausted (no API key) at the time of the check.
- arXiv API (searches "Becker AND Doring", "atomization AND Hopf", newest first): nothing relevant after
  2026-10-18.
- Zenodo API ("Becker-Döring" atomization): no relevant record.
- Web search ("Becker-Döring atomization Hopf bifurcation supercritical stable periodic solutions Pego Velázquez"):
  only the primary source and Niethammer–Pego–Schlichting–Velázquez (SIAM J. Appl. Math. 2022, arXiv:2102.06751,
  a different model, formal argument).
- GitHub user listings (mmislan, ipitchford): not accessible from this session (repository-scoped token); result 8's
  check of 2026-10-18 found nothing relevant there.

## Methods

- Yu. A. Kuznetsov, *Elements of Applied Bifurcation Theory*, 3rd ed., Springer 2004: projection formula (3.20) for
  the first Lyapunov coefficient; Hopf theorem and centre manifold reduction (§§3.4–3.5, 5.4). Used from memory, as
  in result 8 (where the formula was validated numerically). Not re-read; see NEEDS.md.
- R. Krawczyk, *Newton-Algorithmen zur Bestimmung von Nullstellen mit Fehlerschranken*, Computing 4 (1969)
  187–201; R. E. Moore, *A test for existence of solutions to nonlinear systems*, SIAM J. Numer. Anal. 14 (1977)
  611–615. Krawczyk operator (standard; the parametric version is proved in `result.md`, Lemma 2.1). Not re-read.
- S. Gershgorin, *Über die Abgrenzung der Eigenwerte einer Matrix*, Izv. Akad. Nauk SSSR 6 (1931) 749–754
  (Gershgorin discs, standard). Not re-read.
- Argument principle / winding number by interval arithmetic: standard; the version used is proved in `result.md`,
  Lemma 4.1.
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. 66
  (2017) 1281–1292; used through python-flint 0.9. Not re-read.
