# References

## Primary source

- R. L. Pego, J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
  Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605.
  - **Read in full** (arXiv v2, LaTeX source and PDF).
  - Model (1)–(4); matrix \(B\) (27); Lemma 3.1; Theorems 3.2–3.3 and the open sentence (§3.3); Table 1.
  - Continuum approximation (§4); characteristic function \(F\) (71)–(76); Prop. 5.3 and Lemma 5.4.
  - \(Q(z;\kappa)\) (87); Hopf conditions (89)–(91); Lemmas 7.1–7.3; convergence results of §§8–10.
  - The novelty sub-agent checked the accepted manuscript on IOPscience: the open sentence is identical in arXiv
    v1, v2 and the accepted manuscript. The typeset version is paywalled and was not seen.

## Novelty check (citing works, 2020–2026)

Done on 2026-10-18 by a sub-agent: OpenAlex (17 citing works), Semantic Scholar (16), 20 distinct. Every arXiv full
text was grepped for "supercritical", "subcritical", "Lyapunov", "stable periodic" and "atomization". None
determines the criticality.
- **Closest.** B. Niethammer, R. L. Pego, A. Schlichting, J. J. L. Velázquez, *Oscillations in a Becker–Döring
  model with injection and depletion*, SIAM J. Appl. Math. 2022, arXiv:2102.06751. It has a formal
  direction-of-bifurcation computation ("we expect a supercritical Hopf bifurcation") for a *different* model.
  Checked by the sub-agent; I did not read it myself.
- **Other citing works** (title and abstract, or full-text grep). They treat different models or cite the paper in
  passing:
  - Budzinskiy–Matveev–Krapivsky (PRE 2021);
  - Kalinov et al. (JCP 2022);
  - Fortin (J. Phys. A 2022);
  - Fortin–Choi (J. Phys. A 2023);
  - Matveev (review, arXiv:2606.10785, 2026), which states that stability of periodic solutions is still open for
    its own model;
  - Doumic–Fellner–Mezache–Velázquez (Nonlinearity 2025);
  - Fagan–MacKay–Wood (EPJB 2024; EPL 2021);
  - Krapivsky (arXiv:2605.22383);
  - Deuschel–Köppl–Steenbeck–Worschech (arXiv:2609.07622);
  - Krapivsky–Matveev (PRE 2024);
  - Ferreira–Lukkarinen–Nota–Velázquez (SIMA 2024);
  - Franco–Kepka–Velázquez (ARMA 2026);
  - Timokhin et al. (2022);
  - Matveev et al. (CMMP 2022);
  - da Costa–Drmota–Grinfeld (2020);
  - Mezache (thesis 2019);
  - Calvo–Hingant–Yvinec (KRM 2023) and Sekiyama–Ohtsuki–Yamamoto (Physica A 2022), abstracts only.
- **Limits.** Google Scholar was not queried.

## Methods

- Yu. A. Kuznetsov, *Elements of Applied Bifurcation Theory*, 3rd ed., Springer 2004, §3.5: the projection formula
  for the first Lyapunov coefficient in \(\mathbb R^n\). Used from memory and validated numerically (planar test
  against the Guckenheimer–Holmes formula; direct simulation).
- J. Guckenheimer, P. Holmes, *Nonlinear Oscillations, Dynamical Systems, and Bifurcations of Vector Fields*,
  Springer 1983, §3.4: the planar formula for the Hopf coefficient, used in test V1.
- R. Krawczyk, *Newton-Algorithmen zur Bestimmung von Nullstellen mit Fehlerschranken*, Computing 4 (1969) 187–201;
  R. E. Moore, *A test for existence of solutions to nonlinear systems*, SIAM J. Numer. Anal. 14 (1977) 611–615.
  The Krawczyk operator: standard, not re-read.
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. 66
  (2017) 1281–1292. Used through python-flint 0.9, including `acb_mat.eig` (rigorous eigenvalue enclosures), not
  re-read.
- R. Veltz, *BifurcationKit.jl* (v0.8.7). Independent Hopf normal form computation; documentation only.
- O. Diekmann, Ph. Getto, M. Gyllenberg, *Stability and bifurcation analysis of Volterra functional equations in
  the light of suns and stars*, SIAM J. Math. Anal. 39 (2008) 1023–1069. Cited for the rigorous Hopf theory of
  renewal equations (§7.2 of notes.md). **Not read** for this work.

## Completion audit (25–26 October 2026)

The complete branches 1--3 proof and refreshed primary-source/method checks are documented in
[all-M/references.md](all-M/references.md) and [all-M/verification.md](all-M/verification.md).
The earlier source and novelty record is retained above.
