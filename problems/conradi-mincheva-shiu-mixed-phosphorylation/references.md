# Sources and novelty audit

Audit date: 5 November 2026. Solver: Claude Opus 5.5 (effort xhigh).

## Primary source

Carsten Conradi, Maya Mincheva, Anne Shiu, *Emergence of oscillations in a mixed-mechanism phosphorylation system*,
Bull. Math. Biol. 81 (2019), no. 6, 1829–1852, DOI [10.1007/s11538-019-00580-6](https://doi.org/10.1007/s11538-019-00580-6);
arXiv [1809.02886](https://arxiv.org/abs/1809.02886) (v2, 1 February 2019).

Read: the full arXiv v2 PDF (24 pages) and the journal version, which is the publisher PDF deposited at NSF PAR,
<https://par.nsf.gov/servlets/purl/10090542>, also 24 pages. Used: network (1), ODEs (2), conservation laws (3), Table 1,
Prop. 3.1 (uniqueness), Prop. 3.2 (monomial parametrisation), Prop. 3.6, Prop. 4.1 (H1–H4 > 0, b6 > 0 when
k2 = k6 = k9), Thm. 4.2, inequality (13), Prop. 4.3, Thm. 4.5, Section 6.

* The **question** is unnumbered: arXiv v2 p. 20; journal version p. 21 of the PDF, i.e. p. 1849. It reads "Consider the
  mixed-mechanism network, and any choice of rate constants and total amounts. If the unique steady state x\* in 𝒫
  is locally stable, does it always follow that x\* is globally stable?" (checked on the page image: x\* is starred).
* The **conjecture** that an unstable steady state implies a stable periodic orbit is **Conjecture 6.1** in both the
  arXiv v2 (p. 18) and the journal version (p. 20 of the PDF, i.e. p. 1848). **Conjecture 6.2** is a different
  statement (oscillations for arbitrarily large K_tot, orbit in x5 independent of K_tot). The disagreement between
  sources mentioned in the campaign brief does not come from the paper, where the numbering is the same in both
  versions; it probably comes from the caption of Figure 5, which says "See Conjecture 6.2" about large K_tot.

## Other sources used

* T. Suwanmajo, J. Krishnan, J. R. Soc. Interface 12 (2015) 20141405, doi:10.1098/rsif.2014.1405. Cited through the
  source for uniqueness of steady states (not read).
* W. M. Liu, J. Math. Anal. Appl. 182 (1994) 250–256, doi:10.1006/jmaa.1994.1079, and X. Yang, Appl. Math. Lett. 15
  (2002) 615–621, doi:10.1016/S0893-9659(02)80014-3. Hurwitz-determinant criterion for a simple imaginary pair
  (statement as in the source's Prop. 2.3; papers not re-read).
* B. Hassard, N. Kazarinoff, Y.-H. Wan, *Theory and Applications of Hopf Bifurcation* (CUP 1981), Ch. 1; Yu. A.
  Kuznetsov, *Elements of Applied Bifurcation Theory*, 3rd ed. (Springer 2004), §§3.4–3.5, Ch. 5; Yu. A. Kuznetsov,
  SIAM J. Numer. Anal. 36 (1999) 1104–1124, doi:10.1137/S0036142998335005. Hopf theorem, stability of the bifurcating
  orbit, and the projection formula for l1 (standard results, cited from memory, not re-read).
* T. Kapela, M. Mrozek, D. Wilczak, P. Zgliczyński, CAPD::DynSys, Commun. Nonlinear Sci. Numer. Simul. 101 (2021)
  105578, doi:10.1016/j.cnsns.2020.105578; P. Zgliczyński, C¹ Lohner algorithm, Found. Comput. Math. 2 (2002) 429–465,
  doi:10.1007/s102080010025. Software: CAPD 6.1, github.com/CAPDGroup/CAPD (cloned 5 November 2026).
* R. Krawczyk, Computing 4 (1969) 187–201, doi:10.1007/BF02234767; A. Neumaier, *Interval Methods for Systems of
  Equations* (CUP 1990). The Krawczyk test.
* F. Johansson, Arb, IEEE Trans. Comput. 66 (2017) 1281–1292, doi:10.1109/TC.2017.2690633 (via python-flint 0.9).
* S. Rao, Math. Biosci. 301 (2018) 159–166, doi:10.1016/j.mbs.2018.05.001: the source's reference [36]
  (Michaelis–Menten limit). Abstract only.

DOIs were checked against Crossref on 5 November 2026.

## Novelty search

* **Citing works** (OpenAlex, work W2891313905, 23 citing records, all titles screened). The relevant ones are
  Conradi–Feliu–Mincheva (sequential distributive, no Hopf, 2019/2020), Obatake–Shiu–Tang–Torres (ERK model, 2019),
  Pantea–Romanovski (Vestnik SPbU 2020, absence of Hopf in two MAPK models, abstract read), Conradi–Mincheva (cyclic
  distributive, J. Math. Biol. 2024, abstract read), and Vassena (SIAM J. Appl. Dyn. Syst. 2025 and DCDS-B 2023,
  abstracts read). None computes Lyapunov coefficients for the mixed-mechanism network or answers the question.
* **Named in the brief.** Vassena, arXiv:2608.27081 (sequential distributive network, no Hopf under mass action;
  abstract read); Feliu–Kaihnsa, arXiv:2405.16179 (reduced fully distributive networks, absence of Hopf; abstract
  read); Zenodo 10.5281/zenodo.23446981 (Evidence Press, "A reverse reaction can create attracting oscillations",
  18 Oct 2026; record metadata and GitHub README read). That record says it does not resolve the motivating
  dual-phosphorylation question, which is in any case the Feliu–Kaihnsa question on the fully distributive network.
* **Zenodo API** queries "mixed-mechanism phosphorylation", "mixed mechanism phosphorylation oscillations",
  "Conradi Mincheva Shiu", "processive distributive Hopf subcritical", "dual phosphorylation Lyapunov coefficient":
  no relevant record.
* **GitHub.** `gh search repos` ("mixed-mechanism phosphorylation", "Conradi Mincheva Shiu", "mixed mechanism
  oscillations phosphorylation") returned no repositories. `gh search code` ("mixed-mechanism", "Conradi Mincheva
  Shiu", "mixed mechanism network" + Hopf, the ODE monomials) returned no relevant file. The ipitchford account has
  251 repositories; the only phosphorylation-related ones are `gk-shared-phosphatase-connectedness` and
  `reverse-reaction-attracting-oscillations`, neither about this question. The mmislan account (Palomar) has p038,
  p064, p079 and p080, all on sequential distributive phosphorylation, not the mixed mechanism.
* **Web search** (mixed-mechanism phosphorylation, subcritical Hopf, Lyapunov coefficient, global stability) found only
  the source and the related distributive-network papers.

Conclusion: no prior resolution was found in this bounded audit.
