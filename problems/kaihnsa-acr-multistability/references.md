# References and novelty check

| Source | Read | Use |
|---|---|---|
| N. Kaihnsa, T. Nguyen, A. Shiu, *Absolute concentration robustness and multistationarity in reaction networks: conditions for coexistence*, Eur. J. Appl. Math. 35 (2024) 566–600, doi:10.1017/S0956792523000335, arXiv:2307.04186 | full text (arXiv v1) | Open statement (§6, Discussion); Definition 2.6 (ACR); Lemmas 2.3–2.5; Propositions 3.1, 3.5, 3.10; Theorems 4.1, 4.8, 5.1, 5.5; Example 4.10 |
| G. Shinar, M. Feinberg, *Structural sources of robustness in biochemical reaction networks*, Science 327 (2010) 1389–1391, doi:10.1126/science.1183372 | full text | origin of ACR |
| M. Feinberg, *The existence and uniqueness of steady states for a class of chemical reaction networks*, Arch. Rational Mech. Anal. 132 (1995) 311–370, doi:10.1007/BF00375614 | statement of the deficiency-one theorem (checked by the complexes sub-study) | uniqueness of positive steady states for one linkage class of deficiency one with one terminal strong linkage class |
| M. Feinberg, *Foundations of Chemical Reaction Network Theory*, Springer 2019, doi:10.1007/978-3-030-03858-8 | standard reference | Laplacian factorization, terminal strong linkage classes |
| M. Feinberg, F. J. M. Horn, *Chemical mechanism structures and the coincidence of the stoichiometric and kinetic subspaces*, Arch. Rational Mech. Anal. 66 (1977) 83–97, doi:10.1007/BF00250853 | standard reference | kernel of the Laplacian, kinetic subspace |
| B. Joshi, T. D. Nguyen, M. D. Johnston, *Bistability, absolute concentration robustness, and hysteresis in dual-site futile cycles with bifunctional enzymes*, arXiv:2605.11389 (May 2026) | full text (grepped, §4.5 and §5 read) | earlier example of ACR with bistability: 6 species, dimension 4 (network (E,E,C2,C2)); does not address minimal sizes |
| Y. Jiao, X. Tang, X. Zeng, *Multistability of small zero-one reaction networks*, J. Math. Biol. (2025), doi:10.1007/s00285-025-02306-w, arXiv:2406.11586 | arXiv full text (grepped, §2.2 read) | standard definition of (exponentially) stable steady states and multistability; cites KNS but does not treat ACR with multistability |
| X. Si, X. Tang, *Obstruction of absolute concentration robustness by conservation laws in non-redundant zero-one networks*, arXiv:2604.14592 (v5, Sept. 2026) | abstract and introduction | ACR classification for zero-one networks of dimension ≤ 2; no multistability |
| N. Meshkat, A. Shiu, A. Torres, *Absolute concentration robustness in networks with low-dimensional stoichiometric subspace*, Vietnam J. Math. 50 (2022) 623–651, doi:10.1007/s10013-021-00524-5 | abstract and introduction | ACR in one-dimensional and two-reaction networks; not used in the proofs |
| B. Joshi, A. Shiu (SIADS 2017), B. Joshi, G. Craciun (JMB 2022), C. Conradi et al. (PLoS Comput. Biol. 2017), L. D. García Puente et al. (JSC 2025) | local copies obtained, not needed | background on small networks, ACR algebra, determinant criteria |

## Novelty check (27 Oct 2026)

- **Citing works.** OpenAlex lists 4 works citing KNS (W4390503227): García Puente et al. (ACR, algebra and geometry,
  JSC 2024), Jiao–Tang–Zeng (multistability of zero-one networks, JMB 2025), an algorithm for the equivalence of
  zero-one networks (2025), and a paper on the Wnt pathway (2025). None addresses ACR with multistability.
- **Recent ACR literature.** OpenAlex title/abstract search "absolute concentration robustness" since June 2024
  (15 works). The only one combining ACR and bistability is Joshi–Nguyen–Johnston, arXiv:2605.11389, which gives a
  6-species, 4-dimensional example and does not discuss minimal sizes. arXiv abstract searches for "concentration
  robustness" with "bistab*" or "multistab*" return only that paper.
- **Later work of the authors.** Nguyen is a coauthor of arXiv:2605.11389 (above); nothing else on this question was
  found.
- **AI campaigns.** `ipitchford` / Evidence Press repositories: no title on ACR or multistability (only unrelated
  "robustness" and "concentration" titles). `mmislan`: `palomar-p082-acr-zero-divisors` (zero-divisor criteria for
  ACR, not multistability). `gh search repos "absolute concentration robustness"`: only ACR tools (ACRE, ACORNS,
  DynamicACR, asymptoticACR). Zenodo search for the key phrases: no relevant record.
- **Conclusion.** The conjecture appears to be open as of 27 October 2026.
