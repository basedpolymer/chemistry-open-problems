# References and novelty check

## Sources read

| Source | Identifier | What was read |
|---|---|---|
| B. Y. Bisbis, *On the steady states and dynamic behavior of two CSTR's in series*, PhD thesis, Univ. Missouri–Columbia, Dec. 1987 | doi:10.32469/10355/98283 (MOspace, open access) | Full text (70 pages, OCR plus page images of pp. 3–5, 27, 29, 44). Model eqs. (1)–(8), dimensionless groups p. 3, steady-state equations (9)–(14), Chapter II (sequential bifurcations, nine states for n = 2, numerical example with ε = 0, B = 13), concluding remarks p. 27 (the conjecture), Chapter III p. 29 (which parameters are shared), Jacobian p. 31, references p. 45 |
| C. Galanics, I. Szalai, *Dynamics of pH oscillators in continuous stirred tanks in series*, ChemPhysChem 25 (2024) e202400610 | doi:10.1002/cphc.202400610, PMC11614366 | Full text (Europe PMC XML): quotes Varma's bound 2^(N+1) − 1 with N + 1 stable |
| S. Svoronos, R. Aris, G. Stephanopoulos, *On the behaviour of two stirred tanks in series*, Chem. Eng. Sci. 37 (1982) 357–366 | doi:10.1016/0009-2509(82)80087-0 | Abstract only (ASU Pure page): two identical tanks, up to seven steady states for a given Damköhler number, at most three stable |
| R. Li, *Steady-state multiplicity of isothermal reacting systems in a sequence of two CSTRs*, Chem. Eng. Sci. 49 (1994) 3053–3057 | doi:10.1016/0009-2509(94)E0091-4 | Abstract only, through a Scribd description and a web-search summary of the ADS abstract: nine states (four stable) for two isothermal CSTRs, "extended" to 3^N states, 2^N stable, for N CSTRs. Isothermal model (no energy balance), so a different problem. Paywalled; see NEEDS.md |
| A. Varma, *On the number and stability of steady states of a sequence of continuous-flow stirred tank reactors*, Ind. Eng. Chem. Fundam. 19 (1980) 316–319 | doi:10.1021/i160075a016 | Not read (paywalled). Content known from Bisbis (p. 7) and Galanics–Szalai: bound 2^(n+1) − 1, at most n + 1 stable, for adiabatic and isothermal operation |
| M. Kubíček, H. Hofmann, V. Hlaváček, J. Sinkule, *Multiplicity and stability in a sequence of two nonadiabatic nonisothermal CSTR*, Chem. Eng. Sci. 35 (1980) 987–996 | doi:10.1016/0009-2509(80)85088-3 | Not read (paywalled, no abstract in OpenAlex/Semantic Scholar). Bisbis (p. 7) says its plots show at most seven states. One web-search summary attributed a 3^N claim to it, but this turned out to be a conflation with Li (1994) |
| D. G. Retzloff, P. C.-H. Chan, C. C. Chicone, I. J. Papick, *Maximal multiplicity for sequential bifurcations of a first-order reaction occurring in CSTRs coupled in series*, SIAM J. Appl. Math. 52 (1992) 1136–1147 | doi:10.1137/0152067 | Abstract only (OpenAlex): two reactors, nine critical points, winged-cusp singularity. No n-reactor statement in the abstract |
| G. Dangelmayr, I. Stewart, SIAM J. Math. Anal. 15 (1984) 423–445; SIAM J. Appl. Math. 45 (1985) 895–918 | doi:10.1137/0515033, doi:10.1137/0145054 | Not read; cited for context, as in Bisbis |
| A. Uppal, W. H. Ray, A. B. Poore, Chem. Eng. Sci. 31 (1976) 205–214 | doi:10.1016/0009-2509(76)85058-0 | Not read; source of the nondimensionalisation, as stated in Bisbis |
| L. Perko, *Differential Equations and Dynamical Systems*, 3rd ed., Springer 2001 | doi:10.1007/978-1-4613-0003-8 | Standard linearisation stability facts (§2.9), cited from knowledge |
| J. Bremer, K. Sundmacher, *Novel Multiplicity and Stability Criteria for Non-Isothermal Fixed-Bed Reactors*, Front. Energy Res. 8 (2021) 549298 | doi:10.3389/fenrg.2020.549298 | Relevant passages (open access): notes the 3^n combinatorial count for a cascade model of a fixed bed; with its constraints only seven of nine combinations occur for two stages. Different constrained setting; no attainment claim |

## Novelty check (5 November 2026)

- OpenAlex citing works of Varma (1980), 18 items, titles checked: none treats attainment of 3^n for jacketed
  cascades. Citing works of Retzloff et al. (1992): one (AIChE J. 1997, catalytic CSTR), unrelated.
- Web searches (standard and extended) on "3^n steady states CSTRs in series", "cascade of n nonadiabatic CSTRs
  maximum number of steady states", "sequence of N CSTRs steady states of which stable": only Li (1994, isothermal),
  the Frontiers 2021 fixed-bed paper, and the Bisbis thesis.
- Zenodo API (`CSTRs in series steady states 3^n`, `cascade CSTR multiplicity 3^n`, `Bisbis CSTR`): nothing relevant.
- `gh search repos` ("CSTR series steady states", "CSTR cascade multiplicity"): nothing. `gh search code`
  ("3^n steady states", "CSTRs in series" steady states): only teaching/simulation code and this repository's own
  candidate list.
- `gh api users/ipitchford/repos` and `users/mmislan/repos`: no CSTR or cascade repository.

Conclusion: no proof that 3^n steady states (with 2^n stable) are attained for every n, or even for n = 3, in the
jacketed nonisothermal model was found. The closest statement is Li's (1994) abstract for isothermal sequences, a
different model.
