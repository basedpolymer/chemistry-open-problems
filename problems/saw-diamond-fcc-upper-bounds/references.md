# References and novelty search

## Sources

| Reference | Identifier | What was read | Used for |
|---|---|---|---|
| S. E. Alm, "Upper bounds for the connective constant of self-avoiding walks", *Combin. Probab. Comput.* 2 (1993) 115–136 | doi:10.1017/S0963548300000547 | full text | previous best rigorous bounds: §6.3.3 and Table 11, face-centred cubic μ < 10.27777 (Table 2: 10.277761; n = 10, m = 4, 313 × 313 matrix); §6.3.4 and Table 12, tetrahedral (= diamond) lattice μ < 2.90940 (n = 22, m = 9, 1003 × 1003 matrix); earlier bounds quoted there: fcc 10.361 (Ahlberg–Janson, counts n ≤ 12), diamond 2.9175 (Wall–White 1976) |
| N. Madras, G. Slade, *The Self-Avoiding Walk*, Birkhäuser (1993; reprint 2013) | doi:10.1007/978-1-4614-6025-1 | §1.2 and the notes of chapter 1 (full text) | μ ≤ (c_N/c_1)^{1/(N−1)} (attributed to Alm via Ahlberg–Janson) and the Ahlberg–Janson polynomial (1.6.1); Table 1.1 covers hypercubic lattices only |
| A. Pönitz, P. Tittmann, "Improved upper bounds for self-avoiding walks in Z^d", *Electron. J. Combin.* 7 (2000) R21 | doi:10.37236/1499 | full text | loop-size state space (hypercubic lattices only) |
| Result 9 of this repository (`problems/cubic-saw-upper-bound/`) | — | full (notes, sources, checker) | method: forbidden-set automata, candidate children through contacts, Collatz–Wielandt certificates |
| O. Couronné, "New upper bound for the connective constant for square-lattice self-avoiding walks", arXiv:2211.16146 | arXiv:2211.16146 | full text (searched) | eigenvector-optimised choice of children (square lattice) |
| R. D. Schram, G. T. Barkema, R. H. Bisseling, N. Clisby, "Exact enumeration of self-avoiding walks on BCC and FCC lattices", *J. Stat. Mech.* (2017) 083208 | doi:10.1088/1742-5468/aa819f; arXiv:1703.09340 | not read; counts taken from OEIS A001336 (b-file "from Schram et al.", n ≤ 24) | exact c_n, fcc, n ≤ 24 (cross-check a_n ≥ c_n; bounds derivable from the counts) |
| OEIS A001394 (number of n-step SAWs on the diamond lattice), b-file n ≤ 32 | https://oeis.org/A001394 | entry and b-file | exact c_n, diamond, n ≤ 32 (terms to n = 27 from Guttmann 1989); cross-check; bounds derivable from the counts |
| OEIS A001336 (number of n-step SAWs on the fcc lattice), b-file n ≤ 24 | https://oeis.org/A001336 | entry and b-file | as above, fcc |
| A. J. Guttmann, "On the critical behaviour of self-avoiding walks II", *J. Phys. A* 22 (1989) 2807–2813 | doi:10.1088/0305-4470/22/14/027 | not read (cited via OEIS and Alm) | diamond counts to n = 27 |
| N. Clisby, "The growth constant for self-avoiding walks on the fcc and bcc lattices", *J. Phys. A* 55 (2022) 465003 | doi:10.1088/1751-8121/aca189 | abstract (OpenAlex) | estimate μ_fcc = 10.03705785(14) |
| Q. He, "Upper bounds for the connective constant of weighted self-avoiding walks", arXiv:2508.01993 (2025) | arXiv:2508.01993 | full text searched | Alm's method restated; numerical results for square, hexagonal, triangular, simple cubic only |
| S. Kim, F. Pinna, arXiv:2509.04568 (2025) | arXiv:2509.04568 | full text searched for "diamond", "fcc", "face" | no statement about these lattices (quotes Z^2 and Z^3 bounds only) |
| G. R. Grimmett, Z. Li, "Self-avoiding walks and connective constants", arXiv:1704.05884; "Bounds on connective constants of regular graphs", arXiv:1210.6277 | arXiv | full text searched | no numerical bound for the diamond or fcc lattice |
| C. P. Herrero, "Self-avoiding walks and connective constants in zeolitic frameworks", *J. Phys.: Condens. Matter* 7 (1995) | doi:10.1088/0953-8984/7/47/010 | abstract | numerical estimates only |

## Novelty search (25 October 2026)

* Literature: Alm 1993 is the last published rigorous upper bound for both lattices that we found. OpenAlex: all
  39 works citing Alm 1993 (titles scanned): the only later upper-bound papers are Noonan 1998 and Pönitz–Tittmann 2000
  (hypercubic lattices), Alm 2005 (two-dimensional Archimedean and Laves lattices), Alm 2004 (hexagonal lattice) and
  He 2025 (weighted walks, square/cubic/hexagonal/triangular). Madras–Slade (1993) and MathWorld list hypercubic
  lattices only. Web searches for "connective constant" with "diamond lattice", "face-centred/centered cubic", "upper
  bound" returned nothing newer.
* Bounds derivable from the published exact counts (computed here, not published as such as far as we know):
  diamond (c_32/c_1)^{1/31} = 2.919627 and Ahlberg–Janson (N = 32) 2.917190, both weaker than Alm's 2.90940;
  fcc (c_24/c_1)^{1/23} = 10.244119 and Ahlberg–Janson (N = 24) 10.240585, both **stronger** than Alm's 10.277761
  (the exact counts to n = 24 date from 2017). The best rigorous fcc bound available before this work is therefore
  10.240585 (or 10.244119 with the elementary inequality alone), derived from Schram et al. 2017.
* AI campaigns and preprint servers:
  * `gh api "users/ipitchford/repos?per_page=100" --paginate`: 193 repositories; names and descriptions grepped for
    saw, self-avoid, connective, walk, polymer, lattice, diamond, fcc, cubic, growth constant: no relevant hit.
  * `gh api "users/mmislan/repos?per_page=100"`: 27 repositories, reaction networks and binary polymer models only.
  * `gh search repos "connective constant"`: `tejasvi/connective-constant-bound` and `tejasvi-releases/*` (October 2026:
    upper bound 2.6586513 for the **square** lattice, Zenodo 10.5281/zenodo.23599811), `vitplanocka/saw-connective-constant`
    (2025–2026, **lower** bounds for Z^2), `qidong-he/weighted_connective_constant` (He 2025), Union-Jack lattice
    projects, Monte Carlo projects. Nothing on the diamond or fcc lattice. `gh search repos` for "diamond lattice
    self-avoiding", "fcc self-avoiding", "self-avoiding upper bound": no results.
  * Zenodo API (`https://zenodo.org/api/records?q=...`): "connective constant" (10 hits: the square-lattice upper-bound
    records above, a hard-square record, a Z^2 lower-bound archive), "connective constants" (1), "self-avoiding" AND
    "diamond" (3, unrelated), "self-avoiding" AND "face-centered" (0), "self-avoiding" AND "fcc" (0),
    "self-avoiding walk" AND "upper bound" (9, same square-lattice records). Raw outputs in `certificates/sources/`
    (git-ignored).
