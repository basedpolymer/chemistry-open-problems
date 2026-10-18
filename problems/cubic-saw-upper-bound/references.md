# References

| Reference | Identifier | What was read | Used for |
|---|---|---|---|
| A. Pönitz, P. Tittmann, "Improved upper bounds for self-avoiding walks in Z^d", *Electron. J. Combin.* 7 (2000) R21 | doi:10.37236/1499 | full text (PDF from combinatorics.org) | the method; Table 2 (d = 2..5 bounds; d = 3: 4.7599, 4.7476, 4.7387 for k = 10, 12, 14; d = 4: 6.8040 at k = 14; d = 5: 8.8602 at k = 10); reproduced exactly (notes §5) |
| O. Couronné, "New upper bound for the connective constant for square-lattice self-avoiding walks", arXiv:2211.16146 (2022) | arXiv:2211.16146v2 | full text | multiple children + eigenvector-based choice; small bridges; remark that only deletion of old vertices and small bridges carry over to d ≥ 3 |
| N. Clisby, R. Liang, G. Slade, "Self-avoiding walk enumeration via the lace expansion", *J. Phys. A* 40 (2007) 10973 | doi:10.1088/1751-8113/40/36/003; preprint https://personal.math.ubc.ca/~slade/se.pdf | full text (§7.2) | "It is erroneously claimed in [47] that c26 gives µ ≤ 4.7114 for d = 3; in fact c26 gives the weaker bound µ ≤ 4.7626"; their list of best upper bounds 4.7387, 6.8179, 8.8602, 10.8886 (d = 3..6; the d = 4 entry is P–T's k = 10 value, P–T's own table has 6.8040 at k = 14) |
| Clisby–Liang–Slade enumeration data (saw_d3.txt, saw_d4.txt, saw_d5.txt) | https://personal.math.ubc.ca/~slade/lacecounts/ | data files | exact c_n (d = 3: n ≤ 30) for the path-count cross-check; agrees with OEIS A001412 |
| S. Kim, F. Pinna, "Existence and bounds of growth constants for restricted walks, surfaces, and generalisations", arXiv:2509.04568 (2025) | arXiv:2509.04568 | full text searched; quote checked | state of the art as of 2025: "4.57213 ≤ µ_SAW ≤ 4.73871 for the cubic lattice (Hara et al., 1993, Pönitz and Tittmann, 2000)" |
| Q. He, "Upper bounds for the connective constant of weighted self-avoiding walks", arXiv:2508.01993 (2025), *J. Phys. A* (2025) | arXiv:2508.01993; doi:10.1088/1751-8121/ae280e | full text (§2–3) | Alm-type method, small m, n only; no improvement of the unweighted d = 3 bound |
| S. E. Alm, "Upper bounds for the connective constant of self-avoiding walks", *Combin. Probab. Comput.* 2 (1993) 115–136 | doi:10.1017/S0963548300000547 | not read (cited via P–T, He) | alternative method (matrices between walks of lengths m and n) |
| J. Noonan, "New upper bounds for the connective constants of self-avoiding walks", *J. Stat. Phys.* 91 (1998) 871–888 | doi:10.1023/A:1023023831510 | not read (cited via P–T and CLS) | earlier bounds 2.6939 (d = 2), 4.7476 (d = 3) |
| N. Clisby, "Calculation of the connective constant for self-avoiding walks via the pivot algorithm", *J. Phys. A* 46 (2013) 245001 | arXiv:1302.2106 | abstract | numerical estimate µ_3 ≈ 4.684039931(27) |
| T. Hara, G. Slade, A. D. Sokal, "New lower bounds on the self-avoiding-walk connective constant", *J. Stat. Phys.* 72 (1993) 479 | doi:10.1007/BF01048021 | not read | lower bound 4.57213 quoted for context |

## Novelty search (2026-10-18)

* OpenAlex: all 43 works citing Pönitz–Tittmann 2000 (W1721599859) and all 39 citing Alm 1993 (W2005105063) —
  titles scanned; none improves the d = 3 upper bound (the only connective-constant upper-bound papers are Alm 2005
  (Archimedean/Laves lattices), Couronné 2022 (d = 2), He 2025 (weighted, small memory)).
* Semantic Scholar: citations of arXiv:2211.16146 — one (Kim–Pinna 2025), which quotes 4.73871 as current.
* arXiv full-text API query `abs:"connective constant"`, 60 most recent entries up to 2026-09-23: nothing on
  rigorous upper bounds for Z^3.
* Web searches ("upper bound connective constant self-avoiding walk simple cubic lattice 4.7387 improved";
  "connective constant upper bound cubic lattice 2024 2025 2026 automaton"): only the erroneous 4.7114
  (MacDonald et al. 2000, corrected by CLS 2007) and the sources above.

## Additional novelty search requested by the orchestrator (2026-10-18, ~12:45–13:00 local time)

Goal: rule out an unpublished/AI-campaign improvement of the Z^3 upper bound (GitHub, Zenodo, arXiv 2024–2026).
Raw outputs: `~/work/cubic-saw/logs/zenodo_search.log`, `~/work/cubic-saw/logs/arxiv_search.log` (WSL).

1. **GitHub, Evidence Press** — `gh api "users/ipitchford/repos?per_page=100" --paginate`: 166 repositories
   (created 2026); names/descriptions grepped for saw, self-avoid, connective, walk, polymer, lattice, cubic,
   automat: only two hits on "cubic" (`exotic-affine-spheres-quadratic-cubic`, `four-point-inductive-cubical-excision`),
   unrelated. **No SAW / connective-constant repository.**
2. **GitHub search** — `gh search repos "connective constant"` (6 results):
   * `vitplanocka/saw-connective-constant` (created 2026-08-22): two rigorous **lower** bounds for µ(Z^2)
     (2^27/51082230 = 2.627483725…, then 2.627483880…), exact-integer certificates, Lean 4. Its
     `docs/LITERATURE.md` (read) reports two errata: (2a) the abstract of Pönitz–Tittmann prints
     2.679192495 for the d = 2, k = 22 constant whereas the automaton's constant is certified in
     [2.679192641375, 2.679192650414] — only the 9-digit string is wrong, Table 2's four-decimal entries are
     correct; (2b) Couronné's planar pruning rules (§4.1–4.2) are applied as move filters, which is unsound as
     written, repairable without changing the value. **No Z^3 work.** Relevance here: (2a) does not affect the
     d = 3 record (Table 2 value 4.7387; our P–T automaton gives 4.738630775); we independently confirm (2a):
     our P–T-rule automaton for d = 2, K = 22 (8 880 776 states) gives λ ∈ [2.6791926494, 2.6791926495]
     (floating-point Collatz–Wielandt bracket), inside their certified interval and above the printed digits.
     (2b) does not concern our method, which only forbids steps onto remembered visited sites.
   * `qidong-he/weighted_connective_constant` (2025, code of arXiv:2508.01993, Alm method, small memory);
     `alex21347/Upper-bound-to-the-connective-constant-of-Union-Jack-Lattice` (2020, Alm method, Union-Jack
     lattice); three Monte Carlo estimation repositories. None on Z^3 upper bounds.
   * `gh search repos "self-avoiding walk upper bound"`, `"self-avoiding walk cubic"`: no results;
     `gh search repos self-avoiding --sort updated` (40 most recent): simulations, enumerations (taxi walks,
     A252653, NE-prudent walks), nothing on rigorous upper bounds.
   * `gh search code`: `"4.7387" self-avoiding` and `Ponitz Tittmann`, `"connective constant" cubic`: no results;
     `"4.73871"`: only unrelated numeric data files.
3. **Zenodo API** (`https://zenodo.org/api/records?q=...&size=25`): exact phrases `"connective constant"` (2 hits),
   `"connective constants"` (1), `"self-avoiding" AND "upper bound"` (1), `"self-avoiding walk" AND "connective"` (2),
   `"self-avoiding walks" AND "Z^3"` (0), `"finite memory" AND "self-avoiding"` (0), plus author names and free-text
   queries. The only SAW connective-constant record is 10.5281/zenodo.21546041 (Numaro, 2026-07-25):
   certificate archive for the **lower** bound µ(Z^2) ≥ 2.6273856 (Kesten bridges, span 17), submitted to
   `teorth/optimizationproblems` (constant 38a). The other hit (10.5281/zenodo.23164583) is about the hard-square
   lattice gas.
4. **teorth/optimizationproblems** (Tao's list of constants, pushed 2026-10-18): the only SAW entry is 38a
   (square lattice; upper bound listed 2.679193); **no entry for Z^3**.
5. **arXiv API**, submittedDate 2024-01-01 … 2026-12-31: `abs:"connective constant" AND abs:"upper bound"` (2:
   Kim–Pinna 2509.04568, He 2508.01993); `all:` variant (same 2); `abs:"self-avoiding" AND abs:"upper bound" AND
   abs:cubic` (1: Kim–Pinna); `abs:"self-avoiding walk" AND abs:"cubic lattice"` (6: trails Monte Carlo
   2609.21724 — abstract read, numerical only; copolymers; quantum-inspired enumeration; Kim–Pinna; obstacle
   lattices; star polymers); `abs:"self-avoiding walks" AND abs:"connective constant"` (8: Grigorchuk graphs,
   nested fractals, hard-core zero-freeness, cubic *graphs* 2601.12571 — abstract read, exact µ for
   gadget-substituted 3-regular graphs, not Z^3 —, Kim–Pinna, He, particle systems, Cayley graphs);
   `ti:"connective constant"` (6, none on Z^d upper bounds); `abs:"finite memory" AND abs:"self-avoiding"` (1,
   unrelated); `abs:"self-avoiding" AND abs:"Collatz"` (0).

Conclusion of the extended search: no improvement of µ(Z^3) ≤ 4.7387 (nor of the d = 4, 5, 6 bounds) was found in
journals, arXiv, Zenodo or GitHub as of 2026-10-18; the active AI campaigns found (vitplanocka, Numaro) work on
**lower** bounds for Z^2.

| Additional reference | Identifier | What was read |
|---|---|---|
| vitplanocka, *saw-connective-constant* (GitHub, 2026) | github.com/vitplanocka/saw-connective-constant | README, docs/LITERATURE.md §1–2 |
| Numaro, certificate archive µ(Z^2) ≥ 2.6273856 (Zenodo, 2026) | doi:10.5281/zenodo.21546041 | record metadata/description |
| T. Tao et al., *optimizationproblems*, constant 38a | github.com/teorth/optimizationproblems | constants/38a.md, README list |
| X. Su, Z. Zhou, Q. Wu, "Self-avoiding trails in two and three dimensions" (2026) | arXiv:2609.21724 | abstract |
| B. Grant, Z. Li, "Self-avoiding walks on cubic graphs and local transformations" (2026) | arXiv:2601.12571 | abstract |
| G. Slade, "Bounds on the self-avoiding-walk connective constant", J. Fourier Anal. Appl. (Kahane special issue, 1995), reprinted CRC 2020 | doi:10.1201/9780429332838-32; https://personal.math.ubc.ca/~slade/kahane.pdf | full text (6 p.). Table 1: best rigorous upper bounds at the time are Alm's: 4.756, 6.832, 8.881, 10.903 for d = 3..6 (lower bounds from Hara–Slade–Sokal and the lace expansion). Confirms that d ≥ 5 upper bounds come from Alm-type methods, later improved by Noonan 1998 and Pönitz–Tittmann 2000 |
