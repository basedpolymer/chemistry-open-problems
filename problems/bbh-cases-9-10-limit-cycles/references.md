# References and novelty check

## Sources used

| Source | Read | Use |
|---|---|---|
| M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*, Stud. Appl. Math. 152 (2024) 249–278, doi:10.1111/sapm.12639, arXiv:2304.02303 | full text, arXiv v1 (pdftotext) | Open statement (§6, p. 22), notation (13)–(18), Lemma 2 (unique equilibrium), Lemma 4 (Dulac 1/(xy)), §4.2 (generalised LVA), Theorem 9 (cases 9, 10) and its proof |
| Same paper, published version, open access, PMC10953466 (Europe PMC full-text XML) | Theorem 3, §6, theorem/lemma labels | Numbering map: arXiv Lemma 2/3/4 = published Lemma 1/2/3; Theorem 9 = Theorem 3; Theorem 10 = Theorem 4. Theorem 3 and §6 text identical to arXiv |
| This repository, result 2: `problems/bbh-three-reaction-no-limit-cycle/` (README, paper.tex, verification.md, scripts) | full | The example and the subfamily d₃ = 0, d₁ ≥ 2c₁ of 9a were settled there with the Dulac function y^(−b). Not modified |
| B. Boros, GitHub `balazsboros/reaction_networks` (cloned 24 Oct 2026) | file list, commit log | Folder `3reactions/` (focal-value notebook) unchanged; last commit 30 Sep 2026 ("journal info"). Nothing on non-Hopf networks |

**Standard facts used (F4 of `result.md`), not re-read for this work.**
- Poincaré index theory (sum of indices inside a periodic orbit is 1, a saddle has index −1).
- The Bendixson–Dulac criterion.
- The Jordan curve theorem.

They are found, for example, in L. Perko, *Differential Equations and Dynamical Systems*, 3rd ed., Springer (2001),
Ch. 3 (Bendixson's criteria; index theory). Section numbers are quoted from memory, not re-checked. BBH cite
M. Y. Li, Rocky Mountain J. Math. 25 (1995) 351–363 for Bendixson–Dulac.

The comparison lemma (Lemma 3.4 of `result.md`) is proved in full there. It is in the spirit of Filippov-transformation
arguments for Liénard systems. We did not locate a reference with exactly this statement, and none is relied on.

## Novelty check (24 Oct 2026)

**Citing works.**
- Semantic Scholar lists 12 citing works, 2023–2026: arXiv 2610.06130 (Yeung–Erban), 2608.18509, 2608.01116,
  2507.09272, 2503.00008, 2409.06877 (*Positive equilibria in mass action networks: geometry and bounds*),
  2406.13451 (BBH, Nonlinear Dyn. 2024), 2406.09514, 2406.05057 (Craciun–Erban), 2402.18523, 2312.12897 and
  2307.04186.
- OpenAlex (W4386927497) lists 10. Its API budget was exhausted on the shared IP on 24 Oct, so the full OpenAlex list
  was not re-fetched; the 18 Oct check of result 2 used it.
- Full text grepped on 24 Oct for "limit cycle", "case 9", Dulac, Liénard:
  - 2406.13451: mentions cases 8–10 only to place 104 Hopf networks;
  - 2409.06877 and 2610.06130: BBH cited only in bibliography or introduction.
- The other arXiv citing works were grepped by the result-2 verifier on 18 Oct, with the same outcome.

**Later work of the authors.**
- arXiv author listings (Banaji, Boros, Hofbauer), up to 24 Oct 2026: newest relevant items are Banaji 2609.07265
  (fold/cusp via dual equations) and Boros et al. 2510.03621 (disguised toric locus). Neither addresses the question.
- Boros's GitHub repository: see above.

**AI campaigns.**
- `mmislan` (Palomar): 27 repositories, read via GitHub search on 24 Oct. CRNT, autocatalytic sets, D-stability, and
  one on classical Liénard lower bounds (`palomar-lienard-limit-cycles`). None on BBH or three-reaction networks.
- `ipitchford` (Evidence Press): 191 repositories, names and descriptions filtered. The CRN-related ones are
  `reverse-reaction-attracting-oscillations`, `quadratic-mass-action-sharpness`, `two-host-dormancy-oscillations` and
  `wilhelm-heinrich-large-gain`. None is on BBH §6.

**Zenodo.** API queries with required terms gave 0 relevant hits: "three-reaction" + "mass-action", Banaji Boros
Hofbauer, "limit cycle" + "mass-action" + quadratic, "Andronov-Hopf" + "mass-action", Dulac + "reaction network",
Lienard + "mass-action", the arXiv id and the DOI. Free-text queries returned only unrelated records.

**Web search** on the key phrases and on the example network: only the BBH paper itself and unrelated work.

**Conclusion: no prior answer to the general question was found.** The example and one subfamily were answered in
this repository on 18 Oct 2026 (result 2). The general answer here appears to be new.
