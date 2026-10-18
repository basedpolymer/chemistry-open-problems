# References

## Primary source (the open problem)

- B. Boros, J. Hofbauer, *Oscillations in planar deficiency-one mass-action systems*, J. Dyn. Diff. Equat.
  36 (Suppl 1):175–197 (2024), doi:10.1007/s10884-021-10051-z; arXiv:2103.00972 (v1, 2021, the only arXiv
  version); open access PMC10901957.
  **Read: full text** (arXiv PDF and LaTeX source `sections/*.tex`). Used: §5 (network (15), ODE (16), scaling
  (20)–(23)), §5.1 (the conjecture, quoted verbatim in README/notes), §3.2 and §4.1 (closed-form / numerical
  focal values used to validate our normalisation: \(L_1\) of §3.2, §4.1, §5.1; \(L_2\) of §4.1 and §3.2; \(L_3\)
  of §4.1). The published version (PMC) has the same wording of §5.1 (checked by the novelty-check agent and
  by `rapports/verify_crnt.md` §C1).
- B. Boros, GitHub repository `balazsboros/reaction_networks`, folder `dfc1thm_oscillation` (notebook supporting
  §5.1 with numerics only; no proof). **Read: inspected by the novelty-check agent (notebook PDF).**

## Novelty check (2021–2026): works checked, none proves four limit cycles for three reactions

- B. Boros, J. Hofbauer, *Limit cycles in mass-conserving deficiency-one mass-action systems*, EJQTDE
  2022(42):1–18, doi:10.14232/ejqtde.2022.1.42, arXiv:2202.10406. **Read: abstract + introduction** (three limit
  cycles, mass-conserving networks; nothing on four cycles with three reactions).
- B. Boros, J. Hofbauer, *Some minimal bimolecular mass-action systems with limit cycles*, Nonlinear Anal. RWA
  72:103839 (2023), doi:10.1016/j.nonrwa.2023.103839, arXiv:2202.11034. **Read: abstract + full-text search.**
- M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*, Stud. Appl.
  Math. 152(1):249–278 (2024), doi:10.1111/sapm.12639, arXiv:2304.02303. **Read: abstract + full-text search**
  (three-reaction quadratic-source class: no isolated periodic orbit with trimolecular targets; supercritical
  Hopf only — does not cover real exponents).
- M. Banaji, B. Boros, J. Hofbauer, *Bifurcations in planar, quadratic mass-action networks with few reactions and
  low molecularity*, Nonlinear Dyn. 112(23):21425–21448 (2024), doi:10.1007/s11071-024-10068-1,
  arXiv:2406.13451. **Read: abstract + statement of the result "no bifurcations of codimension greater than
  m−2" (networks without conservation laws, m reactions)** — explains why the exponents must be varied.
- G. Craciun, R. Erban, *Planar chemical reaction systems with algebraic and non-algebraic limit cycles*,
  arXiv:2406.05057 (J. Math. Biol. 2025). **Read: abstract + full-text search** (four algebraic cycles, many
  reactions, degree 9; not three reactions).
- R. Erban, H.-W. Kang, *Chemical systems with limit cycles*, arXiv:2211.05755. **Read: abstract.**
- N. Yeung, R. Erban, *Classification of strongly non-hyperbolic critical points of planar second-order
  chemical reaction systems*, arXiv:2610.06130 (Oct. 2026). **Read: abstract + full-text search.**
- Citing-work lists: Semantic Scholar API (13 records) and OpenAlex (W3133534872, 12 records); Google Scholar
  only partially reachable. GitHub `mmislan` (Palomar, 23 repositories): nothing on limit cycles / CRNT
  oscillations. Web searches on "four limit cycles" + "mass-action" + "three reactions": only the original
  paper. (Novelty check performed 2026-10-18 by a sub-agent; see `PROGRESS.md`.)

## Methods (background; cited, not re-derived unless stated)

- N. N. Bautin, *On the number of limit cycles which appear with the variation of coefficients from an
  equilibrium position of focus or center type*, Mat. Sb. 30(72) (1952) 181–196; AMS Transl. 100 (1954).
  (Finite cyclicity of a weak focus; used only in a remark.)
- A. A. Andronov, E. A. Leontovich, I. I. Gordon, A. G. Maier, *Theory of Bifurcations of Dynamic Systems on a
  Plane*, Wiley/Israel Program for Scientific Translations, 1973. (Return map near a focus, multiple foci.)
- R. Roussarie, *Bifurcations of Planar Vector Fields and Hilbert's Sixteenth Problem*, Birkhäuser, 1998.
- M. Han, P. Yu, *Normal Forms, Melnikov Functions and Bifurcations of Limit Cycles*, Springer, 2012.
- V. G. Romanovski, D. S. Shafer, *The Center and Cyclicity Problems: A Computational Algebra Approach*,
  Birkhäuser, 2009. (Complex Lyapunov-function recursion for focus quantities.)
- C. Christopher, C. Li, *Limit Cycles of Differential Equations*, Birkhäuser, 2007.
- R. Krawczyk, *Newton-Algorithmen zur Bestimmung von Nullstellen mit Fehlerschranken*, Computing 4 (1969)
  187–201; R. E. Moore, *A test for existence of solutions to nonlinear systems*, SIAM J. Numer. Anal. 14 (1977)
  611–615; A. Neumaier, *Interval Methods for Systems of Equations*, Cambridge Univ. Press, 1990.
- F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput.
  66(8) (2017) 1281–1292; FLINT / python-flint 0.9.0 (exact multivariate polynomials, Arb balls).
- PARI/GP (cross-check of the gcd over the number field in §7 of notes.md); mpmath 1.3.0; SciPy 1.18 (DOP853).

(The background references are cited from general knowledge for context; none of the proofs in `notes.md`
relies on a specific theorem number from them. Lemmas 4.2, 6.1, 6.2 and the successive-perturbation argument
are proved in `notes.md`.)
