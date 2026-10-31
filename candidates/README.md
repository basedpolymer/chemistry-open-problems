# Candidate problems for the next session

Checked on 31 October 2026. No solver has worked on these, except where stated. They come from a scouting pass (Claude
Sonnet 5.5, effort medium) and a second novelty sweep aimed at unindexed resolutions (same model). The ordering and the
remarks below are the orchestrator's.

**Before starting any of them:**

1. Redo the novelty check of [PROTOCOL.md](../PROTOCOL.md) §1, including `gh search code`. On 31 October two targets that
   looked open turned out to be solved in places a normal literature search misses: a Zenodo working paper (result 30)
   and a GitHub repository plus a personal site (result 31). AI-assisted campaigns react to open statements within days.
2. Read the source PDF and copy the open sentence verbatim into the problem's README. The statements below are
   paraphrases. The scouting agents read the sources, mostly through `pdftotext`; their readings were not re-checked.
3. Fix the exact definitions before any computation. Several statements below need this (see each entry).

Limits of the sweep: GitHub code search is noisy, so "no hit" is weak evidence; arXiv was searched through its metadata
API only; OpenAlex and Semantic Scholar were rate limited; HAL, OSF, Figshare, SSRN and similar hosts were covered only by
restricted web searches.

## Summary

Feasibility means the chance of a *complete* resolution by a very strong AI with a computer algebra system.

| Tier | Problem | Field | Feasibility | Expected form | PDF needed |
|---|---|---|---|---|---|
| A | [Rigorous Hopf bifurcation in the Becker–Döring bubblelator](#a1) | Nucleation kinetics | medium-high | theorem + certified computation | no |
| A | [Quadratic trimolecular (2, 4, 2) networks](#a2) | Reaction networks | medium | per-network theorems, or two limit cycles | no |
| A | [Mixed-mechanism dual phosphorylation](#a3) | Reaction networks | medium | counterexample | no |
| A | [Complexity of the multistage minimum utility cost problem](#a4) | Chemical engineering | medium | NP-hardness or an algorithm | no |
| A | [Maximal energy of bicyclic graphs, two odd cycles](#a5) | Mathematical chemistry | medium, long | proof + interval checks | yes |
| B | [Unimodality of TAC(E)](#b1) | Chemical engineering | medium-high | counterexample or theorem | no |
| B | [n CSTRs in series, 3ⁿ steady states](#b2) | Chemical engineering | medium-high (construction) | construction; upper bound harder | yes |
| B | [Fifteen solutions for three consecutive reactions](#b3) | Chemical engineering | medium | certified parameter set | yes |
| B | [Wiener index of orientations](#b4) | Mathematical chemistry | medium (search) | counterexample | no |
| B | [Clar–Fries question for benzenoids](#b5) | Mathematical chemistry | medium (search) | counterexample | yes |
| B | [Minimal order of explosive endotactic networks](#b6) | Stochastic reaction networks | low-medium | exact value | no |
| B | [Gelation threshold at homogeneity 1](#b7) | Aggregation kinetics | medium-low | new estimate | no |
| B | [Nondegeneracy conjecture](#b8) | Reaction networks | low-medium | counterexample | no |
| C | [Brownian similarity profile uniqueness](#c1) | Aggregation kinetics | very low | proof | no |
| C | [Random-jump heptagons](#c2) | Polymer physics | low | exact value | no |
| C | [Fullerene duals, Conjecture 1(c) and 1(d)](#c3) | Mathematical chemistry | low | counterexample | no |
| C | [Randić index and radius](#d1) (one run, unresolved) | Mathematical chemistry | low | proof | no |

## PDFs to supply

Put them in `candidates/sources/` (git-ignored), named by first author and year. These are the only sources no agent could
read, and the only ones needed before work starts.

| Priority | DOI | Reference | Used for |
|---|---|---|---|
| 1 | [10.1021/ci000164z](https://doi.org/10.1021/ci000164z) | Gutman, Vidović, J. Chem. Inf. Comput. Sci. 41 (2001) | A5, exact statement of the conjecture |
| 2 | [10.1016/0009-2509(94)e0091-4](https://doi.org/10.1016/0009-2509(94)e0091-4) | Li, Chem. Eng. Sci. 49 (1994) | B2, priority check for the isothermal case |
| 2 | [10.1137/0152067](https://doi.org/10.1137/0152067) | Retzloff, Chan, Chicone, Papick, SIAM J. Appl. Math. 52 (1992) | B2, the case n = 2 (abstract read only) |
| 2 | [10.1016/0362-546X(82)90015-3](https://doi.org/10.1016/0362-546X(82)90015-3) | Chicone, Retzloff, Nonlinear Anal. 6 (1982) | B3, the bound 15 |
| 2 | [10.1016/0022-247X(87)90003-5](https://doi.org/10.1016/0022-247X(87)90003-5) | Retzloff, Chan, Mohamed, J. Math. Anal. Appl. 124 (1987) | B3, the bound attained for two reactions |
| 2 | [10.1007/s10910-013-0193-0](https://doi.org/10.1007/s10910-013-0193-0) | Graver, Hartung, Souid, J. Math. Chem. 51 (2013) | B5, original question |
| 3 | [10.1016/0012-365X(88)90199-9](https://doi.org/10.1016/0012-365X(88)90199-9) | Fajtlowicz, Discrete Math. 72 (1988) | D1, history of the statement only |
| 3 | [10.1016/S0012-365X(99)00206-X](https://doi.org/10.1016/S0012-365X(99)00206-X) | Caporossi, Hansen, Discrete Math. 212 (2000) | D1, history of the statement only |

For the verification of result 14, also useful: Grisvard, *Elliptic Problems in Nonsmooth Domains*,
[10.1137/1.9781611972030](https://doi.org/10.1137/1.9781611972030), Theorem 3.2.1.2 (its exact statement is cited, not read).
Details in [the problem's NEEDS.md](../problems/taylor-dispersion-polygon-skewness/NEEDS.md).

Published versions exist for several candidates below, but their arXiv or author copies are free and were read, so they are
not needed.

## Tier A: hard, a complete resolution is plausible

<a id="a1"></a>
### Rigorous Hopf bifurcation in the Becker–Döring bubblelator

- **Statement.** Prove a Hopf bifurcation theorem for the Becker–Döring model with injection of monomers and depletion of
  large clusters, in which the source verifies transversality and nonresonance only numerically (its appendix C).
- **Source.** Niethammer, Pego, Schlichting, Velázquez, *Oscillations in a Becker–Döring model with injection and depletion*,
  SIAM J. Appl. Math., [arXiv:2102.06751](https://arxiv.org/abs/2102.06751) §§1.4, 3, 4.2. The source says that such a theorem
  has not yet been proven for a model of this type, which is weaker wording than "open problem".
- **Status.** Open. No citing work claims rigour. Read arXiv:2508.07884 (injection and fragmentation, no Hopf).
- **Attack.** Reduce to a Volterra-type equation with infinite delay; certify the eigenvalue crossing and the sign of the
  first Lyapunov coefficient in interval arithmetic. The experience of results 8 and 25 helps, but their model is finite in
  a way this one is not, and the infinite-delay functional setting is the main risk.

<a id="a2"></a>
### Quadratic trimolecular (2, 4, 2) networks

- **Statement.** Two questions on planar mass-action networks with four reactions: do the 91 networks without invariant axis
  and without Bautin bifurcation admit at most one limit cycle? Can a homoclinic orbit occur in the seven networks that have
  both fold and Andronov–Hopf bifurcations, with bifurcation sets that do not meet and no Bogdanov–Takens point?
- **Source.** Banaji, Boros, Hofbauer, Nonlinear Dyn. (2024), doi:10.1007/s11071-024-10068-1,
  [arXiv:2406.13451](https://arxiv.org/abs/2406.13451) §7, p. 22. Result 16 settled the bimolecular sibling in the same paper.
- **Status.** Open; no citer addresses them. Read arXiv:2610.27188 (canard points), 2610.06130, 2609.07265, 2610.32739. The
  authors' notebooks are at github.com/balazsboros/reaction_networks.
- **Attack.** Focal values per network, Dulac and Cherkas criteria, validated Poincaré maps; a pair of nested cycles is
  certified as in result 5. The homoclinic question is the smaller of the two and can be done first.

<a id="a3"></a>
### Mixed-mechanism dual phosphorylation

- **Statement.** For the network with processive phosphorylation and distributive dephosphorylation (9 species, unique positive
  steady state per class): if the steady state is locally stable, is it globally stable? A conjecture in the same section
  says that if it is unstable, the class contains a stable periodic orbit.
- **Source.** Conradi, Mincheva, Shiu, Bull. Math. Biol. 81 (2019), doi:10.1007/s11538-019-00580-6,
  [arXiv:1809.02886](https://arxiv.org/abs/1809.02886) §6, pp. 18 and 20. Check the conjecture's number: sources disagree
  between 6.1 and 6.2.
- **Status.** Open. Zenodo 10.5281/zenodo.23446981 says itself that it does not resolve this. Read Vassena arXiv:2608.27081 and
  Feliu, Kaihnsa arXiv:2405.16179.
- **Attack.** The monomial parametrization of the source's Theorem 3.2 gives the first Lyapunov coefficient along the Hopf
  surface. A subcritical point gives a counterexample, certified with an interval Poincaré map. The yes direction looks
  far harder.

<a id="a4"></a>
### Complexity of the multistage minimum utility cost problem

- **Statement.** Is the multistage minimum utility cost problem for heat-exchanger networks polynomial or NP-hard?
- **Source.** Letsios et al., [arXiv:1909.12328](https://arxiv.org/abs/1909.12328) §6.4, pp. 21–22; Comput. Chem. Eng. 2020,
  doi:10.1016/j.compchemeng.2019.106599 (journal version matched by title, not opened).
- **Status.** Open; none of its nine citers settles it.
- **Attack.** Reduction from pooling, 3-Partition or bin packing through the bilinear mixing terms, or a convex
  reformulation. Fix the model first: "multistage" must be read exactly as in the paper.

<a id="a5"></a>
### Maximal energy of bicyclic graphs, two odd cycles

- **Statement.** The Gutman–Vidović conjecture says that the bicyclic graph of maximal energy is P_n^{6,6} for n = 14 and
  n ≥ 16. What remains is the case of two odd cycles, disjoint or sharing an edge.
- **Source.** Gao, Li, Yang, Zheng, Discrete Appl. Math. 2024, doi:10.1016/j.dam.2024.09.014 (free preprint at Nankai);
  original conjecture in the PDF above.
- **Status.** Open after 25 years; the authors expect this case to be much harder than the bipartite one.
- **Attack.** The Coulson-integral quasi-order fails with odd cycles, so compare integrals directly: family inequalities, then
  interval arithmetic for the finite remainder. Long, paper-length work.

## Tier B: promising, with a caveat

<a id="b1"></a>
### Unimodality of TAC(E)

- **Statement.** For a fixed minimal stage-wise heat-exchanger network, the total annualized cost as a function of utility
  energy is claimed to be unimodal between its minimum and maximum energies.
- **Source.** Chang, Peccini, Wang, Costa, Bagajewicz, AIChE J. 66(7) (2020), doi:10.1002/aic.16267, §2; Part II,
  doi:10.1002/aic.16264, §2, where both statements are left as conjectures. Author copies at ou.edu/class/che-design/pub-papers/.
- **Status.** Open: Oliva et al. (2024) still condition on it.
- **Caveat.** Informal statement. A precise model and exponent b are needed, and the answer may depend on b < 1.
- **Attack.** One-parameter calculus per exchanger; symbolic or interval search for two local minima on chains of 2–3
  exchangers with b < 1 (a counterexample is the likely outcome); otherwise quasi-convexity for b = 1.

<a id="b2"></a>
### n CSTRs in series, 3ⁿ steady states

- **Statement.** For n exothermic first-order CSTRs in series, each with its own jacket, the maximum number of steady states is
  3ⁿ, of which at most 2ⁿ are stable. Proved for n = 1, 2.
- **Source.** Bisbis, PhD thesis, Univ. Missouri 1987, doi:10.32469/10355/98283, Ch. II, p. 27 (open access).
- **Caveat.** A remark in a thesis, not a published conjecture. Li (1994) claims the analogue for isothermal sequences; it is a
  different model and has not been read. A review (PMC11614366) quotes a related count with garbled wording.
- **Attack.** Compose the multivalued maps of one stage; compute the images of the fold curves exactly; construct 27 states for
  n = 3 or find an obstruction. The construction is more accessible than the upper bound.

<a id="b3"></a>
### Fifteen solutions for three consecutive reactions

- **Statement.** Is the Chicone–Retzloff bound of 15 steady states attained for three consecutive first-order reactions in one
  CSTR? The analogous bound of 7 is attained for two reactions.
- **Source.** Berdouzi, PhD thesis, Univ. Missouri 1987, doi:10.32469/10355/98342, Conclusion, p. 54 (open access).
- **Caveat.** A single line in a thesis; the scalar equation must be derived again.
- **Attack.** Rolle reduction to a polynomial in x and ln x, find a degenerate point of high multiplicity and perturb it; certify
  15 roots with Sturm sequences or Arb.

<a id="b4"></a>
### Wiener index of orientations

- **Statement.** Conjecture 4.3: among orientations of a graph, a minimizer of the Wiener index can be chosen acyclic. Conjecture
  4.5: for chordal graphs the maximizer is strongly connected. Read the exact hypotheses in the source.
- **Source.** Knor, Škrekovski, Tepeh, Ars Math. Contemp. 24 (2024) #P4.07, doi:10.26493/1855-3974.3077.63a, §4, p. 7.
- **Status.** Open; 18 citers, none on these two. Read Meng, Zhang, Zhang, AMC 2026, doi:10.1016/j.amc.2025.129836.
- **Attack.** Exhaustive orientation search for n ≤ 9, then SAT or ILP. A proof looks out of reach.

<a id="b5"></a>
### Clar–Fries question for benzenoids

- **Statement.** Two questions on the relation between Clar and Fries numbers, answered negatively for fullerenes and open for
  benzenoids.
- **Source.** Fenton, Hartung, Art Discrete Appl. Math. (online January 2026); Fowler et al., ADAM 9 (2026) #P1.02, questions 1 and 2.
- **Attack.** ILP or SAT search over generated benzenoids up to 30–40 hexagons. A negative search result gives little, so this
  only pays off if a counterexample exists.

<a id="b6"></a>
### Minimal order of explosive endotactic networks

- **Statement.** What is the smallest order of an endotactic stochastic mass-action system that explodes?
- **Source.** Xu, [arXiv:2409.05340](https://arxiv.org/abs/2409.05340) Remark 4.8, with an example of order 10. Order 1 is
  excluded by Xu, arXiv:2601.00176 (every first-order endotactic system is exponentially ergodic).
- **Status.** Open; nothing lowers the order 10. Read arXiv:2511.10223, Zenodo 10.5281/zenodo.22728357, AIM problem 272.
- **Attack.** Lower the order of the example step by step and look for an obstruction at order 2 or 3. An exact value needs
  both directions, and only one is clearly accessible.

<a id="b7"></a>
### Gelation threshold at homogeneity 1

- **Statement.** For coagulation kernels of the form (xy)^{1/2} log^{α/2}(e + x) log^{α/2}(e + y), decide gelation for α in (1, 2].
- **Source.** Fournier, C. R. Math. 363 (2025) 583, doi:10.5802/crmath.738, [arXiv:2501.03867](https://arxiv.org/abs/2501.03867) §3, p. 4.
  The critical power depends on the kernel (Props 6, 7).
- **Status.** Partly addressed in a neighbouring setting: Ali, arXiv:2607.26167, proves gelation for all α > 1 in the discrete
  Oort–Hulst–Safronov system, with an additive-type kernel. It does not settle Fournier's continuous kernels, but its method
  must be read first.
- **Attack.** Transfer Ali's argument, or a new moment estimate at the critical logarithmic scale.

<a id="b8"></a>
### Nondegeneracy conjecture

- **Statement.** If the capacity for positive steady states of a network is finite, it equals its capacity for nondegenerate ones.
- **Source.** Joshi, Shiu, SIADS 16 (2017), doi:10.1137/16M1069705, [arXiv:1603.06441](https://arxiv.org/abs/1603.06441), Conj. 2.3,
  restated as open in Feliu, Henriksson, Pascual-Escudero, arXiv:2412.17798.
- **Status.** Open; partial cases known. Read arXiv:2608.18509 and 2507.09272.
- **Attack.** Counterexample search only.

## Tier C: very hard, kept for an idea

<a id="c1"></a>
### Brownian similarity profile uniqueness

See [its folder](../problems/brownian-coagulation-profile-uniqueness/). A complete proof must exclude every other profile at the
physical kernel, where perturbative theory stops. It needs a new global argument, not more estimates.

<a id="c2"></a>
### Random-jump heptagons

Exact knot probabilities for seven uniform points on the sphere joined cyclically (Geisz, Peterson, Shonkwiler,
[arXiv:2609.24860](https://arxiv.org/abs/2609.24860) §5). Unlike six points, the convex hull does not determine the knot type,
and the exact value probably has no closed form accessible by the same method.

<a id="c3"></a>
### Fullerene duals, Conjecture 1(c) and 1(d)

Only a counterexample looks reachable. One run proved part (b) and left these two open; read
[all-n/result.md](../problems/fullerene-dual-cospectral/all-n/result.md) first (not verified).

- For n ≥ 60, read as "one spectrum separates all isomers", (d) is equivalent to (c): both say that σ(T⁶) is injective.
- The exhaustive scan covers n ≤ 170, and exact T⁶-cospectral classes occur only at n = 32, 36, 40, 52.
- A(T⁶) − D is exactly cospectral for colour twins, so T⁶ misses only by a potential on the flanking hexagons. This is the
  place to look for a mechanism.

## Attempted, unresolved

<a id="d1"></a>
### Randić index and radius

- **Statement.** R(G) ≥ r(G) − 1 for every connected graph: Conjecture 1.1 of Cygan, Pilipczuk, Škrekovski, MATCH 67 (2012)
  451–466 (open access), who prove it for maximum degree at most 4. Stronger form: R(G) ≥ r(G) except for even paths, attributed
  to Caporossi and Hansen by Deng, Tang, Zhang, Filomat 29 (2015), doi:10.2298/FIL1506369D,
  [arXiv:1210.2543](https://arxiv.org/abs/1210.2543). The failing paths are exactly P₄, P₆, …, with R − r = √2 − 3/2.
- **Status.** Open on 31 October 2026 (bounded check). Read Doig, arXiv:2309.11613, and chapter 12 of Farmer's 2024 thesis (UNC
  Greensboro), which depends on an unproved conjecture.
- **One GPT-6-Astra run (xhigh), no resolution.** Its checks are in
  [verify_obstructions.py](../problems/randic-radius-conjecture/scripts/verify_obstructions.py) and its full notes in a
  git-ignored `PROGRESS.md` in that folder. None of this is verified.
  - Works: contracting a closed neighbourhood N[v] lowers the radius by at most one, for every v.
  - Missing: a vertex whose contraction always lowers R by at least 1 when the maximum degree is at least 3. With the metric
    lemma, the paths and the cycles, this would give R ≥ r + √2 − 3/2 by induction.
  - Refuted, with explicit graphs: choosing a vertex of maximum degree, even a unique one (R can increase); choosing the largest
    incident weight; spanning trees with at most 2R internal vertices; a block estimate R ≥ diameter.
- **Next attack.** A global argument that selects or averages contractions, or a different weighting. More local rules or
  more small graphs are not a demonstrated route, so a further run is worthwhile only with a new idea.

## Rejected during scouting

Solved or unusable, so not to be proposed again: fully distributive dual futile cycle (Vassena arXiv:2608.27081; Zenodo
10.5281/zenodo.23438810); positive recurrence of weakly reversible networks (analogue of the Global Attractor Conjecture,
AIM problem 272); Randić index against diameter (arXiv:1104.0426); maximal Wiener index of cubic graphs and of theta graphs;
Kaihnsa's attainable-region conjecture (arXiv:1810.03547); instantaneous gelation for (ij)^a; the differential sedimentation
kernel (arXiv:2506.11573); equilateral stick numbers of the trefoil and figure-eight; the Wang–Sontag 2n − 1 bound
(Zenodo 10.5281/zenodo.22770847); Vassena–Stadler Conjecture 17 (palomar-p003); recent Sombor index conjectures.
