# Iosif's network duality: preclusters, invariant polyhedral supports, siphons, and what the duality transfers

**Status: partial result.** Every precise form of the conjectured "duality relation" is false, with exact criteria
and minimal counterexamples. The relation left open by the source is answered in three parts:
- the duality is the reverse dual of Petri-net theory;
- it transfers siphon structure, persistence, uniqueness and stability on an explicit self-dual class of networks;
- outside that class it transfers none of weak reversibility, persistence or convergence.

The conjecture as written is not a decidable statement, and the rank-two persistence dichotomy is open; see
[What remains open](#what-remains-open). Manuscript: [paper.pdf](paper.pdf).

**Attribution.**
- Research: Claude Opus 5.5 (effort xhigh), for the original refutation (25 October 2026) and agents P, E, Q1, Q2,
  Q3, E2 and D2 (27 October 2026). GPT-6-Astra (effort xhigh), for agents G and A2.
- Literature: Claude Opus 5.5 (effort xhigh; effort medium for the second search).
- Verification: GPT-6-Astra (effort xhigh) for the original report ([verification.md](verification.md)); Claude Opus
  5.5 (effort xhigh) for the new results, by three fresh-context verifiers V1–V3, and for the review of the rewritten
  manuscript (R1, R2).
- Writing of the manuscript (28 October 2026): Claude Opus 5.5 (effort xhigh).

**Acknowledgment.** We thank Alexandru Iosif for helpful correspondence, which shaped questions Q1–Q3 and pointed us
to recent work on the Global Attractor Conjecture.

## Problem

A. Iosif, *Duality in mass-action networks*, arXiv:2603.08767 (v1 9 March 2026, v2 15 July 2026; also HAL
hal-03058670 v3). Section 3, p. 7 of v2 (identical in v1):

> "Conjecture 22. For conservative mass-action networks that do not have two species with exactly the same rates and
> have at least one positive steady state, there is a duality relation between the sets of preclusters and of maximal
> invariant polyhedral supports.
>
> Remark 23. The condition ẋᵢ ≠ ẋⱼ for i ≠ j is necessary in order to make sure no further clustering on the left
> kernel is done. We require that there is at least one positive steady state so that the internal cycles are a
> refinement of preclusters.
>
> We conclude with the following question:
>
> Question 24. Is there a duality relation between preclusters and siphons?"

The motivation follows (p. 7–8): "It is well known that the Global Attractor Conjecture holds in small dimensions
[4], and we wonder whether there is a similar argument, through duality, in small codimensions." The dual network
exchanges species and rate constants: N* has one species per reaction of N and one reaction per species, with
stoichiometric matrix Γᵀ. Networks with matrix ±Γᵀ appeared before (Famili–Palsson 2003, Jamshidi–Palsson 2009,
Blanchini–Giordano 2022, Duvall 2024); they transfer linear and local properties (conservation pools, spectra,
structural stability), not siphons or persistence. The statement has changed status over time. The 2020 HAL v1 states it as a theorem, without
proof. A 2024 talk and a 2025 poster abstract say it is proved, and the 2026 preprint calls it a conjecture.

## Result

The manuscript asks three questions that the source supports, and answers them. Below, N is a network, N* its dual,
s the rank of Γ (the same for N and N*), q the number of Petri clusters, and (H1)–(H3) the hypotheses of Conjecture
22 (conservative, consistent, distinct rows of Γ).

**The precise forms fail** (§3.1, condensed from the [archived report](archive/paper-precise-forms.pdf)). Equality,
refinement in either direction, equal cardinality and the siphon versions are all false. Exact criteria say when each
holds. The counterexamples are tiny: X₁ ⇄ X₂, the 3-cycle X₁ → X₂ → X₃ → X₁, the one-site futile cycle.

**1. Preclusters and siphons are related by transversality, and preclusters can be repaired** (§2–3).
- Iosif's dual is the weighted *reverse dual* of Petri nets (Hack 1972): places and transitions are exchanged, every
  arc is reversed, weights are kept.
- In this dictionary, conservation laws of N are flux vectors of N*, and siphons of N are the T-traps of N*.
- A minimal siphon contains at most one species of each class of species with the same consumers. This extends
  Desel–Esparza Thm 4.30(a) for free-choice nets. So a precluster with two elements never lies in a minimal siphon.
- Definition 19 depends on how a choice is split into binary choices. The Conradi–Flockerzi clusters P̂ are
  invariant under duplication and splitting.
- For P̂, the maximal invariant polyhedral supports refine the preclusters for every consistent (respectively
  conservative) network. Equality holds exactly when twin pairs are coupled, in particular on the species side of
  every weakly reversible network satisfying (H1)–(H3).
- Remark 23 for P̂ is the isolation property, which gives toric steady states for all rates.

**2. The duality transfers structure on a self-dual class** (§4, Figure 46).
- **Cluster rank class.** If N is connected and satisfies (H1), (H2) and s = q − 1, the rank equation of free-choice
  and equal-conflict nets, then neither N nor N* has a critical siphon. Both are persistent for every rate vector.
  The siphon part is known for ordinary nets (Desel 1994); the weighted version was not found in the literature.
- **Proportional-conflict networks** (a self-dual subclass). N and N* each have a positive steady state for every
  rate vector, exactly one per class, and the Jacobian is Hurwitz on the stoichiometric subspace at every positive
  point.
- **Well-formed equal-conflict networks** (reactions sharing a reactant have the same reactant complex). N has toric
  steady states and is globally stable; N* is globally stable when s ≤ 2.
- **Boundary of the dynamical statements.** The network NC3 = {D → B+C, C → A, A+B → B+C, 2A+B → C+D} satisfies the
  cluster rank equation but is not proportional-conflict. It has three positive steady states in one class, one of
  them a saddle.

**3. Outside this class, the dynamics do not transfer** (§5–6).
- Weak reversibility, deficiency, endotacticity and the siphon criterion of Angeli–De Leenheer–Sontag are not
  preserved. Weak reversibility was already noted by Blanchini–Giordano 2022.
- **Persistence is not preserved.** The dual of 2X₂ ⇄ X₁+X₂ loses a species for every rate vector. Drain species and
  drain sets give such extinction in every dimension.
- **Convergence is not preserved, even when persistence is.** The 3-cycle 3X₃ → X₁+X₂ → X₂+2X₃+3X₄ → 3X₃ is weakly
  reversible with deficiency zero. Its dual is persistent for every rate vector and has a unique positive
  equilibrium in every class. Yet at the rates (17, 400, 70000, 10⁶) that equilibrium is a source, and the class
  contains a periodic orbit (Theorem 54). The certificate is exact; a figure illustrates the orbit.
- **Rank one.** Persistence of the dual is decided exactly by face types (IN, OUT, ESC, DOM). For weakly reversible N
  it is all-or-none over rate vectors. Only a weaker form of the "if" half is in print (Zhang–Gao–Dochain 2023).
- **Any rank.** For weakly reversible N, the dual has no dominated facet.
- **Rank two.** Among 18 861 networks, 4 873 have duals persistent for all rates and 13 817 for none; 171 are
  undecided. Two exact families refute the conjecture that a facet attracting the flow forces extinction.

**The Global Attractor Conjecture** (§7). The duality preserves the dimension s and exchanges the codimension n − s
with the cycle-space dimension of N*. So it cannot turn small codimension into small dimension. The GAC itself is
now claimed in full by Craciun (arXiv:1501.02860 v3) and by OpenAI's preprint of 25 September 2026; neither proof
was checked here.

**What is computer-assisted.** Every example and finite claim of the manuscript is certified in exact arithmetic
(`scripts/certify_duality.py`, 130 claims; `scripts/certify_examples.py`, 142 claims for the archived report). Counts
on large families (enumerations of agents E, E2, D2, recounted by the verifiers) are checks, not proofs. The figure of
the periodic orbit is a numerical illustration only.

## What remains open

1. The rank-two dichotomy: is the dual of a weakly reversible network with (H1)–(H3) and s = 2 persistent for all
   rate vectors or for none? 171 networks are undecided.
2. Global stability of the dual of equal-conflict networks for s ≥ 3, and global stability for
   proportional-conflict networks.
3. Absolute consistency under the cluster rank equation, and for duals of unbalanced equal-conflict networks with
   s = m.
4. A liveness theory for weighted impure equal-conflict nets.
5. An exact algorithm for P̂ when s ≥ 2, and whether P̂ is the finest invariant repair of Definition 19.
6. Which facets that attract the flow force extinction in higher rank.

## Reproduce

From `scripts/`, with Python 3 (standard library; `figures.py` needs scipy and matplotlib, run under WSL Fedora 44):

```bash
python certify_duality.py       # 130 exact claims of the manuscript                      (~1 s)
python certify_examples.py      # 142 exact claims of the archived report                 (< 1 s)
python figures.py               # Figure 55, numerical illustration of Theorem 54         (~16 s)
python motifs.py                # archived report, Table 14
python enumerate_small.py 3 2 4 12 ../certificates/enumeration_n3_m2_r4.json    # archived enumeration, ~75 s
python summarize_enumeration.py ../certificates/enumeration_*.json
```

Each certificate script prints PASS/FAIL per claim and exits with a nonzero status on failure. The manuscripts are
built from the repository root with `python templates/build_tex.py problems/iosif-siphon-precluster-duality/paper.tex`
(and `.../archive/paper-precise-forms.tex`). The working notes, their scripts and the verification reports are in
`gac-duality/`; each folder's notes state how to re-run its computations.

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript (41 pages): Petri dictionary, precise forms, corrected preclusters, transversality, what the duality transfers and does not transfer, persistence in rank one and two, the GAC, open problems, verification record and credits |
| `scripts/certify_duality.py`, `certificates/certify_duality.log` | Exact certificates for every example and finite claim of the manuscript |
| `scripts/figures.py`, `figures/`, `certificates/figures.log` | Numerical illustration of Theorem 54 (periodic orbit), in logarithmic coordinates |
| `archive/` | Earlier report (25 October 2026, revised 27 October): exact criteria and counterexamples for the precise forms (`paper-precise-forms.tex`, `.pdf`, `references-precise-forms.bib`), verified; certified by `scripts/certify_examples.py` |
| `scripts/crn.py`, `certify_examples.py`, `motifs.py`, `enumerate_small.py`, `summarize_enumeration.py`, `certificates/enumeration_*` | Library and certificates of the archived report |
| `gac-duality/` | Working notes of the agents (P, E, L, G, Q1, Q2, Q3, E2, D2, A2, L2) with their scripts and logs; verification reports `verification/V1–V3/`; prior-art table `literature2/novelty.md`; the writer's list of changes `writer/CHANGES.md`; reviews `review/R1–R2/` |
| `verification.md`, `verification/` | Verification of the archived report, with addenda for the new results |
| `references.md` | Sources read and novelty check |

`certificates/sources/` and the agents' `sources/` folders hold copyrighted PDFs; they are git-ignored.
