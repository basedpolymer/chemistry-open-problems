# Iosif: are preclusters dual to maximal invariant polyhedral supports, and to siphons?

**Status: solved (preliminary), with a negative answer.** Every precise formulation of the conjectured "duality
relation" that we could extract is false, and we prove exact criteria for when each relation holds. The conjecture
itself does not define "duality relation", so a relation of a different kind is not excluded (see "Open" below). Not
yet verified.

**Attribution.** Solver: Claude Opus 5.5 (effort xhigh), 25 Oct 2026, covering the novelty check, proofs, scripts and
manuscript. Verification: pending.

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

The abstract and the introduction state both as conjectures: "we conjecture that there is a duality relation between
preclusters and maximal invariant polyhedral supports. Given the close relation between maximal invariant polyhedral
supports and siphons, we also conjecture that siphons and preclusters are dual objects."

Definitions used: siphons (Def. 13), maximal invariant polyhedral supports (MIPS, Def. 17 and Remark 18), preclusters
(Def. 19; the iterated version is Algorithm 4.4 of the author's thesis, DOI 10.25673/13949), the dual network
(x^Ye →k x^Yp)* = k^(Yeᵀ) →x k^(Ypᵀ) (§3), and internal cycles (Def. 20). The statement has changed status over
time. The 2020 HAL v1 states it as a theorem, without proof and without the positive-steady-state hypothesis. A 2024
talk and a 2025 poster abstract (Ananiadi and Iosif) say it is proved. The 2026 preprint calls it a conjecture.

## Result

**Making the conjecture precise** (manuscript §1):

- The hypotheses are self-dual. N is conservative if and only if its dual N* is consistent, and vice versa (Γ(N*) = Γᵀ).
  "At least one positive steady state" for some rate constants means consistency. "No two species with exactly the
  same rates" means that the rows of Γ are pairwise distinct.
- Definition 17 (each decorated invariant polytope invariant under swapping xᵢ, xⱼ) is equivalent to Remark 18. Both
  say that xᵢ and xⱼ are in one MIPS if and only if their columns in a conservation-law matrix are positively
  proportional, if and only if eᵢ − μeⱼ ∈ im Γ for some μ > 0. So MIPS partition the species. The MIPS of N* are the
  classes of fully coupled reactions (flux coupling analysis).
- The one-step precluster graph of the paper (Def. 19) and the iterated one of the thesis (Alg. 4.4, in any order)
  have the same connected components.

**Exact criteria** (Theorem 7, Corollaries 8–10). For a consistent N:

- the preclusters P(N) refine the MIPS of N* if and only if every pair of reactions with a common source is fully
  coupled;
- P(N) = MIPS(N*) if and only if, in addition, every class of at least two fully coupled reactions contains such a pair.

The species side is obtained by applying this to N*, where "same source" becomes "identical rows of Ye". For weakly
reversible N, P(N*) always consists of singletons. For weakly reversible N with deficiency zero, P(N) ≠ MIPS(N*)
always.

**Counterexamples** (Theorem 1, Table 11). Every network below satisfies all hypotheses of Conjecture 22 and has an
explicit positive steady state.

| Network | Refutes |
|---|---|
| T: X₁ → X₂ → X₃ → X₁ (weakly reversible, deficiency 0, self-dual, rows/columns of Γ pairwise linearly independent) | equality, "MIPS refine preclusters", equal cardinality, on both sides and in all three readings of MIPS; all siphon formulations except "every precluster lies in a minimal siphon"; Remark 23 |
| R: X₁ ⇄ X₂ (smallest) | the same as T |
| W: X₁ → X₂, X₁ → X₃, X₂ → X₁, X₃ → X₁, and its dual W* | "preclusters refine MIPS" (reaction side for W, species side for W*) |
| U: X₄ → X₂+X₃, X₃+X₄ → X₁, X₁+X₂ → 2X₄ | "preclusters refine MIPS", species side, with the strongest hypotheses on N and N* |
| D: X₁ → X₂, X₁ → X₃, X₂+X₃ → 2X₁ (self-dual) | Question 24, "every precluster lies in a minimal siphon" (both sides), although both equalities of Conjecture 22 hold for D |
| F: one-site futile cycle | "preclusters refine MIPS" on both sides; "every precluster lies in a minimal siphon" |
| N₃ = 2X₁ ⇄ X₁+X₂ → 2X₂ and 2X₁ ⇄ X₁+X₂ (thesis Examples 4.1 and 4.12) | equality and refinements on the reaction side; Remark 23 |

The siphon formulations of Question 24 are: precluster = minimal siphon, refinement in both directions, "every
precluster is a siphon", "every precluster is a union of minimal siphons", and equal numbers. All fail on both sides
of the duality.

**What is computer-assisted.** All claims about the explicit networks have hand proofs in the manuscript. They are
re-checked by exact rational computation (132 claims). The exhaustive enumeration of 15 513 small networks satisfying
the hypotheses confirms the criteria and gives the statistics of Table 15; it is a check, not part of the proofs.
Nothing is numerical (no floating point anywhere).

**Open.** "Duality relation" is not defined in the source. Any relation that forces a refinement or a bijection
between blocks is excluded, but a relation of another kind is not. The rate-dependent clusters of the thesis were only
checked on N₃ (where the cluster version also fails). The motivation of the source, a case of the Global Attractor
Conjecture in small codimension, is not addressed.

## Reproduce

From `scripts/`, with Python 3.12 (standard library only; run under WSL Fedora 44):

```bash
python certify_examples.py      # 132 exact claims on the explicit networks            (< 1 s)
python motifs.py                # Table 14: standard biochemical motifs                  (< 10 s)
python enumerate_small.py 2 2 5 12 ../certificates/enumeration_n2_m2_r5.json    # (n, molecularity, r_max, cores)
python enumerate_small.py 3 2 4 12 ../certificates/enumeration_n3_m2_r4.json    # ~75 s on 12 cores
python enumerate_small.py 3 3 3 12 ../certificates/enumeration_n3_m3_r3.json    # ~50 s
python enumerate_small.py 4 1 12 12 ../certificates/enumeration_n4_m1_r12.json  # ~20 s
python enumerate_small.py 4 2 3 12 ../certificates/enumeration_n4_m2_r3.json    # ~5 s
python summarize_enumeration.py ../certificates/enumeration_*.json             # Table 15
```

Each script prints PASS/FAIL per claim and exits with a nonzero status on failure. The manuscript is built with
`python templates/build_tex.py problems/iosif-siphon-precluster-duality/paper.tex` from the repository root.

## Files

| File | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript: definitions, Proposition 3 (MIPS), Lemma 5 (precluster variants), Theorem 7 and Corollaries 8–10 (criteria), Theorem 1 (counterexamples), enumeration |
| `scripts/crn.py` | Library: exact extreme rays, dual network, siphons, three precluster constructions, three readings of MIPS |
| `scripts/certify_examples.py` | Exact certificates for every network of the manuscript, including Examples 15 and 25 of the paper and Examples 1.22, 4.1, 4.5 and 4.12 of the thesis |
| `scripts/motifs.py` | Standard motifs (Michaelis–Menten, futile cycles, thesis network G, a triangle) |
| `scripts/enumerate_small.py` | Exhaustive enumeration with the check of the criteria |
| `scripts/summarize_enumeration.py` | Combined statistics (Table 15) |
| `certificates/*.log`, `certificates/*.json` | Outputs (all PASS), enumeration statistics and smallest failing networks |
| `references.md` | Sources read and novelty check |
| `NEEDS.md` | One optional paywalled source (Conradi–Flockerzi 2012) |

`certificates/sources/` holds the arXiv, HAL, thesis and slide PDFs. It is git-ignored because the files are
copyrighted.
