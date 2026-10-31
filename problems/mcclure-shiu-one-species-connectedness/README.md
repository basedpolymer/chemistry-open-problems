# McClure–Shiu: are multistationarity regions of one-species networks with four or five reactions connected?

<!-- plain-words:start -->
**The question.** A reaction network with one species can settle at several concentrations, depending on its rate constants. McClure and Shiu asked whether the rate constants giving at least two steady states always form one connected piece when the network has four or five reactions. They had proved this for three reactions and found a counterexample with six.

**The answer.** Yes, always, according to a preprint by Patrick White dated 7 September 2026. It settles every size: the region splits only when the network has exactly three source complexes, each with reactions going both ways, and then into exactly two pieces. It is unrefereed. This folder does not contain a new proof: it records the claim and an audit of it, which found no error.

**How.** Rate constants matter only through the coefficients of the steady-state polynomial. Polynomials with two prescribed roots form a convex set, and rescaling the concentration moves the roots, so one connected family sweeps out the whole region. With three two-way sources that convex set is a line, and the zero polynomial cuts it in two.

<p align="center"><img src="figure.svg" alt="Polynomials with roots 1 and 2: one straight family for five reactions, two opposite rays for three two-way sources" width="720"></p>

<p align="center"><sub>With five reactions, polynomials sharing the steady states 1 and 2 pass from one coefficient sign pattern to another along a straight line. With three two-way sources, only a polynomial and its negative remain, and they meet only at zero.</sub></p>

<!-- plain-words:end -->

**Status: already solved (claimed in an unrefereed preprint). No manuscript here, by PROTOCOL §1.3.** The audit below found
no error, but it is the solver's own audit, not the separate verification of PROTOCOL §4.

**Attribution.** Claimed proof: Patrick White (independent researcher), with a large-language-model assistant according
to the paper's Appendix A.3; the commit trailers of his repository name Claude Opus 4.8, Claude Fable 5.1 and Claude
Opus 5. Novelty check and audit in this folder: Claude Opus 5.5 (effort not recorded). Separate verification: pending.

## Problem

A. McClure and A. Shiu, *On the connectedness of multistationarity regions of small reaction networks*, J. Symbolic
Comput. 125 (2024) 102323, doi:[10.1016/j.jsc.2024.102323](https://doi.org/10.1016/j.jsc.2024.102323),
arXiv:2303.03960v2, Section 6, p. 14:

> "Indeed, it is an open question whether, for one-species networks with four or five reactions (e.g., {0 ⇆ A, 2A ⇆ 3A}), the multistationarity regions are always connected."

A network is a finite set of reactions mᵢA → pᵢA. The region is the set of positive rate vectors κ ∈ ℝʳ for which
f_κ(x) = Σᵢ κᵢ(pᵢ − mᵢ)x^{mᵢ} has two or more positive roots (Definitions 3.3 and 4.2). A reversible pair counts as two
reactions, so the quoted example has four. The source proves connectedness for r ≤ 3 (Theorem 5.3) and gives a
disconnected example with r = 6 (Proposition 5.1).

## Who solved it

Patrick White, *Connectedness of multistationarity regions of one-species reaction networks*, draft of 7 September 2026,
14 pages, [PDF](https://read.byclaude.net/crnt-connectivity.pdf), with code and a Lean 4 formalization of Lemma 3.1 at
[github.com/pw/crnt-connectivity](https://github.com/pw/crnt-connectivity) (first commit 25 June 2026, last 9 September
2026). Theorem 1.1:

> "The multistationarity region M(G) is connected unless G has exactly three distinct source complexes and each of them is mixed. In that case M(G) has exactly two connected components."

A source is mixed if it has both a reaction producing A and one consuming A. Three mixed sources need r ≥ 6, so the
region is connected for r = 4 and r = 5, and for r = 6 it is disconnected exactly when there are three sources, each
with one reaction in each direction. Every set is proved path-connected. The paper defines M(G) with f_κ ≢ 0.

The proof has three parts: a linear reduction to coefficient space with convex fibres (Lemma 2.3); path-connectedness of
each sign region, by a continuous section of the fibres "roots 1 and τ" (Lemma 3.1); and an orthant graph that glues the
sign regions (Propositions 4.4 and 5.2). Proposition 6.1 separates the exceptional case by the sign of c₃.

No refutation or erratum was found, and no other proof. The arXiv, Zenodo and Crossref do not list this paper. The
acknowledgements thank Anne Shiu for comments on two earlier versions; no public statement from her was found. See
[references.md](references.md).

## How the answer depends on the reading of "region"

| Meaning of "region" | r = 4, 5 | Three mixed sources | Where |
|---|---|---|---|
| f ≢ 0, two distinct positive roots | path-connected | exactly two components | White, Theorem 1.1 |
| Literal Definition 3.3 (f ≡ 0 counts) | path-connected | path-connected, through f ≡ 0 | Last cell: White, Remark 2.1. First cell: below |
| Two nondegenerate steady states | path-connected | exactly two components | Complement below, not in White's paper |

The literal reading is not idle. At f ≡ 0 every x > 0 is a steady state, and Joshi–Shiu (their Example 2.4 and
Theorem 3.6) count that as multistationary. Then Proposition 5.1 of the source fails as written: on the segment from
κ(U) to κ(V) the polynomial is (1 − 2t)x(x − 1)(x − 2), and the midpoint (2, 2, 5/2, 5/2, 3/2, 3/2) has f ≡ 0, so the
whole segment lies in the literal region. Lemma 3.5(1) has the same mismatch, with {A → 0, A → 2A}. Both are correct once
f ≡ 0 is excluded. White's Remark 2.1 makes this point. The answer for r = 4, 5 does not depend on it: if all sources are
mixed and r ≤ 5, there are at most two sources, M(G) is empty, and the literal region is one convex fibre.

## What was audited here

Nothing below reuses or runs the author's scripts, which were not read. Lemma numbers are White's.

| Part of the paper | Check here | Result |
|---|---|---|
| Lemma 2.3, reduction to coefficients | Read; the same reduction had been derived here before the paper was found | no error found |
| Lemma 3.1, Steps 1–7 | Read line by line; formulas (14)–(18) tested in exact rationals on 2,682 random cases with s ≤ 8 | no error found, PASS |
| Lemmas 4.2, 4.3, Proposition 4.4 | Read line by line | no error found |
| Lemma 5.1, Proposition 5.2 | Read line by line, induction redone; orthant graph reimplemented from Definition 2.6, all 9,840 signatures with s ≤ 8 | no error found, PASS |
| Proposition 6.1, Theorem 1.1, Corollary 1.2 | Read line by line | no error found |
| The answer for r = 4, 5 | No signature realisable with r ≤ 5 is exceptional; all 4,845 and 15,504 networks on complexes 0A to 4A | PASS |
| The classification itself | Second proof, below, without the orthant graph | agrees |
| `CRNTConnectivity.lean` | Read as text: `region` and both theorem statements match Lemma 3.1, and `sorry` occurs only in the comment "sorry-free". **Not built**, as no Mathlib is available here | consistent |

**A shorter route, found during the audit.** The paper works orthant by orthant and then glues. That step is not needed.
Let C be the whole box of allowed coefficient vectors and F(τ) = {d ∈ C : P_d(1) = P_d(τ) = 0}. This set is convex, and
M_G is the image of {(u, τ, d) : u > 0, τ > 1, d ∈ F(τ) \ {0}} under rescaling, as in White's (13). If some source is not
mixed, 0 ∉ C and F(τ) \ {0} = F(τ) is convex. If all sources are mixed and s ≥ 4, F(τ) is a linear subspace of dimension
s − 2 ≥ 2, so removing 0 leaves it connected. If the region is not empty, Descartes' rule gives an achievable pattern
with two sign changes, and White's section s(τ) in that orthant joins the fibres exactly as in his Step 7. For s ≤ 2 the
region is empty. For s = 3 with all sources mixed, F(τ) is a line and F(τ) \ {0} is two rays: the exception.
Degenerate points need no separate treatment, and Sections 4 and 5 of the paper can be dropped. This simplifies White's
proof and gives a check on Proposition 5.2. It is not a new result. Two worked segments are in `audit_record.py`, checks
D1 and D2. For example, for {A → 0, A → 2A, 2A → 3A, 3A → 2A, 4A → 5A}, the polynomials x(x + a)(x − 1)(x − 2),
−3 ≤ a ≤ 1/2, form one straight segment of rate constants crossing from pattern (+, +, −, +) to (−, +, −, +).

**Complement: the nondegenerate reading.** M(G) is not open and contains degenerate points, and White's paper does not
treat the set N of coefficient vectors with two simple positive roots. Let A be the open set of c with an alternation
P_c(x₁)P_c(x₂) < 0, P_c(x₂)P_c(x₃) < 0, split into the types (+, −, +) and (−, +, −). For fixed x₁ < x₂ < x₃ and fixed
type the condition is convex. It is satisfiable for every triple when that sign pattern is achievable, and then that
type is path-connected; by Descartes' rule it is empty otherwise. If both patterns are achievable and G is not
exceptional, an alternating pattern of length four is achievable (check B2, complete for s ≤ 5 and so for r ≤ 5, and by
comparing the two index triples in general), and a four-term polynomial with three simple roots has both types. Finally
N ⊂ A, and A \ N lies in the triple-root locus, of dimension at most s − 2, which cannot disconnect an open subset of ℝˢ.
So N is path-connected for r ≤ 5. This paragraph has had no separate check.

**Limits of this audit.** The hand reading was done by one model, in one session, after it had derived the same
classification itself, so it may share the paper's blind spots. The Lean proof was not compiled. Exhaustive checks up
to s = 8 or s = 12 prove nothing beyond those sizes, although s ≤ 5 is all that r ≤ 5 needs. The preprint has not been
peer reviewed.

**Corrections to the brief.**

* The even-case claim fails for four or more sign changes: (x − 1)²(x − 2)² has two distinct roots but never changes
  sign. Such points lie in the region and need attaching, not only the double-plus-simple points of the odd case.
* The difficulty does not lie in sources with opposite net signs. The fibres are convex across sign patterns.
* The region is not open in general. For {0 ⇆ A, 2A ⇆ 3A} and κ = (2, 5, 4, 1), f = −(x − 1)²(x − 2) has two roots, and
  raising κ₃ by 1/100 leaves one. Remark 1.2 of the source covers only the regions computed there, so path-connectedness
  must be proved directly.
* Two citing papers (arXiv:2403.16556, 2607.25456) give the source as "A. Dennis and A. Shiu". Searches should use both
  names.

## Reproduce

From this directory (Python 3.12, standard library; `py -3` on Windows if `python` is missing):

```bash
python scripts/audit_record.py     # about 5 s, exact arithmetic, 17 PASS/FAIL checks
python scripts/make_figure.py      # figure.svg, with 2 exact checks of what it prints
python scripts/novelty_audit.py    # needs network and gh; results change with time
```

## Files

| File | Content |
|---|---|
| `references.md` | Sources read and novelty audit |
| `figure.svg`, `scripts/make_figure.py` | Figure and its script |
| `scripts/audit_record.py`, `certificates/audit-record.log` | Exact audit: source statements, combinatorics, Lemma 3.1, worked segments |
| `scripts/novelty_audit.py`, `certificates/novelty-api.json` | API part of the novelty search |
