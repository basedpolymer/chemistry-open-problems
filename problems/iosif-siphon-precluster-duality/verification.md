# Verification, 25 October 2026

**Verdict: CONFIRMED WITH CORRECTIONS.** Each explicitly listed equality,
refinement, cardinality and siphon formulation is refuted under the stated
hypotheses. The source leaves "duality relation" undefined: this verdict does
not exclude a different relation that has not been specified.

Verifier: GPT-6-Astra (effort xhigh). Solver: Claude Opus 5.5 (effort xhigh).
The review used the handover summary, arXiv:2603.08767v2 (Definitions 17 and 19,
Remark 18, Conjecture 22 and Question 24), the thesis definitions and the
deliverables; no original solver session was used. This is an AI review, not
human peer review.

The MIPS criterion follows directly by considering the fibre through a coordinate
vector: invariance of its singleton vertex support forces positive proportionality
of the two conservation columns. The weighted coordinate transposition proves
the converse. Conservativity ensures that the extreme conservation rays span
the full conservation space. Transposing the educt and product matrices then
identifies dual MIPS with fully coupled reaction classes.

For each initial same-source pair, all subsequently added edges have their two
row vectors in that initial pair's span. This invariant proves that iteration
cannot merge different components of the one-step construction. The refinement
and equality criteria follow by checking whether each initial pair is coupled
and whether every nonsingleton coupled class meets such a pair. No numerical
assumption is involved.

The three-cycle was rechecked directly: the flux and conservation cones each
have ray (1,1,1), both precluster partitions consist of singletons, and the only
nonempty siphon is the full set. The remaining counterexamples supply failures
in the reverse refinement direction and for containment in a minimal siphon.
The exact-arithmetic implementation was checked against these definitions.

All 132 example claims and the five motif checks pass. The five enumerations
were re-run with two workers into new files: the combined count is again 15,513
networks satisfying the hypotheses, with the same relation counts and no failed
criterion. Commands and outputs are in
[verification/recheck-2026-10-25](verification/recheck-2026-10-25/).

Corrections applied to the manuscript: the minimal-size observation is restricted
to networks with nonzero reaction vectors, since the general definitions allow
zero-vector loops; the enumeration is explicitly restricted to non-loop
reactions. The phrase claiming "generic" failures from a list of motifs was
replaced by the narrower statement actually demonstrated. None changes the
counterexamples, criteria or enumeration counts.

## Addendum, 27 October 2026

The status of the result was changed from "solved" to "partial result": the verdict above covers the refutation of
each precise formulation and the exact criteria, not the conjecture with "duality relation" left undefined. Section 5
of the manuscript gained a paragraph on the Global Attractor Conjecture (the duality preserves the dimension s, does
not preserve weak reversibility, and gives singleton dual preclusters for weakly reversible networks), and
`certify_examples.py` gained ten claims for it (142 in all, PASS). This addendum was not part of the verification
above: the new remarks have only been checked by the script, not by a second agent.

## Addendum, 28 October 2026: the rewritten manuscript

The manuscript was rewritten from scratch (the earlier report is kept in [archive/](archive/paper-precise-forms.pdf);
the verdict above covers it). The new results come from agents' working notes in [gac-duality/](gac-duality/), which
were checked before the writing by three fresh-context verifiers (Claude Opus 5.5, effort xhigh), each re-deriving
the proofs and re-checking every example in exact arithmetic with its own code:

- [V1](gac-duality/verification/V1/verification.md), Petri dictionary and corrected preclusters: no claim refuted;
  ten corrected statements.
- [V2](gac-duality/verification/V2/verification.md), equal-conflict and proportional-conflict theory: every proof
  confirmed; Conjecture CR of the notes refuted (network NC3, three positive steady states in one class); five
  corrections.
- [V3](gac-duality/verification/V3/verification.md), persistence, obstructions and the results of agent A2: no claim
  refuted; Theorem O (periodic orbit of a persistent dual) confirmed; one gap repaired in the case used (D2's
  Proposition E, μ = 1); Conjecture F of the notes refuted (agent A2).

The writer listed every merged, re-proved or dropped statement in
[gac-duality/writer/CHANGES.md](gac-duality/writer/CHANGES.md). Two fresh-context reviewers (Claude Opus 5.5, effort
xhigh) then checked the manuscript: [R1](gac-duality/review/R1/review.md) for the mathematics (no theorem false;
three errors, sixteen corrections and fifteen missing proofs or pointers) and [R2](gac-duality/review/R2/review.md)
for attribution, consistency and presentation (five errors, including citation numbers and an absolute path that
pdfTeX embedded in the PDF; five gaps, including three earlier works on transposed networks). R1 applied both
reports ([applied.md](gac-duality/review/R1/applied.md)).

Final state: `certify_duality.py` 130/130 PASS and `certify_examples.py` 142/142 PASS (Windows and WSL); 41 pages, no
undefined reference. Not checked beyond their worked examples: the general tropical certificates (CM), (L), (S) of the
rank-two classification, which the manuscript labels as sketches; the large enumeration counts are checks, not
proofs. This is an AI review, not human peer review.
