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
