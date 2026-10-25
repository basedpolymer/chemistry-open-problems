# Verification, 25 October 2026

**Verdict: CONFIRMED.** The four closed benchmark values are 37sp-yfyv = 36,
20sp1 = 19, 23sp1 = 22 and unbalanced10 = 25. The lower bound 32 for the
balanced15 family is also confirmed. This does not close the other benchmarks.

Verifier: GPT-6-Astra (effort xhigh). Solver: Claude Opus 5.5 (effort xhigh).
The review used the handover summary, original benchmark data and published model,
the manuscript, and the certificate files; no original solver session was used.
This is an AI review, not human peer review.

The block bound follows by taking the connected components of the bipartite match
graph. Every component is balanced and obeys each cumulative heat inequality;
with c components it needs at least N-c edges. The exact subset enumeration and
the three-block exclusion were re-run. The matching upper bounds are explicit
networks, checked with exact flows and separately as feasible points of the
2018 MILP with each of its three big-M data files.

For unbalanced10 the strengthened MILP was audited row by row. The residual upper
bound is the total heat cascade; direct and match capacities follow from flow
cuts. The cold-window and hot-window inequalities count, respectively, delivered
heat and heat originating in a window. Each partner-count row follows by sorting
nonnegative coefficients. Thus every original feasible network with at most 24
matches maps to a feasible point of the certified model.

The 540,587,890-byte VIPR certificate was rechecked from its recorded SHA-256:

- the LP regenerated identically, and every row and bound was re-derived;
- the certificate problem matched a relaxation of that LP;
- the hardened VIPR checker accepted both the raw and completed certificate;
- the separate Python checker checked all 191,697 derivations and 14,796,916
  multipliers, including assumptions, integer rounding and disjunctions, and
  obtained a contradiction with no remaining assumptions;
- both checkers rejected all four mutated certificates.

All six reproduction steps exited successfully. Commands, exit codes, timings
and full outputs are in [verification/recheck-2026-10-25](verification/recheck-2026-10-25/).
The VIPR file remains outside Git, as documented in the certificate manifest.
The floating-point explorations and the unpatched VIPR checker are not used to
justify the lower bound.
