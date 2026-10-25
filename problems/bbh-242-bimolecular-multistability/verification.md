# Verification, 25 October 2026

**Verdict: CONFIRMED.** Up to the stated dynamical equivalence, exactly the ten
networks identified by Banaji, Boros and Hofbauer admit two asymptotically stable
equilibria in the closed quadrant. The stable pair must be the origin and a
positive equilibrium; no third stable equilibrium is possible.

Verifier: GPT-6-Astra (effort xhigh). Solver: Claude Opus 5.5 (effort xhigh).
The review used the handover summary, the published source (especially Table 1,
Lemma 27 and the question in Section 7), and the deliverables; no original solver
session was used. This is an AI review, not human peer review.

The boundary arguments were checked including nonhyperbolic equilibria. A
positive-axis equilibrium forces invariance of that axis. Its scalar quadratic
has at most one attracting positive root. Two attracting axis equilibria would
consume all four reactions and force positive transverse growth. An attracting
axis equilibrium and an attracting positive equilibrium are excluded by the
remaining one- or two-source cases and the determinant identity b p(x).

The index argument was re-derived. For an origin attracting relative to the
quadrant, the displacement homotopy (flow_t(z)-z)/t is nonzero on a small quarter
circle for every positive t. Its endpoint directions remain in their respective
closed half-planes. At a uniformly large time it points inward, giving the
angle correction -2 pi. On the opposite edge of a large triangle the quadratic
field turns by at most 2 pi: a larger sweep would meet one line three times,
forcing the degree-two determinant polynomial to vanish identically. Adding the
two axis segments gives index sum n-1 with n <= 1. An attracting positive
equilibrium therefore requires another positive equilibrium.

The five original reproduction scripts passed. In addition,
[rederive_four_sources.py](verification/rederive_four_sources.py) enumerates the
reaction sets and reconstructs the four-source proof without importing any
solver module. Positive circuits are obtained from determinants and the affine
source relation from cofactors. It reproduces 27,195 rank-two networks, 838
nontrivial four-source classes, nine constant families and 30 simple-fold
classes. Exact polynomial root counts confirm the required sign relation.

The ten explicit rate vectors were checked at the origin and (1,1). The three
singular origins are covered by the strict local Lyapunov inequality for x+3y;
a Jacobian eigenvalue test alone would be insufficient there. The 30-class table
then excludes the other 20 classes by inflow or instability at the origin.

Commands, exit codes and outputs are in
[verification/recheck-2026-10-25](verification/recheck-2026-10-25/).
The numerical index sanity check passed, but is not part of the proof.
