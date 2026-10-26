# Earlier attempts retained for provenance

These are historical execution summaries, not the final verdict. The initial run had
14 jobs. The revised partial run rechecked the corrected winding/contour implementation
and added onset certificates. Its three coefficient jobs failed only because the requested
upper endpoint 0.0788111 exceeded the exact stored decimal 0.078811099999999995;
all mathematical rows passed. The coverage round corrected that request and all three
jobs passed. Both endpoints exceed 1/sqrt(161), so the integer range is unchanged.
The regression round ran the two newly registered regression/formula jobs.

The early onset-j3 trial failed a chart-transition containment assertion: a whole M-box
was unnecessarily required to fit in the last s-box. The final checker instead certifies
a narrow transition point in both boxes, and all three full onset runs passed.

The parent summary assembles all 21 required successful jobs with their run provenance;
it does not claim they were executed in a single invocation. Source and artifact hashes
are in the parent manifest. Earlier failed checks were not silently relabelled as passes.
