# Certified enclosure of the knotting probability of random equilateral hexagons

*Solver: Claude Opus 5.5 (effort xhigh). Date: 24 October 2026.*
*Status: partial result, preliminary (not yet verified; certificate spot-checked, see §6).*

## 0. Statement

Let **Pol̂(6)** be the space of closed equilateral hexagons in R³ (unit edges) modulo translations and
rotations, with its standard (symplectic) probability measure μ (Cantarella–Shonkwiler, Ann. Appl. Probab. 26
(2016) 549–596, §4 and Theorem 13: the symplectic volume of Pol̂³(n; 1) equals the standard measure; their
Theorem 1, the Duistermaat–Heckman theorem, makes the action-angle map measure preserving). Let p₆ = μ{knotted}.

> **Theorem (certified).**  L ≤ p₆ ≤ U with the exact dyadic rationals L, U of Section 6:
> 3.0675·10⁻⁵ < L = 3.0675613…·10⁻⁵ and U = 8.4944462…·10⁻⁴ < 8.4945·10⁻⁴.

For comparison: Hake (JKTR 28 (2019) 1950040; arXiv:1810.11882, Theorem 4.1) proves p₆ ≤ (14 − 3π)/192 ≈ 0.02383;
Cantarella–Shonkwiler (Proposition 16) prove p₆ ≤ 1/2; no positive lower bound was known. Monte Carlo
estimates: 1.3·10⁻⁴ ± 0.2·10⁻⁴ (Cantarella–Shonkwiler §4.1, 1.3·10⁶ samples), 1.37·10⁻⁴ (Hake §4),
and ours, 1.3639·10⁻⁴ ± 0.0117·10⁻⁴ (10⁸ samples; Section 7; numerics, not part of the proof).

What is *proved* (Sections 1–4, pen and paper, using only the classical facts F1–F3 below) is a reduction of
p₆ to the volume of an explicit semi-algebraic subset of a 6-dimensional box, together with a sufficient
condition for knottedness. What is *certified by computer* (Section 5) is an outer and an inner
approximation of that volume by finitely many boxes, with interval arithmetic (outward rounding,
no transcendental library functions) and exact dyadic summation. Calvo's characterisation of hexagonal
trefoils is **not** used in the proof (it is used only to interpret the two patterns R, L, Section 3.4).

### Classical facts used

* **F1 (triangle moves).** If the closed triangle conv(v_{i−1}, v_i, v_{i+1}) of a polygonal knot meets the
  polygon only in the edges [v_{i−1}, v_i], [v_i, v_{i+1}] and in the endpoints v_{i−1}, v_{i+1} of the two
  neighbouring edges, then deleting v_i gives an equivalent (ambient isotopic) polygonal knot
  (Reidemeister's Δ-moves; e.g. Crowell–Fox, *Introduction to knot theory*, Ch. I; Burde–Zieschang, *Knots*, §1.A).
* **F2 (geometric isotopy).** A continuous path in the space of embedded n-gons consists of equivalent knots
  (Calvo, JKTR 10 (2001) 245–267 = arXiv:math/9904037, §1: path components of Geo(n) are unions of knot types;
  standard by compactness and F1).
* **F3 (Fox 3-colourings).** The number of Fox 3-colourings of a regular diagram is a knot invariant; it is 3 for
  the unknot (and 9 for the trefoil) (Fox 1961; Crowell–Fox Ch. VI; Przytycki, Banach Center Publ. 42 (1998)).
  A generic linear projection of a polygonal knot is a regular diagram of it.

We also use that a planar (or degenerate) triangle is unknotted, trivially.

## 1. Action-angle coordinates and the reduction to a volume

**1.1 Coordinates.** Use the triangulation T135 of the hexagon (diagonals v1v3, v3v5, v5v1 of lengths d1, d2, d3).
The moment polytope is P = {d ∈ [0,2]³ : d_i ≤ d_j + d_k}, Vol(P) = 4. For (d, θ) ∈ int(P) × T³ put (standard position)

* v1 = (0,0,0), v3 = (d1,0,0), v5 = (p, q, 0), p = (d1² − d2² + d3²)/(2d1), q = √(d3² − p²) > 0;
* for k = 1,2,3 let [a_k, b_k] = [v1,v3], [v3,v5], [v5,v1] (oriented), m_k its midpoint, n_k the unit vector in
  the plane z = 0 obtained by rotating b_k − a_k by +90° (it points into the central triangle),
  h_k = √(1 − d_k²/4), and
  v_{2k} = m_k + h_k (cos θ_k n_k + sin θ_k e_z).

Then |v_{i} − v_{i+1}| = 1 and v_{2k} rotates at unit speed about the line a_k b_k when θ_k varies (θ_k = 0: flap folded onto
the central triangle; θ_k = π: planar convex hexagon). This is the action-angle map α of
Cantarella–Shonkwiler (Theorem 13: the angle θ_k acts by bending about the k-th diagonal, the moment map records
the diagonal lengths); the choice of the zero of each angle is a translation of T³, which preserves Haar measure.
By their Theorem 1 (Duistermaat–Heckman) and Theorem 13,

  p₆ = (1/(Vol(P)(2π)³)) ∫_P ∫_{T³} χ_knot(α(d,θ)) dθ dd = (1/(32π³)) Vol₆{(d,θ) : α(d,θ) knotted}.   (1)

**1.2 Genericity.** Call a hexagon *generic* if no four of its vertices are coplanar. Each of the 15 orientation
determinants O(a,b,c,e) = det(v_b − v_a, v_c − v_a, v_e − v_a) is a real-analytic function of (d,θ) on the
connected set int(P) × T³ and is nonzero at an explicit point (`scripts/check_geometry.py`, G1, ball arithmetic);
hence non-generic (d, θ) form a null set and may be ignored in (1). For generic hexagons, "the edge [p,q] meets
the triangle conv(a,b,c)" (vertex-disjoint) is equivalent to "meets its relative interior", and is decided by
the signs of O(a,b,c,p), O(a,b,c,q), O(p,q,a,b), O(p,q,b,c), O(p,q,c,a): it holds iff the first two have
opposite signs and the last three are equal. We say that [p,q] *pierces* the triangle.

## 2. Elementary lemmas

Indices are mod 6, e_j = [v_j, v_{j+1}], T_i = conv(v_{i−1}, v_i, v_{i+1}).

**Lemma 2.1 (shared vertex).** Let T = conv(x,a,b) and T' = conv(x,a',b') be non-degenerate triangles with only
the vertex x in common. If the segment [a',b'] meets the plane Π of T in exactly one point z and z lies in the
relative interior of T, then [a,b] ∩ T' = ∅.

*Proof.* z ∈ [a',b'] ⊂ Π', so z lies on the line ℓ = Π ∩ Π' (Π ≠ Π' since a', b' ∉ Π), and ℓ = line(x, z).
In Π', ℓ passes through the vertex x and the point z of the opposite side, so ℓ ∩ T' = [x, z]. In Π, ℓ passes
through the vertex x and the interior point z, hence leaves T through a point z* of the opposite side [a,b], with z
strictly between x and z*. Since [a,b] ⊂ Π and [a,b] ⊄ ℓ (otherwise x, a, b would be collinear),
[a,b] ∩ Π' ⊂ [a,b] ∩ ℓ = {z*}, and z* ∉ [x,z] = ℓ ∩ T'. ∎

**Lemma 2.2 (pentagons).** Every generic closed pentagon is unknotted.

*Proof.* In a pentagon only e_{i+2} can pierce T_i (all other edges share a vertex with T_i). T1 = conv(v2, v5, v1) and
T3 = conv(v2, v3, v4) share only v2, and the sides opposite v2 are e5 and e3. If e3 pierced T1, Lemma 2.1 would
give e5 ∩ T3 = ∅; so T1 and T3 are not both pierced. An unpierced T_i meets the rest of the polygon only in the
two neighbouring edges at v_{i±1} (genericity), so F1 deletes v_i and leaves a quadrilateral; in a quadrilateral
every triangle T_i is unpierced (all edges touch it), so F1 again reduces it to a triangle, which is unknotted. ∎

**Lemma 2.3 (knotted hexagons).** If a generic hexagon is knotted, every triangle T_i (i = 1..6) is pierced by
exactly one edge.

*Proof.* Only e_{i+2} and e_{i+3} can pierce T_i (e_{i−1}, e_i lie in ∂T_i; e_{i−2}, e_{i+1} meet the plane of T_i
only at v_{i−1}, v_{i+1}). If neither does, F1 deletes v_i and Lemma 2.2 applies. If both do, let Π_i, Π_{i+3} be
the planes of T_i and T_{i+3} = conv(v_{i+2}, v_{i+3}, v_{i+4}) (vertex-disjoint from T_i) and ℓ = Π_i ∩ Π_{i+3}.
The piercing points z, z' of e_{i+2}, e_{i+3} lie on ℓ, are distinct, and lie on ∂T_{i+3}; hence ℓ ∩ T_{i+3} = [z, z']
⊂ int T_i (convexity). Then ∂T_i ∩ T_{i+3} ⊂ ∂T_i ∩ ℓ ∩ T_{i+3} = ∅, so e_{i−1}, e_i (the only edges that could
pierce T_{i+3}) miss it; F1 deletes v_{i+3} and Lemma 2.2 applies. ∎

**Lemma 2.4 (half-spaces).** In standard position, if the numbers sin θ1, sin θ2, sin θ3 are nonzero and not all
of the same sign, the hexagon is unknotted (for generic (d,θ)).

*Proof.* z(v_{2k}) = h_k sin θ_k with h_k > 0. Some v_{2k} is alone on its side of the plane z = 0, say z(v_{2k}) > 0
> z(v_{2k±2}). The open triangle T_{2k} lies in z > 0, while the only edges that can pierce it, e_{2k+2} and e_{2k+3},
join v_{2k+3} (in z = 0) to v_{2k+2}, v_{2k+4} (in z < 0) and lie in z ≤ 0. So T_{2k} is unpierced; Lemma 2.3. ∎

(This is the argument of Hake's Lemma 3.4; it improves Cantarella–Shonkwiler's 1/2 to 1/4 already.)

**2.5 Symmetries.** (`check_geometry.py`, G3.) (a) α(d, −θ) is the mirror image of α(d, θ) in the plane z = 0.
(b) α(d2,d3,d1, θ2,θ3,θ1) is properly congruent to the relabelled hexagon (v3,v4,v5,v6,v1,v2) of α(d,θ).
(c) α(d3,d2,d1, θ3,θ2,θ1) is improperly congruent to (v1,v6,v5,v4,v3,v2). Relabelling, rigid motions and mirror
images preserve knottedness, so χ_knot ∘ α is invariant under θ ↦ −θ and under the simultaneous action of S₃ on
(d1,d2,d3) and (θ1,θ2,θ3). With Lemma 2.4 and (1):

  p₆ = (2·6/(32π³)) · K₀ = 3K₀/(8π³),   K₀ = Vol₆{(d,θ) ∈ C × (0,π)³ : α(d,θ) knotted},   (2)

where C = {2 ≥ d1 ≥ d2 ≥ d3, d1 ≤ d2 + d3} is a fundamental chamber (Vol C = 2/3).

## 3. The two pierce patterns

**3.1 Pair conditions.** Each flap apex v_{2k} depends on θ_k only. Define, on C × (0,π)³,

| name | condition | depends on |
|---|---|---|
| A12 | e4 = [v4,v5] pierces T2 = (v1,v2,v3) | d, θ1, θ2 |
| A23 | e6 = [v6,v1] pierces T4 = (v3,v4,v5) | d, θ2, θ3 |
| A13 | e2 = [v2,v3] pierces T6 = (v5,v6,v1) | d, θ1, θ3 |
| B12 | e1 = [v1,v2] pierces T4 | d, θ1, θ2 |
| B23 | e3 = [v3,v4] pierces T6 | d, θ2, θ3 |
| B13 | e5 = [v5,v6] pierces T2 | d, θ1, θ3 |

and R = A12 ∧ A23 ∧ A13, L = B12 ∧ B23 ∧ B13. Each of the five determinants of each condition is of the form
u(θ_i)ᵀ M(d) u(θ_j) with u(θ) = (1, cos θ, sin θ), because a determinant is affine in each of its rows and only
v_{2k} depends on θ_k.

**Proposition 3.2 (necessary condition).** For generic (d,θ) ∈ C × (0,π)³: knotted ⇒ R or L; and R ∧ L is empty.

*Proof.* By Lemma 2.3, T2 is pierced by exactly one of e4 (A12), e5 (B13); T4 by exactly one of e6 (A23), e1 (B12);
T6 by exactly one of e2 (A13), e3 (B23). Lemma 2.1 applied to (T2, T4) (common vertex v3; opposite sides e1, e4)
excludes A12 ∧ B12; applied to (T4, T6) (vertex v5; sides e3, e6) it excludes A23 ∧ B23; applied to (T6, T2) (vertex
v1; sides e5, e2) it excludes A13 ∧ B13. Of the 2³ choices only (A12, A23, A13) and (B13, B12, B23) avoid these. ∎

**Proposition 3.3 (embeddedness).** If all determinant signs of R (resp. L) hold strictly at (d,θ), the hexagon is
embedded; the same holds on any connected set where they hold strictly.

*Proof (R).* A12 strict: e4 meets the plane of T2 once, inside T2, so e4 ∩ (e1 ∪ e2) = ∅; Lemma 2.1 on (T2,T4) gives
e1 ∩ T4 = ∅, so e1 ∩ (e3 ∪ e4) = ∅. Likewise A23 gives e6 ∩ (e3 ∪ e4) = ∅ and e3 ∩ (e5 ∪ e6) = ∅, and A13 gives
e2 ∩ (e5 ∪ e6) = ∅ and e5 ∩ (e1 ∪ e2) = ∅. These are all nine non-adjacent pairs. Adjacent edges overlap only if
three consecutive vertices are collinear, which would make one of the strict determinants O(v1,v2,v3,v4),
O(v3,v4,v5,v6), O(v5,v6,v1,v2), O(T2,·), O(T4,·), O(T6,·) vanish. (L: symmetric, with B13, B12, B23.) ∎

**3.4 Remark (Calvo).** Calvo (arXiv:math/9904037, §2, citing Lemma 8 of Topology Appl. 112 (2001)) states that a
hexagon is a trefoil iff Δ2, Δ4, Δ6 are all nonzero (then all equal ±1). Hence R and L are, up to null sets,
exactly the right- and left-handed trefoils of positive curl; p₆ = 3(Vol R + Vol L)/(8π³). We do not use this;
`mc_validate.py` confirms it on 10⁸ samples (C2, C5, C5').

**Corollary 3.5.** p₆ ≤ (3/(8π³)) [Vol₆(R) + Vol₆(L)] (Prop. 3.2), and p₆ ≥ (3/(8π³)) Vol₆(S) for any set S ⊂ R ∪ L
on which every connected component, along which the strict signs hold, contains one point with a 9-colourable
diagram (Prop. 3.3, F2, F3).

## 4. Chart of the chamber

C is parametrised by (d1, s, t) ∈ [0,2] × [0,1]² via d2 = d1(1+s)/2, d3 = d1((1−s)/2 + ts), with Jacobian d1² s/2
and d1 − d2 = d1(1−s)/2, d2 − d3 = d1 s(1−t), d2 + d3 − d1 = d1 ts (`check_geometry.py`, G4). In these
coordinates q = d1 √Q / 2 with Q = ts(1 − s²(1−t)²)(2 + ts) ≥ 0 (no division by d1), and
n2 = (−√Q/(2r2), (r3² − r2² − 1)/(2r2), 0), n3 = (√Q/(2r3), −(1 − r2² + r3²)/(2r3), 0), r2 = d2/d1, r3 = d3/d1.
For a box D = [2a₁, 2a₁+2]/2^{L1} × [a₂, a₂+1]/2^{L2} × [a₃, a₃+1]/2^{L3} and θ-cells of side π/2^{Lθ}, the
contribution of N cells of D × (θ-grid) to (2) is exactly

  w(D)·N,  w(D) = ((a₁+1)³ − a₁³)(2a₂+1) / (4 · 2^{3L1 + 2L2 + L3 + 3Lθ})  (a dyadic rational; π cancels).   (3)

## 5. Certified computation (`src/hexcert.cpp`, `scripts/bnb_driver.py`)

For each box D and level Lθ (n = 2^{Lθ} cells per angle on (0,π)):

1. **Coefficients.** For each of the 30 determinants (6 conditions × 5), the 3×3 coefficient matrix M(d) is
   evaluated with interval arithmetic and first-order forward automatic differentiation over D: an enclosure
   M(D), the value M(c) at the centre c of D, and enclosures G_m(D) of ∂M/∂(d1,s,t)_m.
2. **Cell enclosures.** For a θ-cell Θ_i × Θ_j, F = u(θ_i)ᵀ M u(θ_j) is enclosed by the intersection of
   (a) u(Θ_i)ᵀ M(D) u(Θ_j) (interval products), and (b) the centred form
   F(c; θ_c) + H_i ∂_iF(c; Θ) + H_j ∂_jF(c; θ_{c,i}, Θ_j) + Σ_m (u(Θ_i)ᵀ G_m(D) u(Θ_j))·(D_m − c_m)
   (mean value theorem in θ and in d; H = [−π/2^{Lθ+1}, π/2^{Lθ+1}]). The values cos, sin at the grid points
   kπ/2^{14} are enclosures computed with mpmath at 60 digits and rounded outward (`scripts/make_trig_table.py`);
   on a cell of [0,π] cos is monotone and sin is monotone on each half, so cell ranges come from the endpoints.
   Two exact algebraic simplifications are applied before evaluating: (i) a determinant whose four vertices are
   v1, v3, v5 and one apex v_{2k} equals ±d1·q·h_k·sin θ_k (the sign being the parity of the vertex order), hence has
   a known strict sign on int(P) × (0,π)³ and is not evaluated; (ii) a determinant containing a whole flap triangle
   {v_{2k−1}, v_{2k}, v_{2k+1}} equals h_k times the determinant in which the row of v_{2k} is replaced by
   (0, cos θ_k n_k + sin θ_k e_z) (the midpoint part vanishes since m_k lies on the line v_{2k−1}v_{2k+1}); as
   h_k > 0 on int(P), the reduced determinant has the same sign. Both remove spurious degeneracies of the enclosures
   at the faces d_k = 2 and at the flat-central-triangle face d1 = d2 + d3. Finally the components of the unit
   vectors n2, n3 are intersected with [−1, 1].
3. **Status.** A condition is *certified true* on D × cell if the sign pattern of Section 1.2 holds with all five
   enclosures strictly signed, *certified false* if the first two enclosures have the same strict sign or two of
   the last three have strictly opposite signs, *unknown* otherwise. A quadtree on each angle pair refines unknown
   cells down to level Lθ (signs certified on a parent are inherited).
4. **Trace formula.** With bitmaps P12, P23, P13 (inner: certified true; outer: not certified false) of the three
   pairs, the number of 3D cells of (0,π)³ in the corresponding product set is
   Σ_{a,b} P12[a,b] · popcount(P23[b,·] ∧ P13[a,·]) (exact integer).
   outR, outL (outer, for R and L) bound the cells that may contain points of R, L; by Prop. 3.2 the knotted set
   of D × (0,π)³ is covered.
5. **Knottedness of inner cells.** The certified R-cells (resp. L-cells) of D are split into face-connected
   components (union–find on runs). Because of the simplifications (i), (ii), the strict signs are certified for the
   true determinants on (D ∩ int P) × (Θ ∩ (0,π)³) for every certified cell Θ; these sets are convex, and two
   face-adjacent cells give sets whose union is path-connected (the common face meets int P × (0,π)³), so each
   component is a path-connected set of embedded hexagons (Prop. 3.3). For each component, at the centre of D and of one of its cells, the six
   vertices are enclosed by intervals, the xy-projection is checked to be regular with certified signs (all
   segment-intersection orientations, height differences at crossings and orders of crossings along edges
   strictly decided), and the Fox 3-colouring matrix is reduced exactly mod 3. If the count is 9, all cells of
   the component are counted (certR, certL); otherwise none (the program also reports components certified as
   *unknotted*, which would contradict Prop. 3.3 + Calvo; there were none).
6. **Exactness.** Interval operations: IEEE binary64 with the rounding mode set to upward (lower ends computed as
   −(−x ⊕ ...)), compiled with `-frounding-math -ffp-contract=off`, no `-ffast-math`; square roots widened by one ulp;
   division by an interval containing 0, infinities and NaNs give (−∞, ∞). Tested against exact rational
   arithmetic (`scripts/test_interval.py`). Box weights and sums are exact integers over 2^{400} (Python).
7. **Branch and bound.** Starting from 16·8·8 boxes with Lθ = 6, in each round the 2000 leaves with the largest gap
   w(D)(out − cert) are refined: boxes coarser than 2⁻⁶ (in s, t) are split into 8 with Lθ tied to the box size
   (Lθ = level + 2); finer boxes are either re-evaluated with Lθ + 1 or bisected along one d-axis, whichever carries
   the largest share of the remaining uncertainty according to a (non-rigorous) attribution computed by the
   certifier from the widths of the centred-form terms of the undecided cells. Every final leaf, including those
   with zero outer count, is listed in the certificate; `verify_certificate.py` checks that the leaves form a k-d
   tree partition of [0,2] × [0,1]² (recursive bisection check and exact total volume 2), recomputes every leaf
   from scratch, and re-sums the exact bounds.

Then U = Σ_D w(D)(outR + outL) and L = Σ_D w(D)(certR + certL) satisfy L ≤ p₆ ≤ U by (2), (3), Cor. 3.5.

## 6. Results

| | value |
|---|---|
| upper bound U | 8.494446213400·10⁻⁴ (exact: `upper_num`/2⁶⁰⁰ in `certificates/hexknot.summary.json`) |
| lower bound L | 3.067561320913·10⁻⁵ (exact: `lower_num`/2⁶⁰⁰) |
| U/L | 27.69 |
| previous bounds | 0 ≤ p₆ ≤ (14−3π)/192 = 0.02383 (Hake; proof incomplete, §7) — U is 28.05 times smaller |
| leaves / evaluations / rounds | 245 448 / 292 472 / 24 (stopped at the owner's request; resumable) |
| certificate | `certificates/hexknot.leaves.gz` (SHA-256 of the text 799b8710…a7), 2.1 MB |
| anomalies | 0 components certified unknotted; all components with certified cells certified knotted |

The run used the octree phase only (24 of the 24 rounds split boxes coarser than 2⁻⁶ isotropically; 6 θ-refinements
in the last round). Verification performed: tiling (k-d tree + exact volume 2), exact re-summation, recomputation
of 300 random leaves (identical) — `certificates/verify_sample300.log`. **Not performed:** recomputation of all
245 448 leaves (~30 min, one core). An earlier run with an older certifier (no exact-sign rules) reached
U = 2.08·10⁻³, L = 1.1·10⁻⁵ after 75 000 evaluations (not kept as a certificate).

## 7. Numerical validation (not part of the proof)

* **Monte Carlo, 10⁸ samples** (`scripts/mc_validate.py 100000000 2026`, log `certificates/mc_validate.log`):
  d uniform in P (rejection from [0,2]³), θ uniform in [0,2π)³, floating point. Knotted fraction
  1.3639·10⁻⁴ ± 1.17·10⁻⁶ (13 639 knots); the classifier used for "knotted" is the exact Fox 3-colouring count of the projection on a random
  plane (integer arithmetic mod 3 on the Gauss diagram). Checks, all PASS: the count is always 3 or 9 (C1); it agrees
  with Calvo's criterion on every sample (C2); every knotted sample has each T_i pierced exactly once (C3, Lemma 2.3);
  the A/B conditions on the same angle pair never hold together (C4, Lemma 2.1); in (0,π)³ knotted ⇔ R or L (C5),
  with Δ = +1 on R and −1 on L (C5'); knots occur only in (0,π)³ ∪ (π,2π)³ (C6, Lemma 2.4).
* **Box-level consistency** (`scripts/check_boxes_mc.py`, run with the final certifier; PASS on all 16 boxes, log not
  saved — rerun takes ~1 min): for 16 boxes (10 interior, 6 touching degenerate faces) of various sizes, the Monte Carlo frequency
  of R and L inside D × (0,π)³ (4·10⁵ samples each, separate Python implementation of the geometry) lies within
  the certified [inner, outer] interval (±4 s.e.), and every certified component was certified knotted.
* **Hake's lemmas** (`scripts/hake_lemma_check.py`, certified with ball arithmetic, not numerics): knotted
  hexagons with curl +1 such that d1 > d2 > d3 and θ2 > π/2, resp. d1² > d2² + d3² and θ1 < π/2, exist; this
  contradicts Lemmas 3.9/3.12 of arXiv:1810.11882 (the same assertion for J = (±1, 1)), which her proof of
  Theorem 4.1 uses. In the Monte Carlo sample of knots in the chamber C × (0,π)³ (22 090 knots), 0.6 % violate the
  first and 1.4 % the second assertion, all of them in pattern L; Lemma 3.8 (θ_i + θ_j < π) was never violated.

## 8. Assumptions, limitations, what is not proved

* **Used without proof:** F1–F3 (classical knot theory), the measure-preserving property of the action-angle
  map (Cantarella–Shonkwiler, Theorems 1 and 13, resting on Duistermaat–Heckman, Kapovich–Millson,
  Howard–Manon–Millson, Hitchin), correctness of IEEE-754 arithmetic as implemented by the hardware/compiler with
  the stated flags (tested against exact rationals), and of mpmath's cos/sin at 60 digits (used with a 10⁻⁴⁰ margin).
* **Not used:** Calvo's Δ-criterion, Hake's lemmas, any Monte Carlo result.
* **Heuristics** (the refinement order, the attribution of the uncertainty used to choose the split axis, the
  sample cells used for the Fox test) influence only the tightness, never the validity, of the bounds: every leaf
  of the final tiling is re-certifiable from scratch by `verify_certificate.py`.
* **Gap.** The enclosure is limited by first-order blur: a cell is counted as certainly knotted only if it is
  knotted for all d in the box, so the error is proportional to the area of the boundary of the knotted set times the
  box size, in the d-directions (dominant) and in θ. Halving the d-box size costs a factor 8 (or 2 per axis), so
  the gap shrinks slowly with computing time; see README for running with more cores.
