# Conjecture 1 of Bille–Buchstaber–Spodarev: state after the all-n work

Solver: Claude Opus 5.5 (effort not recorded). Verification: pending (nothing here has been checked by a second agent).

| Part | State |
|---|---|
| (a) T_n^6 determined by its spectrum | **false** in the literal sense (n = 34: C6 ∪ K1 and the subdivided claw, re-checked, R1) |
| (b) cospectral T_n-isomers for every feasible n ≥ 60 | **PROVED for every even n ≥ 60** (Theorems A, B below); no use of the scan |
| (c) σ(T_n^6) separates the isomers, n ≥ 54 | **open**; no counterexample for n ≤ 170 (section 4) |
| (d) one of σ(T), σ(T^5), σ(T^6) separates all isomers | **open**; for n ≥ 60 it is equivalent to (c) at the same n (section 4) |

So the resolution is **not complete**: (c) and (d) are neither proved nor refuted. No manuscript was written.

## 1. Theorem A (colour twins): the mechanism

Of the 890 cospectral classes with n ≤ 138, 882 come from a single identity, found by LLL-reducing the fixed space of the
orthogonal matrix Q with Q A(G2) = A(G1) Q: its short vectors are the triangles of a clique partition of T.

Let H be a simple plane triangulation with all degrees in {4, 6} whose degree-4 vertices are pairwise non-adjacent (there
are exactly six by Euler). All degrees being even, the faces have a proper 2-colouring, black and white; every edge
lies on one black and one white face, so both colour classes have m faces, where |V(H)| = m + 2. For a colour c let T_c
be the graph on the faces of colour c, two being adjacent iff they share a vertex of H, and N_c the
(faces of colour c) × V(H) incidence matrix.

**Theorem A.** T_black and T_white are duals of fullerenes with n = 2(m − 2) atoms and have the same characteristic
polynomial.

*Proof.* (i) N_c^T N_c = A(H) + D(H)/2 for both colours: a vertex v lies on d(v)/2 faces of colour c, an edge on exactly
one, and two vertices in a common face are adjacent. So N_black^T N_black = N_white^T N_white.
(ii) N_c N_c^T = 3I + A(T_c): two distinct faces of one colour share at most one vertex, since two common vertices
u, v would put the edge uv (H is simple) on two faces of the same colour.
(iii) N_c has m rows and m + 2 columns, so x^2 det(xI − N_c N_c^T) = det(xI − N_c^T N_c), the same polynomial for both
colours. Hence det(xI − 3I − A(T_black)) = det(xI − 3I − A(T_white)).
(iv) T_c is a fullerene dual. The plane dual Γ of H is cubic, bipartite (colour classes = black, white faces), with faces of
length 4 and 6. Inside each hexagon of Γ join its three vertices of colour c by a triangle, inside each quadrilateral join its
two vertices of colour c by a diagonal. This draws T_c in the plane. Its faces are these triangles and, around each
vertex of Γ of the other colour, the triangle on its three neighbours; so T_c is a triangulation, simple by (ii). A face of
colour c has at most one degree-4 corner, so its degree in T_c is 6 or 5 (the sum over its corners v of d(v)/2 − 1). A simple
triangulation is 3-connected, Euler gives twelve vertices of degree 5, and the dual is a fullerene. ∎

Same proof: N_black W N_black^T and N_white W N_white^T are cospectral for every diagonal W on V(H).
Conversely T = T_c(H) iff E(T) is partitioned into triangular faces and a perfect matching of the pentagons by edges of T
(`scripts/twin.py`); in the pairs, each degree-4 vertex of H is one adjacent pentagon pair of T. A + 3I = N N^T then
gives λ_min(T) ≥ −3.

## 2. Theorem B: one helical family reaches every even n ≥ 60

For L ≥ 15 let H_L have vertices 0, 1, …, L and p0..p4, q1, q2, r0..r3, s1, s2 and these faces (`scripts/family.py`):
- tube: B_t = {t, t+1, t+8} and W_t = {t, t+7, t+8}, 0 ≤ t ≤ L − 8 (a triangular lattice cylinder of chirality (7,1));
- cap 1 (20 triangles, bounded by the 8-cycle 0, 1, …, 7): p0 0 1, p0 1 p4, p0 p4 p1, p0 p1 6, p0 6 7, p0 7 0, p1 p4 q2,
  p1 q2 p2, p1 p2 5, p1 5 6, p2 q2 p3, p2 p3 q1, p2 q1 4, p2 4 5, p3 3 q1, p3 q2 p4, p3 p4 2, p3 2 3, p4 1 2, q1 3 4;
- cap 2 (18 triangles, bounded by the 8-cycle e−7, …, e, e = L; "−k" means e − k): r0 −7 −6, r0 −6 r1, r0 r1 r3, r0 r3 −1,
  r0 −1 e, r0 e −7, r1 r3 r2, r1 r2 s1, r1 s1 −5, r1 −5 −6, s1 −5 −4, s1 −4 r2, r2 −4 −3, r2 −3 s2, r2 s2 r3, s2 −3 −2,
  s2 −2 r3, r3 −2 −1.

**Theorem B.** For every L ≥ 20 the graphs T_b(L), T_w(L) (black = the colour of all B_t) are non-isomorphic cospectral duals of
fullerenes with n = 2L + 20. Hence Conjecture 1(b) holds for every even n ≥ 60.

*Proof.* Every edge lies on two faces and every vertex link is a cycle: for 8 ≤ t ≤ L − 8 the link is t+1, t+8, t+7, t−1, t−8, t−7,
and the other links do not depend on L (cap 2 is cap 1's data shifted with L; finite check, C1). With V = L + 14, F = 2L + 24,
E = 3L + 36 we get χ = 2, so H_L is a triangulated sphere, simple for L ≥ 15. Its degree-4 vertices are 0, q1, q2, L, s1, s2,
pairwise non-adjacent, the others have degree 6. Theorem A applies: m = L + 12.

Non-isomorphism. Let h_j be the number of hexagons of T with exactly j pentagon neighbours, an isomorphism invariant.
Pentagons are the faces containing a degree-4 vertex; those at cap 1 only involve helix vertices ≤ 8, those at cap 2 only
vertices ≥ L − 8. A face of the tube spans at most 8 consecutive labels, so for L ≥ 25 no hexagon touches pentagons of both
caps. Translation t ↦ t + 1 keeps the colours of tube faces, so the two contributions do not depend on L. The computation
at L = 25 (C6) gives for all L ≥ 25: (h_1, h_2, h_3, h_4) = (11, 6, 5, 1) for T_b and (8, 9, 4, 1) for T_w. For
L = 20, …, 24 the vectors (h_0, …, h_4) are, black / white: (0,6,7,6,1) / (1,3,10,5,1), (0,7,8,5,1) / (1,4,11,4,1),
(1,8,6,6,1) / (2,5,9,5,1), (2,8,6,6,1) / (3,5,9,5,1), (2,9,7,5,1) / (3,6,10,4,1). ∎

Certificates: `check_family.py` (exact, L = 20..130, n up to 280): C1–C7 all PASS, with integer identities N_b^T N_b =
N_w^T N_w and N N^T = 3I + A, integer characteristic polynomials, nauty canonical forms. Only L = 20..25 is needed by the proof.
`certificates/family_pairs.json`: H_L, T_b, T_w, polynomial for n = 60..160.

## 3. What the scan now looks like (n ≤ 138)

- `check_classification.py`: **all 890 classes are explained**: 882 colour twins, and the 8 cap rotations of `../result.md`
  (n = 44, 50, 68, 74, 92, 98, 116, 122). K1, K2 PASS.
- `check_converse.py` (all isomers, n ≤ 100): the twin pairs are exactly the classes of the scan apart from those rotations.
  X1, X2 PASS. This also cross-checks the completeness of the old scan.
- The family H_L is the (7,1) tube seen in the data; n = 42 (L = 11) and n = 54 (L = 17) are members too.
- The earlier route is superseded: in the twin pairs sampled, no side of a cut has an invariant Schur complement, and a
  potential on one side breaks cospectrality. The identity is global, not a gadget. Tubes (p,q) of H need p ≡ q (mod 3).

Novelty (three web searches, not a full check): no publication of Theorem A or B was found. Fowler's unpublished 2005
construction with Duncan (`../references.md`) may or may not be the same; nothing about it is available.

## 4. (c) and (d): what is reduced, computed, excluded

**Reduction.** For n ≥ 60, Theorem B says σ(T) is not injective on isomers, and T^5(T_b(L)) ≅ T^5(T_w(L)) (R2, local for
L ≥ 25) says the same for σ(T^5). So (d), read as "one of the three spectra separates all isomers of C_n", is for n ≥ 60
**equivalent to (c) at n**: both say that σ(T^6) is injective on the isomers. Read pairwise (any two isomers differ in
one of the three spectra), (d) is weaker, and is implied by (c).

**Computed.** σ(T^6) is injective on the isomers for every even n with 54 ≤ n ≤ 170, and (d) holds for n ≤ 170. New exhaustive
scan: 353 136 799 isomers for 140 ≤ n ≤ 170, complete hashing, every collision settled by exact T^6 polynomials. The range
n ≤ 138 was redone the same way; it gives again exact classes at n = 32, 36, 40, 52 only (S1, S2 PASS, `certificates/t6_scan.json`).

**Excluded.**
- Twins and cap rotations give no T^6 pair: K3 (exact, all 890 classes). In the family T_b^6, T_w^6 differ in tr A^5 by −20 (R3).
- How close twins get is exact. Let B_hh be the (black hexagons) × (white hexagons) edge-adjacency matrix of H, a square
  matrix. Then B_hh B_hh^T = 3I − D_b + A(T_b^6) and B_hh^T B_hh = 3I − D_w + A(T_w^6), where D_c(f) counts the pentagon pairs of
  which the hexagon f is a common neighbour. So **A(T_b^6) − D_b and A(T_w^6) − D_w are cospectral** (R4): T^6 fails by a
  potential on the 12 flanking hexagons. The potential cannot be removed: Γ is connected, so any deleted vertex set
  leaves boundary rows of smaller norm.
- Pentagon rings disconnecting T^6: bounded size only; the classes n = 32, 36, 40, 52 have isomorphic T^6.
- Theorem A does not transfer (sketch, not a full proof). Two decompositions of one graph into triples would make a
  closed 2-coloured pseudo-surface with χ ≥ (a_2 + 2a_1)/3, where a_k counts lines with k hexagons and a_1 + a_2 = 30 − e, e the
  number of pentagon adjacencies. For large n one has e ≤ 20, so χ > 2, and this forces cut cliques in T^6, which do not exist
  for large n. Holes are the obstruction: at a boundary edge only one colour is present.
- Schur gadgets for T^6, i.e. two fillings of a ring with equal W^6 (sketch, not a proof). Let δ be the least hexagon-distance
  from the ring at which the fillings differ. Then W_1 − W_2 has leading term Σ_X u_s u_s^T − Σ_Y u_s u_s^T, where X, Y are the
  sites at distance δ that are hexagons in one filling only and u_s counts shortest paths from the ring. The extreme entry of a
  site with maximal shadow on the ring can only cancel against a site with the same shadow. So no such gadget exists
  while shadows are proper distinct arcs, which is the case for caps (depth small against the ring). For T the cap I
  escapes this because pentagons stay in the walks.
- Tube twists: the hash collisions include such pairs, e.g. the cap-I rotation at n = 140. Their T^6 polynomials first differ
  at c_17 to c_24 (walks must run along the tube); none is exact.
- High symmetry: the only sporadic coincidences met are the three A + D pairs (n = 44, 84, 92, walk rank 4–8). No analogue
  for T^6 occurs up to n = 170; beyond that nothing was searched.

**Missing.** A proof that σ(T^6) is injective on isomers for all n ≥ 54, or one pair refuting it. I know no method for the
first. The evidence favours (c) being true: twins share the large-scale metric of H and still differ at small k in T^6.

## 5. Files (`all-n/scripts`, run under WSL with the `sci` environment)

| File | Content |
|---|---|
| `check_family.py [LMAX]` | certificate for Theorem B, C1–C7 |
| `check_classification.py` | K1–K3, all 890 classes |
| `check_converse.py [NMAX]` | X1, X2, twins = scan classes |
| `check_reduction.py [LMAX]` | R1 literal (a), R2 T^5, R3 tr A^5, R4 T^6 with potential |
| `make_family_certificate.py`, `certificates/family_pairs.json` | explicit pairs |
| `settle_t6.py`, `t6_summary.py`, `certificates/t6_scan.json` | exact settlement of the T^6 hash collisions, S1, S2 |
| `fdc.py`, `twin.py`, `family.py` | library |
| `explore_Q2.py`, `explore_near.py` | exploration only: how the mechanism was found; anatomy of near-cospectral T^6 pairs |
| `certificates/check_*_log.txt` | stored outputs of the check scripts |

T^6 scan: `buckygen -q Nd | dualspec16 2 16` (`../scripts/dualspec.c`, complete hashing), collisions settled exactly.
