# Result: cospectral fullerene duals and Conjecture 1 of Bille–Buchstaber–Spodarev

Notation (as in the source). C_n: the set of combinatorial isomers of fullerenes with n atoms (n = 20 or even
n ≥ 24, "feasible"); isomers are classes of combinatorially equivalent polytopes, i.e. isomorphism classes of
fullerene graphs (mirror images are the same isomer). T_n: the dual graph of a fullerene (vertices = facets,
adjacent iff the facets share an edge): a plane triangulation with m = n/2 + 2 vertices, twelve of degree 5
(pentagons) and the others of degree 6 (hexagons). T_n^5, T_n^6: the subgraphs of T_n induced by the degree-5,
resp. degree-6 vertices. "Cospectral" (Definition 1 of the source) = same spectrum of the adjacency matrix,
equivalently (Lemma 4 of the source) the same characteristic polynomial. Isomer numbers C_n:k are the ranks of
the lexicographically smallest face spirals (Fowler–Manolopoulos Atlas convention, the convention C_{n,k} of the
source; our implementation reproduces C60:1812 = I_h and the source's pair C44:37/C44:38).

## 1. The statements (verbatim)

### 1.1 arXiv:1910.00846v1 (2 Oct 2019), Section 3, p. 9

> "Conjecture 1. a) For all feasible n ≠ 44 two Cn–isomers are isomorphic iff they are cospectral with respect
> to Tn.
> b) For all feasible n ≥ 54 two Cn–isomers are isomorphic iff they are cospectral with respect to Tn6.
> c) For any feasible n at least one of the spectra σ(Tn), σ(Tn5), σ(Tn6) is unique for all Cn–isomers."

It rests on Table 1 (p. 9), whose T_n row is 1 at n = 44 and 0 for all other 32 ≤ n ≤ 60 ("This table can be
extended to n > 60 with all zeroes in its first and third rows"), computed with "MATLAB functions eig(),
charpoly() and trace() with double precision", for all even 24 ≤ n ≤ 150.

So **Conjecture 1(a) (arXiv)** concerns the **whole dual T_n**, the **adjacency spectrum**, and claims that the
spectrum determines the isomer (except at n = 44).

### 1.2 Published version, J. Math. Chem. 59 (2021) 264–288 (online 18 Nov 2020, CC BY), Section 3

p. 274: "For all even 24 ≤ n ≤ 150, we checked the existence of cospectral pairs of Cn-isomers within Tn, Tn5
and Tn6 graphs by using Lemma 4 and computing all Newton polynomials in Python's int-precision. [...] Considering
the whole dual graph Tn, one finds a lot of non-isomorphic but cospectral graphs. In fact, it seems likely that
infinitely many such pairs exist." Table 1 (p. 274) now has T_n row 1, 1, 1, 3, 2 at n = 42, 44, 50, 54, 60.
p. 275:

> "Conjecture 1
> (a) For all feasible n the graph Tn6 is uniquely determined by its spectrum.
> (b) For all feasible n ≥ 60 at least two Cn–isomers exist having the same spectrum with respect to Tn.
> (c) For all feasible n ≥ 54 two Cn–isomers are isomorphic iff they are cospectral with respect to Tn6.
> (d) For any feasible n at least one of the spectra σ(Tn), σ(Tn5), σ(Tn6) is unique for all Cn–isomers."

So the published Conjecture 1(a) is about **T_n^6**, and the published 1(b) asserts the **opposite** of the
arXiv 1(a): T_n-cospectral non-isomorphic isomers exist for every n ≥ 60.

### 1.3 A later conjecture of the same authors (arXiv:2410.19322v1, Section 5.1, p. 13; J. Chem. Inf. Model. 2025)

> "Conjecture 3. For α, β ∈ (0, 1) and feasible n1, n2, let Tn1, Tn2 be two dual fullerenes. It holds that
> Tn1 ≃ Tn2 ⇔ chα,β(Tn1) = chα,β(Tn2)."   with   chα,β(Tn) := tr(exp(αA + βD))  (Definition 5),

followed by "Open Problem 4. Prove or falsify Conjecture 3." (The JCIM version could not be read: NEEDS.md.)

### 1.4 Fowler's 2005 announcement

P. W. Fowler, CMS Winter Meeting 2005, abstract "Cospectrality and the fullerenes": "As yet, no cospectral
fullerenes are known, but many pairs of fullerenes with cospectral duals have now been found (P. W. Fowler and
M. Duncan, to be published). [...] A construction for (infinite) families of cospectral fullerene duals will be
presented." Fowler's "dual" is the triangulation T_n with its adjacency spectrum: the same notion as Bille et al.
The abstract gives no pair, and we found no publication (references.md). It announces counterexamples to the
arXiv Conjecture 1(a) fourteen years before it was stated, but it is not a refutation in print.

## 2. Results

Throughout, N denotes the largest n of the exhaustive scan: **N = 138** (every even n with 20 ≤ n ≤ 138; the scan was
stopped at the owner's request during n = 140; counts per n in `certificates/counts_per_n.csv`). Matrix types
scanned: A(T), A(T^6), A(T^5) for all n ≤ 138; A + D for n ≤ 122; 2A + D for n ≤ 104.

**R1 (arXiv Conjecture 1(a)) — false; already refuted in print by the authors' own journal version.**
The smallest counterexample is **n = 42**: the isomers **C42:6** and **C42:22** (canonical spirals
(1,2,3,4,7,10,17,18,19,20,21,22) and (1,2,3,5,7,9,16,17,19,20,22,23)) have non-isomorphic duals T42 (23 vertices)
with the same characteristic polynomial

x^23 − 63x^21 − 84x^20 + 1542x^19 + 3876x^18 − 16676x^17 − 67344x^16 + 50826x^15 + 529528x^14 + 408770x^13
− 1738312x^12 − 3290096x^11 + 1137952x^10 + 7772140x^9 + 5023808x^8 − 5192475x^7 − 8179896x^6 − 2199987x^5
+ 1898460x^4 + 1171962x^3 + 1404x^2 − 115992x − 19440.

Non-isomorphism is certified by several invariants (different canonical spirals; different T42^5 and T42^6
spectra; automorphism groups of orders 4 and 2; nauty canonical forms). No two isomers with n < 42 have
cospectral duals (exhaustive, exact). Further counterexamples: n = 50 (C50:4/C50:5), 54 (3 pairs), 60 (2 pairs),
and every n from 60 to N (R2). The authors' own journal version (Table 1, p. 274) lists exactly our counts for
n ≤ 60 (1, 1, 1, 3, 2 pairs at n = 42, 44, 50, 54, 60) and replaced the conjecture; so the arXiv 1(a) was
corrected in print by its authors (Nov 2020). Our exact computation is an independent confirmation.

**R2 (published Conjecture 1(b)) — verified for 60 ≤ n ≤ N, and PROVED for the two infinite residue classes
n ≡ 20 (mod 24), n ≥ 44, and n ≡ 2 (mod 24), n ≥ 50.**
(i) Exhaustive exact scan: for every even n with 60 ≤ n ≤ N there is at least one pair of non-isomorphic
isomers with T_n-cospectral duals (number of classes per n in counts_per_n.csv; e.g. 2, 1, 2, 5, 2, 3 classes for
n = 60, …, 70; 16 for n = 100). (ii) Theorem 2 below (cap rotation on the thinnest zigzag nanotube) gives, for
every odd k ≥ 1, cospectral non-isomorphic pairs with n = 32 + 12k (= 44, 68, 92, …) and n = 38 + 12k
(= 50, 74, 98, …). The residue classes other than 2 and 20 mod 24 remain open beyond N.
These families also give a proof of Fowler's announced "infinite families of cospectral fullerene duals".

**R3 (published Conjecture 1(c)) — verified for 54 ≤ n ≤ N (finite verification only).** No two non-isomorphic
C_n isomers with 54 ≤ n ≤ N have cospectral T_n^6. For n < 54 the only T^6-cospectral classes are at
n = 32, 36, 40, 52, each consisting of two isomers whose T^6 graphs are isomorphic (this agrees with the arXiv
version and with the sentence "For n ∈ {32, 36, 40, 52} a pair of distinct isomers exists with the same
spectrum σ(Tn6), since their graphs Tn6 are isomorphic" of both versions; the journal's Table 1 entries 1, 2, 1, 2
at n = 30, 34, 38, 40 and the sentence "we found cospectral pairs only for n ∈ {30, 32, 34, 36, 38, 40, 52}"
are NOT reproduced: at n = 30, 34, 38 all T^6 spectra are distinct, and n = 40 has one class).

**R4 (published Conjecture 1(a)) — false in its literal sense; true on [24, N] in the restricted sense.**
"Determined by its spectrum" has a standard meaning (no non-isomorphic graph whatsoever is cospectral with it).
In that sense (a) fails first at n = 34: T34^6 of C34:2 is C6 ∪ K1, which is cospectral with the subdivided claw
(spectrum {±2, ±1, ±1, 0}); more examples at n = 36 (C36:15) and n = 38 (C38:4, C38:9, C38:15); none for n ≤ 32
(exhaustive over all graphs with the same number of vertices and edges, nauty geng). In the restricted sense
suggested by the authors' check ("We also checked that these cospectral pairs of Tn6-graphs are isomorphic":
cospectral T^6 graphs of two C_n isomers are isomorphic), it holds for every n ≤ N (finite verification).

**R5 (published Conjecture 1(d), arXiv 1(c)) — verified for n ≤ N (finite verification).** For n ∈ {32, 36, 40, 52}
σ(T_n) separates all isomers; for every other n ≤ N, σ(T_n^6) does. (σ(T_n^5) is non-unique for every n ≥ 40.)

**R6 (Conjecture 3 of arXiv:2410.19322) — false for every α = β ∈ (0, 1).** The duals of **C44:72** and
**C44:80** are non-isomorphic and have the same signless Laplacian A + D characteristic polynomial (exact), hence
ch_{α,α}(T) = tr exp(α(A + D)) coincide for every α. Further A + D-cospectral pairs: one each at n = 84 and n = 92, and none other for
n ≤ 122 (the A + D type was scanned exhaustively up to n = 122 only). For the authors' preferred
parameters (α, β) = (1/2, 1/4), i.e. the matrix 2A + D, no cospectral pair exists for n ≤ 104 (exhaustive); a
coincidence of these transcendental characters without cospectrality was not investigated. Hence: under the
natural (universal) reading "for all α, β ∈ (0, 1)" Conjecture 3 is false; under an existential reading (some
fixed pair (α, β) separates all duals) our pair only excludes the diagonal α = β. Since 2A + D-cospectrality is
implied by "A + cD-cospectral for every c", no two duals with n ≤ 104 have the same bivariate polynomial
det(xI − A − cD) either.

**R7 (side question; Fowler 2005: "no cospectral fullerenes are known").** No two non-isomorphic fullerene graphs
C_n (the carbon skeletons; Hückel adjacency spectrum) with n ≤ 100 are cospectral (same complete-hash method,
K = 24, all key collisions settled exactly; `scripts/scan_fullerene_graphs.py`, finite verification); in particular
the cospectral-dual pairs with n ≤ 100 never have cospectral fullerene graphs.

## 3. Theorem 2 (cap rotation) and its proof

**Lemma 1 (Schur complement).** Let a graph G have vertex set I ⊔ R, adjacency blocks A_I, A_R, and I–R
incidence B (|I| × |R|). Then, as polynomials,
det(xI − A_G) = p_I(x) · det( xI − A_R − W(x)/p_I(x) ),  p_I = det(xI − A_I),  W = Bᵀ adj(xI − A_I) B.
Hence if two graphs G1 = (I1, R, B1), G2 = (I2, R, B2) share R and satisfy p_{I1} = p_{I2} and W_1 = W_2, they are
cospectral, *whatever R is*. (Proof: block LU factorisation, valid for x not an eigenvalue of A_I, hence
identically.)

**The tube and the cap.** In the thinnest zigzag nanotube of fullerene duals ("6-row tube"), a *row* is an
induced 6-cycle of hexagons; consecutive rows satisfy (r+1, j) ~ (r, j−1), (r, j) (indices mod 6). The cap
**I(δ)** (δ ∈ {0, 1}) on top of a row S = (S_0..S_5) consists of a row O_0..O_5 with O_j ~ S_{j−1}, S_j, where
O_j is a pentagon iff j ≡ δ (mod 2), and a triangle of pentagons t_0 t_1 t_2 with t_i ~ O_{p−1}, O_p, O_{p+1},
p = 2i + δ. It has 9 faces, 6 of them pentagons, symmetry C3v; I(1) is I(0) rotated by one step (60°). It is the
cap of the pair C44:37/C44:38 of Bille et al.

**Lemma 2.** For I(0) and I(1) glued to S: p_I(x) = x^9 − 18x^7 − 20x^6 + 72x^5 + 120x^4 − 54x^3 − 180x^2 − 81x
for both, and W_0(x) = W_1(x) (6 × 6 matrices over Z[x]).

*Proof 1 (computation, `scripts/tube_family.py`, claim F2a).* Entries of W have degree ≤ 8; the two matrices are
evaluated exactly (python-flint, over Q) at 12 integer points ≥ 10 and agree, which is an identity test.

*Proof 2 (conceptual).* p_I is the same because I(1) ≅ I(0) as graphs. Index O and S by Z/6 and use the Fourier
basis e_a(j) = ω^{aj}/√6, ω = e^{2πi/6}. The S–O incidence (S_j ~ O_j, O_{j+1}) is circulant with multiplier
β_a = 1 + ω^{±a}, so β_3 = 0. By the Schur complement on the triangle t, the O-block of (xI − A_I)^{-1} is
R = (xI − A_O − Kᵀ(xI − A_t)^{-1}K)^{-1}, where A_O (a 6-cycle) is diagonal in this basis and K (t_i ~ O_{p−1},
O_p, O_{p+1}, p = 2i + δ) sends e_a to κ_a ω^{aδ} times the C3-mode a mod 3 of the triangle, with
κ_a = 1 + 2cos(πa/3) = 3, 2, 0, −1, 0, 2 for a = 0..5. Hence Kᵀ(xI − A_t)^{-1}K, and therefore R, is block diagonal
for the blocks {0, 3}, {1}, {2}, {4}, {5} (the would-be couplings 1–4 and 2–5 carry the factor κ_4 = κ_2 = 0).
Then W = Bᵀ R B has the only possible off-diagonal coupling 0–3 multiplied by β_3 = 0, so W is diagonal in the
Fourier basis, i.e. circulant. Rotating the cap by one step conjugates W by the cyclic shift, which fixes every
circulant: W_1 = W_0. ∎

**Theorem 2.** (A) Let G be any fullerene dual in which a hexagon row S of a 6-row tube carries the cap I(δ)
(all other vertices form the rest R). Replacing I(δ) by I(1 − δ) gives a fullerene dual G' with the same
characteristic polynomial as G. (B) Let G(Z, k, δ) = Z + k rows + I(δ), where Z is one of the caps of the
6-row tube (sizes 7, 8, 9, 10, 12; `certificates/caps_6row_tube_n52.json`). For **Z = I (|Z| = 9) and Z = Z12
(the 12-face cap) and every odd k ≥ 1**, G(Z, k, 0) and G(Z, k, 1) are **non-isomorphic** and cospectral; their
sizes are n = 32 + 12k and n = 38 + 12k respectively. (For Z7, Z8, Z10, and for even k, the two graphs are
isomorphic.)

*Proof.* (A) is Lemmas 1–2 with R = the rest (I has neighbours only in S); G' is a fullerene dual because the
cap is replaced by a congruent disk with the same boundary.
(B) Validity: rows are triangulated annuli, gluing disks along their boundary rows gives a triangulated sphere
with degrees 5/6 and twelve 5's (also re-checked by planarity for k ≤ 12, claim F1). Non-isomorphism: let
φ : G0 → G1 be an isomorphism (G_δ = G(Z, k, δ)).
1. φ maps degree-5 vertices onto degree-5 vertices, hence induces an isomorphism of the pentagon graphs. In G_δ
   the pentagons of I form a component "net" (triangle with three pendant vertices) of T^5; Z12's pentagons form
   K1 ⊔ (5 vertices, 7 edges), so for Z = Z12 φ maps the I-pentagons to the I-pentagons. The vertex set of the cap
   is the net plus the three hexagons with ≥ 2 neighbours in it, so φ(I(0)) = I(1) and φ maps the rest R onto R.
   (For Z = I there are two nets; see step 5.)
2. Inside R, the degree-4 vertices form S, and row r is the set of vertices at distance k − r from S; so
   ψ = φ|_R is an automorphism of R preserving every row and the cap Z. Its restriction to the patch
   H = Z ∪ (row 1) is an automorphism of H, and ψ is determined by it (row r+1 vertex j is the unique common
   neighbour, in row r+1, of row-r vertices j−1, j). Conversely every automorphism of H extends: a rotation
   j ↦ j + t of row 1 acts as the same rotation on every row; a reflection j ↦ c − j of row 1 acts on row r as
   j ↦ c + r − 1 − j.
3. A dihedral map g of S = row k sends the cap I(δ) to I(δ') with δ' − δ ≡ t (rotation by t) or with δ' ≡ δ iff
   c_k is odd (reflection j ↦ c_k − j), because it sends O_j (attached to S_{j−1}S_j) to O_{j+t}, resp. O_{c_k+1−j}.
4. Hence G0 ≅ G1 by an end-preserving map iff Aut(H) contains a rotation with t odd or a reflection with
   c + k − 1 even. The automorphisms of the five cap patches and their actions on the row are computed in
   `tube_family.py` (exact, VF2): Z7: all of D6; Z8: rotations {0, 3}, reflections c ∈ {0, 3}; I: rotations
   {0, 2, 4}, reflections c ∈ {1, 3, 5}; Z10: rotations {0, 3}, c ∈ {2, 5}; Z12: identity and c = 5. So for Z12
   (and for I, as far as end-preserving maps go) the two graphs are non-isomorphic exactly for odd k; for Z7,
   Z8, Z10 never.
5. For Z = I an isomorphism might exchange the two caps. Such a φ maps the rows of G0 to the rows of G1 in
   reverse order; the tube maps reversing the rows are (r, j) ↦ (k + 1 − r, j + b − r) and (r, j) ↦
   (k + 1 − r, a − j) (the only maps compatible with (r+1, j) ~ (r, j−1), (r, j)). Writing β for the parity of
   the bottom cap and δ for the top one, the second kind never changes the relative parity, and the first kind
   maps (β, 0) onto (β, 1) iff k is even. So for odd k no cap-exchanging isomorphism exists either. ∎

The prediction of steps 2–5 (isomorphic iff "parity-swapping symmetry exists") is cross-checked against nauty
canonical forms for all five caps and all k ≤ 12 (claim F3; n up to 182), together with exact characteristic
polynomials (F2b). Consequences: n = 44, 68, 92, 116, 140, 164, … (Z = I) and n = 50, 74, 98, 122, 146, …
(Z = Z12). Of the 6-row-tube caps, sizes 7, 8, 9, 10, 12 (+ multiples of 6 by extra rows) are the only ones
found in all fullerenes up to C64 (no cap of size ≡ 5 mod 6); unless larger caps of this tube exist, the
mechanism does not reach the residue class n ≡ 0 (mod 12); it also does not explain most pairs (e.g. not C42:6/22; `scripts/classify_pairs*.py`).

## 4. Method of the exhaustive scan

1. **Generation.** buckygen 1.1 (Brinkmann–Goedgebeur–McKay) writes every isomer of C_n once, as the dual
   triangulation T_n in planar code (no `-d`), in a deterministic order. The isomer count is compared with OEIS
   A007894 for every n (all equal; consistency check of the generator, not a proof of completeness).
2. **Complete hashing.** For each isomer and each matrix M ∈ {A(T), A(T^6), A(T^5), 2A+D, A+D} the C program
   `dualspec.c` computes the closed-walk numbers tr(M^k), k = 1..K, modulo 2^32 (n ≤ 104, K = 20) or 2^16
   (n ≥ 106, K = 16), and hashes them (128 bits). Reduction mod 2^w is a ring homomorphism, so two cospectral
   matrices always have the same key: **no cospectral pair can be missed** (no false negatives). All key
   collisions are printed.
3. **Exact settlement.** Every graph in a key collision is extracted (`pcselect.c`) and the integer
   characteristic polynomials are computed with python-flint (`fmpz_mat.charpoly`); a key group is split into
   exact spectra, and only groups of ≥ 2 isomers with identical polynomials count. Hash false positives
   ("mixed" groups in scan_log.txt) are thus eliminated. Non-isomorphism inside a class is guaranteed by
   buckygen (one graph per isomorphism class) and certified again by nauty canonical forms; for T^6 classes the
   T^6 graphs themselves are compared with nauty.
4. **Independent re-check** (`verify_pairs.py`, from the graph6 strings alone): validity of every dual, equality
   of exact polynomials, non-isomorphism through canonical face spirals (an isomorphism invariant) and networkx
   VF2, recomputation of the recorded spirals, consistency of the counts table with the class files and OEIS, the
   C34:2 example, and (with `--exhaustive-small`) the minimality n = 42 by an exact exhaustive recount.

Cost: the scan to n = 138 took about 4 hours of one heavily shared core (timings per n in scan_log.txt); C138 has
6 344 698 isomers, C150 would have 14 059 174.

## 5. What is proved, certified, or only numerical; assumptions

- Exact (integer) arithmetic everywhere for the claims; floating point is used nowhere in the certificates.
- **Assumptions:** (i) buckygen generates every isomer exactly once (published, widely used; counts agree with
  A007894); the "no pair for n < 42" and all counts per n depend on it. (ii) Correctness of FLINT's
  characteristic polynomial and of nauty's canonical labelling (the latter is re-checked by independent
  invariants for every certified pair). (iii) The Atlas numbering assumes the canonical-spiral convention
  (lexicographically smallest pentagon vector over all starting faces and both orientations); our code asserts
  that every isomer for n ≤ 70 has a spiral.
- R2(ii)/Theorem 2 is a proof (Lemma 1 is elementary linear algebra; Lemma 2 and the automorphism groups of
  five fixed patches are finite exact computations; the rest is the argument above). The completeness of the
  cap list of the 6-row tube is only used for the negative remark (no residue ≡ 0 mod 12), not for the theorem.
- R1, R4 (literal), R6 are refutations by explicit exact counterexamples.
- R3, R4 (restricted), R5 and R2(i) are **finite verifications** (n ≤ N), not proofs.
- Discrepancies with the journal version's Table 1 (T^5 counts at n = 48, 54, 56, 58 differ by one: ours are
  42, 112, 148, 147 = the arXiv values; T^6 row at n = 30, 34, 38, 40) are reported, not explained.

## 6. State at the stop (24 Oct 2026) and what remains

- Scan complete for 20 ≤ n ≤ 138 (A(T), A(T^6), A(T^5)); 890 A(T)-cospectral classes in total; `verify_pairs.py`
  ALL PASS on these data (`certificates/verify_pairs_log.txt`).
- Published Conjecture 1: (a) literal — false (C34:2); (a) restricted, (c), (d) — no counterexample for n ≤ 138
  (finite verification); (b) — holds for every 60 ≤ n ≤ 138 and, by Theorem 2, for all n ≡ 2, 20 (mod 24); open
  for the other residue classes beyond 138.
- Remaining: finish n = 140..150 (command in README), recompute Atlas numbers up to n = 70
  (`make_certificates.py 70`, ~10 min), a mechanism covering the other residue classes for 1(b), and the JCIM
  wording of Conjecture 3 (NEEDS.md).
