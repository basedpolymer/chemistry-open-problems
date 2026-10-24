# The Kirchhoff Šoltés problem: non-existence theorems and exhaustive searches

Status: **partial result** (not independently verified). No new Kirchhoff Šoltés graph was found; C₅
remains the only known one. We prove non-existence for several natural classes and extend the exhaustive
searches. Labels used below: **[proof]** = complete mathematical proof given here; **[proof + exact
computation]** = proof whose finite part is an exact (integer / rational / modular) computation by the
listed script; **[numerical]** = floating-point exploration only, not a proof.

## 0. Problem, conventions

Gottwald, Majstorović Ergotić, Došlić, *Šoltés problem for the Kirchhoff index of a graph*, Appl. Math.
Comput. 510 (2026) 129694, doi:10.1016/j.amc.2025.129694, arXiv:2609.13751 (read in full, v1), p. 5:

> **Problem 3.** Find all graphs G so that the equality Kf(G) = Kf(G − v) holds for all v ∈ V(G).

p. 9: "Our numerical studies failed to identify any additional Kirchhoff Šoltés graphs. We investigated
connected graphs with up to 9 vertices. [...] Additionally, we examined all vertex-transitive graphs with
up to 23 vertices [...]. None of these graphs were Kirchhoff Šoltés." p. 21: "the existence of other
solutions has not been ruled out"; abstract: "the original problem of determining whether there are
additional solutions beyond C₅".

Kf(G) = Σ_{u<w} r_G(u,w) (effective resistances, unit resistors) = n Σ_{i=1}^{n−1} 1/μ_i, where
μ_1 ≥ … ≥ μ_{n−1} > 0 are the non-zero Laplacian eigenvalues; a(G) = μ_{n−1} is the algebraic
connectivity; Rt_G(v) = Σ_u r_G(u,v); Γ = L(G)^+ (Moore–Penrose inverse); τ(G) = number of spanning trees.

**Convention (as in the source, p. 9: "x is also not a good vertex since it is a cut vertex").** A vertex v
is *good* if G − v is connected and Kf(G − v) = Kf(G). A *Kirchhoff Šoltés graph* (KS graph) is a connected
graph with n ≥ 3 vertices all of which are good. Hence every KS graph is 2-connected. We write
Δ_v(G) = Kf(G) − Kf(G − v).

## 1. Summary of the claims

1. **Thread lemma** (Lemma 2, Cor. 3) [proof]. In a 2-connected graph that is not a cycle, along every
   *thread* (maximal path whose internal vertices have degree 2) the map i ↦ Kf(G − v_i) is a strictly
   convex quadratic polynomial; hence at most two internal vertices of a thread are good. A KS graph
   other than C₅ has all threads with at most two internal vertices.
2. **Bounded cyclomatic number** (Thm. 4) [proof + exact computation]. A KS graph ≠ C₅ with cyclomatic
   number c = m − n + 1 has n ≤ 8(c − 1); there is none with c ≤ 5 (c = 2: all bicyclic graphs, i.e. theta
   graphs; c = 3, 4, 5: 1 038, 76 603 and N₅ candidate graphs, all rejected; see §11 for N₅).
3. **Expander bound** (Thm. 6) [proof]. If G is KS then a(G)·(n − 1)·Σ_v 1/d_v ≤ 2 Kf(G); in particular
   a(G)² ≤ 2n / Σ_v d_v^{−1} ≤ 2·(average degree). For vertex-transitive G violating it, *no* vertex is good.
   Corollary 7: no KS graph has minimum degree δ > (n − 2 + √(2(n − 1)))/2; joins with both sides larger
   than √(2(n − 1)) are not KS.
4. **Dense bound** (Thm. 8) [proof]. If a(G) > 1 then Δ_v(G) ≥ 1 − Σ_i (n − μ_i)/(μ_i(μ_i − 1)) for every v
   (sharp: equality for the source's graphs G_s with half of the vertices good).
5. **Local deletion formula** (Prop. 9) [proof] expressing Kf(G − v) through Γ and Γ² on the closed
   neighbourhood of v, and its distance-regular specialisation (Prop. 10) [proof].
6. **Complete multipartite graphs** (Thm. 11) [proof]: Δ_v = φ(n_j) − S with φ strictly increasing; the good
   vertices form the union of all parts of one size, and no complete multipartite graph is KS.
7. **Strongly regular graphs** (Thm. 12) [proof + exact computation]: C₅ is the only KS strongly regular
   graph. For triangle-free SRGs: Δ_v = N(k,μ)/(μ(k+μ−2)(k²+kμ−k+μ)) > 0 except C₅; for Moore graphs of
   diameter 2: Δ_v = k(k−2)(k+1)/(k²+1) (includes the hypothetical 57-regular Moore graph: 18183/325).
8. **Distance-regular families** (Prop. 13) [exact computation, finite ranges]: hypercubes Q_d (d ≤ 120), odd
   graphs O_{k} (k ≤ 81), incidence graphs of generalized 3-, 4-, 6-gons of order (q,q) (q ≤ 150, existence
   not needed), Hamming H(d,q) (d, q ≤ 30), Johnson J(m,k) (m ≤ 40), the 13 cubic DRGs, cycles C_n (n ≤ 200):
   Δ_v ≠ 0 except for C₅.
9. **Exhaustive searches** [exact computation]: no KS graph on 6 ≤ n ≤ 11 vertices (all 2-connected graphs,
   nauty geng −C, 910 914 360 graphs for n ≤ 11; the published range was n ≤ 9); no vertex-transitive
   KS graph on 6 ≤ n ≤ 47 vertices (full Holt–Royle census; published range n ≤ 23). Counts in §11.

Consequently a KS graph other than C₅, if it exists, has ≥ 12 vertices, is 2-connected, is not
complete multipartite and not strongly regular, has cyclomatic number ≥ 6, threads of length ≤ 3, poor
expansion (a(G)² ≤ 2·HM(degrees)), and if vertex-transitive has ≥ 48 vertices.

## 2. Tools

**Lemma 1 (grounded formula).** Let H be connected on N vertices, g ∈ V(H), A = L(H) with the row and
column of g deleted, M = A^{−1}. Then det A = τ(H), r(u,g) = M_uu, r(u,w) = M_uu + M_ww − 2M_uw (u, w ≠ g),
and

  Kf(H) = N·tr(M) − 1ᵀM1,  so α(H) := τ(H)·Kf(H) ∈ ℤ.

*Proof.* The first identity is the matrix-tree theorem. M is the Green function of H grounded at g, so
the potential of a unit current from u to g (potential 0 at g) is M e_u, giving r(u,g) = M_uu and
r(u,w) = (e_u − e_w)ᵀM(e_u − e_w). Summing over pairs:
Σ_{u<w; u,w≠g} (M_uu + M_ww − 2M_uw) + Σ_{u≠g} M_uu = (N − 2)tr M − (1ᵀM1 − tr M) + tr M. Since
M = adj(A)/det A with adj(A) integral, τ(H)Kf(H) ∈ ℤ. ∎

Further standard facts used: Rayleigh monotonicity (deleting edges does not decrease any effective
resistance); r_H(x,y) ≤ length of any x–y path; Kron reduction (eliminating a vertex v of degree d by
Schur complement replaces the star at v by a complete graph on N(v) with conductances 1/d, without
changing resistances between the other vertices); the Woodbury identity.

## 3. Threads

**Lemma 2 (thread formula).** Let G be connected, P = x v_1 v_2 … v_{ℓ−1} y (ℓ ≥ 2 edges) a path whose
internal vertices have degree 2 in G, x ≠ y, and suppose R = G − {v_1, …, v_{ℓ−1}} is connected, with
n_R vertices. For 1 ≤ i ≤ ℓ − 1 put a = i − 1, b = ℓ − 1 − i. Then

  Kf(G − v_i) = Kf(R) + n_R[a(a+1)/2 + b(b+1)/2] + a·Rt_R(x) + b·Rt_R(y)
        + (a³ − a)/6 + (b³ − b)/6 + ab(a + b + 2)/2 + ab·r_R(x,y).

For fixed s = a + b = ℓ − 2, this is a polynomial of degree ≤ 2 in a with coefficient of a² equal to
n_R − 1 − r_R(x,y).

*Proof.* G − v_i is R with a pendant path p_1 … p_a at x (p_j at distance j from x) and a pendant path
q_1 … q_b at y. Pendant paths are joined by bridges, so r(p_j, u) = j + r_R(x,u) (u ∈ R),
r(p_j,p_k) = |j − k|, r(p_j,q_k) = j + k + r_R(x,y), and symmetrically for q. Summing: pairs in R give
Kf(R); pairs p–R give n_R a(a+1)/2 + a Rt_R(x); pairs p–p give Σ_{j<k}(k − j) = (a³ − a)/6; similarly
for q; pairs p–q give b a(a+1)/2 + a b(b+1)/2 + ab r_R(x,y) = ab(a+b+2)/2 + ab r_R(x,y). For the
leading coefficient substitute b = s − a: n_R[a(a+1) + b(b+1)]/2 contributes n_R a²;
(a³ + b³)/6 = (s³ − 3s²a + 3sa²)/6 contributes (s/2)a²; ab = sa − a² contributes −((s+2)/2 + r_R(x,y))a².
Total: n_R + s/2 − (s+2)/2 − r_R(x,y) = n_R − 1 − r_R(x,y). ∎

**Lemma 2′.** For a connected graph R and x ≠ y, r_R(x,y) ≤ n_R − 1, with equality iff R is a path with
end vertices x and y.

*Proof.* r_R(x,y) ≤ d_R(x,y) ≤ n_R − 1 (keep only a shortest path; Rayleigh). If equality holds, a
shortest x–y path P′ is Hamiltonian. If R had an edge e = p_ip_j ∉ P′ (j ≥ i + 2), then R ⊇ P′ + e, and in
P′ + e the segment p_i…p_j (resistance j − i) is in parallel with e, so r_{P′+e}(x,y) = n_R − 1 −
(j − i) + (j − i)/(j − i + 1) < n_R − 1, and r_R ≤ r_{P′+e} by Rayleigh: contradiction. ∎

**Corollary 3.** Let G be 2-connected, not a cycle, and P a thread (maximal path whose internal vertices
have degree 2 in G) with t ≥ 1 internal vertices. Then (i) at most two internal vertices of P are good;
(ii) if t = 2 and both are good, then Rt_R(x) = Rt_R(y); (iii) consequently every KS graph is either C₅ or
has all threads with at most two internal vertices.

*Proof.* The end vertices x, y of P have degree ≥ 3, and x ≠ y (otherwise x would be a cut vertex). R is
connected: every component of R contains x or y (a component avoiding both would be a component of G),
and if x and y were in different components, the component of x contains a neighbour of x outside P
(deg x ≥ 3), so x would be a cut vertex of G. R is not an x–y path (else G would be a cycle), so by Lemmas
2 and 2′ the function a ↦ Kf(G − v_{a+1}), a = 0, …, t − 1, is the restriction of a real quadratic with
positive leading coefficient; it takes the value Kf(G) at most twice: (i). For t = 2 (s = 1), the formula
gives Kf(G − v_2) − Kf(G − v_1) = Rt_R(x) − Rt_R(y) (all other terms coincide for (a,b) ∈ {(0,1),(1,0)}):
(ii). For (iii): a KS graph is 2-connected; if it is a cycle C_n then (n³ − n)/12 = Kf(C_n) =
Kf(P_{n−1}) = (n − 2)(n − 1)n/6 forces n = 5 (source, Prop. 6); otherwise apply (i). ∎

Check: `scripts/check_thread_lemma.py` (exact, 175 deletions, random R; also the second differences
2(n_R − 1 − r_R(x,y)) and criterion (ii)).

## 4. Bounded cyclomatic number

**Theorem 4.** Let G be a KS graph other than C₅ and c = m − n + 1 its cyclomatic number. Then c ≥ 2,
n ≤ 8(c − 1), and G is obtained from a loopless multigraph K without cut vertex, of minimum degree ≥ 3,
with m_K − n_K = c − 1 (so n_K ≤ 2(c − 1), m_K ≤ 3(c − 1)), by subdividing every edge of K zero, one or two
times. **For c ≤ 5 there is no such graph**; in particular no bicyclic (theta) graph and no tricyclic graph
is KS, and every KS graph other than C₅ has m ≥ n + 5.

*Proof.* c ≥ 2 since a 2-connected graph with c = 1 is a cycle. Suppressing the degree-2 vertices of
a 2-connected graph that is not a cycle gives its kernel K: vertices = vertices of degree ≥ 3 of G, edges
= threads; K is loopless (a thread from x to x would make x a cut vertex) and has no cut vertex (a cut
vertex of K is one of G), and c(K) = c(G). From 2m_K = Σ deg ≥ 3n_K and m_K = n_K + c − 1 we get
n_K ≤ 2(c − 1), m_K ≤ 3(c − 1). By Corollary 3 every thread has t_e ∈ {0,1,2} internal vertices, so
n = n_K + Σ_e t_e ≤ 2(c − 1) + 2·3(c − 1). (G simple forces at most one edge of each parallel class of K
to have t_e = 0.) For the computation, `scripts/kernel_enum.py c` enumerates, for every n_K ≤ 2(c − 1),
all 2-connected simple graphs on n_K vertices (nauty `geng -C`; for n_K = 2 the edge K₂) as underlying
graphs, all multigraphs with exactly m_K = n_K + c − 1 edges on them (nauty `multig -e`), keeps those of
minimum degree ≥ 3, and outputs every admissible subdivision (multisets of t-values per parallel class).
A loopless multigraph on ≥ 3 vertices has no cut vertex iff its underlying simple graph is 2-connected, so
the list of kernels is complete. Every generated graph is tested by `kscheck` (exact criterion of §11).
Results (`certificates/search/kernel_c*.log`): c = 2: 1 kernel, 7 graphs (the theta graphs θ(a,b,c) with
lengths in {1,2,3}, at most one equal to 1; exact values in `certificates/small_cases.txt`); c = 3: 4
kernels, 1 038 graphs; c = 4: 17 kernels, 76 603 graphs; c = 5: see §11. None is KS. ∎

*Remarks.* (1) A 2-connected cactus is a cycle, so C₅ is the only KS cactus (in particular the only KS
unicyclic graph, as in the source). (2) Good vertices of degree 2 do exist in non-cycle graphs (e.g. two
2-connected graphs on 8 vertices, graph6 ``G?qa`o`` and ``G?qa`k``), so the bound "at most two per thread" cannot
be replaced by "none" in general.

## 5. An expander bound

For a vertex v of G and u ≠ w, let φ = Γ(e_u − e_w) be the potential of a unit current from u to w in G
and E_v(u,w) = Σ_{x∼v}(φ_x − φ_v)² the power dissipated in the edges at v.

**Lemma 5.** Let G be connected, v a vertex with H = G − v connected, a_v = a(H). For u, w ∈ V(H):

  E_v(u,w) ≤ r_H(u,w) − r_G(u,w) ≤ (1 + 1/a_v)·E_v(u,w).

*Proof.* Work in ℝ^{V(H)} and on the subspace 1^⊥. Kron reduction of v gives L_{H′} = L_H + S with
S = P_N − (1/d)1_N1_Nᵀ (N = N(v), d = |N|), and r_G(u,w) = bᵀL_{H′}^+b, b = e_u − e_w. S is the orthogonal
projector onto Y_N = {y supported on N, Σy = 0} ⊂ 1^⊥; write S = QQᵀ with Q orthonormal. On 1^⊥ put
X = L_{H′}^{−1}; then L_H = X^{−1} − QQᵀ and by Woodbury
r_H(u,w) − r_G(u,w) = cᵀ(I − K)^{−1}c with c = QᵀXb, K = QᵀXQ. On 1^⊥, L_H ⪰ a_v I, hence
L_H + S ⪰ a_v I + S, and for unit y ∈ Y_N, (a_vI + S)y = (a_v + 1)y; by operator monotonicity of the
inverse, 0 ⪯ K ⪯ (1 + a_v)^{−1}I, so I ⪯ (I − K)^{−1} ⪯ (1 + 1/a_v)I. Finally ‖c‖² = (Xb)ᵀS(Xb), Xb is
the restriction of φ up to an additive constant (killed by S), and φ_v = mean of φ over N (φ is harmonic
at v), so ‖c‖² = Σ_{x∈N}(φ_x − φ_v)² = E_v(u,w). ∎

**Lemma 5′ (sums).** With L*_v = Σ_{x∼v}(e_x − e_v)(e_x − e_v)ᵀ: Σ_{u<w} E_v(u,w) = n·tr(Γ²L*_v) =: nσ_v,
Σ_v nσ_v = 2Kf(G), and ε_v := Σ_{w≠v} E_v(v,w) ≥ (n − 1)/d_v.

*Proof.* φ_x − φ_v = (e_u − e_w)ᵀy with y = Γ(e_x − e_v), and Σ_{u<w}((e_u − e_w)ᵀy)² = n‖y‖² − (1ᵀy)² =
n‖y‖². Summing over x gives n tr(Γ²L*_v). Σ_v L*_v = 2L and Γ²L = Γ on 1^⊥, so Σ_v nσ_v = 2n tr Γ = 2Kf.
For the current from v to w the edge currents φ_v − φ_x (x ∼ v) sum to 1, so their squares sum to ≥ 1/d_v. ∎

**Theorem 6 (expander bound).** Let G be 2-connected, a = a(G), a_* = min_v a(G − v) > 0. Then

  Σ_v Δ_v(G) ≥ (1 + 1/a_*) Σ_v ε_v − 2Kf(G)/a_* ≥ (1/a_*)·[a(G)(n − 1)Σ_v d_v^{−1} − 2Kf(G)].

Hence if a(G)(n − 1)Σ_v d_v^{−1} > 2Kf(G), and in particular if a(G)² > 2n/Σ_v d_v^{−1}, then G is not KS.
If moreover G is vertex-transitive, then Δ_v(G) > 0 for every v (no vertex is good).

*Proof.* Δ_v = Rt_G(v) − Σ_{u<w; u,w≠v}[r_{G−v}(u,w) − r_G(u,w)] (bookkeeping). By Lemma 5 the sum is
≤ (1 + 1/a_v)(nσ_v − ε_v) ≤ (1 + 1/a_*)(nσ_v − ε_v), since nσ_v − ε_v = Σ_{u<w; u,w ≠ v}E_v(u,w) ≥ 0. Sum
over v and use Σ_v Rt(v) = 2Kf and Lemma 5′. For the second inequality use 1 + a_* ≥ a(G) (Fiedler 1973:
a(G − v) ≥ a(G) − 1; proof: let x be a unit Fiedler vector of G − v (x ⊥ 1, xᵀL_{G−v}x = a(G − v)), extended
by x_v = 0; then x ⊥ 1 in ℝ^V and a(G) ≤ xᵀL_Gx = xᵀL_{G−v}x + Σ_{y∼v}x_y² ≤ a(G − v) + 1) and ε_v ≥ (n − 1)/d_v. Finally Kf(G) = nΣ1/μ_i ≤ n(n − 1)/a(G). For vertex-transitive G
all Δ_v are equal. ∎

**Corollary 7.** Since d_v ≤ n − 1, every graph with a(G)² > 2(n − 1) is not KS. In particular:
(i) if δ(G) > (n − 2 + √(2(n − 1)))/2 then G is not KS (because a(G) ≥ 2δ − n + 2: a(G) = n − λ_max(L(Ḡ))
and λ_max(L(Ḡ)) ≤ 2Δ(Ḡ) = 2(n − 1 − δ) by Gershgorin); (ii) a join G = H₁ + H₂ with
min(n₁ + a(H₂), n₂ + a(H₁))² > 2(n − 1), e.g. with min(n₁, n₂) > √(2(n − 1)), is not KS (the Laplacian
spectrum of a join is {0, n} ∪ {λ + n₂ : λ ∈ spec L(H₁) minus one 0} ∪ {μ + n₁ : μ ∈ spec L(H₂) minus one 0}).

Examples: the bound excludes all Johnson graphs J(m,k) (a = m, d = k(m − k) ≤ m²/4), Paley graphs P(q)
for q ≥ 9 (a = (q − √q)/2), and every d-regular graph with second adjacency eigenvalue ≤ 2√(d − 1)
(e.g. Ramanujan graphs) when d ≥ 12, since then a ≥ d − 2√(d − 1) and (d − 2√(d − 1))² > 2d. It is not tight for sparse graphs (C₅: a² ≈ 1.91 < 4).
Check: `scripts/check_expander_dense.py` (exact identities of Lemma 5′ and the lower half of Lemma 5;
floating-point check of the upper half and of Theorems 6 and 8).

## 6. A dense bound

**Theorem 8.** Let G be connected with a(G) > 1 and Laplacian eigenvalues μ_1, …, μ_{n−1} > 0. Then for
every vertex v (G − v is connected since κ(G) ≥ a(G) > 1)

  Δ_v(G) ≥ 1 − Σ_{i=1}^{n−1} (n − μ_i)/(μ_i(μ_i − 1)) ≥ 1 − (n(n − 1) − 2m)/(a(a − 1)).

So if Σ_i (n − μ_i)/(μ_i(μ_i − 1)) < 1, no vertex of G is good. The bound is attained by the source's
graphs G_s = K_s + K_{s+1,1} (spectrum {2s+2 (s+1 times), s+1 (s times), 0}: the sum equals 1 and the
odd-numbered vertices are good).

*Proof.* L(G) + L(Ḡ) = nI − J, so the non-zero eigenvalues of G are n − θ, θ running over the n − 1
eigenvalues of L(Ḡ) on 1^⊥; hence Kf(G) = Σ_θ n/(n − θ) = (n − 1) + Φ(Ḡ, n) with
Φ(F, N) = Σ_{θ∈spec L(F)} θ/(N − θ) (zero eigenvalues contribute 0). Likewise Kf(G − v) = (n − 2) +
Φ(Ḡ − v, n − 1). The Laplacian of F = (Ḡ − v) ∪ {v} is L(Ḡ) minus a positive semidefinite matrix, so by
Weyl θ_i(F) ≤ θ_i(Ḡ); f(θ) = θ/(n − 1 − θ) is increasing on [0, n − 1) and θ_max(Ḡ) = n − a(G) < n − 1,
so Φ(Ḡ − v, n − 1) ≤ Σ_θ θ/(n − 1 − θ) (sum over spec L(Ḡ)). Therefore
Δ_v ≥ 1 − Σ_θ θ/((n − 1 − θ)(n − θ)), and θ = n − μ. The second form uses Σθ = 2|E(Ḡ)| and
(n − 1 − θ)(n − θ) ≥ (a − 1)a. ∎

## 7. A local formula for Kf(G − v)

**Proposition 9.** Let G be connected on n vertices, v a non-cut vertex with N = N(v), d = |N|, Γ = L^+,
Γ₂ = Γ². Put, for a, b ∈ N,
M_ab = Γ_ab − Γ_av − Γ_vb + Γ_vv, (M²)_ab = (Γ₂)_ab − (Γ₂)_av − (Γ₂)_vb + (Γ₂)_vv + n(Γ_av − Γ_vv)(Γ_bv − Γ_vv),
s_a = n(Γ_vv − Γ_av), t_a = −n[(Γ₂)_av − (Γ₂)_vv] − n²Γ_vv(Γ_av − Γ_vv),
Z = [[M_NN − I, s],[sᵀ, n − 1 + n²Γ_vv]], W = [[M²_NN, t],[tᵀ, n²(Γ₂)_vv + n³Γ_vv²]] (size d + 1). Then

  Kf(G − v) = (n − 1)·[Rt_G(v) − 1 − tr(Z^{−1}W)],  Rt_G(v) = tr Γ + nΓ_vv,  Kf(G) = n tr Γ.

*Proof.* Let M = L_v^{−1} (grounded at v; M_uw = Γ_uw − Γ_uv − Γ_vw + Γ_vv). Then L_{G−v} = M^{−1} − P_N
and, with J the all-ones matrix on V∖v, Y = L_{G−v} + J/(n − 1) satisfies tr Y^{−1} = tr L_{G−v}^+ + 1, so
Kf(G − v) = (n − 1)(tr Y^{−1} − 1). Y = M^{−1} + UCUᵀ with U = [E_N, 1], C = diag(−I_d, 1/(n − 1)); Woodbury
gives tr Y^{−1} = tr M − tr((C^{−1} + UᵀMU)^{−1}UᵀM²U), and C^{−1} + UᵀMU = Z, UᵀM²U = W, using
M1 = n(Γ_vv1 − Γe_v)|_{V∖v}, 1ᵀM1 = n²Γ_vv, tr M = Rt_G(v) and, writing M̃ = (I − 1e_vᵀ)Γ(I − e_v1ᵀ),
M̃² = (I − 1e_vᵀ)[Γ² + nΓe_ve_vᵀΓ](I − e_v1ᵀ) (because Γ1 = 0). ∎

Check: `scripts/check_deletion_formula.py` (exact, 381 deletions). Only entries of Γ, Γ² on {v} ∪ N(v)
(and tr Γ) enter; for vertex-transitive graphs with explicit Γ (Cayley graphs of abelian groups via
characters, distance-regular graphs) this gives Kf(G − v) at cost O(n d²) or in closed form.

**Proposition 10 (distance-regular graphs).** Let G be distance-regular with intersection array
{b_0, …, b_{D−1}; c_1, …, c_D}, k = b_0, k_i the sizes of the distance layers. For a radial f with
Σ_l k_l f_l = 0, the radial solution of Lψ = f is k_i b_i(ψ_i − ψ_{i+1}) = Σ_{l≤i} k_l f_l (flux through
the layer boundary); Γ_i (value of Γ at distance i) solves this with f = e_v − 1/n and Γ₂ with f = Γe_v,
normalised by Σ k_l ψ_l = 0. Two neighbours of v are at distance 1 or 2, so M_NN and M²_NN lie in
span{I, A_loc, J − I − A_loc}, A_loc the adjacency matrix of the local graph (a_1-regular), and
tr(Z^{−1}W) = Σ_θ w(θ)/z(θ) + (2skt − z₂₂w_Σ)/(s²k), the sum over the k − 1 non-principal eigenvalues θ of
A_loc (θ = 0 when a_1 = 0), with z(θ) = (m_0 − 1 − m_2) + (m_1 − m_2)θ, w(θ) = (w_0 − w_2) + (w_1 − w_2)θ,
m_0 = 2(Γ_0 − Γ_1), m_1 = Γ_0 − Γ_1, m_2 = Γ_2 − 2Γ_1 + Γ_0, w_j = (Γ₂)_j − 2(Γ₂)_1 + (Γ₂)_0 + n(Γ_1 − Γ_0)²,
w_Σ = w_0 + a_1w_1 + (k − 1 − a_1)w_2, and the 2 × 2 block uses Z_NN1 = 0 (because M1_N = 1). This holds
at every vertex (no vertex-transitivity needed).

*Proof.* Γ is a polynomial in A, hence constant on distance classes; the layer equations are current
conservation. Substitute into Proposition 9 and diagonalise with respect to 1_N and the eigenvectors of
A_loc orthogonal to 1_N. ∎

Check: `scripts/check_drg_formula.py` (C₄–C₉, K₄,₄, K₆, Petersen, Q₃, Q₄, Heawood, Clebsch, O₄, H(2,3),
H(3,3), Paley(13): formula = direct exact computation).

## 8. Complete multipartite graphs

**Theorem 11.** Let G = K_{n_1,…,n_k}, k ≥ 2, n = Σn_i ≥ 3, G not a star. For v in a part of size n_j:

  Δ_v(G) = φ(n_j) − S,  φ(x) = (n + x − 2)/(n − x) + x(x − 1)/((n − x)(n − 1 − x)),
  S = Σ_i n_i(n_i − 1)/((n − n_i)(n − 1 − n_i)).

φ is strictly increasing on [1, n − 2]. Hence the set of good vertices is empty or the union of all parts of
one size; G is never KS: if all parts have size s then Δ_v = [k(k−1)s² − 2ks + 2]/[(k−1)s((k−1)s − 1)] > 0.

*Proof.* Spectrum: n (k − 1 times), n − n_i (n_i − 1 times), 0; so Kf(G) = k − 1 + nΣ_i(n_i − 1)/(n − n_i).
If n_j ≥ 2, G − v = K_{…, n_j − 1, …}; if n_j = 1 the part disappears (G − v is connected unless G is a
star). Subtracting term by term: for i ≠ j, n(n_i − 1)/(n − n_i) − (n − 1)(n_i − 1)/(n − 1 − n_i) =
−n_i(n_i − 1)/((n − n_i)(n − 1 − n_i)); the j-term gives (n + n_j − 2)/(n − n_j) (also when n_j = 1, where it
equals 1 and the part count drops by one). Thus Δ_v = φ(n_j) − S. With n = x + 2 + y, x = 1 + z
(y, z ≥ 0), the numerator of φ′ is (2y + 3)(y + z + 1)(y + z + 2) > 0. For K_{k×s}: φ(s) − kh(s) with
h(x) = x(x − 1)/((n − x)(n − 1 − x)) simplifies to the stated fraction; k(k−1)s² − 2ks + 2 = ks((k−1)s − 2) + 2
> 0 unless (k − 1)s = 1, i.e. G = K₂. ∎

Check: `scripts/check_multipartite.py` (formula vs exact Kf for all 1 177 (partition, part) cases with
n ≤ 12; closed form for K_{k×s}, k, s ≤ 12; monotonicity symbolic). This recovers the source's Prop. 7
structure (in K_{4,3} the part of size 3 is good: φ(3) = S = 5/2).

## 9. Strongly regular graphs

**Theorem 12.** C₅ is the only strongly regular KS graph. More precisely, for an SRG(n,k,λ,μ) and every
vertex v,

  Δ_v = Rt(v) − ((n − 1)/P)·Σ_θ W(θ)/Q(θ),  P = μn = k² − kλ + kμ − k + μ,
  Q(θ) = k² − kλ + kμ − 2k + λ − θ,  W(θ) = k² − 2kλ + 2kμ + k + λ² − 2λμ + μ² − μ + (2k − λ + μ)θ,

θ over the k − 1 non-principal eigenvalues of the local graph of v, Rt(v) = 2Kf/n,
Kf = k(k² − 2kλ + 2kμ − k + λ² − 2λμ + λ + μ² − 2μ)/μ². For triangle-free SRGs (λ = 0):

  Δ_v = N(k,μ) / [μ(k + μ − 2)(k² + kμ − k + μ)],
  N = k⁴ + 3k³μ − 5k³ + 3k²μ² − 9k²μ + 5k² + kμ³ − 3kμ² + 5kμ − k + μ³ − 2μ² + μ,

and N > 0 unless (k, μ) = (2, 1); for Moore graphs (μ = 1), Δ_v = k(k − 2)(k + 1)/(k² + 1).

*Proof.* The formula is Proposition 10 for D = 2 (symbolic computation in `check_srg_theorem.py`).
Triangle-free: for k ≥ 4 every coefficient of N as a polynomial in μ is a polynomial in k − 4 with
non-negative coefficients and positive constant term; N(2, μ) = (μ − 1)(3μ² + 7μ + 6); N(3, μ) =
4μ³ + 16μ² + 16μ − 12 > 0. Now let G be any KS SRG. Imprimitive connected SRGs are complete multipartite
(Theorem 11). For a primitive SRG with eigenvalues k > r > 0 > s: (a) conference graphs (n = 4t + 1,
r = (√n − 1)/2) with n ≥ 9 satisfy a(G)² = ((n − √n)/2)² > 2k = n − 1, excluded by Theorem 6; n = 5 is C₅.
(b) Otherwise r, s are integers with r ≥ 1, s ≤ −2 (r = 0 gives μ = k, s = −1 gives λ = k − 1, both
imprimitive), so 1 ≤ μ = k + rs ≤ k − 2r, i.e. a(G) = k − r ≥ (k + 1)/2, and ((k + 1)/2)² > 2k for k ≥ 6:
excluded by Theorem 6. (c) The feasible primitive parameter sets with k ≤ 5 are (5,2,0,1), (9,4,1,2),
(10,3,0,1), (16,5,0,2) (`check_srg_theorem.py` enumerates them). For λ = 0 the closed form applies; for
(9,4,1,2) the θ's lie in [s, r] = [−2, 1] (Cauchy interlacing; θ_2(A_loc) ≤ λ_2(A) = r and θ_min ≥ s) with
mean −λ/(k − 1), and W/Q = −(2k − λ + μ) + C/(Q(0) − θ) with C > 0 is convex increasing, so Jensen and the
chord bound give Δ ∈ [64/45, 8/5] > 0. ∎

(The interval method of (c) proves Δ_v > 0 at every vertex for all 1 195 feasible primitive parameter
sets with n ≤ 400 except C₅ [exploration, floating point; not needed for the theorem].)

## 10. Distance-regular families

**Proposition 13 [exact computation].** Δ_v ≠ 0 (exact rational arithmetic, Proposition 10) for: cycles
C_4–C_200 except C₅ (Δ < 0 for n ≥ 6); hypercubes Q_2–Q_120 (Δ > 0); odd graphs O_3–O_81 (Δ > 0); incidence
graphs of generalized triangles, quadrangles and hexagons of order (q,q), 2 ≤ q ≤ 150 (formula applied to
the feasible arrays, existence not needed; Δ > 0 except q = 2 for quadrangles and hexagons); Hamming
graphs H(d,q), 2 ≤ d, q ≤ 30 (Δ > 0); Johnson graphs J(m,k), m ≤ 40 (Δ > 0; all are excluded by Theorem 6
anyway); the 13 cubic distance-regular graphs (5 with Δ > 0, 8 with Δ < 0); Moore arrays {k, k−1; 1, 1},
k ≤ 300; triangle-free SRG arrays with k ≤ 120. Log: `certificates/drg_families.log`
(`scripts/drg_families.py`).

## 11. Exhaustive searches

**Exact criterion (`scripts/kscheck.c`).** For each graph G (graph6, n ≤ 64): reject if G or some G − v
is disconnected (cut vertex). Otherwise compute k_p(H) = N tr(A^{−1}) − 1ᵀA^{−1}1 mod p for H = G and
H = G − v (Lemma 1, A grounded Laplacian, LDLᵀ factorisation and inverse of the unit triangular factor
over 𝔽_p), p₁ = 1073741789 (fallback p₂ = 1073741783 if p | τ). Since k_p(H) ≡ α(H)τ(H)^{−1}, k_p(G) ≠
k_p(G − v) implies p ∤ α(G)τ(G − v) − α(G − v)τ(G), hence Kf(G) ≠ Kf(G − v): a rigorous rejection.
Vertices are tested in order until one fails; a graph passing all vertices for both primes is reported
as SURVIVOR and checked with exact rationals; a graph where both primes divide some τ would be reported
INCONCLUSIVE (never happened). In vertex-transitive mode (`-t`) only vertex 0 is tested. Validation:
`scripts/test_kscheck_modular.py` compares the residues with exact rational values of Kf for 119 graphs
(2 126 Kirchhoff indices, n ≤ 40).

**Theorem 14 [exact computation].** C₅ is the only KS graph with at most 11 vertices.

| n | 2-connected graphs (geng −C) | survivors |
|---|---|---|
| 3 | 1 | 0 |
| 4 | 3 | 0 |
| 5 | 10 | 1 (C₅) |
| 6 | 56 | 0 |
| 7 | 468 | 0 |
| 8 | 7 123 | 0 |
| 9 | 194 066 | 0 |
| 10 | 9 743 542 | 0 |
| 11 | N11 | 0 |

(counts = OEIS A002218; logs `certificates/search/geng_C_n*.log`; reproduce with
`scripts/run_geng_search.sh 3 11`.)

**Theorem 15 [exact computation].** There is no vertex-transitive KS graph on 6 ≤ n ≤ 47 vertices.
Every graph of the Holt–Royle census (Zenodo 10.5281/zenodo.4010122, all vertex-transitive graphs on
10–47 vertices, all valencies, including disconnected ones) was tested (`scripts/run_census.sh`); counts per
n agree with OEIS A006799 and are in the table below; n ≤ 9 is covered by Theorem 14.

CENSUS_TABLE

## 12. Numerical explorations (not proofs)

* Square tori C_a × C_a: Δ/Kf changes sign between a = 6 (+0.0049) and a = 10 (−0.0011); exact:
  Kf(C₈ × C₈) = 554632/357, Kf(C₈ × C₈ − v) = 17147406699025/11037120528 (difference ≈ −0.0213).
  Along every family examined (tori, prisms, generalized Petersen graphs GP(n,k), circulants) Δ changes
  sign at small order and then keeps the sign of the "dimension" (negative for 1- and 2-dimensional
  families, positive for expanders), so near-misses are frequent and exact solutions would be coincidences.
* NUMERICAL_EXTRA

## 13. What is not proved

The problem remains open: we do not exclude KS graphs with n ≥ 12 in general, vertex-transitive ones with
n ≥ 48, nor sparse graphs with cyclomatic number ≥ 6. The searches of §11 are exact but finite; the
explorations of §12 are floating point.
