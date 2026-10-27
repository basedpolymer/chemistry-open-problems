# The Global Attractor Conjecture through transpose duality

**Agent G: GPT-6-Astra, effort xhigh. 27 October 2026.**

Status: working mathematical report. Statements marked **proved** have proofs below; exact checks are corroboration, not substitutes for those proofs. Statements attributed to published papers/preprints are explicitly identified. No novelty or independent verification is claimed. No private correspondence is quoted or used as mathematical evidence.

## 1. Definitions and scope

A network N has n species and r labelled reactions y_j -> y'_j, with nonnegative integer complexes. Set Ye=(y_1|...|y_r), Yp=(y'_1|...|y'_r), Gamma=Yp-Ye, S=im Gamma, s=rank Gamma, q=n-s. Its mass-action system is x'=Gamma v(x), v_j(x)=k_j x^{y_j}, k_j>0. A positive compatibility class is (x0+S) intersect R^n_{>0}; its closure is P=(x0+S) intersect R^n_{>=0}, x0>0. Codimension means **n-s**, not the codimension of an individual boundary face.

The Iosif dual N* has r species, n labelled reactions a_i -> b_i, where a_i is row i of Ye and b_i is row i of Yp. Thus Gamma*=Gamma^T. Self-loops and reaction labels are retained when taking duals or subnetworks. This makes (N*)*=N literally. Removing a self-loop afterwards does not affect its ODE but changes the subsequent transpose construction.

H1 (conservativity) means there is w>0 with Gamma^T w=0. H2 (consistency) means there is v>0 with Gamma v=0, equivalently a positive equilibrium for **some** rate vector. H3 means the rows of Gamma are pairwise distinct. A semiflow is a nonzero nonnegative vector in the relevant kernel. An internal cycle is an extreme ray of ker Gamma intersect R^r_{>=0}; its support is a minimal semiflow support. A siphon Sigma is a nonempty species set such that every reaction producing a species in Sigma consumes a species in Sigma. A T-trap J is a nonempty reaction set for which the union of produced species is contained in the union of consumed species.

Let m be the number of distinct complexes, ell the number of linkage classes, I the complex/reaction incidence matrix, and delta=m-ell-s the deficiency. Weak reversibility (WR) means that every reaction belongs to a directed cycle in the complex graph. Complex balance at c>0 means I v(c)=0, stronger than Gamma v(c)=0. Persistence here means liminf_{t->infinity} x_i(t)>0 for every species of every positive trajectory under consideration; the trajectories in our main statements are bounded, so this is equivalent to having no boundary omega-limit point.

For context, preclusters P(N) are reaction components obtained from common-source pairs by the span rule of Iosif's Definition 19, and MIPS(N) are species classes with positively proportional columns in a matrix of extreme conservation laws. Result 17 already proves exact refinement criteria and the discreteness of P(N*) when N is WR and H3 holds. Those partitions are not used as proxies for dynamics in this report.

**Answer in scope.** Transposition gives exact and useful descriptions of boundary tests in terms of reaction sets, and a restricted codimension-one certificate for strong endotacticity when the number of species is prime (Proposition 4.1). It does not give a general transfer of persistence/global attraction, even when all H1--H3 hold and the receiving class contains a positive equilibrium. Small q controls the dual's *stoichiometric flux nullity*, not its dynamical dimension. Sections 3--5 provide precise obstructions; Sections 6--7 give the meaningful boundary dictionary. No new GAC case beyond established dynamical theorems is claimed.

## 2. Elementary invariants and the exact WR characterization

### Proposition 2.1 — rank, bounded classes and deficiency [proved]

1. N* has r species, n reactions and rank s. Consequently dim ker Gamma*=n-s=q and dim ker (Gamma*)^T=r-s. H1 and H2 exchange.
2. A compatibility polyhedron containing a positive point is bounded iff N is conservative. Hence under H1 and H2 both N and N* have bounded classes. Without both assumptions boundedness of classes is not invariant: 0 <-> A has an unbounded class, whereas its dual U -> V has bounded classes.
3. Write beta*=dim ker I*=n-m*+ell*. Then

   **q = beta* + delta*.**

   If N* is WR and nonempty, beta* >= ell*, so ell*+delta* <= q. In particular q=1 and N* WR force ell*=1, delta*=0 and beta*=1: the dual complex graph is one directed cycle (a loop is allowed in the trivial rank-zero case).

**Proof.** The first assertion is transposition and rank-nullity. The recession cone of P is S intersect R^n_{>=0}. Stiemke's alternative says this is {0} iff S^perp contains a strictly positive vector. Equivalently a positive conservation law bounds each coordinate, while a nonzero nonnegative recession direction makes every such P unbounded. For (3), rank I*=m*-ell*, and Gamma*=Y* I*, so rank Gamma*<=rank I* and delta*>=0. Rank-nullity gives the identity. Each WR linkage class with an edge has at least one independent circulation. If the total graph cycle nullity is one, there is only one linkage class; in its strongly connected graph every vertex has at least one incoming and outgoing edge, and the edge/vertex counts agree, forcing exactly one of each. This proves the cycle assertion. This argument includes parallel edges and self-loops with their labels. QED.

The number of complexes, number of linkage classes and deficiency are not preserved. The dual counts are exactly those of the **row graph** with vertices the distinct rows of Ye and Yp and edge i from a_i to b_i. Section 3 changes (m,ell,delta)=(2,1,0) into (4,2,1) without changing Gamma. In contrast A <-> B has the same n, r, Gamma, m, ell and delta as the network in Section 3, but its dual still has (m*,ell*,delta*)=(2,1,0). Thus even these primal counts and Gamma together do not determine the dual counts.

### Proposition 2.2 — general and specialized tests for dual WR [proved]

For any N, N* is WR iff there are u_i>0 satisfying, for every row-vector vertex c,

    sum_{i:a_i=c} u_i = sum_{i:b_i=c} u_i.                  (2.1)

If N itself is WR and satisfies H3, its source rows a_i are pairwise distinct, as are its product rows b_i. In this case the following are equivalent:

- N* is WR;
- the sets of source and product rows are equal;
- there is a unique permutation pi of [n] with b_i=a_{pi(i)} for all i.

The dual graph is then the disjoint union of the cycles of pi, m*=n, ell*=c(pi), and

    delta* = q-c(pi).

Every pi-orbit C gives a primal conservation law 1_C. In particular q=1 forces pi to be one n-cycle.

**Proof.** A finite directed graph has a strictly positive circulation iff every edge lies on a directed cycle. In one direction sum the characteristic circulations of cycles covering all edges. Conversely a positive circulation on the acyclic condensation graph cannot cross between strongly connected components: take a source component and use flow balance, then induct. This proves (2.1). If a_i=a_j in a WR N, every product complex is also a source, so b_i=b_j; hence Gamma_i=Gamma_j, contradicting H3. Apply the same argument to the reversed WR network to get distinct b_i. Thus the dual has at most one edge entering or leaving any vertex. It is WR precisely when every target is a source, equivalently the two finite sets agree. The permutation and counts follow. Summing b_i-a_i on a permutation cycle gives zero, so 1_C is a conservation law. QED.

This is substantially more restrictive than consistency: (2.1) implies Gamma^T u=0, but the latter forgets balance at the individual row complexes. It also shows that the WR-to-WR locus under H3 is specified by finitely many row-permutation equalities, not by a general conservation argument.

One can describe the original complex graph as well: `y'_i = y_{pi(i)}` in species coordinates for every reaction. Thus one fixed coordinate permutation sends every source to its product. Once parallel copies are merged, each source has only one outgoing edge; original WR makes each linkage class a complete orbit of that coordinate permutation. In particular **unweighted total concentration is conserved**. The dual-WR hypothesis is substantially stronger than the existence of an arbitrary weighted positive conservation law.

### Complex balance and endotacticity [proved]

A network admits complex balance for some positive rates iff it is WR. Indeed a complex-balanced flux is a positive incidence circulation; conversely choose such a circulation, take c=1, and use the circulation as the rate vector. Thus (2.1) is also the exact condition that N* admit *some* complex-balanced kinetics. It does not transfer the particular rates or equilibrium of N.

Here is an exact row-geometric description of the dual's endotacticity. Discard self-loop edges for these geometric definitions, as in the standard reaction-network theorems, and let E={i:Gamma_i!=0} be the remaining dual reaction labels. For z in R^r set h=Ye z and d=Gamma z. Using a sweep from large to small scalar product, N* is endotactic iff for every z with d!=0, every index maximizing h_i over {i:d_i!=0} has d_i<0. It is strongly endotactic iff for every such z every index maximizing h_i over **all i in E** has d_i<=0 and at least one maximizing index has d_i<0. These are just the usual supporting-hyperplane tests applied to sources a_i and displacements b_i-a_i; an outward active reaction must have an inward active reaction at strictly greater source height. Algebraic duals still retain the original labels. These tests give a precise image of the hypotheses, but the properties are not invariant: the counterexample below starts with a reversible single-linkage network and ends with a non-endotactic one.

## 3. A dynamical obstruction in one dimension [proved]

Consider the reversible network

    M: 3A -> 2A+B,   2A+B -> 3A.

Its dual is

    D=M*: 3U+2V -> 2U+3V,   V -> U.

Both matrices are

    Gamma = [[-1, 1], [1, -1]].

Both satisfy H1--H3, have s=q=1, positive conservation and flux vectors (1,1), and bounded compatibility segments. M has two complexes, one linkage class and deficiency zero. It is reversible and strongly endotactic. D has four complexes, two linkage classes and deficiency one. It is not WR. With sweep vector (0,1), its first source has height 2, its second source height 1, and the first reaction increases V. Hence D is not endotactic.

Put **both rate constants equal to one in both networks**. For M, with a+b=T>0,

    a' = a^2(b-a) = a^2(T-2a).

The sign is positive for 0<a<T/2 and negative for T/2<a<T. Monotonicity and uniqueness therefore give a(t)->T/2 from every positive initial point. The equilibrium is complex balanced because the two fluxes agree. This proves global attraction directly, independently of any general GAC theorem.

For D, u+v=T and

    v' = u^3 v^2-v = v[(T-v)^3 v-1].                     (3.1)

In the class T=2 there is the positive equilibrium (u,v)=(1,1). Nevertheless start at (15/8,1/8). The function h(v)=(2-v)^3 v has derivative (2-v)^2(2-4v)>0 on [0,1/8], and

    h(1/8)=3375/4096 < 1.

Consequently (0,1/8] is forward invariant for v and

    v' <= -(721/4096) v,
    0<v(t)<= (1/8) exp[-(721/4096)t].

Thus (u(t),v(t))->(2,0). Global existence follows from the conserved bounded total; positivity follows from the factor v and the mass-action vector field. D is nonpersistent **in a class that contains a positive equilibrium**. All numbers and the transpose are checked exactly in `scripts/check_gac_duality.py`.

This refutes the unrestricted implication

    N* satisfies a known GAC case  =>  N is persistent/globally attracting

under H1--H3, even when a positive equilibrium in the tested class is assumed. Take N=D and N*=M. It is not a counterexample to the GAC: D is not complex balanced, nor WR. It does not exclude a transfer theorem with additional hypotheses on both systems. It proves that transposition alone supplies no such dynamical theorem.

### Common-factor cancellation is not a repair [proved]

Cancelling species common to both sides of each reaction in D produces U -> V, V -> U. This is linear, reversible and globally attracting in every positive class, but D itself has the extinction trajectory just proved. The same cancellation in M produces A <-> B. All these networks have the same Gamma. Thus neither the stoichiometric subspace nor a canonical network obtained from the positive/negative parts of Gamma determines persistence. Cancelling common factors is not a conjugacy or a common positive time change of the whole vector field: each reaction loses a different monomial factor.

## 4. What small codimension does and does not buy [proved]

For every n>=2, the monomolecular directed cycle

    X1 -> X2 -> ... -> Xn -> X1

is conservative, consistent, WR and deficiency zero; has distinct rows; and has q=1, s=n-1. Its dual is again a monomolecular directed cycle with rank n-1. Proof: the reaction columns span the sum-zero hyperplane, and the incidence matrix is transposed by the dual. This gives arbitrarily large primal and dual dynamical dimensions at the **smallest possible positive conservative codimension**.

Transposition, reversal of all arrows, relabelling, cancellation of common source/product factors, and invertible diagonal changes of coordinates preserve rank. Any sequence of these operations therefore cannot turn q<=d into rank<=d. A genuine injective dynamical change of coordinates likewise cannot reduce the dimension of a compatibility class. A noninjective projection or lumping needs additional dynamical structure, which is not provided by a transpose identity.

There is a positive but limited statement: if q=1 and N* is WR, Proposition 2.1 puts N* in Anderson's single-linkage case (indeed in deficiency zero). This proves global attraction **for N*** for all positive rates. Section 3 shows why the additional step back to N cannot be inferred from H1--H3, even in this situation. If N is itself complex balanced, its convergence still needs a theorem about N or a separately proved transfer mechanism.

### Proposition 4.1 — a restricted positive route in codimension one [proved]

Suppose N is WR, has pairwise distinct rows of Gamma, has no zero-vector reactions, and q=1. Suppose also that N* is WR and the number n of species is **prime**. Then N is strongly endotactic. Therefore, by Gopalkrishnan--Miller--Shiu, every mass-action system on N is permanent; if its chosen rates are complex balanced, it has the global attractor property.

This uses information extracted from the dual to recognize a known dynamical class of the **original** network. It does not infer convergence of N from convergence of N*. It is a structural sufficient condition, not a claim of a new GAC case beyond the strongly endotactic theorem.

In fact H3 is redundant in this particular proposition: q=1 and dual WR imply by Proposition 2.1 that the dual graph is one deficiency-zero cycle, so its incidence and stoichiometric kernels both equal span(1,...,1). Thus ker Gamma^T=span(1,...,1); equality of two rows of Gamma would put e_i-e_j in that kernel, which is impossible. H3 is retained in the statement to match the project's hypotheses.

**Proof.** Proposition 2.2 gives a single n-cycle pi on the species. After relabelling, its linear action T is the cyclic shift. Every original linkage class, ignoring parallel reaction copies, is the complete T-orbit of a complex y. The reaction is nonzero, so d=Ty-y is a nonzero rational vector in the sum-zero space V.

For prime n, V is an irreducible Q[T]-module. Here is a direct justification. The minimal polynomial of the cyclic shift on Q^n is t^n-1, and on V it is the cyclotomic polynomial Phi_n(t)=1+t+...+t^{n-1}. This polynomial is irreducible over Q: Phi_n(t+1)=((t+1)^n-1)/t is Eisenstein at the prime n. If a polynomial p of degree <n-1 annihilated a nonzero d in V, then gcd(p,Phi_n)=1 and Bezout's identity would force d=0. Hence d,Td,...,T^{n-2}d are independent. All of them are reaction vectors of this linkage class. Thus **each linkage class spans V=S**.

For any w not orthogonal to S, w is nonconstant on the complexes of every linkage class. A linkage class containing a globally maximal-height source has an edge leaving its maximal-height plateau: follow a directed path from a maximum to a lower vertex. That edge starts at the global maximum and points strictly inward. Every reaction starting at the global maximum points weakly inward, since its target is another source. This is strong endotacticity. The quoted permanence/GAC conclusion is exactly the known strongly endotactic theorem. QED.

**A nontrivial five-species example [proved, exactly checked].** On indices modulo 5 take the two disjoint linkage cycles

    Xi+X(i+2) -> X(i+1)+X(i+3),
    2Xi+2X(i+2) -> 2X(i+1)+2X(i+3),       i=1,...,5.

There are ten reactions, ten complexes, two linkage classes, s=4, q=1 and delta=4. The dual is one five-edge cycle with deficiency zero. All-one rates make N complex balanced at the all-one state. Proposition 4.1 proves strong endotacticity and hence permanence for all rates. The proper set {X1,X2,X3} is a critical siphon: it meets every source and every product, while the only conservation-law ray has full support. It is realizable on the total-5 class, for example at (0,0,0,2,3); its face has codimension 3, so it is neither a facet nor a vertex. There are no blind reactions, but dim S0=1, so Wiuf's condition (a) fails. This example lies outside the single-linkage, rank<=3, ADLS, and all-siphon-facets/vertices hypotheses; the operative known theorem is strong endotacticity.

**The primality condition cannot simply be removed [proved, exactly checked].** On four species take the monomolecular cycle A->B->C->D->A together with A+C <-> B+D. Then q=1, H3 holds and the dual is WR (the common permutation is the four-cycle). In direction w=(3,2,1,2), the two bimolecular sources both have height 4 and stationary reactions; all monomolecular sources have height <=3 and some of their reactions have nonzero projection. There is no strictly inward reaction at the globally maximal source height, so N is not strongly endotactic. This refutes only the proposed geometric conclusion without primality, not persistence of this example (which the monomolecular cycle already makes easy by ADLS).

### Gale duals and conservation-coordinate projections [proved]

Let L=S^perp and let B have row space L and rank q. Then

    d(Bx)/dt = B Gamma v(x)=0.

The q conservation coordinates carry the zero dynamics. They cannot encode the nonconstant trajectories in an s-dimensional class.

An orthogonal-complement construction that realizes L itself as a stoichiometric subspace in the same R^n has rank q, but if N is conservative, L contains w>0. Its positive classes are consequently unbounded, containing x+t w for all t>=0. It cannot have a positive conservation law, since such a law would belong to L^perp=S and have strictly positive scalar product with w. This obstructs a **conservative** Gale-dual realization in the same positive orthant. It does not say that no nonconservative complex-balanced realization of L exists, or exclude every other construction called a Gale dual.

For a fixed complex-balanced system with equilibrium c, the positive equilibrium variety is

    E = {c exp(l): l in L},

of dimension q. This is a family of equilibria across classes, not a q-dimensional system carrying the original trajectories. Birch's intersection theorem gives one point of E in each positive class. In the conservative case this can also be seen by minimizing the strictly convex function

    F_c(x)=sum_i [x_i log(x_i/c_i)-x_i+c_i]

on the compact P: a boundary minimizer is ruled out by its directional derivative toward a positive point, and the interior critical equation is log(x/c) in L. Uniqueness of the minimizer is not convergence of all trajectories to it. In particular F_c has a finite continuous extension to boundary points because x log x -> 0; it is not a boundary barrier. A small-dimensional equilibrium variety does not supply the missing boundary argument.

## 5. A noninvariant siphon criterion [proved]

Let

    H: 2A <-> 2B,   2A <-> A+B,

with the forward/backward ordering displayed. It is WR, conservative and consistent and satisfies H3. Its only minimal siphon is {A,B}, which supports the law a+b. Therefore it satisfies the Angeli--De Leenheer--Sontag (ADLS) criterion.

Its dual, with four species U1,...,U4, has reactions

    2U1+2U3+U4 -> 2U2+U3+2U4,
    2U2+U4 -> 2U1+U3.

The singleton {U4} is a siphon, but it supports no conservation law: row 4 of the dual stoichiometric matrix is (1,-1). Thus the ADLS criterion is not self-dual, even starting with a conservative WR network. The dual is also conservative, consistent and has pairwise distinct rows. This is a failure of the criterion; no assertion about this dual's persistence is needed.

## 6. A useful duality: boundary tests become reaction-subnetwork tests

### Proposition 6.1 — siphons, cycle coverage and persistence [proved]

Under the identification of species of N with reactions of D=N*, Sigma is a siphon of N iff it is a T-trap of D. Moreover

    every siphon of N contains a semiflow support
    iff every T-trap of D contains an internal-cycle support.             (6.1)

If N is conservative, either equivalent condition implies persistence for every positive mass-action rate vector. If the chosen kinetics of N are complex balanced, it also implies global attraction to the unique positive equilibrium in each positive class.

**Proof.** Dual species j is produced by some reaction in Sigma iff some i in Sigma has (Yp)_ij>0; it is consumed by Sigma iff some i in Sigma has (Ye)_ij>0. These are exactly the two implications in the siphon definition. The cones of conservation laws of N and internal fluxes of D coincide, since Gamma*=Gamma^T. A nonzero nonnegative vector supported in Sigma decomposes into extreme rays; every ray in its nonnegative decomposition has support contained in Sigma. This proves (6.1).

For completeness, the ADLS boundary argument is as follows. Conservativity bounds a positive trajectory and makes its omega-limit set compact. If a point z of that set has zero set Sigma, Sigma must be a siphon. Otherwise a reaction with all its reactants positive at z produces a missing species i, giving f_i(z)>0: consumption terms of i vanish there. A short backward segment of the complete orbit through z inside the invariant omega-limit set would have negative coordinate i, a contradiction. If w>=0 is a nonzero conservation law supported in Sigma, w.x(t)=w.x(0)>0, while w.z=0, also a contradiction. Thus the omega-limit set avoids the boundary. For complex-balanced kinetics the standard strict free-energy dissipation in the interior, compactness, and the unique positive equilibrium give convergence. QED.

This is a true persistence theorem stated entirely through reaction sets and internal cycles of the dual. It is a re-expression of ADLS, not a new dynamical theorem and not an application of persistence of D.

### Proposition 6.2 — face codimension is dual subnetwork rank [proved]

Fix x0>0. For Sigma subset [n], suppose the face

    F_Sigma={x in P: x_i=0 for i in Sigma}

contains a point with x_i>0 for every i outside Sigma. Then

    codim_P F_Sigma = rank Gamma[Sigma,:]
                   = rank Gamma_D[:,Sigma].                            (6.2)

The last quantity is the stoichiometric rank of the subnetwork of D retaining **only the reactions in Sigma**, with all dual species retained.

**Proof.** P has dimension s because x0>0. Near a point in the relative interior of this face, the only local linear restrictions in addition to x-x0 in S are x_Sigma=0. Thus its tangent space is ker(pi_Sigma|_S), with dimension s-rank(pi_Sigma|_S). The image of this projection is im Gamma[Sigma,:], yielding (6.2); transpose gives the final equality. QED.

For a realizable siphon, the corresponding set of dual reactions is a T-trap. Hence **facets correspond to rank-one T-traps; vertices to rank-s T-traps**, provided realizability in the specified primal class is retained. The ADLS noncritical siphons give empty faces. Face codimension cannot be inferred merely from the cardinality of Sigma or from q.

Consequently the Anderson--Shiu GAC criterion has this exact equivalent formulation: for every nonempty realizable T-trap Sigma of D, the reaction subnetwork D|Sigma has rank either 1 or s. This is a useful class of dual certificates, but the theorem being applied concerns the original network and its faces.

### Proposition 6.3 — Wiuf's conditions in dual matrices [proved translation]

Assume N is WR. Fix a siphon Sigma, and let J0 consist of reactions of N whose source has no species in Sigma. By WR and the siphon property their targets also avoid Sigma. Put M=Gamma_D=Gamma^T. Wiuf's spaces, extended by zero on Sigma when necessary, are

    S0(Sigma)=S intersect {x_Sigma=0},
    Sbl(Sigma)=im Gamma[:,J0] subset S0(Sigma).

Then his condition (a) is exactly

    rank M[J0,:] + rank M[:,Sigma] = s.                                (6.3)

For condition (b), project the original complexes to coordinates Sigma, retain labelled reactions and form the embedded linkage classes. A nonzero embedded complex is Sigma-minimal if no distinct active embedded complex is componentwise <= it; an active linkage class is Sigma-minimal if it contains such a complex. For each such class L, let J_L be its original reaction labels. Wiuf's condition (diamond) is exactly

    ker M[J_L,Sigma] intersect R^{Sigma}_{>0} = empty.                  (6.4)

Thus the indicated **species-and-reaction restricted dual network is not consistent**.

**Proof.** The first equality is rank-nullity as in Proposition 6.2, since dim S0=s-rank Gamma[Sigma,:] and dim Sbl=rank Gamma[:,J0]. The inclusion makes dimension equality sufficient. For (b), the embedded class has stoichiometric subspace im A with A=Gamma[Sigma,J_L]. Wiuf asks that im A contain a nonzero nonnegative vector. By Stiemke's alternative this holds iff ker A^T contains no strictly positive vector. Since A^T=M[J_L,Sigma], (6.4) follows. All restrictions are ordinary submatrices, with parallel/self-loop reactions retained. QED.

Projection and transpose obey `(projection of N to Sigma)* = (D restricted to reactions Sigma)` before removal of loops or duplicate labels. However embedded **complex linkage classes** of N do not become the usual linkage classes of D: they must be reconstructed from the projected primal complexes. Likewise the ordering that defines Sigma-minimality is an order on those complexes. This is the remaining geometry in (6.4).

Applied to every realizable siphon of a specified positive class, (6.3)--(6.4) are an equivalent way to check the hypotheses of Wiuf's Theorem 7.2. The conclusion of that theorem is used **as stated in his preprint**, not proved here. In particular this translation does not claim the theorem holds for arbitrary T-traps independent of class realizability, nor that its hypotheses are self-dual.

## 7. Known cases: hypothesis, dual image, verdict

In the first five rows, the middle column explains what applying the named theorem to **D=N*** would require of N. In the last three rows it rewrites a theorem **about N** as a certificate in D. This distinction prevents an unsupported dynamical transfer.

| Known result and actual hypothesis | Exact image under transpose | Verdict |
|---|---|---|
| Anderson (2011): a complex-balanced system with one linkage class has the global attractor property. | The row-complex graph a_i -> b_i is strongly connected, and the chosen dual rates are complex balanced. Under original WR and H3 this means pi is one cycle. | Natural and decidable row-graph class, but proves convergence only in D. If q=1 and D is WR, this case always applies to D. Section 3 blocks an unconditional return to N. |
| Pantea (2012): GAC for complex-balanced systems with stoichiometric dimension <=3; bounded WR trajectories in dimension <=2 are persistent. | rank Gamma^T=s is unchanged; WR/complex balance must still be checked in D. | No new small-codimension case. |
| Craciun--Nazarov--Pantea (2013): two-species endotactic systems with rates bounded above and away from zero are permanent; GAC for three-species systems. | D has r species; imposing r<=2 or r<=3 also imposes s<=2 or s<=3 on N. | Ambient-species bounds should not be confused with Pantea's stoichiometric-rank bound. No q-to-s conversion. |
| Gopalkrishnan--Miller--Shiu (2014): strongly endotactic systems are permanent; in particular complex-balanced strongly endotactic systems satisfy GAC. | The row-geometric strict supporting-hyperplane test in Section 2; complex balance for the chosen dual kinetics if a GAC conclusion is wanted. | Not invariant; q alone gives no test. However q=1, both orientations WR, and prime n give the original-network certificate of Proposition 4.1. M in Section 3 is strongly endotactic while D is not endotactic. |
| Xu (2026; arXiv v2 2025): first-order endotactic mass-action systems have a unique positive equilibrium in each positive class and exponential global attraction. | Every source row of Ye has sum <=1, together with the dual row-endotacticity test. | A severe row-sparsity condition. Under original WR and H1 it already forces a class covered by the siphon criterion; under H3 it reduces to disjoint monomolecular cycles (Proposition 7.1). |
| Anderson--Shiu (2010), Theorem 4.6: complex-balanced N whose siphon faces are empty, facets or vertices. Their WR persistence theorem uses empty/facet siphon faces and bounded trajectories. | Every realizable dual T-trap subnetwork has rank 1 or s (GAC); rank 1 in the bounded WR persistence criterion. | Exact and useful boundary certificates by (6.2), not persistence of D. Ambient q does not bound these ranks. |
| Wiuf (2026), Theorem 7.2: for every realizable siphon in the class, S0=Sbl and each Sigma-minimal active embedded linkage class satisfies (diamond). | Rank identity (6.3) and absence of a positive kernel vector in the submatrices (6.4), with the same realizability and embedded-complex ordering. | Meaningful exact reformulation. No claim of self-duality or a new GAC case; the preprint's theorem has not been independently audited here. |
| Angeli--De Leenheer--Sontag (2007), Theorem 2: conservative N and every siphon contains a nonnegative conservation-law support. | Every T-trap of D contains an internal-cycle support. | A proved dual structural certificate for persistence of N, and GAC if N has complex-balanced kinetics. It is not the ADLS criterion of D; Section 5 refutes self-duality. |

### Proposition 7.1 — the first-order-dual class already has persistence [proved]

Suppose N is WR and conservative, has no zero reaction vectors, and every row of Ye has sum <=1. Ignore unused species, which have constant positive coordinates. Then distinct nonzero source complexes have disjoint supports, each coefficient is one, and each source occurs in exactly one labelled reaction. The original reaction graph is a disjoint union of directed cycles on disjoint nonempty species groups. Every siphon contains a conservation-law support, so ADLS proves persistence. Every positive rate choice is complex balanced. If H3 also holds, each complex contains exactly one species, so N is a union of monomolecular cycles.

**Proof.** An integer nonnegative row with sum <=1 has at most one nonzero entry, equal to one. Thus each participating species appears in exactly one source. There is no zero complex in a nontrivial conservative WR component: its positive conserved mass would have to equal the mass of every other complex in that component. Sources are therefore nonempty, disjoint groups, and cannot repeat. In a WR graph every product is a source and every vertex has outdegree one, so its components are cycles. Distinct source supports are disjoint, hence their complex vectors are linearly independent. This gives deficiency zero. Alternatively complex balance can be constructed on each cycle by choosing one common positive flux and setting the product of concentrations in each group to that flux divided by its outgoing rate constant.

If a siphon contains a species in one group of a cycle, the reaction producing that group forces it to meet the predecessor group; induction forces it to meet every group of that cycle. Pick one species from each group and put coefficient one on these chosen species and zero elsewhere. Every reaction on that cycle consumes and produces exactly one chosen species, and every other reaction consumes and produces none. This is a conservation law supported in the siphon. Finally two species in one group have identical rows of Gamma: both are consumed at its outgoing edge and produced at its incoming edge. H3 excludes such pairs. QED.

Self-loops can be handled separately by retaining their species as constant or deleting that trivial component before the statement. The theorem is not claimed as a new persistence class; its purpose is to show why the simple-looking first-order dual route supplies no additional GAC case under the working hypotheses.

## 8. The Lyapunov route: what is transferable

Blanchini--Giordano (2022), DOI 10.1109/LCSYS.2021.3081369, work with a class of monotone systems determined by the **net** matrix Gamma. Their Assumption 3 requires a reaction to depend on a species iff its net stoichiometric coefficient is negative; their Assumption 4 makes those partial derivatives strictly positive. Their Theorem 1 transfers **structural local stability** across Gamma -> Gamma^T. Their Theorem 3 transfers a common polyhedral Lyapunov norm in concentrations to its polar norm in the dual's **rates**, and conversely. Theorem 3/Remark 2 and the formulas were checked on PDF page 5; the kinetic assumptions are on page 2.

For mass action, the derivative factor is

    G(x)=Dv(x)=diag(v(x)) Ye^T diag(x)^{-1}.

The dependence part of their assumption is equivalent to `(Ye)_ij>0 iff Gamma_ij<0`. It fails for catalysts and for autocatalytically produced reactants. This **support condition is itself invariant under transposition** and allows partial net consumption, for example a reactant coefficient 3 and product coefficient 2 in one entry. Pure networks, with Ye and Yp disjoint entrywise, are a sufficient subclass: there Ye=(-Gamma)_+, Yp=(Gamma)_+, so even the canonical source/product matrices agree. Purity is not necessary for the monotone dependence condition. For general mass-action networks the support condition fails; the extinction example has a reactant that is net produced, so it lies outside this class.

Even at the algebraic Jacobian level, matching a transpose does not match equilibria. At any positive x and flux v, choosing dual concentration z_j=1/v_j and dual flux w_i=1/x_i yields

    H=diag(w) Ye diag(z)^{-1}=G^T.

Then Gamma G and Gamma^T H have the same nonzero eigenvalues (use AB/BA and transpose). But `Gamma v=0` does **not** imply `Gamma^T(1/x)=0`. The matched dual point is generally not an equilibrium. This identity is valid without the monotone support restriction; that restriction concerns the scope of the structural theorem. This observation was also recorded independently by agent L and has been re-derived here.

The pseudo-Helmholtz F_c is strictly convex, smooth in the interior, not homogeneous and not a norm. Its Hessian is diag(1/x_i), and its Fenchel conjugate on R^n is sum_i c_i(exp(p_i)-1), not a polyhedral dual norm. Thus applying the norm-duality theorem to F_c is not justified. The fact that their Remark 4 extends to arbitrary **norms** does not remove this issue.

There is a second, separate issue in obtaining a persistence conclusion. Their Remark 2 requires coefficients bounded above and **away from zero**, plus a nonsingular inclusion, to infer asymptotic convergence from a nonincreasing Lyapunov norm. In the **rate representation**, the coefficients are instantaneous derivatives of monomials; these can vanish as a boundary is approached. Boundedness gives upper bounds, not these lower bounds.

The **concentration representation about a positive equilibrium** is different. Its mean-value derivative coefficients do have uniform positive lower bounds on a bounded nonnegative class. Indeed, if c>0 is the equilibrium and y>=0, then c+t(y-c)>=c/2 componentwise for 0<=t<=1/2. For every genuine reactant derivative of a mass-action monomial, integrating its positive derivative over that half-interval gives a positive lower bound depending only on c and the rate constants. Boundedness gives a uniform upper bound as well. Thus, when a common concentration-space polyhedral Lyapunov norm and the required nonsingularity on S are available, the sufficient global-convergence result cited by Blanchini--Giordano can apply without first assuming persistence. This is a legitimate conditional Lyapunov route. Neither existence of that common norm nor nonsingularity of the whole structural inclusion follows merely from complex balance of a particular mass-action rate vector. Their polar-norm transfer must also be used with the correct concentration/rate orientation; rate-space lower bounds cannot be silently replaced by the mean-value bounds.

### Example 8.1 — complex balance does not imply structural local stability [proved]

Take the three-cycle of complexes

    3A+2C -> 2A+3B -> 2B+3C -> 3A+2C.

Its stoichiometric matrix is

    Gamma = [[-1,-2, 3], [3,-1,-2], [-2,3,-1]].

It has rank 2, conservation law a+b+c, positive flux (1,1,1), and deficiency zero. With rates one it is complex balanced at (1,1,1). Every reactant is strictly net consumed, so the monotone dependence support assumption holds, despite the network not being pure. The rows of Yp are a permutation of the rows of Ye, so its Iosif dual is also WR with deficiency zero.

The same matrix and the same reactant-dependence pattern admit the monotone rates

    g1=a^10 c^2,   g2=a^2 b,   g3=b^10 c.

All three equal one at (1,1,1), so this point remains an equilibrium. These rates are also ordinary mass-action rates of another nonnegative-integer realization of the same Gamma:

    10A+2C -> 9A+3B,   2A+B -> 3C,   10B+C -> 3A+8B.

The derivative factor at that equilibrium is

    G = [[10,0,2], [2,1,0], [0,10,1]],

and the Jacobian Gamma G, restricted to S={a+b+c=0} in the basis (e1-e3,e2-e3), is

    J_S = [[-15,27], [24,-25]],
    det J_S = -273 < 0.

It has a positive eigenvalue. Thus the structural Jacobian class of the original complex-balanced network is not Hurwitz, and cannot possess any common nonincreasing concentration-space Lyapunov norm: a positive-eigenvalue eigenvector would make every norm grow. A GAC Lyapunov argument that requires a common norm for *all* admissible monotone derivatives therefore cannot be inferred from complex balance, even within the monotone support class and even when both Iosif orientations are WR. This does not contradict stability of the original mass-action kinetics: their exponents constrain G to a much smaller family.

**Verdict:** this literature supplies genuine structural stability transfer on its stated class; it does not supply a general entropy, boundary barrier or persistence transfer for Iosif duals. Additional Lyapunov certificates might prove individual cases. Section 3 excludes extending the conclusion to all mass-action systems merely because Gamma has a stable canonical realization.

## 9. Craciun v3 and the status of the general problem

The arXiv record of 1501.02860 lists v3 on **23 September 2026**; the PDF title page says version 3.0, 24 September 2026. Its comment identifies a major revision of the extension from low dimensions to general dimension. We read the introductory theorem statements, proof strategy, the relevant construction statements in Sections 7--8, and compared the earlier v2 outline. We did **not** check the full proof.

The v3 organization separates embedding variable-rate WR systems into toric differential inclusions (Theorem 3.1), construction of zero-separating hypersurfaces, and general-dimensional constructions. Sections 5 and 6 treat dimensions 2 and 3. Section 7 (pp. 49--76) develops subdivisions, pre-blueprints and faithful blueprints; Theorem 7.1 and the construction in Section 7.4 proceed inductively on dimension. Section 8 (pp. 77--81) turns the restricted blueprints into an exhaustive family of zero-separating hypersurfaces in arbitrary dimension. These dimensions refer to the ambient logarithmic/concentration space and its geometric subdivisions, not dim ker Gamma. The older v2 already claimed the general theorem, but its higher-dimensional extension was much shorter; v3 develops that construction much more explicitly. This is a description of the revision, not certification of its validity or a diagnosis of a particular gap in v2.

If v3's general theorem is accepted, every complex-balanced system in the proposed small-codimension class is already covered, with no transpose argument needed. A useful remaining contribution would then be a simpler proof, an effective boundary certificate, or a structural explanation of duality, not a new set of convergent complex-balanced systems.

The required AI-campaign search also found materially broader recent results. Agent L pointed to OpenAI's *Boundedness and persistence of weakly reversible mass-action systems* (25 September 2026). We independently fetched and read its Theorem 1.1 and Remark 3.5: it states global boundedness and persistence for every finite WR mass-action system with fixed positive rates, and derives complex-balanced convergence. The repository supplies a Lean development under `lean/OAI/Analysis/MassAction/`, with the top-level theorem `OAI.Problem326.global_bounded_persistent_solution` in `Main.lean`. Its catalogue scope and the top-level proof were read. The `ComparatorChallenges/MassAction.lean` file is the **target statement**, ending in a placeholder; it is not the proof. We did not run Lean or audit the imported proof dependency tree, so this report records an available formalization, not our own validation of it.

The same repository's *Uniform Permanence in Weakly Reversible Mass-Action Systems* (5 October 2026), Theorem 1.1, states a common compact convex absorbing set for each positive compatibility class. Its first two pages and theorem were read. The formalization scope file inspected expressly covers the earlier boundedness/persistence statement; the stronger class-uniform conclusion should not be called independently formalized based on that file alone. Finally `mmislan/palomar-endotactic-permanence` claims uniform permanence for arbitrary finite endotactic systems with bounded measurable rate coefficients and reports a Lean verification. Its README and claim scope were read; its proof/formalization were not audited.

These sources make a novelty claim of a new GAC *case* inappropriate. The contribution here is the structural certificate in Proposition 4.1, the explicit obstructions, and the boundary dictionary; none depends on the recent general proofs. See `references.md` for exact source links, search scope and limitations.

## 10. Computation, refutations and remaining work

`scripts/check_gac_duality.py` uses only integers and `fractions.Fraction` and imports the repository's `crn.py` without modifying it or writing its bytecode. It certifies explicit matrices, conservation/flux cones, graph counts, siphons, row matching, boundary-subnetwork ranks and the rational extinction inequalities. See the log for the current count.

`scripts/enumerate_gac_checks.py` checks the 705 labelled networks formed from three species and the six total-degree-two complexes by taking all unions of 1--3 reversible pairs (575 networks) and all simple directed 3/4-cycles (130 networks, rotation fixed but orientations retained). Exactly 702 satisfy H3. For those 702 it checks the row-permutation characterization and q=beta*+delta*. The ADLS criteria hold on both sides for 18, only the original for 177, only the dual for 0, and neither for 507; 13 have WR duals. These are **finite-family observations**, not general theorems or evidence that persistence itself fails in all the criterion failures. The first asymmetry witness has an unused third species; Section 5 removes it to give a cleaner two-species hand proof.

**Refuted:** preservation of complex counts, linkage counts, deficiency, WR, possible complex balance, endotacticity, strong endotacticity, the ADLS criterion, and persistence/global attraction under H1--H3; the inference that a persistent/complex-balanced canonical cancellation implies persistence of the uncancelled system; a q-to-rank dimension swap through transposition or the listed rank-preserving operations; and the inference from complex balance to structural local stability/common concentration Lyapunov norms over all monotone kinetic derivatives. Boundedness of classes is preserved under H1--H2 but not unconditionally.

**Proved positive statements:** exact row-graph and row-permutation WR criteria; q=beta*+delta* and the q=1 cycle consequence; the prime-species codimension-one strong-endotacticity certificate (Proposition 4.1); the siphon/T-trap and semiflow/cycle persistence certificate; the face-codimension/subnetwork-rank identity; the two exact Wiuf matrix translations; and the first-order-dual reduction.

**Conjectured:** no unproved universal mathematical assertion is advanced as a theorem. No new GAC case has been established. Whether a more restrictive, dynamically justified dual construction gives a simpler proof for a genuinely broader class remains open in this report. The possibility is not excluded by these obstructions. Independent verification and a broader novelty check of the derived characterizations would be needed before a publication claim.
