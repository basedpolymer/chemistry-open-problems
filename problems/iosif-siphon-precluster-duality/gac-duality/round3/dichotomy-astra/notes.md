# A2: rate independence and a counterexample to Conjecture F

Agent A2, **GPT-6-Astra (effort xhigh)**, 27 October 2026.
All paths below are relative to this folder. This note contains no private
correspondence. Status words distinguish proofs, exact finite computations,
conjectures, refutations, and open questions.

## 1. Setting and current conclusions

A network N has species X_1,...,X_n, reactions y_j -> y'_j, source and product
matrices Ye and Yp, and Gamma=Yp-Ye. Its dual has species k_1,...,k_r and one
reaction `row_i(Ye) -> row_i(Yp)` for each original species, with independent
positive rate constant kappa_i. Its stoichiometric matrix is Gamma^T.

(H1) means a strictly positive conservation vector for N; (H2) means a strictly
positive vector in ker(Gamma), equivalently a positive steady state for some
mass-action rate vector; (H3) means pairwise distinct rows of Gamma. Persistence
means that every positive trajectory, in every compatibility class, has no
boundary omega-limit point. ALL/SOME/NONE quantify over the positive dual rate
vectors, with all classes included for each rate vector.

| Claim | Status |
|---|---|
| Q1-7: SOME is impossible for WR N satisfying (H1)--(H3), rank two | **OPEN** |
| F: a realizable OUT facet forces NONE | **REFUTED**, by Theorem A below |
| A: the explicit four-species three-cycle below has an ALL dual and an OUT facet | **PROVED**, exact checks being finalized |
| AF: the family N_(u,R) is ALL for R>=3 and NONE for R=1,2 | **PROVED**, u>=3 |
| E: WR + (H1) implies absolute consistency of the dual, in every rank | **PROVED** |
| E+: the positive dual equilibrium set is a smooth global graph in logarithmic coordinates | **PROVED** |
| L: a single-linkage rank-two N cannot give a two-way siphon-edge flow in N* | **PROVED** |
| O: the ALL dual in Theorem A has a unique positive equilibrium in each class, but has a periodic orbit at explicit rational rates | **PROVED**, by exact source certificate and planar dynamics |

Theorem A is an obstruction to a proposed route toward the dichotomy. It is not
a counterexample to the dichotomy itself.

## 2. An ALL dual with an OUT facet

**Theorem A (proved).** Let N be the directed three-cycle

```
3X3 -> X1+X2 -> X2+2X3+3X4 -> 3X3.
```

Then N is weakly reversible, satisfies (H1)--(H3), has stoichiometric rank two
and deficiency zero. Its dual is persistent for every positive rate vector,
although every positive class has an OUT facet.

### 2.1 Exact network certificates

The matrices, with reactions in the displayed cyclic order, are

```
Ye = [[0,1,0],       Gamma = [[ 1,-1, 0],
      [0,1,1],                [ 1, 0,-1],
      [3,0,2],                [-3, 2, 1],
      [0,0,3]]                [ 0, 3,-3]].
```

The strictly positive conservation vector is `(11,1,4,1)`; its value on each
complex is 12. A positive flux vector is `(1,1,1)`. All four rows of Gamma are
distinct. The minor in rows 1,2 and columns 1,2 has determinant 1, while the
column sum is zero, so the rank is exactly two. The graph is a directed
three-cycle, with three complexes and one linkage class, giving deficiency
`3-1-2=0`.

Write the dual species as a,b,c and its rate constants as alpha,beta,gamma,delta.
The dual is

```
b       -> a          (alpha)
b+c     -> a+b        (beta)
3a+2c   -> 2b+3c      (gamma)
3c      -> 3b         (delta).
```

Its equations are exactly

```
a' = alpha*b + beta*b*c - 3*gamma*a^3*c^2,
b' = -alpha*b + 2*gamma*a^3*c^2 + 3*delta*c^3,
c' = -beta*b*c + gamma*a^3*c^2 - 3*delta*c^3.
```

Every positive compatibility class is the interior of the triangle
`P_T={a,b,c>=0: a+b+c=T}`, for T>0. It is compact and forward invariant. Positive
initial data stay positive for every finite time, since every negative term
in a coordinate equation contains that coordinate as a factor. Consequently
solutions exist for all positive times and have nonempty compact omega-limit
sets.

The edge c=0 is a siphon edge. Its minimal transverse order is 1, contributed
only by `-beta*b*c`; the inward term has order 2. Thus it is OUT. Along the edge
`b'=-alpha*b`, so flow is one-way toward e_T=(T,0,0). All points of its open edge
belong to positive classes, proving realizability.

### 2.2 The only possible boundary omega-limit point

Let Omega be an omega-limit set of an interior trajectory. Every point of Omega
has a total trajectory inside the compact invariant set Omega. At an interior
point of a=0 the a component of the vector field is strictly positive. At an
interior point of b=0 the b component is strictly positive. Such points cannot
have a backward trajectory in P_T, hence cannot belong to Omega.

At (0,T,0) and (0,0,T) the vector field is nonzero. At a vertex of a convex
polygon a nonzero vector cannot be tangent in both forward and backward time,
so these vertices are likewise excluded. At e_T the field vanishes.

An open point of c=0 has `b(t)=b(0)exp(-alpha*t)` on its total edge trajectory.
This exceeds T in finite backward time, so no such total trajectory lies in
P_T. It follows that

```
Omega intersect boundary(P_T) is contained in {e_T}.
```

### 2.3 An explicit integrating factor excludes e_T

Set `V=c*exp(-beta*b/alpha)`. For c>0,

```
(log V)' = c*[ gamma*a^3
              -(3*delta+2*beta*gamma*a^3/alpha)*c
              -(3*beta*delta/alpha)*c^2 ].                 (A.1)
```

The cancellation of the term `-beta*b` is exact. Define

```
E = 3*delta + 2*beta*gamma*T^3/alpha + 3*beta*delta/alpha,
epsilon = min(T/2, 1, gamma*T^3/(16*E)).
```

On the neighborhood `U={b+c<epsilon}` in P_T, one has a>T/2, c<epsilon<=1,
and the bracket in (A.1) is at least

```
gamma*T^3/8 - E*epsilon >= gamma*T^3/16 > 0.
```

Therefore V is strictly increasing along every interior trajectory while it
stays in U. For example, when all four rates and T equal 1, one can take
`epsilon=1/128`, with the exact bound `(log V)'>=c/16`.

For completeness, the local inequality excludes visits accumulating at e_T,
not just convergence to e_T. Suppose x(t_j)->e_T. If x(t) eventually stays in
U, monotonicity of V contradicts V(x(t_j))->0. Otherwise take the last entry
times s_j before t_j, so x(s_j) is on the relative boundary of U and the whole
segment (s_j,t_j) stays in U. Passing to a subsequence, s_j tends to infinity
and x(s_j) tends to a point q in Omega on that boundary. Monotonicity gives
`0<=V(x(s_j))<=V(x(t_j))->0`, hence c(q)=0. Section 2.2 forces q=e_T, whereas
e_T is not on the relative boundary of U. This contradiction proves
`Omega intersect boundary(P_T)=empty`.

The argument applies to every T>0 and every alpha,beta,gamma,delta>0. This
proves ALL and refutes Conjecture F. In particular, an OUT open edge may direct
trajectories toward a vertex whose center direction repels the interior; its
transverse attraction alone does not imply extinction.

## 3. A sharp infinite family

**Theorem AF (proved).** For integers u>=3 and R>=1, let

```
N_(u,R): uX3 -> X1+X2 -> X2+2X3+R X4 -> uX3.
```

These networks satisfy all the hypotheses of the question. Their duals are
NONE for R=1 or R=2 and ALL for R>=3. All have an OUT facet. Theorem A is u=R=3.

*Proof.* A positive conservation vector is
`(uR+2,u-2,R+1,u-2)`, with value `u(R+1)` on all three complexes. The positive
flux `(1,1,1)`, the determinant-1 minor in the first two rows, and the four
distinct stoichiometric rows certify the remaining hypotheses. The dual is

```
b -> a;  b+c -> a+b;  u a+2c -> 2b+u c;  R c -> R b.
```

In particular

```
b' = -alpha*b + 2*gamma*a^u*c^2 + R*delta*c^R,
c' = -beta*b*c + (u-2)*gamma*a^u*c^2 - R*delta*c^R.
```

For R>=3 the only possible boundary omega-limit point is again e_T=(T,0,0),
by the exact same edge and vertex argument in Section 2.2. With
`V=c exp(-beta*b/alpha)`, one obtains

```
(log V)' = c*[ (u-2)*gamma*a^u - R*delta*c^(R-2)
              -(2*beta*gamma/alpha)*a^u*c
              -(R*beta*delta/alpha)*c^(R-1) ].            (AF.1)
```

Let `L=(u-2)*gamma*(T/2)^u`,
`E=R*delta+2*beta*gamma*T^u/alpha+R*beta*delta/alpha`, and
`epsilon=min(T/2,1,L/(2E))`. On `b+c<epsilon`, (AF.1) is at least `Lc/2>0`.
The last-entry argument of Section 2.3 proves ALL.

For R=2 choose `T=min(1,delta/((u-2)*gamma))`. Then
`c' <= -delta*c^2`, so every positive trajectory of that class has c(t)->0.
For R=1 choose `T=min(1,delta/(2*(u-2)*gamma))`. Then
`c' <= -(delta/2)*c`, again giving extinction. These choices use
`T^u<=T` and `T^(u+1)<=T` for 0<T<=1. This proves NONE for every rate vector.

For R>=2 the unique least-order transverse term on c=0 is `-beta*b*c`; for
R=1 the additional least-order term `-delta*c` is also outward. Thus c=0 is
OUT in every case. The dichotomy within this family is structural, not a
change of rates on one fixed network. QED.

## 4. Absolute consistency of every conservative WR dual

The following statement is separate from persistence. It exploits the source
span consequence of WR-Σ directly.

**Theorem E (proved, any rank).** Suppose N is conservative and
`im(Gamma) subset im(Ye)`. Then for every positive dual rate vector kappa,
N* has at least one positive equilibrium. In particular this holds for every
conservative weakly reversible N. Neither (H2) nor (H3) is needed for this
existence conclusion.

Moreover, if Gamma has rank s, the positive equilibrium set of N* is a smooth
manifold of dimension r-s. In logarithmic coordinates it is a global graph
over a fixed linear space of dimension r-s, independent of kappa.

*Proof.* Choose a full-column-rank matrix Pi (n by s) with image im(Gamma),
and let pi_i be its rows. Choose a matrix R (r by s) with `Ye R=Pi`; the source
span assumption guarantees this is possible, and R has full column rank.
Fix a linear complement L to im(R) in R^r. Every logarithmic concentration
vector is uniquely `z=q+Rt`, where q belongs to L and t belongs to R^s.

For fixed q define positive coefficients
`c_i=kappa_i exp((Ye q)_i)` and the strictly convex function

```
F_q(t) = sum_i c_i exp(pi_i dot t).
```

Its Hessian is `Pi^T diag(c_i exp(pi_i dot t)) Pi`, positive definite. Let
`w>0` satisfy `Pi^T w=0`, as supplied by (H1). For every nonzero v some
`pi_i dot v` is positive: otherwise their strictly positively weighted sum
could be zero only if all were zero, contrary to the full rank of Pi. By
compactness of the unit sphere, `max_i pi_i dot v` has a positive lower bound
on that sphere. Consequently F_q is coercive and has exactly one minimizer
t(q). At that point

```
Pi^T (kappa_i exp((Ye(q+Rt(q)))_i))_i = 0.
```

Because Pi and Gamma have the same image, this is precisely the dual
equilibrium equation `Gamma^T (kappa_i k^(row_i Ye))_i=0`. The positive definite
Hessian and the implicit function theorem show that t(q) depends smoothly on
q. Conversely every positive equilibrium must satisfy this critical-point
equation, so uniqueness of the minimizer gives exactly this graph. Exponentiation
is a diffeomorphism, proving the manifold statement.

For WR N every product complex is some source complex, hence
`Yp=Ye Σ` and `Gamma=Ye(Σ-I)` for a suitable selector matrix Σ. This proves
the required source-span inclusion. QED.

The minimization argument is the standard convex exponential-sum/moment-map
argument underlying Birch-type theorems; no novelty is claimed for that
convex-analytic fact. The application here is to the existence and smoothness
of equilibria of the transposed, generally non-WR network. For related primary
references see Müller--Regensburger,
[arXiv:1209.6488](https://arxiv.org/abs/1209.6488) and
[arXiv:1406.6587](https://arxiv.org/abs/1406.6587) (abstracts consulted).

**Important quantifier.** The graph sections above are affine sections in
logarithmic concentration space, not stoichiometric compatibility classes.
Theorem E does not assert a positive equilibrium in every compatibility class,
or attraction, or persistence. In particular, Theorem AF with R=2 is a family
of counterexamples to any such inference.

## 5. Boundary flow for a single linkage class

**Proposition L (proved).** Suppose N has one linkage class and rank two. On a
siphon edge of N*, the nonzero active reaction vectors cannot point in both
directions. Thus its edge flow is either identically zero or one-way,
independently of the positive rate vector.

*Proof.* Factor Gamma=Pi G with Pi of full column rank two. Active dual
reactions on the edge have original species indices i whose source loads on
its zero set Z vanish; their pi_i are parallel to the edge. If two nonzero
vectors pi_i and pi_h had opposite signs, some lambda>0 would give
`pi_i+lambda*pi_h=0`, hence the original conservation vector
`e_i+lambda*e_h`. At the source of any reaction in Z both species coefficients
are zero. This nonnegative conservation law therefore has value zero there,
and hence on every complex in the linkage class. Both species have coefficient
zero in every complex, so pi_i=pi_h=0, a contradiction. QED.

More generally, in a single-linkage N any boundary equilibrium of N* has all
active nontrivial reactions absent. Indeed, its flux vector v is a nonnegative
conservation law of N. At a source indexed by a zero dual coordinate,
`v dot y_j=0` because every active dual source avoids that coordinate. This
value is constant on the linkage class. Nonnegativity then forces every
species with v_i>0 to be absent from every complex, hence to have zero
stoichiometric row. Thus equilibrium balances between distinct nonzero active
dual reactions cannot occur on the boundary in this setting.

## 6. A persistent dual with an unstable equilibrium and a periodic orbit

**Theorem O (proved).** For the dual in Theorem A, every positive rate vector
and every class have exactly one positive equilibrium. Nevertheless, at

```
(alpha,beta,gamma,delta) = (17,400,70000,1000000),
T = 201/100,
```

that equilibrium is `(a,b,c)=(1,1,1/100)` and is a hyperbolic source in its
two-dimensional class. The class contains a nonconstant periodic orbit.

*Uniqueness proof.* At equilibrium, the identity

```
b' - 2c' = (2*beta*c-alpha)*b + 9*delta*c^3
```

gives

```
0<c<alpha/(2*beta),
b(c) = 9*delta*c^3/(alpha-2*beta*c),
a(c)^3 = 3*delta*c*(alpha+beta*c)/(gamma*(alpha-2*beta*c)).
```

Both a(c) and b(c) are strictly increasing on this interval. Their sum with c
increases continuously from 0 to infinity, so each T>0 gives exactly one
positive equilibrium.

*Exact source certificate.* At the displayed point, the four reaction rates
are `(17,4,7,1)`, whose product with the dual stoichiometric matrix is zero.
With local coordinates (b,c), eliminating a by `a=T-b-c`, the Jacobian is

```
J = [[-59,3658],
     [-25,  79]],
trace(J)=20,  det(J)=86789,
eigenvalues = 10 +/- i*sqrt(86689).
```

An integer quadratic-form certificate is

```
H = [[  93655, -290457],
     [-290457,13471234]],
J^T H + H J = 3471560 I,
```

with H positive definite. Thus the quadratic form centered at this equilibrium
strictly increases in a sufficiently small punctured neighborhood. The
equilibrium cannot belong to the omega-limit set of any other interior
trajectory.

By Theorem A such an omega-limit set is a compact subset of the interior, and
by uniqueness it contains no equilibrium. The planar Poincare--Bendixson
theorem therefore makes it a nonconstant periodic orbit. For example one may
start at `(1,99/100,1/50)`, which is a different positive point in the same
class. No numerical integration is used in this existence proof.

This is not a claim of hyperbolicity or orbital stability of a particular
periodic orbit. It proves existence, and shows that the primal's weak
reversibility and deficiency zero do not transfer convergence to the unique
positive equilibrium, even when the dual is persistent for all rates.

## 7. A rejected boundary-cycle candidate

An attempted SOME construction used six species and nine reactions, arranged
as three directed three-cycles. Its exact complexes and dual are generated by
`scripts/cycle_candidate.py`; the source profiles were chosen to create a
cyclic boundary flow with rate-dependent local eigenvalue ratios.

It is **NONE**, not SOME. In the script's cyclic reaction order take

```
p=(0,1,0,1,1,1,1,1,1).
```

This is a realizable boundary vertex. For
`z=(0,-1,0,-2,0,0)`, the vector `p+Gamma^T z/4` is strictly positive, so p lies
in a positive compatibility class. In local coordinates `(k_1,k_3)` its exact
Jacobian is

```
[[-kappa_2,        0],
 [ kappa_2, -kappa_4]].
```

It is Hurwitz for every positive rate vector. Local asymptotic stability gives
interior trajectories converging to p. This exact obstruction invalidates the
candidate regardless of what its boundary cycle does in a different class.

## 8. Work remaining

The all-or-none dichotomy Q1-7 is not settled by this note. A SOME example still
requires a persistence proof for one complete rate vector in every class and
an extinction certificate for another vector. No such example is claimed.
