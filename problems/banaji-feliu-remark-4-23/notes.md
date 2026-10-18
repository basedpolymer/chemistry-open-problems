# Sharp equilibrium bounds for (n, n, n−1) mass-action networks: an answer to Banaji and Feliu's Remark 4.23

## The question {#sec:question .unnumbered}

Banaji and Feliu [1] study positive equilibria of mass-action networks with \(n\) species, \(m\) reactions and
stoichiometric rank \(r\), called \((n,m,r)\) networks. Their Corollary 4.21(i) states that
\((n,n,n-1)\) networks admit no more than \(n\) positive nondegenerate equilibria on any stoichiometric class.
In Remark 4.23 they write (arXiv:2409.06877v5, p. 30):

> "In the case n = 2, both bounds in Corollary 4.21 are achieved by the (2, 2, 1) network X + 2Y → 3Y, Y → X
> […]. We remark, however, that quadratic (2, 2, 1) networks do not admit multiple positive nondegenerate
> equilibria […]. We leave open whether the bound on positive nondegenerate equilibria in Corollary 4.21(i) is
> achieved for general n and, if so, the minimum source molecularity of a network needed to achieve the bound."

We answer both parts.

**Theorem A.** For every \(n\ge 2\), there is an \((n,n,n-1)\) mass-action network with source molecularity 3 which,
for suitable rate constants, has exactly \(n\) positive equilibria on some stoichiometric class, all of them
nondegenerate. The bound of Corollary 4.21(i) is therefore attained for every \(n\).

**Corollary B.** For every \(n\ge 2\), the minimum source molecularity of an \((n,n,n-1)\) network attaining the
bound \(n\) is exactly 3.

*Proof of Corollary B from Theorem A.* For \(n=2\) this is the content of Remark 4.23 itself (the network
\(X+2Y\to 3Y,\ Y\to X\) has source molecularity 3, and quadratic \((2,2,1)\) networks admit at most one positive
nondegenerate equilibrium per class). For \(n\ge 3\), Corollary 5.16 of [1] states that quadratic
\((n,n,n-1)\) networks (all sources of molecularity at most 2) admit no more than two positive nondegenerate
equilibria on each stoichiometric class, which is fewer than \(n\). Theorem A provides networks of source
molecularity 3. \(\square\)

## Notation {#sec:notation}

Species are \(X_1,\dots,X_n\) with concentrations \(x\in\mathbb R^n_{>0}\). Reaction \(R_j\) has source complex
\(a_j\in\mathbb Z^n_{\ge0}\), product complex \(b_j\in\mathbb Z^n_{\ge 0}\) and reaction vector \(\gamma_j=b_j-a_j\).
We write \(\Gamma=[\gamma_1|\cdots|\gamma_m]\in\mathbb Z^{n\times m}\) for the stoichiometric matrix and
\(S\in\mathbb Z^{m\times n}\) for the matrix whose \(j\)-th row is \(a_j\). Mass-action kinetics with rate constants
\(\kappa\in\mathbb R^m_{>0}\) gives
\[
\dot x=\Gamma\, v(x),\qquad v_j(x)=\kappa_j\,x^{a_j}=\kappa_j\prod_i x_i^{a_{ji}},
\]
with Jacobian \(J(x)=\Gamma\,D_{v(x)}\,S\,D_{1/x}\), where \(D_u\) is the diagonal matrix with diagonal \(u\).
Following [1, §2], a positive equilibrium \(x\) is *nondegenerate* if \(J(x)\) acts as a nonsingular
transformation on \(\operatorname{im}\Gamma\), i.e. \(\ker J(x)\cap\operatorname{im}\Gamma=\{0\}\).
Stoichiometric classes are the sets \((x_0+\operatorname{im}\Gamma)\cap\mathbb R^n_{\ge0}\).

For \(\lambda>0\) and \(q\in\mathbb R^n\) we write \(\lambda^q=(\lambda^{q_1},\dots,\lambda^{q_n})\), and \(\circ\)
denotes the componentwise product.

## The family of networks {#sec:family}

Fix \(n\ge3\). Define the signs
\[
\varepsilon_i=(-1)^{\lfloor (i-1)/2\rfloor}\qquad(i=1,\dots,n),\qquad\text{i.e. } \varepsilon=(+,+,-,-,+,+,-,-,\dots),
\]
and \(s=-\varepsilon_n\in\{\pm1\}\). The network \(\mathcal N_n\) has the \(n\) reactions
\[
\begin{aligned}
R_1&:\ X_1\ \longrightarrow\ X_1+\gamma_1,\\
R_{2k}&:\ 2X_{2k-1}+X_{2k}\ \longrightarrow\ 2X_{2k-1}+X_{2k}+\gamma_{2k}\qquad(2k\le n),\\
R_{2k+1}&:\ 2X_{2k}+X_{2k+1}\ \longrightarrow\ 2X_{2k}+X_{2k+1}+\gamma_{2k+1}\qquad(2k+1\le n),
\end{aligned}
\]
with reaction vectors (\(e_i\) denotes the \(i\)-th unit vector)
\[
\begin{aligned}
\gamma_1&=s\,e_1+e_n,\\
\gamma_{2k}&=c_k\,(e_{2k}-e_{2k-1}), & c_k&=(-1)^{k+1}s,\\
\gamma_{2k+1}&=d_k\,(e_{2k}+e_{2k+1}), & d_k&=(-1)^{k}s .
\end{aligned}
\]
Written out, the reactions are

| reaction | if the coefficient is \(+1\) | if the coefficient is \(-1\) |
|---|---|---|
| \(R_1\) (\(s\)) | \(X_1\to 2X_1+X_n\) | \(X_1\to X_n\) |
| \(R_{2k}\) (\(c_k\)) | \(2X_{2k-1}+X_{2k}\to X_{2k-1}+2X_{2k}\) | \(2X_{2k-1}+X_{2k}\to 3X_{2k-1}\) |
| \(R_{2k+1}\) (\(d_k\)) | \(2X_{2k}+X_{2k+1}\to 3X_{2k}+2X_{2k+1}\) | \(2X_{2k}+X_{2k+1}\to X_{2k}\) |

All products are nonnegative, the sources are \(a_1=e_1\) and \(a_j=2e_{j-1}+e_j\) (\(2\le j\le n\)), so the
source molecularity is 3 and the product molecularity is at most 5. For instance,

- \(\mathcal N_3\): \(X_1\to 2X_1+X_3,\quad 2X_1+X_2\to X_1+2X_2,\quad 2X_2+X_3\to X_2\);
- \(\mathcal N_4\): \(X_1\to X_4,\quad 2X_1+X_2\to 3X_1,\quad 2X_2+X_3\to 3X_2+2X_3,\quad 2X_3+X_4\to X_3+2X_4\).

**Lemma 1 (stoichiometry).** \(\mathcal N_n\) is an \((n,n,n-1)\) network with
\(\operatorname{im}\Gamma=\varepsilon^\perp\) (so \(\ker\Gamma^{\mathsf T}=\operatorname{span}(\varepsilon)\)) and
\(\ker\Gamma=\operatorname{span}(\mathbf 1)\), \(\mathbf 1=(1,\dots,1)\). In particular it is dynamically nontrivial and
its stoichiometric classes are the sets \(\{x\ge0:\varepsilon\cdot x=T\}\), \(T\in\mathbb R\).

*Proof.* For \(t=1,\dots,n-1\) let \(u_t=e_{t+1}-e_t\) if \(t\) is odd and \(u_t=e_t+e_{t+1}\) if \(t\) is even. Then
\(\gamma_{t+1}\) is a nonzero multiple of \(u_t\) (\(\gamma_{2k}=c_k u_{2k-1}\), \(\gamma_{2k+1}=d_k u_{2k}\)).
If \(t\) is odd, \(\varepsilon_t=\varepsilon_{t+1}\); if \(t\) is even, \(\varepsilon_{t+1}=-\varepsilon_t\). Hence
\(\varepsilon\cdot u_t=0\) for all \(t\). Also \(\varepsilon\cdot\gamma_1=s+\varepsilon_n=0\). The vectors
\(u_1,\dots,u_{n-1}\) are linearly independent: deleting the first coordinate leaves an upper bidiagonal
\((n-1)\times(n-1)\) matrix with unit diagonal (\(u_t\) has coefficient \(1\) at position \(t+1\) and is supported on
\(\{t,t+1\}\)). Therefore \(\operatorname{im}\Gamma=\operatorname{span}(u_1,\dots,u_{n-1})=\varepsilon^\perp\) and
\(\operatorname{rank}\Gamma=n-1\).

We check \(\sum_j\gamma_j=0\) coordinatewise. Coordinate 1: \(s-c_1=s-s=0\). Coordinate \(2k\): \(c_k+d_k=0\) if
\(2k<n\); if \(2k=n\), the coordinate is \(c_k+1=(-1)^{k+1}s+1\), and \(s=-\varepsilon_{2k}=-(-1)^{k-1}\) gives
\(c_k=-1\). Coordinate \(2k+1\ge3\): \(d_k-c_{k+1}=(-1)^ks-(-1)^{k+2}s=0\) if \(2k+1<n\); if \(2k+1=n\), the
coordinate is \(d_k+1\), and \(s=-\varepsilon_{2k+1}=-(-1)^k\) gives \(d_k=-1\). So \(\Gamma\mathbf 1=0\), and since
\(\operatorname{rank}\Gamma=n-1=m-1\), \(\ker\Gamma=\operatorname{span}(\mathbf 1)\). \(\square\)

**Lemma 2 (exponents).** The source matrix \(S\) is lower bidiagonal with unit diagonal, hence invertible, and
\(q:=S^{-1}\mathbf 1\) has entries
\[
q_i=\frac{1+2(-2)^{i-1}}{3}\qquad(i=1,\dots,n),\qquad q=(1,-1,3,-5,11,-21,43,\dots).
\]
These are pairwise distinct nonzero integers, with
\[
q_{2\lfloor n/2\rfloor}<\cdots<q_4<q_2<0<q_1<q_3<q_5<\cdots .
\]

*Proof.* \(Sq=\mathbf 1\) reads \(q_1=1\) and \(2q_{i-1}+q_i=1\), i.e. \(q_i-\tfrac13=-2(q_{i-1}-\tfrac13)\). For odd \(i\),
\(q_i=(1+2^i)/3>0\) increases with \(i\); for even \(i\), \(q_i=(1-2^i)/3<0\) decreases with \(i\). \(\square\)

## Equilibria lie on a monomial curve {#sec:curve}

**Lemma 3.** Fix \(\kappa\in\mathbb R^n_{>0}\) and put \(C=\exp(-S^{-1}\log\kappa)\in\mathbb R^n_{>0}\)
(componentwise logarithm and exponential). The positive equilibria of \(\mathcal N_n\) are exactly the points
\[
x(\lambda)=C\circ\lambda^{q},\qquad\lambda>0,
\]
and at \(x(\lambda)\) every reaction has rate \(v_j=\lambda\). Conversely, every \(C\in\mathbb R^n_{>0}\) arises from
the rate constants \(\kappa_j=C^{-a_j}\).

*Proof.* By Lemma 1, \(x>0\) is an equilibrium iff \(v(x)\in\ker\Gamma=\operatorname{span}(\mathbf 1)\), i.e.
\(\kappa_jx^{a_j}=\lambda\) for all \(j\) and some \(\lambda>0\). Taking logarithms,
\(S\log x=(\log\lambda)\mathbf 1-\log\kappa\), so \(\log x=(\log\lambda)\,q+\log C\). The converse follows from
\(\log\kappa=-S\log C\). \(\square\)

On the class \(\varepsilon\cdot x=T\), positive equilibria therefore correspond bijectively (the map
\(\lambda\mapsto x(\lambda)\) is injective because \(q\ne0\)) to the positive zeros of
\[
F(\lambda)=\sum_{i=1}^n\varepsilon_iC_i\,\lambda^{q_i}-T .
\]

**Lemma 4 (nondegeneracy).** The equilibrium \(x(\lambda)\) is nondegenerate if and only if \(F'(\lambda)\neq0\).

*Proof.* At \(x=x(\lambda)\) we have \(v=\lambda\mathbf 1\), so \(J=\lambda\,\Gamma S D_{1/x}\). Then \(Jy=0\) iff
\(SD_{1/x}y\in\ker\Gamma=\operatorname{span}(\mathbf 1)\) iff \(D_{1/x}y\in\operatorname{span}(S^{-1}\mathbf 1)=\operatorname{span}(q)\)
iff \(y\in\operatorname{span}(x\circ q)\). Since \(x'(\lambda)=\lambda^{-1}\,x(\lambda)\circ q\), we get
\(\ker J=\operatorname{span}(x'(\lambda))\). Because \(J(\mathbb R^n)\subseteq\operatorname{im}\Gamma=\varepsilon^\perp\),
\(\ker J\cap\operatorname{im}\Gamma=\{0\}\) iff \(x'(\lambda)\notin\varepsilon^\perp\) iff
\(\varepsilon\cdot x'(\lambda)=F'(\lambda)\neq0\). \(\square\)

## Signomials with prescribed roots {#sec:signomials}

**Lemma 5 (Descartes' rule for real exponents).** Let \(e_0<e_1<\dots<e_N\) be real and \(c\in\mathbb R^{N+1}\setminus\{0\}\).
The number of zeros of \(f(t)=\sum_k c_kt^{e_k}\) in \((0,\infty)\), counted with multiplicity, is at most the number
of sign changes of \((c_0,\dots,c_N)\) (zero entries being ignored).

This is classical (Laguerre 1883; see Pólya and Szegő [3, Part V, Ch. 1] or Jameson [4]). *Proof sketch.* Induct on
the number \(V\) of sign changes. If \(V=0\), \(f\) has no positive zero. Otherwise pick \(\theta\) strictly between
two consecutive exponents where the sign changes, and consider \(g(t)=(t^{-\theta}f(t))'=\sum_kc_k(e_k-\theta)t^{e_k-\theta-1}\).
Multiplication by \(e_k-\theta\) flips exactly the signs on one side of the chosen change, so \(g\) has \(V-1\) sign
changes, and by Rolle's theorem \(g\) has at least \(Z(f)-1\) positive zeros with multiplicity. \(\square\)

**Lemma 6.** Let \(e_0<\dots<e_N\) be real and \(0<\lambda_1<\dots<\lambda_N\). There is a nonzero
\(c\in\mathbb R^{N+1}\) with \(f(\lambda_r)=\sum_kc_k\lambda_r^{e_k}=0\) for all \(r\). Every such \(c\) has all entries
nonzero and strictly alternating in sign, and \(\lambda_1,\dots,\lambda_N\) are then all the positive zeros of \(f\),
each of them simple.

*Proof.* Existence: \(N\) homogeneous linear equations in \(N+1\) unknowns. By Lemma 5, \(f\) has at least \(N\)
positive zeros, so \((c_0,\dots,c_N)\) has at least \(N\) sign changes; with \(N+1\) entries this forces all entries
to be nonzero and to alternate. Lemma 5 again shows that \(f\) has at most \(N\) positive zeros counted with
multiplicity, so the \(\lambda_r\) are all of them and each is simple. \(\square\)

## Proof of Theorem A {#sec:proof}

The case \(n=2\) is in [1, Remark 4.23]. Let \(n\ge3\), choose \(0<\lambda_1<\dots<\lambda_n\), and list the
\(n+1\) distinct exponents \(\{q_1,\dots,q_n,0\}\) in increasing order (Lemma 2):
\[
q_{2\lfloor n/2\rfloor}<\dots<q_4<q_2<\ 0\ <q_1<q_3<\dots<q_{2\lceil n/2\rceil-1}.
\]
Lemma 6 gives a strictly alternating coefficient vector \(c\) whose signomial vanishes exactly at the \(\lambda_r\), all
simple. Replacing \(c\) by \(-c\) if necessary, the coefficient of \(\lambda^{q_1}\) is positive. Now read the signs
\(\varepsilon_i\) in the same exponent order: \(\varepsilon_{2k}=(-1)^{k-1}\) alternates along the negative exponents
\(q_2,q_4,\dots\) and \(\varepsilon_{2k+1}=(-1)^k\) alternates along the positive exponents \(q_1,q_3,\dots\), with
\(\varepsilon_1=\varepsilon_2=+1\). So the sequence of \(\varepsilon\)'s, with an entry \(-1\) inserted at the position of
exponent \(0\), alternates and starts at \(\lambda^{q_1}\) with \(+\), exactly like \(c\). Hence the coefficient of
\(\lambda^{q_i}\) in \(c\) has sign \(\varepsilon_i\), and the constant coefficient is negative. Define
\[
C_i=\varepsilon_i\,c_{(q_i)}>0,\qquad T=-c_{(0)}>0,\qquad \kappa_j=C^{-a_j},
\]
where \(c_{(e)}\) is the entry of \(c\) at exponent \(e\). By Lemma 3, the positive equilibria of \(\mathcal N_n\) with
these rate constants on the class \(\varepsilon\cdot x=T\) correspond to the positive zeros of \(F=f\), namely
\(\lambda_1,\dots,\lambda_n\). They are simple, so the \(n\) equilibria \(x(\lambda_r)\) are nondegenerate by Lemma 4,
and there are no others on that class. \(\square\)

**Remarks.**

1. If the \(\lambda_r\) are rational, then \(c\) can be chosen rational, so \(C\), \(T\), \(\kappa\) and the equilibria
   \(x(\lambda_r)=C\circ\lambda_r^q\) are rational (the exponents are integers).
2. The stoichiometric classes of \(\mathcal N_n\) are unbounded (\(\varepsilon\) has entries of both signs), as they
   must be by Corollary 4.21(ii) of [1] once \(n\ge3\).
3. \(\mathcal N_n\) is an \((n,n,n-1)\) network with trivial partition, so it also shows that the bound in
   Theorem 4.20(i) of [1] is attained, for every \(n\), in the case \(m=n\).
4. The exponents \(q_i\) grow like \(2^i\). This is harmless for the theorem but makes the rational witnesses large;
   the computer search during exploration found other families with linearly growing exponents, but without a
   uniform proof.

## A small explicit instance {#sec:example}

For \(n=3\) take \(\lambda_r=r\). The exponents are \(q=(1,-1,3)\) and
\(\lambda\,F(\lambda)=-C_3\lambda^4+C_1\lambda^2-T\lambda+C_2\) must equal \(-C_3(\lambda-1)(\lambda-2)(\lambda-3)(\lambda+6)
=-C_3(\lambda^4-25\lambda^2+60\lambda-36)\). With \(C_3=1\): \(C=(25,36,1)\), \(T=60\). The network
\[
X_1\xrightarrow{\ 1/25\ }2X_1+X_3,\qquad 2X_1+X_2\xrightarrow{\ 1/22500\ }X_1+2X_2,\qquad 2X_2+X_3\xrightarrow{\ 1/1296\ }X_2
\]
has, on the class \(x_1+x_2-x_3=60\), exactly the three positive equilibria
\[
(25,36,1),\qquad(50,18,8),\qquad(75,12,27).
\]
The Jacobian restricted to \(\operatorname{im}\Gamma\) has determinant \(7/450\), \(-1/225\) and \(1/150\) respectively,
with eigenvalues \(\{-1.053,-0.0148\}\), \(\{-0.412,+0.0108\}\) and \(\{-0.384,-0.0174\}\): two stable equilibria and a
saddle (script `example_n3.py`).

## Computer verification {#sec:verification}

The proof above is complete without computer assistance. As an independent consistency check,
`scripts/certify_family.py` builds \(\mathcal N_n\) for \(n=3,\dots,12\) and verifies in exact rational arithmetic
(python-flint): rank and kernels of \(\Gamma\); the exponent formula; the strict sign alternation of the
generalized-Vandermonde kernel vector for \(\lambda_r=r\); positivity of \(C\), \(T\), \(\kappa\); that the \(n\) points
\(x(r)\) satisfy \(\Gamma v(x)=0\) and \(\varepsilon\cdot x=T\) exactly; that the restricted Jacobian has nonzero
determinant at each of them; and that \(F\) has exactly \(n\) sign changes. The log is
`certificates/certify_family_n3-12.log` (301 checks, all PASS).

## References {#sec:references .unnumbered}

1. M. Banaji, E. Feliu, *Positive equilibria in mass action networks: geometry and bounds*, J. Math. Biol. 93
   (2026), doi:10.1007/s00285-026-02429-8; arXiv:2409.06877v5 (19 Aug 2026).
2. M. Banaji, C. Pantea, *The inheritance of nondegenerate multistationarity in chemical reaction networks*,
   SIAM J. Appl. Math. 78 (2018) 1105–1130 — cited by [1] (its reference [11]) for the \(n=2\) computations; not
   needed here.
3. G. Pólya, G. Szegő, *Problems and Theorems in Analysis II*, Springer, Part V, Chapter 1.
4. G. J. O. Jameson, *Counting zeros of generalised polynomials: Descartes' rule of signs and Laguerre's
   extensions*, Math. Gazette 90 (2006) 223–234.
