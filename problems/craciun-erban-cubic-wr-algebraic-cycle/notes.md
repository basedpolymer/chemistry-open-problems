# A cubic weakly reversible mass-action system with an algebraic limit cycle

## The question {#sec:question .unnumbered}

Craciun and Erban [1] introduce chemical analogues of Hilbert's 16th problem. For \(n\in\mathbb N\), \(W_n\) is the
set of two-species mass-action systems whose reactions have at most \(n\) reactant and at most \(n\) product
molecules (the set \(M_n\) of [1, Definition 2]) and whose reaction graph is weakly reversible
[1, Definition 3]. \(W^a(n)\) is the maximum number of algebraic limit cycles of systems in \(W_n\)
[1, Definition 5]. Their Table 2 records \(W^a(2)=0\), \(W^a(3)\ge0\) and \(W^a(4)\ge1\), and they write
(arXiv:2406.05057v2, §4 and §6):

> "The question of whether cubic weakly reversible two-species systems can give rise to algebraic limit cycles
> remains open, leaving us with inequality \(W^{(a)}(3)\ge 0\) in Table 2."
>
> "[…] it remains an open question whether a cubic weakly reversible system can have an algebraic limit cycle."

**Theorem.** \(W^a(3)\ge1\). Explicitly, the weakly reversible network with integer rate constants
\[
\begin{aligned}
&0\xrightarrow{434}2Y,\qquad 0\xrightarrow{448}X+Y,\qquad Y\xrightarrow{3374}0,\qquad 2Y\xrightarrow{1728}3Y,\\
&3Y\xrightarrow{194}Y,\qquad 3Y\xrightarrow{11}2X+Y,\qquad X+Y\xrightarrow{328}X+2Y,\\
&X+2Y\xrightarrow{250}3Y,\qquad 2X+Y\xrightarrow{216}3X,\qquad 3X\xrightarrow{28}X+2Y
\end{aligned}
\]
belongs to \(W_3\), and its mass-action system
\[
\begin{aligned}
\dot x&=f(x,y)=-56x^3+216x^2y-250xy^2+22y^3+448,\\
\dot y&=g(x,y)=56x^3-216x^2y+250xy^2+328xy-410y^3+1728y^2-3374y+1316
\end{aligned}
\]
has the circle \(\mathcal C=\{(x-2)^2+(y-2)^2=1\}\subset\mathbb R^2_{>0}\) as a hyperbolic, orbitally asymptotically
stable limit cycle.

## Proof {#sec:proof}

**1. The network is in \(W_3\).** Its eight complexes \(0, Y, 2Y, 3Y, X+Y, X+2Y, 2X+Y, 3X\) have molecularity at most
3, so every reaction has at most 3 reactant and 3 product molecules. The reaction graph is strongly connected:
\(0\to2Y\to3Y\to Y\to0\), \(0\to X+Y\to X+2Y\to3Y\), and \(3Y\to2X+Y\to3X\to X+2Y\to3Y\) are directed cycles whose union
contains every reaction and passes through every complex. Hence every reaction lies on a directed cycle: the
network is weakly reversible (it has a single linkage class). The mass-action equations follow from
\(\dot z=\sum_{\text{reactions}}\kappa\,(\text{product}-\text{source})\,x^{a}y^{b}\).

**2. The circle is invariant.** With \(h(x,y)=(x-2)^2+(y-2)^2-1\), direct expansion gives the polynomial identity
\[
f\,h_x+g\,h_y=s\,h,\qquad s(x,y)=-112x^2+544xy-448x-820y^2+1728y-1008 .
\]
So \(\dot h=s\,h\) along solutions, and \(\mathcal C=\{h=0\}\) is invariant. It lies in \([1,3]^2\), inside the open
positive quadrant.

**3. There is no equilibrium on \(\mathcal C\).** Parametrise \(\mathcal C\) by \(x=2+\cos\theta\), \(y=2+\sin\theta\). On
\(\mathcal C\), the field is tangent, so \(F=(f,g)=P\,(-\sin\theta,\cos\theta)\) with
\(P=-(y-2)f+(x-2)g\), and \(\dot\theta=P\). We show \(P>0\) on \(\mathcal C\). The resultant
\(\operatorname{Res}_y(P,h)\in\mathbb Q[x]\) (degree 4; the leading coefficient of \(h\) in \(y\) is 1) has no real root in
\([1,3]\) (Sturm count), so \(P\) and \(h\) have no common real zero, i.e. \(P\neq0\) on \(\mathcal C\); since
\(P(3,2)=792>0\) and \(\mathcal C\) is connected, \(P>0\) on \(\mathcal C\). Independently, a ball-arithmetic evaluation on
4096 subintervals of \([0,2\pi]\) gives \(P\ge 348\) on \(\mathcal C\).

Consequently \(\mathcal C\) is a periodic orbit, traversed counter-clockwise, of period
\(T=\int_0^{2\pi}d\theta/P(\theta)\in[0.0103320434082921\pm4.3\times10^{-17}]\).

**4. The periodic orbit is hyperbolic and stable.** For a periodic orbit \(\gamma\) of a planar \(C^1\) vector field, the
derivative of the Poincaré return map at \(\gamma\) equals \(\exp\!\left(\oint_\gamma\operatorname{div}F\,dt\right)\)
(e.g. Perko [3, §3.4, Theorem 2]). Here
\[
\oint_{\mathcal C}\operatorname{div}F\,dt=\int_0^{2\pi}\frac{\operatorname{div}F(2+\cos\theta,2+\sin\theta)}{P(\theta)}\,d\theta
\in[-3.68475192932676\pm3.3\times10^{-15}],
\]
computed by rigorous ball-arithmetic integration (Arb's `acb_calc_integrate` [4] via python-flint, which returns a
guaranteed enclosure for integrands analytic near the path; \(P\) has no zero near \([0,2\pi]\) by step 3). The
integral is negative, so the return map is a contraction at \(\gamma\): \(\mathcal C\) is a hyperbolic, orbitally
asymptotically stable periodic orbit, in particular isolated among periodic orbits, i.e. a limit cycle. It is
algebraic, being a component of \(\{h=0\}\). As a consistency check, \(\oint_{\mathcal C}s\,dt\) has the same enclosure
(the coordinate \(h\) is transversal to \(\mathcal C\) and satisfies \(\dot h=sh\), so both integrals give the
transversal Floquet exponent). \(\square\)

## Remarks {#sec:remarks}

1. **Construction.** With \(h\) fixed, the conditions "\(fh_x+gh_y=sh\)" and "the realisation uses only complexes of
   molecularity \(\le3\)" are linear in the rate constants and in the coefficients of \(s\); weak reversibility is
   enforced by requiring a strictly positive circulation on the support of the rate vector. A linear program found
   a 17-reaction network during the exploration phase (also certified, see below); a mixed-integer program
   minimising the number of reactions (SCIP) then returned the 10-reaction support above, on which the solution of
   the linear identity is unique up to a time rescaling. This gives the integer rate constants.
2. **A second example.** The 17-reaction network listed in `scripts/certify_cycle.py` (network N17; 10 complexes,
   rational rate constants) also has \(\mathcal C\) as a stable hyperbolic limit cycle, with
   \(\oint\operatorname{div}F\,dt\in[-0.636539668381478\pm2\times10^{-16}]\) and period \(1.91598\ldots\).
3. **Consequences for the table of [1].** Since \(W_3\subseteq M_3\) and \(W_3\subseteq W_4\), the example also re-proves
   \(M^a(3)\ge1\) and \(W^a(4)\ge1\) with a cubic system.
4. **Deficiency.** The network has 8 complexes, one linkage class and rank 2, hence deficiency 5; it is not covered
   by the deficiency-zero or deficiency-one theorems, consistent with the presence of a limit cycle.

## Verification {#sec:verification}

`scripts/certify_cycle.py` checks, for both networks, in exact rational arithmetic (sympy): membership in \(W_3\)
(molecularities, strong connectivity), the invariance identity with its cofactor, the resultant and the sign of
\(P\); and in rigorous ball arithmetic (python-flint/Arb): the positivity of \(P\) on 4096 subintervals, and enclosures
of the period, of \(\oint s\,dt\) and of \(\oint\operatorname{div}F\,dt\). Log: `certificates/certify_cycle.log`
(all PASS, about 3 s).

## References {#sec:references .unnumbered}

1. G. Craciun, R. Erban, *Planar chemical reaction systems with algebraic and non-algebraic limit cycles*, J. Math.
   Biol. (2025), doi:10.1007/s00285-025-02221-0; arXiv:2406.05057v2.
2. M. Feinberg, *Foundations of Chemical Reaction Network Theory*, Springer (2019) — weak reversibility, deficiency.
3. L. Perko, *Differential Equations and Dynamical Systems*, 3rd ed., Springer (2001), §3.4.
4. F. Johansson, *Numerical integration in arbitrary-precision ball arithmetic*, in: Mathematical Software – ICMS
   2018, LNCS 10931, Springer (2018) 255–263, arXiv:1802.07942.
