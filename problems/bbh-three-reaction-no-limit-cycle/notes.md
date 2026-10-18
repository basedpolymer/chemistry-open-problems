# No limit cycle in the Banaji–Boros–Hofbauer three-reaction example

## The question {#sec:question .unnumbered}

Banaji, Boros and Hofbauer [1] classify the three-reaction, planar, quadratic mass-action systems that admit an
Andronov–Hopf bifurcation (their Theorem 9, cases 9 and 10). In their Conclusions (§6, p. 22 of arXiv:2304.02303)
they write:

> "However, we leave it open whether the networks in cases 9 and 10 admitting no Andronov–Hopf bifurcations can
> have a limit cycle. For example, we do not know whether the system
> 2X → 3X + 2Y (κ₁), X + Y → 0 (κ₂), Y → X + Y (κ₃),
> ẋ = κ₁x² − κ₂xy + κ₃y, ẏ = 2κ₁x² − κ₂xy
> admits a limit cycle."

We settle the example (it has no periodic orbit at all, for any rate constants) and, with the same argument, the
whole subfamily of case-9 networks to which it belongs.

**Theorem 1.** For all \(\kappa_1,\kappa_2,\kappa_3>0\), the system
\[
\dot x=\kappa_1x^2-\kappa_2xy+\kappa_3y,\qquad \dot y=2\kappa_1x^2-\kappa_2xy
\]
has no periodic orbit in the closed positive quadrant. In particular it has no limit cycle.

**Theorem 2.** Let \(c_1,c_3,d_1\) be positive integers and consider the networks
\[
2X\xrightarrow{\kappa_1}(2+c_1)X+d_1Y,\qquad X+Y\xrightarrow{\kappa_2}0,\qquad Y\xrightarrow{\kappa_3}(1+c_3)X+Y,
\]
i.e. case 9 of [1, Theorem 9] with reaction vectors \(c=(c_1,-1,c_3)\), \(d=(d_1,-1,0)\). The mass-action system
\[
\dot x=c_1\kappa_1x^2-\kappa_2xy+c_3\kappa_3y,\qquad \dot y=d_1\kappa_1x^2-\kappa_2xy
\]
has a periodic orbit for some \(\kappa\) only if \(c_1<d_1<2c_1\), which is exactly the region where it undergoes an
Andronov–Hopf bifurcation [1, Theorem 9]. For \(d_1\le c_1\) or \(d_1\ge 2c_1\) it has no periodic orbit for any
\(\kappa>0\). The example of Theorem 1 is \(c_1=c_3=1\), \(d_1=2\), on the boundary \(d_1=2c_1\) of the Hopf region.

## Proof {#sec:proof}

**Bendixson–Dulac on the open quadrant.** Write \(F=(f,g)\) for the vector field of Theorem 2 and take the Dulac
function \(B(x,y)=y^{-b}\) on \(\mathbb R^2_{>0}\), with
\[
b=1+\frac{2\kappa_1(d_1-c_1)}{\kappa_2}.
\]
A direct computation (script `certify_dulac.py`) gives
\[
y^{\,b+1}\,\operatorname{div}(BF)=Q(x,y):=\frac{d_1\kappa_1\bigl(2(c_1-d_1)\kappa_1-\kappa_2\bigr)}{\kappa_2}\,x^2+2d_1\kappa_1\,xy-\kappa_2\,y^2 ,
\]
a quadratic form whose \(y^2\) coefficient is \(-\kappa_2<0\) and whose discriminant is
\[
\Delta=4d_1\kappa_1\bigl(\kappa_1(2c_1-d_1)-\kappa_2\bigr).
\]
(The exponent \(b\) minimises \(\Delta\) over all Dulac functions \(y^{-b}\).) If \(d_1\ge 2c_1\), then
\(\kappa_1(2c_1-d_1)-\kappa_2\le-\kappa_2<0\), so \(\Delta<0\) and \(Q\) is negative definite: \(\operatorname{div}(BF)<0\) on
the open quadrant. Since \(\mathbb R^2_{>0}\) is simply connected, the Bendixson–Dulac criterion excludes periodic
orbits contained in it.

**The boundary.** On \(\{y=0,x>0\}\) we have \(\dot y=d_1\kappa_1x^2>0\), and on \(\{x=0,y>0\}\) we have
\(\dot x=c_3\kappa_3y>0\); the origin is an equilibrium. A periodic orbit in the closed quadrant cannot contain an
equilibrium, and it cannot touch a boundary point where the field points strictly into the open quadrant (it would
have to reach that point from inside). Hence every periodic orbit would lie in the open quadrant, which is
impossible when \(d_1\ge2c_1\).

**The case \(d_1\le c_1\).** The stoichiometric matrix \(\Gamma=\begin{pmatrix}c_1&-1&c_3\\ d_1&-1&0\end{pmatrix}\) has rank 2
(its rows are not proportional since \(c_3\neq0\)), so \(\ker\Gamma\) is spanned by \(u=c\times d=(c_3,\;c_3d_1,\;d_1-c_1)\). A
positive equilibrium \(x\) requires a positive flux vector \(v(x)\in\ker\Gamma\cap\mathbb R^3_{>0}\), i.e. \(u\) or \(-u\) positive.
If \(d_1\le c_1\) this fails, so the system has no positive equilibrium. A periodic orbit in the open quadrant
would enclose an equilibrium in its interior (index theory), which would be positive; hence there is none, and the
boundary argument above excludes orbits touching the axes.

**Theorem 1** is the case \(c_1=c_3=1\), \(d_1=2\): \(b=1+2\kappa_1/\kappa_2\), \(Q=-\frac{2\kappa_1(2\kappa_1+\kappa_2)}{\kappa_2}x^2+4\kappa_1xy-\kappa_2y^2\),
\(\Delta=-8\kappa_1\kappa_2<0\). For the record, the unique positive equilibrium is
\((2\kappa_3/\kappa_2,\ 4\kappa_1\kappa_3/\kappa_2^2)\) with \(\det J=4\kappa_1\kappa_3^2/\kappa_2>0\) and \(\operatorname{tr}J=-2\kappa_3<0\): it is
globally the only recurrent object in the open quadrant. \(\square\)

**Why the example sits on the Hopf boundary.** For case 9 with \(c_3>0\), Theorem 9 of [1] gives an Andronov–Hopf
bifurcation iff \(\tfrac12(d_3/c_3+d_1/c_1)<d_2/c_2<d_1/c_1\). In the subfamily, \(d_2/c_2=1\) and \(d_3=0\), so the condition
reads \(c_1<d_1<2c_1\). For the example, \(\tfrac12(0+2)=1<1<2\) fails by equality. Theorem 2 shows that in this
subfamily, "no Andronov–Hopf bifurcation" implies "no periodic orbit".

## What remains open {#sec:open}

The general question of [1] (all case-9 and case-10 networks without an Andronov–Hopf bifurcation: case 9 with
\(d_2\ne-1\) or \(d_3\ne0\), including the branch \(c_3=0\), and all of case 10) is not settled by this note. The
same Dulac strategy with \(B=x^\alpha y^\beta\) reduces it to sign conditions on a cubic-plus-quadratic form; this is
being attacked separately.

## Verification {#sec:verification}

`scripts/certify_dulac.py` checks every identity above symbolically (sympy): the Dulac computation, the
coefficients and the discriminant in the general subfamily and in the example, the vector \(u=c\times d\), the
equilibrium and its Jacobian, and the boundary behaviour. Log: `certificates/certify_dulac.log`.

## References {#sec:references .unnumbered}

1. M. Banaji, B. Boros, J. Hofbauer, *Oscillations in three-reaction quadratic mass-action systems*, Stud. Appl.
   Math. 152 (2024) 249–278, doi:10.1111/sapm.12639, arXiv:2304.02303.
