# Hopf branches j ≥ 4 in Becker–Döring with atomization: onsets, criticality, spectra

Extension of result 8 (`../README.md`, `../paper.tex`) and of its all-M extension (`../all-M/result.md`, cited as
"all-M Lemma x").  Primary source: R. L. Pego, J. J. L. Velázquez, *Temporal oscillations in Becker–Döring
equations with atomization*, Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605
("PV (n)" = their equation (n)).

Everything marked **[certificate]** is a computer-assisted step in Arb ball arithmetic (python-flint 0.9, 192–256
bits), done by a script that prints PASS/FAIL per claim and exits non-zero on failure (Section 9).  Floating point
is used only to find starting points.

---------------------------------------------------------------------------------------------------------------

## 0. Notation and statements

Model, matrix $B=B(K,M)$ (PV (27)), $A=1+K$, mass hyperplane $H$, the exact identity $f(A\mathbb 1+v)=Bv+v_1Sv$,
first Lyapunov coefficient $\ell_1$ of the system restricted to $H$ (Kuznetsov's projection formula): as in
result 8 §§2–3 and all-M §0.  $\kappa=K\sqrt M$, $\varepsilon=M^{-1/2}$, $z=M(\varphi-1)$, $\lambda=(A-\varphi^{-1})
(\varphi-1)$, $\Lambda=M\lambda/K$.  $F(\varphi)$ is PV's characteristic function (PV (71)–(76)); eigenvalues of $B$
correspond to zeros of $F$ (PV Prop. 5.3; all-M Lemma 1.2).  $t_j$ is the $j$-th positive root of $\tan t=t$ with
$\cos t<0$, $\kappa_j^0=(\sqrt{1+t_j^2}-1)^{1/2}$ (PV (89)–(91)).

For real $M>3$ we use $\varphi^M=\exp(M\log\varphi)$ (principal branch; all points below have $\operatorname{Re}
\varphi>0$).  The *Hopf set* is
$$\mathcal H=\{(M,K,\varphi):\ M>3,\ K>0,\ \operatorname{Im}\varphi>0,\ F(\varphi;K,M)=0,\ \operatorname{Re}\lambda(\varphi)=0\}.$$
For integer $M$ its points are exactly the pairs of eigenvalues $\pm i\omega$, $\omega>0$, of $B(K,M)$
(PV Prop. 5.3, with $\varphi$ the root of modulus $\ge A^{-1/2}$).  Non-integer $M$ is a device to connect the
integers.

**Definition (branch $j$).**  $\Gamma_j$ is the connected component of $\mathcal H$ that contains the points
$(M,K,\varphi)$ with $M\to\infty$, $\kappa\to\kappa_j^0$, $z\to it_j$ — PV's $j$-th Hopf curve (their Lemma 7.3 and
Theorem 3.2 construct it for $M>M_{0,k}$ by the implicit function theorem from the limit root $(\kappa_j^0,it_j)$;
all-M Theorem 7.2 continues it in $\varepsilon$).  For integer $M$, a point of $\Gamma_j$ is written
$\kappa_j(M)$.

**Theorem A (branches $4\le j\le13$).**  Let $j\in\{4,\dots,13\}$  and let
$M_j$ be as in the table ($M_4=142$, $M_5=210$, $M_6=289$, …, $M_{13}=1135$).
1. *(onset)* $\Gamma_j$ is a $C^1$ arc which is the graph of a function $M\mapsto(K,\varphi)$ over the interval
   $(M^*_j,\infty)$, where $M^*_j\in(M_j-1,M_j)$ is enclosed in the table; $K\to\infty$ as $M\downarrow M^*_j$.
   Hence, for integer $M$, $\Gamma_j$ contains exactly one point for every $M\ge M_j$ and none for $M<M_j$:
   **branch $j$ first appears at $M=M_j$**, at very large $K$ ($K\approx14963$ for $j=4$, $M=142$;
   $K\approx4.53\times10^5$ for $j=13$, $M=1135$).
2. *(Hopf hypotheses)* For every integer $M\ge M_j$, at $\kappa=\kappa_j(M)$: $\pm i\omega$ are simple eigenvalues of
   $B$, $\operatorname{Re}d\lambda/d\kappa>0$, and $0,\pm i\omega$ are the only eigenvalues of $B$ on $i\mathbb R$
   (hypothesis H2; in particular no resonance $\pm ik\omega$, $k\ge2$, and $0$ is simple and removed on $H$).
3. *(instability of the equilibrium)* Exactly $2j-2$ zeros of $F$ (with multiplicity, $|\varphi|\ge1$) correspond to
   eigenvalues with $\operatorname{Re}\lambda>0$; in particular $B|_H$ has an eigenvalue with
   $\operatorname{Re}\lambda>0$ at $\kappa_j(M)$.
4. *(criticality)* $\ell_1<0$ at $\kappa_j(M)$ for every integer $M\ge M_j$, **except $(j,M)=(4,142)$ and
   $(13,1135)$, where $\ell_1>0$** ($\ell_1=+5.8962\ldots\times10^{-15}$ and $+2.3159\ldots\times10^{-20}$,
   normalization $q^{\mathsf H}q=1$).

Consequently, for $j\in\{4,\dots,13\}$ and every integer $M\ge M_j$ the Hopf bifurcation at $\kappa_j(M)$ of the
system restricted to $H$ is, on its two-dimensional centre manifold, **supercritical** (a unique small periodic
orbit for $\kappa>\kappa_j(M)$, of amplitude $\asymp(\kappa-\kappa_j(M))^{1/2}$, attracting within the centre
manifold), except at $(j,M)=(4,142)$ and $(13,1135)$ where it is **subcritical** (the small orbits exist for
$\kappa<\kappa_j(M)$ and repel within the centre manifold).  In the full mass hyperplane **all these periodic orbits are unstable**
(Section 6).  So the answer to PV's question for the higher branches is: the bifurcating solutions are not stable,
and on the centre manifold the bifurcation is supercritical with the two exceptions $M=142$ (branch 4) and
$M=1135$ (branch 13), the first integers of their branches, in the range covered.

**Theorem B ($K\to\infty$; every $M\le1150$).**  Let $u(M)$ be the number of eigenvalues, counted with algebraic
multiplicity, of the limit matrix $C(M)=\lim_{K\to\infty}B(K,M)/A$ in the open right half-plane, divided by 2.
For every integer $3\le M\le1150$: $C(M)$ has no eigenvalue on $i\mathbb R\setminus\{0\}$, $0$ is simple,
$u(M)=\#\{j:M_j\le M\}$ with $(M_j)_{j\le13}=(15,44,87,142,210,289,379,480,591,713,844,985,1135)$, and for every such
$M$ there is $K_0(M)$ with: for $K\ge K_0(M)$, $B(K,M)|_H$ has exactly $2u(M)$ eigenvalues with
$\operatorname{Re}\lambda>0$, all non-real.  The $j$-th pair of $C(M)$ crosses the imaginary axis at the real value
$M=M^*_j\in(M_j-1,M_j)$, transversally from left to right ($j\le13$, enclosures in Section 8).
($M_1=15$, $M_2=44$, $M_3=87$ are result 8's onsets of branches 1–3.)

---------------------------------------------------------------------------------------------------------------

**Theorem C (every branch, large $M$).** For every $j\ge1$,
$$\lim_{M\to\infty}\frac{\ell_1(q_1=1)}M=L(t_j)
=-\frac{R(R-1)(7R^2-34R+23)}{t_j(R+1)(5R-3)((R-1)^3+R+1)}<0,
\qquad R=\sqrt{1+t_j^2}.$$
Consequently the branch is supercritical on the centre manifold for every integer $M$ beyond a threshold
$\widetilde M_j$ that is not explicit. The bifurcating orbits are stable for $j=1$ and unstable for $j\ge2$.

*Proof.* Evaluate all-M's regular formula at $\varepsilon=0$, using
$e^{it}=-(1+it)/R$ and the rational parametrization
$R=(4\tau^4+1)/(4\tau^2)$, $t=(4\tau^4-1)/(4\tau^2)$,
$\kappa=(2\tau^2-1)/(2\tau)$. The symbolic certificate `limit_symbolic.py` proves that the result
is identically $L$ (S2), that every divisor is nonzero for $\tau\ge2$ (S3), and that $L<0$ there (S4).
Every Hopf root has $\tau>2.12$. The simple limit root and transversal crossing of PV Lemmas 7.1–7.3
give a local curve in $\varepsilon$; the regular formula is continuous and equals $\ell_1(q_1=1)/M$
for $\varepsilon=M^{-1/2}$. PV Theorem 3.2 supplies the spectral hypotheses for all sufficiently large $M$.
This argument uses the preliminary all-M regular formula. $\square$

## 1. What is reused

From all-M (`../all-M/result.md`, scripts copied unchanged into `scripts/`): the $\varepsilon$-regular characteristic
function $Q^\varepsilon$ and the bound for $T=A^{-M}$ (all-M §1), the parametric Krawczyk test (Lemma 2.1), the
localization lemmas 3.1, 3.3, 3.3′, the argument principle in ball arithmetic (Lemma 4.1), the region lemma 4.2′,
the monotonicity lemma 4.4 and **Theorem 4.5** (spectrum at the $j$-th Hopf point, any $j$), the regular formula for
$\widehat L=\ell_1(q_1=1)/M$ and its mean value form (all-M §5, Lemma 5.2).  From result 8 (`scripts/bdlib.py`, copied
unchanged): the closed form of $\ell_1$ at integer $M$ and the Krawczyk test in $(\kappa,z)$ at integer $M$ (result 8
§§4–5; verified there).  None of these needs $j\le3$: all-M's routines take $j$ as a parameter, and Theorem 4.5 is
stated for any $j$.

---------------------------------------------------------------------------------------------------------------

## 2. The limit $K\to\infty$

**Lemma 2.1 (splitting).**  $B=AC+E$ with $C,E$ independent of $K$:
$(Cv)_1=-(M+2)v_1-\sum_{\ell=2}^{N}v_\ell+Mv_M$, $(Cv)_\ell=v_{\ell-1}-v_\ell$ ($2\le\ell\le N$),
$(Cv)_M=v_1+v_N-v_M$; $(Ev)_1=2v_2+\sum_{\ell=3}^Nv_\ell+(1-M)v_M$, $(Ev)_\ell=v_{\ell+1}-v_\ell$, $(Ev)_M=0$.
*Proof.* Insert $K=A-1$ in PV (27): $1-K=-A+2$, $-K=-A+1$, $MK+1=AM+1-M$. $\square$

**Lemma 2.2 (characteristic polynomial of $C$).**  With
$$G(\varphi)=-\varphi^{M-2}(\varphi^3+M\varphi^2-2M\varphi+M)+(M+1)\varphi-M,$$
$\det(\nu I-C)=-G(1+\nu)/\nu$ for all $\nu$.  $G(1)=G'(1)=0$, so $R:=G/(\varphi-1)^2$ is a polynomial of degree $M-1$
with integer coefficients, and the eigenvalues of $C$ with algebraic multiplicity are $\nu=0$ and $\nu=\varphi-1$,
$\varphi$ the roots of $R$ (with multiplicity).
*Proof.* For $\varphi=1+\nu\ne0$ let $w_\ell=\varphi^{1-\ell}$ ($\ell\le N$), $w_M=(1+\varphi^{2-M})/\varphi$.  Rows
$2,\dots,M$ of $(C-\nu)w$ vanish, and row 1 equals $g(\nu)=-(M+2)-\nu-\sum_{k=1}^{M-2}\varphi^{-k}+M w_M$.  Let $W$ be
the identity with first column replaced by $w$ ($\det W=w_1=1$).  The first column of $(C-\nu)W$ is $g\,e_1$, so
$\det(C-\nu)=g\cdot\det D_{11}$, where $D_{11}$ ($C-\nu$ without row and column 1) is lower triangular with diagonal
$-\varphi$: $\det(C-\nu)=(-\varphi)^{M-1}g$.  Summing the geometric series, $\varphi^{M-1}(\varphi-1)g=G(\varphi)$.
Hence $\det(\nu-C)=(-1)^M\det(C-\nu)=-G(1+\nu)/\nu$ for $\nu\notin\{0,-1\}$, and everywhere by continuity.
(Also checked symbolically for $M=3..12$, claim K1.) $\square$

**Lemma 2.3 (where the unstable roots are).**  (a) $G$ has no root with $\operatorname{Re}\varphi\ge1$, $|\varphi|\ge3$;
(b) $G<0$ on $(1,\infty)$.  Hence every root of $R$ with $\operatorname{Re}\varphi\ge1$ lies in $|\varphi-1|<4$, and the
roots with $\operatorname{Re}\varphi>1$ are non-real (they come in conjugate pairs).
*Proof.* (a) For $r=|\varphi|\ge3$, $\operatorname{Re}\varphi\ge1$: $|\varphi+M|\ge M+1$, so
$|\varphi^3+M\varphi^2-2M\varphi+M|\ge(M+1)r^2-2Mr-M\ge r^2(2M/9+1)$, and
$|G|\ge r^M(2M/9+1)-(2M+1)r>0$ since $r^{M-1}(2M/9+1)\ge9(2M/9+1)>2M+1$.  (b) $-G(x)=h(x):=x^{M-2}(x^3+M(x-1)^2)-(M+1)x+M$
has $h(1)=h'(1)=0$ and $h''>0$ on $[1,\infty)$ (every term of $h''$ is $\ge0$ and $x^{M-2}(6x+2M)>0$). $\square$

**Lemma 2.4 (the limit of $F$).**  $\lim_{A\to\infty}F(\varphi)/A^3=-\varphi^3G(\varphi)$, locally uniformly in
$\varphi$ with $|\varphi|>0$, for $M\ge4$ (the $(A\varphi)^{-M}$ terms are $O(A^{3-M})$) [sympy, K1].  The onset chart
of Section 3 makes this uniform with derivatives.

**Proposition 2.5 (count at large $K$).**  If $R$ has no root on $\operatorname{Re}\varphi=1$ and $2u$ roots in
$\operatorname{Re}\varphi>1$, then for all $K$ large enough $B(K,M)|_H$ has exactly $2u$ eigenvalues (with
multiplicity) in $\operatorname{Re}\lambda>0$, all non-real.
*Proof.* $B/A=C+E/A$; the eigenvalues depend continuously on $1/A\in[0,1)$ (with multiplicity).  At $1/A=0$ they
are $0$ (simple) and $\varphi-1$ for the roots of $R$, none on $i\mathbb R$ except 0.  The simple eigenvalue 0 of
$B$ persists for all $K$ (PV Lemma 3.1).  For $1/A$ small the other eigenvalues stay in the same open half-plane
as their limits; those near the $2u$ non-real unstable ones are non-real. $\square$

**[certificate] (kinf.py).**  K2: for every integer $3\le M\le1150$ the winding number of $G$ along the boundary of
$\{\operatorname{Re}\varphi>1,|\varphi-1|<4\}$, indented to the left around the double zero $\varphi=1$ by a half-disc
of radius $0.02/M^2$ (on which $R\ne0$ is certified from $R(1)$ and an enclosure of $R'$ on the disc), is computed
by all-M Lemma 4.1 (each edge subdivided until the mean-value enclosure of $G$ excludes 0): it equals
$2+2u(M)$, so $R$ has no root on $\operatorname{Re}\varphi=1$ and $2u(M)$ roots in $\operatorname{Re}\varphi>1$.
K3: $u(3)=0$, $u$ is nondecreasing with unit jumps exactly at the $M_j$ listed in Theorem B.
K4: for $j\le13$, the system $G(1+i\eta;M)=0$ (real $M$) has a unique zero $(M^*_j,\eta_j)$ in a Krawczyk box, with
$M_j-1<M^*_j<M_j$, and $d\operatorname{Re}\varphi/dM=\operatorname{Re}(-\partial_MG/G')>0$ there.

---------------------------------------------------------------------------------------------------------------

## 3. The onset chart (regular at $K=\infty$)

With $s=1/A=1/(1+K)\in(0,1)$, dividing PV (71)–(75) by $A^3$ gives exactly (a polynomial identity in $s$,
$\varphi$, $M$, $\varphi^{\pm M}$, $s^M$; claim K5 of `kinf.py`):
$$\frac F{A^3}=f_1\big(H_1-s^{M-1}\varphi^{-M}H_2\big)-\big(G_1-\varphi^MG_2\big)f_2,$$
$f_1=\varphi^M+s(\varphi-1)$, $f_2=s-\varphi+s^{M-1}\varphi^{1-M}$,
$G_1=\varphi(s-\varphi)\big((1-s)\varphi+(M(1-s)+s)(\varphi-1)\big)$,
$G_2=(s-\varphi)\big(M(\varphi-1)+(1-s)\varphi+\varphi(\varphi-1)\big)$,
$H_1=(\varphi-1)\big(s(1-s)+(M(1-s)+s)(s-\varphi)\big)$,
$H_2=(\varphi-1)\big((s+M\varphi)(s-\varphi)+s(1-s)\varphi\big)$,
and $\lambda/A=(1-s/\varphi)(\varphi-1)$.  The *onset system*
$$\mathcal O(s,M,u,v)=\big(\operatorname{Re}F/A^3,\ \operatorname{Im}F/A^3,\ \operatorname{Re}\lambda/A\big),\qquad
\varphi=1+u+iv,$$
has the same zeros as the Hopf system for $s\in(0,1)$ ($A>0$), is $C^1$ on $[0,1)\times\{\operatorname{Re}\varphi>0\}$
for $M\ge3$ ($s^{M-1}$ is $C^1$ up to $s=0$), and at $s=0$ reduces to $\varphi\,G(\varphi)\cdot(-\varphi^2)=0$,
$u=0$: the zero $(M^*_j,0,\eta_j)$ of K4.  For $s$-balls containing 0, $s^{p}$ ($p=M\pm\dots\ge2$) and its
derivatives are enclosed by monotonicity ($t^p\in[0,b^{p_{\min}}]$ for $t\in[0,b]$, etc.;
`onset_sweep.s_pow`).

**[certificate] (onset_sweep.py).**  Part A ($s$ as parameter, $s\in[0,2\cdot10^{-4}]$): on every $s$-box, the
parametric Krawczyk test (all-M Lemma 2.1 with parameter $s$) proves a unique zero $(M,u,v)(s)$ in a box, and
$dM/ds>0$, $v>0$ on it.  Part B ($M$ as parameter, from inside the last $A$-box up to $M_{HI}>\varepsilon_c^{-2}$):
on every $M$-box, a unique zero $(s,u,v)(M)$ in a box, with $s>0$, $v>0$ (the curve is a graph over $M$).
Junction A/B: the first $B$-box lies inside the last $A$-box.  Junction with the $\varepsilon$-chart: an $M$-sub-box
whose image $(\varepsilon,\kappa,z)=(M^{-1/2},(1/s-1)M^{1/2},M(\varphi-1))$ lies inside the last Krawczyk box of the
$\varepsilon$-sweep of Section 4 (for its $\varepsilon$-range) and inside a box of part B.

---------------------------------------------------------------------------------------------------------------

## 4. The branch as one arc; proof of Theorem A.1

Three certified pieces, each the graph of a $C^1$ function over an interval with uniqueness in its boxes:
(I) the $\varepsilon$-chart, $\varepsilon\in[0,\varepsilon_c]$ (`l1_sweep.py`, all-M Theorem 7.2's method),
which near $\varepsilon=0$ is PV's curve from $(\kappa_j^0,it_j)$; (II) part B of the onset chart, a graph over
$M\in[M_{B0},M_{HI}]$; (III) part A, a graph over $s\in[0,2\cdot10^{-4}]$ with $M$ increasing in $s$.  The junctions
of Section 3 give a common point of (I) and (II) and a common sub-arc of (II) and (III); by uniqueness in the boxes
and since each piece is a 1-manifold, the union $\Gamma$ is one $C^1$ arc on which $M$ is a strictly monotone
parameter, $M\in(M^*_j,\infty)$ ($M\to\infty$ at the $\varepsilon=0$ end, $M\to M^*_j$ and $K\to\infty$ at the $s=0$
end).  $\Gamma$ is open in $\mathcal H$ (each of its points has a neighbourhood in which $\mathcal H$ is $\Gamma$: the
point lies in the interior of the union of at most two consecutive boxes, where the zero set is the graph) and
closed in $\mathcal H$ (its closure adds only the two ends, which are not in $\mathcal H$: $M=\infty$, resp.
$K=\infty$).  So $\Gamma$ is a connected component of $\mathcal H$; it contains PV's points, hence $\Gamma=\Gamma_j$.
For integer $M$: $M\ge M_j=\lceil M^*_j\rceil$ gives exactly one point, $M<M_j$ none. $\square$

---------------------------------------------------------------------------------------------------------------

## 5. Proof of Theorem A.2–A.4

*Integers $M\ge M_c$* ($M_c^{-1/2}\le\varepsilon_c$): **[certificate]** `l1_sweep.py` on $[0,\varepsilon_c]$ (all-M
Lemma 2.1 and the checks of all-M §2: simplicity, transversality, $\operatorname{Im}z>0$, $f(\varphi)\ne0$; and
$\widehat L<0$ on every box by all-M Lemma 5.2), re-proved line by line by `verify_sweep.py`; `h2cover.py`: all-M
Theorem 4.5 on merged parameter boxes (every Krawczyk box re-proved first; for a group of consecutive boxes the
checks are made for all $(\varepsilon,\kappa)$ in the union of their $\varepsilon$-ranges times the hull of their
$\kappa$-ranges, which contains $(\varepsilon,\kappa_j(\varepsilon))$; counts $N(P)=2j+2$, $N(P_-)=4$).  Contour
parameters: $R_0=2\pi j+7$ (rounded), $Y=R_0+0.6$, level values $(c_L,c_R)=(-0.1,0.25)$, and where these fail
(near $\varepsilon=0$ for $j\ge5$, where the neighbouring zeros are at distance $\approx\log(t_{j\pm1}/t_j)$ from the
axis) the narrow pair $c_L=-\tfrac12\log(t_{j+1}/t_j)+0.02$, $c_R=\tfrac12\log(t_j/t_{j-1})-0.02$.

*Integers $M_j\le M<M_c$*: **[certificate]** `pointwise.py` at each $M$: (a) result 8's Krawczyk test in
$(\kappa,z)$ at integer $M$, simplicity, transversality, $\omega>0$; (b) result 8's closed form of $\ell_1$ in both
normalizations, sign certified; (c) all-M's Krawczyk test at the ball $\varepsilon=M^{-1/2}$, whose box contains the
box of (a) (same zero), $\operatorname{Im}z>0$, $|AE+\varepsilon^2z|>0$; (d) all-M Theorem 4.5 at that point
(localization by Lemma 3.3′ part (II′) on $[R_0\varepsilon^2,a_{big}]$ and the tail (III), so no condition $K<1$;
Lemma 4.2′ needs $\varepsilon^2R_0\le2/5$, true here).  `pointwise.py --verify` re-proves every line from the stored
box.  `link_check.py`: the box of (a), mapped to the onset chart, lies in a re-proved box of `onset_sweep.py`, so the
point is $\Gamma_j$'s point at that $M$ (and by Theorem A.1 the only one).

In Theorem 4.5 the spurious zeros $\varphi=1/A,\pm A^{-1/2}$ have $\operatorname{Re}z\le-1/(2\varepsilon(\kappa^{-1}+
\varepsilon))$, which is $<-5$ for $\varepsilon\le0.09$, $\kappa\ge3$, far left of $x_{\min}\ge-0.13$; this is the only
place where all-M used $\varepsilon\le0.0789$, and $M\ge142$ gives $\varepsilon\le0.0840$.

**Separate cross-checks.**  `check_overlap.py`: on $M_c\le M\le3M_c$ (branch 4: $229\le M\le1000$) result
8's closed form (integer $M$) lies inside the $\varepsilon$-chart enclosure of $\widehat L$ and its Hopf box inside
the $\varepsilon$-chart Krawczyk box.  `xcheck_sparse.py`: result 8's verifier's $O(M)$ sparse code (reduced
coordinates, 320-bit midpoint iteration, no shared code) re-locates the Hopf point at $M=142,143$ ($j=4$) and
$M=1135,1136$ ($j=13$) and reproduces $\ell_1$ to $6\times10^{-17}$ relative, including the positive signs at $M=142$
and $M=1135$ (logs `xcheck_j4.log`, `xcheck_j13.log`).

---------------------------------------------------------------------------------------------------------------

## 6. Instability of the bifurcating orbits ($j\ge2$)

By Theorem A.2 the Hopf theorem with centre-manifold reduction applies at $\kappa_j(M)$ (system on $H$): the
two-dimensional centre manifold carries the bifurcating orbits $\Gamma_\kappa$, and their direction and stability
within it are given by the sign of $\ell_1$.  By Theorem A.3, $B|_H$ has an eigenvalue $\lambda_u$ with
$\operatorname{Re}\lambda_u>0$ at $\kappa_j(M)$.  As $\kappa\to\kappa_j(M)$ along the family, $\Gamma_\kappa$
shrinks to the equilibrium and its period tends to $2\pi/\omega$, so its monodromy matrix (in $H$) tends to
$\exp(2\pi B|_H/\omega)$, which has the eigenvalue $e^{2\pi\lambda_u/\omega}$ of modulus $>1$.  By continuity of the
spectrum, $\Gamma_\kappa$ has a Floquet multiplier outside the unit circle for $\kappa$ close to $\kappa_j(M)$: the
orbits are unstable (all-M Theorem 2.3, same argument).

---------------------------------------------------------------------------------------------------------------

## 7. Behaviour near the onsets (observations, floating point / midpoint balls)

Along $\Gamma_j$ near $M^*_j$ (where $\kappa\to\infty$), the computed values follow
$\widehat L\,\kappa^3\approx c_j+a_j\kappa$ with $c_j\approx-5$ and a small $a_j>0$
($a_4\approx7.9\times10^{-5}$, $a_5\approx3.4\times10^{-5}$, $a_7\approx1.1\times10^{-5}$), while
$M-M^*_j\approx C_j/\kappa$ ($C_4\approx3.3\times10^3$, $C_5\approx5.9\times10^3$, $C_7\approx1.45\times10^4$).
So the sign of $\ell_1$ changes once on each branch, at $M^*_j+w_j$ with a narrow window $w_j\approx C_ja_j/|c_j|$
($w_4\approx0.046$, $w_5\approx0.035$, $w_7\approx0.027$): an integer onset is subcritical exactly when
$M_j-M^*_j<w_j$.  This explains $M=15$ (branch 1, result 8) and $M=142$ (branch 4: $M_4-M^*_4=0.0182$), and
predicted, from the fractional parts of $M^*_j$, the cases $j=13$ ($M_{13}-M^*_{13}=0.0049$) and $j=40$
($M_{40}-M^*_{40}\approx0.0004$, floating point) among $j\le40$.  The prediction for $j=13$ is confirmed by Theorem A:
$\ell_1(q^\mathsf Hq=1)=+2.316\times10^{-20}>0$ at $(13,1135)$, and $\ell_1<0$ from $M=1136$ on (the sparse cross-check
`xcheck_sparse.py` reproduces both signs at $M=1135,1136$).  The case $j=40$ is not computed.  The empirical rule
$w_j\approx0.18/j$ ($jw_j=0.184,0.175,0.189$ for $j=4,5,7$) makes $\sum_jw_j$ diverge, so infinitely many subcritical
onsets are to be expected if the fractional parts of $M^*_j$ are spread out; this is a heuristic, not a claim.

---------------------------------------------------------------------------------------------------------------

## 8. Results (filled from the certificate logs)

| $j$ | $M^*_j$ (approx.) | $M_j$ | $\ell_1(M_j)$ (unit norm, approx.) | $\varepsilon_c$ | $M_c$ | boxes / H2 groups | integer points |
|---|---|---|---|---|---|---|---|
| 4 | 141.9817605587 | 142 | +5.896225e-15 | 0.06610 | 229 | 1028 / 252 | 142–228 (87) |
| 5 | 209.6289461098 | 210 | -5.627984e-12 | 0.05456 | 336 | 1032 / 268 | 210–335 (126) |
| 6 | 288.7385783342 | 289 | -4.811125e-13 | 0.04653 | 462 | 1084 / 272 | 289–461 (173) |
| 7 | 378.8658344869 | 379 | -1.821088e-14 | 0.04059 | 607 | 1180 / 296 | 379–606 (228) |
| 8 | 479.6599235850 | 480 | -1.233926e-13 | 0.03609 | 768 | 1225 / 319 | 480–767 (288) |
| 9 | 590.8342776983 | 591 | -5.413867e-15 | 0.03252 | 946 | 1277 / 341 | 591–945 (355) |
| 10 | 712.1483888900 | 713 | -3.661047e-13 | 0.02962 | 1140 | 1340 / 362 | 713–1139 (427) |
| 11 | 843.3961316756 | 844 | -6.227954e-14 | 0.02721 | 1351 | 1406 / 375 | 844–1350 (507) |
| 12 | 984.3978885172 | 985 | -3.164222e-14 | 0.02519 | 1576 | 1479 / 399 | 985–1575 (591) |
| 13 | 1134.9950223691 | 1135 | +2.315987e-20 | 0.02346 | 1817 | 1557 / 449 | 1135–1816 (682) |

All displayed decimals are approximations to certified values. Full onset enclosures and signs are in
`certificates/logs/kinf.log` and the pointwise verification logs. Theorem A covers branches 4–13; the
$K=\infty$ onset computation (Theorem B) covers branches 1–13.  The parameters of each run are in
`certificates/logs/run_branch_params.txt`.

---------------------------------------------------------------------------------------------------------------

## 9. Scripts

(see README)

---------------------------------------------------------------------------------------------------------------

## 10. What is proved, what is not

**Proved (computer-assisted, this work):** Theorems A (for $j=4,\dots,13$), B and C.
**Taken from result 8 and all-M (verified, resp. not yet verified):** the closed form of $\ell_1$ (result 8, verified),
all-M's lemmas and Theorem 4.5 (all-M is not yet verified; its scripts are reused unchanged).
**Standard results used:** Hopf theorem with centre-manifold reduction, Kuznetsov's formula, Krawczyk–Moore, Arb.
**Not proved:**
* branches $j\ge14$ for every $M$: Theorem C covers only sufficiently large $M$, with a threshold that is not
  explicit; no argument uniform in $j$ with explicit constants;
* that the onset $M_j$ is the first $M$ at which $B$ has $j$ crossings of the imaginary axis for *some* $K$: Theorem A.1
  concerns the branch $\Gamma_j$ (PV's continuation); other components of $\mathcal H$ (crossings that do not belong to
  $\Gamma_1,\dots,\Gamma_j$) are not excluded for $M<M_j$, nor in general; Theorem B gives the count only for
  $K\ge K_0(M)$, with $K_0(M)$ not explicit;
* that the $2j-2$ unstable zeros are simple (only "at least one unstable eigenvalue" is used);
* global statements about the periodic orbits away from the bifurcation.
