# Hopf criticality in Becker–Döring with atomization for every M: proofs, constants, certificates

Extension of result 8 (`../README.md`, `../paper.tex`, `../verification.md`) to all $M$.
Primary source: R. L. Pego, J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605 (v2, read in full). We write
"PV (n)" for equation (n) and "PV Lemma/Prop. x" for their numbered statements; numbering is that of the arXiv v2 PDF.

Everything below marked **[certificate]** is a computer-assisted step in Arb ball arithmetic (python-flint 0.9,
192 bits); every such step is re-done by a script that prints PASS/FAIL (Section 9). Steps marked **[result 8]** or
**[verifier of result 8]** are taken from the independently verified result 8 and are not redone here, except
where said.

---------------------------------------------------------------------------------------------------------------

## 0. Statements

Notation (as in PV and result 8): $n_\ell$, $1\le\ell\le M$, $N=M-1$, $J_\ell=n_\ell n_1-n_{\ell+1}$,
$$\dot n_\ell=J_{\ell-1}-J_\ell\ (2\le\ell\le N),\qquad \dot n_M=J_N-Kn_M,\qquad \dot n_1=-J_1-\textstyle\sum_{\ell=1}^NJ_\ell+MKn_M .$$
$A=1+K$; the constant state $n\equiv A$ is an equilibrium; $B=B(K,M)$ is the Jacobian there (PV (27));
$H=\{v:\sum_\ell \ell v_\ell=0\}$ is the invariant mass hyperplane (in the coordinate $v=n-A\mathbb 1$);
$\kappa=K\sqrt M$ is the bifurcation parameter; $\varepsilon=M^{-1/2}$. The exact identity
$f(A\mathbb 1+v)=Bv+v_1Sv$ (result 8, §2) defines the linear map $S$; $\ell_1$ denotes the first Lyapunov
coefficient of the system restricted to $H$ (Kuznetsov's projection formula, result 8 §3). Its sign does not depend
on the normalization of the eigenvector; we use $q_1=1$ and the scale-free quantity
$$\widehat L:=\ell_1(q_1=1)/M=\varepsilon^2\,\ell_1(q_1=1).$$

**Theorem 1 (first branch, every $M$).**
1. *(verifier of result 8, exact algebra)* For $3\le M\le14$ the constant equilibrium is linearly stable on $H$ for
   every $K>0$ and no eigenvalue of $B|_H$ reaches $i\mathbb R$: there is no Hopf bifurcation.
2. *(result 8, Theorem A, and its verifier)* For $M=15$ there is a Hopf point $\kappa_1(15)$ with $\pm i\omega$ simple,
   transversal crossing, hypothesis H2 (no other eigenvalue on $i\mathbb R$), and $\ell_1>0$: **subcritical**.
3. For every integer $M\ge16$ there is $\kappa_1(M)>0$ such that at $\kappa=\kappa_1(M)$:
   * $B$ has the simple eigenvalues $\pm i\omega$, $\omega>0$, and $\operatorname{Re}d\lambda/d\kappa>0$;
   * every other eigenvalue of $B|_H$ has negative real part (in particular H2 holds);
   * $\ell_1<0$.

   Hence the Hopf bifurcation of the restricted system at $\kappa_1(M)$ is **supercritical**: for $\kappa>\kappa_1(M)$
   close to $\kappa_1(M)$ there is a unique small periodic orbit near the constant equilibrium in the mass
   hyperplane, of amplitude $\asymp(\kappa-\kappa_1(M))^{1/2}$, and it is orbitally asymptotically stable within the
   hyperplane; for $\kappa<\kappa_1(M)$ close to $\kappa_1(M)$ there is none.

   *Sources:* $16\le M\le160$: result 8 (Theorem A: Hopf box, simplicity, transversality, $\ell_1<0$) and its
   verifier (H2 and the count "no unstable pair", `acb_mat.eig`). $M\ge161$: **this work**, certificate
   `certificates/l1h2_j1.csv` (Theorem 7.1 below), which is independent of result 8; for $161\le M\le10^5$ it concerns
   the same Hopf points as result 8 (all 99 840 points checked, Section 8).
4. *(asymptotics)* $\widehat L(\varepsilon)$ is given by an expression that is continuous on $[0,\varepsilon_0]$
   (Section 6); the limit $\lim_{M\to\infty}\ell_1(q_1=1)/M$ exists, and the certified enclosures give
   $\widehat L\in[-0.0093670,-0.0092545]$ for every $M\ge10^{10}$ and
   $\widehat L(0)\in[-0.00931144,-0.00931137]$, which contains result 8's formal limit value $L_1=-0.0093114$.

**Theorem 2 (second and third branches).** Let $j\in\{2,3\}$ and $\kappa_j(M)$ be the $j$-th Hopf point
(result 8 for $M\le10^5$; Krawczyk boxes of this work for $M\ge161$).
1. $\ell_1<0$ at $\kappa_j(M)$ for every $M\ge44$ ($j=2$) and every $M\ge87$ ($j=3$): result 8 for $M\le10^5$,
   this work (`certificates/l1_j2.csv`, `l1_j3.csv`) for $M\ge161$ [ranges: see Section 7].
2. Spectral hypothesis: for $j=2$ and every $M\ge161$, and for $j=3$ and every $M$ in the range of Section 7, at
   $\kappa_j(M)$ the only eigenvalues of $B|_H$ on $i\mathbb R$ are $\pm i\omega$, and exactly $2(j-1)$ zeros of $F$
   (counted with multiplicity) correspond to eigenvalues with $\operatorname{Re}\lambda>0$. For $M\le160$ the same
   was certified by the verifier of result 8.
3. Meaning of "supercritical" for $j\ge2$: the restriction of the flow to the 2-dimensional centre manifold at
   $\kappa_j(M)$ undergoes a supercritical Hopf bifurcation (the small periodic orbits exist for $\kappa>\kappa_j(M)$
   and attract within the centre manifold). In the full mass hyperplane these orbits are **unstable**: for $\kappa$
   near $\kappa_j(M)$ the equilibrium has $j-1$ unstable conjugate pairs, which persist as $2(j-1)$ Floquet
   multipliers of the small orbits outside the unit circle. They are saddle-type cycles; the question "are the
   bifurcating solutions stable" has the answer *no* for $j\ge2$, and *yes* (for $M\ge16$) for $j=1$.

What is **not** proved is listed in Section 10.

---------------------------------------------------------------------------------------------------------------

## 1. The $\varepsilon$-regular characteristic function

PV reduce the eigenvalue problem $BV=\lambda V$ to the zeros of
$F(\varphi)=f_1g_2-g_1f_2$ (PV (71)–(75)), with $\lambda=(A-\varphi^{-1})(\varphi-1)$ (PV (61)) and the sorted
representation (PV (76)) $F=-P_1+\varphi^MP_2+A^{-M}R_1+(A\varphi)^{-M}R_2$, where (PV (102), (103), (140), (141))
$$P_1=K[(A\varphi-1)^2\varphi^2-(\varphi-1)^2]+(KM+1)S(\varphi),\quad S=(\varphi-1)(A\varphi-1)(A\varphi^2-1),$$
$$P_2=\varphi(A\varphi-1)^3-A^2(\varphi-1)^2+M(A\varphi-1)A^2(\varphi-1)^2,$$
$$R_1=-A^2[(A\varphi-1)^2\varphi^2-(\varphi-1)^2],\qquad R_2=A(\varphi-1)^3+A^2\varphi(A\varphi-1)^2(M(\varphi-1)+\varphi).$$
*Check F1* (exact, symbolic $M$): these closed forms reproduce $f_1g_2-g_1f_2$.

**Scaled variables.** $\varepsilon=M^{-1/2}$, $K=\kappa\varepsilon$, $\varphi=1+\varepsilon^2z$. Exact identities
(*check F2*):
$$A\varphi-1=\varepsilon g,\ g=\kappa+A\varepsilon z;\qquad A\varphi^2-1=\varepsilon h,\ h=g\varphi+\varepsilon z;$$
$$\lambda=\varepsilon^3gz/\varphi=:\varepsilon^3\hat\Lambda,\qquad \Lambda:=M\lambda/K=z+\frac{\varepsilon z^2}{\kappa\varphi};$$
$$E:=\varphi^M=\exp\!\big(z\,\mathsf L(\varepsilon^2z)\big),\quad \mathsf L(u)=\frac{\log(1+u)}u;\qquad T:=A^{-M}=\exp\!\big(-\tfrac\kappa\varepsilon\mathsf L(\kappa\varepsilon)\big);$$
$$Q^\varepsilon(z;\kappa):=K^{-3}F(1+z/M)=E\,\tilde P_2-\tilde P_1+\frac T\varepsilon\,\frac W{\kappa^3},$$
$$\tilde P_2=\frac{\varphi g^3-\varepsilon A^2z^2+A^2z^2g}{\kappa^3},\quad \tilde P_1=\frac{\kappa(g^2\varphi^2-\varepsilon^2z^2)+(\kappa+\varepsilon)zgh}{\kappa^3},\quad
W=-A^2(g^2\varphi^2-\varepsilon^2z^2)+\frac{A\varepsilon^4z^3+A^2\varphi g^2(z+\varphi)}{E}.$$
$\tilde P_1,\tilde P_2,W$ are polynomials in $(z,\kappa,\varepsilon)$ (and $1/E$); $E$ is analytic in
$(z,\varepsilon)$ for $|\varepsilon^2z|<1$; at $\varepsilon=0$ one gets PV's limit function
$Q(z;\kappa)=e^z(1+z^2/\kappa^2)-(1+z)$ (PV (87)). The only non-analytic term is $T/\varepsilon$, which is
exponentially small:

**Lemma 1.1 (bound for $T$).** For $\varepsilon>0$, $\kappa>0$, on a box $\varepsilon\in[\varepsilon_{lo},\varepsilon_{hi}]$,
$\kappa\in[\kappa_{lo},\kappa_{hi}]$, $\sup T/\varepsilon^k$ is bounded by the smaller of
(a) $\exp(-\kappa_{lo}/\varepsilon_{hi}+\kappa_{hi}^2/2)/\varepsilon_{hi}^k$, valid if $\varepsilon_{hi}<\kappa_{lo}/k$
(from $\log(1+x)\ge x-x^2/2$ for $x\ge0$, and since $\varepsilon\mapsto e^{-c/\varepsilon}\varepsilon^{-k}$ is increasing on
$(0,c/k)$), and (b) $T(\varepsilon_{hi},\kappa_{lo})/\varepsilon_{lo}^k$ if $\varepsilon_{lo}>0$ (since
$\partial_\varepsilon[\log(1+\kappa\varepsilon)/\varepsilon^2]=[x/(1+x)-2\log(1+x)]/\varepsilon^3<0$, $x=\kappa\varepsilon$,
$T$ is increasing in $\varepsilon$ and decreasing in $\kappa$) (`allm_core.T_over_eps_bound`; (b) matters when
$K=\kappa\varepsilon$ is not small, e.g. on the third branch near $M=180$, where $\kappa_3\approx26$). Derivatives:
$\partial_\kappa(T/\varepsilon)=-T/(\varepsilon^2A)$, $|\partial_\varepsilon(T/\varepsilon)|\le(3\kappa+\varepsilon)T/\varepsilon^3$.

In all certificates $T$-dependent terms enter as disks of radius given by Lemma 1.1, so every enclosure is valid for
all $\varepsilon$ in a box $[\varepsilon_{lo},\varepsilon_{hi}]$ **including $\varepsilon_{lo}=0$** (the value at
$\varepsilon=0$ being the limit $M\to\infty$). Formulas are evaluated on real balls for $\varepsilon$; where a ball
dips below 0 by the 30-bit radius rounding of Arb, the analytic parts are still enclosed and $T$ is bounded from the
upper endpoint only. Functions $\mathsf L$, $\mathsf M_1=\mathsf L'$, $\mathrm{Ex}(b)=(e^b-1)/b$,
$\mathrm{Ex}_1=\mathrm{Ex}'$, $\mathrm{Ex}_2=\mathrm{Ex}_1'$, $\mathsf N_1(u)=(1-\mathsf L(u))/u$ are evaluated by
their power series with explicit tail bounds (`ad.py`), so they are valid on balls containing $0$.

**Lemma 1.2 (eigenvalues ↔ zeros).** (PV Prop. 5.3) If $\lambda\ne-1-A\pm2\sqrt A$ is an eigenvalue of $B$, then
$F(\varphi)=0$ for both roots $\varphi$, $1/(A\varphi)$ of $A\varphi^2-(\lambda+A+1)\varphi+1=0$. (PV Lemma 5.4,
with the correction of result 8's verifier) if $S(\varphi)\ne0$, $\varphi$ is a simple zero of $F$ and
$f(\varphi)\ne0$ or $f(1/(A\varphi))\ne0$, where $\varphi f(\varphi)=A\varphi^M+\varphi-1=AE+\varepsilon^2z$, then
$\lambda$ is a simple eigenvalue. (PV Lemma 3.1) $0$ is a simple eigenvalue of $B$, removed by restriction to $H$.
Since $\lambda=-1-A\pm2\sqrt A<0$, Lemma 1.2 applies to every eigenvalue with $\operatorname{Re}\lambda\ge0$.

---------------------------------------------------------------------------------------------------------------

## 2. Hopf points, uniformly over $\varepsilon$-boxes

Unknowns $X=(\kappa,x,y)$, $z=x+iy$; equations $G(X;\varepsilon)=(\operatorname{Re}Q^\varepsilon,
\operatorname{Im}Q^\varepsilon,\operatorname{Re}\Lambda)=0$ (the same system as result 8, written in regular form).

**Lemma 2.1 (parametric Krawczyk).** Let $E=[\varepsilon_{lo},\varepsilon_{hi}]$ with midpoint $\varepsilon_m$, $X$ a
box with midpoint $m$, $Y$ a real matrix, and
$$\mathcal K=m-Y\,G(m;\varepsilon_m)-\big(Y\,\partial_\varepsilon G(m;E)\big)(E-\varepsilon_m)+\big(I-Y\,DG(X;E)\big)(X-m),$$
where $\partial_\varepsilon G(m;E)$ and $DG(X;E)$ are interval enclosures of $\{\partial_\varepsilon G(m;\varepsilon)\}$
and $\{D_XG(\xi;\varepsilon):\xi\in X\}$ over $\varepsilon\in E$ (the product $Y\,\partial_\varepsilon G$ is formed
before multiplying by the scalar interval $E-\varepsilon_m$). If $\mathcal K\subset\operatorname{int}X$, then for
every $\varepsilon\in E$, $G(\cdot;\varepsilon)$ has exactly one zero in $X$.

*Proof.* Fix $\varepsilon\in E$. By the mean value theorem applied to each component,
$G(m;\varepsilon)=G(m;\varepsilon_m)+g'\,(\varepsilon-\varepsilon_m)$ with a vector $g'\in\partial_\varepsilon G(m;E)$ (the
components may use different intermediate points), so $Y\,G(m;\varepsilon)\in YG(m;\varepsilon_m)+(Y\partial_\varepsilon
G(m;E))(E-\varepsilon_m)$. Hence the Krawczyk operator of $G(\cdot;\varepsilon)$ (R. Krawczyk, Computing 4 (1969);
R. E. Moore, SIAM J. Numer. Anal. 14 (1977)) satisfies $K_\varepsilon(X)\subset\mathcal K\subset\operatorname{int}X$,
which gives existence and uniqueness. $G$ is $C^1$ in $X$, and $C^1$ in $\varepsilon\in[0,\varepsilon_{hi}]$ up to
$\varepsilon=0$ (the $T$-part and its derivatives tend to $0$). $\square$
(`allm_core.krawczyk`; forming $Y\partial_\varepsilon G$ first avoids the wrapping effect, which on the second branch
near $M\approx230$ made the boxes $\approx40$ times too large. The box is found by $\epsilon$-inflation from the
variation of Newton solutions at the ends of $E$, `h2_lib.hopf_box`; only the final inclusion test matters for the
proof, and it is re-done from the stored box by `verify_sweep.py`.)

Derivatives are computed by forward dual numbers in $(z,\kappa,\varepsilon)$ (`allm_core.DN`), with
$dE/dz=E/\varphi$ and $dE/d\varepsilon=2\varepsilon z^2E\,\mathsf M_1(\varepsilon^2z)$ (exact).

**Checks on the whole box** (`allm_core.hopf_checks`), all certified:
* *simplicity*: $0\notin\partial_zQ^\varepsilon(X;E)$, so $\varphi$ is a simple zero of $F$;
* $|AE+\varepsilon^2z|>0$ ($f(\varphi)\ne0$) and $\operatorname{Im}z>0$; also $S(\varphi)\ne0$ (for $\varepsilon>0$,
  $|\varphi-1|=\varepsilon^2|z|>0$ and $|\varphi-1/A|,|\varphi\mp A^{-1/2}|\ge K/A-\varepsilon^2|z|>0$ on the boxes);
  by Lemma 1.2, $\pm i\omega$ ($\omega=\kappa\varepsilon^3\operatorname{Im}\Lambda$) are simple eigenvalues of $B$;
* *transversality*: $d\operatorname{Re}\Lambda/d\kappa>0$ along the root curve $z(\kappa)$
  ($dz/d\kappa=-\partial_\kappa Q/\partial_zQ$); since $\lambda=\kappa\varepsilon^3\Lambda$,
  $\operatorname{Re}d\lambda/d\kappa=\kappa\varepsilon^3\,d\operatorname{Re}\Lambda/d\kappa>0$ at
  $\operatorname{Re}\Lambda=0$.

---------------------------------------------------------------------------------------------------------------

## 3. Localization of the unstable spectrum (explicit form of PV Props. 8.1, 8.4, 9.1, 9.2)

Write $w=\varphi-1$, $a=|w|$, $u=\operatorname{Re}w$.

**Lemma 3.1 (region).** (*check F3*) $\lambda=w(K+Aw)/(1+w)$ and
$\operatorname{Re}\lambda=Au-(u+a^2)/|1+w|^2$, $|1+w|^2=1+2u+a^2$. Hence, for $|\varphi|\ge1$:
$\operatorname{Re}\lambda\ge0\iff a^2(1-Au)\le u(K+2Au)$; and then $u\ge0$ (so $\operatorname{Re}\varphi\ge1$,
$|K+Aw|^2=K^2+2KAu+A^2a^2\ge\max(K,Aa)^2$) and $u\ge u_{\min}(a)$, the positive root of
$2Au^2+(K+Aa^2)u-a^2=0$; $u_{\min}$ is increasing in $a$ and decreasing in $K$.

*Proof of $u\ge0$:* with $\varphi=re^{i\theta}$, $r\ge1$: $\operatorname{Re}\lambda=\cos\theta(Ar+1/r)-(A+1)\ge0$ gives
$\operatorname{Re}\varphi=r\cos\theta\ge(A+1)r^2/(Ar^2+1)\ge1$ (PV Lemma 8.3). The rest is algebra. $\square$

**Lemma 3.2 (Gershgorin).** If $K<1$ and $\operatorname{Re}\lambda\ge0$ then $|\lambda|\le2A+2$, hence (with
$|\lambda|\ge a\,Aa/(1+a)$) $a\le a_{\max}:=2+2\sqrt2$.

*Proof.* Row 1 of $B$: centre $-A(M+2)$, radius $2+2MK-4K$; $A(M+2)-(2+2MK-4K)=M(1-K)+6K>0$, so that disc lies in
$\operatorname{Re}<0$. Rows $2..N$: centre $-(A+1)$, radius $A+1$; row $M$: centre $-A$, radius $2A$. $\square$

**Lemma 3.3 (localization, explicit).** Let $R_0>0$, $0<a_1<a_{\max}$, and a parameter box
$\varepsilon\in E$, $\kappa\in\mathcal K$ with $K<1$. Suppose (I) and (II) hold:

(I) $\rho_a+\rho_b+\rho_c+\rho_d+\rho_e+\rho_f<1$ for all $(\varepsilon,\kappa)$ in the box, where
$$\rho_a=(1+a_1)\Big(\frac\kappa{AR_0}+\varepsilon\Big)^2,\quad\rho_b=\frac\varepsilon\kappa,\quad
\rho_c=(1+a_1)^2\Big(\frac{\kappa^2}{A^2R_0^2}+\frac{\kappa\varepsilon}{AR_0}\Big),\quad\rho_d=\frac{\varepsilon^2}{A^2},$$
$$\rho_e=\frac{\kappa^2+\kappa\varepsilon}{A^2R_0}+\frac{(\kappa\varepsilon+\varepsilon^2)(2+a_1)}A,\quad
\rho_f=\frac{N_{f0}\,T\varepsilon^{-3}+N_{f1}\,T\varepsilon^{-5}}{A^2R_0^2\kappa},$$
$N_{f0}=A^2((K+Aa_1)^2(1+a_1)^2+a_1^2)+Aa_1^3+A^2(1+a_1)^2(K+Aa_1)^2$, $N_{f1}=A^2(1+a_1)(K+Aa_1)^2a_1$;

(II) on every subinterval $[a_-,a_+]$ of a cover of $[a_1,a_{\max}]$:
$$\Big(1+2u_{\min}(a_-)+a_-^2\Big)^{M_{lo}/2}\,\underline{\tfrac{|P_2|}M}-\overline{\tfrac{|P_1|}M}-\overline{\tfrac{A^{-M}(|R_1|+|R_2|)}M}>0,$$
with $M_{lo}=\varepsilon_{hi}^{-2}$ and the monotone bounds (all in terms of $u\in[u_{\min}(a_-),a_+]$,
$a\in[a_-,a_+]$, $1/M\in(0,\varepsilon_{hi}^2]$):
$|P_2|\ge MA^2a^2|K+Aw|-|\varphi||K+Aw|^3-A^2a^2$,
$|P_1|\le K(|K+Aw|^2|\varphi|^2+a^2)+(KM+1)a|K+Aw|(K+2Aa+Aa^2)$,
$A^{-M}(|R_1|+|R_2|)\le T[A^2(|K+Aw|^2|\varphi|^2+a^2)+Aa^3+A^2|\varphi||K+Aw|^2(Ma+|\varphi|)]$.

Then for every integer $M$ with $M^{-1/2}\in E$ and $\kappa\in\mathcal K$: every zero $\varphi$ of $F$ with $|\varphi|\ge1$,
$\operatorname{Re}\lambda(\varphi)\ge0$ and $|\varphi-1|\le a_{\max}$ satisfies $|z|=M|\varphi-1|<R_0$; by Lemmas 1.2
and 3.2 this applies to every eigenvalue $\lambda\ne0$ with $\operatorname{Re}\lambda\ge0$ (taking the root
$\varphi$ with $|\varphi|\ge A^{-1/2}$, which then has $|\varphi|\ge1$ by Lemma 3.1's proof). If $\varepsilon_{lo}>0$, (I) may be
replaced by (II) on $[R_0\varepsilon_{lo}^2,a_1]$.

*Proof.* Let $\varphi$ be such a zero with $a\ge R_0/M$. By Lemma 3.2, $a\le a_{\max}$. Put
$\mathcal T=MA^2a^2|K+Aw|$. Using $|\varphi|\le1+a$, $|K+Aw|\in[\max(K,Aa),K+Aa]$, $1/a\le M/R_0$, $K\sqrt M=\kappa$:
the six terms of $|P_2|-\mathcal T$ and of $|P_1|$ and of the exponentially small part (with
$|(A\varphi)^{-M}|\le A^{-M}$ since $|\varphi|\ge1$) are bounded by $\rho_a\mathcal T,\dots,\rho_f\mathcal T$ for
$a\le a_1$ (each bound is the elementary computation recorded in `allm_core.loc_check_I`; e.g.
$|\varphi||K+Aw|^3/\mathcal T\le(1+a)(K/(Aa)+1)^2/M\le(1+a_1)(\kappa/(AR_0)+\varepsilon)^2$ and
$(KM+1)a|K+Aw||A\varphi^2-1|/\mathcal T\le(K+M^{-1})(K+2Aa+Aa^2)/(A^2a)\le\rho_e$). Then
$|\varphi^MP_2|\ge|P_2|>|P_1|+A^{-M}|R_1|+|(A\varphi)^{-M}R_2|$, so $F(\varphi)\ne0$. For $a\in[a_1,a_{\max}]$ use
$|\varphi|^M=(1+2u+a^2)^{M/2}\ge(1+2u_{\min}(a_-)+a_-^2)^{M_{lo}/2}$ and (II). $\square$

**Lemma 3.3′ (localization without Gershgorin, any $K$).** Replace (II) by (II′): on every subinterval of a cover of
$[a_1,a_{big}]$ ($a_{big}=12$) the inequality of (II) holds with $\underline{|P_2|/M}$ replaced by
$\max(\mathrm{LB}_a,\mathrm{LB}_b)$, where $\mathrm{LB}_a$ is the bound of (II) and, with $r=|\varphi|\in[(1+2u_{\min}(a_-)+a_-^2)^{1/2},1+a_+]$,
$$\mathrm{LB}_b=A^3r^3-\frac{3A^2r^3+(3A+A^2)r^2}M-A^2(1+2A)r^2-A^2(2+A)r-\frac{(1+2A^2)r}M-A^2\Big(1+\frac1M\Big),$$
which follows from $P_2=A^3\varphi^3(\varphi+M)-3A^2\varphi^3+c_2\varphi^2+c_1\varphi+c_0$ (*check F7*),
$c_2=3A-MA^2(1+2A)-A^2$, $c_1=-1+MA^2(2+A)+2A^2$, $c_0=-(MA^2+A^2)$, and $|\varphi+M|\ge M+1$ for
$\operatorname{Re}\varphi\ge1$; and (III) on $a\ge a_{big}$: with $r\ge r_b=a_{big}-1\ge2$, $a\le1.5r$,
$|K+Aw|\le\alpha r$, $K+2Aa+Aa^2\le\beta r^2$ ($\alpha=K/2+1.5A$, $\beta=K/4+3.75A$, *check F7*):
$|P_1|/M\le C_1r^4$, $A^{-M}(|R_1|+|R_2|)/M\le C_2r^4$ and $|P_2|/M\ge c\,r^3$ with explicit $C_1,C_2,c$
(`allm_core.loc_tail`), and $r_b^{M_{lo}-1}c>C_1+C_2$. Then the conclusion of Lemma 3.3 holds for **every** zero
$\varphi$ of $F$ with $|\varphi|\ge1$ and $\operatorname{Re}\lambda(\varphi)\ge0$, without the condition $K<1$.
*Proof.* As for Lemma 3.3 on $[R_0/M,a_1]$ and $[a_1,a_{big}]$; on $a\ge a_{big}$:
$|\varphi^MP_2|\ge r^M c\,Mr^3>(C_1+C_2)Mr^4\ge|P_1|+A^{-M}(|R_1|+|R_2|)$ because $r^{M-1}\ge r_b^{M_{lo}-1}$. $\square$
(`allm_core.loc_check_II2`, `loc_check_II2_range`, `loc_tail`; used for $j\ge2$, `h2j_lib`.)

Implemented: `allm_core.loc_check_I`, `loc_check_II`, `loc_check_II_range` (adaptive subdivision); the
parameters used: $R_0=8$ ($j=1$), $20$ ($j=2$), $26$ ($j=3$), $a_1=0.1$. For boxes near $\varepsilon=0$ the
sum in (I) is $0.58$ ($j=1$), $0.56$ ($j=2$), $0.68$ ($j=3$); for moderate $M$ (I) fails and (II) on
$[R_0\varepsilon_{lo}^2,a_{\max}]$ is used (the growth factor $|\varphi|^M$ is then large).

---------------------------------------------------------------------------------------------------------------

## 4. Counting zeros (argument principle in ball arithmetic)

**Lemma 4.1 (winding number).** Let $\Gamma$ be a closed counterclockwise polygon and $(\varepsilon,\kappa)$ range over
a box. Subdivide every edge until the mean-value enclosure of $Q^\varepsilon$ over the edge and the parameter box
excludes $0$ (`allm_core.Q_centered`; the $T$-part as a disk). On each such sub-edge the image lies in a convex set
not containing $0$, so the argument increment is in $(-\pi,\pi)$ and equals the principal difference of the arguments
at the endpoints; summing these (computed at the midpoint parameters from balls of relative radius $<10^{-6}$) gives
$2\pi$ times an integer, which by homotopy invariance is the number of zeros of $Q^\varepsilon(\cdot;\kappa)$ inside
$\Gamma$ for **every** parameter in the box (`allm_core.winding_count`).

**First branch** (`h2_lib`). $D=\{|z|\le R_0,\ \operatorname{Re}z\ge-\sigma\}$, $R_0=8$, $\sigma=0.3$; $\Gamma$ is the
convex polygon with vertices on the circle of radius $R_1=R_0/\cos(\Delta/2)$ (96 chords, each at distance $\ge R_0$
from $0$) closed by the chord $\operatorname{Re}z=-\sigma$, so that $D\subset\operatorname{int}\Gamma$
(`allm_core.contour_D`).

**Lemma 4.2 (region).** If $|z|\le R_0$ and $\operatorname{Re}z\le-\sigma$, then
$\operatorname{Re}\Lambda\le\max(-\sigma+\varepsilon\sigma^2/\kappa,\,-R_0+\varepsilon R_0^2/\kappa)+
\varepsilon^3R_0^3/(\kappa(1-\varepsilon^2R_0))<0$ (certified on each box; `allm_core.region_check`).

*Proof.* $\operatorname{Re}\Lambda=x+\frac\varepsilon\kappa\operatorname{Re}z^2+\frac\varepsilon\kappa
\operatorname{Re}(z^2(1/\varphi-1))$, $|1/\varphi-1|\le\varepsilon^2R_0/(1-\varepsilon^2R_0)$; if $|y|\ge|x|$,
$\operatorname{Re}z^2\le0$; otherwise $\operatorname{Re}z^2\le x^2$ and $x+\varepsilon x^2/\kappa$ is convex on
$[-R_0,-\sigma]$. $\square$

**Theorem 4.3 (H2, first branch).** Suppose that on a box $E\times X$ the checks of Section 2, Lemma 3.3
($R_0=8$), Lemma 4.2 ($\sigma=0.3$) hold, that $|z_H|<R_0$, $\operatorname{Re}z_H>-\sigma$ on $X$, and that the winding
count of Lemma 4.1 on $\Gamma$ equals $4$. Then for every integer $M$ with $M^{-1/2}\in E$, at $\kappa=\kappa_1(M)$
(the zero in $X$): the eigenvalues of $B$ with $\operatorname{Re}\lambda\ge0$ are exactly $0,\pm i\omega$, all
simple; all other eigenvalues of $B|_H$ have negative real part.

*Proof.* $z=0$ is a zero of multiplicity $\ge2$ ($F(1)=F'(1)=0$, *check F3*, PV Prop. 9.1(i)), and $z_H\ne\bar z_H$
are zeros in $D$: these are 4 zeros, so there are no others in $\Gamma$. Let $\lambda\ne0$ be an eigenvalue with
$\operatorname{Re}\lambda\ge0$; by Lemma 1.2 and Lemma 3.1 there is a zero $\varphi$ of $F$ with $|\varphi|\ge1$ and
$\lambda=\lambda(\varphi)$; by Lemma 3.3 $|z|<R_0$; since $\operatorname{Re}\Lambda=M\operatorname{Re}\lambda/K\ge0$,
Lemma 4.2 gives $\operatorname{Re}z>-\sigma$; so $z\in D$, $z\ne0$, hence $z\in\{z_H,\bar z_H\}$ and
$\lambda=\pm i\omega$. Simplicity: Section 2 and PV Lemma 3.1. $\square$

**Branches $j\ge2$** (`h2j_lib`). The curve $\operatorname{Re}\Lambda=0$ bends to the right
($x\approx\varepsilon y^2/\kappa$), so contours follow level curves of $\operatorname{Re}\Lambda$.

**Lemma 4.4 (monotonicity).** $\Lambda'(z)=1+\frac\varepsilon\kappa\frac{z(2+\varepsilon^2z)}{\varphi^2}$; if
$\frac\varepsilon\kappa\rho(2+\varepsilon^2\rho)/(1-\varepsilon^2\rho)^2<1$ with $\rho$ the largest $|z|$ on the
rectangle $\mathcal R=[x_{\min},R_0]\times[-Y,Y]$, then $\partial_x\operatorname{Re}\Lambda=\operatorname{Re}\Lambda'>0$
on $\mathcal R$ (certified: `h2j_lib.mono_check`).

**Theorem 4.5 (spectrum at the $j$-th Hopf point).** Let $\ell(y)$ and $r(y)$, $|y|\le Y$, be piecewise linear graphs
with $\operatorname{Re}\Lambda<0$ on $x=\ell(y)$ and $\operatorname{Re}\Lambda>0$ on $x=r(y)$ (interval checks on
every edge, for all parameters of the box), $R_0\le Y$, $x_{\min}=\min\ell<0$, and $P=\{\ell(y)<x<R_0,|y|<Y\}$,
$P_-=\{\ell(y)<x<r(y),|y|<Y\}$. If Lemma 3.3′ (with $R_0$), Lemma 4.4 on $\mathcal R=[x_{\min},R_0]\times[-Y,Y]$,
Lemma 4.2 with $\sigma=-x_{\min}$, $|y_H|<Y$, and the counts $N(P)=2j+2$, $N(P_-)=4$ hold, then at $\kappa_j(M)$: the only eigenvalues of $B$ on $i\mathbb R$ are $0,\pm i\omega$,
and exactly $2j-2$ zeros of $F$ with $|\varphi|\ge1$ (with multiplicity) correspond to eigenvalues with
$\operatorname{Re}\lambda>0$.

*Proof.* Let $\lambda\ne0$ be an eigenvalue with $\operatorname{Re}\lambda\ge0$ and $\varphi$ its root with
$|\varphi|\ge A^{-1/2}$ (so $|\varphi|\ge1$, Lemma 3.1). By Lemma 3.3′, $|z|<R_0$; hence $|y|<Y$, $x<R_0$, and
$x>x_{\min}$ (Lemma 4.2), so $z\in\mathcal R$; by monotonicity along the horizontal line, $x>\ell(y)$: $z\in P$.
If moreover $\operatorname{Re}\lambda=0$, then $\operatorname{Re}\Lambda(z)=0$ and monotonicity gives
$\ell(y)<x<r(y)$: $z\in P_-$. Now $P_-$ contains the four zeros $0,0,z_H,\bar z_H$, hence no other; so
$z\in\{z_H,\bar z_H\}$ and $\lambda=\pm i\omega$. Conversely, a zero of $Q^\varepsilon$ in $P\setminus P_-$ lies in
$\mathcal R$ with $x\ge r(y)$, so $\operatorname{Re}\Lambda>0$; it is not one of the spurious zeros of $F$
($\varphi=1/A,\pm A^{-1/2}$ have $\operatorname{Re}z\le M(A^{-1/2}-1)<x_{\min}$, and $\varphi=1$ is $z=0\in P_-$), so
by PV Prop. 5.3 it gives an eigenvalue $\lambda(\varphi)$ with $\operatorname{Re}\lambda>0$; and every eigenvalue
with $\operatorname{Re}\lambda>0$ arises in this way from a zero in $P\setminus P_-$ (first part of the proof). These
zeros have total multiplicity $N(P)-N(P_-)=2j-2$. $\square$

Parameters: $j=2$: $R_0=20$, $Y=20.4$, $\ell$ on $\operatorname{Re}\Lambda\approx-0.3$, $r$ on
$\operatorname{Re}\Lambda\approx0.25$; $j=3$: $R_0=26$, $Y=26.6$, levels $-0.15$ and $0.25$ (48 chain vertices,
shifted outwards by $0.02$, refined until the interval check holds).

---------------------------------------------------------------------------------------------------------------

## 5. The first Lyapunov coefficient in transfer-function form

Kuznetsov's formula with $\mathcal C=0$, $\mathcal B(x,y)=x_1Sy+y_1Sx$ (result 8, §3), $q_1=1$, $p$ with
$pB=i\omega p$, $pq=1$, gives, with $r(\mu):=(\mu-B)^{-1}Sq$ ($r(0)$ the solution in $H$) and the linear functional
$\sigma(y)=pSy$:
$$\ell_1=\frac1\omega\operatorname{Re}\Big[\sigma(r(0))+\sigma(\overline{r(0)})+2\operatorname{Re}(r(0)_1)\,\sigma(q)
+\sigma(r(2i\omega))+r(2i\omega)_1\,\sigma(\bar q)\Big].\tag{5.1}$$
*Derivation.* $\mathcal B(q,\bar q)=Sq+S\bar q$, so $w_{11}=(B|_H)^{-1}\mathcal B(q,\bar q)=-(r(0)+\overline{r(0)})$;
$w_{20}=(2i\omega-B)^{-1}2Sq=2r(2i\omega)$; $\langle p,\mathcal B(q,w)\rangle=\sigma(w)+w_1\sigma(q)$; insert into
$\ell_1=\frac1{2\omega}\operatorname{Re}[-2\langle p,\mathcal B(q,w_{11})\rangle+\langle p,\mathcal B(\bar q,w_{20})\rangle]$.
(Checked: dense Kuznetsov formula at $M=40,100$, *check F5*; result 8's closed form at $M=10^3..10^6$, 20 digits.)

**Lemma 5.1 (row 1 = mass).** If $r$ satisfies rows $2,\dots,M$ of $(\mu-B)r=Sq$, then row 1 holds iff
$\mu\sum_\ell\ell r_\ell=0$. *Proof.* $\boldsymbol\ell B=\boldsymbol\ell S=0$ give
$(Br)_1=-\sum_{\ell\ge2}\ell(Br)_\ell=-\sum_{\ell\ge2}\ell(\mu r_\ell-(Sq)_\ell)$ and
$\sum_{\ell\ge2}\ell(Sq)_\ell=-(Sq)_1$. $\square$
So $r(\mu)$ is determined by the bulk rows, row $M$, and $\sum\ell r_\ell=0$ (for $\mu=0$ the last condition selects
$r\in H$). This replaces the monomer row, whose terms are $O(\varepsilon^{-3})$ and cancel.

### 5.1 Atoms and closed forms

All vectors are finite sums of *atoms* in $\ell=1..M$ (`l1reg`):
* $\mathrm{Slow}(c,\tilde b)$: $c\,\beta^\ell$, $\beta=1+\varepsilon^2\tilde b$, $b:=M\log\beta=\tilde b\,\mathsf L(\varepsilon^2\tilde b)$;
* $\mathrm{BM}(c,\tilde g)$: $c\,\gamma^{\ell-M}$, $\gamma=A(1+\varepsilon^2\tilde g)$ (boundary layer at $\ell=M$; $|\gamma^{-M}|=T e^{-\operatorname{Re}\tilde g\mathsf L(\varepsilon^2\tilde g)}$);
* $\mathrm{Dd}(c,\tilde g_1,\tilde g_2)$: $c(\gamma_1^{\ell-M}-\gamma_2^{\ell-M})/(\gamma_1-\gamma_2)$ (two nearly equal layers);
* $\mathrm B1(c,\tilde p)$: $c\,\psi'^\ell$, $\psi'=(1+\varepsilon^2\tilde p)/A$ (layer at $\ell=1$).

Coefficients are *scaled numbers* $c=v\varepsilon^p$ ($v$ regular, $p\in\mathbb Z$ tracked exactly). Closed forms
(all exact; regular part + exponentially small part, the latter enclosed in a disk by Lemma 1.1), with
$\sum'=\sum_{k=1}^{N}$:
* $\sum'\beta_1^k\beta_2^k=\varepsilon^{-2}[\mathrm{Ex}(b_1+b_2)\mathsf L(u)-\varepsilon^2]$, $u=\beta_1\beta_2-1$ (no
  exponentially small part);
* $\sum'\beta^k\gamma^{k-M}=\dfrac{\beta^M-\beta\gamma^{1-M}}{\beta\gamma-1}$, $\beta\gamma-1=\varepsilon d$,
  $d=\kappa+A\varepsilon(\tilde b+\tilde g+\varepsilon^2\tilde b\tilde g)$: regular part $e^{b}/(\varepsilon d)$;
* $\sum'\beta^k\,\mathrm{Dd}_k$: regular part $-\beta e^b/(\varepsilon^2d_1d_2)$ (the divided difference done
  analytically); the exponentially small part by $\sup|f'|$ on the segment, $f(\gamma)=\gamma^{1-M}/(\beta\gamma-1)$;
* $\sum'\psi'^k\beta^k$: regular part $\psi'\beta A/(\varepsilon d')$, $d'=\kappa-\varepsilon(\tilde p+\tilde b+\varepsilon^2\tilde p\tilde b)$;
  $\sum'\psi'^k\gamma^{k-M}$ and $\sum'\psi'^k\mathrm{Dd}_k$ are exponentially small;
* $\sum_{\ell=1}^M\ell\beta^\ell=\varepsilon^{-4}\beta[\mathsf L(u)^2\mathrm{Ex}_1(b)+\varepsilon^2e^b\mathsf N_1(u)]$, $u=\varepsilon^2\tilde b$;
* $\sum_{\ell=1}^M\ell\gamma^{\ell-M}=\varepsilon^{-3}[\gamma/d_\gamma-\varepsilon\gamma/d_\gamma^2]+\gamma^{1-M}/(\gamma-1)^2$, $\gamma-1=\varepsilon d_\gamma$;
* $\sum_\ell\ell\,\mathrm{Dd}_\ell=\varepsilon^{-4}[-(1-\varepsilon^2)/(d_1d_2)+\varepsilon(d_1+d_2)/(d_1^2d_2^2)]$ + exp. small;
* values at $\ell=1,N,M$.
*Check F6*: all closed forms agree with direct summation for random atoms at $M=37,200$ (enclosures overlap).

### 5.2 The blocks

With $g,h,\hat\Lambda,E$ as in Section 1, $T_q=T/\varepsilon$ (a disk):
* eigenvector ($q_1=1$): $q=C_1E\,(1/\varphi)^\ell+\varepsilon^{-1}\hat C_2(A\varphi)^{\ell-M}$,
  $C_1=\dfrac{g-A^2\varphi T_q/E}{gE/\varphi+\varepsilon^2zAT_q/E}$, $\hat C_2=\dfrac{AE+\varepsilon^2z}{gE+\varepsilon^2zA\varphi T_q/E}$
  (from PV's $(c_1,c_2)=(f(\psi),-f(\varphi))$, $\psi=1/(A\varphi)$);
* adjoint: $p_k=a_1\varphi^k+a_2\psi^k+p_c$, $\theta=\hat\Lambda/(\kappa+\varepsilon^2\hat\Lambda)$,
  $c_A=(\kappa+\varepsilon)\theta+\kappa$, $a_1=c_A\psi+\varepsilon^2T_q(Az-\varepsilon\hat\Lambda)/E$,
  $a_2=-(c_A\varphi-E(g/\varphi+\varepsilon^2\hat\Lambda))$, $p_c=-(1-\varepsilon^2\theta)(a_1\varphi+a_2\psi)$
  (columns $2..N$ of $pB=\lambda p$ give the form; column $M$, divided by $\varepsilon$, gives $(a_1,a_2)$; the
  normalization $\langle p,q\rangle=\sum'p_kq_k+p_Mq_M$, of order $\varepsilon^{-2}$);
* $(pS)_k=p_{k+1}-p_k-p_1=a_1\varepsilon^2z\varphi^k-a_2\frac{\varepsilon g}{A\varphi}\psi^k-\varepsilon^2\theta(a_1\varphi+a_2\psi)$
  for $1\le k\le N$, $(pS)_M=0$ (divided by $\langle p,q\rangle$): orders $\varepsilon^4,\varepsilon^3,\varepsilon^4$;
* resolvent, $\mu=\varepsilon^3m$ ($m=0$, or $m=2\hat\Lambda$ for $\mu=2\lambda=2i\omega$), slow root
  $s_1=1+\varepsilon^2\zeta$ of $s^2-(A+1+\mu)s+A=0$ ($\zeta=-2m/(c_m+\sqrt{c_m^2+4\varepsilon m})$,
  $c_m=\kappa+\varepsilon^2m$; $\zeta=0$ for $\mu=0$), $s_2=A/s_1$:
  $$r=\frac{zC_1E}{\varepsilon(m-\hat\Lambda)}(1/\varphi)^\ell+\frac{\hat C_2g}{\varepsilon\Xi}\,\mathrm{Dd}_\ell(A\varphi,s_2)
  +\frac{\hat h_1}\varepsilon s_1^\ell+\frac{\hat H_2}{\varepsilon^2}s_2^{\ell-M},\qquad
  \Xi=\frac{\kappa-\varepsilon\zeta}{1+\varepsilon^2\zeta}+A\varepsilon z-\varepsilon^2m;$$
  the particular part solves the bulk rows termwise (both bases of $q$ have $D_\mu(\beta)=\lambda-\mu$); the two
  boundary layers $c_p(A\varphi)^{\ell-M}$ and $h_2s_2^{\ell-M}$, which are each $O(\varepsilon^{-3})$ and cancel to
  $O(\varepsilon^{-2})$, are rewritten exactly as $c_p(w-s_2)\mathrm{Dd}_\ell+(h_2+c_p)s_2^{\ell-M}$ using
  $(w-s_2)(w+s_2-A-1-\mu)=-(\mu-\lambda)w$, $w=A\varphi$ (*check F4*); $(\hat h_1,\hat H_2)$ solve the $2\times2$
  system (row $M$, mass), whose rows are scaled by $\varepsilon$ and $\varepsilon^5$; the row-$M$ coefficient of
  $s_2^{\ell-M}$ is written exactly as $-\varepsilon(\kappa-\varepsilon\zeta+\varepsilon^2m)$ + exp. small.

Every block then has its natural order (verified: $\sigma(q),\sigma(\bar q)=O(\varepsilon^2)$,
$r(\cdot)_1=O(\varepsilon^{-1})$, the five terms of (5.1) $=O(\varepsilon)$), and with $\omega=\varepsilon^3\operatorname{Im}\hat\Lambda$:
$$\widehat L=\varepsilon^2\ell_1(q_1=1)=\frac{\operatorname{Re}\,b}{\operatorname{Im}\hat\Lambda},\qquad
b:=\varepsilon^{-1}\big[\text{bracket of }(5.1)\big],\tag{5.2}$$
where $b$ is computed as a *regular* expression: no power $\varepsilon^{p}$ with $p<0$ remains, the
$T$-terms are disks. So (5.2) can be evaluated on $\varepsilon$-boxes containing $0$.

**Lemma 5.2 (mean value form).** Let $\Phi(\kappa,x,y,\varepsilon;\nu)$ be the function computed by `l1reg.blocks`
when each exponentially small expression is replaced by an independent variable $\nu_i$ ($\Phi$ is smooth in
$(\kappa,x,y,\varepsilon)$ for fixed $\nu$, and $\widehat L=\Phi(\cdot;\nu(\cdot))$). For a box $\mathcal B$ with
midpoint $c$ and the disks $\nu_i\in\mathcal D_i$ valid on $\mathcal B$:
$\widehat L(\xi)\in\Phi(c;\mathcal D)+\nabla\Phi(\mathcal B;\mathcal D)\cdot(\mathcal B-c)$ for all $\xi\in\mathcal B$.
*Proof.* $\widehat L(\xi)=\Phi(c;\nu(\xi))+[\Phi(\xi;\nu(\xi))-\Phi(c;\nu(\xi))]$ and the mean value theorem in the
first four variables at fixed $\nu$. $\square$ The gradient is computed by forward AD in 4 real directions (`ad.py`,
`l1reg.Lhat_box`). The disks $\mathcal D_i$ are computed once, from the box evaluation (all inputs are balls over
$\mathcal B$), and **replayed** in the midpoint evaluation, so that $\nu(\xi)\in\mathcal D$ for every
$\xi\in\mathcal B$ (`l1reg.expsmall`, record/replay). Example: on $\varepsilon\in[0,10^{-8}]$ this gives $\widehat L=-0.0093114\pm3.6\cdot10^{-8}$
(naive interval evaluation: $\pm2.9\cdot10^{-5}$).

---------------------------------------------------------------------------------------------------------------

## 6. Regularity at $\varepsilon=0$; the limit

Because (5.2) is evaluated with $\varepsilon$-balls containing $0$, the certificates prove at once that
$\widehat L$ is bounded on every box $[0,\delta]$ and give its range. Continuity at $\varepsilon=0$ follows from the
construction (every operation is continuous on the closed box, the $T$-terms and their derivatives tend to $0$), so
$\lim_{M\to\infty}\ell_1(q_1=1)/M=\widehat L(0)$, the value of (5.2) at $\varepsilon=0$, $(\kappa,z)=(\kappa_1^0,it_1)$.
Certified: $\widehat L(0)\in[-0.00931144,-0.00931137]$ (box $[0,10^{-8}]$); result 8's formal value $L_1=-0.0093114$
(Proposition B of result 8) lies inside. The identity between (5.2) at $\varepsilon=0$ and result 8's closed form
$L_j$ is **not** proved symbolically here; it is only confirmed to the accuracy of the enclosures.

---------------------------------------------------------------------------------------------------------------

## 7. Certificates (sweeps over $\varepsilon$)

Each certificate file has one line per box $[\varepsilon_{lo},\varepsilon_{hi}]$; consecutive boxes share endpoints
(identical decimal strings; `allm_core.eps_interval` rounds outwards), the first starts at $\varepsilon=0$.

**Theorem 7.1 (first branch, $M\ge161$).** `certificates/l1h2_j1.csv` (`l1_sweep.py 1 0 0.0788111 --h2`): 654
boxes covering $[0,0.0788111]\supset[0,161^{-1/2}]$ ($161^{-1/2}=0.0788110406\ldots$), all PASS: on each box the
Krawczyk test of Lemma 2.1, the checks of Section 2, $\widehat L<0$ on the whole box (Lemma 5.2), and the hypotheses of
Theorem 4.3 ($R_0=8$, $\sigma=0.3$, $a_1=0.1$; winding count 4). Runtime 528 s. Ranges: $\widehat L\in
[-0.0337,-0.0081]$; $\kappa_1\in[1.8979,3.3285]$.
An earlier H2-only sweep (`h2_sweep.py`, `certificates/h2_j1.csv`, 201 boxes) is redundant with it.

**Theorem 7.2 (branches 2 and 3, $\ell_1$, $M\ge161$).** `certificates/l1_j2.csv` (589 boxes) and
`certificates/l1_j3.csv` (834 boxes), `l1_sweep.py J 0 0.0788111`, cover $[0,0.0788111]$; all PASS: Krawczyk box of
the $j$-th Hopf point (continued from $(\kappa_j^0,it_j)$ at $\varepsilon=0$), simplicity, transversality,
$\operatorname{Im}z>0$, $f(\varphi)\ne0$, and $\widehat L<0$ on the whole box. Ranges: $j=2$:
$\kappa_2\in[3.1547,10.538]$, $\widehat L\in[-0.01066,-0.00169]$; $j=3$: $\kappa_3\in[4.0319,30.040]$,
$\widehat L\in[-0.00471,-0.000109]$. Runtimes 169 s and 262 s.

**Theorem 7.3 (branches 2 and 3, spectrum).** `certificates/h2j_j2.csv`, `h2j_j3.csv` (`h2j_sweep.py`, on the
boxes of Theorem 7.2, each Krawczyk box re-proved from the stored data): Theorem 4.5 with $R_0=20$, $Y=20.4$ and
levels $-0.3$, $0.25$ ($j=2$), $R_0=26$, $Y=26.6$, levels $-0.15$, $0.25$ ($j=3$). [RESULTS]

## 8. Consistency with result 8

`check_vs_result8.py 1 certificates/l1h2_j1.csv ../certificates/sweep_j1*.csv.gz`: for each of the 99 840 Hopf points
of result 8 with $161\le M\le10^5$, result 8's box (midpoint $\pm$ stored radius) lies inside our Krawczyk box for
the $\varepsilon$-box containing $M^{-1/2}$ (so the two certificates concern the same zero, by uniqueness), and result
8's value $\ell_1(q_1=1)/M$ lies inside our enclosure of $\widehat L$. PASS. (Same for $j=2,3$, Section 7.)

---------------------------------------------------------------------------------------------------------------

## 9. Scripts (PASS/FAIL)

| script | claims | runtime |
|---|---|---|
| `check_formulas.py` | F1–F6: identities (sympy), dense entrywise enclosure test of every block at $M=40,100,300$, sums vs direct summation | 2 min |
| `h2_sweep.py OUT [161]` | H2 for $j=1$, all $M\ge161$ (Theorem 4.3) | 1 min |
| `l1_sweep.py J EPS_LO EPS_HI OUT [--h2]` | Hopf boxes, nondegeneracy, $\widehat L<0$ (and H2 for $j=1$) | 9 min ($j=1$, `--h2`) |
| `h2j_sweep.py J IN OUT R0 Y cL cR` | Theorem 4.5 on the boxes of an `l1_sweep` file | ~2 s/box |
| `verify_sweep.py FILE EPS_LO EPS_HI [--h2] [STRIDE]` | re-proves every line of a certificate from the stored data, and coverage | as the sweep |
| `check_vs_result8.py J OURS R8FILES` | same Hopf points as result 8, $\ell_1$ values inside our enclosures | 1 min |

Libraries: `allm_core.py` (regular $Q^\varepsilon$, Krawczyk, localization, region, winding), `l1reg.py` (regular
$\widehat L$), `ad.py` (AD and series), `h2_lib.py`, `h2j_lib.py`.

---------------------------------------------------------------------------------------------------------------

## 10. What is proved, what is assumed, what is not proved

**Proved (computer-assisted, this work).** Theorem 7.1 (first branch, all $M\ge161$: Hopf point, simplicity,
transversality, H2, no unstable eigenvalue, $\ell_1<0$); the asymptotic enclosure of Section 6; the branch-2/3
statements of Section 7 on the stated ranges.

**Taken from result 8 and its verifier (independently verified there).** $3\le M\le14$ (no crossing), $M=15$
(subcritical), $16\le M\le160$ (supercritical, H2), branches 2, 3 for $M\le10^5$ ($\ell_1<0$) and $M\le160$ (H2).

**Standard results used, not re-proved.** The Hopf bifurcation theorem in $\mathbb R^n$ with the centre manifold
reduction and the sign rule of the first Lyapunov coefficient (Kuznetsov, *Elements of Applied Bifurcation
Theory*, 3rd ed., §3.5 and §5.4; the formula was validated in result 8 by a planar test, simulations and
BifurcationKit); the Krawczyk–Moore theorem; correctness of Arb (python-flint 0.9). PV Prop. 5.3, Lemma 3.1 and Lemma
5.4 (with the verifier's correction of its first step, re-checked here as $f(\varphi)\ne0$).

**Not proved.**
* For $j\ge4$ nothing beyond result 8's floating-point checks and the formal limit coefficient.
* That $\kappa_1(M)$ is the *first* loss of stability (for all $0<\kappa<\kappa_1(M)$) for $M>57$; the stability of
  the bifurcating orbits only needs the spectrum at $\kappa_1(M)$, which is proved.
* The identity of the limit $\widehat L(0)$ with result 8's closed-form $L_1$ (equal to $\sim10^{-7}$ relative, Section 6).
* Global statements (existence of periodic orbits far from the bifurcation, their stability for $\kappa$ not close
  to $\kappa_1$).
