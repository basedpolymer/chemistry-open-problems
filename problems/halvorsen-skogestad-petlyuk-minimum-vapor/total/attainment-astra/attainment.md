# Finite-stage attainment: a constructive prefractionator theorem and the remaining network problem

Attainment: GPT-6-Astra (effort xhigh).

Status: **partial result, not independently verified**. Full FTC attainment for
three or more components is **not proved here**. The new result is Theorem 1:
for every number of components and every real feed quality, an isolated
directly coupled prefractionator admits explicit positive finite-stage profiles
arbitrarily close to its preferred-split vapor flows. The unwanted extreme
components in its return streams decay exponentially. The return streams are
boundary data of this theorem; closing them through the other FTC columns remains
an additional problem. Section 6 distinguishes this result from full attainment.

The physical model and the meaning of the objective are those of
[`../foundations/class.md`](../foundations/class.md). The already established
stage-level lower bound is used only to state the outstanding target. Binary
attainment is Agent S's Proposition 10 in
[`../sections/sections.md`](../sections/sections.md), not a new result claimed here.

## 1. Target, provenance, and notation

For a normalized process feed \(z>0\), distinct positive volatilities
\(\alpha_1>\cdots>\alpha_N\), and real \(q\), let

\[
\sum_i\frac{\alpha_i z_i}{\alpha_i-\theta_j}=1-q,\qquad
\theta_j\in(\alpha_{j+1},\alpha_j),\qquad
H_j=\sum_{i>j}\frac{\alpha_i z_i}{\theta_j-\alpha_i}.
\]

The requested upper bound is: for every \(\delta>0\), a **closed** finite-stage
FTC steady state has every assigned-product purity at least \(1-\delta\) and
total positive external duty at most \(\max_j H_j+\delta\).
Theorem 1 below does not assert this closed-network statement.

The original open statement is: “We have not presented a complete proof, so the
above conclusion is a conjecture.” Halvorsen–Skogestad, Part 3, Section 6.3,
p. 626, [doi:10.1021/ie0108651](https://doi.org/10.1021/ie0108651).
The local text of this passage was checked. The algebraic preferred-split
construction used below is Theorem `achieve` of `../../paper.tex`.
The projective section identity is also recorded by Agents N and S; it is
derived again below because its signs matter.

In Sections 2–5 there are \(m\ge3\) local components. Write their volatilities as
\(a_1>\cdots>a_m>0\), their net feed as \(f_i>0\), and its vapor equivalent as
\(U\in\mathbb R\). The feed need not be normalized; its quality is
\(1-U/\sum_i f_i\). Its \(m-1\) intervening roots are denoted by
\(\theta_1,\ldots,\theta_{m-1}\).

A section has positive flows \(V,L\), normalized liquid composition \(x\), and
\[
S(x)=\sum_i a_i x_i,\quad y_i(x)=a_i x_i/S(x),\quad
G_x(t)=\sum_i\frac{a_i x_i}{a_i-t}.
\]
The constant net component flow is positive upwards:
\[
w=Vy(x_k)-Lx_{k+1},\qquad \sum_i w_i=V-L.
\]
Thus, with \(M=V\operatorname{diag}(a)-wa^{\mathsf T}\),
\[
x_{k+1}=\frac{Mx_k}{L S(x_k)}.                                      \tag{1}
\]
Indeed the numerator is \(Vy(x_k)S(x_k)-wS(x_k)\), and its component sum is
\((V-\sum_i w_i)S(x_k)=L S(x_k)\).
Conversely, a positive normalized sequence satisfying (1) satisfies every
cross-section component balance. Taking differences gives the interior tray
balances.

## 2. A slack deformation preserving the internal common roots

At the preferred split write the top and bottom net products as \(d^0,b^0\):
\[
d^0_1=f_1,\quad d^0_m=0,\quad b^0=f-d^0,\quad
d^0_i>0\ (i<m),\quad b^0_i>0\ (i>1).
\]
Write its positive vapor flows as \(V_T^0,V_B^0\), where
\(V_T^0-V_B^0=U\). The preferred-split rational identities are
\[
\begin{aligned}
\Psi_T^0(t)&=V_T^0-\sum_i\frac{a_i d_i^0}{a_i-t}
 =V_T^0\,\frac{\prod_{j=1}^{m-1}(t-\theta_j)}
                   {\prod_{i=1}^{m-1}(t-a_i)},\\
\Psi_B^0(t)&=V_B^0+\sum_i\frac{a_i b_i^0}{a_i-t}
 =V_B^0\,\frac{\prod_{j=1}^{m-1}(t-\theta_j)}
                   {\prod_{i=2}^{m}(t-a_i)}.                       \tag{2}
\end{aligned}
\]
These identities also give \(L_T^0=V_T^0-\sum d^0>0\) and
\(L_B^0=V_B^0+\sum b^0>0\). For completeness, their residues give the top formula
in Theorem `achieve`; its normalization is \(d^0_1=f_1\). Subtracting the feed
equation gives the second identity. Counting signs in these interlacing products
gives positivity of \(d^0,b^0,V_T^0,V_B^0\); evaluating at zero gives positive
liquid flows. Thus no finite-stage assertion is assumed in (2).

Empty products in what follows equal one. Set
\[
\beta=\frac{\prod_{i=2}^{m-1}a_i}{\prod_{j=2}^{m-2}\theta_j},\qquad
K(t)=
\frac{(t-\beta)\prod_{j=2}^{m-2}(t-\theta_j)}
     {\prod_{i=2}^{m-1}(t-a_i)}.                                   \tag{3}
\]
For \(m=3\), \(K=1\), after cancellation. For \(m>3\),
\(a_{m-1}<\beta<a_2\): divide the displayed product by \(a_2\) and compare
\(a_{j+1}<\theta_j\); divide by \(a_{m-1}\) and compare \(a_j>\theta_j\).

Define \(h_1=h_m=0\) and, for \(2\le i\le m-1\),
\[
h_i=
\frac{(a_i-\beta)\prod_{j=2}^{m-2}(a_i-\theta_j)}
 {a_i\prod_{\substack{2\le k\le m-1\\k\ne i}}(a_i-a_k)}.             \tag{4}
\]
Partial fractions, including removable poles, give
\[
K(t)=1-\sum_i\frac{a_i h_i}{a_i-t},\qquad \sum_i h_i=0.              \tag{5}
\]
For the last equality use \(K(0)=1\), which follows from (3).

For \(\eta>0\), prescribe
\[
\begin{gathered}
d=d^0+\eta h,\quad b=b^0-\eta h,\quad
V_T=V_T^0+\eta,\quad V_B=V_B^0+\eta,\\
L_T=L_T^0+\eta,\qquad L_B=L_B^0+\eta.                               \tag{6}
\end{gathered}
\]
All total net flows and the feed are unchanged. For sufficiently small \(\eta\)
the positive entries of \(d,b\) and all four phase flows remain positive. Moreover
\[
\Psi_T=\Psi_T^0+\eta K,\qquad \Psi_B=\Psi_B^0+\eta K.                \tag{7}
\]
Consequently \(\theta_2,\ldots,\theta_{m-2}\) remain common roots **exactly**.
This is more restrictive than adding the same vapor slack with every component
split fixed; that simpler prescription generally fails to preserve the internal
common roots when \(m\ge4\).

Denote the two remaining top roots by \(r,u\), and the two remaining bottom roots
by \(v,s\), the first being the larger in each case. For sufficiently small
\(\eta>0\),
\[
a_1>r>v>a_2,\qquad a_{m-1}>u>s>a_m.                               \tag{8}
\]
To prove this, the roots at \(\eta=0\) are simple. At either extreme feed root,
\(K(\theta_j)>0\), because the evaluation point lies outside all its poles and
zeros and numerator and denominator have the same sign. At a top root
\(\partial_t\Psi_T^0=-\sum_i a_i d_i^0/(a_i-t)^2<0\); at a bottom root
\(\partial_t\Psi_B^0=\sum_i a_i b_i^0/(a_i-t)^2>0\).
Implicit differentiation in (7) therefore moves both top roots upwards and
both bottom roots downwards. The root ordering and simplicity persist by
continuity. The middle roots cannot account for any additional eigenvalues:
each numerator has degree \(m-1\).

## 3. Positive feed matching

**Lemma 1 (interlacing interpolation).** If
\(a_i>\chi_i>a_{i+1}\), there is exactly one normalized vector \(X\) with
\(G_X(\chi_i)=0\), and it is strictly positive. Explicitly,
\[
\omega_i=\frac{\prod_{j=1}^{m-1}(a_i-\chi_j)}
 {a_i\prod_{k\ne i}(a_i-a_k)},\qquad
X_i=\frac{\omega_i}{\sum_k\omega_k}.                               \tag{9}
\]

*Proof.* The numerator and denominator of each \(\omega_i\) have \(i-1\)
negative factors, so \(\omega_i>0\). Partial fractions give the claimed zeros.
For uniqueness, multiply \(G_X(t)\) by \(\prod_i(t-a_i)\). Its numerator has
degree at most \(m-1\) and the \(m-1\) specified distinct zeros. It is therefore
a scalar multiple of their product. Residues determine every \(X_i\);
normalization fixes the scalar. The same argument for a homogeneous perturbation
shows the square matrix with rows \(\mathbf1^{\mathsf T}\) and
\((a_i/(a_i-\chi_j))_i\) is nonsingular. ∎

Apply the lemma to
\[
(\chi_1,\ldots,\chi_{m-1})=(r,\theta_2,\ldots,\theta_{m-2},s).        \tag{10}
\]
Thus the feed matching in this construction is a nonsingular **linear** system,
not an unproved transversality assumption.

We need four normalized pinch vectors. Put \(c=a_m\), \(a=a_1\). The top vectors
\(e_c,e_u\), and bottom vectors \(e_v,e_a\), are as follows:
\[
\begin{array}{lll}
(e_c)_i=\dfrac{c d_i}{L_T(a_i-c)}&(i<m),&
(e_c)_m=1-\sum_{i<m}(e_c)_i,\\
(e_u)_i=\dfrac{u d_i}{L_T(a_i-u)}&(i<m),&(e_u)_m=0,\\
(e_v)_i=\dfrac{v b_i}{L_B(v-a_i)}&(i>1),&(e_v)_1=0,\\
(e_a)_i=\dfrac{a b_i}{L_B(a-a_i)}&(i>1),&
(e_a)_1=1-\sum_{i>1}(e_a)_i .
\end{array}                                                       \tag{11}
\]
The root equations normalize \(e_u,e_v\). Direct substitution shows
\[
M_T e_c=V_Tc e_c,\quad M_T e_u=V_Tu e_u,\quad
M_B e_v=V_Bv e_v,\quad M_B e_a=V_Ba e_a,                            \tag{12}
\]
where \(M_T=V_T\operatorname{diag}(a)-da^{\mathsf T}\) and
\(M_B=V_B\operatorname{diag}(a)+ba^{\mathsf T}\).
The vectors \(e_u,e_v\) are positive on their indicated faces.
At \(\eta=0\), \(e_c=e_a\) is the positive vector in Lemma 1 with
\(\chi=\theta\), as follows from (2) and (13) below.
Hence \(e_c,e_a\) are strictly positive for small positive \(\eta\).

For any normalized pinch \(e\), including a pinch at the volatility of a
zero-net-flow component, the identity is
\[
\Psi(t)=\left(L-\frac{Vt}{S(e)}\right)G_e(t)
       =-\frac{V}{S(e)}(t-\phi)G_e(t),\quad \phi=LS(e)/V.           \tag{13}
\]
Indeed \(w_i=(Va_i/S(e)-L)e_i\), and
\(\sum_i a_i^2 e_i/(a_i-t)=S(e)+tG_e(t)\). This proof does not divide by
\(a_i-\phi\); there is no extra term when \(\phi=a_h,w_h=0,e_h>0\).
Evaluations at poles mean rational identities or their continuous extensions.

**Lemma 2 (strict positive cones).** The vector \(X\) from (9)–(10) has
representations
\[
X=Ae_c+Be_u=Ce_v+De_a,\quad A+B=C+D=1,\quad A,B,C,D>0.              \tag{14}
\]

*Proof.* Both top vectors vanish under \(G(\,r\,)\) and every
\(G(\theta_j)\), \(2\le j\le m-2\), by (7), (12), (13).
On the interval \((c,a_{m-1})\), \(G_{e_c}\) is strictly increasing with zero
\(u\), so \(G_{e_c}(s)<0\). All components of \(e_u\) have volatility larger
than \(s\), so \(G_{e_u}(s)>0\). Their unique convex combination with \(G(s)=0\)
has strictly positive weights and, by Lemma 1, equals \(X\).
For the bottom vectors, both vanish under \(G(s)\) and the middle-root
functionals. On \((a_2,a)\), \(G_{e_a}\) is increasing and has zero \(v<r\);
hence \(G_{e_a}(r)>0\). All nonzero components of \(e_v\) have volatility below
\(r\), so \(G_{e_v}(r)<0\). The same argument gives the second combination. ∎

## 4. Explicit finite-stage profiles

**Theorem 1 (physical isolated prefractionator, all \(m\), all real quality).**
Fix \(m\ge3,a,f,U\) as in Section 1. There exist \(\eta_0,C_0>0\) and
\(0<\rho<1\), depending on these fixed data, with the following properties.
For every \(0<\eta\le\eta_0\), flows (6) are positive and have an adiabatic
feed stage and finite equilibrium-stage sections of arbitrary lengths. Their
end-return compositions can be chosen so that every stage composition is
strictly positive, their net products are exactly \(d,b\), and the heavy
component in the upper liquid return and light component in the lower return's
equivalent liquid composition obey
\[
x_{n,m}\le C_0\eta^{-1}\rho^n,\qquad
x_{-k,1}\le C_0\eta^{-1}\rho^k.                                   \tag{15}
\]
Here \(n,k\ge1\) count cross-section steps from the feed-stage liquid \(X\).
Each section then has \(n-1\), respectively \(k-1\), interior trays between
the feed stage and its boundary return. Adding or deleting this fixed indexing
offset has no asymptotic effect.
The actual lower vapor return has composition \(y(x_{-k})\);
its light fraction is at most \((a_1/a_m)x_{-k,1}\), so the same exponential
bound holds after enlarging \(C_0\).

*Proof.* Use the four positive weights of Lemma 2 and define, for every
integer \(j\ge0\),
\[
\begin{aligned}
x_j&=\frac{A(c/u)^j e_c+B e_u}{A(c/u)^j+B},\\
x_{-j}&=\frac{C e_v+D(v/a)^j e_a}{C+D(v/a)^j}.                     \tag{16}
\end{aligned}
\]
Both expressions equal \(X\) at zero. They are convex combinations of normalized
vectors; the full-support vector has positive weight at every finite \(j\).
All coordinates are therefore strictly positive. Equations (12) show
that the first sequence is the normalized forward orbit of \(M_T\), and the
second the normalized backward orbit of \(M_B\). Consequently (1) gives
\[
V_T y(x_j)-L_Tx_{j+1}=d,\qquad
V_B y(x_{-(j+1)})-L_Bx_{-j}=-b.                                  \tag{17}
\]
Every interior tray balance follows by taking adjacent differences.
At the feed stage,
\[
L_Tx_1+V_By(x_{-1})+f-V_Ty(X)-L_BX=-d-b+f=0 .
\]
Its CMO vapor balance is \(V_T-V_B=U\); it requires no external heat.
Thus these are exact finite-stage physical profiles, with boundary streams
provided by their indicated positive compositions.

The two extreme traces satisfy the more explicit bounds
\[
x_{j,m}\le (A/B)(c/u)^j,\qquad x_{-j,1}\le(D/C)(v/a)^j.             \tag{18}
\]
All parameters are analytic near \(\eta=0\); the internal roots are fixed and
the other roots remain simple. At zero \(A=D=1,B=C=0\).
Furthermore
\[
B=\frac{-G_{e_c}(s)}{G_{e_u}(s)-G_{e_c}(s)}.
\]
Since \(G_{e_c}(u)=0\), \(u'(0)>0>s'(0)\), and \(G'_{e_c}(\theta_{m-1})>0\),
the numerator has a strictly positive first derivative at zero, and its
denominator has a strictly positive limit. Hence \(B'(0)>0\).
Likewise
\[
C=\frac{G_{e_a}(r)}{G_{e_a}(r)-G_{e_v}(r)}
\]
has \(C'(0)>0\), since \(r'(0)>0>v'(0)\).
After decreasing \(\eta_0\), \(A/B,D/C\le C_0/\eta\).
At zero \(c/u=c/\theta_{m-1}<1\) and \(v/a=\theta_1/a<1\);
continuity bounds both by a single \(\rho<1\). This proves (15). ∎

In particular, for fixed data, a sufficient section length is
\[
\max(n,k)\ \hbox{chosen with both }n,k\ge
\frac{\log(C_0/(\eta\varepsilon))}{-\log\rho},
\]
rounded upwards and replaced by one if the right side is smaller. This controls
the double limit \(\eta\downarrow0,\varepsilon\downarrow0\) for the isolated
prefractionator. Constants are not claimed uniform over coalescing volatilities,
vanishing feed fractions, or unbounded feed quality.

### 4.1 Ternary matching of prescribed return ratios

For \(m=3\), more is possible without changing any of the prescribed net flows.
Keep a fixed \(0<\eta\le\eta_0\). Write \(a>b>c\), and fix any
\(0<\xi,\zeta<1\).

**Theorem 2 (ternary return-ratio matching).** For all sufficiently large
integers \(n,k\) there is an exact, strictly positive finite-stage
prefractionator with the flows of (6), upper-return ratio
\(x_{n,1}/(x_{n,1}+x_{n,2})=\xi\), lower-return ratio
\(x_{-k,2}/(x_{-k,2}+x_{-k,3})=\zeta\), and
\[
x_{n,3}=O((c/u)^n),\qquad x_{-k,1}=O((v/a)^k).
\]
The feed-stage composition tends to \(X\) from Lemma 2. Bounds and the required
stage threshold may depend on \(\eta,\xi,\zeta\); for fixed \(\eta\) they are
uniform when \(\xi,\zeta\) range in compact subsets of \((0,1)\).

*Proof.* Put \(T=M_T/V_T\), \(B_0=M_B/V_B\), and let
\(f_U=(1-\xi,-\xi,0)\), \(f_D=(0,1-\zeta,-\zeta)\).
Solve the three-by-three **linear** system
\[
\mathbf1^{\mathsf T}X_{n,k}=1,\quad
f_UT^nX_{n,k}=0,\quad f_DB_0^{-k}X_{n,k}=0.                        \tag{21}
\]
The eigenvalues of \(T\) are \(r>u>c>0\); those of \(B_0\) are \(a>v>s>0\).
The top eigenvector at \(r\) is supported on \(AB\) and has opposite signs
in those two coordinates. Thus \(f_Ue_r\ne0\) for every \(0<\xi<1\).
Likewise the bottom eigenvector at \(s\) is supported on \(BC\), has opposite
signs there, and \(f_De_s\ne0\).

Scale the second row of (21) by \(r^{-n}\) and the third by \(s^k\).
By their spectral expansions these rows converge respectively to nonzero
multiples of the left eigenvectors at \(r\) and \(s\).
Those left eigenvectors are \(G(\,r\,)\) and \(G(\,s\,)\):
direct multiplication of the row \((a_i/(a_i-t))_i\) by \(M\), using
\(\Psi(t)=0\), gives eigenvalue \(Vt\). The limiting normalized system is
therefore Lemma 1. It is invertible, so (21) has a unique solution for all
sufficiently large \(n,k\), and
\[
X_{n,k}=X+O((u/r)^n+(s/v)^k).                                     \tag{22}
\]
In particular \(X_{n,k}>0\). Decompose it in the top eigenbasis with
coefficient of \(e_u\) tending to \(B>0\), and coefficient of \(e_c\)
tending to \(A>0\).
Equation (21) determines the \(e_r\) coefficient after multiplication by
\((r/u)^n\); hence \(u^{-n}T^nX_{n,k}\) has a limit \(W_U\) on the \(AB\)
face, with \(f_UW_U=0\).
Its value under the left functional \(G(\,u\,)\) is
\(B G_{e_u}(u)>0\). That functional is positive on every nonzero positive
vector on the \(AB\) face, so
\(W_U\) is a **positive** multiple of \((\xi,1-\xi,0)\).
The third coordinate before normalization is exactly the positive
\(e_c\) coefficient times \(c^n(e_c)_3\). Therefore the normalized upper
endpoint is strictly positive for large \(n,k\) and its third coordinate is
\(O((c/u)^n)\).

In the bottom eigenbasis, the \(e_v\) and \(e_a\) coefficients tend to
\(C,D>0\). The same argument applied to \(v^kB_0^{-k}X_{n,k}\) gives a
positive multiple of \((0,\zeta,1-\zeta)\): here \(G(\,v\,)\) is negative
on every nonzero positive vector on the \(BC\) face, and its limiting value
is \(C G_{e_v}(v)<0\), so the multiplier is again positive.
The remaining light coordinate is positive and \(O((v/a)^k)\).

Finally, (19) maps the positive upper endpoint backwards through every top
stage, and \(M_B\ge0\), with positive diagonal, maps the positive lower
endpoint forwards through every bottom stage. Thus **every** intermediate
composition is positive, not just the endpoints. The balance verification is
exactly (17) and the feed balance following it, with \(X_{n,k}\) in place
of \(X\). All limiting coefficients and inverse matrices are uniformly
bounded on the stated compact sets, which proves the uniformity assertion. ∎

This theorem fixes two major-component ratios and determines the small return
traces. It does **not** allow four independent return-composition coordinates
to be specified. The latter distinction is exactly relevant to closing the
ternary Petlyuk network.

## 5. Two useful exact cautions for a network proof

**No positivity from the forward matrix alone.** In general \(M_T\) has negative
entries. The positive-cone proof above is essential. One elementary safe
direction is also available: if \(w=d\ge0\), \(L=V-\sum d>0\), then
\[
M_T^{-1}=\frac1V\operatorname{diag}(a)^{-1}
 +\frac1{VL}\operatorname{diag}(a)^{-1}d\mathbf1^{\mathsf T}\ge0.    \tag{19}
\]
This follows by direct multiplication or the rank-one inverse formula.
For \(w=-b\le0\), \(M_B\ge0\) directly. These facts allow safe continuation
from the ends in the appropriate direction; they do not solve junction matching.

**A specific rank deficiency if net flows are frozen.** For the ternary
topology let \(F,U,D,A,B,C\) denote the feed, upper coupling, lower coupling,
light-product, middle-product and heavy-product liquid compositions.
Consider a putative sharp limit with fixed ideal net component flows and
simple spectra, and select the face-separating pinch branches. The twelve
endpoint equations, two from each of the six sections, read
\[
\begin{array}{c|cc}
\text{section}&\text{bottom-end equation}&\text{top-end equation}\\ \hline
\text{prefractionator top}&G_F(r)=0&U_C=0\\
\text{prefractionator bottom}&D_A=0&G_F(s)=0\\
\text{outer top}&\text{none}&A_B=A_C=0\\
\text{upper middle}&B_A=0&U_C=0\\
\text{lower middle}&D_A=0&B_C=0\\
\text{outer bottom}&C_A=C_B=0&\text{none}.
\end{array}                                                       \tag{20}
\]
On six two-dimensional composition simplices this matrix has rank ten:
the two \(F\) rows are independent by Lemma 1, there is one row on each
of \(U,D\), and two on each of \(A,B,C\). The duplicate rows are explicit.
The free directions are the major-component compositions within the \(AB\)
and \(BC\) coupling faces.

This calculation concerns **frozen net component flows and the indicated
branch only**. It is not a proof that every enlarged, rescaled network
formulation is singular. It does disprove the claim that these particular
endpoint equations automatically have an invertible junction-composition
Jacobian merely because the individual section eigenvalues are simple.
Some choices of stage-length ratios instead approach other pinch branches.

Agent S's Lemma 12 supplies a smooth exponential section chart only when its
pivot coefficients stay nonzero and the other chart coefficients stay bounded.
That assumption fails at a transition between two pinches. For example, (14)
has \(B,C\to0\) as slack tends to zero; the constants in (18) grow as \(1/\eta\).
More generally, spectral elimination gives a product of two small coefficients
equal to an exponentially small parameter times two bounded coefficients.
Its zero-parameter limit is a union of branches, not automatically a smooth
manifold. This is consistent with S's later warning in Section 8.4.

## 6. What remains unproved

Theorem 1 gives an **open boundary-value subassembly**. The two returning
streams have been chosen by construction. In a closed FTC they must be the
actual streams from the adjoining equilibrium junctions. Also, the trace
components returning into a binary finishing column cannot simply be discarded:
they affect the small component fluxes to its product ends and to the shared
middle-product stage.

One possible completion is a proof for suitably scaled trace-flow variables
that simultaneously establishes:

1. equality of the physical returning compositions at every thermal coupling;
2. all finite-stage component flows and compositions are nonnegative;
3. impurity tending to zero in **every external product**, including middle
   products;
4. global vapor and liquid balances with only the allowed FTC external duties,
   and total vapor at most \(\max H_j+o(1)\).

No uniform inverse, contraction, degree argument, or invariant region proving
these four assertions has been established here, even for the full ternary
network. In particular, it would be incorrect to replace this step by
“apply the implicit function theorem” or by the existence of some steady state
at prescribed total flows. The latter need not control purity.

Combining Agent S's binary theorem with the lower bound proves the requested
attainment statement for \(N=2\), all real \(q\). This work adds Theorem 1 for
every local \(m\), but **does not extend the closed-network theorem to \(N=3\)**.
None of the numerical observations below changes that logical status.

## 7. Reproduction and numerical evidence

Run under the project's WSL scientific environment, from this folder:

```text
python scripts/check_exact.py
python scripts/check_prefractionator.py
python scripts/check_ternary_returns.py
python scripts/check_ftc_numerical.py
```

Set `OMP_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=2`, `MKL_NUM_THREADS=2`,
and `PYTHONDONTWRITEBYTECODE=1`. The runner `scripts/run_checks.sh` also sets
these and a 3 GiB virtual-memory limit. All output files go in this folder.

`check_exact.py` checks rational partial fractions, positivity and eigenvector
identities, the zero-total redistribution, the pinch flux identity, Cauchy
interpolation, and the rank-ten example. Exact sample checks run for
\(m=3,\ldots,8\). The proofs above are general; the finite collection of checks
is not their replacement.

`check_prefractionator.py` checks 25 explicit-profile examples at 70-digit
precision: \(m=3,\ldots,7\), qualities \(-0.5,0,0.5,1,1.5\), up to 65 section
steps. It checks positivity, normalization, every cross-section balance, the
feed-stage balance, and (18). These are **numerical** checks, not interval
certificates. The largest balance residual in the recorded run is below
\(2.0\,10^{-70}\). Logs are in `runs/`.

`check_ternary_returns.py` checks (21) at 120-digit precision for several
different prescribed return ratios and section lengths. It checks every
interior stage, rather than treating endpoint positivity as sufficient.

The closed ternary experiment in `check_ftc_numerical.py` is separately labelled
numerical. It uses Agent P's finite-stage solver read-only, supplies explicit
phase flows with positive vapor slack, and rechecks all component balances,
normalizations, product impurities and the stage-level lower bound.
The quantitative results and limitations are in `runs/check_ftc_numerical.log`.
For \(R=H_{\max}+0.06\), the maximum impurity decreases from
0.00684513 at 20 stages per section to 0.000121556 at 40 stages. For
\(R=H_{\max}+0.03\), the corresponding values are 0.00926765 and 0.000257189.
All four cases pass the stated numerical acceptance checks.

**Failed extended attempts are retained.** The initial 80-stage runs did not
pass: at excess 0.06 the absolute balance residual was \(1.54\,10^{-9}\) and
the solver reported failure; at excess 0.03 the iteration diverged
(absolute residual 2.35). Their apparent product purities are not accepted
evidence. Both records, the numerical warnings, and the FAIL outcomes remain
in `runs/ftc_exploratory.json` and `runs/check_ftc_exploratory.log`.
Run `python scripts/check_ftc_numerical.py --extended` to attempt these cases;
that optional command may exit nonzero. The default runner checks the four
accepted 20/40-stage cases and does not claim to validate the failed cases.
Automatic warning paths were removed from the saved log; warning categories
and messages were retained.

For comparison, Agent P's `../simulation/runs/attain_T61_FTC.json`, as read in
this session, gives \(H_{\max}=1.365722513569909\). At prescribed impurity
\(0.01\), its 20- and 40-stage runs have respectively
\(R=1.3725713067027712\) and \(1.343556140716847\), with successful independent
re-solves. The 80-stage row explicitly reports `resolve_ok=false`; it must not
be silently described as independently re-solved. Values below the pure-product
limit at fixed nonzero impurity are compatible with the established finite-purity
lower bound and do not refute attainment.

## 8. Sources and novelty scope

- Halvorsen–Skogestad, Parts 1–3: the existing project source review and
  preferred-split theorem were read; Part 2 Section 6.2 and the Part 3 open
  passage were checked in the local full texts. Part 2,
  [doi:10.1021/ie0108649](https://doi.org/10.1021/ie0108649), explicitly bases
  its minimum-energy result on the infinite-stage analytical model.
- Halvorsen, *Minimum energy requirements in complex distillation arrangements*,
  NTNU thesis 2001:43: local Section 3.4.4 read for pinch compositions and the
  zero-net-flow limit. The singular formula there is avoided by (11)–(13).
- Project sources: `../../paper.tex` (preferred-split algebra),
  `../foundations/class.md`, `../sections/sections.md`,
  `../network/PROGRESS.md`, and `../simulation/` (as identified above).
- The project's existing 2026-10-20 novelty searches, including the independent
  verifier and Agent F, were read rather than represented as new searches of
  every database. Additional web queries in this session on finite-stage
  attainment, the conjecture, and steady-state existence/uniqueness found
  classical and numerical discussions, not a complete finite-stage FTC
  attainment proof. This is a limited search conclusion, not a proof of novelty.
- Halvorsen–Skogestad, *Energy efficient distillation*,
  [doi:10.1016/j.jngse.2011.06.002](https://doi.org/10.1016/j.jngse.2011.06.002):
  author-copy search extract read; it discusses the infinite-stage target and
  practical finite-stage energy. No finite-stage existence theorem is invoked
  from it.

No new paper requiring a paywalled PDF was used in a proof. No publication,
commit, push, or contact with authors was made. All claims await independent
verification under the project protocol.
