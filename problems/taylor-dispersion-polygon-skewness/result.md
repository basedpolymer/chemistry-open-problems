# Sign of the short- and long-time skewness of Taylor dispersion in regular polygonal ducts

Solver: Claude Opus 5.5 (effort xhigh). Status: **solved (preliminary)** for the
long-time and short-time signs (computer-assisted proofs, not yet verified); **partial result**
for the "all times" claim.

Throughout, `P_n` is the regular n-gon with **apothem (inradius) 1**, centred at 0, with side 0 equal to
`{x = 1, |y| <= t}`, where

    t = tan(pi/n),   R = sec(pi/n) = sqrt(1+t^2) (circumradius),   |P_n| = n t,   |dP_n| = 2 n t.

`<f> = |P_n|^{-1} int_{P_n} f dA` is the cross-sectional average. All signs below are invariant under dilation
(the paper's Table 1 uses circumradius 1; its section 3 uses inradius 1).

## 0. Results

| | statement | proof |
|---|---|---|
| Prop. 1 | long-time skewness `Sk(tau) = S_inf tau^{-1/2} + O(tau^{-3/2})`, `S_inf = 3 Pe^3 <u g2> / (2 + 2 Pe^2 <u g1>)^{3/2}` (the paper's eq. (17)) | analytic (Sec. 2) |
| Lemma 1 | `<u g2> = 2 <u g1^2>`; hence `sign S_inf = sign I(Omega)`, `I(Omega) := int u g1^2`, for every Pe > 0 | analytic |
| **Thm A** (n = 3) | `<u g1> = 9/1540`, `<u g1^2> = -7101/476476000 < 0`: the long-time skewness of the equilateral triangle is **negative** for every Pe; `3<u g2>/(2<u g1>)^{3/2} = -789 sqrt(770)/309400 = -0.0707623...`; `g1` is an explicit polynomial of degree 6 | exact (sympy) |
| **Thm B** (n >= 4) | `I(P_n) > 0` for **every** n >= 4: the long-time skewness of every other regular polygon is **positive** for every Pe | computer-assisted: ball arithmetic for 4 <= n <= 40, uniform radial perturbation argument for all n >= 30 |
| **Thm C** (short time) | geometric skewness `S^G_n = <u^3>/<u^2>^{3/2} > 0` for **every** n >= 3 (n = 3: `4 sqrt(21)/99`, exact); moreover `<u^3> = t^4/360 + O(t^5)` (explicit two-sided bounds), so `S^G_n ~ (8 sqrt(3)/15) tan^4(pi/n) ~ 90.0/n^4` | computer-assisted: 4 <= n <= 150 and all n >= 140 |
| Thm D (finite times, partial) | explicit certified windows `(0, t0]` and `[T1, inf)` on which the sign of `Sk(tau)` is proved, for n = 3, 4, 5, 6, 8 (Section 7) | computer-assisted |

Intermediate times `t0 < tau < T1` were left open at this stage; they are settled in `all-times/` (positive at all times for
every n >= 4, exactly one sign change for n = 3).

## 1. The source and the open statement

M. Aminian, R. Camassa, R. M. McLaughlin, *Mass distribution and skewness for passive scalar transport in pipes
with polygonal and smooth cross sections*, Stud. Appl. Math. 141(3) (2018) 399–417, doi:10.1111/sapm.12230;
arXiv:1708.05097 (v1 16 Aug 2017, v2 18 Aug 2017, v3 5 Feb 2018). Full texts of v1–v3 and of the published
version were read.

* Published abstract (p. 399): "For the case of the triangle and localized, cross-wise uniform initial
  distributions, short-time skewness is calculated exactly to be positive, while long-time asymptotics shows it to
  be negative. [...] The equilateral triangle appears to be the only regular polygon with this property—all others
  possess positive skewness at all times."
* arXiv abstract (identical in v1, v2, v3): "Interestingly, the equilateral triangle is the only regular polygon
  with this property, all other polygons possess positive skewness at all times, although this cannot cannot be
  proved on finite times due to the lack of closed form flow solutions for such geometries."
* Introduction (p. 401; arXiv v2/v3 p. 2): "We apply this analysis, using finite element solvers, and establish
  remarkably that the only polygonal domain yielding negative long-time solute negative skewness is the triangle."
  (v1: "...the only polygonal domain yielding long time solute skewness is the triangle.")
* Section 4 (p. 406): "We use a finite element package to calculate the asymptotic coefficients, and show that in
  fact the triangle is the only one of the regular polygons for which this occurs." The computation is
  Mathematica 10 FEM with `MaxCellMeasure -> 10^-3`, `MeshOrder -> 2` (p. 406); Table 1 (p. 407) lists
  `n = 3, 4, 5, 6, 8, 20, inf` (arXiv v3 also `n = 7, 100`).
* Conclusion (p. 416): "the triangle stands apart from all other regular polygons in that the concentration
  distribution's skewness approaches zero through negative values, despite its short-time sign behavior agreeing
  with all other regular polygons."

So the sign claims for n >= 4 (long and short time) rest on finite-element numerics for finitely many n, and
the finite-time claim is explicitly flagged as unproved. Differences between versions: the published abstract
replaces "is the only" by "appears to be the only" and drops the "cannot be proved" clause; the published
eq. (4) has the correct normalisation 1/6 (arXiv eq. (3) has `1/(12a)`, which gives `Lap u = -1`); the published
version adds exact triangle asymptotics (eq. (14): `M3 ~ (27/19250) Pe^3 t^3 + (9/175) Pe^3 t^4`) and a
point-source section; Table 1 is rounded and the columns n = 7, 100 are dropped. The long-time formula
(eqs. (15)–(17) published = (11)–(13) arXiv) is unchanged. Eq. (12) (published) still prints
`<u^3> = 1/19250`, `<u^2> = 3/700`, which correspond to `u/3`; the ratio 0.185 is unaffected.

Novelty check (references.md): no later work proves or disproves these signs; the only later paper on regular
polygonal ducts (Feng–Chu, arXiv:2605.22982, 2026) treats the Taylor coefficient of Brownian rods, not skewness.

## 2. The skewness functionals

**Setting.** `Omega` a bounded convex polygon, `ut` the Poiseuille profile `-Lap ut = 2`, `ut = 0` on `dOmega`
(any positive multiple of `ut` gives the same signs; the paper uses this normalisation), `u = ut - <ut>` (frame
moving with the mean speed). The tracer solves `T_tau + Pe u(y,z) T_x = Lap_3 T` in `R x Omega`, `d_n T = 0`,
with cross-sectionally uniform initial data `T(x,y,z,0) = rho(x)`, `rho >= 0`, `int rho = 1`, with finite moments.
As in the source we use the Aris moment hierarchy: `c_p(y,z,tau) = int x^p T dx` solves

    d_tau c_p = Lap c_p + p(p-1) c_{p-2} + p Pe u c_{p-1},   d_n c_p = 0,   c_p(0) = rho_p := int x^p rho,

(`Lap` = Neumann Laplacian of `Omega` from now on), and `M_p = <c_p>`; the skewness of the cross-sectionally
averaged distribution is `Sk = mu_3/mu_2^{3/2}` with central moments `mu_k`. Translating x we may assume
`rho_1 = 0`. Let `mu_2(Omega) > 0` be the first nonzero Neumann eigenvalue, `e^{tau Lap}` the Neumann heat
semigroup; it is self-adjoint, preserves means, and `||e^{tau Lap} v|| <= e^{-mu_2 tau} ||v||` if `<v> = 0`.

Cell problems (paper eqs. (15), (16)): `g1, g2` are the mean-zero weak solutions of the Neumann problems
`-Lap g1 = u`, `-Lap g2 = 2(u g1 - <u g1>)` (they exist since the data have mean zero).

**Proposition 1 (long time).** `M_1 = 0`, `mu_2(tau) = (2 + 2 Pe^2 <u g1>) tau + O(1)` and
`mu_3(tau) = 3 Pe^3 <u g2> tau + O(1)`; hence `tau^{1/2} Sk(tau) -> S_inf = 3 Pe^3 <u g2>/(2 + 2 Pe^2 <u g1>)^{3/2}`.

*Proof.* `c_0 = 1`. `c_1 = Pe A(tau)` with `A(tau) = int_0^tau e^{s Lap} u ds = g1 - e^{tau Lap} g1` (spectral
calculus: `u = sum u_i phi_i`, `g1 = sum (u_i/lambda_i) phi_i`); `M_1 = 0` since `<g1> = 0`. Then
`M_2' = 2 + 2 Pe <u c_1> = 2 + 2 Pe^2 <u g1> - 2 Pe^2 <u e^{tau Lap} g1>`, and the last term is
`O(e^{-mu_2 tau})`. Next `c_2 = rho_2 + 2 tau + Pe^2 S` with `S' = Lap S + 2 u A`, `S(0) = 0`, and
`M_3' = 6 M_1 + 3 Pe <u c_2> = 3 Pe^3 <u S>` (because `<u> = 0`). By Duhamel and self-adjointness,
`<u S(sigma)> = 2 int_0^sigma <(e^{(sigma-s)Lap} u) u A(s)> ds`, and inserting `A = g1 - e^{s Lap} g1`,

    <u S(sigma)> = <u g2> - <u e^{sigma Lap} g2> - 2 int_0^sigma <(e^{(sigma-s)Lap} u) u (e^{s Lap} g1)> ds,     (2.1)

whose last two terms are bounded by `e^{-mu_2 sigma}(||u|| ||g2|| + 2 sigma ||u|| ||u||_inf ||g1||)` (norms
w.r.t. the average). Integrating, `M_3 = rho_3 + 3 Pe^3 <u g2> tau + O(1)`. Since `M_1 = 0`, `mu_2 = M_2`,
`mu_3 = M_3`. []

**Lemma 1.** `<u g2> = 2 <u g1^2>` and `<u g1> = <|grad g1|^2> > 0`.

*Proof.* Weak forms: `int grad g1 . grad phi = int u phi`, `int grad g2 . grad phi = 2 int (u g1 - <u g1>) phi`
for all `phi in H^1`. Taking `phi = g2`, then `phi = g1`: `int u g2 = int grad g1 . grad g2 = 2 int u g1^2 - 2<u g1> int g1
= 2 int u g1^2`. []

Hence, **for every Pe > 0, the long-time sign of the skewness is the sign of `I(Omega) := int_Omega u g1^2`**, and
the paper's long-time coefficient is `kappa := 3<u g2>/(2<u g1>)^{3/2} = 6<u g1^2>/(2<u g1>)^{3/2}` (the
`Pe -> inf` limit of `S_inf`; it is invariant under `u -> c u` and scales like length). Check on the unit disk:
`u = 1/4 - r^2/2`, `g1 = r^4/32 - r^2/16 + 1/48`, `<u g1> = 1/768`, `<u g1^2> = 1/184320`, `kappa = sqrt(6)/10
= 0.24495` (Table 1: 0.245).

**Proposition 2 / 3 (short and all times).** With `rho_1 = 0`: `mu_3(tau) = kappa_3(rho) + 3 Pe^3 K(tau)`,

    K(tau) := int_0^tau <u S(sigma)> dsigma = 2 int_0^tau int_0^sigma <(e^{(sigma-s)Lap} u) u A(s)> ds dsigma,

and `mu_2(tau) > 0`. So for a symmetric longitudinal profile (e.g. a delta release, `kappa_3(rho) = 0`) the sign
of `Sk(tau)` equals the sign of `K(tau)` **for every Pe and every tau**. Lemma 11 (Section 7) gives
`K(tau) = <u^3> tau^3/3 + O(tau^{7/2})`, so the short-time sign is the sign of `<u^3>`, i.e. of the geometric
skewness `S^G = <u^3>/<u^2>^{3/2}` of the PRL 2015 paper (paper eq. (12)).

## 3. The equilateral triangle (Theorem A)

Triangle with inradius 1 (paper eq. (4)): sides `y = -1`, `y = 2 +- sqrt(3) z`,
`ut = (1+y)((2-y)^2 - 3z^2)/6`, `<ut> = 3/10`. The cell problem has the **polynomial** solution

    g1 = y^6/960 - y^5/96 - y^4 z^2/64 + y^4/32 + y^3 z^2/48 + y^2 z^4/64 + y^2 z^2/16 - 11 y^2/120
         + y z^4/32 - z^6/960 + z^4/32 - 11 z^2/120 - <.>,

(`-Lap g1 = u` identically and `d_n g1 = 0` on the three sides, checked symbolically; a degree-<= 5 ansatz has no
solution, so 6 is minimal). Then `<u g1> = 9/1540`, `<u g1^2> = -7101/476476000`, so by Lemma 1
`<u g2> = -7101/238238000 < 0` and `kappa_3 = -789 sqrt(770)/309400 = -0.0707623189...` (paper eq. (18):
`-0.07076`; circumradius 1: `-0.0353812`, Table 1: `-0.035`). Also `<u^2> = 27/700`, `<u^3> = 27/19250`,
`S^G = 4 sqrt(21)/99 = 0.185155`. (`g2` itself is not a polynomial of degree <= 12; it is not needed.)
Script: `scripts/triangle_exact.py`.

## 4. A-posteriori toolkit (n >= 4)

Notation: `||.||` is the L2(P) norm unless stated; `||.||_dP` the L2(dP) norm; `sd(phi) = ||phi - <phi>||`.

**Lemma 2 (Poiseuille error, maximum principle).** If `ut_N` is a polynomial with `-Lap ut_N = 2`, then
`E := ut - ut_N` is harmonic with `E = -ut_N` on `dP`; hence `min_dP(-ut_N) <= E <= max_dP(-ut_N)`,
`|E| <= delta := max_dP |ut_N|`, `|E - <E>| <= osc := osc_dP(ut_N)`.
Moreover `0 <= ut <= R^2/2` (comparison with the circumscribed disk).

**Lemma 3 (L2 bound for harmonic functions).** For `E` harmonic in `P_n`, `E in H^2(P_n)`:
`||E||^2 <= ((1+R)/2) ||E||_dP^2`.

*Proof.* `ut in H^2(P_n)` (convex polygon; Grisvard 1985, Thm 3.2.1.2). Green's identities give
`2 int E^2 = int E^2 (-Lap ut) = -2 int ut |grad E|^2 + int_dP E^2 (-d_n ut) <= int_dP E^2 (-d_n ut)`, using
`ut >= 0`, `Lap(E^2) = 2|grad E|^2`, `ut = 0` on `dP`. On side 0, `ut <= v := (1-x)(x+R)` (maximum principle,
`P_n subset {-R <= x <= 1}`, `-Lap v = 2`); `ut - v` is harmonic, vanishes on the open side, so (Schwarz
reflection) `ut` is smooth up to the open side and `0 <= -d_n ut <= -d_n v = 1 + R` there. By symmetry this holds
on every side. []

In applications `E = ut - ut_N` (or `E = ut - Y`, Section 6), which is in `H^2` because `ut` is.

**Lemma 4 (Neumann error).** Let `g_N` be a polynomial with `-Lap g_N = u_N := ut_N - <ut_N>` exactly,
`<g_N> = 0`, and `e := d_n g_N` on `dP`. Then `w := g1 - g_N` satisfies
`||grad w|| <= eta := p ||E|| + tau_tr ||e||_dP` and `||w|| <= omega := p eta`, where
`p = 2R/pi` and `tau_tr = sqrt(2p^2 + 2Rp)`.

*Proof.* Weak form: `int grad w . grad phi = int f phi - int_dP e phi` for all `phi in H^1`, `f := u - u_N = E - <E>`,
`<w> = 0`. Poincaré–Wirtinger for convex domains (Payne–Weinberger 1960; Bebendorf 2003, Thm 3.2):
`||w|| <= (diam/pi) ||grad w|| <= p ||grad w||` since `diam P_n <= 2R`. Trace identity (apothem 1, so
`x . n = 1` on `dP`): `||w||_dP^2 = int_dP w^2 (x.n) = int_P (2 w^2 + 2 w x.grad w) <= 2||w||^2 + 2R ||w|| ||grad w||
<= tau_tr^2 ||grad w||^2`. With `phi = w`: `||grad w||^2 <= ||f|| ||w|| + ||e||_dP ||w||_dP <= (p||f|| + tau_tr||e||_dP)||grad w||`,
and `||f|| <= ||E||`. []

**Lemma 5 (perturbation of I).** With `I_N := int u_N g_N^2`, `q := u_N g_N - <u_N g_N>`,
`Uinf := max(<ut_N> + delta, R^2/2 + delta - <ut_N>) >= sup_P |u_N|`:

    |I - I_N| <= ||E|| sd(g_N^2) + osc (2 ||g_N|| omega + omega^2) + 2 sd(u_N g_N) omega + Uinf omega^2 =: Err.

*Proof.* `I - I_N = int f g^2 + int u_N (g^2 - g_N^2)`, `g = g_N + w`. Since `int f = 0`,
`int f g_N^2 = int f (g_N^2 - <g_N^2>)`; `|f| <= osc`; `int u_N (g^2 - g_N^2) = 2 int q w + int u_N w^2`
(`<w> = 0`). []

**Lemma 6 (exact integrals).** For `F_d` homogeneous of degree d, on the half-sector `T = {0 <= y <= t x <= t}`:
`int_T F_d dA = (d+2)^{-1} int_0^t F_d(1, tau) dtau` (cone coordinates `(x,y) = (lam, lam tau)`, Jacobian `lam`).
For `D_n`-invariant `F`, `int_P F = 2n int_T F`; the restriction to `dP` is determined by `F(1, tau)`,
`0 <= tau <= t`, and the outward normal derivative on side 0 is `d_x F_d(1,tau) = d F_d(1,tau) - tau d_tau F_d(1,tau)`
(Euler). For radial `phi(r^2)` this gives `<r^{2k}>_P = (k+1)^{-1} sum_j C(k,j) t^{2j}/(2j+1)`: averages of radial
polynomials are polynomials in t alone.

**Lemma 7 (Rellich bound).** For `E` harmonic in `P_n`, `E in H^2`: `||d_n E||_dP <= (t + R) ||d_s E||_dP`.

*Proof.* The field `Z = (x.grad E) grad E - |grad E|^2 x/2` is divergence-free in 2D when `Lap E = 0`
(`div Z = (x.grad E) Lap E`). On a side, `x.n = 1` and `x.grad E = d_n E + tau d_s E` with `|tau| <= t`
(`tau` = arclength from the side midpoint). So `0 = int_dP Z.n = int_dP ((d_nE)^2/2 + tau d_sE d_nE - (d_sE)^2/2)`,
i.e. `X^2 <= Y^2 + 2t XY` with `X = ||d_nE||_dP`, `Y = ||d_sE||_dP`, whence `X <= (t + sqrt(1+t^2)) Y`. []

## 5. Long-time sign for n >= 4 (Theorem B)

**Finite part, 4 <= n <= 40** (`scripts/lt_certify.py`). Approximants
`ut_N = (1-r^2)/2 + sum_{k=0}^K a_k Re((w/S)^{nk})` and
`g_N = particular solution + sum_{j=1}^J b_j Re((w/S)^{nj}) - const` (w = x + iy; `S` = binary approximation of
R; K = J = 16 (n <= 5), 10 (n <= 8), 6 (n <= 14), 4 (n >= 15)); the coefficients come from floating-point least
squares on the side and are then exact binary numbers (stored in `certificates/lt_certificates_4_40.json`).
`-Lap ut_N = 2` and `-Lap g_N = u_N` hold *exactly* (`Lap(r^{2a} Re w^m) = 4a(a+m) r^{2a-2} Re w^m`; checked
term by term). Every quantity of Lemmas 2–5 is computed in arb ball arithmetic (python-flint, 248–500 bits here and 352–660 bits in Section 6, chosen from the degree):
boundary sup/oscillation by Taylor enclosures on 240 subintervals of `[0, t]`, all integrals exactly through
Lemma 6. Result: `I_N - Err > 0` for every `4 <= n <= 40`, with `Err/I_N` between `0.1%` (n = 4) and `7.9%`
(n = 15). Example (apothem 1): n = 4: `I >= 5.0306e-5` (I_N = 5.0352e-5, Err <= 4.6e-8); n = 40:
`I >= 1.6953e-5`.

**Uniform part, all n >= 30** (`scripts/lt_uniform.py`). Take the radial approximants
`ut_N = (1-r^2)/2 + t^2/6`, `u_N = 1/4 + t^2/12 - r^2/2`, `g_N = r^4/32 - (1/16 + t^2/48) r^2 + b_0(t)`. Then
(exactly, sympy): on the side `ut_N = t^2/6 - tau^2/2`, so `delta = t^2/3`, `osc = t^2/2`,
`||ut_N||_dP^2 = 2nt^5/45`; `e = (3 tau^2 - t^2)/24`, `||e||_dP^2 = n t^5/360`; `sup_P |u_N| = 1/4 + 5t^2/12`;
and by Lemma 6 every term of Lemma 5 divided by `|P| = nt` is an explicit function of t only (polynomials and
square roots), e.g. `<u_N g_N^2> = 1/184320 + t^2/110592 + 17 t^4/4147200 + O(t^6)`. The inequality
`F(t) := <u_N g_N^2> - Err/(nt) > 0` is checked in arb interval arithmetic on 400 subintervals covering
`[0, tan(pi/30)]` (minimum of the lower bound `8.1957e-7`, against `F(0) = 1/184320 = 5.43e-6`). Since
`t = tan(pi/n) <= tan(pi/30)` for n >= 30, `I(P_n) > 0` for all n >= 30. Together: **Theorem B** for all n >= 4.

Enclosures of the paper's coefficient (circumradius 1, `scripts/table1.py`), see Section 8.

## 6. Short-time sign (Theorem C)

**6.1 Finite part, 4 <= n <= 150** (`scripts/sg_certify.py`). With `ut_N`, `u_N`, `E`, `f = E - <E>` as above,
`phi = u_N^2 - <u_N^2>`, and a polynomial `Phi` with `-Lap Phi = phi` exactly plus a fitted harmonic part
(so that `rho := Phi|_dP` is small):

    int u^3 = int u_N^3 + 3 int phi E + 3 int u_N f^2 + int f^3,
    int phi E = int_dP rho d_nE - int_dP E d_nPhi          (Green; E = -ut_N on dP; int_dP d_nE = 0),

so with Lemmas 2, 3, 7:

    int u^3 >= int u_N^3 + 3 int_dP ut_N d_nPhi - 3 ||rho||_dP (t+R) ||d_s ut_N||_dP - (3 Uinf + osc) ((1+R)/2) ||ut_N||_dP^2.

(Lemma 8.) All terms are computed in ball arithmetic (K = 16 for n <= 6, 10 for n <= 12, 6 beyond). Result: the
lower bound is positive for every `4 <= n <= 150`; the relative error `(int u_N^3 - lower bound)/int u_N^3` is at most 1.1% (n = 150) (log
`certificates/run_all.log`, data `certificates/sg_certificates_4_150.json`).

**6.2 Uniform part, all n >= 140** (`scripts/sg_uniform.py`). Write `ut = Y + H` with the radial
`Y = (1-r^2)/2 + t^2/6` and `H` harmonic, `H = b(tau) := tau^2/2 - t^2/6` on each half side;
`Y_c = Y - <Y> = 1/4 + t^2/12 - r^2/2`, `H_c = H - <H>`, `phi = Y_c^2 - <Y_c^2>`, `Phi` radial with
`-Lap Phi = phi`, `rho = Phi|_dP - mean`. Exactly,

    int u^3 = [int Y_c^3 - 3 int_dP b d_nPhi] + 3 int_dP rho d_nH + 3 int Y_c H_c^2 + int H_c^3,

and the bracket equals `n t M(t)` with (exact) **`M(t) = t^4/360 + 19 t^6/7560`**. Bounds:
`|3 int_dP rho d_nH| <= 3 ||rho||_dP (t+R) ||d_s b||_dP` (Lemma 7; `||d_s b||_dP^2 = 2nt^3/3`, `int_0^t rho^2 =
t^9/4082400 + ...`), `|3 int Y_c H_c^2 + int H_c^3| <= (3(1/4 + 5t^2/12) + t^2/2) int H^2`, and (Lemma 10)
`int H^2 <= n t (vhat + ehat)^2` where, with the explicit harmonic function `V = sum_{k<=6} beta_k Re(w^{nk})`,
`beta_k = 2t^2(-1)^k/(k^2 pi^2)` (the Fourier coefficients of b):

* (10a) `||V|| <= sum |beta_k| (int_P r^{2nk})^{1/2}`, `int_P r^{2nk}/(nt) <= (t/(pi k)) exp(k pi t/(1 - t^2/3))`
  (Lemma 6, `log(1+t^2) <= t^2`, `t - t^3/3 <= arctan t = pi/n <= t`);
* (10b) `||H - V|| <= sqrt((1+R)/2) ||b - V||_dP` (Lemma 3), with
  `||b - V||_{L2(0,t)} <= t^{5/2} sqrt(2/(3 pi^4 K^3)) + sqrt(t) sum_k |beta_k| d_k(t)`,
  `d_k(t) = k pi t^2/(3 - t^2) + exp(k pi t/(2(1 - t^2/3))) - 1 >= sup |cos(k pi tau/t) - Re((1+i tau)^{nk})|`
  (Fourier tail `sum_{k>K} k^{-4} <= 1/(3K^3)`; `|tau/t - arctan(tau)/arctan(t)| <= t^2/(3 - t^2)`).

Dividing by `n t^5`, `G(t) := M(t)/t^4 - (error terms)/t^4` is an explicit continuous function on `[0, 1]`
with `G(0) = 1/360`; arb interval arithmetic on 400 subintervals gives `G > 2.45e-4 > 0` on `(0, tan(pi/140)]`.
Hence `<u^3> > 0` for all n >= 140, and with 6.1, **`S^G_n > 0` for all n >= 4** (n = 3: Section 3).

*Asymptotics.* The bounds are two-sided. With K = 6 fixed the remainder is `<= 5e-5 t^4 + O(t^5)` (about 2% of the main term `t^4/360`); choosing
`K = K(t) ~ t^{-2/3}` gives `ehat = O(t^3 log(1/t))`, `vhat = O(t^{5/2})`, hence rigorously
`<u^3>_{P_n} = t^4/360 + O(t^5)`. Similarly `<u^2> = <Y_c^2> + O(t^{5/2}) = 1/48 + O(t^2)`, so
`S^G_n = (8 sqrt(3)/15) tan^4(pi/n) (1 + O(1/n))`, i.e. `S^G_n ~ 90.0/n^4` (arXiv Table 1, n = 100: `8.6e-7`;
certified enclosure `[8.348e-7, 8.397e-7]`, Section 8).

## 7. Finite times: certified windows (Theorem D, partial result)

`scripts/alltime_windows.py`. For symmetric cross-sectionally uniform initial data and every Pe,
`sign Sk(tau) = sign K(tau)` (Proposition 3). Let `c0 = sup_{x>0}(1-e^{-x})/sqrt(x) <= 0.64` (checked in ball arithmetic; the supremum is 0.6382 at x = 1.2564).

**Lemma 11 (short times).** `|K(tau) - <u^3> tau^3/3| <= C_a tau^{7/2} + C_b tau^4`,
`C_a = (32/105) c0 ||grad u|| ||u^2||`, `C_b = (pi/48) c0^2 ||u||_inf ||grad u||^2` (average norms).
*Proof.* Spectrally, `||(e^{r Lap} - 1)u||^2 = sum (1 - e^{-r lambda_i})^2 u_i^2 <= c0^2 r sum lambda_i u_i^2
= c0^2 r ||grad u||^2`. Write `e^{r Lap}u = u + rho_r` and `A(s) = s u + R_s`, `R_s = int_0^s rho_r dr`,
`||R_s|| <= (2/3) c0 s^{3/2} ||grad u||`. Then
`K(tau) = <u^3> tau^3/3 + 2 int_0^tau int_0^sigma [ s<rho_{sigma-s} u^2> + <u^2 R_s> + <rho_{sigma-s} u R_s> ] ds dsigma`
and, by Cauchy–Schwarz, the three terms are bounded by `s c0 (sigma-s)^{1/2} ||grad u|| ||u^2||`,
`(2/3) c0 s^{3/2} ||grad u|| ||u^2||` and `(2/3) c0^2 (sigma-s)^{1/2} s^{3/2} ||u||_inf ||grad u||^2`. The double
integrals are `2 int int s (sigma-s)^{1/2} = 16 tau^{7/2}/105`, `2 int int (2/3) s^{3/2} = 16 tau^{7/2}/105`,
`2 int int (2/3)(sigma-s)^{1/2} s^{3/2} = (4/3)(pi/64) tau^4` (Beta integrals `B(2,3/2) = 4/15`, `B(3/2,5/2) = pi/16`). []

**Lemma 12 (long times).** From (2.1): `|K(tau) - <u g2> tau| <= int_0^inf e^{-mu_2 sigma}(a1 + 2 sigma b1) dsigma
= a1/mu_2 + 2 b1/mu_2^2`, `a1 = ||u|| ||g2||` (spectrally `||g2|| <= ||2(u g1 - <u g1>)||/mu_2 = 2 sd(u g1)/mu_2`),
`b1 = ||u|| ||u||_inf ||g1||`. Hence `sign K(tau) = sign <u g2>` for `tau > T1 := (a1/mu_2 + 2 b1/mu_2^2)/|<u g2>|`.

With certified enclosures of all norms (exact for n = 3) and `mu_2 = pi^2/4` (square, exact) or the
Payne–Weinberger bound `pi^2/diam^2` (others), the script gives (t0 rounded down, T1 rounded up):

| n | `Sk(tau) > 0` for `0 < tau <= t0` | sign of `Sk(tau)` for `tau >= T1` |
|---|---|---|
| 3 | `t0 = 3.279e-3` | **negative**, `T1 = 371` |
| 4 | `t0 = 5.299e-4` | positive, `T1 = 32.0` |
| 5 | `t0 = 1.377e-4` | positive, `T1 = 41.5` |
| 6 | `t0 = 4.383e-5` | positive, `T1 = 34.2` |
| 8 | `t0 = 6.725e-6` | positive, `T1 = 24.3` |

(`t0` shrinks with n because `<u^3>` does, `S^G_n ~ 90/n^4`; `T1` is driven by the crude constant of Lemma 12.)

So for the triangle the sign change predicted by the paper is certified to occur inside `(t0, T1)`; for n >= 4
positivity is certified on `(0, t0] U [T1, inf)` only. Times are in units of `a^2/kappa` with `a` = apothem.

## 8. Table 1 of the paper, rigorous enclosures

Circumradius-1 normalisation as in the paper (`kappa` scales like length: `kappa(R=1) = cos(pi/n) kappa(apothem 1)`);
printed values from arXiv v3 (the published Table 1 rounds them and drops n = 7, 100). Enclosures from
`scripts/table1.py` (outward-rounded).

| n | LT coefficient, paper | certified enclosure | `S^G`, paper | certified enclosure |
|---|---|---|---|---|
| 3 | -0.035 | `-789 sqrt(770)/618800 = -0.0353812` (exact) | 0.185 | `4 sqrt(21)/99 = 0.1851546` (exact) |
| 4 | 0.139 | [0.13849, 0.13856] | 0.0813 | [0.081294, 0.081297] |
| 5 | 0.194 | [0.19390, 0.19429] | 0.042 | [0.042262, 0.042270] |
| 6 | 0.215 | [0.21474, 0.21603] | 0.024 | [0.024264, 0.024274] |
| 7 | 0.225 | [0.22384, 0.22699] | 0.0149 | [0.014939, 0.014955] |
| 8 | 0.231 | [0.22903, 0.23281] | 0.0097 | [0.0097036, 0.0097168] |
| 20 | 0.243 | [0.23937, 0.24655] | 3.97e-4 | [3.9611e-4, 3.9832e-4] |
| 100 | 0.245 | [0.24384, 0.24589] | 8.6e-7 | **[8.3475e-7, 8.3971e-7]** |
| inf | 0.245 | `sqrt(6)/10 = 0.2449490` (exact) | 0 | 0 |

All printed values agree with the enclosures up to the printed rounding, except the arXiv-only value
`S^G_100 = 8.6e-7`, which is about 3% too large (it does not appear in the published table). The asymptotic law
`(8 sqrt(3)/15) tan^4(pi/100) (1 + O(1/n))` gives `9.0e-7 (1 - O(0.03))`, consistent with the enclosure. The
LT enclosures for n >= 6 are wide only because the bound on `||grad(g1 - g_N)||` is crude; they are not needed
for the sign.

## 9. Certificates

* `certificates/lt_certificates_4_40.json`: for each n, K, J, the coefficients `a_k`, `b_j` and the scale `S` as
  hexadecimal floats (exact binary numbers), the certified `delta`, `||ut_N||_dP`, `||E||`, `||e||_dP`, `eta`,
  `omega`, `I_N`, the four terms of `Err`, the lower bound of `I`, and enclosures of `kappa`.
* `certificates/sg_certificates_4_150.json`: for each n, K, J, coefficients `a_k` (Poiseuille) and `c_j`
  (harmonic part of `Phi`), the boundary term, the Rellich and second-order terms, the lower bound of `int u^3`
  and an enclosure of `S^G`.
* `certificates/run_all.log`: complete output of all scripts (the run used for this report).
Re-running a script recomputes everything from scratch (the least-squares fits are deterministic); a verifier
may also feed the stored coefficients directly into the ball-arithmetic part.

## 10. Assumptions, trusted base, limitations

* **Model.** Results are about the functionals defined by the Aris hierarchy for cross-sectionally uniform
  initial data, exactly as in the source; Proposition 1 makes the long-time asymptotics rigorous at that level
  (standard semigroup facts). Unimodality/median statements of the paper are not addressed.
* **Cited theorems.** H^2 regularity of the Dirichlet problem on convex polygons (Grisvard 1985, Thm 3.2.1.2,
  statement recalled, not re-read); Payne–Weinberger/Bebendorf Poincaré constant `diam/pi` (Bebendorf 2003,
  Thm 3.2, read); Schwarz reflection; Green's identities for `H^2` functions on Lipschitz domains; spectral
  theorem for the Neumann Laplacian. The Rellich identity is used only for `H^2` harmonic functions, where it
  follows from the divergence theorem.
* **Computations.** Trusted: python-flint 0.9.0 / Arb (ball arithmetic, polynomial arithmetic), sympy 1.14
  (exact rational/polynomial identities). Floating point (numpy) is used only to *choose* coefficients.
* **Not proved.** The sign of the skewness for intermediate times (`t0 < tau < T1`) for any n >= 4, i.e. the
  paper's "positive skewness at all times"; the monotonicity in n claimed for Table 1.
