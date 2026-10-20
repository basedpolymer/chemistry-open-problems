# The Halvorsen–Skogestad conjecture for general N (Agent G)

> **Verifier's note (independent verification, see `../verification.md`).** (1) The N = 5 counts of §7–8 (1456
> sign-hard pairs per j, 28 open (basic, j) pairs, 1152 (variant, j) pairs) come from an earlier version of
> `theory/scripts/signcheck.py`; the final version gives 1312 sign-hard pairs per j and an open set of 20 (basic, j)
> pairs = 864 (variant, j) pairs (`theory/scripts/opencount.py`, equal to Agent C's path enumeration), as in the
> paper. (2) The open case is not true in general: for the pair (basic 161, j = 3) the permissive model has a feasible
> point with R < V_FTC (all flows positive, so (N) holds), certified by interval arithmetic
> (`../verification/krawczyk_certificate.py`). The proofs of G1, G3, G4, G5, G8 and of the strict-model statement G2
> were checked and are correct (G2 does not actually need (N)).

Status of this file: version 1 (2026-10-20). Model, notation and the cut identity are those of
`theory/setting.md` and `theory/proofs.md` (Agent T) and of `computation/PROGRESS.md` (Agent C). Everything
stated as a Lemma/Theorem below has a complete proof here, under the explicit assumption (N) of Section 1.
Floating-point scripts in `general/scripts/` only test the statements; they are not part of the proofs
(except the enumeration counts of Section 7, which are exact combinatorics).

Summary of what is proved for every number of components N >= 3:

* **Lemma G1 (recovery ordering).** In every column of every regular-column configuration, at every
  TG-feasible point (excess vapor allowed), the top recoveries d_p/f_p are non-increasing in p.
* **Theorem G2 (STRICT model).** If every column for which j is an interior Underwood root operates at its
  threshold (no excess vapor), then R >= H_j. In the strict infinite-stage model (Agent C) this holds for every
  j, i.e. the conjecture is true for all N in that model. (Same statement as Agent C's sign-graph theorem;
  the proof below is independent, in "flux" form, and shows that only the excess of the *type-(d)* columns of
  the split j matters.)
* **Theorem G3 (PERMISSIVE/TG model, merge-free splits).** If the mixed sub-network M_j has no merge node
  (no mixed stream produced by two splits), then R >= H_j at every TG-feasible point, excess vapor allowed.
  Hence the conjecture holds, for all N, in the full TG model, for every configuration in which M_{j*} is
  merge-free for the index j* of the highest peak; and for every configuration, for the peaks j whose M_j is
  merge-free (always the case for j = 1 and j = N-1).
* **Theorem G6 (hybrid).** At a TG-feasible point, R >= H_j holds as soon as no merge node of M_j descends
  from a type-(d) column that carries excess vapor. (Contains G2 and G3.)
* **Theorem G8 (simple diamonds, TG model).** R >= H_j at every TG-feasible point if every merge node of M_j
  is the bottom of a simple FTC-type diamond (type-(d) column -> type-(b) and type-(c) children -> merge).
  With G3 this gives, for every N, the conjecture in the full TG model for all configurations whose mixed
  networks are merge-free or made of simple diamonds, and an alternative proof of the N = 4 case.

What remains open for N >= 5 in the TG (permissive) model: splits j whose mixed network contains a merge node
that is not the bottom of a simple diamond, below a type-(d) column with positive excess, when the sign
relaxation of Agent T does not already close the case (Section 8; for N = 5 these are 28 (basic
configuration, j) pairs, 1152 (variant, j) pairs; adversarial numerics find only ties there). Degenerate points
(a stream lacking one of its components) are excluded throughout (assumption (N), Section 9).

---------------------------------------------------------------------------------------------------------

## 1. Model, assumption (N), regions

Data: N components, a_1 > ... > a_N > 0 (constant relative volatilities), feed F = 1, z_p > 0, any real q.
I_r = (a_{r+1}, a_r). theta_j = feed root in I_j. P_j = sum_{p<=j} a_p z_p/(a_p - theta_j),
H_j = P_j - (1-q) = sum_{p>j} a_p z_p/(theta_j - a_p). V_FTC = max_j H_j (total reboiler vapor of the FTC).

Configurations, transfer streams (TC / HX / side draw), and the TG constraints are those of
`theory/setting.md` §2–3. For a column c (stream [i,m], split [i,k]/[l,m], l <= k+1): top net flows d_p >= 0
(p <= k), bottom net flows b_p >= 0 (p >= l), feed flows f = d + b, net vapor feed U_c, vapors V_T - V_B = U_c,
threshold vapors Y_T - Y_B = U_c, 0 <= Y_T <= V_T; for each r in {l-1,...,k} a root theta_{c,r} in I_r of
sum_p a_p f_p/(a_p - theta) = U_c with g_T(theta_{c,r}) <= Y_T, and equality for the *interior* roots
l <= r <= k-1. The excess e_c := V_T - Y_T = V_B - Y_B >= 0.

* PERMISSIVE model = TG literal (excess allowed). STRICT model (Agent C) = PERMISSIVE plus e_c = 0 for every
  column having at least one interior root. STRICT ⊂ PERMISSIVE.

**Assumption (N) (non-degeneracy).** Every stream carries a positive net flow of each of its components:
for every column c and every p in its stream [i,m], f_{c,p} > 0.
(Under (N) every feed root of every column lies in the open interval I_r, and every function used below
has its poles exactly where stated.) The degenerate points (some f_{c,p} = 0) are not treated here; see §9.

**Fluxes** (Agent T's Lam, Agent C's Psi). For a section with vapor V and upward net flows w,
Lam(theta) = V - sum_p a_p w_p/(a_p - theta). For a column c and theta in I_j:
  T_c(theta) := Lam_{c,T}(theta) = V_T - sum_{p<=k} a_p d_p/(a_p - theta),
  B_c(theta) := -Lam_{c,B}(theta) = -V_B - sum_{p>=l} a_p b_p/(a_p - theta),
  a_c(theta) := T_c + B_c = U_c - sum_p a_p f_p/(a_p - theta)   (the column's feed flux).
On I_j all three are non-increasing (strictly if a flow is positive). Write every such function as
  y(theta) = U_y - G_y(theta) + K_y(theta),
  G_y(theta) = sum_{p<=j} a_p w_p/(a_p - theta)   (light part: >= 0, non-decreasing on I_j),
  K_y(theta) = sum_{p>j}  a_p w_p/(theta - a_p)   (heavy part: >= 0, non-increasing on I_j),
with w the net flows carried (w = d for T_c, w = b for B_c, w = f for a_c).

**Regions for the split j.** L = streams inside [1,j] and light products; H = streams inside [j+1,N] and heavy
products; M = M_j = mixed streams [i,m], i <= j < m (always contains [1,N]). Each mixed column c is of one of
four types:
  (a) k = j, l = j+1: both products leave M (sharp j-splitter);
  (b) k <= j, l <= j: distillate in L ("light exit"), residue mixed;
  (c) k >= j+1, l >= j+1: distillate mixed, residue in H ("heavy exit");
  (d) k >= j+1, l <= j: both products mixed. Then l <= j <= k-1: **j is an interior root of c**.
Links of M: from a mixed column to a mixed product (TC, HX or side draw). A **merge node** is a mixed stream
produced by two splits (a mixed side draw; both links TC).

**Propagation.** For a mixed column c with producers, for every theta (TG definitions; `theory/proofs.md`):
  a_c(theta) = sum over its links e of y_e(theta),
  y_e = T_P - kappa_e for a distillate link from P (kappa_e = condenser duty L_T(P) if HX, else 0),
  y_e = B_P + rho_e  for a residue link from P   (rho_e = reboiler duty V_B(P) if HX, else 0),
and a_{[1,N]}(theta) = Lam_F(theta) = (1-q) - sum_p a_p z_p/(a_p - theta), which vanishes at theta_j.
Light flows are conserved: G_{a_c} = sum_e G_{y_e}, and G_{a_c} = G_{T_c} + G_{B_c}.

**Cut identity** (Agent T Prop. 1 = Agent C's identity for the canonical cut; proof in `theory/proofs.md`,
exact symbolic check in `computation/scripts/verify_identity_exact.py`):
  R - H_j = R(L) + Q_cond(H) + Q_cond(M) + sum_{light exits c} T_c(theta_j),               (*)
light exits = mixed columns whose distillate lies in L (types (a), (b)); Q_cond(M) = sum of the condenser
duties kappa_e on HX links into mixed streams. Hence

  R - H_j >= Tot_j := sum_{light exits c} T_c(theta_j) + Q_cond(M),                         (**)

and every result below is a proof of Tot_j >= 0.

## 2. Single-column facts

Let c be a mixed column. Under (N), a_c has exactly one zero theta^c in I_j (a_c -> +inf at a_{j+1}^+,
-inf at a_j^-). If j is in the root range {l-1,...,k} of c, theta^c = theta_{c,j}.

**Lemma G0.**
(i) (Underwood at theta^c) T_c(theta^c) = -B_c(theta^c) =: s_c >= 0 (the Underwood slack), and s_c = e_c if c
    is of type (d). [type (a),(b) with k = j, (c) with l = j+1: TG at the root theta_{c,j}; type (b) with
    k < j: T_c >= 0 on I_j because T_c is non-increasing on (-inf, a_k) and T_c(theta_{c,k}) >= 0 with
    theta_{c,k} > a_{k+1} >= a_j; type (c) with l > j+1: symmetric; type (d): equality at the interior root.]
(ii) (sign rule at every level) for every theta in I_j: T_c(theta) >= 0 or B_c(theta) <= 0. [T_c non-increasing,
    -B_c non-decreasing, both equal s_c >= 0 at theta^c.] In particular T_c(theta) >= min(0, a_c(theta)).
(iii) For a type-(d) column, T_c(theta) and B_c(theta) have opposite strict signs exactly on the "window"
    (z_2, z_1), where z_1 = zero of T_c, z_2 = zero of B_c, z_2 <= theta^c <= z_1, and z_1 = z_2 = theta^c iff
    e_c = 0.

## 3. Lemma G1 — recovery ordering (every column, TG model)

**Lemma G1.** Let c be any column (split [i,k]/[l,m]) at a TG-feasible point with f_{c,r} > 0 for all r in
[i,m]. Then for all p < q in [i,m]: r_p >= r_q, where r_p := d_p/f_p (with d_p := f_p for p < l and d_p := 0
for p > k). Equivalently d_p b_q >= d_q b_p.

*Proof.* Only p < q both in [l,k] needs a proof. Let D(theta) = Y_T - g_T(theta) (top threshold flux),
Btil(theta) = Y_B + sum_{r>=l} a_r b_r/(a_r - theta), F(theta) = U - sum_r a_r f_r/(a_r - theta) = D - Btil.
Facts: F vanishes at every feed root theta_r (r in [i,m-1]); D(theta_r) = Btil(theta_r) = 0 for interior r;
D(theta_r) = Btil(theta_r) >= 0 for r = l-1 and r = k (TG); also for r < l-1 and r > k (D non-increasing on
(-inf,a_k) and Btil non-decreasing on (a_l, +inf), compared with the outer roots).
Suppose r_p < r_q and pick lam with r_p < lam < r_q, lam distinct from the finitely many r_s (s in [l,k]).
Then 0 < lam < 1. Put H := (1-lam) D + lam Btil = D - lam F = K0 - sum_{s in [i,m]} a_s c_s/(a_s - theta),
c_s := d_s - lam f_s. Signs: c_s = (1-lam) f_s > 0 for s < l, c_s = -lam f_s < 0 for s > k, c_p < 0 < c_q,
all c_s != 0. So H = P(theta)/prod_s(a_s - theta) with deg P <= n := m - i + 1, P not identically 0, and
P(a_s) != 0 for all s.
Near a_s, H ~ -a_s c_s/(a_s - theta): on I_r = (a_{r+1}, a_r) the sign of H at the right end is -sgn c_r and at
the left end +sgn c_{r+1}. Let sigma_s = sgn c_s. The sequence (sigma_i,...,sigma_m) starts with +, ends
with -, and contains - (at p) before + (at q): it has C >= 3 sign changes, and every change sits at an
interval I_r with r in [l-1,k] (outside, the signs are constant).
* If sigma_r = sigma_{r+1}, the ends of I_r have opposite signs: >= 1 zero of P in I_r.
* If sigma_r != sigma_{r+1} (so r in [l-1,k]), both ends of I_r have the same sign, and I_r contains the root
  theta_r where H = 0 (r interior) or H >= 0 with both ends negative (r = l-1: ends -sigma_{l-1} = -,
  sigma_l = -; r = k: ends -sigma_k = -, sigma_{k+1} = -). In both cases P has >= 2 zeros in I_r counted with
  multiplicity (same sign at both ends and a zero, or a non-negative interior value between two negative ends).
Total: >= (n-1-C) + 2C = n-1+C >= n+2 zeros of a nonzero polynomial of degree <= n. Contradiction. QED

Notes. (1) No positivity of d, b, Y_T, Y_B is needed, only f > 0 on [i,m]. (2) For the root column of K this
is Agent T's Lemma 4. (3) Exploration check: LP minimum of r_p - r_{p+1} over the TG polytope of random
columns (feed roots fixed, the constraints are linear in (d, Y_T)): `general/scripts/recovery_order_lp.py`,
1587 LPs, minimum +9.7e-5, PASS. Also checked at NLP points by `gamma_check.py`.

## 4. Credit functionals

For a function y = U_y - G_y + K_y as above, a level mu in I_j, and the fixed root theta_j, define
  F_mu(y)  := G_y(mu) - G_y(theta_j) + min(0, y(mu)) = min(G_y(mu), U_y + K_y(mu)) - G_y(theta_j),
  Gamma(y) := G_y(zeta_y) - G_y(theta_j) = max_{mu in I_j} F_mu(y),   zeta_y := zero of y in I_j.
(The max: F_mu(y) = G_y(mu) - G_y(theta_j) is non-decreasing for mu <= zeta_y and F_mu(y) = U_y + K_y(mu) -
G_y(theta_j) is non-increasing for mu >= zeta_y.) Properties:
* F_mu is superadditive for a fixed mu: F_mu(y'+y'') >= F_mu(y') + F_mu(y'') (G linear, min(0,.) superadditive).
* F_{theta_j}(y) = min(0, y(theta_j)); Gamma(y) >= min(0, y(theta_j)).
* For the main column, Gamma(Lam_F) = F_{theta_j}(Lam_F) = 0 (zeta = theta_j).

**Lemma G4 (local inequalities).** Let c be a mixed column with feed flux a = a_c, links e1 (distillate, if
mixed) and e2 (residue, if mixed) carrying y1 = T_c - kappa, y2 = B_c + rho, and let ell = T_c(theta_j) if c
is a light exit (types a, b), ell = 0 otherwise. Then, with kappa counted only if e1 exists:
 (G4.1) [min form, level theta_j] ell + kappa + sum_{e} min(0, y_e(theta_j)) >= min(0, a(theta_j))
        for types (a), (b), (c), and for type (d) if e_c = 0 or theta_j is outside the window of c;
 (G4.2) [Gamma form] ell + kappa + sum_{e} Gamma(y_e) >= Gamma(a) for all four types, any excess.

*Proof.* Write G = G_a, G^T = G_{T_c}, G^B = G_{B_c}; G = G^T + G^B (light flows).
(G4.1). Since min(0, t-kappa) + kappa >= min(0,t) and min(0, b+rho) >= min(0,b), it suffices to treat
kappa = rho = 0 with t = T_c(theta_j), b = B_c(theta_j), a = t + b.
 (a) ell = t >= min(0,a) by G0(ii).
 (b) t + min(0,b) = min(t, a) >= min(0,a) if t >= 0; if t < 0 then b <= 0 (G0(ii)) and t + b = a.
 (c) min(0,t) >= min(0,a): if t >= 0 clear; if t < 0 then b <= 0 and t >= a.
 (d) min(0,t) + min(0,b) = min(0,a) when t, b have the same sign; by G0(iii) this holds when theta_j is not in
     the open window, in particular when e_c = 0 (empty window).
(G4.2). (a) Top = light flows only: T_c(theta) = V_T - G(theta), so ell = T_c(theta^c) + G(theta^c) - G(theta_j)
 = s_c + Gamma(a) >= Gamma(a).
 (b) Top light only, all heavy go down. Let z := zeta_{y2} (zero of B_c + rho). Then
 ell + Gamma(y2) - Gamma(a) = s_c - [G^B(theta^c) - G^B(z)]. If z >= theta^c the bracket is <= 0. If
 z < theta^c: Lam_B(theta^c) - Lam_B(z) = s_c - rho = [G^B(theta^c) - G^B(z)] + [K^B(z) - K^B(theta^c)], the
 second bracket >= 0, so s_c >= rho + [G^B(theta^c)-G^B(z)] >= the bracket.
 (c) All light flows go up: G_{y1} = G. With z := zeta_{y1}: Gamma(y1) + kappa - Gamma(a) = kappa + G(z) -
 G(theta^c). If z >= theta^c done. If z < theta^c: T_c(z) - T_c(theta^c) = kappa - s_c = [G(theta^c) - G(z)] +
 [K^T(z) - K^T(theta^c)] >= G(theta^c) - G(z), so kappa >= s_c + G(theta^c) - G(z).
 (d) Here T_c(theta^c) = e, B_c(theta^c) = -e (e = e_c). Let z1 = zeta_{y1} (T_c(z1) = kappa),
 z2 = zeta_{y2} (Lam_B(z2) = rho). The claim is
     X := [G^T(z1) - G^T(theta^c)] + kappa - [G^B(theta^c) - G^B(z2)] >= 0.
 Case kappa <= e, rho <= e: z1 >= theta^c >= z2 and
   e - kappa = T_c(theta^c) - T_c(z1) = (z1 - theta^c) * sum_p x_p,  x_p = a_p d_p /(|a_p - theta^c| |a_p - z1|),
   e - rho  = Lam_B(theta^c) - Lam_B(z2) = (theta^c - z2) * sum_p y_p, y_p = a_p b_p /(|a_p - theta^c| |a_p - z2|),
 and G^T(z1) - G^T(theta^c) = sigma1 (e-kappa), G^B(theta^c) - G^B(z2) = sigma2 (e-rho), where
 sigma1 = X_L/(X_L+X_H), sigma2 = Y_L/(Y_L+Y_H) (L: p <= j, H: p > j). By Lemma G5 below sigma1 >= sigma2,
 hence X = sigma1 e + (1-sigma1) kappa - sigma2 (e - rho) >= (sigma1 - sigma2) e >= 0.
 Case kappa > e: z1 < theta^c and T_c(z1) - T_c(theta^c) = kappa - e >= G^T(theta^c) - G^T(z1), so the first
 two terms of X are >= e; the last bracket is <= e - rho <= e if rho <= e and < 0 if rho > e.
 Case kappa <= e < rho: first bracket >= 0, last bracket < 0. QED

**Lemma G5 (share comparison).** In the setting of (G4.2)(d) with z1 >= theta^c >= z2: X_L Y_H >= X_H Y_L.
*Proof.* For p light, q heavy: X_p Y_q - X_q Y_p = (a_p a_q/(|a_p-theta^c||a_q-theta^c|)) * [v1_p v2_q d_p b_q -
v2_p v1_q d_q b_p] with v1_s = 1/|a_s - z1|, v2_s = 1/|a_s - z2|. For light p (a_p > z1 >= z2):
v1_p/v2_p = (a_p - z2)/(a_p - z1) >= 1; for heavy q (a_q < z2 <= z1): v1_q/v2_q = (z2-a_q)/(z1-a_q) <= 1. Hence
v1_p v2_q >= v2_p v1_q, and with d_p b_q >= d_q b_p (Lemma G1) the bracket is >= v2_p v1_q (d_p b_q - d_q b_p)
>= 0. Summing over p <= j < q gives the claim. QED

## 5. Theorem G2 (no excess at type-(d) columns; STRICT model)

**Theorem G2.** Under (N), at any TG-feasible point such that every type-(d) column of M_j has e_c = 0
(more generally: theta_j is outside its window), Tot_j >= 0, hence R >= H_j. In the STRICT model this holds
for every j: R >= V_FTC for every configuration, every N, every q.

*Proof.* Sum (G4.1) over all mixed columns. Each link into a mixed stream appears once on the left (at its
producer) as min(0, y_e(theta_j)); each mixed stream s != [1,N] has on the right min(0, a_s(theta_j)) =
min(0, sum_{e into s} y_e(theta_j)) >= sum_{e into s} min(0, y_e(theta_j)) (superadditivity; covers merge
nodes). Telescoping leaves sum_c ell_c + Q_cond(M) >= min(0, a_{[1,N]}(theta_j)) = 0. QED
(Agent C proved the STRICT statement independently by a sign-graph/cut argument; the two proofs certify
each other. The proof shows that the PERMISSIVE model differs from STRICT, for the split j, only through the
excess of the type-(d) columns of M_j.)

## 6. Theorem G3 (merge-free mixed networks, PERMISSIVE model) and Theorem G6

**Theorem G3.** Under (N), if M_j has no merge node, then Tot_j >= 0 (so R >= H_j) at every TG-feasible point,
for every TC/HX choice and every q.

*Proof.* M_j is a tree rooted at [1,N]; let L(c) := (light-exit fluxes T(theta_j) of the columns of the subtree
of c) + (condenser duties on links inside the subtree). By (G4.2) and induction from the leaves,
L(c) >= Gamma(a_c) for every mixed c (each child c' of c has a_{c'} = y_e for its unique link e). At the root,
Tot_j = L([1,N]) >= Gamma(Lam_F) = 0. QED

**Theorem G6 (hybrid).** Under (N), let D_+ be the set of type-(d) columns of M_j whose window contains
theta_j (e.g. those with e_c > 0 and a_c(theta_j) small). If no descendant of a column of D_+ is a merge node,
then Tot_j >= 0.
*Proof.* Let W = descendants of D_+ (excluded D_+ itself), U0 = M_j \ W. No column of W has a producer in U0
(it would be a merge node), so W is a forest whose roots are children of D_+ columns. Use (G4.1) on U0 \ D_+,
(G4.2) on W, and for c in D_+: ell + kappa + Gamma(y1) + Gamma(y2) >= Gamma(a_c) >= min(0, a_c(theta_j)). The
links from U0 \ D_+ to U0 are credited min(0, y(theta_j)) (merges inside U0 are handled by superadditivity as
in G2), the links from D_+ into W are credited Gamma(y). Telescoping as in G2 and G3. QED

### 6.1 Blocks, levels, and the diamond theorem G8

**Block identity.** Let B be a set of mixed columns with an entry column e such that every link into B goes to
e, and let its "child blocks" be the mixed streams C outside B fed by links from B, each of which receives *all*
its links from B (so a_C = sum of the links from B into C). Put
  L(B) := sum_{light exits in B} T(theta_j) + (condenser duties on links inside B and into its children)
          + sum_C Gamma(a_C).
For any level mu in I_j one has the two exact expressions (light form; heavy form, obtained from the light form
with the flux balance at level mu inside B and G(zeta)-G(mu) = a(mu) + K(zeta) - K(mu) for zeta = zero of a)
  L(B) - Gamma(a_e) = sum_{light exits} T(mu) + kappa_B + sum_C [G_C(zeta_C) - G_C(mu)] + [G_e(mu) - G_e(theta^e)]
                    = sum_{heavy exits} Lam_B(mu) + rho_B + sum_C [K_C(zeta_C) - K_C(mu)] + [K_e(mu) - K_e(theta^e)],
(zeta_C = theta^C, kappa_B / rho_B = condenser / reboiler duties on the links of B and into its children).
If the mixed network is partitioned into blocks forming a tree (each block's entry fed only by its parent
block) and every block satisfies L(B) >= Gamma(a_e), then Tot_j = L(root block) >= Gamma(Lam_F) = 0.

**Theorem G8 (simple diamonds).** Under (N), suppose every merge node s of M_j is the bottom of a simple
diamond: s is the residue of a type-(b) column c1 and the distillate of a type-(c) column c2, c1 and c2 are not
merge nodes, and c1, c2 are the distillate and the residue of one type-(d) column d. Then R >= H_j at every
TG-feasible point (excess vapor and all TC/HX choices allowed).
*Proof.* Blocks: {d, c1, c2} for each diamond (entry d, single child block s), singletons otherwise; this is a
block tree (an entry d that is itself a merge node is the bottom of the diamond above, whose block feeds it
entirely). Singletons satisfy the block inequality by (G4.2). For a diamond take mu = theta^d (so the last
brackets vanish). With kappa on d->c1 and rho on d->c2: a_{c1}(mu) = e_d - kappa and a_{c2}(mu) = rho - e_d.
* If theta^s >= mu (light form): T_{c1}(mu) + kappa >= 0 [if kappa <= e_d, theta^{c1} >= mu and T_{c1}(mu) >=
  s_{c1} >= 0; if kappa > e_d, G0(ii) gives T_{c1}(mu) >= a_{c1}(mu) = e_d - kappa]; G_s(theta^s) - G_s(mu) >= 0;
  d and c2 have no light exit. Hence L(B) >= Gamma(a_d).
* If theta^s <= mu (heavy form): Lam_{B,c2}(mu) + rho >= 0 [symmetric argument]; K_s(theta^s) - K_s(mu) >= 0;
  d and c1 have no heavy exit. Hence L(B) >= Gamma(a_d). QED

**Corollaries.** (1) j = 1 and j = N-1: M_j is a chain (all mixed streams are [1,m], resp. [i,N]), merge-free:
R >= H_1, R >= H_{N-1} in the TG model (= Agent T's Theorem A). (2) N = 3: the conjecture in the TG model.
(3) N = 4: at j = 2 the merge-free sign-hard family K is covered by G3, and the four basic configurations with
a merge at j = 2 (those containing ABC, BCD, BC) are simple diamonds, covered by G8. Hence G1 + G3 + G8 give an
alternative proof of Agent T's Theorem B (N = 4, TG model, all 152 configurations), independent of the
explicit computations of Agent T's Lemmas 4-5. (4) Any N: every configuration whose mixed networks M_j are, for
all j, merge-free or made of simple diamonds satisfies the conjecture in the TG model.

## 7. Structure of M_j and counts

**Proposition G7 (grid picture).** Identify the mixed stream [i,m] with the grid point (i,m) in
[1,j] x [j+1,N], and let P_j be the set of points of the configuration's mixed streams. By the regularity rule
(distillate of [i,m] = largest present [i,k], k < m; residue = smallest present [l,m], l > i), a mixed column's
mixed distillate is the next point of P_j below it in its grid column, and its mixed residue the next point
of P_j to its right in its grid row. Consequently: a point is a merge node iff it has a point of P_j above it in
its column and to its left in its row; type (d) iff it has a point below and a point to the right; light exits
are the lowest points of the grid columns, heavy exits the rightmost points of the rows. M_j is a planar DAG
(links go down or right), every non-root node has 1 or 2 producers, so its cycle rank is the number of merge
nodes. Every merge node s has a type-(d) ancestor: follow two root-to-producer paths, one through each
producer of s; at the first node where they part, both products are on the paths, hence mixed (type (d)).
(Informally: bounded faces of the planar drawing <-> merge nodes as bottom-right corners, with a type-(d)
top-left corner; not used in the proofs.)

Counts (`general/scripts/classify_mixed.py`, exact enumeration on top of `computation/scripts/hsmodel.py`):

| N | basic configs | (basic, j) pairs | merge-free (G3) | with a merge | of which simple diamonds (G8) |
|---|---|---|---|---|---|
| 4 | 18 | 54 | 50 | 4 | 4 |
| 5 | 203 | 812 | 640 | 172 | 88 |
| 6 | 4373 | 21865 | 13616 | 8249 | 2113 |

Agent T's sign relaxation (`theory/scripts/signcheck.py`, LP over sign patterns) applied to all 6128 N = 5
configurations: sign-hard (variant, j) pairs: 1456 at j = 2 and 1456 at j = 3, none at j = 1, 4. Of these,
880 + 880 are merge-free (closed by Theorem G3) and 576 + 576 have a merge node that is not a simple diamond
(none of the sign-hard pairs is covered by G8: the simple-diamond cases are already sign-easy at N = 5); they
come from 28 (basic configuration, j) pairs, listed by `classify_mixed.py 5 --sign` (e.g. basic 31, j = 2:
ABCDE->ABCD/BCDE; ABCD->ABC/BCD; BCDE->BCD/E; ABC->AB/C; BCD->B/CD; AB->A/B; CD->C/D, merge node BCD, whose
residue-producer ABCD is itself of type (d)). For N = 4 every sign-hard pair is merge-free (the 16 K variants).
(Agent C's permissive sign-path enumeration, a different relaxation, leaves 20 of these basic pairs.)

## 8. What is not proved (PERMISSIVE model with merges)

The two local mechanisms fail in complementary places:
* (G4.1) (sign rule at theta_j) fails at a type-(d) column whose window contains theta_j (positive excess), by
  at most e_c: ell + kappa + sum min(0,y_e) >= min(0,a) - e_c.
* (G4.2) (Gamma) is exact at every column but Gamma is NOT superadditive: at a merge node s with in-links y', y''
  the defect Gamma(y') + Gamma(y'') - Gamma(y'+y'') can be positive (physically: mixing a vapor-rich and a
  liquid-rich feed lowers the threshold vapor). Numerically 21 of 23 merge defects at N = 4 PERMISSIVE points
  were positive, up to 0.93 (`gamma_check.py`).
In simple diamonds the block identity of §6.1 at the level of the diamond's top closes the case (Theorem G8:
light form if the merge root is above that level, heavy form otherwise). The open cases are merges whose
two producer paths are longer or contain another type-(d) column (e.g. basic 31/72, j = 2: the residue
producer ABCD of the merge BCD is itself a dipole). There, a single level per block fails: in basic 72, j = 2,
with the root and ABCD carrying excess, the only admissible block is {ABCDE, ABCD, BCDE} and at its level both
the light and the heavy form have one term of the wrong sign when the merge root lies below theta_j. A
first-order expansion in the excess e of ABCD (root excess -> 0, all Underwood slacks 0) gives
  Tot_j ~ e * [ sigma(top of ABCD) - sigma(feed of BCD) ],
sigma = sensitivity-weighted light share; G1 only orders sigma(top of ABCD) >= sigma(bottom of ABCD), whereas
the feed of BCD also contains the distillate of BCDE. So a proof must use more than G1 and the sign rule:
presumably the composition constraint imposed on the distillate of BCDE by its outer Underwood root at j (it
forces enough C into that distillate). Numerically there is no counterexample: targeted adversarial
minimisation of R/H_2 jointly over (alpha, z, q) and the operating point (`general/scripts/adv_merge.py`,
12 starts each, PERMISSIVE model) gives 0.9999999974 (basic 72, all TC), 1.0000000130 (basic 31), 1.0000000022
(basic 124), 0.9999999988 (basic 194), 1.0000009 and 1.0000000088 (two HX variants of basic 72): ties to
solver accuracy, never below 1 - 1e-8 (log `general/runs/adv_merge_j2_a.log`). Precisely open:

> For N >= 5, the PERMISSIVE (TG literal) model: configurations and splits j such that M_j contains a merge node
> that is not the bottom of a simple diamond, below a type-(d) column carrying excess vapor, and the sign
> relaxation does not close the case. For N = 5: the 28 (basic, j) pairs of §7 (1152 (variant, j) pairs).
> For N = 6 not yet enumerated against the sign relaxation (8249 (basic, j) pairs have a merge, 2113 of them
> are simple diamonds).

In the STRICT model nothing remains open (Theorem G2), apart from the degenerate points of §9.

## 9. Degenerate points

All statements assume (N). At a point where some stream lacks one of its components, roots may sit at the
interval ends and some zeros used above may not exist. Agent T's Lemma 1 treats endpoint roots by noting that
the feed equation then forces the corresponding flow to vanish; the Gamma-machinery would need zeros defined
as endpoints (Gamma(y) = sup over the closed interval). Not written.

## 10. Files

* `general/scripts/recovery_order_lp.py` — LP test of Lemma G1 (PASS).
* `general/scripts/classify_mixed.py` — merge/type counts of M_j for N = 4, 5, 6; with `--sign`, intersection
  with Agent T's sign relaxation (exact enumeration + LP).
* `general/scripts/gamma_check.py` — NLP stress test (PERMISSIVE / STRICT points) of G1, the sign rule, (G4.1),
  (G4.2), Theorem G3 and the conjecture itself; prints PASS/FAIL. Results: N = 4, 40 trials (74 points), PASS;
  N = 5, 60 trials on basic 21, 31, 64, 194 (90 points), PASS (min local Gamma slack -1e-9; (G4.1) at type-(d)
  columns fails by up to 0.42, as expected with excess; merge defects positive in 103/120 cases, up to 6.2).
  STRICT model, N = 4, 30 trials (57 points): PASS, and (G4.1) holds also at type-(d) columns (min -4e-10),
  as Theorem G2 predicts.
* `general/scripts/nlp_g.py` — frozen copy of Agent C's NLP model (so that the checks do not depend on files
  being edited).
* `general/scripts/adv_merge.py` — adversarial minimisation of R/H_j for the open merge cases (exploration).
