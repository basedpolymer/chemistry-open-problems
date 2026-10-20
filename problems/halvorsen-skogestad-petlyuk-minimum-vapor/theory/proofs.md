# Proofs (Agent T) — lower bound R >= H_j in the TG model

**Superseded by paper.tex (10:45), which contains the corrected and extended versions (lemma numbering differs).**

Notation of theory/setting.md. Fix a split index j in {1,...,N-1} and write theta := theta_j (process feed
root in I_j), Lam_s := Lam_s(theta). "Mixed" = stream [i,m] with i <= j < m. All flows (d, b, f) are >= 0.

## Lemma 0 (monotonicity)
Let J be an open interval containing no a_p of the components present in the section. On J,
Lam_T(x) = V_T - sum_{p<=k} a_p d_p/(a_p - x) is nonincreasing, Lam_B(x) = V_B + sum_{p>=l} a_p b_p/(a_p - x)
is nondecreasing, Lam_F(x) = U - sum_p a_p f_p/(a_p - x) is nonincreasing, and strictly so as soon as one
flow is positive (derivatives -sum a d/(a-x)^2, +sum a b/(a-x)^2, -sum a f/(a-x)^2). The same holds on a
closed interval [a_{r+1}, a_r] at the endpoints where the expressions are finite. Lam_T - Lam_B = Lam_F.

## Lemma 1 (column sign rule)
Let c be a mixed column, a_c := Lam_{c,F}(theta), t_c := Lam_{c,T}(theta), b_c := -Lam_{c,B}(theta)
(so a_c = t_c + b_c). Then  a_c >= 0 => t_c >= 0  and  a_c <= 0 => b_c <= 0.
Proof. Split [i,k]/[l,m].
(S) l-1 <= j <= k: TG gives theta_{c,j} in [a_{j+1}, a_j] with Lam_{c,F}(theta_{c,j}) = 0 and
Lam_{c,T}(theta_{c,j}) >= 0, Lam_{c,B}(theta_{c,j}) >= 0 (g_T <= Y_T <= V_T, g_B <= Y_B <= V_B). The values are
finite (if theta_{c,j} is an endpoint a_r, the feed equation forces f_r = 0, hence d_r = b_r = 0). If the
column has no positive flow the claim is trivial; otherwise Lam_{c,F} is strictly decreasing on I_j, so
a_c >= 0 forces theta <= theta_{c,j} and then t_c = Lam_T(theta) >= Lam_T(theta_{c,j}) >= 0; a_c <= 0 forces
theta >= theta_{c,j} and Lam_B(theta) >= Lam_B(theta_{c,j}) >= 0.
(T) k < j (distillate inside [1,j]): the top section contains only p <= k, so Lam_T is nonincreasing on
(-inf, theta_{c,k}] (all its poles are >= a_k >= theta_{c,k}); theta < a_j <= a_{k+1} <= theta_{c,k}, hence
t_c >= Lam_T(theta_{c,k}) >= 0 unconditionally; if a_c <= 0 then Lam_B(theta) = t_c - a_c >= 0.
(B) l > j+1: symmetric (Lam_B(theta) >= Lam_B(theta_{c,l-1}) >= 0 unconditionally). QED

Inflow propagation (definitions of TG): for a mixed child c,
a_c = sum over its producers of [ t_P - kappa_P ] (distillate link; kappa_P = L_T(P) if HX, 0 if TC/side draw)
    + [ b_P' + rho_P' ] (residue link; rho_P' = V_B(P') if HX, 0 if TC/side draw),
and a_{[1,N]} = Lam_F^{feed}(theta_j) = 0. (Linear: Lam of the vapor transfer stream is D - sum a d/(a - x)
= Lam_T(P) - L_T(P); Lam of the liquid residue transfer is -sum a b/(a - x) = -Lam_B(P') + V_B(P').)

## Proposition 1 (cut identity; = Agent C's identity for the canonical cut)
Node regions L, M, H as in setting.md sec. 5. Every node conserves Lam(x) up to: +Q at a reboiler, -Q at a
condenser, -(-a_p z'_p/(a_p - x)) for a liquid product withdrawal of z'_p of component p, + Lam_F^{feed}(x) at
the feed stage (checked for: column feed stage, partial condenser and reboiler on a transfer stream, total
condenser / reboiler on a product, liquid side draw between two sections). Summing over M u H at x = theta_j:
   R - H_j = R(L) + Q_cond(H) + Q_cond(M) + sum_{light exits s} Lam_s.                      (*)
Light exits = top sections of mixed columns with distillate inside [1,j]. Q_cond(M) = sum of kappa over HX
distillate links into mixed streams.

## Theorem A (extreme peaks). For all N, every configuration, every TG-feasible point: R >= H_1, R >= H_{N-1}.
Proof (j=1). Mixed streams are [1,m], each the distillate of exactly one stream: a path
c_0 = [1,N] -> c_1 -> ... -> c_n, one light exit (top of c_n, distillate A). Lemma 1: t >= min(0,a).
a_{c_{i+1}} = t_{c_i} - kappa_i, so min(0,a_{c_{i+1}}) >= min(0,a_{c_i}) - kappa_i; from a_{c_0}=0,
t_{c_n} >= -sum kappa_i = -Q_cond(M). Insert in (*). j = N-1: mirror image (heavy side, residues). QED

## Lemma 2 (path). If every mixed column has at most one mixed child (then M is a path from the root,
links top or bottom, TC or HX), R >= H_j.
Proof. Backward induction: T_i := sum of t over light exits among c_i..c_n, K_i := sum of kappa on links
inside c_i..c_n; claim T_i >= min(0,a_i) - K_i. Top link: a_{i+1} = t_i - kappa_i >= min(0,a_i) - kappa_i.
Bottom link (c_i's top is a light exit): T_i = t_i + T_{i+1} >= t_i + min(0, b_i + rho_i) - K_{i+1}
>= min(0, a_i) - K_i in both cases b_i+rho_i >= 0 (use t_i >= min(0,a_i)) and < 0 (use t_i + b_i = a_i).
Then (*) with T_0 >= -K_0 = -Q_cond(M). QED

## Lemma 3 (N = 4, j = 2, diamond). Configurations containing ABC, BCD and BC: R >= H_2.
Mixed: C0=ABCD (a=0, so t_0 >= 0 >= b_0... precisely t_0 = -b_0 >= 0), C1=ABC (a_1 = t_0 - kappa_1),
C2=BCD (a_2 = b_0 + rho_2), C3=BC (side draw: a_3 = b_1 + t_2). Light exits: top of C1, top of C3.
t_1 >= min(0,a_1) >= -kappa_1. If a_3 >= 0: t_3 >= 0, done. Else t_3 >= a_3 and t_1 + t_3 >= a_1 + t_2;
if a_2 >= 0, t_2 >= 0 and the sum is >= -kappa_1; else t_2 >= a_2 and the sum is
>= t_0 - kappa_1 + b_0 + rho_2 = rho_2 - kappa_1 >= -kappa_1. Insert in (*). QED

## Lemma 4 (recovery ordering in the root column of K). If [1,4] splits as ABC/BCD (TG: theta_2 interior,
theta_1, theta_3 boundary), then its top recoveries satisfy d_B/z_B > d_C/z_C.
Proof. Phi_1 := g_T(theta_2) - g_T(theta_1) >= 0 and Phi_3 := g_T(theta_2) - g_T(theta_3) >= 0 (TG).
On the diagonal d_B = r z_B, d_C = r z_C, using sum_p a_p z_p/(a_p - theta_r) = 1-q:
g_T(theta_r) = r(1-q) + (1-r) A(theta_r) + r D(theta_r), A(x) = a_A z_A/(a_A - x) (increasing),
D(x) = a_D z_D/(x - a_D) (decreasing). So Phi_1 = -(1-r)A12 + r D12, Phi_3 = (1-r)A32 - r D32 with
A12 = A(theta_1)-A(theta_2), D12 = D(theta_2)-D(theta_1), A32 = A(theta_2)-A(theta_3),
D32 = D(theta_3)-D(theta_2), all > 0, and
A12 D32 / (A32 D12) = (a_A - theta_3)(theta_1 - a_D) / ((a_A - theta_1)(theta_3 - a_D)) > 1.
Hence Phi_1 >= 0 and Phi_3 >= 0 cannot both hold on the diagonal (r = 0, r = 1 and 0<r<1 separately).
dPhi_1/dd_B = a_B/(a_B - theta_2) + a_B/(theta_1 - a_B) > 0 and dPhi_3/dd_B = a_B/(a_B-theta_2) - a_B/(a_B-theta_3)
> 0: lowering d_B below r z_C... lowers both, so no feasible point has d_B/z_B <= d_C/z_C. QED
(Consequence: d_B > 0 and b_C = z_C - d_C > 0.)

## Lemma 5 (family K, all 16 TC/HX variants). R >= H_2.
C0 = ABCD -> ABC/BCD: interior root theta_2 gives t_0 = Lam_{0T}(theta) = V_{0T} - Y_{0T} =: e >= 0 and
b_0 = -e. C1 = ABC -> AB/C fed with d = (z_A, d_B, d_C): Lam_{1F} = Lam_{0T} - kappa_1; TG root rho_1 in I_2,
Lam_{1T}(rho_1) >= 0. C2 = BCD -> B/CD fed with b = (b_B, b_C, z_D): Lam_{2F} = -Lam_{0B} + rho_2; root rho_2,
Lam_{2T}(rho_2) >= 0. With Delta_p(x) := 1/(a_p - x) - 1/(a_p - theta) (sign of x - theta for every p):
  t_1 >= sum_{A,B} a_p d_p Delta_p(rho_1),   e - kappa_1 = sum_{A,B,C} a_p d_p Delta_p(rho_1),
  t_2 >= a_B b_B Delta_B(rho_2),             -e + rho_2 = sum_{B,C,D} a_p b_p Delta_p(rho_2).
Need t_1 + t_2 + kappa_1 >= 0 (then (*)). Cases: rho_1, rho_2 >= theta: both t >= 0. rho_1, rho_2 <= theta:
t_1 + kappa_1 >= e and t_2 >= -e. rho_1 <= theta <= rho_2: t_2 >= 0, t_1 + kappa_1 >= e.
rho_1 >= theta >= rho_2: t_1 + kappa_1 + t_2 >= e - X - Z with X = a_C d_C Delta_C(rho_1), Z = a_B b_B |Delta_B(rho_2)|,
Y = a_B d_B Delta_B(rho_1), W = a_C b_C |Delta_C(rho_2)|, X + Y <= e, Z + W <= e. If XZ <= YW then
X + Z <= e (X/(X+Y) + Z/(Z+W) <= 1). XZ <= YW <=> (d_C b_B)/(d_B b_C) <= Rat(rho_2)/Rat(rho_1) with
Rat(x) = Delta_C(x)/Delta_B(x) = (a_B - x)(a_B - theta)/((x - a_C)(theta - a_C)), decreasing on I_2; the left side
is <= 1 by Lemma 4 (d_B z_C >= d_C z_B <=> d_B b_C >= d_C b_B) and the right side is >= 1. QED

## Theorem B (N = 4). For every one of the 152 configurations, every TG-feasible operating point, every
a_1>a_2>a_3>a_4>0, z>0, q real: R >= max(H_1, H_2, H_3) = V_FTC.
Proof. j = 1, 3: Theorem A. j = 2: the mixed streams are among ABCD, ABC, BCD, BC; ABC (resp. BCD) is present
iff it is the distillate (resp. residue) of ABCD; BC, if present, is the residue of ABC and the distillate of
BCD whenever those are present. Path cases (Lemma 2), diamond (Lemma 3), family K (Lemma 5). QED
Verified independently by the sign-relaxation LP (theory/scripts/signcheck.py): for N = 4 the only
configurations not settled by Lemma 1 alone are the 16 variants of K, j = 2.

## Achievability (FTC attains V_FTC) — status
HS Part 3 construction: every FTC column at its preferred split (all carried-over roots active, e = 0),
except the binary bottom-row columns [i,i+1], which carry excess e_i := max_r P_r^T - P_i^T >= 0; the liquid
side-draw balances V_T[i+1,i+2 col] = V_B[i,i+1 col] reduce to e_{i+1} - e_i = P_i^T - P_{i+1}^T (cut values),
and R = H_{N-1} + e_{N-1} = V_FTC. Nonnegativity of all preferred-split flows: to be checked (HS's Vmin
diagram) — Agent C verifies numerically that the FTC reaches V_FTC in the model.
