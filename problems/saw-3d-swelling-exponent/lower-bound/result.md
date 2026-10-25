# Lower-bound solver: result

Solver: Claude Opus 5.5 (effort xhigh). Date: 2026-10-24.

## 0. Verdict

**No new unconditional lower bound on E_N|omega_N|^2 was obtained.** The best rigorous lower bound for strict SAW on
Z^3 that I know of is still Madras's (CITED below, L0): E_N||omega(N)||_2^2 >= N^{4/9}/36 for every N >= 1. Every rung
L1-L4 is still open, and so are the weakened versions L1' (subsequence, Cesaro, generating function).

This file contains the following.

1. The frontier: Madras 2014, stated exactly (CITED), with the parts of its proof that I reuse re-proved here (Madras
   omits the proof of his Lemma 2.1).
2. **Theorem 3 (PROVED, conditional).** A reorganisation of Madras's argument into one inequality with two inputs:
   (H1) a lower bound R on the span, and (H2) a bound C_0 m^{1-alpha} on the *average* multiplicity of Madras's
   reflection map. The output is E_N|omega_N|^2 >= (R^2/(4 C_0))^{2/(3-alpha)}.
   - With the inputs that are known unconditionally (R = N^{1/d}/3, alpha = 0) it gives back Madras's exponent.
   - Any power saving alpha > 0 in the average multiplicity improves the exponent to 4/(d(3-alpha)).
   - alpha = 1 gives the open bound E >= c N^{2/d}.
   - Heuristically, with the true inputs (R ~ N^nu, multiplicity m^{o(1)}) the inequality gives back the correct
     exponent 2 nu. The method therefore loses nothing beyond these two inputs.
3. **Proposition 2 (PROVED).** The exact structure of the multiplicity. Two near walks with the same Madras image
   differ by a "fold" at a renewal depth j <= m of the head or of the tail, measured from the walk's top. After the
   fold the end piece of the walk is translated by 2j e_I, and it must still avoid the rest of the walk. This turns
   (H2) into a precise statement about renewal depths near the top and about head-tail collisions.
4. Heuristic Monte Carlo evidence (not a proof) that (H2) holds with a large alpha. For N = 4000 the average
   multiplicity is 2.7, 3.1 and 3.9 for m = 16, 32 and 64, against Madras's worst case m + 1 = 17, 33 and 65.

Nothing in Sections 2 and 3 is claimed to be deep. Theorem 3 and Proposition 2 are elementary consequences of
Madras's construction. Their only purpose is to isolate exactly what is missing. I found neither in print, but they are
the kind of observation an expert would make in passing, so I claim no novelty for them as results.

## 1. Setting and the cited frontier

d >= 2, Z^d nearest-neighbour. S_N = N-step self-avoiding walks omega = (omega(0) = 0, ..., omega(N)) (strict
self-avoidance). P_N and E_N denote the uniform measure on S_N and its expectation. ||x||_inf = max_i |x_i| and
Lambda_r = {x : ||x||_inf <= r}.

**CITED (L0).** Madras, *A lower bound for the end-to-end distance of the self-avoiding walk*, Canad. Math. Bull. 57
(2014) 113-118, doi:10.4153/CMB-2012-022-6, Proposition 1.1 (p. 114). Read in full: the Cambridge open-access PDF and
the publisher PDF. "Let N be a positive integer and let p be a real number in [1, inf). Then
(E_N ||omega(N)||_p^p)^{1/p} >= (E_N ||omega(N)||_inf^p)^{1/p} >= (1/6) N^{p/((p+1)d)}."
With p = 2: E_N||omega(N)||_2^2 >= N^{4/(3d)}/36 for every N >= 1, which is N^{4/9}/36 in d = 3.

Open statements, verbatim:
- Madras 2014, p. 114: "Even more simply, one could reasonably expect a lower bound with nu_L = 1/d, since the
  average distance of all points of an N-step self-avoiding walk from omega(0) must be at least of order N^{1/d}.
  Alas, no rigorous proof is known even for this bound."
- Madras 2014, p. 118 (Acknowledgments): "I was hoping that the argument could be improved, for example, by showing
  that the restriction of T to S_{N,<=m} is typically much better than (m+1)-to-one. However, no such improvement was
  apparent."
- Bauerschmidt, Duminil-Copin, Goodman, Slade, *Lectures on self-avoiding walks*, arXiv:1206.2092, p. 6, around
  (1.30): "it has not been proved that cn <= E_n|omega(n)|^2 <= Cn^{2-eps}, or even that the endpoint is typically as
  far away as the surface of a ball of volume n, i.e., cn^{2/d} <= E_n|omega(n)|^2."
- Duminil-Copin, Hammond, *Self-avoiding walk is sub-ballistic*, arXiv:1205.0401 (CMP 324 (2013)), Section 1.5,
  Question 4: "Show that lim_{n->inf} n^{-2/d} <||gamma_n||^2> = inf."
- Madras-Slade, *The Self-Avoiding Walk* (2013 reprint), p. 14: "no general proof is known of the mean-field bound
  nu >= 1/2"; p. 18: "it has not yet been proved rigorously that G_{z_c}(0,x) is even finite for d = 2, 3 or 4, for any
  value of x != 0." The statements describe the state of knowledge in 1993/2013.
  The literature agent's sweep (literature/quick-facts.md) found nothing newer that changes the lower side.

### Madras's map (definitions from Madras 2014, Section 2)

For omega in S_N define the following.
- span(omega) = max_{0<=k<=N} ||omega(k)||_inf (Madras writes 9omega9).
- L = min{k > 0 : ||omega(k)||_inf = span(omega)}.
- I = min{i : |omega_i(L)| = span(omega)}.
- H = omega_I(L). Then |H| = span(omega).
- U = U_H is the reflection through the hyperplane {x_I = H}: (U x)_I = 2H - x_I and (U x)_i = x_i for i != I.
- T omega = psi, with psi(k) = omega(k) for k <= L and psi(k) = U omega(k) for k > L.

**Lemma 1 (Madras 2014, Lemmas 2.1 and 2.2; re-proved here).** Let omega in S_N, rho = span(omega), and let m be an
integer with ||omega(N)||_inf <= m <= rho. Then:
- (0) psi = T omega is in S_N.
- (i) psi_i(N) = omega_i(N) in [-m, m] for i != I.
- (ii) psi_I(N) has the sign of H.
- (iii) 2|H| - m <= |psi_I(N)| <= 2|H| + m, and |H| <= 2|H| - m.
- (iv) ||psi(N)||_inf = |psi_I(N)|.
- (v) if H > 0, then psi_I(k) < H for k < L, psi_I(L) = H and psi_I(k) >= H for k >= L. The case H < 0 is the same
  with the inequalities reversed.
- (vi) [Madras Lemma 2.2] If N > 3^d, then span(omega) > N^{1/d}/3.

*Proof.* Take H > 0; the case H < 0 is symmetric.
- (v): for k < L we have ||omega(k)||_inf < rho = H. This includes k = 0, since rho >= 1. Hence omega_I(k) < H. For
  every k, omega_I(k) <= rho = H, so for k >= L, psi_I(k) = 2H - omega_I(k) >= H.
- (0): psi(L) = omega(L) is fixed by U, and U is an isometry of Z^d, so psi is a nearest-neighbour path. Points with
  indices k < L <= l are distinct by (v). Points with indices both < L, or both >= L, are distinct because omega and
  U are injective.
- (i): immediate from the definition of T.
- (ii) and (iii): psi_I(N) = 2H - omega_I(N) and |omega_I(N)| <= m <= H.
- (iv): |psi_I(N)| >= 2H - m >= H >= m >= |psi_i(N)| for i != I.
- (vi): if span(omega) <= N^{1/d}/3, the N+1 distinct sites of omega lie in a cube with at most
  (2 N^{1/d}/3 + 1)^d < N sites when N > 3^d. This is Madras's computation, p. 116. []

## 2. Multiplicity: definition, average version, exact sibling structure

For real R > 0 and an integer 1 <= m < R, define A_{m,R} = {omega in S_N : ||omega(N)||_inf <= m, span(omega) >= R}.
For psi in S_N, define K(psi) = |{theta in A_{m,R} : T theta = psi}|. Define the **average multiplicity**

    kappa_N(m,R) = (1/|A_{m,R}|) * sum_{omega in A_{m,R}} K(T omega)   (when A_{m,R} is non-empty).

**Proposition 1 (PROVED; for R = N^{1/d}/3 this is Madras's Lemma 2.3 with average multiplicity in place of the worst
case).** Let 1 <= m < R.
- (a) If omega in A_{m,R}, then ||T omega(N)||_inf >= 2R - m > R.
- (b) If theta, phi in A_{m,R} with T theta = T phi, then I(theta) = I(phi) and sign H(theta) = sign H(phi). Moreover,
  H(theta) = H(phi) implies theta = phi. Consequently K(psi) <= m + 1 for all psi.
- (c) |T(A_{m,R})| >= |A_{m,R}| / kappa_N(m,R), so P_N(||omega(N)||_inf >= R) >= P_N(A_{m,R}) / kappa_N(m,R).

*Proof.*
- (a): Lemma 1(iii),(iv) with |H| = span >= R.
- (b): Let psi = T theta = T phi. Suppose J := I(phi) != I(theta). Lemma 1(i) for theta gives |psi_J(N)| <= m. Lemma
  1(iii),(iv) for phi give |psi_J(N)| >= 2R - m > m, a contradiction. The signs agree by Lemma 1(ii). If
  H(theta) = H(phi) = H, then by Lemma 1(v), L(theta) = L(phi) = min{k : psi_I(k) = H}. Hence
  theta = psi[0,L] followed by U_H(psi[L,N]) = phi. By Lemma 1(iii), |H| lies in an interval of length m. With the
  sign fixed, H takes at most m + 1 integer values.
- (c): |A| = sum_{psi in T(A)} K(psi) and sum_{omega in A} K(T omega) = sum_{psi in T(A)} K(psi)^2. By Cauchy-Schwarz,
  (sum K)^2 <= |T(A)| * sum K^2, hence |T(A)| >= |A|^2 / sum_{omega} K(T omega) = |A|/kappa. Finally
  T(A) is contained in {||psi(N)||_inf >= R} by (a). Divide by |S_N|. []

**Definitions (renewal depths at the top).** Let omega in S_N with Madras data (L, I, H). Suppose H > 0; for H < 0
replace x_I by -x_I throughout. For an integer j >= 1:
- the **tail has a renewal at depth j** if tau_j := min{k >= L : omega_I(k) = H - j} exists and omega_I(k) <= H - j
  for all k >= tau_j;
- the **head has a renewal at depth j** if, with sigma_j := min{k : omega_I(k) = H - j} (which exists when
  H - j >= 0), omega_I(k) >= H - j for all sigma_j <= k <= L.

**Proposition 2 (exact sibling structure; PROVED).** Let 1 <= m < R, omega in A_{m,R} with H > 0, psi = T omega,
and theta in A_{m,R} \ {omega} with T theta = psi. Then I(theta) = I and H(theta) = H + j for an integer j with
1 <= |j| <= m. Moreover:
- (a) If j > 0, the tail of omega has a renewal at depth j, and
  theta = omega[0, L] ∘ U_H(omega[L, tau_j]) ∘ (omega[tau_j, N] + 2j e_I).
  In particular (omega[tau_j, N] + 2j e_I) ∩ omega[0, L] is empty, and theta(N) = omega(N) + 2j e_I, so
  |omega_I(N) + 2j| <= m.
- (b) If j = -i < 0, the head of omega has a renewal at depth i, and
  theta = omega[0, sigma_i] ∘ U_{H-i}(omega[sigma_i, L]) ∘ (omega[L, N] - 2i e_I).
  In particular (omega[L, N] - 2i e_I) ∩ omega[0, sigma_i] is empty, and theta(N) = omega(N) - 2i e_I, so
  |omega_I(N) - 2i| <= m.

(Here a∘b denotes the concatenation of two pieces that share their junction point.) Different theta correspond to
different j. Consequently

    K(T omega) <= 1 + #{j in [1,m] : the tail renews at depth j, (omega[tau_j,N] + 2j e_I) ∩ omega[0,L] is empty}
                    + #{i in [1,m] : the head renews at depth i, (omega[L,N] - 2i e_I) ∩ omega[0,sigma_i] is empty}.

*Proof.* By Proposition 1(b), I(theta) = I, H(theta) > 0 and H(theta) != H. Lemma 1(iii) applied to theta and to
omega puts both H and H(theta) in the same window, so |j| <= m. Write H' = H(theta) and L' = L(theta). By Lemma 1(v)
applied to theta, L' = min{k : psi_I(k) = H'}, and theta = psi[0, L'] ∘ U_{H'}(psi[L', N]). Since theta is in
Lambda_{H'} (its span is |H'|), theta_I(k) <= H' for all k.

(a) j > 0. By Lemma 1(v) for omega, psi_I(k) < H < H' for k < L, and psi_I(k) = 2H - omega_I(k) for k >= L.
Hence L' = min{k >= L : omega_I(k) = H - j} = tau_j, and tau_j exists. For k >= tau_j,
theta_I(k) = 2H' - psi_I(k) = omega_I(k) + 2j, so theta_I <= H' is equivalent to omega_I(k) <= H - j: the tail
renews at depth j. For k >= tau_j, theta(k) = U_{H'} U_H omega(k) = omega(k) + 2(H' - H) e_I, because the product of
reflections in two parallel hyperplanes is a translation. For L <= k <= tau_j, theta(k) = psi(k) = U_H omega(k). For
k <= L, theta(k) = omega(k). Self-avoidance of theta gives the disjointness. The endpoint identity is the case k = N,
and ||theta(N)||_inf <= m gives the last inequality.

(b) j = -i < 0, so H' = H - i. For k <= L, psi_I(k) = omega_I(k), with omega_I(0) = 0, omega_I(L) = H and unit
steps. Since 0 <= H' < H (H' > 0 by Proposition 1(b)), the level H' is first reached at sigma_i < L, so L' = sigma_i.
For k >= sigma_i, theta_I(k) = 2H' - psi_I(k) <= H' iff psi_I(k) >= H'. For sigma_i <= k <= L this is
omega_I(k) >= H - i: the head renews at depth i. For k > L, psi_I(k) >= H > H' holds automatically. For k >= L,
theta(k) = U_{H'} U_H omega(k) = omega(k) - 2i e_I. The rest is as in (a).

The final bound: theta is determined by H(theta) (Proposition 1(b)), hence by j. Each j or i that occurs satisfies
the stated necessary conditions. []

## 3. Conditional lower bound

**Theorem 3 (PROVED, conditional).** Let d >= 2 and N >= 1. Let R >= 2 be real, alpha in [0,1] and C_0 >= 1, and
assume:
- (H1) P_N(span(omega) < R) <= 1/4;
- (H2) kappa_N(m, R) <= C_0 m^{1-alpha} for every integer m with 1 <= m < R and A_{m,R} non-empty.

Then

    E_N ||omega(N)||_2^2  >=  E_N ||omega(N)||_inf^2  >=  ( R^2 / (4 C_0) )^{2/(3-alpha)}.

*Proof.* Let B = (E_N||omega(N)||_inf^2)^{1/2}. Then B >= 1, because ||omega(N)||_inf >= 1 for N >= 1. By Markov's
inequality, P_N(||omega(N)||_inf >= 2B) <= 1/4. Let m = ceil(2B) - 1. This is an integer with 1 <= m < 2B, and since
||omega(N)||_inf is an integer, {||omega(N)||_inf < 2B} = {||omega(N)||_inf <= m}. Hence
P_N(||omega(N)||_inf <= m) >= 3/4.
- Case m >= R. Then B > m/2 >= R/2, so B^2 >= R^2/4 >= (R^2/(4C_0))^{2/(3-alpha)}. The last inequality holds
  because R^2/(4C_0) <= R^2/4, because R^2/4 >= 1, and because 2/(3-alpha) <= 1.
- Case m < R. By (H1), P_N(A_{m,R}) >= 3/4 - 1/4 = 1/2. Proposition 1(c) and (H2) give
  P_N(||omega(N)||_inf >= R) >= 1/(2 C_0 m^{1-alpha}) >= 1/(2 C_0 (2B)^{1-alpha}). Markov's inequality gives
  B^2 >= R^2 P_N(||omega(N)||_inf >= R) >= R^2 / (2^{2-alpha} C_0 B^{1-alpha}). Hence
  B^{3-alpha} >= R^2/(2^{2-alpha} C_0) >= R^2/(4 C_0). []

**Corollary 4 (Madras's exponent recovered; PROVED, not new).** Take N >= 6^d and R = N^{1/d}/3 >= 2.
- (H1) holds with probability 0 by Lemma 1(vi).
- (H2) holds with alpha = 0 and C_0 = 2, since kappa <= m + 1 <= 2m.

So E_N||omega(N)||_2^2 >= (N^{2/d}/72)^{2/3} = 72^{-2/3} N^{4/(3d)}. In d = 3 and for N >= 216 this is
E_N|omega_N|^2 >= N^{4/9}/17.4. Same exponent as Madras; only the constant differs, which is of no interest.

**What the theorem isolates.** With R = N^{1/d}/3 (the only span bound known), (H2) with a fixed alpha > 0 would give
the new exponent 4/(d(3-alpha)), which is > 4/(3d). In d = 3 that is 4/(3(3-alpha)), e.g. 0.533 for alpha = 1/2 and
2/3 for alpha = 1. alpha = 1 (bounded average multiplicity) gives c N^{2/d}, the bound stated open by Madras, BDGS and
DCH (quotes in Section 1). Beyond 2/d the theorem needs (H1) with R = N^{1/d+eta}. Then the exponent is
4(1/d+eta)/(3-alpha).

HEURISTIC remark (not used anywhere): if R ~ N^nu and kappa(m) = m^{o(1)} (alpha -> 1), the bound reads
N^{2nu - o(1)}, which is the true order. So Madras's construction is not intrinsically lossy. The whole loss lies in
the two inputs: (H1) beyond the volume scale, and (H2) beyond the worst case.

**What (H2) requires (from Proposition 2).** An alternative preimage at relative level j requires all three of:
- a renewal of the head or the tail at depth |j| <= m below the top in the direction I;
- that the end piece of the walk, translated by 2j e_I, still avoids the head;
- that the translated endpoint stays in Lambda_m.

So (H2) is a statement about (1) the density of renewal depths in the first m levels below the walk's maximum, and
(2) head-tail interpenetration at scale m. Neither is controlled by any result I found.
- For bridges under the critical renewal measure P_iSAB, zero density of renewal *levels* is equivalent to
  B_T(z_c) -> 0. This is not known on Z^d: see literature/quick-facts.md (proved only on the hexagonal lattice,
  Krachun-Panagiotis). Duminil-Copin-Hammond Thm 2.5 (E_iSAB|gamma| = inf) gives zero density of renewal *times*,
  not levels.
- No estimate on head-tail collisions is available.

## 4. Heuristic evidence on (H2) (Monte Carlo; NOT a proof, floating point)

Script: `scripts/madras_multiplicity_mc.cpp`. It samples uniform N-step SAWs on Z^3 with the pivot algorithm, 10
pivot attempts between samples; equilibration is 40N attempts from a straight rod, from a single seed per run. For
every sample with ||omega(N)||_inf <= m it computes K(T omega) exactly: every candidate level is decoded and checked
for self-avoidance, Madras data and endpoint.
- The near set used is S_{N,<=m}, with no span restriction.
- Measured mean span of near walks: 39-43 at N = 1000 and 88-95 at N = 4000. The fraction of near walks with
  span <= m is 0 for m <= 16, 11% for (N, m) = (1000, 32), and 3.4% for (4000, 64).
- Rows where this fraction is large (N = 200, and N = 1000 with m = 64) are left out of the table, because there the
  multiplicity is not Madras's.
- Candidate preimages with a different reflection coordinate are not searched. They can exist only when span <= m.

Raw output: `data/mc_multiplicity_N{200,1000,4000}.txt`.

| N | m | near samples | mean K (= kappa estimate) | m+1 | mean # levels passing renewal test | mean # passing also self-avoidance |
|---|---|---|---|---|---|---|
| 1000 | 8 | 166 | 1.95 | 9 | 5.6 | 2.0 |
| 1000 | 16 | 1319 | 2.59 | 17 | 8.8 | 2.9 |
| 1000 | 32 | 10114 | 3.19 | 33 | 14.1 | 4.2 |
| 4000 | 16 | 87 | 2.72 | 17 | 8.0 | 2.8 |
| 4000 | 32 | 681 | 3.07 | 33 | 12.3 | 3.4 |
| 4000 | 64 | 5733 | 3.88 | 65 | 19.9 | 4.8 |

Reading (HEURISTIC):
- kappa grows like a small power of m, roughly m^{0.25-0.35}, or logarithmically. It looks almost independent of N
  at fixed m. If this persists, (H2) would hold with alpha around 0.65-0.75, giving an exponent around 0.57-0.59 in
  d = 3. Not proved.
- Both mechanisms of Proposition 2 matter. Renewal depths near the top are frequent: about 40% of depths for m = 16,
  about 25% for m = 64. Head-tail collisions after folding remove roughly three quarters of the renewal-admissible
  levels.
- Caveats: autocorrelation was not estimated; samples at small m are few; a single seed per run.

## 5. Status summary

| Item | Status | Where |
|---|---|---|
| E_N||omega(N)||_2^2 >= N^{4/(3d)}/36, all N (L0) | CITED (Madras 2014, Prop. 1.1) | Section 1 |
| Lemma 1 (Madras Lemmas 2.1, 2.2) | PROVED here (Madras omits the proof of 2.1) | Section 1 |
| Proposition 1 (average-multiplicity form of Madras Lemma 2.3) | PROVED (elementary) | Section 2 |
| Proposition 2 (exact sibling structure) | PROVED (elementary) | Section 2 |
| Theorem 3 (conditional bound from (H1), (H2)) | PROVED (conditional) | Section 3 |
| (H2) with some alpha > 0 | CONJECTURED (MC evidence only) | Section 4 |
| L1, L1', L2, L3, L4 | FAILED (not obtained) | see attempts.md |
