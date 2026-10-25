# Upper-bound solver: result

Solver: Claude Opus 5.5 (effort xhigh). Date: 2026-10-24. Status: **no unconditional progress on U1/U2.**
What is proved below is a *conditional* quantitative sub-ballisticity theorem for strict SAW on Z^3, an
unconditional inequality behind it, and an elementary lemma that reduces the two critical inputs of the
Krachun–Panagiotis method to a single one in dimension three. None of this is U1 or U2; see Sections 6-7.

## 0. Summary of claims

| # | Statement | Status |
|---|---|---|
| L1 | Arch lemma: A(z) <= sum_{M>=1} B_M(z)^2 for every z >= 0 | PROVED (Section 2) |
| P2 | Walk-to-bridge reduction with explicit loss (Hammersley–Welsh, BDGS form, partitions bounded elementarily) | PROVED from CITED BDGS steps (Section 3) |
| P3 | Bridge decomposition bound (Krachun–Panagiotis Lemma 4.1 transplanted to Z^3, with L1 replacing their Lemma 2.2) | PROVED (Section 4) |
| TA | **Theorem A** (unconditional): explicit bound on P_n(max_k ||gamma_k||_inf >= h) in terms of B_m(z_c) and sum_M e^{-2tM} B_M(z_c)^2 | PROVED (Section 4) |
| CB | **Corollary B** (conditional U1): if S := sum_{T>=1} B_T(z_c)^2 < infinity then P_n(max_k ||gamma_k||_inf >= 24 log(2K) n/log n) <= exp(-(log(2K)/4) n^{2/3}) for n >= n_0(S), K = max(1, 4 z_c + 2S); hence E_n||gamma_n||^2 <= (1+o(1)) 1728 (log 2K)^2 n^2/(log n)^2 | PROVED, conditional on hypothesis (H) |
| CB2 | Same with the weaker hypothesis sum_{M<=L} B_M(z_c)^2 <= C_2 (1+log L)^p: rate n log log n / log n | PROVED, conditional (Section 5) |
| P5 | GF version: m(z_c e^{-t}) >= (1/m) log(1/B_m(z_c)) whenever z_c e^{-t} theta^{1/m} < z_c, theta >= 1 + 4z' + 2A(z') | PROVED (Section 5) |
| CC | **Corollary C** (conditional): under (H), the mass satisfies m(z_c e^{-t}) >= t log(1/t)/(9 log(1 + 4z_c + 2S)) for small t (correlation length o((z_c - z)^{-1}) with a log rate) | PROVED, conditional on (H) |
| P4 | B_T(z_c) is the renewal sequence of the heights of irreducible bridges; B_T(z_c) -> 1/E_iSAB[height]; (H) implies E_iSAB[height] = infinity | PROVED from CITED (Section 6) |
| H | (H) holds in d = 3 (B_T(z_c) ~ T^{-a}, a = 1 - gamma_b/nu ~ 0.66 > 1/2) and fails in d = 2 (a = 1/4) | HEURISTIC (scaling + numerics, Section 6) |

Everything is written for d = 3; every proof works verbatim on Z^d, d >= 2, with 4 = 2(d-1) and 6 = 2d.

## 1. Notation and cited facts

Z^3 with coordinates x = (x_1, x_2, x_3); the "height" of a point is x_1. ||x||_inf = max_i |x_i|,
||x|| = ||x||_2 <= sqrt(3) ||x||_inf. SAW_n = n-step self-avoiding walks gamma = (gamma(0) = 0, ..., gamma(n)),
c_n = |SAW_n|, P_n, E_n uniform. mu = lim c_n^{1/n}, z_c = 1/mu. For a walk omega we write omega_1(i) for the
height of omega(i), |omega| for its number of steps.

* **Bridge** (BDGS Definition 2.1, the same as DCH Section 2.2 up to the choice of coordinate): an n-step SAW with
  omega_1(0) < omega_1(i) <= omega_1(n) for i = 1, ..., n; the 0-step walk is a bridge. Its height is
  H(omega) = omega_1(n) - omega_1(0) (= its span). b_{n,T} = number of n-step bridges from 0 of height T;
  b_n = sum_T b_{n,T}; b^{>=h}_n = sum_{T>=h} b_{n,T}. (b_{n,T} = 0 unless T <= n, and b_{0,0} = 1.)
* **Critical bridge partition function at height T**: B_T(z) = sum_{n>=0} b_{n,T} z^n in [0, infinity]; B_0 = 1.
* **Half-space walk** (BDGS Definition 2.2): omega_1(0) < omega_1(i) for i = 1..n; h_n their number (h_0 = 1).
* **Upper arch**: an l-step SAW, l >= 2, with omega_1(l) = omega_1(0) and omega_1(i) > omega_1(0) for 0 < i < l.
  **Lower arch**: same with "<". A(z) = sum over upper arches from 0 of z^{|omega|}; by the reflection
  x_1 -> -x_1 the lower arches from 0 have the same generating function.
* P_D(A) = number of ways to write A = A_1 + ... + A_k with A_1 > ... > A_k >= 1 (BDGS Theorem 2.4 notation);
  P*_D(N) = max(1, max_{1<=A<=N} P_D(A)) for N >= 0 (so P*_D(0) = 1).
* **Hypothesis (H)**: S := sum_{T>=1} B_T(z_c)^2 < infinity.

Cited facts (each read in the source):

* (F1) c_n >= mu^n for all n >= 0. [BDGS lectures, arXiv:1206.2092, eq. (1.12) (from the submultiplicativity
  (1.10) and Lemma 1.1).] Also d <= mu <= 2d - 1, so 1/5 <= z_c <= 1/3 in d = 3 [BDGS (1.13)].
* (F2) Splitting at the last minimum of the height: for omega in SAW_n, let x_1 = min_i omega_1(i) and
  m = max{i : omega_1(i) = x_1}. Then (omega(m), ..., omega(n)) translated is an (n-m)-step half-space walk and
  (omega(m) - e_1, omega(m), omega(m-1), ..., omega(0)) translated is an (m+1)-step half-space walk, and
  omega -> (m, the two half-space walks) is injective. [BDGS, proof of Theorem 2.6, eq. (2.14) and (2.15).]
* (F3) For a half-space walk, define n_0 = 0 and the alternating last-extremum times n_1 < n_2 < ... < n_K = n
  and the decreasing values A_1 > A_2 > ... > A_K >= 1 as in BDGS (2.8)-(2.9); A_1 = max_j omega_1(j) - omega_1(0)
  is the span. If h_n[a_1, ..., a_k] counts n-step half-space walks with K = k and A_i = a_i, then
  h_n[a_1, ..., a_k] <= b_{n, a_1 + ... + a_k}. [BDGS, proof of Proposition 2.5, eqs. (2.10)-(2.11): reflections
  of the part after n_1, injective for fixed (a_1, ..., a_k).]
* (F4) The radius of convergence of sum_n c_n z^n is z_c [BDGS (1.12): c_n^{1/n} -> mu].

## 2. Elementary lemmas

**Lemma 0.** (a) For a, c >= 0 and h >= 1: b^{>=h}_a b_c <= b^{>=h}_{a+c} and b_a b^{>=h}_c <= b^{>=h}_{a+c};
also B_{T+T'}(z) >= B_T(z) B_{T'}(z).
(b) B_T(z_c) <= 1 for every T >= 0.
(c) For t >= 0 and z = z_c e^{-t}: B_T(z) <= e^{-tT} B_T(z_c).
(d) P*_D(N) <= (N+1)^{sqrt(2N)+1} for N >= 1.

*Proof.* (a) Concatenation: if beta (a steps, height H) and beta' (c steps, height H') are bridges, the walk
beta o beta' (beta' translated to start at beta(a)) has its beta-part at heights <= beta_1(a) and its beta'-part,
except the junction, at heights > beta_1(a); so it is self-avoiding, its start is the unique minimum, its end is a
maximum: a bridge of height H + H' with a + c steps. For fixed a the map (beta, beta') -> beta o beta' is injective,
and H + H' >= h if H >= h or H' >= h. For B_{T+T'} >= B_T B_{T'}: the junction is the last index of height
<= T, so the map is injective on (height T) x (height T') bridges of all lengths, and lengths add.
(b) This is Hutchcroft, ECP 2018 (arXiv:1708.09460v3), Lemma 2.3 with (2.2), and also Madras–Slade (2013
reprint) (4.1.12) with Theorem 4.1.13 (B_z(L) <= e^{-L Mbar(z)}, Mbar(z_c) = 0); we repeat the argument. For
0 <= z < z_c, sum_T B_T(z) = sum_n b_n z^n <= sum_n c_n z^n < infinity by (F4). If B_{T0}(z) > 1 for some T0,
then by (a) B_{k T0}(z) >= B_{T0}(z)^k -> infinity, contradicting summability. So B_T(z) <= 1 for z < z_c, and
B_T(z_c) = lim_{z up z_c} B_T(z) <= 1 by monotone convergence.
(c) A bridge of height T has at least T steps, and z^n = z_c^n e^{-tn} <= z_c^n e^{-tT} for n >= T.
(d) A partition of A into k distinct parts has A >= 1 + 2 + ... + k, so k <= sqrt(2A); it is a k-subset of
{1, ..., A}. Hence P_D(A) <= sum_{k=1}^{floor(sqrt(2A))} binom(A, k) <= sqrt(2A) A^{sqrt(2A)}
<= (A+1)^{sqrt(2A)+1}, which is nondecreasing in A. QED

**Lemma 1 (arches are bounded by pairs of bridges).** For every z >= 0, A(z) <= sum_{M>=1} B_M(z)^2.

*Proof.* Let alpha be an upper arch from 0 with l steps. Let M = max_i alpha_1(i) >= 1 and
k = min{i : alpha_1(i) = M}; then 1 <= k <= l - 1. Put beta = (alpha(0), ..., alpha(k)) and
rho = (alpha(l) , alpha(l-1), ..., alpha(k)) - alpha(l) (the reversed second part, translated to start at 0).
* beta is a bridge of height M: alpha_1(0) = 0 < alpha_1(i) for 0 < i <= k (interior points of an arch are
  strictly above 0) and alpha_1(i) <= M = alpha_1(k).
* rho is a bridge of height M: rho_1(0) = 0 and rho_1(j) = alpha_1(l-j) > 0 for 1 <= j <= l-k (these are
  interior points of the arch, or alpha(k) whose height is M >= 1), and rho_1(l-k) = M is the maximum.
* The map alpha -> (beta, rho) is injective: alpha is beta followed by the reversal of rho translated so that
  its endpoint is beta(k); and |alpha| = |beta| + |rho|.
Summing z^{|alpha|} = z^{|beta|} z^{|rho|} over arches with maximal height M gives at most B_M(z)^2; sum over M.
QED

(Remark: dropping the mutual avoidance of beta and rho is the only loss. In d = 2 the right-hand side is expected
to be infinite at z_c while A(z_c) is finite, so the lemma is useless there; see Section 6.)

## 3. From walks to bridges (Hammersley–Welsh with an explicit loss)

**Proposition 2.** For n >= 1 and h >= 1,
  #{gamma in SAW_n : max_k gamma_1(k) >= h} <= 2 (n+1) P*_D(n+1)^2 b^{>=h}_{n+1}.
Consequently, for every h >= 1,
  P_n( max_{0<=k<=n} ||gamma(k)||_inf >= h ) <= 12 mu (n+1) P*_D(n+1)^2 z_c^{n+1} b^{>=h}_{n+1}.     (3.1)

*Proof.* Let h^{>=h}_j be the number of j-step half-space walks with span >= h (h^{>=h}_0 = 0). By (F3) the span
of a half-space walk is A_1 <= a_1 + ... + a_k, so summing (F3) over all k and a_1 > ... > a_k with a_1 >= h,
  h^{>=h}_j <= sum_{A >= h} P_D(A) b_{j,A} <= P*_D(j) b^{>=h}_j,  and likewise  h_j <= P*_D(j) b_j.
Let gamma in SAW_n with max_k gamma_1(k) >= h and split it as in (F2). With x_1 = min gamma_1 <= gamma_1(0) = 0,
the first half-space walk has span max_{i>=m} gamma_1(i) - x_1 and the second has span
max_{i<=m} gamma_1(i) - x_1 + 1; since max_k gamma_1(k) - x_1 >= h, one of the two spans is >= h. By injectivity
of (F2),
  #{...} <= sum_{m=0}^{n} [ h^{>=h}_{n-m} h_{m+1} + h_{n-m} h^{>=h}_{m+1} ]
         <= sum_{m=0}^{n} P*_D(n+1)^2 [ b^{>=h}_{n-m} b_{m+1} + b_{n-m} b^{>=h}_{m+1} ]
         <= 2 (n+1) P*_D(n+1)^2 b^{>=h}_{n+1},
using that P*_D is nondecreasing and Lemma 0(a). For (3.1): the event {max_k ||gamma(k)||_inf >= h} is the union
over the 6 signed coordinate directions u of {max_k u . gamma(k) >= h}; each has the probability of
{max_k gamma_1(k) >= h} by the lattice symmetry sending u to e_1 (a bijection of SAW_n). Divide by
c_n >= mu^n (F1) and write mu^{-n} = mu z_c^{n+1}. QED

(The same reduction, with the Hardy–Ramanujan bound instead of Lemma 0(d), is used in DCH Section 2.4 and in
Krachun–Panagiotis, proof of Theorem 1; the loss e^{O(sqrt n log n)} here is harmless.)

## 4. The bridge decomposition and Theorem A

**Proposition 3.** Let N >= 1, m >= 1, h >= 2m, t >= 0, s := floor(h/m) - 1 (>= 1), and
  K_t := max{ 1, 4 z_c + 2 sum_{M>=1} e^{-2tM} B_M(z_c)^2 }   (possibly +infinity).
Then
  z_c^N b^{>=h}_N <= m^2 N e^{tN} (2 K_t)^{(N+1)/m} ( e^{-tm} B_m(z_c) )^{s}.

*Proof.* Put z = z_c e^{-t}; then z_c^N b^{>=h}_N = e^{tN} sum_{gamma} z^{|gamma|}, the sum over N-step bridges
gamma from 0 with H(gamma) >= h. Fix such a gamma, with height H >= h.

*Pigeonhole.* For i in {0, ..., m-1} let V_i(gamma) = #{u in [0, N] : gamma_1(u) = i mod m}. Since
sum_i V_i = N + 1, there is i with V_i(gamma) <= (N+1)/m. We bound the sum over gamma by the sum over i of the
sum over those gamma with V_i(gamma) <= (N+1)/m.

*Decomposition (i fixed).* Levels l_j = i + jm, j = 0, ..., s. Then 0 <= l_0 and
l_s <= m - 1 + (floor(h/m) - 1) m <= h - 1 <= H - 1. Since gamma goes from height 0 to H by steps changing the
height by at most 1, it visits every level in [0, H]. Define
  t+_j = max{u : gamma_1(u) = l_j},  t-_0 = min{u : gamma_1(u) = l_0},  t-_j = min{u > t+_{j-1} : gamma_1(u) = l_j}.
After t+_{j-1} the walk never returns to l_{j-1} and ends at H > l_{j-1}, so gamma_1(u) > l_{j-1} for all
u > t+_{j-1}; as it reaches H >= l_j, t-_j exists and t-_j <= t+_j. Hence
0 <= t-_0 <= t+_0 < t-_1 <= t+_1 < ... < t-_s <= t+_s < N (the last because gamma_1(N) = H > l_s). Pieces:
  alpha = gamma[0, t-_0],  c_j = gamma[t-_j, t+_j] (0 <= j <= s),  b_j = gamma[t+_j, t-_{j+1}] (0 <= j <= s-1),
  tau = gamma[t+_s, N],
each translated to start at 0, so that gamma = alpha o c_0 o b_0 o c_1 o ... o b_{s-1} o c_s o tau.
* alpha is a bridge of height l_0 = i in [0, m-1]: its start is the unique minimum of the bridge gamma, and its
  points before the first visit of level i lie strictly below i.
* b_j is a bridge of height exactly m: its points after the start are > l_j (after the last visit of l_j) and its
  points before the end are < l_{j+1} (before the first visit of l_{j+1} after t+_j); the end is at l_{j+1}.
* tau is a bridge of height H - l_s in [1, N]: points after t+_s are > l_s, and the end has height H, the maximum.
* c_j starts and ends at level l_j. Let I_j >= 1 be its number of visits to level l_j, at times
  t-_j = u_1 < ... < u_{I_j} = t+_j. Each subwalk gamma[u_r, u_{r+1}] has no visit to l_j in its interior: either
  it is one step inside the plane {x_1 = l_j} (4 possibilities) or it has length >= 2 and lies strictly above or
  strictly below l_j in its interior, i.e. it is an upper or a lower arch. All visits counted by the I_j are
  visits of gamma to levels = i mod m at distinct times, so sum_j I_j <= V_i(gamma) <= (N+1)/m.
The map gamma -> (alpha, (c_j), (b_j), tau), with each c_j further cut into its I_j - 1 segments, is injective
(concatenation recovers gamma), and z^{|gamma|} is the product of the z^{|piece|}. Dropping all mutual
avoidance constraints and the length constraint, and writing V = floor((N+1)/m),
  sum_{gamma : V_i(gamma) <= (N+1)/m} z^{|gamma|}
    <= [ sum_{T=0}^{m-1} B_T(z) ] . B_m(z)^s . [ sum_{T=1}^{N} B_T(z) ] .
       sum_{(I_0,...,I_s) in Z_{>=1}^{s+1}, sum I_j <= V}  prod_j (4z + 2A(z))^{I_j - 1}.
By Lemma 0(b),(c) the first two brackets are <= m and <= N, and B_m(z) <= e^{-tm} B_m(z_c). By Lemma 1 and
Lemma 0(c), 4z + 2A(z) <= 4 z_c + 2 sum_M e^{-2tM} B_M(z_c)^2 <= K_t. The number of tuples of positive integers
(I_0, ..., I_s) with sum <= V is binom(V, s+1) <= 2^V, and each product is <= K_t^V since K_t >= 1. Hence the
last factor is <= (2K_t)^V <= (2K_t)^{(N+1)/m}. Summing over the m values of i and multiplying by e^{tN} gives
the claim. QED

**Theorem A (unconditional).** For all integers n >= 1, m >= 1, h >= 2m and all real t >= 0,
  P_n( max_{0<=k<=n} ||gamma(k)||_inf >= h )
    <= 12 mu (n+1)^2 m^2 P*_D(n+1)^2 e^{t(n+1)} (2 K_t)^{(n+2)/m} ( e^{-tm} B_m(z_c) )^{floor(h/m) - 1},
with K_t as in Proposition 3 and P*_D(n+1) <= (n+2)^{sqrt(2n+2)+1}.

*Proof.* (3.1) and Proposition 3 with N = n + 1. QED

Theorem A is vacuous unless B_m(z_c) is small and K_t is finite; it converts information on *critical* bridges
into a displacement bound. The only lattice-specific inputs of Krachun–Panagiotis' Lemma 4.1 (their Theorem 2,
B_T <= 100 T^{-1e-10}, and their Lemma 2.2, sum_k G_k < infinity, both from the parafermionic observable on the
hexagonal lattice) are replaced by the quantities B_m(z_c) and sum_M e^{-2tM} B_M(z_c)^2.

## 5. Conditional corollaries

**Corollary B (hypothesis (H)).** Assume S = sum_{T>=1} B_T(z_c)^2 < infinity and put K = max{1, 4z_c + 2S},
lambda = 24 log(2K), c = log(2K)/4 (> 0.17). There is n_0 = n_0(S) such that for all n >= n_0,
  P_n( max_{0<=k<=n} ||gamma(k)||_inf >= lambda n / log n ) <= exp(-c n^{2/3}),
and therefore
  E_n ||gamma(n)||^2 <= 3 (lambda n/log n + 1)^2 + 3 n^2 exp(-c n^{2/3}) = (1 + o(1)) 1728 (log 2K)^2 n^2/(log n)^2.

*Proof.* Take t = 0, so K_0 <= K. Let n >= 8, M = ceil(n^{1/3}) (so n^{1/3} <= M <= 2 n^{1/3}). Since
sum_{M<=k<2M} B_k(z_c)^2 <= S, some m in [M, 2M) has B_m(z_c)^2 <= S/M. Let h = ceil(lambda n/log n); for n
large, h >= 4M >= 2m. Then s = floor(h/m) - 1 >= lambda n/(2M log n) - 2 and, with g := (1/2) log(M/S)
(positive for n large), log(1/B_m(z_c)) >= g >= (1/6) log n - (1/2) log S. Theorem A gives
  log P <= E(n) + ((n+2)/M) log(2K) - s g,
  E(n) := log(12 mu) + 2 log(n+1) + 2 log(2M) + 2 (sqrt(2n+2) + 1) log(n+2) = O(sqrt(n) log n).
Now s g >= lambda n g/(2M log n) - 2g >= lambda n/(12 M) - lambda n log^+(S)/(4 M log n) - log(M/S), and
lambda n/(12 M) = 2 log(2K) n/M. Hence
  log P <= - log(2K) n/M + [ 2 log(2K)/M + 6 log(2K) log^+(S) n/(M log n) + log(M/S) + E(n) ].
For n >= n_0(S) the bracket is <= (1/2) log(2K) n/M (it is o(n/M) since n/M >= n^{2/3}/2 and E(n) = O(sqrt n log n);
explicitly it suffices that log n >= 48 log^+ S and that the remaining terms are <= (3/8) log(2K) n/M). Thus
log P <= -(1/2) log(2K) n/M <= -(1/4) log(2K) n^{2/3}. For the second moment: ||x||^2 <= 3 ||x||_inf^2,
||gamma(n)||_inf <= n, and E||gamma(n)||_inf^2 <= h^2 + n^2 P(||gamma(n)||_inf >= h). QED

Since B_T(z_c) <= C_1 T^{-a} with a > 1/2 implies S <= C_1^2 zeta(2a) < infinity, the same conclusion holds under
that polynomial-decay hypothesis.

**Corollary B2 (borderline hypothesis).** Assume that sum_{M=1}^{L} B_M(z_c)^2 <= C_2 (1 + log L)^p for all
L >= 1, some C_2 > 0, p > 0 (this includes B_T(z_c) <= C_1 T^{-1/2}, with p = 1). Then there are C, c > 0 and n_0
such that P_n(max_k ||gamma(k)||_inf >= C n log log n / log n) <= exp(-c n^{2/3}) for n >= n_0; hence
E_n||gamma(n)||^2 = O(n^2 (log log n)^2/(log n)^2).

*Proof.* Let n >= 8, M = ceil(n^{1/3}), and choose m in [M, 2M) with
B_m(z_c)^2 <= (1/M) sum_{k<2M} B_k(z_c)^2 <= C_2 (1 + log 2M)^p / M. Then
  log(1/B_m(z_c)) >= g0 := (1/6) log n - (p/2) log(1 + log(4 n^{1/3})) - (1/2) log C_2 = (1/6 - eps_n) log n,
with eps_n = O(log log n / log n). Apply Theorem A with t = 1/m: e^{t(n+1)} (2K_t)^{(n+2)/m} <= (2e K_t)^{(n+2)/m},
and the factor e^{-tms} <= 1 is dropped. Splitting the summation index k into blocks [jm, (j+1)m), on which
e^{-2k/m} <= e^{-2j}, and using log((j+1)m) <= log m + j and 1 + a + b <= (1+a)(1+b),
  sum_{k>=1} e^{-2k/m} B_k(z_c)^2 <= sum_{j>=0} e^{-2j} C_2 (1 + log m + j)^p <= C_2 c_p (1 + log m)^p,
  c_p := sum_{j>=0} e^{-2j} (1+j)^p,
so Lambda_n := log(2e K_{1/m}) <= p log(1 + log(2M)) + C_3 = O(log log n), and Lambda_n >= log(2e) > 1.
Let h = ceil(24 Lambda_n n / log n) (so h = O(n log log n/log n) and h >= 2m for n large). Then
s = floor(h/m) - 1 >= 12 Lambda_n n/(M log n) - 2 and s log(1/B_m) >= s g0, so by Theorem A
  log P <= E(n) + Lambda_n (n+2)/M - (12 Lambda_n n/(M log n) - 2) (1/6 - eps_n) log n
        <= - Lambda_n (n/M) (1 - 12 eps_n) + 2 Lambda_n/M + (1/3) log n + E(n),
with E(n) = O(sqrt(n) log n) as in Corollary B. For n large, 12 eps_n <= 1/4 and the last three terms are
<= (1/4) n/M <= (1/4) Lambda_n n/M, so log P <= -(1/2) Lambda_n n/M <= -(1/4) n^{2/3}. The second-moment bound
follows as in Corollary B. QED

**Proposition 5 (generating-function version: mass).** Let m(z) be the mass of the two-point function; for
0 < z < z_c, m(z) = lim_{L->inf} -(1/L) log B_L(z) [Madras–Slade Proposition 4.1.8 (existence of the limit
Mbar(z) for the point-to-plane bridge generating function B_z(L) = B_L(z)) and Theorem 4.1.14 (m = Mbar on
(0, z_c))]. (Notation clash: m(.) is the mass, the integer m is the block height.) Let t > 0, theta >= 1 and
m >= 1 be such that z' := z_c e^{-t} theta^{1/m} < z_c and
theta >= 1 + 4z' + 2A(z'). Then
  m(z_c e^{-t}) >= (1/m) log(1/B_m(z')) >= (1/m) log(1/B_m(z_c)).

*Proof.* Put z = z_c e^{-t} and fix L >= 2m. For a bridge gamma of height exactly L (any length) choose i with
V_i(gamma) <= (|gamma|+1)/m (pigeonhole as in Proposition 3); then theta^{(|gamma|+1)/m - V_i(gamma)} >= 1, so
  B_L(z) <= sum_i sum_gamma z^{|gamma|} theta^{(|gamma|+1)/m - V_i(gamma)}
         = theta^{1/m} sum_i sum_gamma (z')^{|gamma|} theta^{-V_i(gamma)}.
Decompose gamma exactly as in Proposition 3 (with h = L, s = floor(L/m) - 1) and use theta^{-V_i} <=
theta^{-(I_0 + ... + I_s)}. The c_j-factor becomes, for each j, sum_{I>=1} theta^{-I} (4z' + 2A(z'))^{I-1}
= 1/(theta - 4z' - 2A(z')) <= 1. Hence
  B_L(z) <= theta^{1/m} m . m . B_m(z')^{floor(L/m)-1} . B(z'),   B(z') = sum_T B_T(z') < infinity (z' < z_c).
Take -(1/L) log and let L -> infinity. The last inequality is B_m(z') <= B_m(z_c). QED

**Corollary C (mass under (H)).** Assume (H) and put theta_0 = 1 + 4z_c + 2S. There is t_0(S) > 0 such that for
0 < t <= t_0,
  m(z_c e^{-t}) >= t log(1/t) / (9 log theta_0),
i.e. the correlation length xi(z) = 1/m(z) satisfies xi(z) (z_c - z) -> 0 at rate O(1/log(1/(z_c - z))).

*Proof.* By Lemma 1, 1 + 4z' + 2A(z') <= theta_0 for every z' <= z_c, so theta = theta_0 is admissible as soon as
m > log(theta_0)/t. Let M = ceil(2 log(theta_0)/t); every m in [M, 2M) qualifies, and one of them has
B_m(z_c)^2 <= S/M. Proposition 5 gives m(z_c e^{-t}) >= (1/(2M)) (1/2) log(M/S) = log(M/S)/(4M). As t -> 0,
M ~ 2 log(theta_0)/t and log(M/S) = log(1/t) + O(1), so the right side is (1 + o(1)) t log(1/t)/(8 log theta_0),
which is >= t log(1/t)/(9 log theta_0) for t small. QED

(Consistency check, not needed elsewhere: Corollary C with the elementary inequality of `attempts.md` C1,
P_n(||gamma(n)||_inf >= L) <= 6 z^{-2} e^{tn + 8/t} e^{-m(z) L}/(1 - e^{-m(z)}), at t = n^{-1/2} again gives an
endpoint bound at scale L ~ C n/log n.)

**Remark (what the method cannot give; heuristic).** (i) The gain per bridge block of height m is
log(1/B_m(z_c)); a polynomial rate (U2) from Theorem A would need log(1/B_m(z_c)) >> log m, i.e. B_m(z_c) decaying
faster than any power along usable m. Section 6, item 4 shows that B_T(z_c) is not o(1/T) on average (Abel sense);
isolated very small values are not excluded rigorously, but nothing suggests them. So the method is expected to
stop at n/log n, as in Krachun–Panagiotis. (ii) If only
B_T(z_c) <= C T^{-a} with a < 1/2 is known, Theorem A with t = 1/m gives (by the computation of Corollary B2, with
K_t <= C' m^{1-2a}) only a linear-speed bound with v > (1 - 2a)/a + o(1) — nontrivial only for a > 1/3, and not
sub-ballisticity. This remark is a computation, not a verified claim.

## 6. What is known about hypothesis (H), and why it is plausible in d = 3 only

1. (PROVED from CITED.) Renewal structure: every bridge is uniquely a concatenation of irreducible bridges and
   sum_{gamma in iSAB} mu^{-|gamma|} = 1 (DCH Lemma 2.2, Kesten; Madras–Slade (4.2.3)-(4.2.4); the height
   version of the renewal equation is Madras–Slade (4.2.9), B_z(L) = sum_{j=1}^L A_z(j) B_z(L-j) + delta_{L,0},
   valid for every z > 0 as an identity of power series with nonnegative coefficients, in particular at z_c). Under the law P_iSAB(gamma) = mu^{-|gamma|} the
   height H of an irreducible bridge is an a.s. finite positive integer, and expanding a bridge of height T into
   its irreducible pieces gives
     B_T(z_c) = sum_k P(H_1 + ... + H_k = T) =: u_T,   H_i i.i.d. ~ P_iSAB(H in .),
   i.e. B_T(z_c) is the renewal sequence of the heights. Since P_iSAB(H = 1) >= P_iSAB(gamma = one step) =
   mu^{-1} > 0, the Renewal Theorem in the form of Madras–Slade Theorem 4.2.2(b) (hypothesis f_1 > 0, f = 1;
   here f_k = P_iSAB(H = k), g_n = delta_{n,0}) gives B_T(z_c) -> 1/E_iSAB[H] (= 0 if the mean is infinite). Hence **B_T(z_c) -> 0 iff E_iSAB[H] = infinity**, and
   (H) implies both. (DCH Theorem 2.5 gives E_iSAB[|gamma|] = infinity, the statement about *lengths*; the
   statement about *heights* is stronger, since H <= |gamma|. For lengths, Madras–Slade (1993, p. 92) still wrote
   "it is believed that lim_N b_N/mu^N = 0, but there is no known proof of this"; DCH settled that, without rate.)
2. (Elementary.) sum_T u_T^2 is the expected number of common renewal heights of two independent copies of the
   height renewal process; so (H) implies that two independent critical infinite bridges (DCH's
   P_iSAB^{tensor N}) a.s. have finitely many common cut-plane heights (and conversely, since the common renewal
   heights form a renewal process, which is a.s. finite only if its expected number of points is finite). Equivalently (Lemma 0(a) is an identity here) sum_T B_T(z_c)^2 is the critical weight of
   bridges of height 2T having a renewal point at height T, summed over T.
3. (Open as far as I could find.) I did not find a proof of B_T(z_c) -> 0 on Z^d, d = 2, 3, 4, in the sources
   read (DCH 2013, DGHM 2016, DCGHM 2020, Hutchcroft 2018, Hammond 2019 x2, BDGS lectures, Glazman–Manolescu 2020,
   Krachun–Panagiotis 2023). On the hexagonal lattice, B_T -> 0 is [BBMdG+14, Thm 10] / Glazman–Manolescu, and
   B_T <= 100 T^{-1e-10} is Krachun–Panagiotis Theorem 2 (parafermionic observable). Krachun–Panagiotis (p. 2,
   arXiv v1): "A major obstacle in that direction seems to be the lack of a quantitative bound on the partition
   function of bridges of a given height at the critical point x_c = 1/µ that goes to 0 as the height goes to
   infinity."
4. (CITED, constraint.) Madras–Slade (3.1.14), as restated in Hutchcroft Remark 2.5: B((1-e) z_c) >=
   1 + (1/2) log[(1-e) z_c] + (1/2) log(1/e). With Lemma 0(c) (z = z_c e^{-t}, e = 1 - e^{-t} <= t),
   sum_T e^{-tT} B_T(z_c) >= B(z) >= (1/2) log(1/t) - O(1). Hence B_T(z_c) <= eps/T for all large T is impossible
   for eps < 1/2. So a power law B_T(z_c) ~ T^{-a} compatible with (H) must have a in (1/2, 1].
5. (HEURISTIC.) Scaling: if B_T(z) ~ T^{-a} g(T/xi(z)) with xi ~ (z_c - z)^{-nu}, then
   B(z) = sum_T B_T(z) ~ (z_c - z)^{-nu(1-a)}, so gamma_b = nu (1 - a), a = 1 - gamma_b/nu, where
   b_n ~ mu^n n^{gamma_b - 1}; equivalently a = -gamma_11/nu with the (physics) scaling relation
   gamma_b = gamma_11 + nu of Duplantier–Guttmann (arXiv:1908.03872, abstract). In d = 2, gamma_b = 9/16, nu = 3/4 give a = 1/4, matching Krachun–Panagiotis'
   expectation B_T ~ T^{-1/4}; then sum B_T^2 = infinity. In d = 3, Clisby–Conway–Guttmann (J. Phys. A 49 (2016)
   015004, arXiv:1504.02085; they also give the scaling relation gamma_b = gamma_11 + nu, their eq. (6)) estimate
   gamma_b = 0.199 +- 0.002 (and gamma_11 = -0.3874 +- 0.0024); with nu = 0.587597 this gives
   a ~ 0.661, 2a ~ 1.32 > 1: **(H) is expected to hold in d = 3**, with some margin, and the conditional theorem
   is expected to be "true but unproved input" in d = 3, and false-input in d = 2. In mean field (d >= 5) a = 1.
6. My own exact enumeration (n <= 17, `attempts.md`, Section G) confirms Lemma 1 coefficientwise but is far too
   short to estimate a.

## 7. Novelty assessment (honest)

* Theorem A / Corollary B: the architecture is that of Krachun–Panagiotis, Lemma 4.1 and the proof of their
  Theorem 1 (hexagonal lattice), transplanted to Z^3. What is new, as far as I can tell from the sources listed
  in Section 6.3: Lemma 1, which bounds the half-space arch generating function (the analogue of KP's
  sum_k G_k, which KP obtain from the parafermionic observable and say is "only known in the case of the
  hexagonal lattice among low dimensional lattices") by sum_M B_M(z)^2; and the resulting observation that in
  d = 3 the KP route needs **one** critical input, (H), which is heuristically true in d = 3 and false in d = 2.
  Lemma 1 itself is elementary and its technique is standard: Glazman–Manolescu (arXiv:1708.00395v3, Section 5.4,
  proof of Theorem 3) split a walk at its first visit to its extreme column into two pieces that are (after local
  modification) bridges and bound the result by delta sum_T B_{T}(x, y)^2; Madras–Slade Lemma 4.1.12 and
  Hammersley–Welsh use related reflections. What I have not seen is its use as a substitute for KP's Lemma 2.2 and
  the resulting single-hypothesis statement on Z^3. Madras–Slade Ch. 3-4 (owner-supplied excerpts) contain no
  statement on B_{z_c}(L) as L -> infinity. Not read in full: Beaton et al. 2014, Glazman–Manolescu.
* This is **not** U1 or U2: no explicit rate for E_n|omega_n|^2 on Z^3 is proved unconditionally. DCH's
  Question 1 and the literature's "it is not proved that R_n <= C n^{0.99999}" (Slade 2019, Sec. 2.1.3, via
  literature/quick-facts.md A5) stand.
* No progress on E1.
