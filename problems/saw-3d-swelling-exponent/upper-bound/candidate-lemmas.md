# Upper-bound solver: candidate finite-n statements for the computation lab

Context: `result.md` proves (Corollary B) that the single hypothesis
  (H)  S := sum_{T>=1} B_T(z_c)^2 < infinity,  B_T(z) = sum_n b_{n,T} z^n,
implies E_n|omega_n|^2 = O(n^2/(log n)^2) on Z^3. Heuristically B_T(z_c) ~ T^{-a} with a = 1 - gamma_b/nu ~ 0.66.
Here b_{n,T} = number of n-step bridges from 0 (BDGS Def. 2.1: 0 = x1(w_0) < x1(w_i) <= x1(w_n)) with
x1(w_n) = T, as in `computation/data/bridge_span.tsv`; arches a_{n,M} = number of n-step SAWs from 0 with
x1(w_n) = 0, x1(w_i) > 0 for 0 < i < n, and max_i x1(w_i) = M. Floating point is fine for these tests: none of
them can be a certificate (they are statements about all n or about infinite sums).

## CL1 (exact sanity check of result.md Lemma 1; should HOLD for every n)

  a_{n,M} <= sum_{k=M}^{n-M} b_{k,M} b_{n-k,M}   for all n >= 2, M >= 1.

Informative range: all n <= 19 (I checked n <= 17 with my own enumerator; the lab's tables to n = 19 would make the
check independent of my code). Also report the ratio sum_M a_{n,M} / sum_M sum_k b_{k,M} b_{n-k,M} (mine: 0.29 at
n = 4 decreasing to 0.117 at n = 17). A failure would refute Lemma 1 and Corollary B.

## CL2 (exact identity, sanity check of the renewal structure used in result.md Sec. 6)

Let r_{n,T} = number of n-step bridges of height 2T having a renewal point at height T (an index i with
w_1(j) <= w_1(i) = T < w_1(k) for j <= i < k, Madras–Slade (4.2.1)). Then

  r_{n,T} = sum_k b_{k,T} b_{n-k,T}   for all n, T.

Informative range: n <= 19, T <= 9. Should HOLD exactly.

## CL3 (numerical estimate of the decay exponent of B_T(z_c); HEURISTIC input for (H))

For each T = 1, ..., 8: from the exact b_{n,T}, n <= 19 (or further with transfer matrices / series in the lab),
extrapolate the tail of sum_n b_{n,T} z_c^n (e.g. fit b_{n,T} z_c^n ~ C_T n^{g_T} for large n and add the tail),
z_c = 1/4.684039931 (Clisby 2013 estimate of mu). Then fit B_T(z_c) ~ C T^{-a}.
Question: is a > 1/2 supported (prediction a = 1 - 0.199/0.5876 ~ 0.661, from Clisby–Conway–Guttmann 2016)?
Also report the truncated sums at N = 19 (mine at N = 17: 0.697, 0.544, 0.405, 0.271, 0.159, 0.080 for T = 1..6).
Informative only if the extrapolated values are stable in N; expect large uncertainty for T >= 4.

## CL4 (monotonicity of the critical bridge partition function in the height; structural candidate)

  B_{T+1}(z) <= B_T(z)  for all T >= 0 and 0 < z <= z_c.

On the hexagonal lattice this is Krachun–Panagiotis Lemma 2.3 (parafermionic observable). On Z^3 it is not
known to me. If true, (H) would reduce to a Cesaro bound sum_{T<=L} B_T(z_c) <= C L^{1/2-delta}, i.e. (by the
renewal/Wald identity, attempts.md A3, applied to heights) to E_iSAB[H wedge L] >= c L^{1/2+delta}.
Test: with the exact polynomials truncated at n <= 19, compare B_T^{(19)}(z) and B_{T+1}^{(19)}(z) at
z in {0.10, 0.15, 0.18, 0.20}. For z <= 0.15 the truncation error is at most
sum_{n>19} c_n z^n <= sum_{n>19} 6 * 5^{n-1} z^n (c_n <= 2d(2d-1)^{n-1}), small enough to make the comparison
meaningful for T <= 6; near z_c it is not controlled. A clear violation at small T and z <= 0.15 would refute CL4.

## CL5 (exponential-scale large deviation of bridges; context for the DCH route)

phi_n(v) := -(1/n) log( #{n-step bridges with x1(w_n) >= v n} / b_n ), n <= 19, v in {0.3, 0.4, ..., 0.9}.
Heuristic (Fisher's law) phi(v) ~ c v^{1/(1-nu)} ~ c v^{2.43}; Hutchcroft Remark 2.5 rigorously forces
phi(v) <= (1+o(1)) 2v/log(1/v) as v -> 0. Only of interest if the curves phi_n(v) stabilise in n.
