# Lower-bound solver: candidate finite-n statements for the computation lab

Solver: Claude Opus 5.5 (effort xhigh), 2026-10-24.
Context: `result.md` (Theorem 3: the inputs (H1) and (H2) give E_N|omega_N|^2 >= (R^2/(4C_0))^{2/(3-alpha)}) and
`attempts.md`. Notation: Madras's map T, span, near set S_{N,<=m} = {||omega_N||_inf <= m}, multiplicity
K(psi) = |T^{-1}(psi) ∩ S_{N,<=m}|.

**General caveat on scales.** Every statement whose content is asymptotic in m or in N^{1/3} is uninformative at
exact-enumeration lengths. For n <= 36 we have N^{1/3}/3 <= 1.1, so Madras's "far" threshold is trivial and the
relevant m is 1. Such statements can only be probed by Monte Carlo (pivot algorithm, N = 10^3 to 10^5), which is
HEURISTIC. My own MC is `scripts/madras_multiplicity_mc.cpp` (outputs in `data/`). Exact tests are useful only for
CL3 and as sanity checks.

## CL1. Sublinear average multiplicity of Madras's map (input (H2))

Statement. There are C, alpha > 0 such that for all N and all integers 1 <= m < span-threshold,

    kappa_N(m) := (1/|S_{N,<=m}|) * sum_{omega in S_{N,<=m}} K(T omega) <= C m^{1-alpha}.

- If true: exponent 4/(d(3-alpha)) > 4/(3d) (`result.md` Thm 3).
- If kappa_N(m) <= C, then E >= c N^{2/d}.

Informative range: N >= 10^3, 8 <= m <= N^{0.5}, with m well below the typical span (check the fraction of near walks
with span <= m; the MC prints it).

My MC (HEURISTIC): kappa = 2.72 / 3.07 / 3.88 at m = 16 / 32 / 64 for N = 4000, and 1.95 / 2.59 / 3.19 at
m = 8 / 16 / 32 for N = 1000. Requested: an independent pivot run with autocorrelation estimates, N up to about
3*10^4 and m up to about 256, and a fit of log kappa against log m. Also report the two partial counts:
- levels that pass the renewal test;
- levels that also pass the self-avoidance test.
These separate the two mechanisms of Proposition 2.

Exact sanity check (n <= 19): with m = 1 and span >= 2, the fraction of near walks with K = 2, as a function of n.

## CL2. Renewal-depth density at the top (ingredient of CL1)

Statement. Fix the Madras data (L, I, H > 0). Let r_N(j) = P_N(the tail renews at depth j) and
r'_N(j) = P_N(the head renews at depth j), with the definitions of `result.md` Section 2. Conjecture:
(1/m) * sum_{j<=m} (r_N(j) + r'_N(j)) -> 0 as m -> inf, uniformly in N >= N_0(m).

- If true (with any rate) the multiplicity is o(m), giving a non-quantitative improvement E/N^{4/(3d)} -> inf.
- A rate m^{-alpha} gives CL1 directly, even ignoring collisions.

My MC (HEURISTIC): the density is about 0.4 at m = 16 and about 0.25 at m = 64 (N = 4000).
Requested: r_N(j) as a function of j at N = 10^4, on a log-log plot. Is it ~ j^{-beta}, and with what beta?

Relation to the literature: for bridges this is the height-density of renewals, i.e. B_T(z_c) -> 0. That is open on
Z^d and known on the hexagonal lattice (Krachun-Panagiotis). A numerical estimate of B_T(z_c) decay on Z^3 from the
lab's bridge tables would also be useful.

## CL3. One-step monotonicity: two inequalities whose conjunction gives E_{n+1} - E_n >= 1, hence L3

PROVED in `attempts.md` A5: E_n - E_{n-1} >= 1 follows from
- (I1) E_n|omega_{n-1}|^2 >= E_{n-1}|omega_{n-1}|^2, i.e. Cov_{P_{n-1}}(|eta_end|^2, #free neighbours of eta_end) >= 0;
- (I2) E_n[omega_n · omega_1] >= 1, i.e. E_n[omega^{(1)}(n) | omega(1) = e_1] >= 1 (mean reflection principle).

The lab already reports that E_{n+1} - E_n >= 1 holds for n <= 35, that E_n|omega_k|^2 >= E_k|omega_k|^2 holds for
n <= 19 (which includes (I1) with k = n-1), and that E_n[omega_k·(omega_n - omega_k)] >= 0 holds for n <= 19 (which
includes (I2): with k = n-1 this is E_n[omega_{n-1}·s_n] >= 0, equivalent to (I2) by reversal).

Requested:
- the margins, as exact rationals, for n <= 19 or more:
  - Delta_1(n) = E_n|omega_{n-1}|^2 - E_{n-1}|omega_{n-1}|^2;
  - Delta_2(n) = E_n[omega_n·omega_1] - 1;
- whether each margin is increasing in n (a margin growing like n^{2nu-1} would be the expected behaviour);
- whether the stronger per-walk forms fail, e.g. whether the conditional version of (I2) given the first k steps fails
  for small k. Failures tell us which averaging a proof must use.

## CL4. Critical two-point function at a neighbour (GF route)

Statement: G_{z_c}(e_1) = (mu/3) * sum_m m p_m mu^{-m} < inf in d = 3.
- If true, together with sup_x G_{z_c}(x) <= C, it gives a subsequence version of nu >= 1/d (`attempts.md` A4).
- Open (Madras-Slade p. 18).

Requested (HEURISTIC): partial sums S(M) = sum_{m <= M} m p_m mu^{-m} with a numerical estimate of mu for Z^3 (about 4.684; take the lab's reference value),
using the polygon counts in certificates/sources/lacecounts/poly_d3.txt (or the lab's own), and the local exponent of
the increments m p_m mu^{-m} (expected ~ m^{-1.76}, so convergent). Also: is max_x G_{z}(x) attained at a neighbour of
0 for z close to z_c? (The lab found that for fixed n the argmax of c_n(x) drifts outward, so the analogous GF
question is non-trivial.)

## CL5. Span lower bound (input (H1))

Statement: P_N(span(omega) < N^{1/3 + eta}) <= 1/4 for some eta > 0 and all large N. True value of the span exponent
is nu, about 0.588. The statement is certainly true numerically. What matters is a proof.

Requested (HEURISTIC, cheap): the distribution of span/N^{1/3} at N = 10^3 and 10^4, to document how far the volume
bound is from the truth. No proof strategy is attached; see `attempts.md` A3.
