# Quick facts: rigorous bounds for strict SAW on Z^d, d = 2, 3, 4 (state as of 2026-10-24)

Literature agent. Every statement below was read in the source named (version stated). "Read in full" = whole
paper read; "read: statements + proof overview" = theorem statements, definitions and proof architecture read,
not every line of every proof; "secondary" = taken from another source that I read, primary not read.

Notation: Z^d nearest-neighbour lattice, strict (not weakly) self-avoiding walks from 0, c_n their number,
E_n / P_n uniform measure on n-step SAWs, mu the connective constant (mu^n <= c_n by submultiplicativity + Fekete),
z_c = 1/mu, b_n number of n-step bridges, "closing" = ||omega_n|| = 1.

## A. Mean-square end-to-end distance E_n |omega_n|^2

| # | Statement (exact) | Hypotheses | Source | Status of my reading |
|---|---|---|---|---|
| A1 | For every positive integer N and every p in [1, inf): (E_N ||omega(N)||_p^p)^{1/p} >= (E_N ||omega(N)||_inf^p)^{1/p} >= (1/6) N^{p/((p+1)d)}. With p = 2: **E_N ||omega(N)||_2^2 >= N^{4/(3d)} / 36 for ALL N >= 1** (d = 3: N^{4/9}/36; d = 2: N^{2/3}/36; d = 4: N^{1/3}/36). | Z^d, d >= 2 fixed, strict SAW, uniform measure, every N (no subsequence), explicit constant | Madras, Canad. Math. Bull. 57 (2014) 113-118, **Proposition 1.1** (p. 114), doi:10.4153/CMB-2012-022-6 | Read in full (publisher PDF, 6 pages; owner also supplied the same publisher PDF) |
| A1' | Same Proposition, large p: (E_N ||omega(N)||_inf^p)^{1/p} >= N^{p/((p+1)d)}/6, exponent -> 1/d as p -> inf. Gives nothing better than n^{4/(3d)} for the SECOND moment. | as A1 | idem | idem |
| A2 | Sub-ballisticity: for each v > 0 there is eps > 0 such that, for each n, P_n( max_{0<=k<=n} ||gamma_k|| >= v n ) <= e^{-eps n}. Corollary 1.2: lim_{n} n^{-2} E_n ||gamma_n||^2 = 0. **No rate**: eps = eps(v) is not explicit, so no o(n^2) rate follows. | Z^d, every d >= 2, strict SAW, all n | Duminil-Copin, Hammond, CMP 324 (2013) 401-423, **Theorem 1.1, Corollary 1.2**, doi:10.1007/s00220-013-1811-1 | Read in full, arXiv:1205.0401v1 (only arXiv version; published version not compared) |
| A2' | Bridge form (as stated by Hutchcroft): there is an increasing phi_d : (0,1] -> (0,inf) with #{n-step bridges reaching height >= m} <= b_n exp(-phi(m/n) n) for all n >= m > 0. "The proof of Theorem 1.3 is not quantitative and does not give any estimates on the function phi." | Z^d, d >= 2 | Hutchcroft, ECP 23 (2018), **Theorem 1.3** (restating DCH Cor. 2.4 + Thm 2.5), arXiv:1708.09460v3 | Read in full |
| A3 | Quantitative sub-ballisticity, **hexagonal lattice only**: there is c > 0 such that P_n( max_k ||gamma_k|| >= 10^11 n / log(n+1) ) <= exp(-c n^{2/3}) for every n >= 1. Immediate consequence (my one-line deduction, not stated in the paper): E_n||gamma_n||^2 <= C n^2/(log n)^2 on H. | hexagonal lattice H (2D), strict SAW, all n | Krachun, Panagiotis, Ann. Probab. 54 (2026) 1109-1125, **Theorem 1**, doi:10.1214/24-AOP1730 | Read in full, arXiv:2310.17299v1 (published version not compared) |
| A4 | Critical bridges on H: B_T <= 100 T^{-eps}, eps = 10^{-10}, for every T >= 1 (B_T = partition function at x_c of bridges of height T). Previously only B_T -> 0 (Glazman-Manolescu 2020) and B_T >= c/T (Duminil-Copin-Smirnov). | H only (parafermionic observable) | Krachun-Panagiotis, **Theorem 2** | idem |
| A5 | Nothing quantitative is known on Z^d: "The upper bound states that R_n/n -> 0 but there is no bound on the rate. In particular, it is not proved that there is a constant C such that R_n <= C n^{0.99999}." (R_n^2 = E_n||omega(n)||_2^2) | d >= 2 | Slade, Proc. R. Soc. A 475 (2019) 20180549, Section 2.1.3 (p. 5 of arXiv:1808.04476v2) | Read (Section 2 in full) |
| A6 | d >= 5 only (context): R_n ~ D n^{1/2}, c_n ~ A mu^n (Hara-Slade 1992). | d >= 5 | Slade 2019 Theorem 2.1; BDGS lectures Section 1.6.5 | Secondary (Hara-Slade not read) |

**Nothing newer found** (as of this file): no improvement of Madras' exponent 4/(3d) on Z^d, and no quantitative
sub-ballisticity on Z^d for any d in {2, 3, 4}. Searches: Semantic Scholar citation lists of Madras 2014 (14 entries)
and DCH 2013 (45 entries, through 2026); arXiv listing of all math.PR and math-ph abstracts containing
"self-avoiding" since 2024-01-01; arXiv phrase searches 2023-2026 (details in `novelty.md`). Absence of results is not
proof of absence. OpenAlex was rate-limited (shared IP budget exhausted) at the time of writing.

Open (verbatim sources in `open-statements.md`): E_n|omega_n|^2 >= c n (L3), E_n|omega_n|^2 >= c n^{2/d} (even this
is open, Madras 2014 p. 114, BDGS (1.30), DCH Question 4), E_n|omega_n|^2 <= C n^{2-eps} (U2; DCH Question 1, BDGS
(1.30), Slade 2019), existence of nu.

## B. Probability that the walk ends at a given point

| # | Statement | Hypotheses | Source | Reading |
|---|---|---|---|---|
| B1 | For any eps > 0 and n large enough, P_n(||Gamma_n|| = 1) <= n^{-1/4+eps}. | Z^d, d >= 2, all large n | Duminil-Copin, Glazman, Hammond, Manolescu, Ann. Probab. 44 (2016) 955-983, **Theorem 1.1**, doi:10.1214/14-AOP993 | Read: statements, Sections 1-3 in full, proof overviews of Sections 4-5 (arXiv:1305.1257v2) |
| B2 | sup_{x in Z^d} P_n(Gamma_n = x) -> 0 as n -> inf (non-quantitative). | Z^d, d >= 2 | DGHM 2016, **Theorem 1.2** | idem |
| B3 | Midpoint: there is C with sup_x P_n(Gamma_{floor(n/2)} = x) <= C n^{-1/2} for all n. | Z^d, d >= 2 | DGHM 2016, **Proposition 1.3** | idem |
| B4 | Elementary (Madras-type single unfolding): P_n(||Gamma_n|| = 1) <= 2/3 for all n >= 3d+1. | Z^d | DGHM 2016, Section 2.3 | idem |
| B5 | **Best closing bound in d = 3**: for any eps > 0 and odd n sufficiently high, P_n(Gamma closes) <= n^{-1/2+eps}. ("snake method via Gaussian pattern fluctuation") | Z^d, d >= 2, all large odd n | Hammond, Trans. AMS 372 (2019) 2335-2356, **Theorem 1.2**, doi:10.1090/tran/7494 | Read: abstract, intro, statements (arXiv:1808.09597v1) |
| B6 | d = 2 only: for any eps > 0, P_n(Gamma closes) <= n^{-4/7+eps} on a set of odd n of limit supremum density at least 1/1250 (Thm 1.5(1)); conditional Thm 1.5(2): if theta = lim theta_n and xi = lim xi_n exist then theta + xi >= 5/3. | Z^2 | Hammond, EJP 24 (2019) 1-43, **Theorem 1.5**, doi:10.1214/18-EJP249 | Read: abstract, Section 1 (arXiv:1808.10500v3) |
| B7 | Identity: P_n(Gamma closes) = 2(n+1) p_{n+1} / c_n (p_m = polygons of length m up to translation). | Z^d | Hammond EJP 2019 eq. (1.1) | idem |
| B8 | Not proved (expected): P_n(Gamma_n = x) <= n^{-1/4+eps} for a fixed general x ("For higher dimension, significant difficulties occur and it is unclear how to prove this result.") | | DGHM 2016, after Theorem 1.1 | idem |

## C. c_n / mu^n and related

| # | Statement | Hypotheses | Source | Reading |
|---|---|---|---|---|
| C1 | mu^n <= c_n for all n. | Z^d | submultiplicativity; e.g. BDGS (1.25), Slade 2019 | standard, read |
| C2 | Hammersley-Welsh: if gamma > pi sqrt(2/3) there is n_0(gamma) independent of d with 0 <= kappa(n) - kappa < gamma n^{-1/2} + n^{-1} log d (n >= n_0), where c_n = e^{n kappa(n)}, kappa = log mu; i.e. c_n <= d mu^n e^{gamma sqrt n}. (BDGS Thm 2.6 form: for B > pi(2/3)^{1/2}, n >= n_0(B) independent of d: c_n <= b_{n+1} e^{B sqrt n} <= mu^{n+1} e^{B sqrt n}.) | Z^d, d >= 2 | Hammersley-Welsh, Quart. J. Math. 13 (1962) 108-110, **Theorem** (eq. (3)); BDGS lectures Thm 2.6 | Read in full (owner-supplied scan, 3 pages) and BDGS proof read |
| C3 | **Kesten (best in d = 3)**: there exist alpha_6(d), alpha_7(d) with beta_d^n <= chi_n(d) <= exp(alpha_6 n^{2/(d+2)} log n) beta_d^n (1.8), and exp(-alpha_7 n^{2/(d+1)} log n) beta_d^{2n-1} <= gamma_{2n-1}(d) <= 4(d-1) n beta_d^{2n-1} (1.9) [chi_n = c_n, beta_d = mu, gamma_{2n-1} = number of closing walks of 2n-1 steps]. Kesten: "The right-hand side of (1.8) improves (1.6) [Hammersley-Welsh] for all d >= 3." In d = 3: **c_n <= mu^n exp(C n^{2/5} log n)**. Kesten only sketched the proof ("it does not seem worthwhile to give more than an indication of the proof of Theorem 2"); a **full proof is in Madras-Slade, Theorem 3.3.1** (p. 68: "Let d >= 2. Then there exists a constant Q, depending only on d, such that for every N >= 2, c_N <= mu^N exp[Q N^{2/(d+2)} log N]"; proof due to Kesten, private communication; read: statement and pp. 68-70 of the owner-supplied excerpt). | Z^d; useful for d >= 3 | Kesten, J. Math. Phys. 5 (1964) 1128-1137, **Theorem 2**, doi:10.1063/1.1704216 | Read: statements and intro from owner-supplied scan (page images checked) |
| C4 | Kesten ratio limit: |c_{n+2}/c_n - mu^2| <= A_3 n^{-1/3} (Theorem 4, eq. (3.8)/(1.6)); -A_4 n^{-1/3} <= gamma_{2n+3}/gamma_{2n+1} - mu^2 <= A_5 n^{-1/4} (3.9). Pattern theorem = Theorem 1 (lim sup [chi_n(an,(C,D))/chi_n]^{1/n} < 1). Also e^{-alpha_4 n^{1/2}} beta^n <= gamma_n (1.3). | Z^d, d >= 2 | Kesten, J. Math. Phys. 4 (1963) 960-969, doi:10.1063/1.1704022 | Read: Sections 1-3 statements from owner-supplied scan (page images checked for exponents) |
| C5 | Kesten 1964 Corollary 2 (stated without proof): -A_6 (log n)^{(d+1)/2d} n^{-(d-1)/2d} <= c_{n+2}/c_n - mu^2 <= A_7 (log n)^{(d+2)/(2d+2)} n^{-d/(2d+2)} (1.10). | Z^d | Kesten 1964, Corollary 2 | as C3 |
| C6 | Hutchcroft: c_n <= exp(o(n^{1/2})) mu^n (Theorem 1.2, a corollary of DCH). Theorem 1.4/Corollary 1.5: **any quantitative DCH bridge estimate improves HW**: if phi(eps) = C eps^kappa (kappa > 1; Hutchcroft writes nu, not the Flory exponent) satisfies A2', then c_n <= exp(O(n^{(kappa-1)/(2kappa-1)})) mu^n; also B((1-eps) z_c) = O(eps^{(1-kappa)/kappa}) (Remark 1.6). In d = 3 this beats Kesten iff kappa < 3 (my arithmetic). | Z^d, d >= 2 | Hutchcroft, Electron. Commun. Probab. 23 (2018) 1-8, doi:10.1214/17-ECP94 | Read in full (arXiv:1708.09460v3) |
| C7 | DGHM 2020 (2D only): Z^2: for any eps < 1/466 there are infinitely many n with c_n <= exp(n^{1/2-eps}) mu^n; indeed infinitely many j for which all n in [j, j^{4/3}] satisfy it (Thm 1.2). H: for eps in (0, 1/42), c_n(H) <= exp(n^{1/2-eps}) mu(H)^n for all n large (Thm 1.3). H: power law c_n <= C_0 n^C mu^n **conditional** on Conjecture 1 (Thm 1.4). Nothing for d = 3. | Z^2, H | Duminil-Copin, Ganguly, Hammond, Manolescu, Ann. Probab. 48 (2020) 1644-1692, doi:10.1214/19-AOP1400 | Read: Section 1 and 2.1 (arXiv:1809.00760v2) |
| C7b | Madras, End patterns (J. Stat. Phys. 53 (1988) 689-701): for a pattern P that can start arbitrarily long SAWs, lim inf_N c_N(P)/c_N > 0, and the ratio of these fractions for N and N+2 converges to 1; convergence of c_N(P)/c_N itself is called "a difficult open problem" (abstract). | Z^d, d >= 2 | doi:10.1007/BF01014220 | Read: abstract + Section 1 (owner-supplied) |
| C8 | Open: lim c_{n+1}/c_n = mu for d = 2, 3, 4 (BDGS (1.18)); c_n mu^{-n} -> inf "appears to be open when d equals two or three" (Hammond EJP 2019 Sec. 1.4); existence of gamma open in d = 2, 3, 4; improving Kesten in d = 3, 4 open (BDGS after Thm 2.6). | | see `open-statements.md` | read |
| C9 | Polygons (p_n = n^{-theta_n} mu^n): d = 2 only: for any delta > 0 the set of even n with theta_n >= 3/2 - delta has limiting density one (Hammond AoP 46 (2018) 175-206, Thm 1.3); Madras (J. Stat. Phys. 78, 1995): theta_n >= 1/2 - o(1) in d = 2; Madras 1991 (Progr. Probab. 28): a bound on the generating function of p_n which in d = 3 "would assert lim theta_n >= 1 were this limit known to exist" (Hammond's description; Madras 1991 not read). d >= 2: Thm 1.5 only for 3-edge self-avoiding polygons (a variant, not strict SAW): theta-hat_n >= 2 - 1/d - delta on a density-one set. Classical: lim inf theta_n <= d + 5 (bridges, MS Thm 3.2.4 + Cor 3.1.8, as described by Hammond EJP 2019 Sec. 1.4); lim sup theta_n < inf is open. | | Hammond, "An upper bound on the number of self-avoiding polygons via joining", arXiv:1808.09032v1 | Read: abstract + Section 1 |

## D. Correlation length (very short; details in state-of-the-art.md)

- Chayes-Chayes, CMP 105 (1986) 221-238 (owner-supplied, read Sections 1-3): for z < z_c, point-to-plane function
  G_L(z) satisfies e^{-m(z)L} <= G_L(z) <= chi(z)^2 e^{-m(z)L} (Prop. 3.2), Ornstein-Zernike decay, m real analytic on
  z < z_c; m(z) -> 0 as z -> z_c (remark after Prop. 3.2, via their ref. [23]). **No rate for m(z) -> 0 is given.**
- Madras-Slade Ch. 4 (owner-supplied, read pp. 77-117 statements): m(z) concave in log z, strictly decreasing on
  (0, z_c), m(z) = bridge mass M(z) for z < z_c (Thm 4.1.14), lim_{z->z_c} m(z) = 0 (Cor 4.1.15), and the bridge mass at
  criticality is 0: M(z_c) = 0 (Thm 4.1.13) with lim_L -(1/L) log B_z(L) = M(z) for every z > 0 (Cor 4.4.3), so
  **e^{-o(L)} <= B_{z_c}(L) <= 1** for the critical point-to-plane bridge generating function; irreducible-bridge mass
  m_A(z_c) = 0 (Cor 4.4.4); Cor 4.1.17 = the Chayes-Chayes bound above. B_{z_c}(L) -> 0 is NOT stated (and not found
  elsewhere for Z^d).
- Open per Madras-Slade (p. 15, p. 18), no later proof found: m(z_c) = 0 in d = 2, 3, 4 ("a negative mass at the critical
  point has not yet been ruled out"); finiteness of the critical two-point function G_{z_c}(0, x) for any x != 0.
- No rigorous power-law bound on the correlation-length exponent for strict SAW in d = 2, 3, 4 found. Madras-Slade
  Chapters 3-4 (all theorem statements checked) contain none: p. 78 calls m(z) ~ const.(z_c - z)^{nu-bar} "believed (and
  known for d >= 5)" and "Proving that m(z) -> 0 as z -> z_c is a first step towards proving (4.1.2)"; nu-bar = nu is only
  argued heuristically (Sec. 1.3 for nu_2 under assumed power laws, Sec. 2.1 scaling theory, "no rigorous results in this
  section"). Answer to lower-bound/NEEDS.md: no "nu >= 1/d"-type inequality for m(z) or xi(z) is in these chapters.

## One-line summary for solvers

d = 3 today: N^{4/9}/36 <= E_N|omega_N|^2 = o(N^2) (no rate); P(close) <= n^{-1/2+o(1)}; sup_x P(Gamma_n = x) -> 0
(no rate); mu^n <= c_n <= mu^n exp(C n^{2/5} log n); c_{n+2}/c_n - mu^2 = O(n^{-1/3}). Quantitative sub-ballisticity
exists only on the hexagonal lattice (n/log n), via the parafermionic observable.
