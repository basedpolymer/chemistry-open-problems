# Rigorous state of the art for the SAW swelling exponent (strict SAW on Z^d, focus d = 3)

Literature agent, 2026-10-24. Companion files: `quick-facts.md` (exact statements, short), `open-statements.md`
(verbatim "open" sentences), `novelty.md` (2024-2026 sweep, AI campaigns), `references.md` (what was read).

Conventions. "Read in full" / "read: statements + proof overview" / "secondary" as in `quick-facts.md`. Everything is
for the nearest-neighbour **strict** SAW on Z^d unless a row says otherwise. Labels: CITED (source, theorem) /
OPEN (with the source that says so) / NOT FOUND (I did not find it in the sources read; this is not a proof that it is
open) / HEURISTIC.

## 1. Status table for the ladder of BRIEF-common.md

| Item | Statement | Status in d = 3 (as of 2026-10-24) | Source |
|---|---|---|---|
| L0 | E_N\|omega_N\|^2 >= c N^{4/(3d)} | **CITED, all N, c = 1/36**: E_N\|\|omega(N)\|\|_2^2 >= N^{4/9}/36 in d = 3 | Madras, CMB 57 (2014), Prop. 1.1 (read in full) |
| L1 | exponent > 4/9 for all large n | OPEN / NOT FOUND. Even exponent 2/3 = 2/d ("endpoint as far as the surface of a ball of volume n") is stated open. | Madras 2014 p. 114; BDGS lectures (1.30); DCH Question 4 |
| L1' | > 4/9 along a subsequence / Cesaro / generating-function sense | NOT FOUND. No subsequential or averaged improvement of 4/(3d) found. | sweep in `novelty.md` |
| L2 | n^{1-o(1)} | OPEN (implied by the openness of L1 at exponent 2/3 < 1) | as L1 |
| L3 | >= c n (diffusive) | OPEN in d = 2, 3, 4 | BDGS (1.30); Slade 2019 Sec. 2.1.3; DCH Question 3; Guttmann-Jensen 2022 p. 5 |
| L4 | >= n^{1+delta} | OPEN (stronger than L3) | idem |
| U0 | E_n\|omega_n\|^2 = o(n^2) | **CITED**, all d >= 2 | Duminil-Copin-Hammond, CMP 324 (2013), Thm 1.1 + Cor 1.2 (read in full) |
| U1 | explicit rate, e.g. n^2/(log n)^c | OPEN / NOT FOUND on Z^d for any d. Known only on the hexagonal lattice H: n^2/(log n)^2 follows from Krachun-Panagiotis Thm 1. | Slade 2019 ("there is no bound on the rate"); KP 2026 |
| U2 | n^{2-delta} | OPEN on Z^d (DCH Question 1; BDGS (1.30); Slade 2019: "it is not proved that ... R_n <= C n^{0.99999}"). On H: claimed only in an unrefereed AI catalogue (OpenAI entry 237, diameter n^{3/4+o(1)}), not verified. | see `open-statements.md`, `novelty.md` |
| E1 | lim log E_n\|omega_n\|^2 / log n exists | OPEN in d = 2, 3, 4 ("existence of nu"; Slade 2019: "For d = 2, 3, 4, none of the entries in Table 1 have been proved") | Slade 2019 Sec. 2.1.3; Madras 2014 p. 114 |
| (T) | E_n\|omega_n\|^2 = n^{2 nu + o(1)}, determine nu | OPEN. BDGS 1.6.3: "For d = 3, there are no rigorous results for critical exponents, and no mathematically well-defined candidate has been proposed for the scaling limit"; MS p. 8: "There is no analogous exact solution in three dimensions, and the d = 3 values [...] are from numerical results and field-theoretic calculations"; DCH p. 2: "no significant physical predictions have been made regarding the scaling limit of self-avoiding walk in three dimensions"; Guttmann-Jensen p. 4 (for gamma): "In three dimensions there is no conjectured exact value" | BDGS; MS; DCH; GJ22 |

Related quantities:

| Quantity | Best rigorous statement for Z^3 found | Source |
|---|---|---|
| Endpoint at a neighbour of 0 (closing probability) | <= n^{-1/2+eps}, all large odd n, all d >= 2 | Hammond, Trans. AMS 372 (2019), Thm 1.2 |
| Endpoint at a general point | sup_x P_n(Gamma_n = x) -> 0, no rate, all d >= 2 | DGHM, AoP 44 (2016), Thm 1.2 |
| Midpoint at a point | sup_x P_n(Gamma_{floor(n/2)} = x) <= C n^{-1/2}, all n | DGHM 2016, Prop. 1.3 |
| c_n / mu^n | 1 <= c_n mu^{-n} <= exp(C n^{2/5} log n) in d = 3 | Kesten, JMP 5 (1964), Thm 2 (proof only sketched there) |
| c_n / mu^n, d = 2 | exp(o(n^{1/2})) all n (Hutchcroft 2018); exp(n^{1/2-eps}), eps < 1/466, infinitely many n (DGHM 2020) | |
| c_n mu^{-n} -> inf | OPEN in d = 2, 3 ("appears to be open") | Hammond EJP 2019 Sec. 1.4 |
| Ratio limits | c_{n+2}/c_n -> mu^2 with \|c_{n+2}/c_n - mu^2\| <= A n^{-1/3}; lim c_{n+1}/c_n = mu OPEN in d = 2, 3, 4 | Kesten, JMP 4 (1963), Thm 4; BDGS (1.17)-(1.18) |
| Polygons p_n (up to translation), d = 3 | Kesten 1964 (1.9): exp(-alpha_7 n^{2/(d+1)} log n) mu^{2n-1} <= gamma_{2n-1} <= 4(d-1) n mu^{2n-1} (gamma = closing walks); q_m q_n <= (d-1) q_{m+n} (BDGS (2.26)); lim inf theta_n <= d+5; lim sup theta_n < inf OPEN. Madras 1991 bound "would assert lim theta_n >= 1 [in d = 3] were this limit known to exist" (Hammond's description; Madras 1991 not read). | Kesten 1964; BDGS Sec. 2.2; Hammond AoP 2018 Sec. 1.2, EJP 2019 Sec. 1.4 |
| Bridges | b_n <= mu^n (supermultiplicativity); b_n >= mu^{n-1} e^{-B sqrt(n-1)} (BDGS Cor 2.7); bridge generating function B(1/mu) = inf (BDGS Cor 2.8, via chi(z) <= z^{-1} exp(2(B(z)-1))); bridges of span >= vn are exponentially rare among bridges of length n (DCH, in Hutchcroft's Thm 1.3 form) | BDGS Sec. 2.1; DCH 2013; Hutchcroft 2018 |
| Critical bridge partition function per height, B_T(z_c) | On H: B_T <= 100 T^{-1e-10} (KP Thm 2), B_T >= c/T (Duminil-Copin-Smirnov). **On Z^d: only e^{-o(T)} <= B_T(z_c) <= 1** (bridge mass M(z_c) = 0: MS Thm 4.1.13, Cor 4.4.3; upper bound from supermultiplicativity in the span, MS (4.1.12)). B_T(z_c) -> 0 NOT FOUND for Z^d (by Kesten's renewal structure this should be equivalent to infinite mean span of an irreducible bridge under P_iSAB -- my reading, to be checked). | KP 2026 Sec. 1; MS Ch. 4 |
| Irreducible bridges | Kesten: sum over irreducible bridges of mu^{-\|gamma\|} = 1 (DCH Lemma 2.2); E_iSAB(\|gamma\|) = inf (DCH Thm 2.5); conjecturally stable law (DCH Remark 4.3) | DCH 2013 |
| Correlation length / mass m(z), z < z_c | Ornstein-Zernike decay for all z < z_c; e^{-m(z)L} <= G_L(z) <= chi(z)^2 e^{-m(z)L} for the point-to-plane function (Chayes-Chayes Prop. 3.2 = MS Cor 4.1.17); m real analytic on z < z_c (MS Thm 4.2.8); m(z) -> 0 as z -> z_c (MS Cor 4.1.15). Open per MS p. 15 and p. 18 (no later proof found): m(z_c) = 0 and even finiteness of G_{z_c}(0, x) in d = 2, 3, 4. **No rate and no bound on the correlation-length exponent found** for d = 2, 3, 4. Elementary consequence I derived (check before use): bridges are supermultiplicative in height, so B_T(z) <= e^{-m_B(z) T} with m_B the bridge mass, and BDGS (2.20) with chi(z) >= z_c/(z_c - z) gives sum_T B_T(z) = B(z) >= (1/2) log(z_c/(z_c - z)) - O(1), hence m_B(z) <= C / log(1/(z_c - z)): only a logarithmic divergence of the correlation length is obtained this way. | Chayes-Chayes, CMP 105 (1986) (read Secs. 1-3); BDGS (1.35)-(1.38) |
| Existence of gamma, eta, nu (correlation length) | OPEN in d = 2, 3, 4 | Slade 2019; BDGS Sec. 1.5-1.6 |
| d >= 5 (context only) | nu = 1/2, gamma = 1, Brownian scaling limit (Hara-Slade 1992) | Slade 2019 Thm 2.1 (secondary) |

## 2. Key papers: method, and what blocks (a) Z^3 and (b) a quantitative rate

### 2.1 Madras 2014 (CMB 57, 113-118) -- L0. Read in full.

*Statement.* Prop. 1.1: for every N >= 1 and p in [1, inf), (E_N\|\|omega(N)\|\|_inf^p)^{1/p} >= (1/6) N^{p/((p+1)d)}.
*Method.* Let 9omega9 = max_k \|\|omega(k)\|\|_inf, L the first time this maximum is reached, I the first coordinate where it is
attained, H = omega_I(L). The map T reflects omega after time L through the hyperplane x_I = H (a face of the minimal
cube containing omega). Lemma 2.2 (volume): N + 1 distinct sites cannot fit in the cube of radius N^{1/d}/3, so
9omega9 > N^{1/d}/3. Lemma 2.3: T maps {\|\|omega(N)\|\|_inf <= m} into {\|\|omega(N)\|\|_inf >= N^{1/d}/3} and is at most
(m+1)-to-one (from the image one recovers I, and H up to m+1 values), so
P(\|\|omega(N)\|\|_inf <= m) <= (m+1) P(\|\|omega(N)\|\|_inf >= N^{1/d}/3). Then Markov's inequality twice.
*Why 4/(3d) and what blocks improvement.* The (m+1) multiplicity loss is the whole weakness; Madras' acknowledgment:
"I was hoping that the argument could be improved, for example, by showing that the restriction of T to S_{N,<=m} is
typically much better than (m+1)-to-one. However, no such improvement was apparent." The only geometric input is the
deterministic volume bound N^{1/d}; no information on how far *beyond* N^{1/d} the walk goes is used.
*(a) Z^3:* no obstruction (valid for all d >= 2). *(b) rate:* already explicit.

### 2.2 Duminil-Copin and Hammond 2013 (CMP 324, 401-423) -- U0. Read in full (arXiv v1).

*Statement.* Thm 1.1: for each v > 0 there is eps > 0 with P_n(max_k \|\|gamma_k\|\| >= vn) <= e^{-eps n} for every n;
Cor 1.2: E_n\|\|gamma_n\|\|^2 = o(n^2). Valid on Z^d for every d >= 2.
*Method (proof by contradiction).* Reduce to bridges by Hammersley-Welsh unfolding (cost e^{C sqrt n}). Kesten's
renewal decomposition into irreducible bridges gives a probability measure P_iSAB (weights mu^{-\|gamma\|}, total mass 1,
Lemma 2.2) and the infinite i.i.d. concatenation P_iSAB^{tensor N}. Two incompatible facts:
(i) Thm 2.3 / Cor 2.4: under the *ballistic assumption* (bridges reach height vn with subexponential probability), at
subexponential cost and **along a subsequence**, bridges have a positive density of renewal points (unfolding of
"zigzags", multi-valued map principle, three-case analysis), which forces E_iSAB\|gamma\| < inf.
(ii) Thm 2.5: E_iSAB\|gamma\| = inf. If the mean were finite, the law of large numbers and ergodicity of the shift on
P_iSAB^{tensor Z} give a positive density of *diamond points* (Prop. 4.2); the "stickbreak" surgery (rotate the piece
between two diamond points by a right angle, insert two east-west edges) produces order n^2 alternative bridges of
large width from each typical bridge, with at most (2 nu n)^2 preimages, contradicting the LLN (width o(n)).
*(a) Z^3:* none, the proof is written for Z^d, d >= 2 (rotation in a coordinate 2-plane).
*(b) rate:* every step is qualitative: proof by contradiction along subsequences with o(n) exponents (Thm 3.1,
Prop. 3.2), law of large numbers and ergodic theorem in Thm 2.5 and Prop. 4.2, no control of eps(v). KP 2026 say: "their
result relies on the ergodicity of the latter measure, and for this reason it does not provide any quantitative bound on
the displacement of self-avoiding walk. A major obstacle in that direction seems to be the lack of a quantitative bound
on the partition function of bridges of a given height at the critical point x_c = 1/mu that goes to 0 as the height
goes to infinity." Hutchcroft 2018: "The proof of Theorem 1.3 is not quantitative and does not give any estimates on the
function phi."

### 2.3 Krachun and Panagiotis (AoP 54 (2026) 1109-1125) -- quantitative U1 on the hexagonal lattice. Read in full (arXiv v1).

*Statements.* Thm 1 (H): P_n(max_k \|\|gamma_k\|\| >= 10^{11} n/log(n+1)) <= exp(-c n^{2/3}) for every n >= 1.
Thm 2 (H): B_T <= 100 T^{-eps}, eps = 10^{-10}.
*Method.* (1) Parafermionic observable at x_c = (2+sqrt 2)^{-1/2}, sigma = 5/8 (Duminil-Copin-Smirnov): contour identities
in strips and triangles give cos(3pi/8) A_T + B_T = 1, monotonicity of B_T and of the triangle-walk function D_T,
B_k <= cos(pi/8) D_k, and the summability sum_{k>=1} G_k <= 1/(2 cos(3pi/8)) of the half-plane boundary two-point
function (Lemma 2.2). (2) Glazman-Manolescu's concatenation of three triangle walks into a half-plane U-walk, refined by
restricting to walks with not too many renewal times (Markov's inequality on E[N], Lemma 3.1) to control preimage
multiplicity; this gives the recurrence T^4 D_{18T}^5 <= 2^{17} (sum_{i<=3T} D_i)^4 sum_{k=T}^{21T} G_k (Cor. 3.1),
and an analytic lemma (Lemma 3.4) turns it into polynomial decay. (3) Quantitative sub-ballisticity (Section 4,
Lemma 4.1): a bridge of length n and height >= w = 5 eps^{-1} n/log(n+1) is cut along m = n^{1/3} interleaved families of
horizontal lines; for a good family (at most (n+1)/m line visits) the bridge is a concatenation of s = w/m bridges of
height m (weight <= B_m^s <= (100 m^{-eps})^s = exp(-(5/3+o(1)) n/m)) and of pieces between first and last visits of a
line, each a concatenation of U-walks (weight factor <= C + 1 per visit, C = 2 sum G_k < 3). Then Hammersley-Welsh
unfolding from walks to bridges.
*(a) Z^3:* both analytic inputs come from the parafermionic observable, which exists only on H ("The latter is only
known in the case of the hexagonal lattice among low dimensional lattices and its proof uses the parafermionic
observable"): (i) polynomial decay of B_T at z_c and (ii) finiteness of the critical half-plane arch sum. On Z^3
neither is known (NOT FOUND in sources read).
*(b) rate:* the rate n/log n comes from B_T <= T^{-eps} with tiny eps; a better bridge exponent would improve the log
only if combined with a different decomposition.

### 2.4 Duminil-Copin, Glazman, Hammond, Manolescu (AoP 44 (2016) 955-983) -- endpoint. Read: Secs. 1-3, overviews of 4-5 (arXiv v2).

*Statements.* Thm 1.1 (d >= 2): P_n(\|\|Gamma_n\|\| = 1) <= n^{-1/4+eps} for n large. Thm 1.2 (d >= 2): sup_x P_n(Gamma_n = x)
-> 0. Prop. 1.3: sup_x P_n(Gamma_{floor(n/2)} = x) <= C n^{-1/2}.
*Methods.* Kesten's pattern theorem recast with *shells* (walks modulo swapping type-I/type-II patterns, which differ in
length by 2): Lemma 3.3 (shell probabilities stable when the length changes by 2k, k <= n^{1/5}, using
Hammersley-Welsh), Lemma 3.5 (conditionally on a shell, type-I patterns are allocated uniformly: Gaussian fluctuations
of order sqrt n, giving Prop. 1.3). Thm 1.2: first for bridges, via z-renewal times (a renewal followed by a transverse
edge that can be rotated in 2d-2 ways, Prop. 4.3), the number of which tends to infinity because the bridge
generating function diverges at z_c (Prop. 4.2, using MS Cor 3.1.8); then half-space walks, then walks (Lemmas 4.5-4.6).
Thm 1.1: by contradiction, with polygonal invariance (the law of a closing walk is invariant under cyclic shift, the
hanging point depends only on the polygon), pattern resampling (Lemma 3.5) to move the hanging time by up to sqrt n, and
a nested family of "good" first parts which forces the remaining walk to avoid them with stretched-exponentially small
probability.
*(a) Z^3:* valid for all d >= 2. The authors note that extending Thm 1.1 to a general fixed x is open: "For higher
dimension, significant difficulties occur and it is unclear how to prove this result."
*(b) rate:* Thm 1.2 has no rate because Prop. 4.2 only uses B(z_c) = inf (no divergence speed) and Prop. 4.3 is "routine"
but qualitative; Thm 1.1's exponent 1/4 is limited by the sqrt n pattern fluctuation; Hammond (TAMS 2019) reworks the
method to n^{-1/2+o(1)} and says that without significant modification the DGHM method "cannot prove an upper bound on
this quantity that decays more rapidly than n^{-1/2+o(1)}".

### 2.5 Hammond's snake method: TAMS 372 (2019) 2335-2356 (pattern fluctuation) and EJP 24 (2019) 1-43 (polygon joining). Read: intros, statements, Sec. 4.1 outline.

*Statements.* TAMS Thm 1.2 (d >= 2): closing probability <= n^{-1/2+eps} for all large odd n. EJP Thm 1.5 (d = 2):
<= n^{-4/7+eps} on a set of odd n of upper density >= 1/1250; conditional on existence of theta, xi: theta + xi >= 5/3.
*Method.* A general sufficient condition ("snake method"): assuming a closing probability lower bound n^{-beta},
build a "snake" -- a nested sequence of first parts gamma[0, l] each of which has high conditional closing probability --
and show that a "charming" snake (many terms with high closing probability) forces the second part to avoid many
reflections, contradicting the assumption. The snake is produced by "square-root fuzziness": Gaussian fluctuation of
pattern counts moves the index l ~ n/2 over a window of order n^{1/2} (TAMS), or by Madras' polygon joining (EJP, d = 2).
*(a) Z^3:* the pattern-fluctuation version holds in all d >= 2; the polygon-joining version uses planarity (Madras'
joining).
*(b) limits:* the n^{1/2} window caps the method at n^{-1/2}. Relation to nu: only through the conjectural hyperscaling
theta = 1 + d nu (Hammond TAMS Sec. 1.4, HEURISTIC).

### 2.6 Duminil-Copin, Ganguly, Hammond, Manolescu (AoP 48 (2020) 1644-1692) -- Hammersley-Welsh with polygon insertion. Read: Secs. 1, 2.1-2.3 outline, 5.5.

*Statements.* Z^2: c_n <= exp(n^{1/2-eps}) mu^n for infinitely many n (eps < 1/466); H: for all large n (eps < 1/42); H:
c_n <= C_0 n^C mu^n conditional on Conjecture 1 (connective constant of the universal cover of H minus a ball).
*Method.* Restrict the Hammersley-Welsh map Psi (half-space walk -> decreasing list of bridges) to half-space walks with
many branches (>= n^{1/2-eps}) confined to a strip of width n^{1/2+eps}; such walks have many "near self-touches";
insert self-avoiding polygons of a given length at near self-touches to make Psi a multi-valued map with many outgoing
arrows, while keeping incoming arrows controlled. Needs *polygon abundance* (power-law lower bound on suitably wide
polygons; subsequential on Z^2, all n on H via the observable) and Madras' polygon joining.
*(a) Z^3:* "Lemma 5.5 depends on planarity"; and in d >= 3 the Hammersley-Welsh exponent 1/2 is not the relevant
benchmark since Kesten already has n^{2/(d+2)} log n. Not adapted to d = 3.
*(b)* explicit eps but tiny; Sec. 5.5 says the method seems unlikely to go beyond eps close to 1/12.

### 2.7 Hutchcroft 2018 (ECP 23, 1-8) -- quantitative bridge: DCH => HW. Read in full.

*Statements.* Thm 1.2 (d >= 2): c_n <= exp(o(n^{1/2})) mu^n. Thm 1.4: an explicit Psi built from any increasing phi in
the bridge estimate #{n-step bridges reaching height >= m} <= b_n exp(-phi(m/n) n) gives c_n <= exp(Psi(n) - 2) mu^{n+1}.
Cor 1.5: phi(eps) = C eps^nu (nu > 1) gives c_n <= exp(O(n^{(nu-1)/(2nu-1)})) mu^n. Remark 1.6: B((1-eps) z_c) <=
[1 - (1-eps)^{psi(eps)}]^{-1}, = O(eps^{(1-nu)/nu}) in the power case.
*Method.* A new proof of Hammersley-Welsh without unfolding combinatorics, through the bridge generating function.
*Relevance:* any quantitative DCH (U1/U2 in bridge form) immediately improves c_n bounds in all d, d = 3 included
(compare with Kesten's n^{2/5} log n).

### 2.8 Kesten 1963, 1964 (JMP 4, 960-969; JMP 5, 1128-1137). Owner-supplied scans; statements and Sec. 4 of 1964 read.

*Statements.* 1963: pattern theorem (Thm 1); |c_{n+2}/c_n - mu^2| <= A_3 n^{-1/3} (Thm 4); half-space walks nu_n satisfy
nu_{n+1}/nu_n -> mu (1.8); polygon lower bound e^{-alpha_4 n^{1/2}} mu^n <= gamma_n (1.3). 1964: Thm 2 (1.8)-(1.9)
(c_n <= exp(alpha_6 n^{2/(d+2)} log n) mu^n); Cor 2 (ratio rates in d dimensions, stated without proof).
*Method of the d >= 3 bound (Sec. 4, sketched; only the weaker version with n^{2/(d+1)} is sketched, "(1.8) which is a
refinement ... we do not want to give"):* pigeonhole a common endpoint Z for (2n+1)^{-d} c_n walks, glue two such walks
(one reversed, after a translation Y chosen by averaging so that the two copies meet at most n^beta times), close up
into a closed path v of length <= 2n + 2d n^{(2-beta)/d} - 1 with few double points, then count such paths by attaching
polygons; optimise beta. The d-dependence comes only from the averaging over translates (volume n^{2-beta} spread
over a box of side n^{(2-beta)/d}).
*Proof status:* Kesten writes "it does not seem worthwhile to give more than an indication of the proof of Theorem 2";
the full proof of c_N <= mu^N exp[Q N^{2/(d+2)} log N] (all d >= 2, every N >= 2) is Madras-Slade Theorem 3.3.1
(pp. 68-74, "due to Kesten (private communication)"), built on Lemma 3.3.3 and the half-space decomposition (3.3.4).
MS p. 69: "The full proof was not given in Kesten (1964) because it was hoped that someone would find a better bound,
but almost thirty years later this has still not come to pass for d = 3, 4." Still true in 2026 as far as I found.

### 2.9 Chayes and Chayes 1986 (CMP 105, 221-238) -- correlation length. Owner-supplied, Secs. 1-3 read.

Ornstein-Zernike decay at every z < z_c (beta > beta_c), point-to-plane function with pure exponential decay
(2.5), Gaussian hitting distribution (2.6), mass real analytic; divergence of the correlation length m(z)^{-1} as
z -> z_c is quoted from their ref. [23] via a function P(m) with G_L <= P(m) e^{-mL}. Nothing near criticality
quantitatively. *(a)* all d. *(b)* no rate.

### 2.10 Other recent results checked (not on the ladder)

- Lehner, Lindorfer, Panagiotis, "Self-avoiding walk is ballistic on graphs with more than one end" (Forum Math. Sigma,
  arXiv:2403.13121): ballistic on graphs with > 1 end; irrelevant to Z^3 (one end) except as contrast. Abstract only.
- Grimmett, "On counting polygons in a crystal" (arXiv:2601.01128): polygon growth rate = connective constant on
  quasi-transitive graphs with a square graph height function. Abstract only. Not about nu.
- Liu-Slade (arXiv:2605.15545), Duminil-Copin-Markar-Panis-Slade (arXiv:2605.21438), Duminil-Copin-Panis
  (arXiv:2410.03649): near-critical / mean-field results above the upper critical dimension only (d >= 5 for SAW).
  Abstracts only.
- Lazar (arXiv:2504.13093): formal lattice-point/Fourier "asymptotic" formulas for c_n; no exponent bound; not
  relevant (read: abstract, statements).

## 3. Leads for the solvers (each with the source that suggests it; none is a result)

1. **Better multiplicity for Madras' reflection** (source: Madras 2014, Acknowledgments). Show that T restricted to
   {\|\|omega(N)\|\|_inf <= m} is typically much better than (m+1)-to-one. Madras tried and found "no such improvement".
   Any gain in Lemma 2.3 of the form P(\|\|omega_N\|\|_inf <= m) <= m^{1-a} P(\|\|omega_N\|\|_inf >= N^{1/d}/3) raises the
   exponent. The computation lab can measure the actual multiplicity of T on exact enumerations / samples.
2. **Anti-concentration => displacement** (elementary, my observation; check before use). If sup_x P_n(Gamma_n = x) <= p_n
   then P(\|\|Gamma_n\|\|_inf <= r) <= (2r+1)^d p_n, hence E_n\|Gamma_n\|^2 >= c p_n^{-2/d}. So a *quantitative* DGHM
   Thm 1.2 with p_n = n^{-a} gives exponent 2a/d: a > 2/3 beats Madras (4/(3d)) in every d; L3 in d = 3 needs a >= 3/2
   (local-CLT scale). DGHM's midpoint bound C n^{-1/2} (Prop. 1.3) gives only E\|Gamma_{n/2}\|^2 >= c n^{1/d} for the
   midpoint. Source of the mechanism: DGHM 2016 Prop. 1.3 / Lemma 3.5 (pattern fluctuation of order sqrt n).
3. **Logarithmic divergence of the bridge generating function** (elementary consequence of BDGS (2.20) and
   chi(z) >= z_c/(z_c - z); check before use): B(z) >= (1/2) log(z_c/(z_c - z)) - O(1) as z -> z_c; since b_j z_c^j <= 1,
   this yields partial sums sum_{j<=k} b_j z_c^j >= (1/2) log k - O(log log k). Plugged into DGHM Prop. 4.2 this may give
   a (logarithmic) rate for "bridges have many z-renewal times" and hence a log-rate version of Thm 1.2 for bridges.
   Too weak for L1 by itself, but it is a quantitative handle on the one place where DGHM use divergence at z_c.
4. **KP's route on Z^3** (source: KP 2026 Sec. 1 and Lemma 4.1). Their Section 4 deduction is lattice-agnostic given two
   inputs: (i) B_T(z_c) <= C T^{-eps} for the critical bridge partition function of height T; (ii) finiteness of the
   critical generating function of half-space arches (walks from the boundary hyperplane back to it, KP's sum G_k).
   On Z^3 both are NOT FOUND. Heuristically (ii) should hold (arch count ~ mu^n n^{gamma_11 - 1} with surface exponent
   gamma_11 < 0; HEURISTIC), and Duplantier-Guttmann (arXiv:1908.03872, physics, abstract only) give the scaling
   relation gamma_b = gamma_{1,1} + nu for the bridge exponent, "with compelling numerical evidence" in d = 2, 3; combined
   with b_n ~ mu^n n^{gamma_b - 1} this suggests B_T(z_c) ~ T^{-a} with a = -gamma_{1,1}/nu (= 1/4 in d = 2, matching
   KP's prediction B_T ~ T^{-1/4}). HEURISTIC only. The upper-bound agent's `upper-bound/result.md` develops this route
   (conditional theorem under sum_T B_T(z_c)^2 < inf). Even B_T(z_c) -> 0 on Z^d was not found; a qualitative proof of B_T(z_c) -> 0 on Z^3 would
   itself be a new, checkable sub-goal (equivalently, by Kesten's renewal theory, infinite mean height of irreducible
   bridges -- my reading). KP explicitly call (i) "a major obstacle".
5. **Quantitative DCH in bridge form gives quantitative Hammersley-Welsh** (source: Hutchcroft 2018 Thm 1.4, Cor 1.5,
   Remark 1.6): an explicit phi in "#bridges reaching height >= m <= b_n exp(-phi(m/n) n)" is the natural U1 target.
   Writing Hutchcroft's exponent as kappa (his "nu", not the Flory exponent): phi(eps) = C eps^kappa (kappa > 1) gives
   c_n <= exp(O(n^{(kappa-1)/(2kappa-1)})) mu^n, which beats Kesten's n^{2/5} log n in d = 3 when
   (kappa-1)/(2kappa-1) < 2/5, i.e. kappa < 3 (my arithmetic). Cross-check for solvers: a quantitative sub-ballisticity
   claim in this bridge form with kappa < 3 would also improve Kesten's 1964 bound, a 60-year-old record.
   HEURISTIC: Fisher's tail law P(\|omega_n\| ~ r n^nu) ~ exp(-r^{1/(1-nu)}) suggests kappa = 1/(1-nu) ~ 2.42 in d = 3,
   so the expected truth does lie in the range kappa < 3.
6. **Where DCH loses quantitativity** (source: DCH 2013 Secs. 3-4): Thm 2.5 (E_iSAB\|gamma\| = inf) uses LLN/ergodicity;
   a quantitative replacement would be a lower bound on the tail P_iSAB(\|gamma\| > k) (DCH Remark 4.3 conjectures a
   stable law). Thm 2.3 works along subsequences with subexponential costs. Both must be made quantitative.
7. **Hyperscaling theta = 1 + d nu** (source: Hammond TAMS Sec. 1.4, HEURISTIC): polygon-count exponents and nu are
   conjecturally linked; no rigorous inequality in either direction was found for d = 3. Not a proof route as such.
8. **AI-claimed honeycomb results** (OpenAI catalogue entry 237, unrefereed): a "geometric passage" from strip/bridge
   inputs to diameter exponents on H. All analytic inputs are integrable (parafermionic observable, Yang-Baxter cylinder
   transfer). If someone wants to mine it, the part to look at is Secs. 3-9 of "Mass and covering exponents for
   fixed-length honeycomb walks", which claims to need only (a) B_h ~ h^{-1/4}, (b) strip first-length mass
   <= h^{13/12+o(1)}, (c) polygon length-square mass, (d) a one-arc bound. None of these is available on Z^3.
