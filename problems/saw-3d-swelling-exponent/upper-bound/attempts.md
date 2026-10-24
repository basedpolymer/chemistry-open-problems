# Upper-bound solver: attempts, obstructions, and where each route stops

Solver: Claude Opus 5.5 (effort not recorded), 2026-10-24. Companion of `result.md` (which holds the only
claims marked PROVED). Labels: PROVED / CITED (source, number) / HEURISTIC / OBSTRUCTION (a precise reason why a
specific proof architecture cannot give a rate; not a theorem that no other proof exists).

Notation as in `result.md`: height = first coordinate, B_T(z) = critical-bridge partition function at height T,
z_c = 1/mu, P_iSAB(gamma) = mu^{-|gamma|} on irreducible bridges (DCH Lemma 2.2; Madras–Slade (4.2.4)),
L = |gamma| and H = height of an irreducible bridge under P_iSAB.

## A. Making Duminil-Copin–Hammond (DCH 2013) quantitative

**A1. Architecture (CITED, DCH arXiv:1205.0401v1, read in full).** Contradiction between
(i) Theorem 2.3 / Theorem 3.1 + Proposition 3.2: if P_SAB_n(y(gamma_n) >= vn) >= e^{-o(n)} (ballistic assumption),
then P_SAB_n(|R_gamma| >= delta n) >= e^{-o(n)} for some delta > 0 (renewal points at positive density at
subexponential cost; proved by unfolding zigzags, along subsequences, with epsilon_n -> 0);
(ii) Corollary 2.4: then E_iSAB|gamma| < infinity (large deviations of a renewal process);
(iii) Theorem 2.5: E_iSAB|gamma| = infinity (law of large numbers + ergodicity of P_iSAB^{tensor Z}, diamond points,
"stickbreak" surgery).
Then Section 2.4 transfers bridges to walks with the Hammersley–Welsh loss e^{C sqrt n}.

**A2. What a rate needs (PROVED-level bookkeeping).** Write the bridge large deviation as
#{n-step bridges of height >= vn} <= b_n e^{-phi(v) n} (Hutchcroft 2018, Theorem 1.3 form). With the HW loss
e^{C sqrt n} (Proposition 2 of result.md has e^{O(sqrt n log n)}), any rate requires phi(v_n) n >> sqrt(n) log n,
i.e. an explicit phi with phi(v) >> n^{-1/2} at v = v_n -> 0. Hutchcroft (Thm 1.4, Cor. 1.5, Remark 1.6) shows
that any explicit phi also improves the loss itself (B((1-e) z_c) <= [1 - (1-e)^{psi(e)}]^{-1}); with
phi(v) = C v^nu (nu > 1) one gets E_n|omega_n|^2 <= n^{2 - 2/(2nu-1) + o(1)}. Hutchcroft Remark 2.5: necessarily
phi(v) <= (1+o(1)) * 2v/log(1/v) (so nu <= 1 is impossible).
(Remark, d = 3: Kesten's bound h_n <= mu^n exp(Q n^{2/5} log n) (Madras–Slade Thm 3.3.1, (3.3.3)) can replace
HW in the reduction, via the first-bridge decomposition (3.3.4): a half-space walk of span >= h is a bridge of
span >= h followed by a half-space walk, so #{span >= h} <= n^2 max_i b^{>=h}_i h_{n-i}. The loss becomes
e^{O(n^{2/5} log n)}, at the price of needing the bridge bound at all lengths i <= n. Not needed in result.md.)

**A3. Theorem 2.5 quantitatively = Cesaro decay of b_n/mu^n.** Under P_iSAB^{tensor N}, b_n mu^{-n} = u_n =
P(n is a renewal time) (DCH Lemma 3.9; Madras–Slade (4.2.6)). Wald's identity applied to L wedge (n+1) (the renewal
times in [0,n] are unchanged by this truncation) gives, PROVED (elementary):
  (n+1)/E[L wedge (n+1)] <= sum_{j=0}^n u_j <= (2n+1)/E[L wedge (n+1)].
So a quantitative Theorem 2.5 ("E[L wedge k] >= f(k) -> infinity explicitly") is *equivalent* to an explicit
Cesaro bound (1/n) sum_{j<=n} b_j mu^{-j} <= C/f(n). No such bound is known on Z^3 (Madras–Slade 1993, p. 92:
"it is believed that lim_N b_N/mu^N = 0, but there is no known proof of this"; DCH proved it without rate).
Quantifying the stickbreak proof is possible in principle but extremely lossy: the LLN is used with no moment
(one needs a pigeonhole over scales k < k^A where E[L wedge k] is nearly flat, costing a tower/iterated-log), and
the diamond-point density is only bounded by delta >= mu^{-2K} rho_K (DCH (4.2)) with K the scale at which the
infinite bridge stays in a cone with probability rho_K >= 1/2: exponential in K. One would obtain at best
f(k) of iterated-logarithmic type. I did not carry this out, because of A4.

**A4. OBSTRUCTION in Section 3 (the decisive one).** Quantitatively, Proposition 3.2 would have to show:
P(height >= vn) >= e^{-eps n} implies P(|R| >= delta(v) n) >= e^{-g(eps, v) n} with g -> 0, and then contradict
the renewal large deviation P_SAB_n(|R| >= rho n) <= exp(-I_R(rho) n + C sqrt n). Two facts:
* (PROVED, elementary) I_R(rho) <= rho log mu + o(1): concatenating rho n one-step bridges with a bridge of length
  (1 - rho) n gives P_SAB_n(|R| >= rho n) >= b_{(1-rho)n} mu^{-rho n}/b_n >= mu^{-rho n} e^{-C sqrt n}
  (b_n <= mu^n, b_k >= mu^k e^{-C sqrt k}). Moreover I_R(rho) = rho Lambda*(1/rho) (Cramer for the i.i.d. lengths)
  and Lambda*(x) -> 0 as x -> infinity, so I_R(rho) = o(rho) (heuristically ~ rho^{1/(1-gamma_b)} ~ rho^{1.25} in
  d = 3).
* (Bookkeeping of DCH's proof) The three cases of Prop. 3.2 with explicit parameters: zigzag density threshold eta,
  short zigzags of length <= 1/delta' with delta' = eta/2, number of unfolded zigzags delta'' n with
  delta'' < delta'^3/12, renewal threshold rho < delta'^2/4. Case 2 (few zigzags) costs
  exp(2 eta n log(e/(2 eta))) per iteration and is iterated ceil(2/v) times; Cases 1 and 3 land at renewal
  density min(rho, delta'') = O(eta^3). So one needs (2/v) eta log(1/eta) < I_R(c eta^3) <= c eta^3 log mu: false
  for every small eta. More robustly: every DCH-type step reconstructs the preimage by choosing a subset of the
  renewal points of the output, which costs up to |R_phi| log 2 = Theta(rho) n when the chosen subset is a
  positive fraction, while the gain is I_R(rho) n = o(rho) n. **The unfold-to-create-renewals strategy cannot
  give an exponential-scale rate.** In DCH this is invisible because all costs are e^{o(n)} along a subsequence.

Conclusion A: a rate cannot come from quantifying DCH as written; a mechanism whose entropy cost is per *visit*
rather than per *created renewal* is needed. Krachun–Panagiotis provide one (Route B).

## B. The Krachun–Panagiotis route in Z^3

**B1. KP mechanism (CITED, arXiv:2310.17299v1, Lemma 4.1 and proof of Theorem 1, read in full).** Cut a bridge of
height >= w at the last visit of lines l_{i,j} = {y = jm + i} and the first visit of the next line; pick the
residue i for which the total number of visits to these lines is <= (n+1)/m. The pieces are bridges of height m
(weight B_m each) and returns to a line (weight <= (C+1) per visit, C = 2 sum_k G_k, G_k = half-plane walks from 0
to (k,0)). Result: weight <= n^2 B_m^{w/m} (C+1)^{(n+1)/m}. Inputs: Theorem 2 (B_T <= 100 T^{-1e-10}) and
Lemma 2.2 (sum_k G_k <= 1/(2 cos(3pi/8))), **both from the parafermionic observable** (planar, hexagonal-specific;
identity cos(3pi/8) A_T + B_T = 1 of [GM20, Cor. 2.2]). KP p. 3: the implication "relies on the summability of the
half-plane boundary two-point function ... only known in the case of the hexagonal lattice among low dimensional
lattices".

**B2. Transplanting to Z^3 (PROVED, result.md Theorem A).** Everything except the two inputs is dimension-free
(Proposition 3). New elementary input, result.md Lemma 1: the half-space arch generating function satisfies
A(z) <= sum_M B_M(z)^2 (split the arch at its first highest point). So in Z^d the second input follows from the
first as soon as sum_T B_T(z_c)^2 < infinity. This gives the conditional U1 (Corollary B): (H) implies
E_n|omega_n|^2 = O(n^2/(log n)^2); Corollary B2: the borderline B_T <= C T^{-1/2} gives O(n^2 (log log n)^2/(log n)^2).

**B3. Why d = 3 and not d = 2 (HEURISTIC).** Scaling gives B_T(z_c) ~ T^{-a} with a = 1 - gamma_b/nu
(result.md Sec. 6.5): a = 1/4 in d = 2 (matches KP's "B_T ~ T^{-1/4}"), a ~ 0.661 in d = 3 (gamma_b = 0.199(2),
Clisby–Conway–Guttmann 2016; scaling relation gamma_b = gamma_11 + nu, Duplantier–Guttmann arXiv:1908.03872),
a = 1 in mean field. So (H) (2a > 1) is expected in d = 3, false in d = 2. Lemma 1 loses exactly the mutual
avoidance of the two halves of an arch, which in d = 2 is decisive (A(z_c) < infinity but sum B_M^2 = infinity)
and in d = 3 apparently not.

**B4. Where B stops.** (a) (H) is not proved; not even B_T(z_c) -> 0 is known on Z^d, d = 3, 4, as far as I found
(Section D). (b) Even with (H), the method gives n/log n, never n^{1-delta}: the gain per block is log(1/B_m) =
O(log m) because B_T(z_c) cannot be o(1/T) on average (CITED: Hutchcroft Remark 2.5 / Madras–Slade (3.1.14) give
sum_T e^{-tT} B_T(z_c) >= (1/2) log(1/t) - O(1)). (c) With only B_T <= C T^{-a}, a < 1/2: K_t ~ m^{1-2a} and one
gets only a speed bound v > (1-2a)/a (computation; nontrivial only for a > 1/3).

**B5. Variant tried: crossings instead of visits.** Using edge-planes E_h (crossings) instead of vertex planes:
ballistic bridges have many planes crossed <= 2/v times (DCH (3.2)), but the pieces between consecutive crossings
are again half-space excursions (arches), so the same arch generating function appears. No gain.

**B6. Variant tried: subcritical fugacity for excursions only.** Charging the arches with z = z_c e^{-t},
t ~ 1/m, using the trivial A(z) <= sum_M B_{2M}(z) <= B(z) <= 1/t: per-visit cost log K ~ log m, same order as the
gain a log m; fails (gives at best a speed bound). With Lemma 1 instead: K_t ~ m^{1-2a}, see B4(c). A bootstrap
(feeding a mass bound back into K_t) only improves constants: the mass bound obtained is linear in t.

## C. Reformulations (useful for the next wave)

**C1. Mass formulation (CITED + elementary).** Madras–Slade Prop. 4.1.8, Thm 4.1.14, Cor. 4.1.17: for z < z_c,
B_z(L) <= e^{-m(z) L}, m(z) = bridge mass = two-point mass, e^{-m L} <= G_z(L) <= chi(z)^2 e^{-m L}, m(z) -> 0
(Cor. 4.1.15, no rate). Hutchcroft Lemma 2.4: phi controls xi((1-e)z_c) = m. Conversely (elementary from
c_n >= mu^n, G_z(L) <= chi(z)^2 e^{-m(z)L} (MS (4.1.21)), chi(z) <= z^{-1} e^{2B(z)-2} (BDGS (2.20)) and
B(z_c e^{-t}) <= 1/(1-e^{-t}) <= 2/t for t <= 1 (Hutchcroft, proof of Prop. 2.2)): with z = z_c e^{-t},
  P_n(||gamma(n)||_inf >= L) <= 6 sum_{L'>=L} G_z(L') z^{-n}/c_n <= 6 z^{-2} e^{tn + 8/t} e^{-m(z) L}/(1 - e^{-m(z)}).
So, on the exponential scale, a rate for sub-ballisticity is the same problem as a bound
m(z_c e^{-t}) >= t psi(t) with explicit psi -> infinity ("correlation length << (z_c - z)^{-1}"). Under (H) the GF version of Proposition 3 (sum over all lengths, weight
theta^{(|gamma|+1)/m - V_i}) gives m(z_c e^{-t}) >= c t log(1/t) (sketch, not written up).

**C2. Irreducible-bridge formulation (elementary).** B(z, s) = sum_T B_T(z) e^{sT} = 1/(1 - Lambda(z, s)),
Lambda(z_c e^{-t}, s) = E_iSAB[e^{sH - tL}]; m(z_c e^{-t}) = sup{s : Lambda <= 1}. Sub-ballisticity with rate psi
needs E_iSAB[e^{-t(L - psi H)}] < 1: both a quantitative heavy tail of L (A3) and an exponential-scale
sub-ballisticity of irreducible bridges. Circular as a strategy, but it shows exactly what must be controlled.

**C3. (H) in renewal language (PROVED, result.md Sec. 6).** B_T(z_c) is the renewal sequence of the heights H of
irreducible bridges; (H) iff two independent critical infinite bridges a.s. share finitely many cut-plane heights;
(H) implies E_iSAB[H] = infinity, which is strictly stronger than DCH Theorem 2.5 (E_iSAB[L] = infinity).

## D. Attempts at (H) or at B_T(z_c) -> 0 on Z^3 (all FAILED)

* D1. Stickbreak for heights: suppose E_iSAB[H] < infinity (so cut planes have positive density in height).
  DCH's argument needs the walk to be thin relative to its height to find diamond points (cone condition); with
  E L = infinity the width over n renewal blocks is not controlled by their height. Fails.
* D2. Polygon/Kesten shifting: glue two bridges of height T (one shifted transversally by Y) into an arch or a
  polygon. Making them disjoint needs |Y| comparable to their lengths, and the connecting segment then costs
  z_c^{|Y|}, exponentially small. Only stretched-exponential information results (this is the source of
  Kesten's n^{2/(d+2)}).
* D3. HW-type inequalities: chi(z) <= z^{-1} exp(2B(z) - 2) and c_n >= mu^n give only *lower* bounds on
  sum_T B_T(z) (Madras–Slade (3.1.14)); the needed bound is an *upper* bound on sum_T B_T(z_c)^2.
* D4. Via arches: A(z_c) < infinity would follow from (H) by Lemma 1; the reverse is not available (a pair of
  bridges with different endpoints is not an arch). Circular.
* D5. BDGS Theorem 2.9 trick (pairs of bridges with the same endpoint -> walks to a neighbour): controls
  sum_x b_n(x)^2, not sum_T B_T^2 (different endpoints). No upper bound on polygons at z_c is known that would help
  (polygon GF at z_c finite is itself unproved; best rigorous q_{2n} <= mu^{2n}, BDGS (2.27)).

## E. Existence (E1) and comparison statements

Tried: E_{2n}|omega|^2 <= 4 (c_n^2/c_{2n}) E_n|omega|^2 (split at n, submultiplicativity) — useless since
c_n^2/c_{2n} is only known to be <= exp(O(n^{2/5} log n)) in d = 3 (Kesten). One- and two-step comparisons give
factors (2d-1)^k c_n/c_{n+k} > 1, no polynomial regularity. Lower comparisons need control of the number of
extensions of a walk jointly with its displacement (trapped walks exist with positive frequency). **Nothing
obtained on E1.** The computation lab reports (exact, n <= 36) that E_n/n is nondecreasing and E_n/n^2
nonincreasing; I have no proof of either.

## F. Pieces of literature checked for an existing explicit rate on Z^d

DCH 2013 (Question 1 asks for n^{2-eps}); Hutchcroft 2018 (φ "not quantitative"); Krachun–Panagiotis 2023 (first
quantitative bound, hexagonal only: "Previously, only a non-quantitative o(n) bound was known from the work of
Duminil-Copin and Hammond"); Madras–Slade 1993 Ch. 3-4; BDGS lectures; Glazman–Manolescu 2020 (hexagonal);
Semantic Scholar citation lists of DCH (through 2026: Grimmett 2026, Lazar 2025, Lehner–Lindorfer–Panagiotis 2024,
none quantitative on Z^d) and of KP (2 citing works, neither on Z^d sub-ballisticity). Literature agent:
Slade 2019 states no bound R_n <= C n^{0.99999} is proved. I found no explicit rate on Z^3.

## G. Numerics run by me (single core, < 1 min; heuristic only)

Exact enumeration of strict half-space walks on Z^3 to n = 17 (scratch code; bridges by height agree with the
computation lab's `computation/data/bridge_span.tsv` on all 153 entries n <= 17). Results:
* Lemma 1 holds coefficientwise: a_{n,M} <= sum_k b_{k,M} b_{n-k,M} for all n <= 17 (as it must). Ratio
  (total arches)/(total convolution) = 0.29, 0.22, 0.18, 0.16, 0.15, 0.15, 0.14, 0.14, 0.13, 0.13, 0.13, 0.12,
  0.12, 0.12 for n = 4..17: slowly decreasing, consistent with arches ~ n^{gamma_11 - 1} ~ n^{-1.39} versus the
  convolution ~ n^{-1 - nu(2a-1)} ~ n^{-1.19} (both heuristic).
* Truncated sums sum_{n<=17} b_{n,T} z_c^n: 0.697, 0.544, 0.405, 0.271, 0.159, 0.080 (T = 1..6): lower bounds for
  B_T(z_c), far from converged for T >= 3; useless for estimating a.
