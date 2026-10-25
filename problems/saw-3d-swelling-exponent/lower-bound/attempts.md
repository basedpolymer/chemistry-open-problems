# Lower-bound solver: attempts and exact obstructions

Solver: Claude Opus 5.5 (effort xhigh), 2026-10-24.
Notation: as in `result.md` (strict SAW on Z^d, d = 3 unless said otherwise; span, L, I, H, U, T are Madras's).
Target: beat Madras's E_N|omega_N|^2 >= N^{4/9}/36 (CITED, Madras 2014 Prop. 1.1). **Nothing below beats it.**
Labels: PROVED / CITED / HEURISTIC / BLOCKED (the precise missing input is named).

## 0. Where the exponent 4/(3d) comes from, and a "worst case" consistent with everything known

Madras's only inputs are the following two.
- (V) Volume: span(omega) > N^{1/d}/3 for every walk (Lemma 2.2).
- (M) Reflection: T maps {||omega_N||_inf <= m} into {||omega_N||_inf >= N^{1/d}/3} and is at most (m+1)-to-one.

Write B^2 = E||omega_N||_inf^2. Markov gives P(near(2B)) >= 3/4. (M) then gives P(far) >= c/B, and
B^2 >= N^{2/d} P(far)/9 gives B^3 >~ N^{2/d}.

HEURISTIC "worst case" (WC). Put mass 1 - O(1/B) on endpoints within distance B = N^{2/(3d)} of the origin, roughly
uniformly, and mass of order 1/B at distance of order N^{1/d}.
- WC satisfies (M) for every m, all lattice symmetries, sup_x P(omega_N = x) -> 0 (DGHM Thm 1.2), the closing bound
  P(closing) <= n^{-1/2+eps} (Hammond 2019; in WC P(closing) ~ B^{-3} = n^{-2/3} in d = 3), and sub-ballisticity
  (DCH). I know of no published constraint on the endpoint law that excludes WC.
- PROVED (trivial): the one-dimensional projected constraint system behind (M), namely
  P(max_k omega_1(k) = h, omega_1(N) = a) <= P(omega_1(N) = 2h - a) together with (V), is satisfied by WC.

So any improvement must use the full geometry of the walk, not only its projection on one axis. This is why every
attempt below ends up needing information about cut levels, renewals or collisions.

## A1. Endpoint anti-concentration route — BLOCKED

Reduction (PROVED, elementary): if sup_x P_N(omega_N = x) <= C N^{-a}, then
P(||omega_N||_inf <= r) <= C(2r+1)^3 N^{-a}. Hence E >= c N^{2a/3}. **Beating 4/9 needs a > 2/3.**

Known (CITED):
- sup_x P_n(Gamma_n = x) -> 0 with no rate (DGHM 2016, Thm 1.2).
- P_n(||Gamma_n|| = 1) <= n^{-1/2+eps}, for neighbours of 0 only (Hammond, Trans. AMS 2019, Thm 1.2).
- DGHM 2016, after Thm 1.1, on extending the closing bound to a fixed general x: "For higher dimension, significant
  difficulties occur and it is unclear how to prove this result."

Exact data (computation lab, n <= 19): max_x P_n(omega_n = x) is attained at |x| ~ 0.5 sqrt(E_n) and decays like
n^{-1.7} (HEURISTIC fit), so the truth is a ~ 3 nu ~ 1.76.

Variant tried. Madras's per-endpoint injectivity gives (PROVED) P(omega_N = x) <= sum_{sigma} P(omega_N in ray_sigma(x)
at distance >= 2(N^{1/d}/3 - ||x||_inf) from x), the sum running over the 6 signed axis directions sigma.
(For each walk ending at x, take sigma with M_sigma = max_k sigma·omega(k) = span. Reflecting at the first time this
maximum is attained is valid, and it is injective for fixed (x, sigma), since M_sigma = (sigma·b + sigma·x)/2 where b
is the new endpoint.) This needs a *line* bound
P(omega_N in a given axis line) <= N^{-b}. Beating 4/9 again needs b > 2/3, and nothing of the kind is known.
Obstruction: no quantitative uniform-in-x endpoint bound on Z^3. The known methods (snake method, polygonal
invariance) rely on the endpoint being adjacent to the origin.

## A2. Improving Madras's multiplicity — reduced to two named inputs, not closed

Exact structure (PROVED, `result.md` Prop. 2). Two near walks with the same image differ by a fold at a renewal
depth j <= m of the head or of the tail below the top in direction I. After the fold the end piece of the walk is
translated by 2j e_I, and it must avoid the rest of the walk.

Conditional gain (PROVED, `result.md` Thm 3): an average multiplicity <= C m^{1-alpha} gives the exponent
4/(d(3-alpha)). In d = 3: 0.533 for alpha = 1/2, and 2/3 for alpha = 1.

HEURISTIC MC (`result.md` Section 4): the average multiplicity is 2.7 / 3.1 / 3.9 at m = 16 / 32 / 64 (N = 4000),
growing roughly like m^{0.3}.

Sub-attempts and where each stops:

- **A2.1 Per-slice injectivity / tubes.** For fixed omega_I(N) = a, T is injective (PROVED). Images of different
  slices are disjoint exactly when the images have a single admissible fold level in the window. The images land in
  the tube {|x_j| <= m for j != I} far out. To use this one needs P(endpoint in a far tube) << P(far endpoint), which
  is an anti-concentration statement of the same difficulty as A1 in dimension d-1. Circular.

- **A2.2 Renewal sparsity near the top.** For a fold at depth j to exist, the tail (or the head) must renew at
  depth j.
  - MC: about 25-40% of the depths <= m are renewal depths, so renewal sparsity alone gives at most a constant-factor
    or slowly growing gain in the accessible range.
  - Rigorous zero density of renewal *levels* would need, at least for bridges under P_iSAB, E_iSAB[height] = inf,
    equivalently B_T(z_c) -> 0. This is **not known on Z^d** (literature/quick-facts.md: known only on the hexagonal
    lattice, Krachun-Panagiotis Thm 2; secondary, per the literature agent).
  - DCH Thm 2.5 (E_iSAB|gamma| = inf) concerns renewal *times* and does not transfer.
  - Even granted, transferring from P_iSAB to the top of a uniform near walk needs DGHM-type comparison machinery
    (Lemma 3.3, shells), and the result would be non-quantitative: E/N^{4/(3d)} -> inf, not an exponent.

- **A2.3 Head-tail collisions.** MC: collisions remove about 3/4 of the renewal-admissible levels. A rigorous
  version needs "the end of the walk translated vertically by 2j hits the head for most j <= m", i.e. a lower bound
  on head-tail interpenetration at scale m. No tool is available.
  - The only deterministic collision: when omega(N) lies on the vertical line through 0, one value of j is forbidden.
  - "Anchors" (end patterns, Madras 1988) cannot force collisions, because the transverse offset omega_⊥(N) in
    [-m,m]^2 is not changed by local modifications.

- **A2.4 D_4 freedom at cut levels.** In the image psi, every admissible level is a cut plane, and the 8 transverse
  lattice symmetries can be applied independently at each (always valid). If psi has K admissible levels, the orbit
  has size up to 8^K, and only orbit elements whose transverse endpoint stays in the tube can be images of near walks.
  - The transverse endpoint of a random orbit element is d_0 + sum_i h_i d_i with h_i i.i.d. uniform in D_4, where
    the d_i are the transverse displacements of the pieces between cut levels.
  - The fraction of the orbit in the tube is about min(1, m^2 / sum |d_i|^2). The pieces lie inside a window of
    height m, so typically sum |d_i|^2 ~ K * O(1) <= m << m^2.
  - Result: at best a constant factor (the last piece alone gives <= 2 of the 8 values of h_K, i.e. a factor 4).

- **A2.5 Transverse pivots at the top point.** Rotating the tail about the vertical axis through omega(L) by
  g in D_4 preserves the e_I-profile.
  - The map is therefore injective with no multiplicity (PROVED: L and H can be read off the image), and the
    endpoint moves by (1 - g)(omega_⊥(L)).
  - But it is valid only if g(tail) avoids the head. MC: valid for about 20-50% of near walks, g = 180 degrees.
  - Validity cannot be bounded below, and |omega_⊥(L)| can be small (the top point can sit above the origin).
  - Combined with Madras's reflection (always valid), the 8 maps psi_g have mostly disjoint images. This gains a
    factor of at most 8: a constant.

- **A2.6 Point reflection through convex-hull vertices.** For every vertex p of conv(omega), reflecting the tail
  through p (y -> 2p - y) gives a SAW (PROVED: a generic supporting plane at p meets omega only at p). The endpoint
  goes to 2p - x.
  - Decoding needs p = (b + x)/2 with x in Lambda_m unknown in all d coordinates, so up to (m+1)^d candidates.
  - The number of images equals the number of far hull vertices: at most O(R^{3/2}) by Andrews' bound for lattice
    polytopes (classical, recalled from memory, not checked in a source; not used), and typically small.
  - No gain.

- **A2.7 Encoding the lost coordinate.** Madras loses the information omega_I(N) in [-m,m], about log m bits.
  - Inserting a code segment at the junction: the image has length N + k, and the probability ratio
    c_{N+k}/c_N ~ mu^k exceeds the number of codes that fit in the free layer above the cut. A code in a 2D layer
    needs mu_2^k >= m, so the loss is about m^{log mu / log mu_2} ~ m^{1.6}.
  - Encoding by length (straight vertical insertion of length j = omega_I(N) + m): this makes the map injective
    across slices, but the slices then land in different S_{N+j}. Probabilities are normalised per length, so the sum
    over slices still costs (2m+1) * P(far).
  - Deleting steps loses exactly the entropy it gains.
  - PROVED-level accounting: no encoding gains, because the per-slice map is already lossless. The loss is entirely
    the overlap of the images of different slices, i.e. the multiplicity of A2.2-A2.3.

- **A2.8 More reflection planes.** Diagonal lattice planes x_i ± x_j = const give 18 directions in total, each with
  a valid "reflect at the maximum" map. This is a constant number of maps, each with multiplicity ~ m: a constant gain.
  Composing reflections multiplies multiplicities:
  - top then bottom in direction I gives a translation of the endpoint by 2 * span, but decoding then needs the
    bottom cut level, with no window;
  - reflecting in I then J costs (m+1)^2 and gains at most a factor 2 in |endpoint|^2.

## A3. Span (non-compactness) route — BLOCKED

PROVED (`result.md` Thm 3 with alpha = 0): if P_N(span < N^{1/d+eta}) <= 1/4, then E >= c N^{4/(3d) + 4eta/3}.

What is needed is quantitative non-compactness: walks in a box of side N^{1/d+eta} must be a small fraction of S_N.
Obstructions:
- every comparison with mu^N loses exp(C n^{2/5} log n) in d = 3 (Kesten 1964, Thm 2, as reported in literature/quick-facts.md C3, not read by me; the proof there is only
  sketched) or exp(C sqrt n) (Hammersley-Welsh). The number of confined walks cannot be compared with c_N to
  precision 1/2 that way;
- Kesten's pattern theorem is local, and a positive density of finite patterns is compatible with a compact walk;
- HEURISTIC sketch: HW unfolding of a confined walk produces a bridge of height <= (side)^2/2 at cost 2^{O(side)}.
  One then needs bridges of small height to be rare, which is again a renewal-tail statement (irreducible bridge
  lengths), unknown;
- slab connective constants mu_W < mu are known only qualitatively (no rate in W).

Even P(span <= K N^{1/d}) -> 0 for fixed K would improve only constants.

## A4. Generating-function / Cesaro route — BLOCKED

- PROVED (bathtub argument): xi_2(z)^2 := sum_x |x|^2 G_z(x) / chi(z) >= c_d (chi(z) / sup_x G_z(x))^{2/d}.
- PROVED: with weights w_n = c_n mu^{-n} >= 1, sum_{n<=N} w_n P_n(omega_n = x) <= G_{z_c}(x).
- So **sup_x G_{z_c}(x) < inf would give, for every N, some n <= N with E_n|omega_n|^2 >= c N^{2/d}** (a
  subsequence version of nu >= 1/d), and a GF version.
- Obstruction (CITED, Madras-Slade p. 18): "it has not yet been proved rigorously that G_{z_c}(0,x) is even finite
  for d = 2, 3 or 4, for any value of x != 0." For x = e_1 we have G_{z_c}(e_1) = (mu/d) sum_m m p_m mu^{-m}
  (PROVED: for n >= 3, an n-step SAW from 0 to e_1 plus the missing edge is a polygon through {0, e_1}; summing over
  edge directions gives c_n(e_1) = (n+1) p_{n+1}/d). Finiteness needs
  p_m mu^{-m} = m^{-theta} with theta > 2 in a summable sense.
- What is known: theta >= 3/2 on a density-one set in d = 2 only (Hammond 2018); in d = 3, Madras 1991 "would assert
  lim theta_n >= 1" (Hammond's description; not read).
- A weaker hypothesis sup_x G_z(x) <= chi(z)^{1-delta} also fails to follow from known bounds:
  - Kesten 1964 (1.9) gives only p_n <~ mu^n (secondary, per literature/quick-facts.md C3);
  - closing <= n^{-1/2} needs control of the z-law of the length, and c_n mu^{-n} is uncontrolled up to
    exp(C n^{2/5} log n), so GF statements do not transfer to fixed n anyway.

## A5. Local monotonicity (folklore route to L3) — BLOCKED, exact reduction recorded

PROVED identities:
- E_{n} - E_{n-1} >= 1 for all n implies E_n >= n (L3).
- E_n|omega_n|^2 - E_n|omega_{n-1}|^2 = 1 + 2 E_n[omega_{n-1}·s_n] = 2 E_n[omega_n·omega_1] - 1. The last equality
  uses reversal invariance: E_n|omega_{n-1}|^2 = E_n|omega_n - omega_1|^2.

So E_n - E_{n-1} >= 1 follows from two inequalities:
- (I1) E_n|omega_{n-1}|^2 >= E_{n-1}|omega_{n-1}|^2. This is equivalent to Cov_{P_{n-1}}(|eta_end|^2, ext(eta)) >= 0,
  where ext(eta) is the number of free neighbours of the endpoint.
- (I2) E_n[omega_{n-1}·s_n] >= 0. By lattice symmetry this is equivalent to E_n[omega^{(1)}(n) | omega(1) = e_1] >= 1.
  This is a **mean version of the "reflection principle" that Madras 1988 (p. 691) calls open**: "Is it true that at
  least half of the SAWs in F_N(P) have omega^{(1)}_N >= p^{(1)}_n? Such a 'reflection principle' would be a most
  useful result, but nobody can prove it."
- Likewise, superadditivity E_{k+l} >= E_k + E_l follows from E_n|omega_k|^2 >= E_k|omega_k|^2 together with
  E_n[omega_k·(omega_n - omega_k)] >= 0 (PROVED, using reversal).

Computation lab (exact, finite n only): E_{n+1} - E_n >= 1 holds for n <= 35. Superadditivity holds for n <= 36.
Both auxiliary inequalities hold for n <= 19.

Partial exact fact (PROVED). Let F = {omega in S_n : omega(1) = e_1} and G = {omega in F : 2e_1 not in omega}. The
map omega -> omega[0,1] ∘ U_1(omega[1,n]) is an involution of G, where U_1 is the reflection through x_1 = 1. It is
valid because U_1(omega(k)) = 0 iff omega(k) = 2e_1, and it maps G to G. It sends omega^{(1)}(n) to
2 - omega^{(1)}(n). Hence E_n[omega^{(1)}(n) | G] = 1 exactly, and **(I2) is equivalent to
E_n[omega^{(1)}(n) | omega(1) = e_1, 2e_1 in omega] >= 1.**

Iterating with planes x_1 = 2, 3, ... gives the same exact statement on further "valid reflection" subsets. But every
step leaves a residual set (reflection blocked by the head) whose conditional mean is uncontrolled. Making the choice
of reflection plane canonical *and* involutive on all of F is exactly Madras's open problem: a first-valid-time rule
is not preserved by the reflection.

Attempt at (I2) by reflection at the maximum: the head (containing the first step) is preserved, and the image ends at
2 M_1 - omega_1(N) >= M_1 >= 1. But the map has the same window multiplicity as A2, so it bounds
P(omega_1(N) <= 0 | s_1 = e_1) only by (range) x P(omega_1(N) >= 2), not by 1/2. Blocked for the same reason.

## A6. Pivot detailed balance — no inequality extracted

PROVED identity: for each k, E_N[ sum_{g valid at k} ( |omega(k) + g(omega_N - omega(k))|^2 - |omega_N|^2 ) ] = 0
(pivots are bijections between valid configurations).
- Averaging over the full hyperoctahedral group would give the Pythagorean value
  |omega(k)|^2 + |omega_N - omega(k)|^2, but the acceptance indicator is correlated with geometry.
- Turning the identity into a bound needs that near walks have at least as many valid "outward" pivots as far walks
  have valid "inward" ones. That is an acceptance-rate comparison, with no tool for it.

## A7. Bridges, renewal theory, Hammersley-Welsh — BLOCKED

- Unfolding into a bridge costs exp(C sqrt N) and does not force displacement: a bridge can have small height.
- Kesten's renewal identity (DCH Lemma 2.2, read) and DCH Thm 2.5 control renewal *times*. Displacement needs renewal *heights*:
  E_iSAB[height] = inf, i.e. B_T(z_c) -> 0, which is open on Z^d. Even "irreducible heights have infinite mean"
  gives no rate.
- End patterns (Madras 1988, Thm 2.1, lim inf c_N(P)/c_N > 0) give O(1)-scale smoothness of the endpoint law only.

## A8. Sub-walk / midpoint comparisons — weaker than Madras

- DGHM Prop. 1.3 (sup_x P(Gamma_{n/2} = x) <= C n^{-1/2}) and the volume argument give only
  E_N|omega(N/2)|^2 >= c N^{1/d} (N^{1/3} in d = 3), below 4/9.
- Comparing the law of omega[0,k] under P_N with P_k involves c_k c_{N-k}/c_N, which lies in [1, e^{C sqrt N}] and is
  uncontrolled.

## A9. Moments and interpolation — BLOCKED

- Madras's higher-moment bounds, (E||omega||_inf^p)^{1/p} >= N^{p/((p+1)d)}/6, tend to N^{1/d} but cannot be
  pushed down to p = 2 by Holder without an upper bound on high moments.
- The only upper bound is DCH (no rate) or the trivial |omega_N| <= N, both useless for this.

## Summary of exact obstructions (for the next wave)

1. **(H2): sublinear average multiplicity of Madras's map.** Any power saving improves the exponent. Bounded
   multiplicity gives N^{2/d}. MC suggests multiplicity ~ m^{0.3}. It needs (a) renewal-level sparsity near the top
   of the walk (relates to B_T(z_c) -> 0, open on Z^d) and/or (b) head-tail interpenetration at scale m.
2. **(H1): span >= N^{1/d+eta} with probability >= 3/4.** Needed to go beyond 2/d with this method. Equivalent in
   spirit to quantitative non-compactness.
3. **G_{z_c}(x) < inf uniformly in x.** It would give subsequence and GF versions of nu >= 1/d. Open: Madras-Slade
   p. 18, and polygon exponent theta > 2.
4. **Mean reflection principle E_n[omega^{(1)}(n) | omega(1) = e_1] >= 1, plus Cov(|end|^2, #free neighbours) >= 0.**
   Together they give L3. The first is a mean version of Madras's 1988 open reflection principle.
