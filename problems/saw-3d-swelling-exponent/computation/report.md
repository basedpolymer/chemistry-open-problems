# Computation lab: report (24 October 2026)

Written by the orchestrator from the computation lab's final hand-back (the harness did not let the sub-agent write
this file). Nothing here proves the target statement or any rung of the ladder. "Exact" means integer or rational
arithmetic at finite n; "numerics" is evidence only. Every verdict is regenerated into `data/*_output.txt` by the
scripts in `scripts/`.

## 1. Exact data

All cross-checks pass (`check_enum.py`, `check_aux.py`, `halfspace_tests.py`).

- `src/saw_enum.cpp` (48-fold symmetry reduction) covers n ≤ 19 in 330 s on 3 threads: c_n, c_n(x) per orbit, S_n and
  E_n as exact rationals; bridges b_n (convention 0 < x₁(ω_i) ≤ x₁(ω_n) for i ≥ 1, b₀ = 1); half-space walks; walks
  ending on x₁ = 0; the law of ω_k under n-step walks for every k < n.
- The n ≤ 21 run was aborted at 22 % (machine load 8–15 from other agents put it past the 3-hour limit).
- Further exact data: Gram sums Σ ω_i·ω_j (n ≤ 17), first/last-step counts (n ≤ 18), bridges by height, renewals and
  arches by maximal height (n ≤ 19).
- Cross-checked against a Python brute force, a naive C++ enumerator, OEIS A001412, Schram–Barkema–Bisseling 2011 P_N
  (equal to Clisby–Liang–Slade ρ_n), CLS polygon counts, Clisby–Conway–Guttmann 2016 bridges and irreducible bridges
  (a digit typo found in their arXiv Table 7) and OEIS A116904.
- The longest exact series used is c_n and P_n to n = 36 (Schram–Barkema–Bisseling 2011), recorded with its source in
  `data/reference_series.tsv`.

## 2. Testbed (exact)

**Holds for n ≤ 36:** E_n ≥ n; E_n/n nondecreasing; E_n/n² nonincreasing; E_{n+1} − E_n ≥ 1; superadditivity
(smallest margin 2/5).

**Holds for n ≤ 19:** E_n|ω_k|² ≥ E_k|ω_k|² and E_n[ω_k·(ω_n − ω_k)] ≥ 0; P(x₁(ω_n) = 0), b_n/c_n and
max_x P_n(ω_n = x) nonincreasing in n; c_n(k·d) unimodal in k along the axis, the face diagonal and the body diagonal.

**Fails:**
- c_n(x) is not maximised at the neighbours of 0 (c₃(e₁) = 4 < c₃(1,1,1) = 6; the maximiser sits near 0.5·√E_n).
- Monotonicity of c_n(x) fails along the axis (n = 9), coordinatewise (n = 5) and radially (n = 3).
- n ↦ E_n|ω_k|² is not monotone (n = 4).
- c_n is not log-concave and E_n is not convex (both at even n).
- **Refuted (two separate programs):** nonnegativity of all step correlations, a natural route to the diffusive lower
  bound. Σ over 16-step walks of s₁·s₁₆ = −325626, so E₁₆[s₁·s₁₆] = −54271/16686979329; negative again at n = 18. The
  row sums E_n[s_i·ω_n] ≥ 1 still hold for n ≤ 17.

## 3. Candidate lemmas of the solvers

Upper-bound file (`upper-bound/candidate-lemmas.md`):
- CL1 (arches bounded by pairs of bridges) and CL2 (renewal identity) hold for n ≤ 19.
- CL4, B_{T+1}(z) ≤ B_T(z), verified by exact rational enclosures: with the lab's counts only, z = 1/10 for T ≤ 7 and
  z = 3/20 for T ≤ 3; adding cited inputs (CCG b_n to n = 28; b_n ≤ μ^n, BDGS eq. (2.4); μ ≤ 4.7388 from
  Pönitz–Tittmann 2000), z = 1/10 for T ≤ 9, z = 3/20 for T ≤ 6, z = 9/50 for T ≤ 2. z = 1/5 cannot be decided from this
  data.
- CL3 and CL5 inconclusive: the extrapolation of the tail of B_T(z_c) is uncontrolled for T ≥ 3.

Lower-bound file (`lower-bound/candidate-lemmas.md`):
- CL3: margins Δ₁, Δ₂ ≥ 0 for n ≤ 19, as exact rationals (Δ₂ = 0.3788 at n = 19; both increase within each parity
  class). The per-prefix version of (I2) holds for prefixes of at most 2 steps and fails for every prefix length ≥ 3
  (e.g. +x +y −x): a proof has to average over prefixes.
- CL1, exact (n ≤ 12, m = 1): P(K = 2 | near, span ≥ 2) between 0.06 and 0.20.
- CL1, Monte Carlo (own implementation, validated against exhaustive enumeration): N = 1000: κ = 2.07(4), 2.58(2),
  3.21(2) at m = 8, 16, 32; N = 4000: κ = 2.59(9), 3.09(8), 3.83(8) at m = 16, 32, 64 (agrees with the lower-bound
  solver's own Monte Carlo). Slope of log κ against log m ≈ 0.31.
- CL2 at N = 10⁴: r + r′ ≈ 0.43 at j = 16 and 0.21 at j = 64, i.e. β ≈ 0.5.
- CL5 at N = 10⁴: median span/N^{1/3} = 10.5; P(span < N^{0.53}) ≈ 4 %.
- CL4: partial sum S(32) = 0.0654, increments decaying like m^{−1.78}; the argmax of G_z(x) is e₁, certified for
  z ≤ 3/20.

## 4. Numerics

The pivot code (checked against exact values up to n = 36), run over n = 125 … 32000, gives ν = 0.5876 ± 0.0005 (fit
with an n^{−0.528} correction, n ≥ 250; 0.5870–0.5879 as n_min varies), consistent with Clisby–Dünweg's 0.58759700(40).
The effective exponent from the exact series is still 0.598 at n = 35.

## Best rigorous statements, gaps

- Exact finite-n facts only, chiefly E₁₆[s₁·s₁₆] < 0 and the certified monotonicity of B_T(z) for z ≤ 3/20 and small T.
  No rung of the ladder was reached.
- Not run because of CPU contention: n = 20–21, Gram sums at n = 18. κ at N ≥ 10⁴ is statistically weak.

## Files

- `src/`: saw_enum.cpp, saw_naive.cpp, saw_gram.cpp, end_steps.cpp, halfspace_enum.cpp, cond_prefix.cpp, madras.hpp,
  madras_exact.cpp, madras_mc.cpp, pivot.cpp
- `scripts/`: check_enum.py, check_aux.py, lemma_testbed.py, gf_tests.py, halfspace_tests.py, madras_report.py,
  export_tables.py, fit_nu.py, sawdata.py, brute_force.py, run_pivot.sh
- `data/`: exact tables, reference series with sources, all script outputs

Credits: computation lab, Claude Opus 5.5 (effort xhigh).
