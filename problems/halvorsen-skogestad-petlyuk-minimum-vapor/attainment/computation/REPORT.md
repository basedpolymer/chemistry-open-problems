# Agent C — computation and red team: numerical tests of the finite-stage attainment architecture

Claude Opus 5.5 (effort xhigh), Agent C (computation, red team), attainment wave, 2026-10-20.
Paths are relative to the problem folder `problems/halvorsen-skogestad-petlyuk-minimum-vapor/`; scripts in
`attainment/computation/scripts/`, logs and JSON in `attainment/computation/runs/`.

**Status of this file.** Everything below is **numerical evidence, not proof**: high-precision floating point
(mpmath, 60 to about 1500 significant digits), no interval enclosures. It tests the architecture of
`attainment/BRIEF-common.md` §3 and the statements of Agents D (`design/hypotheses.md`, `design/design.md`), L
(`limit/limit.md`) and A (`analysis/analysis.md`, `theorem.md`). **No counterexample to any claimed statement was
found.** Every discrepancy, caution and failure found is listed in §7.

---

## 1. What was computed

1. **An independent implementation of the forced design** F(δ, η, π) (`scripts/petnet.py`, class `Design`): HS's
   preferred split (paper.tex Theorem `achieve`) relative to each column's actual feed, Astra's slack deformation
   (attainment.md (3)–(6)) of every column with ≥ 3 components, forced binary columns (B-pinned: top total f_p − η_p;
   T-pinned: f_p + η_p; any pattern π ∈ {B, T}^{N−1}), vapor network R = H* + δ at the reboiler, side draws and
   junctions without duty. Checked in every case: total and CMO balances at every unit (≤ 1e-60), all flows > 0.
2. **The infinite-stage limit configuration, independently** (`limit_config`): limit net flows; non-binary junctions
   = Astra's Lemma 1 vector for (r, θ_{i+1..m−2}, s); binary junctions = pinch eigenvector of the forcing section;
   sharp/forced products; all section modes (roots of P_w(φ) = Π_i(α_i − φ)·Ψ(φ)) and the pinch index of every
   section from the mode tables. **Compared with D's JSON exports** (`scripts/compare_D.py`, `runs/compare_D.log`):
   12 configurations, N = 3..6, several patterns and D's random data — agreement ≤ 1.2e-54 in every quantity
   (flows, net flows, compositions, modes, pinch indices, pivots, κ, G_k, Ψ_k, Astra's cone weights A, B, C, D);
   D's limit points solve my G_∞ to ≤ 1.2e-54.
3. **An exact finite-stage solver in junction / net-flow variables.** Unknowns: the composition of every unit
   (junctions and products) and the net flow of every section (N·(#units + #sections) numbers: 36 for N = 3, 306
   for N = 6); equations: normalizations, Σw = V − L per section, unit component balances, and per section the
   N − 1 **exact** chart equations (product form, A's Lemma A.3 / L's Lemma 2.4)

       r above the pinch:  (ℓ_r·x_bot)(ℓ_*·x_top) − (φ_*/φ_r)^{n+1} (ℓ_r·x_top)(ℓ_*·x_bot) = 0,
       r below the pinch:  (ℓ_r·x_top)(ℓ_*·x_bot) − (φ_r/φ_*)^{n+1} (ℓ_r·x_bot)(ℓ_*·x_top) = 0,

   with the polynomial left eigenvectors ℓ_{r,i} = Π_{j≠i}(1 − φ_r/α_j) (smooth through zero trace flows; for a
   trace mode a multiple of L's regularized ℓ̃). Newton with an analytic Jacobian (validated against finite
   differences to 1e-27), linear solves in python-flint `arb_mat`, **two precision audits** (the working precision
   is raised until every nonzero unknown and every interior stage mole fraction exceeds 10^−(digits−30) in
   magnitude), continuation in n from the largest n down. Setting the factors (φ_*/φ_r)^{n+1} to 0 gives G_∞.
4. **Stage-level verification of every reported steady state** (`reconstruct`, `stage_check`): all stage profiles
   are reconstructed by a stable two-sided eigen-expansion, and every cross-section balance V y_k − L x_{k+1} = w,
   every stage balance and every unit balance (with the actual stream compositions, condenser and reboiler
   included) is re-checked componentwise in relative terms; every mole fraction > 0 and S(x) > 0 at every stage and
   unit are required. A zero of G_n is called a steady state only if it passes this check (§7.1 explains why).
5. **Cross-validations with independent codes.** (a) Agent N's double-precision stage solver `stagenet.py`:
   my N = 3 numbers reproduce the network.md §11.3 table (n = 16, 32, 64: impurities, ξ, ζ, ε_C, ε_A) and N's
   unforced numbers (§10.5, n = 16, 32). (b) Agent L's `limitsys.py` (ratio form, its own eigen-solver and trace
   normalization), its Design filled from my specification: L's |G_n| at my steady states ≤ 1e-95 … 1e-228
   (N = 3..6, n = 24..72), L's v_∞ = mine (`scripts/crosscheck_L.py`, `runs/crosscheck_L.log`). (c) D's exports
   (item 2).
6. **Jacobians.** DG_∞(v_∞) and DG_n(v_n) in my variables: singular values; major/trace classification (a column
   is trace iff its limit value is 0, a row is trace iff all its derivatives with respect to major columns vanish at
   v_∞); block sizes; split of the trace block by component; for each component k, the matrix A_k obtained by
   eliminating the trace net flows through the closed-end rows (node balances in the trace contents, oriented
   outflow − inflow) and a test of the nonsingular M-matrix property (Z pattern and entrywise nonnegative inverse);
   L's ratio-form scaling (rows divided by the two pivots).
7. **Hypotheses**, at the limit configuration: (H5) real positive distinct modes and ρ; (H6)/(A5) limit
   half-profiles ξ^bot_k, ξ^top_k (A's formula) componentwise ≥ 0 for k ≤ 600; (H8)/(A3) pivots > 0 with
   1·e_* = 1 = ℓ_*·e_*; (A4) closed coefficients = 0; (H7) Ψ^(k)_σ(α_k) > 0 and the sign of G^(k)_{x_J}(α_k) at every
   closed trace end, (iii) closed ends are trace nodes, (iv) every trace node is closed by some section, (v) nonzero
   limit net flows are major at both ends.

## 2. Main evidence (tables for the paper)

`scripts/paper_tables.py` → `runs/paper_tables.md` and `runs/paper_tables.json` (ALL PASS, about 40 s). Every row
is a steady state of the full stage-by-stage model (all stage and unit balances re-checked, max relative residual
≤ 3e-82; all mole fractions > 0; S > 0 everywhere).

**Table 1.** N = 3, α = (4, 2, 1), z = (1/3, 1/3, 1/3), q = 1, δ = 0.05 (R = H* + δ = 1.415722514,
H* = 1.365722514), prefractionator slack η_[1,3] = 0.025, pattern (B, T), η_1 = η_2 = 0.003 (the data of
network.md §11.3). Limit: J[1,2] = (0.69615425, 0.30384575, 0), J[1,3] = (0.32711348, 0.34246961, 0.33041690),
J[2,3] = (0, 0.55393069, 0.44606931); σ_min(DG_∞) = 0.0072915.

| n | impurities of products 1, 2, 3 | ‖x_J[1,2] − lim‖∞ | ‖x_J[1,3] − lim‖∞ | ‖x_J[2,3] − lim‖∞ | ε_C = w_C(T[1,3]), ε_A = −w_A(B[1,3]) | σ_min(DG_n) | min mole fraction | digits |
|---|---|---|---|---|---|---|---|---|
| 16 | 2.12e-4, 2.11e-2, 3.25e-3 | 1.17e-2 | 3.47e-3 | 4.67e-2 | 1.48e-3, 2.30e-3 | 0.0049615 | 1.8e-11 | 91 |
| 24 | 6.58e-6, 1.77e-2, 2.54e-5 | 4.87e-4 | 2.90e-4 | 8.22e-3 | 2.66e-4, 1.67e-4 | 0.006924 | 2.3e-16 | 96 |
| 32 | 2.15e-7, 1.77e-2, 4.68e-7 | 8.01e-5 | 4.53e-5 | 9.87e-4 | 3.93e-5, 8.07e-6 | 0.0072625 | 8.0e-22 | 101 |
| 48 | 2.31e-10, 1.77e-2, 2.01e-10 | 1.90e-6 | 8.92e-7 | 1.48e-5 | 7.16e-7, 1.77e-8 | 0.0072911 | 1.2e-32 | 112 |
| 64 | 2.47e-13, 1.77e-2, 8.80e-14 | 3.47e-8 | 1.60e-8 | 2.25e-7 | 1.27e-8, 3.89e-11 | 0.0072915 | 1.7e-43 | 123 |
| 96 | 2.83e-19, 1.77e-2, 1.68e-20 | 1.10e-11 | 5.06e-12 | 5.48e-11 | 4.02e-12, 1.88e-16 | 0.0072915 | 3.6e-65 | 144 |
| 128 | 3.25e-25, 1.77e-2, 3.22e-27 | 3.47e-15 | 1.60e-15 | 1.44e-14 | 1.27e-15, 9.11e-22 | 0.0072915 | 7.9e-87 | 166 |
| 192 | 4.28e-37, 1.77e-2, 1.18e-40 | 3.47e-22 | 1.60e-22 | 1.18e-21 | 1.27e-22, 2.13e-32 | 0.0072915 | 3.6e-130 | 209 |
| 256 | 5.63e-49, 1.77e-2, 4.30e-54 | 3.46e-29 | 1.60e-29 | 1.10e-28 | 1.27e-29, 4.98e-43 | 0.0072915 | 1.7e-173 | 253 |

The side product tends to its designed impurity (η_1 + η_2)/P_2 = 0.017682; ‖v_n − v_∞‖ ∝ 0.7765^n while the
largest adjacent mode ratio at the pinches is ρ_∞ = 0.7773.

**Table 2.** N = 4, α = (8, 4, 2, 1), z = 1/4, q = 1, δ = 0.05, η_nb = 0.005 (all three non-binary columns),
pattern (B, B, B), η_p = 0.003 (D's preset n4_1). σ_min(DG_∞) = 9.904e-4. Limit products: 1 pure, 2–4 with
designed impurities 0.0120, 0.0120, 0.01186.

| n | impurities of products 1, 2, 3, 4 | max over the 6 junctions of ‖x_J − lim‖∞ | σ_min(DG_n) | min mole fraction | digits |
|---|---|---|---|---|---|
| 16 | 5.10e-5, 1.55e-2, 2.31e-2, 1.94e-2 | 1.50e-2 | 0.0044799 | 1.3e-23 | 103 |
| 32 | 2.04e-8, 1.22e-2, 1.26e-2, 1.23e-2 | 1.07e-3 | 0.00095545 | 5.9e-44 | 123 |
| 64 | 3.55e-15, 1.20e-2, 1.20e-2, 1.19e-2 | 8.86e-7 | 0.00099035 | 7.5e-86 | 165 |
| 96 | 6.16e-22, 1.20e-2, 1.20e-2, 1.19e-2 | 1.37e-9 | 0.00099042 | 7.7e-128 | 207 |
| 128 | 1.07e-28, 1.20e-2, 1.20e-2, 1.19e-2 | 2.13e-12 | 0.00099042 | 7.8e-170 | 249 |
| 192 | 3.22e-42, 1.20e-2, 1.20e-2, 1.19e-2 | 5.11e-18 | 0.00099042 | 8.1e-254 | 333 |

‖v_n − v_∞‖ ∝ 0.8168^n; ρ_∞ = 0.8169.

**Table 3 (approach to the Petlyuk minimum H*).** N = 3 data of Table 1, pattern (B, T), R − H* = δ exactly.

| δ | η_[1,3] | η_1 = η_2 | n | impurities of products 1, 2, 3 | designed limit impurity of product 2 | σ_min(DG_∞) |
|---|---|---|---|---|---|---|
| 0.05 | 0.025 | 3e-3 | 128 | 3.25e-25, 1.77e-2, 3.22e-27 | 1.768e-2 | 7.29e-3 |
| 0.02 | 0.01 | 1e-3 | 192 | 2.29e-36, 5.96e-3, 8.21e-40 | 5.964e-3 | 2.69e-3 |
| 0.005 | 0.0025 | 2e-4 | 256 | 1.57e-47, 1.20e-3, 3.38e-52 | 1.199e-3 | 5.76e-4 |
| 0.001 | 5e-4 | 3e-5 | 384 | 2.14e-70, 1.80e-4, 8.43e-78 | 1.800e-4 | 9.04e-5 |

This is the order of limits of the theorem: n → ∞ at fixed (δ, η), then η → 0, then δ → 0: purity → 1 with
R − H* = δ → 0.

**η-dependence** (`scripts/eta_dependence.py`, `runs/eta_dependence.log`, ALL PASS): N = 3 data of Table 1 with
η_1 = η_2 = η, n = 16 … 512. n0 obs. = smallest tested n from which every larger tested n is physical;
C_n = ‖v_n − v_∞‖/ρ^n.

| η | n0 obs. | C_64 | C_128 | C_256 | C_512 | impurities of products 1, 2, 3 at n = 512 | ρ |
|---|---|---|---|---|---|---|---|
| 1e-2 | 16 | 0.885 | 0.873 | 0.873 | 0.873 | 6.1e-98, 5.7e-2, 6.2e-109 | 0.7773 |
| 3e-3 | 16 | 2.25 | 1.45 | 1.10 | 1.06 | 1.7e-96, 1.8e-2, 7.7e-108 | 0.7773 |
| 1e-3 | 16 | 11.3 | 11.2 | 11.0 | 10.5 | 4.4e-96, 6.0e-3, 1.7e-107 | 0.7773 |
| 1e-4 | 48 | 99.3 | 99.0 | 98.7 | 98.6 | 6.7e-96, 6.0e-4, 2.5e-107 | 0.7836 |
| 1e-5 | 48 | 983 | 982 | 982 | 982 | 7.0e-96, 6.0e-5, 2.6e-107 | 0.7843 |
| 1e-6 | 64 | 9.87e3 | 9.81e3 | 9.81e3 | 9.81e3 | 7.1e-96, 6.0e-6, 2.6e-107 | 0.7844 |

So ‖v_n − v_∞‖ ≈ (0.0098/η) ρ^n for small η (the constant of the inverse-function step blows up like 1/η, as the
product-form σ_min(DG_∞) ≈ 3.03 η predicts), and the observed threshold n0 grows slowly (16 → 64 as η goes from
1e-2 to 1e-6), consistent with n_0(η) ~ log(1/η)/log(1/ρ).

## 3. The limit system (tests of Theorem L and of the coordinator's Step 2)

In every case of §6 (62 cases, N = 3..8, all four patterns for N = 3, all eight for N = 4, random patterns up to N = 6), except the deliberately invalid `n3_invalid` where only (H6)/(H8) fail:
- the independently computed limit configuration solves G_∞ to the working precision (≤ 1e-54);
- **DG_∞(v_∞) is invertible**: σ_min between 7.9e-7 (D's rand_6_7, α spread 43) and 1.3e-2 in my product-form
  variables (my N-entry variables: 36 for N = 3 … 306 for N = 6);
- **block triangular [[Maj, ∗], [0, Tr]] with square diagonal blocks** (N = 3: 22 + 14; N = 4: 46 + 42; N = 5:
  83 + 92; N = 6: 136 + 170; N = 7: 208 + 282; N = 8: 302 + 434): every trace row has zero derivative with respect to every major variable at v_∞;
- the trace block splits into **one invertible block per component**; after elimination of the trace net flows
  through the closed-end rows, **every A_k is a nonsingular M-matrix** (off-diagonal ≤ 0, inverse ≥ 0; column sums
  ≥ 0 — zero up to rounding at some nodes — while row sums can be negative);
- **(H7)**: Ψ^(k)_σ(α_k) > 0 and G^(k) of the right sign (> 0 at a top end, < 0 at a bottom end) at every closed trace
  end (8 ends for N = 3, 25 for N = 4, …); (iii)–(v) hold;
- **(H6)/(A5) and (H8)/(A3)**: all limit half-profiles ≥ 0 (k ≤ 600), all pivots > 0.

**Role of the forcing** (`scripts/forcing_study.py`, `runs/forcing_study.log`; N = 3 data of Table 1):

| η_1 = η_2 | σ_min(DG_∞), product form | σ_min, L's ratio form | major block | trace block | junction ratio ξ at J[1,2] |
|---|---|---|---|---|---|
| 0.1 | 0.01166 | 0.0620 | 0.0937 | 0.0393 | 0.76716 |
| 3e-3 | 0.00729 | 0.0366 | 0.0147 | 0.0169 | 0.69615 |
| 1e-4 | 3.03e-4 | 0.0334 | 5.09e-4 | 5.76e-4 | 0.69319 |
| 1e-6 | 3.03e-6 | 0.0333 | 5.09e-6 | 5.76e-6 | 0.69309 |
| 1e-8 | 3.03e-8 | 0.0333 | 5.09e-8 | 5.76e-8 | 0.69309 |

In the product form σ_min(DG_∞) ≈ 3.03 η (the pivot of the forcing section at its product end is O(η)); in L's ratio
form it stays ≈ 0.0333 uniformly (L's remark), the η-dependence moving into ‖G_n − G_∞‖_{C¹} = O(ρ^n/η). Either way
n_0(η) grows like log(1/η)/log(1/ρ). **Without forcing** (η = 0; physical steady states from N's `stagenet.py`,
polished in high precision and stage-checked): σ_min(DG_n) = 9.6e-3, 8.2e-5, 9.8e-7, 1.4e-8, 2.0e-10, 2.9e-12 at
n = 8, 16, 24, 32, 40, 48 (geometric decay ≈ 0.59 per stage: the neutral mode of network.md §10.4), and the junction
ratio drifts (ξ = 0.449, 0.364, 0.337, 0.328, 0.323, 0.319), whereas with forcing σ_min(DG_n) → 0.0072915 and
ξ → 0.69615.

## 4. Convergence to the limit configuration

In every case ‖v_n − v_∞‖ (majors) decays like ρ_∞^n, ρ_∞ = the largest ratio between the pinch mode and an adjacent
mode at v_∞ (fitted rate / ρ_∞): 0.7765 / 0.7773 (Table 1), 0.8168 / 0.8169 (Table 2), 0.8289 / 0.8320 (N = 5 base),
0.9320 / 0.9295 (N = 5 mixed), 0.8383 / 0.8544 and 0.8902 / 0.9010 (N = 6), 0.9452 / 0.9452 (D's rand_5_5),
0.9422 / 0.9422 (D's rand_4_3). The ratio ‖v_n − v_∞‖/ρ^n stays bounded (no factor n), consistent with A's
Theorem A.5 (O(ρ^{n+1})); σ_min(DG_n) → σ_min(DG_∞). All traces are resolved: e.g. at n = 256 (Table 1) the smallest
stage mole fraction is 1.7e-173 and every mole fraction is positive.

**Unequal section lengths** (`scripts/unequal_lengths.py`, `runs/unequal_lengths.log`, ALL PASS): N = 3 and N = 4,
random n_σ ∈ [n, 4n] (5 seeds) and one section at n = 24 with all others at 192 (each section in turn): physical
steady states, σ_min(DG_n) ≈ σ_min(DG_∞), distance to v_∞ governed by min n_σ — supports A's statement that only
min_σ n_σ ≥ n_0 matters.

## 5. Tests of D's design claims

`scripts/redteam_design.py` (`runs/redteam_design.log`, ALL PASS) and `scripts/compare_D.py`:
- κ_p(η_nb, δ) − κ_p(η_nb, 0) = δ and κ_p(0, 0) = e_p to 1e-49 (N = 3..6, several q);
  dκ_p/dη_nb ∈ [−4.2, −0.69] (C_1 ≳ 4 for these data, as D reports); the largest common η_nb keeping
  min_p κ_p ≥ δ/2 at δ = 0.05 is 0.008–0.028 for these data.
- The forcing only increases the (H4) margin: Ψ_B(φ^T) (B-pinned) and Ψ_T(φ^B) (T-pinned) increase with η_p
  (η_p = 0 … 0.06; N = 3 and 4; 8 patterns).
- D's root-crossing bounds on η_p were never binding in the data tried (the cap f_p was always smaller); stress
  test at η_p = 0.5 f_p and 0.95 f_p: every check passes and physical steady states exist at n = 64.
- Outside D's region the architecture fails, as it should (case `n3_invalid`, §7.4).

## 6. All cases

`scripts/run_case.py <case>` (cases defined in `scripts/cases.py`; batches with `scripts/run_batch.sh`; table made
by `scripts/summarize.py` → `runs/summary.md`; one log and one JSON per case in `runs/`). Every case runs the full
battery of §1 items 1–7 and the finite-n solves. Columns: ρ = largest adjacent mode ratio at the pinches of v_∞;
"majors / traces" = sizes of the two diagonal blocks of DG_∞; "n0 obs." = smallest tested n from which every
larger tested n gave a physical steady state (Newton from v_∞, or continuation from the next larger n); impurities
of products 1..N at the largest n (side products and, depending on the pattern, product N carry the designed O(η)
impurities). Cases `D_*` use Agent D's own exported data and parameters (D's random generator: q < 0, q > 1,
volatility spreads up to 129; D's export uses η_p = 0.25 × cap, so with tiny feed fractions some designed impurities are large, e.g. product 1 of D_rand_5_6_TTBB is 98 % component 2 by design — the architecture is tested there, purity comes only with η → 0). `n3_invalid` is deliberately outside D's parameter region (§7.4): its FAIL is the
expected outcome.

| case | N | pattern | q | δ | ρ | σ_min(DG_∞) | majors / traces | A_k M-matrices | n0 obs. | largest n (digits) | impurities at largest n | ‖v_n − v_∞‖ | rate | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D_n3_B | 3 | BB | 1 | 0.05 | 0.7784 | 0.00972 | 22 / 14 | yes | 34 | 138 (163) | 4.5e-27, 9.0e-03, 8.9e-03 | 6.5e-16 | 0.7780 | PASS |
| D_n3_N | 3 | BT | 1 | 0.05 | 0.7773 | 0.00729 | 22 / 14 | yes | 34 | 138 (173) | 4.5e-27, 1.8e-02, 2.6e-29 | 1.1e-15 | 0.7724 | PASS |
| D_n4_1 | 4 | BBB | 1 | 0.05 | 0.8169 | 0.00099 | 46 / 42 | yes | 42 | 171 (306) | 8.8e-38, 1.2e-02, 1.2e-02, 1.2e-02 | 2.5e-16 | 0.8166 | PASS |
| D_n5_1 | 5 | BBBB | 0.5 | 0.05 | 0.7899 | 0.000334 | 83 / 92 | yes | 36 | 147 (438) | 2.8e-35, 1.0e-02, 1.0e-02, 1.0e-02, 9.9e-03 | 8.9e-16 | 0.7891 | PASS |
| D_rand_3_1 | 3 | BB | 0.3485 | 0.04148 | 0.7616 | 0.00311 | 22 / 14 | yes | 31 | 127 (195) | 1.1e-53, 6.0e-01, 7.0e-02 | 7.1e-16 | 0.7597 | PASS |
| D_rand_3_2_TB | 3 | TB | 1.208 | 0.02271 | 0.6791 | 0.000936 | 22 / 14 | yes | 22 | 90 (181) | 1.7e-02, 8.1e-36, 1.2e-02 | 3.1e-16 | 0.6791 | PASS |
| D_rand_4_3 | 4 | BBB | 1.512 | 0.09269 | 0.9422 | 8.93e-05 | 46 / 42 | yes | 145 | 580 (1134) | 1.0e-137, 2.2e-01, 7.7e-01, 5.3e-01 | 4.0e-18 | 0.9422 | PASS |
| D_rand_4_4_TBT | 4 | TBT | 1.401 | 0.06928 | 0.9306 | 3.27e-05 | 46 / 42 | yes | 120 | 481 (726) | 9.2e-03, 2.7e-79, 5.0e-01, 1.3e-61 | 2.8e-16 | 0.9306 | PASS |
| D_rand_5_5 | 5 | BBBB | 0.9469 | 0.03379 | 0.9452 | 0.000176 | 83 / 92 | yes | 75 | 300 (1453) | 9.9e-152, 1.4e-01, 9.5e-01, 1.0e-02, 5.0e-02 | 4.1e-11 | 0.9452 | PASS |
| D_rand_5_6_TTBB | 5 | TTBB | 1.31 | 0.0697 | 0.9947 | 6.76e-05 | 83 / 92 | yes | 75 | 300 (1197) | 9.8e-01, 3.0e-02, 6.4e-49, 1.0e-01, 4.1e-01 | 6.1e-05 | 0.9930 | PASS |
| D_rand_6_7 | 6 | BBBBB | -0.7279 | 0.04619 | 0.9738 | 7.94e-07 | 136 / 170 | yes | 75 | 300 (1406) | 3.6e-124, 6.7e-01, 6.9e-03, 1.3e-01, 4.9e-03, 7.3e-02 | 1.1e-06 | 0.9744 | PASS |
| D_rand_6_8_BTBTB | 6 | BTBTB | 0.4849 | 0.07875 | 0.9814 | 1.64e-05 | 136 / 170 | yes | 75 | 300 (1384) | 3.9e-50, 1.2e-02, 7.1e-26, 5.1e-02, 1.3e-24, 1.9e-01 | 1.4e-05 | 0.9831 | PASS |
| n3_BB | 3 | BB | 1 | 0.05 | 0.7784 | 0.00972 | 22 / 14 | yes | 16 | 128 (157) | 3.2e-25, 9.0e-03, 8.9e-03 | 8.1e-15 | 0.7780 | PASS |
| n3_BT | 3 | BT | 1 | 0.05 | 0.7773 | 0.00729 | 22 / 14 | yes | 16 | 128 (166) | 3.2e-25, 1.8e-02, 3.2e-27 | 1.4e-14 | 0.7720 | PASS |
| n3_TB | 3 | TB | 1 | 0.05 | 0.7784 | 0.0132 | 22 / 14 | yes | 16 | 128 (141) | 8.9e-03, 8.5e-15, 8.9e-03 | 8.5e-15 | 0.7780 | PASS |
| n3_TT | 3 | TT | 1 | 0.05 | 0.7773 | 0.0108 | 22 / 14 | yes | 16 | 128 (166) | 8.9e-03, 9.0e-03, 3.2e-27 | 1.5e-14 | 0.7723 | PASS |
| n3_base | 3 | BT | 1 | 0.05 | 0.7773 | 0.00729 | 22 / 14 | yes | 16 | 128 (166) | 3.2e-25, 1.8e-02, 3.2e-27 | 1.4e-14 | 0.7729 | PASS |
| n3_bigeta | 3 | BT | 1 | 0.05 | 0.7773 | 0.0117 | 22 / 14 | yes | 16 | 128 (170) | 4.0e-30, 3.8e-01, 3.2e-31 | 3.5e-15 | 0.7774 | PASS |
| n3_close | 3 | BT | 1 | 0.05 | 0.9617 | 0.00211 | 22 / 14 | yes | 32 | 384 (110) | 1.0e-09, 1.8e-02, 2.3e-09 | 9.7e-06 | 0.9639 | PASS |
| n3_dsmall | 3 | BT | 1 | 0.002 | 0.8027 | 0.000273 | 22 / 14 | yes | 64 | 256 (248) | 1.9e-47, 6.0e-04, 9.0e-52 | 3.2e-23 | 0.8023 | PASS |
| n3_equal_BB | 3 | BB | 1 | 0.02 | 0.8341 | 0.0024 | 22 / 14 | yes | 32 | 192 (185) | 1.8e-17, 8.4e-03, 8.4e-03 | 3.9e-16 | 0.8311 | PASS |
| n3_equal_BT | 3 | BT | 1 | 0.02 | 0.8341 | 0.002 | 22 / 14 | yes | 32 | 192 (223) | 1.8e-17, 1.7e-02, 1.6e-44 | 1.4e-15 | 0.8340 | PASS |
| n3_invalid | 3 | BT | 1 | 0.002 | 0.8106 | 0.00651 | 22 / 14 | yes | None | 128 (166) | 1.6e-24, 1.8e-02, -1.1e-26 | 3.3e-12 | 0.8105 | FAIL: H6/H7/H8, largest n not physical |
| n3_q-0.5 | 3 | BT | -0.5 | 0.05 | 0.8022 | 0.00745 | 22 / 14 | yes | 16 | 128 (172) | 1.9e-29, 1.8e-02, 8.9e-18 | 4.3e-13 | 0.7969 | PASS |
| n3_q-10 | 3 | BT | -10 | 0.05 | 0.8680 | 0.000859 | 22 / 14 | yes | 32 | 256 (302) | 8.3e-74, 1.8e-02, 2.9e-25 | 4.5e-17 | 0.8680 | PASS |
| n3_q-2 | 3 | BT | -2 | 0.05 | 0.8379 | 0.00375 | 22 / 14 | yes | 16 | 128 (181) | 3.8e-33, 1.8e-02, 7.7e-16 | 1.0e-10 | 0.8376 | PASS |
| n3_q0.5 | 3 | BT | 0.5 | 0.05 | 0.7358 | 0.00869 | 22 / 14 | yes | 16 | 128 (161) | 3.1e-26, 1.8e-02, 7.6e-24 | 2.4e-17 | 0.7328 | PASS |
| n3_q0 | 3 | BT | 0 | 0.05 | 0.7773 | 0.00892 | 22 / 14 | yes | 16 | 128 (166) | 1.3e-27, 1.8e-02, 7.1e-20 | 7.4e-15 | 0.7750 | PASS |
| n3_q1.5 | 3 | BT | 1.5 | 0.05 | 0.8022 | 0.00548 | 22 / 14 | yes | 16 | 128 (172) | 1.4e-24, 1.8e-02, 1.3e-29 | 7.6e-13 | 0.7977 | PASS |
| n3_q10 | 3 | BT | 10 | 0.05 | 0.8682 | 0.000549 | 22 / 14 | yes | 32 | 256 (302) | 9.8e-44, 1.8e-02, 9.2e-74 | 8.7e-17 | 0.8682 | PASS |
| n3_q3 | 3 | BT | 3 | 0.05 | 0.8396 | 0.00236 | 22 / 14 | yes | 16 | 128 (182) | 9.1e-24, 1.8e-02, 2.3e-33 | 2.2e-10 | 0.8393 | PASS |
| n3_veryclose | 3 | BT | 1 | 0.05 | 0.9914 | 0.000308 | 22 / 14 | yes | 500 | 4000 (167) | 1.2e-36, 1.2e-02, 4.6e-19 | 2.3e-13 | 0.9914 | PASS |
| n3_wide | 3 | BT | 1 | 0.05 | 0.4495 | 0.0023 | 22 / 14 | yes | 8 | 64 (248) | 1.8e-51, 1.8e-02, 1.1e-55 | 1.1e-22 | 0.4495 | PASS |
| n3_zA_1e-4 | 3 | BT | 1 | 0.02 | 0.9976 | 2.3e-05 | 22 / 14 | yes | 64 | 512 (456) | 1.5e-155, 4.0e-05, 7.6e-91 | 2.7e-07 | 0.9960 | PASS |
| n3_zB_1e-4 | 3 | BT | 1 | 0.02 | 0.9804 | 4.36e-07 | 22 / 14 | yes | 64 | 512 (452) | 5.6e-10, 9.1e-02, 9.6e-70 | 1.3e-06 | 0.9804 | PASS |
| n3_zlight | 3 | BT | 1 | 0.05 | 0.9158 | 0.000754 | 22 / 14 | yes | 16 | 128 (189) | 3.4e-07, 5.7e-02, 2.2e-37 | 2.5e-06 | 0.9166 | PASS |
| n3_zmid | 3 | BT | 1 | 0.05 | 0.9218 | 0.00526 | 22 / 14 | yes | 32 | 128 (162) | 7.6e-37, 7.4e-03, 7.8e-36 | 2.6e-06 | 0.9172 | PASS |
| n3_ztiny | 3 | BT | 1 | 0.05 | 0.9660 | 0.000251 | 22 / 14 | yes | 128 | 768 (768) | 2.3e-231, 1.7e-01, 2.5e-15 | 2.1e-13 | 0.9660 | PASS |
| n4_BBB | 4 | BBB | 1 | 0.05 | 0.7993 | 0.00196 | 46 / 42 | yes | 16 | 96 (209) | 6.3e-22, 1.2e-02, 1.2e-02, 1.2e-02 | 1.4e-10 | 0.7963 | PASS |
| n4_BBT | 4 | BBT | 1 | 0.05 | 0.7993 | 0.00157 | 46 / 42 | yes | 16 | 96 (222) | 6.3e-22, 1.2e-02, 2.3e-02, 3.4e-22 | 3.0e-10 | 0.7925 | PASS |
| n4_BTB | 4 | BTB | 1 | 0.05 | 0.7993 | 0.00238 | 46 / 42 | yes | 16 | 96 (209) | 6.3e-22, 2.3e-02, 1.5e-10, 1.2e-02 | 1.5e-10 | 0.7977 | PASS |
| n4_BTT | 4 | BTT | 1 | 0.05 | 0.7993 | 0.00185 | 46 / 42 | yes | 16 | 96 (222) | 6.3e-22, 2.3e-02, 1.2e-02, 3.4e-22 | 3.1e-10 | 0.7928 | PASS |
| n4_TBB | 4 | TBB | 1 | 0.05 | 0.7993 | 0.00198 | 46 / 42 | yes | 16 | 96 (202) | 1.2e-02, 4.7e-12, 1.2e-02, 1.2e-02 | 1.4e-10 | 0.7961 | PASS |
| n4_TBT | 4 | TBT | 1 | 0.05 | 0.7993 | 0.00159 | 46 / 42 | yes | 16 | 96 (222) | 1.2e-02, 4.7e-12, 2.3e-02, 3.4e-22 | 3.0e-10 | 0.7925 | PASS |
| n4_TTB | 4 | TTB | 1 | 0.05 | 0.7993 | 0.0024 | 46 / 42 | yes | 16 | 96 (202) | 1.2e-02, 1.2e-02, 1.5e-10, 1.2e-02 | 1.5e-10 | 0.7977 | PASS |
| n4_TTT | 4 | TTT | 1 | 0.05 | 0.7993 | 0.00186 | 46 / 42 | yes | 16 | 96 (222) | 1.2e-02, 1.2e-02, 1.2e-02, 3.4e-22 | 3.1e-10 | 0.7928 | PASS |
| n4_base | 4 | BBT | 1 | 0.05 | 0.7993 | 0.00157 | 46 / 42 | yes | 16 | 96 (222) | 6.3e-22, 1.2e-02, 2.3e-02, 3.4e-22 | 3.0e-10 | 0.7925 | PASS |
| n4_bigeta | 4 | BBT | 1 | 0.05 | 0.7993 | 0.00161 | 46 / 42 | yes | 16 | 128 (273) | 1.3e-32, 3.2e-01, 3.9e-01, 5.8e-33 | 5.7e-14 | 0.7929 | PASS |
| n4_close | 4 | BTB | 1 | 0.05 | 0.9449 | 0.000747 | 46 / 42 | yes | 64 | 384 (198) | 5.1e-13, 1.6e-02, 2.2e-10, 7.9e-03 | 1.2e-08 | 0.9448 | PASS |
| n4_equal_BBB | 4 | BBB | 1 | 0.02 | 0.8865 | 0.000391 | 46 / 42 | yes | 32 | 192 (342) | 1.2e-19, 8.8e-03, 1.9e-02, 7.9e-03 | 4.8e-11 | 0.8865 | PASS |
| n4_equal_TBT | 4 | TBT | 1 | 0.02 | 0.8865 | 0.000342 | 46 / 42 | yes | 32 | 192 (381) | 4.8e-03, 1.6e-11, 3.7e-02, 1.7e-41 | 4.8e-11 | 0.8865 | PASS |
| n4_q-5 | 4 | BBT | -5 | 0.05 | 0.8954 | 0.000318 | 46 / 42 | yes | 32 | 128 (300) | 1.8e-36, 1.2e-02, 2.3e-02, 8.9e-17 | 6.3e-09 | 0.8860 | PASS |
| n4_q0 | 4 | BBB | 0 | 0.03 | 0.8173 | 0.00195 | 46 / 42 | yes | 24 | 144 (267) | 4.7e-22, 1.0e-02, 1.3e-02, 5.7e-03 | 1.0e-12 | 0.8173 | PASS |
| n4_q13 | 4 | BTB | 1.3 | 0.03 | 0.9200 | 0.000848 | 46 / 42 | yes | 32 | 192 (213) | 3.4e-25, 1.3e-02, 6.5e-08, 9.9e-03 | 7.9e-08 | 0.9188 | PASS |
| n4_q6 | 4 | BBT | 6 | 0.05 | 0.8954 | 0.000108 | 46 / 42 | yes | 32 | 128 (301) | 1.7e-28, 1.2e-02, 2.3e-02, 9.9e-37 | 1.6e-08 | 0.8913 | PASS |
| n4_smallz2 | 4 | BBT | 1 | 0.05 | 0.8527 | 0.000154 | 46 / 42 | yes | 32 | 192 (358) | 8.1e-41, 5.0e-02, 5.2e-03, 4.6e-42 | 3.8e-13 | 0.8527 | PASS |
| n5_base | 5 | BBBB | 1 | 0.05 | 0.8320 | 0.000329 | 83 / 92 | yes | 16 | 96 (308) | 1.5e-23, 1.0e-02, 1.0e-02, 1.0e-02, 9.9e-03 | 1.8e-09 | 0.8289 | PASS |
| n5_mixed | 5 | BTBT | 0.6 | 0.03 | 0.9295 | 0.000489 | 83 / 92 | yes | 24 | 144 (288) | 5.8e-16, 1.6e-02, 2.5e-12, 1.3e-02, 7.2e-26 | 1.8e-05 | 0.9320 | PASS |
| n6_base | 6 | BBBBB | 1 | 0.05 | 0.8544 | 7.4e-05 | 136 / 170 | yes | 16 | 64 (323) | 4.3e-17, 1.2e-02, 1.2e-02, 1.2e-02, 1.2e-02, 1.2e-02 | 3.4e-06 | 0.8383 | PASS |
| n6_mixed | 6 | TBTBT | 0.8 | 0.03 | 0.9010 | 8.25e-05 | 136 / 170 | yes | 24 | 96 (280) | 5.0e-03, 2.6e-07, 7.9e-03, 2.3e-05, 2.0e-02, 1.9e-12 | 6.9e-05 | 0.8902 | PASS |
| n7_base | 7 | BBBBBB | 1 | 0.05 | 0.8712 | 1.97e-05 | 208 / 282 | yes | 16 | 32 (257) | 8.0e-10, 1.4e-02, 1.4e-02, 1.4e-02, 1.4e-02, 1.5e-02, 1.4e-02 | 1.4e-03 | 0.8567 | PASS |
| n8_base | 8 | BBBBBBB | 1 | 0.05 | 0.8888 | 4.07e-06 | 302 / 434 | yes | 16 | 16 (203) | 4.8e-06, 1.8e-02, 2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02, 2.2e-02, 2.0e-02 | 1.8e-02 | nan | PASS |

## 7. Discrepancies, cautions and failures found (all reported in the PROGRESS files concerned)

1. **G_n has non-physical zeros** (`scripts/spurious_zeros.py`, `runs/spurious_zeros.log`). N = 3, n = 24: the
   forced solution continued to η = 0 in the chart variables ends at a zero with |G_n| < 1e-61 whose stage profiles
   have min x = −2.35 and S = −2.05 at some stages, although junction and product compositions look reasonable.
   N = 5 (`n5_mixed`), n = 48: Newton started from the n = 24 solution converges to a zero at distance 0.17 from v_∞
   with min x = −31.5, coexisting with the physical zero at distance 0.016 found from v_∞. So hypothesis (b) of A's
   Lemma A.1 / L's Theorem 1.4 (b) (S > 0 along every reconstructed profile) is essential and is delivered only
   near v_∞; a numerical zero must always be stage-checked. Reported to A and L (their texts already require it).
2. **Pre-asymptotic n.** For small n (n3_zmid n = 16, n3_dsmall n = 32, n3_ztiny n = 64, N = 6 mixed n = 24; cases
   with ρ = 0.90–0.97 or small η) Newton from v_∞ converges to a non-physical zero or does not converge; for the next
   tested n everything is physical. Not a counterexample (the theorem is for n ≥ n_0), but n_0 can be large:
   it grows with 1/log(1/ρ) and with log(1/η) (§2, η-dependence table).
3. **Trace net flows take both signs for N ≥ 4** (e.g. N = 4, pattern BBT: component 3 flows down T[2,3]; component 1
   flows up B[1,3] at n = 32, sourced by the designed impurity of product 2; signs may change with n); all mole
   fractions stay positive. network.md §11.2's "every net flow has its natural sign" holds for N = 3 only. No step
   of A's or L's proofs uses trace-flow signs (positivity comes from the M-matrix argument of Lemma A.1).
4. **(A5)/(H6) is an independent, essential hypothesis** (case `n3_invalid`: δ = 0.002, prefractionator slack
   0.05 > η̄(δ), binary margin κ_[2,3] = −0.042). Every algebraic check passes (G_∞ = 0, DG_∞ invertible with
   σ_min 6.5e-3, block structure, M-matrices, κ signs), but the non-pinning section B[2,3] has pivot π_top = −0.093
   and a limit half-profile entry −1.25, and the zeros of G_n converging to v_∞ (n = 32, 64, 128) are non-physical
   (min x ≈ −7, S < 0, product-3 impurity −1e-26). So D's parameter region is genuinely needed and A's Lemma 6.4
   must (and does) use (A5). Reported to A and D.
5. **σ_min(DG_∞) ∝ η in the product form** (3.03 η, N = 3); uniform (≈ 0.0333) in L's ratio form (§3). Not a flaw:
   it fixes n_0(η) ~ log(1/η)/log(1/ρ) and the order of limits n → ∞, then η → 0, then δ → 0.
6. **Numerical artifacts of my own code, fixed** (not results): rounding residues ~1e-62 in b_{c,i} flagged by my
   (H7)(v) check (fixed by setting the exact zeros d_i = f_i, b_i = 0); a 0/0 in the right eigenvector when a
   trace root equals α_k to working precision (fixed with L's regularized ẽ); precision exhaustion for the widest
   data (interior stage mole fractions ~1e-1355; fixed by the second precision audit, up to about 1500 digits).
   FAIL lines of these kinds in superseded logs were removed by the final batch (`runs/batch_final.log`).
7. **D's scripts** (reported 18:28, fixed by D at 18:36): `export_config.py` crashed with CLI arguments (argv parsed
   at import of `check_design`), and `check_design`'s tolerance was fixed at import (false FAILs at dps < 60).
8. **Not tested / limits of this evidence:** no interval arithmetic (every number is high-precision floating
   point); finite sets of data (62 cases, N ≤ 8; N = 7, 8 only at small n); D's constants η_1, C_1 are not explicit, so "inside D's region"
   was checked through the explicit conditions (positive splits, κ_p > 0, (H4)), not through η̄(δ) itself.

## 8. Reproduction

From the project root (Git Bash):
`outils/wsl "cd chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/attainment/computation/scripts && <command>"`
with `<command>` one of the following (one thread each; runtimes on the project machine):

| command | what | runtime |
|---|---|---|
| `python run_case.py n3_base` | one case of cases.py (optional n-list, e.g. `16,32`) | 10 s (N = 3) to 10 min (D_rand_6_8) |
| `bash run_batch.sh 4 <case> …` | several cases in parallel, logs and JSON in `../runs/` | — |
| `python paper_tables.py` | Tables 1–3 → `../runs/paper_tables.md` | 40 s |
| `python eta_dependence.py` | threshold n_0 and constants versus η | a few min |
| `python forcing_study.py` | σ_min versus η; forced versus unforced | 1 min |
| `python unequal_lengths.py` | unequal section lengths, N = 3, 4 | 2 min |
| `python spurious_zeros.py` | non-physical zeros of G_n | 1 min |
| `python redteam_design.py` | D's margin identities, forcing monotonicity, stress tests | 1 min |
| `python run_D_export.py CASE ../runs/D_exports/CASE.json 60`, then `python compare_D.py` | D's exports versus my recomputation | 10 s |
| `python crosscheck_L.py` | L's G_n at my steady states | 2 min |
| `python equal_peaks.py` | feeds with equal highest peaks | 5 s |
| `python summarize.py` | summary table of all run_case JSONs | < 1 s |

Environment: WSL Fedora, Python 3.12, mpmath 1.3, numpy 2.x, python-flint 0.9 (`arb_mat`), scipy (only through
Agent N's `total/network/scripts/stagenet.py`, imported read-only by `forcing_study.py`). Every script prints
PASS/FAIL lines and exits non-zero on failure. Library: `scripts/petnet.py`.
