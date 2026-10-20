## Numerical demonstration (Agent C): finite-stage steady states of the forced generalized Petlyuk arrangement at fixed flows

All numbers below are numerical evidence (high-precision floating point, mpmath), not a proof. Every row is a steady state of the full stage-by-stage model (every stage and unit component balance re-checked with the actual stream compositions; all mole fractions > 0; S > 0 everywhere).

### Table 1. N = 3, α = (4, 2, 1), z = (1/3, 1/3, 1/3), q = 1, δ = 0.05, η_[1,3] = 0.025, pattern (B, T), η_1 = η_2 = 0.003

R = H* + δ = 1.415722514 (H* = 1.365722514); product flows ['0.330333', '0.339333', '0.330333']; limit product impurities ['0', '1.7682e-2', '0'] (designed: O(η)).

σ_min(DG_∞(v_∞)) = 0.0072915; limit junction compositions: J[1,2] = (0.69615425, 0.30384575, 0.0); J[1,3] = (0.32711348, 0.34246961, 0.3304169); J[2,3] = (0.0, 0.55393069, 0.44606931)

| n | impurity of products 1..N | ‖x_J[1,2] − limit‖∞ | ‖x_J[1,3] − limit‖∞ | ‖x_J[2,3] − limit‖∞ | ‖v_n − v_∞‖ (majors) | prefractionator trace flows | σ_min(DG_n) | max rel. stage residual | min mole fraction | min S | digits |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 2.12e-4, 2.11e-2, 3.25e-3 | 1.17e-2 | 3.47e-3 | 4.67e-2 | 4.67e-2 | ε_C = w_C(T[1,3]) 1.48e-3, ε_A = −w_A(B[1,3]) 2.30e-3 | 0.0049615 | 6.5e-83 | 1.8e-11 | 1.003 | 91 |
| 24 | 6.58e-6, 1.77e-2, 2.54e-5 | 4.87e-4 | 2.90e-4 | 8.22e-3 | 8.22e-3 | ε_C = w_C(T[1,3]) 2.66e-4, ε_A = −w_A(B[1,3]) 1.67e-4 | 0.006924 | 7.0e-83 | 2.3e-16 | 1.0 | 96 |
| 32 | 2.15e-7, 1.77e-2, 4.68e-7 | 8.01e-5 | 4.53e-5 | 9.87e-4 | 9.87e-4 | ε_C = w_C(T[1,3]) 3.93e-5, ε_A = −w_A(B[1,3]) 8.07e-6 | 0.0072625 | 6.8e-83 | 8.0e-22 | 1.0 | 101 |
| 48 | 2.31e-10, 1.77e-2, 2.01e-10 | 1.90e-6 | 8.92e-7 | 1.48e-5 | 1.48e-5 | ε_C = w_C(T[1,3]) 7.16e-7, ε_A = −w_A(B[1,3]) 1.77e-8 | 0.0072911 | 1.7e-82 | 1.2e-32 | 1.0 | 112 |
| 64 | 2.47e-13, 1.77e-2, 8.80e-14 | 3.47e-8 | 1.60e-8 | 2.25e-7 | 2.25e-7 | ε_C = w_C(T[1,3]) 1.27e-8, ε_A = −w_A(B[1,3]) 3.89e-11 | 0.0072915 | 6.7e-83 | 1.7e-43 | 1.0 | 123 |
| 96 | 2.83e-19, 1.77e-2, 1.68e-20 | 1.10e-11 | 5.06e-12 | 5.48e-11 | 5.48e-11 | ε_C = w_C(T[1,3]) 4.02e-12, ε_A = −w_A(B[1,3]) 1.88e-16 | 0.0072915 | 1.8e-82 | 3.6e-65 | 1.0 | 144 |
| 128 | 3.25e-25, 1.77e-2, 3.22e-27 | 3.47e-15 | 1.60e-15 | 1.44e-14 | 1.44e-14 | ε_C = w_C(T[1,3]) 1.27e-15, ε_A = −w_A(B[1,3]) 9.11e-22 | 0.0072915 | 1.2e-82 | 7.9e-87 | 1.0 | 166 |
| 192 | 4.28e-37, 1.77e-2, 1.18e-40 | 3.47e-22 | 1.60e-22 | 1.18e-21 | 1.18e-21 | ε_C = w_C(T[1,3]) 1.27e-22, ε_A = −w_A(B[1,3]) 2.13e-32 | 0.0072915 | 2.5e-82 | 3.6e-130 | 1.0 | 209 |
| 256 | 5.63e-49, 1.77e-2, 4.30e-54 | 3.46e-29 | 1.60e-29 | 1.10e-28 | 1.10e-28 | ε_C = w_C(T[1,3]) 1.27e-29, ε_A = −w_A(B[1,3]) 4.98e-43 | 0.0072915 | 4.6e-83 | 1.7e-173 | 1.0 | 253 |

Observed rate ‖v_n − v_∞‖ ∝ 0.7765^n; predicted ρ_∞ = max adjacent mode ratio at the pinch = 0.7773.

### Table 2. N = 4, α = (8, 4, 2, 1), z = 1/4, q = 1, δ = 0.05, η_nb = 0.005, pattern (B, B, B), η_p = 0.003

R = H* + δ = 1.322854785 (H* = 1.272854785); product flows ['0.247', '0.25', '0.25', '0.253']; limit product impurities ['0', '1.2000e-2', '1.2000e-2', '1.1858e-2'] (designed: O(η)).

σ_min(DG_∞(v_∞)) = 0.00099042; limit junction compositions: J[1,2] = (0.78097233, 0.21902767, 0.0, 0.0); J[2,3] = (0.0, 0.76719885, 0.23280115, 0.0); J[3,4] = (0.0, 0.0, 0.68359304, 0.31640696); J[1,3] = (0.30681046, 0.32986067, 0.36332887, 0.0); J[2,4] = (0.0, 0.43703406, 0.29747493, 0.26549101); J[1,4] = (0.24833644, 0.25166474, 0.25073373, 0.24926509)

| n | impurity of products 1..N | ‖x_J[1,2] − limit‖∞ | ‖x_J[2,3] − limit‖∞ | ‖x_J[3,4] − limit‖∞ | ‖x_J[1,3] − limit‖∞ | ‖x_J[2,4] − limit‖∞ | ‖x_J[1,4] − limit‖∞ | ‖v_n − v_∞‖ (majors) | σ_min(DG_n) | max rel. stage residual | min mole fraction | min S | digits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 5.10e-5, 1.55e-2, 2.31e-2, 1.94e-2 | 8.06e-3 | 6.20e-3 | 2.17e-3 | 1.50e-2 | 9.56e-3 | 2.96e-3 | 1.50e-2 | 0.0044799 | 7.3e-83 | 1.3e-23 | 1.019 | 103 |
| 32 | 2.04e-8, 1.22e-2, 1.26e-2, 1.23e-2 | 1.18e-4 | 5.06e-4 | 1.61e-4 | 1.07e-3 | 3.10e-4 | 5.69e-5 | 1.07e-3 | 0.00095545 | 2.6e-82 | 5.9e-44 | 1.012 | 123 |
| 64 | 3.55e-15, 1.20e-2, 1.20e-2, 1.19e-2 | 5.94e-8 | 6.94e-7 | 2.39e-7 | 8.86e-7 | 6.36e-7 | 1.07e-7 | 7.37e-7 | 0.00099035 | 3.3e-83 | 7.5e-86 | 1.012 | 165 |
| 96 | 6.16e-22, 1.20e-2, 1.20e-2, 1.19e-2 | 1.80e-11 | 1.00e-9 | 3.52e-10 | 1.37e-9 | 9.82e-10 | 1.65e-10 | 1.00e-9 | 0.00099042 | 2.3e-82 | 7.7e-128 | 1.012 | 207 |
| 128 | 1.07e-28, 1.20e-2, 1.20e-2, 1.19e-2 | 5.40e-15 | 1.53e-12 | 5.35e-13 | 2.13e-12 | 1.51e-12 | 2.56e-13 | 1.53e-12 | 0.00099042 | 2.6e-82 | 7.8e-170 | 1.012 | 249 |
| 192 | 3.22e-42, 1.20e-2, 1.20e-2, 1.19e-2 | 4.87e-22 | 3.63e-18 | 1.27e-18 | 5.11e-18 | 3.61e-18 | 6.15e-19 | 3.63e-18 | 0.00099042 | 1.5e-82 | 8.1e-254 | 1.012 | 333 |

Observed rate ‖v_n − v_∞‖ ∝ 0.8168^n; predicted ρ_∞ = max adjacent mode ratio at the pinch = 0.8169.

### Table 3. Approach to the Petlyuk minimum H* (N = 3 data of Table 1, pattern (B, T)): (δ, η_[1,3], η_p) → 0 jointly, R − H* = δ exactly

| δ | η_[1,3] | η_1 = η_2 | n | R − H* | impurity of products 1, 2, 3 | designed limit impurity of product 2 | σ_min(DG_∞) |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.025 | 0.003 | 128 | 0.05 | 3.25e-25, 1.77e-2, 3.22e-27 | 1.768e-2 | 0.007291 |
| 0.02 | 0.01 | 1e-3 | 192 | 0.02 | 2.29e-36, 5.96e-3, 8.21e-40 | 5.964e-3 | 0.002692 |
| 0.005 | 0.0025 | 2e-4 | 256 | 0.005 | 1.57e-47, 1.20e-3, 3.38e-52 | 1.199e-3 | 0.0005757 |
| 0.001 | 5e-4 | 3e-5 | 384 | 0.001 | 2.14e-70, 1.80e-4, 8.43e-78 | 1.800e-4 | 9.037e-5 |

Reading: at fixed (δ, η) the products 1 and N purify geometrically in n while the side product tends to its designed impurity (η_1 + η_2)/P_2; letting n → ∞, then η → 0, then δ → 0 gives purity → 1 with R − H* = δ → 0 (the order of limits of the theorem). The stage numbers in Table 3 grow because σ_min(DG_∞) ∝ η in the product-form normalization (n_0(η) ~ log(1/η)/log(1/ρ)).
