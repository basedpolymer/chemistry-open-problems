# Verification status (solver's internal audit, superseded by `../verification.md`)

Date: 28 October 2026. Internal reviewer: the solver, GPT-6-Astra.

**Verdict: not yet performed.** This file is not a verifier report. The repository protocol reserves `solved` for a result checked by a separate verifier. The current status is `solved (preliminary)`.

## Internal mathematical audit

1. **Original model:** the ordered gain sum and the factor 2 in the loss term agree with equation (37). Time is the extent-of-reaction time of the paper, not the usual physical-time multiplicative kernel convention.
2. **Support:** subtracting `ell` makes the reaction additive. Initial size `k` generates exactly sizes `ell + (k-ell)j`. Off-support equations have zero initial values and zero gain by induction.
3. **No circular use of mass:** the triangular construction initially uses a candidate denominator. The positive generating series is then summed to prove its actual first moment equals that denominator.
4. **All indices, all parameters:** Lagrange inversion proves the coefficient formula for every `j` and every integer `k > ell >= 1`, including non-integer `p=(2k-ell)/ell`. The finite exact checks are not used as a universal proof.
5. **Convergence:** the real inverse branch turns at `H=1/p`; positive coefficients and their Stirling asymptotic give the convergence radius. An exponential generating-series bound yields uniform moment convergence on every compact subcritical time interval, including zero by removal of the apparent singularity.
6. **Critical endpoint:** the exact `j^(-5/2)` concentration tail has finite mass and infinite second moment. A summable `j^(-3/2)` bound proves uniform mass convergence up to the endpoint. No termwise moment differentiation is attempted where the second moment is infinite.
7. **Uniqueness and maximality:** for any competing pre-gelation solution, uniformly summable second moments justify the mass balance; the system is then triangular and fixes every component. Thus an earlier singularity and a pre-gelation extension through the threshold are excluded.
8. **Scope:** the result does not assert a giant cluster in the stochastic system or a post-gelation rule. `t_g` means the deterministic second-moment threshold, matching the definition in the primary source.

## Exact computation

`scripts/verify_exact.py` completed successfully with **85 PASS, 0 FAIL**. Saved output: `certificates/exact-checks.json`. No floating-point arithmetic occurs in these checks.

## Outstanding gate

A fresh-context verifier must inspect the original statement, rederive the proof, rerun the exact checks, and supply its own verdict before the status can be changed to `solved`.
