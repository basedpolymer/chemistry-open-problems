**Attribution.** All-M solver: Claude Opus 5.5 (effort xhigh). Completion, verification, corrections and
edition: GPT-6-Astra (effort xhigh), 26 October 2026. The finite-range solver and verifier are credited
in [../README.md](../README.md).

# Becker–Döring Hopf branches 1–3 for every M

## Status: solved on the three specified branches, verified with corrections

The original question and full reference are in [../README.md](../README.md). The complete argument is
[Theorem D of the integrated manuscript](../paper.pdf), with the detailed formulas in [result.md](result.md).
The audit and corrections are in [verification.md](verification.md).

Branch 1 starts at M = 15, is subcritical there, and is supercritical with stable small periodic orbits for
**every M ≥ 16**. Branches 2 and 3 start exactly at M = 44 and M = 87 and are supercritical on the centre
manifold for every subsequent integer M, with unstable periodic orbits in the full mass hyperplane.
There is no Hopf for 3 ≤ M ≤ 14 at any K > 0. All statements are local near the specified Hopf points.

For M ≤ 160 the recheck reconstructs the reduced matrices from the vector field and re-proves the coefficient
and full spectrum at all 337 points. For M ≥ 161, 562, 589 and 834 interval boxes certify the Hopf points,
coefficient sign and spectrum. Their **exact stored endpoint** is
`0.078811099999999995`, strictly greater than 1/sqrt(161). The previous decimal `0.0788111` was a rounded
label; exact rational coverage comparisons now enforce the actual bound.
Every endpoint is linked by a separate unique-zero certificate. A chart regular at K = infinity proves the
onsets and links every small-M point to these curves. The symbolic limit coefficient is proved exactly.

## Reproduce

Python 3.12 or later, python-flint 0.9, SymPy, mpmath, NumPy; three processes for the complete suite:

```bash
python scripts/recheck_all.py
```

Run the following from `scripts/` to select an individual check:

| Command | Claim |
|---|---|
| `verify_sweep.py ../certificates/l1h2_j1.csv 0 0.078811099999999995 --h2` | All 562 boxes, coefficient and spectrum on branch 1 |
| `verify_sweep.py ../certificates/l1_j2.csv 0 0.078811099999999995` | All 589 coefficient boxes on branch 2 |
| `verify_sweep.py ../certificates/l1_j3.csv 0 0.078811099999999995` | All 834 coefficient boxes on branch 3 |
| `h2j_sweep.py 2 ../certificates/l1_j2.csv OUT.csv 20 20.4 -0.3 0.25` | Entire spectral certificate for branch 2 |
| `h2j_sweep.py 3 ../certificates/l1_j3.csv OUT.csv 26 26.6 -0.15 0.25` | Entire spectral certificate for branch 3 |
| `check_branch_links.py` | Exact coverage, labels at infinity and every adjacent junction |
| `check_onsets.py J OUTDIR` | Onset, every onset-chart junction, all small-M points linked |
| `check_small_stability.py` | Exact Hurwitz proof for every K > 0, M = 3,…,14 |
| `check_formulas.py` | Symbolic characteristic identities, direct sums and dense resolvents |
| `limit_symbolic.py` | Exact limit formula and nonvanishing divisors |
| `asymptotics.py` | Certified limiting values and bounds for M ≥ 10^10 |
| `check_winding.py` | Counts with multiplicity and orientation; rejection of boundary zeros |

`recheck_all.py` also invokes the earlier verifier's dense code for every M ≤ 160 and checks all 99,840
integer points per branch in the overlap 161 ≤ M ≤ 100,000. Logs and generated certificates are in
`verification/recheck-2026-10-25/`. `run_all.sh` remains the generator for the original compact cover;
the complete completion audit is `recheck_all.py`.

## Scope

This closes both the spectral and large-M gaps for branches 1–3. The first branch is identified by
continuation from its limiting root; a global first-loss-of-stability assertion for every smaller K is not
needed and is not claimed. No simplicity assertion about the extra unstable eigenvalues is used.
[Higher branches](../higher-branches/README.md) retain their separate status.
