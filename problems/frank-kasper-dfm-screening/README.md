# Diblock foam model screen of the Frank–Kasper catalogue

**Status: partial result (numerical screen; not independently verified).** Curved-cell DFM minimiser validated
against Reddy et al. (2018); first DFM evaluation of the Frank–Kasper catalogue of Dutour Sikirić et al. (84 + 37
structures). Exact certificates only for *upper bounds* (no certified ranking).

**Attribution.** Method, code, runs, certificates and write-up: Claude Opus 5.5 (effort not recorded), solver
sub-agent of the orchestrator, 24 Oct 2026. No verification yet.

## Problem

Reddy, Buckley, Arora, Bates, Dorfman, Grason, *Stable Frank–Kasper phases of self-assembled, soft matter spheres*,
PNAS 115 (2018) 10233 (arXiv:1806.01726) minimise the **diblock foam model** energy
𝓕 = (𝒜² ℐ)^{1/3} (𝒜: mean cell area, ℐ: mean second moment of each cell about its centroid, both normalised by
the sphere of the mean cell volume) over curved periodic partitions with fixed cell topology and free cell volumes,
for 11 Frank–Kasper (FK) phases + BCC, and find σ minimal:

> "Among these nearly degenerate, fully unconstrained DFM structures, the σ phase overtakes A15 (minimal for fixed,
> equal volume) as the minimal energy phase (with next lowest energy for P)" — [1], Sec. II, p. 7 of arXiv v1.

The catalogue of FK structures of Dutour Sikirić, Delgado-Friedrichs, Deza (Acta Cryst. A 66 (2010) 602; 84
structures, 9 of them among those of [1]) and its part II (Struct. Chem. 23 (2012) 1103; 37 structures) had never
been evaluated in the DFM. **No paper states this as an open problem**: the question "is σ the DFM minimiser over the
whole catalogue?" is our framing (see `references.md` for the novelty check: none found).

## Result in one paragraph

TO BE COMPLETED (screen running).

## Milestone 1 — reproduction of [1] (`scripts/check_milestone1.py`: 0 FAIL)

See `result.md` §3. Fixed paper boxes, levels 0–4 (up to 4.9×10⁵ triangles), Richardson-extrapolated:
A15 1.0690374 (published 1.069037), C15 1.0692434 (1.069243), Z 1.0692501 (1.069250), C14 1.0690704 (1.069071),
δ 1.0696755 (1.069677): agreement ≤ 1.6×10⁻⁶ (≤ 4.3×10⁻⁷ for the cubic A15, C15). Kelvin foam 1.0972512
(1.097251) and Weaire–Phelan 1.0935413 (1.093541) for the equal-volume area problem. Our values are **lower** than
published for σ (−4×10⁻⁶; box ambiguous), pσ (−2×10⁻⁶), H (−1.6×10⁻⁵), P (−2.8×10⁻⁵), M (−4.0×10⁻⁵), μ (−1.3×10⁻⁴),
BCC (−8.7×10⁻⁶); for BCC, P, M, H, μ an **exact certificate** (integer arithmetic) proves that the published value
is not the infimum for the stated topology and box.

## How to reproduce

```bash
cd problems/frank-kasper-dfm-screening/scripts
pip install gemmi numba scipy numpy            # Python 3.12+; single thread is enforced
# catalogue data (Zenodo 10.5281/zenodo.13691342) -> ../certificates/sources/dutour2010/SpaceFullerene/
python3 check_certificates.py                  # exact certificates (1-3 min): PASS x5
python3 check_milestone1.py                    # Milestone 1 (5 min with stored runs): 0 FAIL
python3 run_dfm.py A15 sigma --box fixed --maxlev 4    # paper boxes (minutes)
python3 screen_pipeline.py ../work/screen_order_dedup.txt --win3 5e-4 --win4 3e-4   # the screen (hours)
python3 make_tables.py                         # tables
```

## Files

TO BE COMPLETED.
