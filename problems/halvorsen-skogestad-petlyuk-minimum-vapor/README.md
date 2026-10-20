# Halvorsen–Skogestad conjecture: the generalized Petlyuk arrangement minimizes total vapor

**Attribution.** Theory: Claude Opus 5.5 (effort xhigh) (Agent T: setting, cut identity, sign rule, extreme peaks,
FTC formula for all N, four-component case, attainment, manuscript). Computation: Claude Opus 5.5 (effort xhigh)
(Agent C: enumeration, Underwood model, counterexample search, exact certificates, and the general-N proof in the
strict model). General N: Claude Opus 5.5 (effort xhigh) (Agent G: recovery ordering for every column,
merge-free theorem in the permissive model). Independent verification: Claude Opus 5.5 (effort xhigh) (verdicts,
corrections and the five-component permissive counterexample: `verification.md`, `verification/`).

**Status: solved in a precisely defined class and model — networks of two-section columns with interval streams, in
the strict infinite-stage Underwood model (TG's formulation in which a column with two or more nominally
distributing components runs at its Underwood vapor), without reverse flows, for every number of components N;
independently verified, with corrections. In the permissive (TG-literal) model: partial result, and the bound is
false for N = 5 (interval-certified counterexample found by the verifier). Open outside the class and model:
several feed stages per column, stream splitters, recycles, reverse net flows, and degenerate operating points at
which a nominally distributing component is absent from one product (there the strict model is more restrictive than
infinite-stage physics).** The Halvorsen–Skogestad conjecture as stated ("any adiabatic distillation arrangement") is
therefore settled only within this class and model.

## Problem

I. J. Halvorsen, S. Skogestad, *Minimum energy consumption in multicomponent distillation. 3. More than three
products and generalized Petlyuk arrangements*, Ind. Eng. Chem. Res. 42 (2003) 616–629, doi:10.1021/ie0108651.

Verbatim open statements:
- Abstract: "We also conjecture that the minimum energy requirement for the generalized extended Petlyuk arrangement
  is lower than the minimum-energy requirement for any distillation configuration when we consider conventional
  adiabatic sections and no internal heat exchange."
- §6.3 (p. 626): "We assume constant relative volatilities, constant molar flows, constant pressure, and no internal
  heat integration. Then the generalized adiabatic Petlyuk arrangement has the lowest need for vaporization compared
  to any other adiabatic distillation arrangement for separation of an arbitrary feed mixture into its pure
  components. We have not presented a complete proof, so the above conclusion is a conjecture."
- §5.3 (p. 624), on the Petlyuk arrangement itself: "We have not carried out a detailed proof for the general
  N-component M-product case for other than the far left and right peaks."

Setting: ideal N-component feed (constant relative volatilities α₁ > … > α_N, composition z, liquid fraction q),
constant molar overflow, infinitely many stages, sharp separation into pure products; objective = total vapor
generated (heat input). With θ_j the roots of Σ α_p z_p/(α_p − θ) = 1 − q, the Petlyuk value is
V_FTC = max_j H_j, H_j = Σ_{p>j} α_p z_p/(θ_j − α_p) (the highest peak of the V_min diagram, referred to the reboiler).

## Result

Competitor class: **two-section column networks** — columns with one feed stage, interval streams [i,m], split
[i,k]/[l,m] with l ≤ k+1, side draws merging any number of links, thermal couplings or heat exchangers on any link,
heaters/coolers at any node, saturated-liquid pure products, single feed. This contains Agrawal's regular-column
configurations (8 / 152 / 6128 for N = 3 / 4 / 5). Each column obeys the Underwood model of Tumbalam Gooty, Agrawal &
Tawarmalani (arXiv:2010.12113), with nonnegative net flows (no reverse flows). Two versions: *strict* (every column
with two or more nominally distributing components runs exactly at its Underwood vapor — the infinite-stage
behaviour whenever those components actually distribute; at degenerate points, where one of them is absent from a
product, the strict model forbids excess vapor that physics would allow) and *permissive* (TG literal: excess vapor
above the threshold allowed in every column).

1. **Strict model, every N:** every network satisfies R ≥ V_FTC at every feasible point (Agent C's sign-graph proof,
   checked by T; independent flux proof by G).
2. **Permissive model:** R ≥ V_FTC for all 152 regular-column configurations with N = 4 (T); for every network and
   every N when the highest peak is the first or the last (T); for every regular-column configuration and every N
   whose mixed sub-network for the highest peak has no merge node, under a non-degeneracy assumption (G); and for
   the FTC itself for every N (T) — i.e. HS's formula "minimum vapor of the generalized Petlyuk arrangement = highest
   peak" is proved for every N in this model (without reverse flows). For N = 5 the remaining cases are 20 (basic
   configuration, j) pairs = 864 (configuration, j) pairs (720 of the 6128 configurations have such a j); the
   classification "settled by the sign rule" of the others is a floating-point LP computation (integer data).
   **The bound fails for the pair (161, 3):** the verifier found a permissive-model point of basic configuration 161
   (all thermal couplings, q = 0) with R/V_FTC ≈ 0.99946, all flows ≥ 1e-6, all inequality slacks ≥ 1e-6, and proved
   the existence of the exact feasible point by an interval Krawczyk test (`verification/krawczyk_certificate.py`,
   data in `verification/permissive_counterexample_N5.json`). The excess vapor sits in a column whose root at the
   split is active, which an infinite-stage column cannot do; no adiabatic arrangement with pure products realizing it
   is known (finitely many stages would leave impurities, and a vapor bypass moves components; see the paper,
   Section "A five-component counterexample", and `verification.md`). So it is an artifact of the permissive
   (TG-literal) formulation, which is not a faithful model of the conjecture; the strict model is the relevant one. The other 19 pairs are undecided. Numerical searches
   (computation/): multistart ipopt with a-posteriori verification; SCIP only as a cross-check at two points (N = 4).
3. **Attainment, every N:** the FTC has a point feasible in both models with R = V_FTC (explicit construction:
   preferred split everywhere, excess vapor only in the binary columns).

Hence, in the strict model, the generalized Petlyuk arrangement minimizes total vapor among all two-section column
networks, for every N and every q. The key tools: the exact cut identity for the Underwood flux
Λ_s(θ_j) = V_s − Σ α_p w_{s,p}/(α_p − θ_j) (conserved at junctions, the feed term vanishing at θ_j), a per-column sign
rule, and the ordering of recoveries.

## Reproduction

From the project root (Git Bash), with `outils/wsl` running WSL Fedora with the scientific Python environment:

```
outils/wsl "bash chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/theory/scripts/run_theory_checks.sh"        # ~2 min
outils/wsl "bash chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/computation/scripts/run_certificates.sh"    # ~40 s
py -3 templates/build_tex.py problems/halvorsen-skogestad-petlyuk-minimum-vapor/paper.tex   # from chemistry-open-problems/
```

No proof of the main results depends on floating point: the scripts check identities, the case analysis for N = 4, the
attainment construction, the cut identity (exact rationals) and the combinatorial core of the strict proof; the
numerical searches (computation/) are evidence only. Exceptions, stated as such: the N = 5 classification of the
permissive cases settled by the sign rule (HiGHS LP, `theory/scripts/signcheck.py`), and the permissive
counterexample, which is a computer-assisted proof in interval arithmetic.

Independent verification (from `chemistry-open-problems/`; WSL as above, paths relative to the problem folder):

```
cd verification && python indep_checks.py                     # exact cut identity, FTC point, TG optima, Table 6.1 (~10 s)
python indep_nlp.py run 4 P2 8 permissive                     # verifier's own NLP (N = 3, 4), a few minutes per run
python krawczyk_certificate.py logs/candidate_161_3_t0_refined.json   # permissive counterexample, interval proof (~5 s)
```

## Files

| Path | Content |
|---|---|
| `paper.tex`, `paper.pdf`, `references.bib` | Manuscript (statements, complete proofs, scope) |
| `references.md` | Sources and novelty check (queries and results) |
| `theory/setting.md`, `theory/proofs.md` | Agent T's working definitions and proofs (superseded by the paper) |
| `theory/scripts/check_algebra.py` | Exact (sympy) identities of the four-component and attainment proofs |
| `theory/scripts/check_ftc_achievability.py` | Attainment construction, N = 3..8, 80-digit check of every constraint |
| `theory/scripts/signcheck.py` | Sign-relaxation LP; `certify4`: the N = 4 case analysis (only family K needs more) |
| `theory/scripts/opencount.py`, `theory/data/` | Open pairs for N = 5 (sign-hard ∩ merge) |
| `theory/scripts/stressK.py`, `probeK.py` | Numerical stress tests on family K (exploration) |
| `computation/` | Agent C: enumeration, model (`nlp.py`), searches, exact certificates (`run_certificates.sh`) |
| `general/` | Agent G: general-N write-up (`general-N.md`) and test scripts |
| `certificates/sources/` | Third-party PDFs (not distributed) |
| `verification.md` | Independent verification report (verdicts, re-derivations, re-runs, corrections) |
| `verification/` | Verifier's independent scripts (`indep_checks.py`, `indep_nlp.py`, `tie_check.py`, `open_n5_probe.py`, `refine_candidate.py`, `krawczyk_certificate.py`, ...), logs and the counterexample data |
