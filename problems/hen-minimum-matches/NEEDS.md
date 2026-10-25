# Needs (for the repository owner)

Supplied on 25 Oct 2026 and read (thank you): Cerdá–Cafaro–Cafaro, IECR 61 (2022) 464–486, and Furman–Sahinidis,
IECR 43 (2004) 3554–3565 (findings in `references.md`). Neither reports a proven optimum for any instance closed here.

Still missing; none is needed for the correctness of the results, they would only complete the novelty check:

- doi:10.1016/j.ces.2024.119763 — Zhang et al., Chem. Eng. Sci. (2024), cites Letsios et al. 2018; the title
  suggests multi-objective HEN design; check for minimum-number-of-matches values on the benchmark instances.
- doi:10.1016/j.compchemeng.2015.05.015 — Chen, Grossmann, Miller, Comput. Chem. Eng. 82 (2015) 68–83, origin of the
  balanced/unbalanced instances; check whether it reports a *proven* optimum for unbalanced10 (or unbalanced15,
  balanced12, balanced15). The CMU preprint (kilthub.cmu.edu, figshare article 6466850) could not be fetched here.
- Optional: Supporting Information of Cerdá–Cafaro–Cafaro 2022 (per-example unit counts; a different, detailed
  minimum-units problem, so probably not comparable).

Suggestion (no action taken; external action needs your approval): the two soundness gaps found in VIPR's
reference checker `viprchk` (rounding over continuous variables, unsplit with a fractional right-hand side; test
certificates in `certificates/vipr/viprchk_tests/`, patch `certificates/viprchk_harden.patch`) could be reported
upstream to github.com/scipopt/vipr. The SCIP 10.0.2 certificate defects (Gomory-cut derivations that cannot be
verified; conflict rows marked as derivable with global bounds) could be reported to github.com/scipopt/scip.
