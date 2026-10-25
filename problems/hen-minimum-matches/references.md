# References and novelty check

## Primary source (read in full: arXiv v2 PDF, 11 Apr 2018)

- D. Letsios, G. Kouyialis, R. Misener, *Heuristics with performance guarantees for the minimum number of matches
  problem in heat recovery network design*, Comput. Chem. Eng. 113 (2018) 57–85, doi:10.1016/j.compchemeng.2018.03.002,
  arXiv:1709.04688 (v2). Read: Sections 2 (transshipment / transportation MILPs), 9 (numerical results), Appendix G
  (Tables G.5, G.9, G.11). Local copy (not committed): `certificates/sources/letsios2018/`.
  - Open statement, Sec. 9.2, p. 30 of arXiv v2 (verbatim): *"Observe that state-of-the-art approaches cannot, in
    general, solve moderately-sized problems with 30-40 streams to global optimality. For example, none of the test
    cases in the Grossmann (2017) or large scale test sets is solved to global optimality within the specified time
    limit."*
  - Solver settings, Sec. 9.2 (verbatim): *"In each solver run, we set absolute gap 0.99, relative gap 4%, and
    maximum number of threads 1."* (The repository code sets only `mipgap = 0.04`; the CPLEX logs show
    "tolerance (0.04/1e-06)", i.e. absolute gap 1e-6.)  Consequently an instance reported as solved may stop with a
    one-unit gap: unbalanced10 (25, bound 24) and 37sp-yfyv (36, bound 35) are such cases (logs re-read here).
  - Sec. 9.4: *"Although CPLEX and Gurobi do not converge to global optimality for many of the Furman (2000),
    Chen et al. (2015a,b), and Grossmann (2017) instances, ..."*
  - Table G.9 also quotes Furman–Sahinidis (2004) CPLEX 7.0 runs (7 h limit, timeouts) for 20sp1 (19), 22sp1 (25),
    23sp1 (23), 37sp-yfyv (36).
- Data, code and logs: github.com/dimletsios/min_matches_heuristics (identical data in
  github.com/cog-imperial/min_matches_heuristics), last commit 2019-02-15; cloned into
  `certificates/sources/min_matches_heuristics/` (git-ignored). Every bound in our README table marked "LKM18" was
  re-derived by us from the `.sol`/`.log` files of `data/mip_solutions/` (script in PROGRESS notes; summary in
  README).

## Other sources used

- K. C. Furman, N. V. Sahinidis, *Approximation algorithms for the minimum number of matches problem in heat
  exchanger network synthesis*, Ind. Eng. Chem. Res. 43 (2004) 3554–3565, doi:10.1021/ie049978l. **Read in full**
  (publisher PDF, 25 Oct 2026; local copy in `certificates/sources/`, not committed). Table 3 gives
  CPLEX 7.0 values 20sp1 = 19\*, 22sp1 = 25\*, 23sp1 = 23\*, 37sp-yfyy = 36\*, the asterisk meaning (footnote a,
  p. 3562) that "an exact algorithm did not converge for these problems after 7 h of CPU time"; the text (p. 3561)
  says that the lower bounding stalled after 19, 41, 13 and 12 CPU h for these four problems. No lower bounds are
  reported. So none of 20sp1, 23sp1, 37sp-yfyv was proved optimal there, and 22 < 23 for 23sp1 is a new upper bound.
  (Their Table 1 lists 37sp-yfyy with 22 hot and 17 cold streams; the data file of the 2018 repository has n = 21,
  m = 17, utilities included; we use the 2018 data file throughout.)
- S. A. Papoulias, I. E. Grossmann, *A structural optimization approach in process synthesis—II: Heat recovery
  networks*, Comput. Chem. Eng. 7 (1983) 707–721, doi:10.1016/0098-1354(83)85023-6 (transshipment model). Not read;
  the model is taken from Letsios et al. Sec. 2.
- J. Cerdá, A. W. Westerberg, *Synthesizing heat exchanger networks having restricted stream/stream matches using
  transportation problem formulations*, Chem. Eng. Sci. 38 (1983) 1723–1740, doi:10.1016/0009-2509(83)85029-5 (transportation model). Not read.
- Y. Chen, I. E. Grossmann, D. C. Miller, *Computational strategies for large-scale MILP transshipment models for
  heat exchanger network synthesis*, Comput. Chem. Eng. 82 (2015) 68–83, doi:10.1016/j.compchemeng.2015.05.015
  (origin of the balanced/unbalanced instances). Not read.
- The "N − 1" (Euler) rule for the minimum number of units, B. Linnhoff, D. R. Mason, I. Wardle, *Understanding heat
  exchanger networks*, Comput. Chem. Eng. 3 (1979) 295–302, doi:10.1016/0098-1354(79)80049-6, and E. C. Hohmann,
  PhD thesis, USC (1971): U_min = N − s, s = number of independent subsystems. Not read; our Lemma 1 (result.md) is a
  self-contained proof of the version we use (components must be balanced *and cascade feasible*; s is computed
  exactly by enumeration).
- K. K. H. Cheung, A. Gleixner, D. E. Steffy, *Verifying integer programming results*, IPCO 2017, LNCS 10328,
  148–160, doi:10.1007/978-3-319-59250-3_13 (VIPR format); L. Eifler, A. Gleixner, *Safe and verified Gomory mixed
  integer cuts in a rational MIP framework*, SIAM J. Optim. 34 (2024) 742–763, doi:10.1137/23M156046X (VIPR 1.1,
  weak derivations). Read: VIPR README and `cert_spec_v1_1.md` only.
- SCIP Optimization Suite 10.0 (SCIP 10.0.2, SoPlex 8.0.2) source from github.com/scipopt/scip (tag v10.0.2,
  git b8eaf98) and github.com/scipopt/soplex (tag v8.0.2); VIPR tools from github.com/scipopt/vipr (commit
  30f2951, 29 Oct 2025).

## Novelty check (24 Oct 2026)

- Citing works of Letsios et al. 2018 (Semantic Scholar API, 11 citing works; OpenAlex reports 13 but its shared
  daily budget was exhausted, so its list could not be retrieved): Zhang et al. CES 2024 (doi:10.1016/j.ces.2024.119763, not read),
  Cerdá–Cafaro–Cafaro IECR 2022 (doi:10.1021/acs.iecr.1c02372, **read in full** on 25 Oct 2026: a different problem, the detailed structure of HENs with substreams and stream splitting at minimum
  utility usage and minimum number of units, or minimum total cost; 45 literature examples with 4–21 process
  streams, codes xSP.w; the main text reports total costs for 24 examples (Table 1) and no minimum-number-of-matches
  value for any data file of the 2018 benchmark; per-example unit counts are only in its Supporting Information,
  not read), Chang et al. AIChE J 2020 Parts I/II (globally optimal HENS,
  different problem), Kouyialis–Wang–Misener Processes 2019 (symmetry detection), Letsios et al. CACE 2020
  (arXiv:1909.12328, survey; read the min-matches section: it calls the problem "a major bottleneck" and reports no
  new values), Goderbauer et al. CACE 2019 (complexity), Leitold et al. 2019/2020 (HEN structural analysis), Nair–
  Karimi IECR 2019 (stageless superstructure). None reports new bounds on the benchmark instances (titles/abstracts
  checked; Cerdá et al. 2022 full text not accessible, listed in NEEDS.md).
- Kouyialis & Misener, *Detecting symmetry in designing heat exchanger networks*, arXiv:1711.02564 (read): symmetry
  of the transshipment model, one small test case; no benchmark bounds.
- Misener group page wp.doc.ic.ac.uk/rmisener/?p=1313: lists the 2018 paper, the 2020 survey and the code; no
  newer results.
- AI campaigns: the GitHub user-listing API is blocked in this sandbox for users `mmislan` and `ipitchford`
  ("sessions are bound to their configured repositories"), so these could not be checked directly; Zenodo searches
  ("minimum number of matches", "heat exchanger network" matches) returned no record on these benchmarks.
- Web searches (several phrasings: instance names "37sp-yfyv", "unbalanced15", "balanced12"; "minimum number of
  matches" + exact/valid inequalities/branch-and-price 2021–2026): nothing beyond the 2018 paper and its repository.
- Conclusion: to our knowledge the instance values have not been improved or certified since 2018.

## Novelty check, second pass (25 Oct 2026, VIPR closure session)

- OpenAlex citing works of Letsios et al. 2018 (`filter=cites:W2754455405`): now retrieved, 13 works: the 11 above
  plus G. Kouyialis, PhD thesis, Imperial College London (2019), doi:10.25560/67897 (read: Chapter 5 is the 2018
  study with the same tables, e.g. unbalanced10 25 solved within tolerance, 37sp-yfyv 36, 20sp1 19\*, 23sp1 23\*),
  and an unrelated 2025 paper on ammonia reactors (doi:10.1016/j.cep.2025.110498). Also D. Letsios et al.,
  ESCAPE 2018 proceedings (doi:10.1016/B978-0-444-64235-6.50101-7) and the Leitold et al. book chapter (2020)
  (titles only).
- GitHub: `gh search repos` on "minimum number of matches", "min_matches_heuristics", "heat exchanger network
  matches", "heat recovery network MILP", "transshipment heat exchanger", "unbalanced10", "37sp-yfyv": only
  dimletsios/min_matches_heuristics and cog-imperial/min_matches_heuristics (no push since Feb 2019) and
  lechugaa/HENS (2019–2020 student replication of the big-M study, no bound tables).
- AI campaigns: `gh api users/ipitchford/repos --paginate` (193 repositories) and `users/mmislan/repos`: no
  heat-exchanger-network or minimum-matches entry.
- Web searches (several phrasings, 2021–2026): nothing beyond the 2018 paper, its repositories and the thesis.
- Y. Chen, I. E. Grossmann, D. C. Miller (2015): the CMU KiltHub preprint (figshare article 6466850) could not be
  fetched (signed download link); only its abstract was read ("find the exact global optimal solution with reduced
  solution times"). Whether it reports a proven optimum for unbalanced10 is therefore not checked (listed in
  NEEDS.md); Letsios et al. 2018 quote no such value.
- Gemini (agy) was not used (quota exhausted).
- Conclusion: to our knowledge none of the optima 37sp-yfyv = 36, 20sp1 = 19, 23sp1 = 22, unbalanced10 = 25 had
  been proved before.

## Software read for the certification (25 Oct 2026)

- VIPR `code/viprchk.cpp` at commit 30f2951 and at master 99eec19 (24 Oct 2026): read the derivation checks
  (`readMultipliers`, `readLinComb`, `Constraint::round`, `canUnsplit`, `dominates`, `_isFalsehood`); found the two
  gaps described in result.md Sec. 4 (rounding over continuous variables; unsplit with fractional right-hand side).
- VIPR specification `cert_spec_v1_1.md` (read in full) for the independent checker.
- SCIP 10.0.2 / SoPlex 8.0.2 sources: not read beyond the build options.
