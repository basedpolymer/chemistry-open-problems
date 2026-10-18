# References and novelty check

| Source | Read | Use |
|---|---|---|
| M. Banaji, E. Feliu, J. Math. Biol. 93 (2026), doi:10.1007/s00285-026-02429-8, arXiv:2409.06877v5 | full text | Theorem 5.5, Remarks 5.4, 5.6, 5.8, definitions |
| M. Banaji, B. Boros, J. Hofbauer, *Bifurcations in planar, quadratic mass-action networks with few reactions and low molecularity*, Nonlinear Dyn. 112 (2024) 21425–21448, doi:10.1007/s11071-024-10068-1, arXiv:2406.13451 | Lemma 20 and Corollary 19 (text) | Correct statement behind Remark 5.6 |
| K. Phillipson, J. M. Rojas, *Fewnomial systems with many roots, and an adelic tau conjecture*, Contemp. Math. 605 (2013) 45–71, doi:10.1090/conm/605/12111, arXiv:1011.4128v5 | Theorem 1.6, Lemma 1.8 (p. 3) | Same polynomial system as the equilibrium equations of 𝒩ₙ (added by the verifier) |
| F. Bihan, *Polynomial systems supported on circuits and dessins d'enfants*, J. London Math. Soc. 75 (2007) 116–132, doi:10.1112/jlms/jdl013 | abstract (arXiv:math/0509219) | Sharp bound n + 1 for circuits (added by the verifier) |
| Anonymous (maintained by I. Pitchford, Evidence Press), *Sharp positive-equilibrium counts in quadratic mass-action networks with two more reactions than species*, v1.0.0-candidate, unrefereed, 4 Oct 2026, doi:10.5281/zenodo.23146115, github.com/ipitchford/quadratic-mass-action-sharpness | README, PRIOR_ART.md, LaTeX source (as data); explicit examples re-checked | Earlier answer to Remark 5.8 for all n ≥ 2 (added by the verifier) |
| M. Banaji, *On the computation of fold and cusp bifurcations in mass action networks via dual equations*, arXiv:2609.07265v1 (7 Sep 2026) | text search | Mentions only the three-equilibrium (2,4,2) example; does not address Remark 5.8 |

## Novelty check

**Solver (18 Oct 2026).** Remark 5.8 is unchanged in v5 (Aug 2026) and in the published version. Semantic Scholar lists
two citing works of the preprint (arXiv:2412.17798, arXiv:2411.15134, both 2024) and OpenAlex none of the journal
version; neither addresses (n, n+2, n) networks. No Palomar (mmislan) repository on this question.

**Verifier (18 Oct 2026).** These citation-database checks are confirmed, but a web search on the problem's key
phrases finds the unrefereed Zenodo/GitHub manuscript above, published on 4 Oct 2026 (Zenodo API: status
"published"). It answers Remark 5.8 affirmatively for every n ≥ 2. Its explicit (3,5,3) example (products of
molecularity ≤ 10) and its all-n family for n = 2..6 were verified exactly with the verifier's own code. Moreover, the
equilibrium system of 𝒩ₙ is the Phillipson–Rojas system G_ε (R = 16 for ε = 1/4; the R = 64 networks of the proof
have the same Gale data as G_{1/8}). The answer to Remark 5.8 is therefore **not new**. The status is "already solved",
with an independent, verified proof here. The Remark 5.6 discrepancy is also noted in the earlier manuscript's
prior-art notes. The mmislan repositories (23, checked 18 Oct 2026) contain nothing on this question.
