Verifier: Claude Opus 5.5 (effort not recorded), fresh context, no access to the solver's reasoning.

# Verification: the exact minimum number of reactant complexes for ACR with bistability

Date: 31 October 2026. Solver of the new result: GPT-6-Astra (effort xhigh).

## Verdict: **CONFIRMED WITH CORRECTIONS**

The minimum is exactly five. Both halves hold: no network with at most four reactant complexes has ACR with two
exponentially stable positive equilibria in one class, and the explicit five-source network does.

No theorem statement, table, count, normal form, script or certificate value had to change. The corrections are one
gap in the exposition (a formula used but not derived), one load-bearing convention that was left unstated, a more
exact citation, and status and wording. They are listed in section 8. The author's own notes are kept in
[self-check.md](self-check.md).

The strongest evidence, beyond the agreement of two codes, is structural (section 4): the exclusion that removes
the dangerous cases has an elementary explanation, and it is exact on all 371 spans concerned.

## 1. Inputs and order of work

I read `PROTOCOL.md`, then `paper.tex`, `README.md`, `references.md` and the parent manuscript, and audited the proof
by hand. I then wrote `verification/enum_sources.py` and `enum_spans.py` from the mathematical specification,
**before opening the solver's scripts**, and only afterwards read and re-ran those. The primary source was read
in full text (arXiv:2307.04186v1, PDF, with pdftotext).

Limit of this separation: my enumeration implements the logic stated in the manuscript, so two matching codes do not
by themselves prove that logic. That is what the hand audit in section 2 and the structural checks in section 4
are for. The methods do differ: exact RREF keys instead of Pluecker coordinates, orbits by union-find on the
labelled family instead of a canonical form, a coefficient tensor instead of a symbolic polynomial, and an LP as a
second opinion on the cone test.

## 2. Hand audit of the lower bound

**Definitions.** KNS Definition 2.6 is the ACR used here, global across classes and requiring a positive steady
state. Complexes have the form 0, X_i, X_i+X_j or 2X_i, and each species appears in some complex. The zero complex
counts as a reactant: Example 4.13 of KNS gives {2Z→Z, X+Y→Z→Y+Z, 0→X} as having "4 reactant complexes". The result
depends on this convention, since the example has only four nonzero reactants. The manuscript did not state it
(correction C2). The Section 6 sentence is quoted verbatim. The exact-minimum question is a refinement, as the
manuscript says; it is not a conjecture of the source.

**(a) Toric reduction.** Write f = N m(x) with one column per source. At a nondegenerate positive equilibrium
im J ⊆ im N ⊆ S and J|_S is invertible, so im N = S. With r ≤ 4 sources and m(x*) > 0 in ker N this gives dim S ≤ r−1 ≤ 3,
which handles fewer than four sources: r ≤ 3 gives dim S ≤ 2. N is a constant matrix, so rank N = 3 holds or fails
for the whole system. If rank N = 3 with four columns, ker N is the line through m(x*), any three columns are
independent, and all positive equilibria satisfy A log x = A log x*, so E = x* ∘ exp(ker A). ACR is then exactly
e_i ∈ row(A). None of this uses genericity of κ, or how many reactions leave each source.

Checked case by case:
- *Arbitrary reactions per source, non-generic κ:* only the aggregated columns enter. p, q, T, U stay free parameters
  in every normal form.
- *A column N_a = 0, or cancelling reactions:* excluded when rank N = 3, since all triples of columns are then independent.
  If the columns have rank 2 then dim S = 2 at a nondegenerate equilibrium, and the parent theorem applies, because it
  covers any number of sources.
- *Inert species:* a species absent from all sources has nonnegative derivative, so at a positive equilibrium no
  reaction produces it. It appears in no complex at all, and removal changes nothing. It cannot be the ACR species.
- *Species count:* n ≥ 3 from rank A = 3. For n = 3 the set E is a point. n ≥ 9 is impossible after removing inert
  species, and n = 7, 8 give no type.
- *Products:* molecularity of products enters only through the list of allowed reaction vectors.

**(b) Stability relative to the class.** J = N₀ diag(m₁,m₂,m₃) A diag(1/x) at equilibrium, so in the basis N₀ of S the
matrix of J|_S is m₄ diag(p) A diag(1/x) N₀, with determinant m₄³p₁p₂p₃ det(A diag(1/x) N₀). This gives formula (6) of the
manuscript, which I re-derived with its constant. A real 3×3 Hurwitz matrix has negative determinant, whether or not it has a complex
pair. So two exponentially stable equilibria in one system have the same determinant sign, and the argument needs
only sign comparisons within one fixed (S, κ). If rank N < dim S, J|_S is singular and nothing is stable. This is
stability in the stoichiometric class, as claimed; for a system with a nondegenerate equilibrium the kinetic and
stoichiometric subspaces coincide anyway.

**(c) At most two nondegenerate equilibria; connectedness.**
- *n = 4:* three exponentials have at most two real zeros with multiplicity, by Rolle. With G the class function along
  the curve, the augmented determinant is a fixed constant times G'(t)/∏x_i, which I derived directly.
- *S4, S5, S6, S8:* I checked each monomial reduction by hand: S4 {a},{b},{d,e}; S5 {a},{b},{e}; S6 {a},{b},{c,e};
  S8 {a,b},{c},{e}. At equilibrium det(W′ diag(x) V′) is a fixed constant times g₁g₂ det[Q; ∇F], so rank Q < 2 gives
  degeneracy, and on the line the sign is a fixed sign times that of F′.
- *The 19 forms:* each scalar reduction was re-derived. Every parameter domain is an interval: it is a line or a
  superlevel set of concave functions meeting a half-line, and the scalar curve covers the whole positive solution set
  of the remaining equations. For S3 forms 1–6 and S7 forms 1–2 strict log-concavity was checked from the signs of
  α, β, γ, δ. For S7 forms 1–2 the substitution u = 1/c is needed: in c itself the function is not concave.
- *Orientation:* `check_forms.py` works with the quantity that matters, D = det(W′ diag(x) V′), not with binomial
  rows. Along each curve D = μ·L′ identically, L being the logarithm of the remaining binomial. I factored each μ by
  hand into a fixed sign times positive quantities:

| Forms | μ | | Forms | μ |
|---|---|---|---|---|
| S3/0 | + b c d | | S7/0 | + b c e |
| S3/1–6 | − a b c | | S7/1–4 | − a b c |
| S3/7, 8 | − a c d | | S9/0 | + b c d f |
| S3/9 | − c d (a+b)/2 | | S9/1 | − a b d f |
| S3/10 | − c d (a+b) | | S9/2 | + a b c d |

**(d) Degenerate equilibria.** Theorem 2 is a statement about nondegenerate equilibria only, and exponentially stable
ones are nondegenerate, so degenerate equilibria cannot help. In each scalar case the count is by multiplicity or by strict
concavity, so a degenerate root leaves no room for two simple ones, and an identically vanishing reduction makes every
point degenerate. The injective case gives at most one equilibrium of any kind, by the logarithmic mean argument, which
I re-derived. The identically vanishing case is not vacuous: 260 of the systems sampled in section 5 have a
positive-dimensional solution set.

**The cited parent theorem.** The proof uses exactly its statement: a two-dimensional at-most-bimolecular system with
ACR has at most one exponentially stable positive steady state per class. It has no restriction on species or sources
and the same definitions of ACR and stability. I re-derived the proof lemma by lemma (the field on the ACR
hyperplane, α₁ = 0, the normal form and determinant, connectedness of V⁺, global ACR, tangency, the final sign
argument) and found no gap. The parent's source list is the S-table here. KNS Proposition 3.10 covers
dimension one as cited, and dimension one also follows from a two-line argument with no ACR hypothesis (added to the proof).

**Lemma 8 and the two exclusions.** The logic of each is right: maximal network, Gordan–Stiemke test by extreme rays,
Cauchy–Binet expansion, and the kinetic cubic, where φ_r = κ_r m_a(x*) gives L_a = m_a(x*) N_a. The one gap: the
manuscript used det(W diag(x) V) to decide degeneracy without connecting it to the Jacobian (correction C1).

## 3. Reimplementation of the enumeration

All runs reproduce the manuscript exactly.

| Claim | Manuscript | Mine |
|---|---|---|
| Source types, n = 4, 5, 6, 7, 8 | 26, 8, 1, 0, 0 | 26, 8, 1, 0, 0 |
| S1–S9 sources and ACR species | table | identical, one orbit each |
| n = 4: distinct kernel entries | ≤ 3, incl. 0 | ≤ 3, on all 358 labelled sets |

| | S1 | S2 | S3 | S7 | S9 |
|---|---|---|---|---|---|
| All rank-three spans | 4333 | 4287 | 4333 | 4397 | 23035 |
| Missing a source | 2965 | 3161 | 2965 | 3050 | 17734 |
| No positive maximal flux | 1126 | 782 | 1126 | 1042 | 4442 |
| Injective | 169 | 301 | 129 | 231 | 569 |
| Class-degenerate | 34 | 6 | 6 | 32 | 144 |
| Kinetically degenerate | 39 | 37 | 43 | 20 | 92 |
| Remaining | 0 | 0 | 64 | 22 | 54 |
| Normal forms | | | 11 | 5 | 3 |

Total 40,385 spans, 140 remaining. Beyond equal totals:
- The label of **every one of the 1,992 spans** reaching the injectivity stage is identical in the two codes, compared
  by exact RREF of S^⊥ (`compare_solver.py`). The tables of forms in the manuscript, typed in by hand, equal my computed
  orbits (`check_forms.py`). I used the full automorphism group of each source set, found by brute force.
- The exact cone test and the HiGHS LP agree on all spans that reach it.
- S1 and S3 having equal first three rows is not a coincidence or a copying error: exchanging E with the slack variable
  2 − |y| permutes the complexes and swaps the two source sets. It acts linearly on reaction vectors, which fixes those
  three rows, but it does not preserve mass-action kinetics, so the later rows differ.

My classification also covers S4, S5, S6, S8 and the 26 four-species types, which the proof does not need.

## 4. The exclusions the theorem depends on

The control in section 5 shows why this deserves attention: spans removed by the kinetic test **do** have classes with
four positive solutions of the toric system. If that exclusion were wrong, the theorem would very likely fall.

- **An elementary proof of the kinetic exclusion, without the cubic** (`kinetic_structure.py`, exact). Each column N_a
  lies in V_a, the span of the allowed vectors of source a in S, and any three columns must be independent. By Rado's
  theorem that needs dim(Σ_{j∈J} V_j) ≥ |J| for every J inside a triple. **All 231 kinetically degenerate spans fail this
  condition, and all 140 surviving spans satisfy it.** So on these spans the cubic test coincides exactly with an
  elementary criterion, and neither over- nor under-excludes. Example: for S3 with conservation of a − b and c + d,
  sources E, 2E, A+B all act in a plane, so their columns are dependent.
- **The positive-flux exclusion is not needed** (`flux_redundancy.py`). The one mutant of section 6 that survives
  disables this test. I measured why: all 8,518 spans it rejects are also labelled injective, class-degenerate or
  kinetically degenerate, and none gives a new form. Lemma 8 is correct but not load-bearing, which removes one
  possible failure point.

## 5. Tests on actual systems (floating point or sampling: tests, not proofs)

- `stress_count.py`: the full polynomial system with no scalar reduction, solved by msolve (certified real root
  isolation), on 100 random classes for each of the 253 surviving spans of S3, S4, S6, S7, S8, S9. Result: never more than
  two positive solutions, opposite signs of D whenever two, and the seed point always recovered. **Control:** the same
  harness on 155 kinetically degenerate mixed-sign spans finds 236 classes with **four** positive solutions, so the
  harness can see a violation.
- `stress_networks.py`: 797 random networks built from reaction lists, on 403 surviving spans including all four-species
  types, with an equilibrium planted through a positive flux. Result on 2,380 classes: no class with more than two equilibria or
  with two stable ones. Each solution is a true zero of that network's mass-action field, classified from the eigenvalues of J|_S.
  A blind Newton multi-start gave 4,075 interior zeros, and the ACR species never moved, supporting the toric description.
  A first run reported 28 "ACR moved" events; all had coordinates of 1e-13 or less, usually exactly 0. They are boundary
  limits reached as log x → −∞, not positive equilibria, and are now counted apart.
- Kinetically degenerate spans: rank of the aggregated matrix below 3 under all 1,300 sampled fluxes.

## 6. The solver's certificates

`scripts/run_all.py` re-run: all 148 checks PASS, 19 s on 3 cores. `sources.json`, `subspaces.json`, `normal-forms.json` and
`five-example.json` came out byte-identical (SHA-256) to those supplied. I read every script: they do what the manuscript
says. The membership adjugate is correct, sources stay separate in the kinetic test, the bounds claimed for integer
overflow are valid (crude but true), and no numerical file is read.

**Corrupted inputs** (`mutate_solver.py`). The `.json` certificates are outputs that no certificate script reads, so
corrupting them proves nothing; the real inputs are the data and tests inside the scripts. On temporary copies I applied 24
corruptions: dropped, altered and spurious normal forms, wrong sources, weakened tests, restricted products, altered expected
values, wrong curves and determinant factors, changed rates, a trimolecular product, a wrong equilibrium and trace, a saddle
called stable. **All 24 were caught.** A 25th, disabling the dual-cone test, survives as expected (section 4). A flipped
label in a copy of `subspaces.json` is caught by my comparator.

## 7. The five-source network; novelty

`check_example.py` rebuilds everything from the reaction list with its own mass-action code.
- Five reactant complexes 0, A+B, A+C, B+C, R; all complexes at most bimolecular; dim S = 3; S^⊥ spanned by (1,1,0,0).
- ṙ = 1 − r exactly, so r = 1 at every positive equilibrium of every class: global ACR. A positive equilibrium exists
  in every class a + b = T: the scalar function is 64T at a = 0 and −5T at a = T.
- Class a + b = 9: exactly (4,5,1/5,1), (6,3,1/3,1), (8,1,1,1), by a lex Gröbner basis and again by msolve.
- Stability from the restricted Jacobian in a basis of S computed from the reaction vectors, not the manuscript's
  coordinates, by exact Routh–Hurwitz. Characteristic coefficients (a₂,a₁,a₀) are (1689/5, 9364/5, 1536) and (113, 1648, 1536)
  at the stable points, with exact real roots ≈ −332.18, −4.62, −1 and −96, −16, −1. The middle point has (207, −562, −768)
  and a positive eigenvalue ≈ 3.66.
- The coupled variant was checked the same way. 200 of 200 rate vectors within 2% keep the picture (floating point).

The regulator R is decoupled in the main example. That is allowed by the definition, and the coupled variant answers any
objection. Minimality is stated correctly: five is the minimum for reactant complexes alone, and the manuscript rightly makes
no simultaneous claim about species or reactions.

**Novelty.** No statement of the exact minimum found. Screened: the arXiv listing for the phrase, 11 citing works on Semantic
Scholar, 4 on OpenAlex, Zenodo, the `ipitchford` and `mmislan` repositories, GitHub and web searches. The closest recent
works (arXiv:2604.14592, 2608.18509, 2605.11389) address zero-one networks or one specific family. Details in `references.md`.
This is a bounded search of titles and abstracts, not a proof of priority.

## 8. Corrections applied (all by the verifier)

| # | Severity | File | Correction |
|---|---|---|---|
| C1 | Minor gap | `paper.tex` | Section 2 lacked the link between det(W diag(x) V) and the Jacobian, yet Section 4 decides degeneracy with it. Added equation (7) with its proof, and cited it in Section 4. Checked in exact arithmetic on a sample (c = det(AC)/det[V C]). |
| C2 | Unstated convention | `paper.tex`, `README.md` | The zero complex counts as a reactant complex; now stated, with KNS Example 4.13. |
| C3 | Citation | `paper.tex` | Corollary 3: exact pointers (parent's dimension-two theorem, KNS Proposition 3.10) and a direct argument for dimension one. |
| C4 | Factual, trivial | `paper.tex`, `README.md` | Runtime "about one minute" changed to under one minute (19 s measured). |
| C5 | Addition | `paper.tex` | Short paragraph recording the two facts of section 4. |
| C6 | Status | `paper.tex`, `README.md`, `references.md` | "Preliminary" and "verification pending" removed, verifier credited, abstract and closing paragraph updated, "solved" per PROTOCOL.md section 5. |
| C7 | Wording | `paper.tex` | "independently generates" changed to "generates from scratch". |
| C8 | Files | folder | `verification.md` renamed `self-check.md` with its content unchanged; README link fixed; this report and `verification/` added. |

`paper.pdf` was rebuilt (10 pages). `certificates/*.log` and `execution.json` were refreshed by re-running `run_all.py`,
because that file hashes `paper.tex`. Harmless points left as they are: the header of `certify_subspaces.py` says "a few
minutes" (it takes 17 s), and the field `independent_verification: false` in `execution.json` is hard-coded by the runner
and describes that runner, not the status of the result.

## 9. What was not checked

- The parent's scripts were not re-run; its theorem was re-derived by reading the proof, and rests on that and on its own
  verifier. No separate stress test of dimension two was done.
- `scripts/novelty_audit.py` was not re-run; I ran my own queries. Citing works were screened by title and abstract.
- Exponential sums (n = 4) and the group argument (S4, S5, S6, S8) were checked by hand and by sampling on spans;
  there is no symbolic certificate of them, and none is needed for arguments this short.
- No formal proof assistant was used. No human has reviewed the result.
- Larger products and non-hyperbolic stability are outside the claim. The root `README.md` was left alone, being outside
  this folder.

## 10. Reproduce

Under WSL, from this folder: `python verification/run_verification.py` (about 12 minutes on 3 cores). It runs every script below,
prints one PASS/FAIL per step and writes full logs to `verification/logs/`.

| Script | Purpose |
|---|---|
| `enum_sources.py`, `enum_spans.py` | Reimplementation, written before reading the solver's code |
| `compare_solver.py` | Span-by-span comparison with `certificates/subspaces.json` |
| `check_forms.py` | Tables of forms, D = μ L′, concavity |
| `kinetic_structure.py`, `flux_redundancy.py` | The structural facts of section 4 |
| `mutate_solver.py` | The 25 corruptions |
| `check_example.py` | The five-source network and its variant |
| `stress_count.py`, `stress_networks.py` | Sampling tests with controls |
