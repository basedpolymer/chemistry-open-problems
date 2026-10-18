Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Independent verification: Banaji–Feliu, Remark 4.23

Date: 18 October 2026. Inputs used: the claimed result, the deliverables in this folder (`paper.tex`/`paper.pdf`,
`README.md`, `scripts/`, `certificates/`), `PROTOCOL.md`, and the primary source arXiv:2409.06877v5 (local copy
`travail/bf/bf.pdf`, stamped "arXiv:2409.06877v5 [q-bio.MN] 19 Aug 2026", 48 pages, re-extracted by the verifier with
`pdftotext`). I did not read `PROGRESS.md`, the solver's logs, `travail/` (other than `bf.pdf`) or `rapports/`.

## Verdict: **CONFIRMED WITH CORRECTIONS**

The main result is correct, and it answers the question as it was asked rather than a weaker variant. For every
n ≥ 2 the bound n of Corollary 4.21(i) is attained by an (n, n, n−1) network, with exactly n positive equilibria on one
stoichiometric class, all nondegenerate. The minimum source molecularity is exactly 3 for every n ≥ 2, given
Banaji–Feliu's own results on quadratic networks. The proof in `paper.tex` is complete and needs no computer. I
checked every step by hand and confirmed the conclusion with an independent computation that does not use the
solver's monomial-curve argument.

One substantive error had to be corrected. The network displayed as 𝒩₄ in `paper.tex` was wrong: it had rank 4, so it
is not a (4,4,3) network. The definitions, the table and the solver's scripts were correct; the error was only in
that displayed example. The other corrections are minor.

## 1. Statement and definitions (checked against the source)

- **Verbatim quote** (Remark 4.23, p. 30 of v5): "We leave open whether the bound on positive nondegenerate equilibria
  in Corollary 4.21(i) is achieved for general n and, if so, the minimum source molecularity of a network needed to
  achieve the bound." The quote in `paper.tex`/`README.md` is exact, and the omissions marked [...] are faithful.
- **Corollary 4.21(i)**: "(n, n, n − 1) networks admit no more than n positive nondegenerate equilibria on any
  stoichiometric class." It carries no further hypotheses. Correctly cited.
- **(n, m, r) network** (§2): n species, m reactions, rank r = rank Γ. A reaction is an ordered pair of *distinct*
  complexes, so every reaction vector is nonzero. Consequence: no (1,1,0) network exists, and n = 1 is vacuous.
- **Nondegenerate** (§2): x ∈ E is nondegenerate if J_g(x) acts as a nonsingular transformation on im Γ, i.e.
  ker J_g(x) ∩ im Γ = {0}. This is relative to the stoichiometric class, and it is exactly the definition used in
  `paper.tex`.
- **Source molecularity**: not formally defined in Banaji–Feliu. "Quadratic" means that all sources have
  molecularity ≤ 2, so the natural reading of "source molecularity of a network" is the maximum molecularity of its
  sources. `paper.tex` used this reading implicitly; a definition has now been added.
- **Corollary 5.16**: "Quadratic (n, n, n−1) networks admit no more than two positive nondegenerate equilibria on each
  stochiometric class. The smallest n for which this bound is achieved is n = 3. [...]". It applies to all quadratic
  (n,n,n−1) networks, with the same notion of nondegeneracy, so it says what the manuscript claims.
- **Case n = 2** (Remark 4.23): "quadratic (2, 2, 1) networks do not admit multiple positive nondegenerate equilibria",
  citing [5, Lemma 20]. Together with Cor. 5.16, this gives source molecularity ≥ 3 for every n ≥ 2.

**Conclusion.** Both parts of the question are answered: "achieved for general n" (yes, for all n ≥ 2) and "minimum
source molecularity" (3, uniformly in n). The construction gives *exactly* n equilibria on the class, all
nondegenerate, which is stronger than "at least n nondegenerate".

## 2. What I recomputed myself, and how

Scripts in `verification/`, written from the reaction table of `paper.tex`. I did not reuse the solver's code. Runs:
WSL Fedora, 1 thread, `ulimit -v 3 GB`.

1. **`verification/indep_groebner.py`** (155 PASS, 0 FAIL, about 4 s; log `indep_groebner.log`). This script does
   not use the monomial curve. Cases: n = 3, 4, 5, 6 with λ_r = r; n = 3, 4, 5 with non-integer λ; n = 5 with clustered
   λ (0.9, 1, 1.1, 1.2, 1.3). For each case:
   - It checks the network from scratch: nonnegative products, max source molecularity 3, rank Γ = n−1, ker Γ spanned
     by a positive vector (computed as (1,…,1)), left kernel z (computed as ε), and q = S⁻¹1.
   - It builds the rate constants with the recipe of the proof.
   - It saturates the **full** system {Γ v(x) = 0 (all n rows), z·x = T} by x₁⋯xₙ in Singular 4.4.1, which keeps
     only torus solutions and so all positive ones. It then computes a reduced lex Gröbner basis and checks that the
     basis has shape-lemma form {g(x₁), xᵢ − hᵢ(x₁)}.
   - It checks that deg g = vdim = 2ⁿ⁻¹ and that g is squarefree.
   - It counts the positive real roots of g with exact Sturm sequences and finds **exactly n**. Those n roots give n
     rational points with all coordinates > 0, which coincide with the manuscript's x(λ_r) = C∘λ_r^q. Since every
     positive equilibrium has x₁ > 0, there are exactly n positive equilibria on the class.
   - It establishes nondegeneracy of each equilibrium by two criteria different from the solver's:
     (a) e_{n−1}(J) = Σ principal (n−1)-minors ≠ 0. Since J(ℝⁿ) ⊆ im Γ, which is J-invariant of codimension 1,
     charpoly(J)(t) = t·charpoly(J|im Γ)(t), so e_{n−1}(J) = det(J|im Γ).
     (b) rank [J; zᵀ] = n.
   - For n = 3 it computes the trace and determinant exactly. At (25,36,1): tr = −961/900, det = 7/450 (stable). At
     (50,18,8): det = −1/225 (saddle). At (75,12,27): tr = −361/900, det = 1/150 (stable). The manuscript's
     determinants and eigenvalues {−1.053, −0.0148}, {−0.412, +0.0108}, {−0.384, −0.0174} are confirmed.
2. **`verification/indep_extra.py`** (31 PASS, 0 FAIL, about 3 min, 110 MB; log `indep_extra.log`).
   - (A) Structural lemmas for n = 2…80 (python-flint, exact): rank n−1, Γ1 = 0, εᵀΓ = 0, S lower bidiagonal with unit
     diagonal, q = (1 + 2(−2)^{i−1})/3 with the ordering used in the proof, source molecularity 3, nonnegative products
     of molecularity ≤ 5.
   - (B) Gröbner/Sturm count for n = 2. 𝒩₂ = {X₁ → X₂, 2X₁ + X₂ → 3X₁} is the authors' example with X ↔ Y, and it has
     exactly 2 positive nondegenerate equilibria. The script also shows that the network displayed as 𝒩₄ in the
     original `paper.tex` has rank 4 (see §4).
   - (C) The proof's construction with λ_r = (2r+1)/5, for n = 2…12. Checked: strict alternation of the kernel
     vector, C > 0, T > 0, n positive equilibria on ε·x = T, det(J|im Γ) ≠ 0 at each (signs alternate +,−,+,…), and,
     for n ≤ 7, a Sturm count of exactly n positive roots of λ^{−min q}F(λ).
3. **By hand**, line by line (§3 below):
   - Lemma 4 (stoichiometry): ε·γ_j = 0, Σγ_j = 0 coordinate by coordinate, including the boundary cases 2k = n and 2k+1 = n, and the
     bidiagonal rank argument.
   - Lemma 5 (exponents): the q formula and ordering.
   - Lemma 6 (curve): equilibria ⇔ v(x) ∈ ker Γ ∩ ℝⁿ₊ = ℝ₊·1. S is unimodular, so x = C∘λ^q.
   - Lemma 7 (nondegeneracy): ker J = span(x∘q) = span(x′(λ)), hence nondegenerate ⇔ F′(λ) ≠ 0.
   - Lemma 8 (Descartes–Laguerre) and its sketch: the derivative trick with θ between two exponents across a sign change
     removes exactly one sign change, and Rolle loses at most one zero.
   - Lemma 9: N zeros force N sign changes among N+1 coefficients, hence all nonzero, alternating, and exactly these
     zeros, all simple.
   - The sign matching in the proof of Theorem 1: in increasing exponent order the coefficient signs are
     …, ε₆ = +, ε₄ = −, ε₂ = +, [constant −T] −, ε₁ = +, ε₃ = −, ε₅ = +, …, which alternate.
   - The n = 3 instance: λF(λ) = −λ⁴ + 25λ² − 60λ + 36 = −(λ−1)(λ−2)(λ−3)(λ+6), with κ_j = C^{−a_j} giving
     (1/25, 1/22500, 1/1296).

## 3. Review of the proof (`paper.tex`)

Every step is correct. Remarks:

- The reduction "positive equilibria on the class ↔ positive zeros of F" is a bijection because λ ↦ x(λ) is
  injective (q₁ = 1), and every positive equilibrium lies on the curve (Lemma 6). "Exactly n" therefore follows
  rigorously from Lemma 9, for every n, not just numerically.
- Nondegeneracy is **proved for all n**: Lemma 7 reduces it to F′(λ_r) ≠ 0, which holds because the zeros are
  simple (Descartes bound = number of prescribed zeros). It is not merely checked numerically.
- In Lemma 7, the phrase "Because J(ℝⁿ) ⊆ im Γ" is not needed for the equivalence; it is harmless.
- Remarks (ii) and (iii) are correct. Bounded classes would need a positive vector in ker Γᵀ = span(ε), which is
  impossible for n ≥ 3, consistent with Cor. 4.21(ii). For m = n, Thm 4.20(i) is Cor. 4.21(i).
- The Descartes–Laguerre lemma is given as a sketch with references (Pólya–Szegő, Jameson). It is classical and
  standard, and the sketch is correct.

## 4. Re-run of the solver's scripts

- `scripts/certify_family.py 3 12`: exit 0, output identical to `certificates/certify_family_n3-12.log`. There are
  **300** checks (15 + 2n per n), not 301 as announced: the 301st "PASS" counted is the final "ALL PASS" line.
- What it actually checks: the network structure, the q formula, the alternation of the generalized-Vandermonde
  kernel, that the n points x(r) satisfy Γv = 0 and ε·x = T exactly, and the restricted Jacobian determinant.
- Its check "no other positive equilibrium" counts the sign changes of the same kernel vector. It is a restatement of
  the alternation check. The completeness claim therefore rests on the Lemma 6 + Descartes argument, not on an
  independent count. This is legitimate as a certificate of the construction, but it is not an independent count.
  My Gröbner/Sturm recount provides one for n ≤ 6.
- `scripts/example_n3.py`: exit 0, all PASS. It uses `sympy.solve` with `positive=True` symbols, which is not by itself
  a rigorous proof of completeness; it is superseded by `indep_groebner.py`.

## 5. Corrections applied

1. **`paper.tex`, displayed examples (substantive).** The original showed
   "𝒩₄: X₁ → X₄, 2X₁+X₂ → 3X₁, 2X₂+X₃ → 3X₂+2X₃, 2X₃+X₄ → X₃+2X₄". With s = −ε₄ = +1, c₁ = +1, d₁ = −1, c₂ = −1, the
   correct 𝒩₄ is "X₁ → 2X₁+X₄, 2X₁+X₂ → X₁+2X₂, 2X₂+X₃ → X₂, 2X₃+X₄ → 3X₃", which is what the reaction table (Table 3 of the PDF) and the solver's
   script produce. The displayed network has a stoichiometric matrix of **rank 4** (no positive flux vector), so it is
   not even a (4,4,3) network. A reader building it would have found no equilibria at all. I corrected it, added 𝒩₅
   (checked) and a sentence warning that the signs matter.
2. `paper.tex` and `README.md`: "301 checks" corrected to "300 checks".
3. Abstract: "Together with their bound of two equilibria for quadratic networks" did not cover n = 2. It now refers
   to their results on quadratic networks (at most two, and at most one for n = 2).
4. Notation: added the definition of source molecularity as the maximum molecularity of a source, consistent with
   Banaji–Feliu's definition of "quadratic".
5. Proof of Theorem 1: added that the formulas also give 𝒩₂ (the authors' example up to renaming), to which the
   argument applies verbatim. Proof of Corollary 2: added that n = 1 is vacuous.
6. "Computer verification" section: one sentence on the independent Gröbner recount.
7. Credits and status: `\PaperCredits` and `\PaperSubtitle` updated, `README.md` status/attribution/reproduction/file
   table updated.
8. `paper.pdf` recompiled with `templates/build_tex.py`: 6 pages, no overfull lines.

The root `README.md` results table, which is outside this folder, has not been modified. It should be updated by
the orchestrator.

## 6. What I only read (not recomputed)

- Banaji–Feliu's upper bounds, used for minimality: Cor. 5.16 (quadratic (n,n,n−1): ≤ 2) and the n = 2 claim in
  Remark 4.23 (quadratic (2,2,1): ≤ 1, via [5, Lemma 20]). These are published results (J. Math. Biol. 2026) and I
  did not re-prove them. The minimality part of the answer depends on them.
- The classical Descartes rule for real exponents: sketch checked, references not consulted.
- Novelty:
  - Semantic Scholar lists 2 works citing arXiv:2409.06877 (2412.17798, 2411.15134, both by Feliu et al., on
    steady-state varieties and toric invariance).
  - OpenAlex lists 0 citations of doi:10.1007/s00285-026-02429-8.
  - Banaji's 2025–2026 arXiv papers (2601.13730 on epidemics; 2609.07265 on fold/cusp bifurcations, abstract only)
    do not address the question.
  - No repository of GitHub user `mmislan` (Palomar) concerns it.
  - Two web searches found nothing.
  - I did not read the two citing papers in full.

## 7. Reservations

- Minimality (= 3) is conditional on Banaji–Feliu's Cor. 5.16 and on their n = 2 remark. This is the reading
  intended by the question, but it is not re-proved here.
- "Source molecularity" is read as the maximum source molecularity. This is the only reading consistent with the
  authors' use of "quadratic". Under another reading (e.g. requiring *all* sources to have molecularity ≥ 3), 𝒩ₙ,
  whose source X₁ has molecularity 1, would need modification; I see no reason to adopt such a reading.
- The novelty check is limited to citation databases, arXiv author listings, GitHub and web search. Work that is
  unpublished or not yet indexed cannot be excluded.
- The witnesses grow quickly: exponents ~2ⁿ, and rational data with thousands of digits for n = 12. This does not
  affect the existence proof.
