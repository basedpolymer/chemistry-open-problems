Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning or logs.

# Verification: ACR with multistability (Kaihnsa–Nguyen–Shiu conjecture)

Date: 27 October 2026.

**Inputs used.**
- The primary source, downloaded by me: arXiv:2307.04186v1 (Kaihnsa, Nguyen, Shiu), read with pdftotext. I used the
  §6 conjecture, Definition 2.6, Lemmas 2.3–2.5, Propositions 3.1, 3.5, 3.9, 3.10, Theorems 1.2, 4.1, 4.9, 5.1, 5.5
  and Example 4.10.
- The deliverables in this folder: `README.md`, `paper.tex` / `paper.pdf`, `scripts/`, `references.md`,
  `references.bib`, and the two sub-study reports `examples/report.md`, `complexes/report.md` with their scripts.
- For the deficiency-one theorem: the statement in Joshi–Shiu, *A survey of methods for deciding whether a reaction
  network is multistationary* (arXiv:1412.5257v3, Theorem 3.2). For the notion of stability: Jiao–Tang–Zeng
  (arXiv:2406.11586).

I checked every proof of the dimension-two section by hand before running anything. My own scripts are in
`verification/`.

## Verdict: **CONFIRMED WITH CORRECTIONS**

The main result is correct. Theorem 15 (dimension two) holds: a two-dimensional at-most-bimolecular mass-action
system with ACR has at most one exponentially stable positive steady state per compatibility class. Its proof is
complete. I found no missing case, no wrong sign and no hidden assumption. The proof in fact shows more: two
nondegenerate positive steady states in one class always have Jacobian determinants (on S) of opposite signs.

The consequences are also correct, and the cited KNS results are used within their hypotheses:
- dimension ≥ 3;
- at least 4 reactant complexes;
- at least 6 complexes, by the kinetic-subspace lemma and the deficiency-one theorem.

The three-species example N1 is correct. I re-derived it from the reaction list with my own exact code (Gröbner basis,
exact root isolation, Routh–Hurwitz, arb eigenvalue enclosures). So the conjecture is settled as claimed:
- true for dimension, reactant complexes and complexes;
- false for the number of species.

Numerically, I attacked the dimension-two proof in two ways, and neither found a counterexample:
- a linear-programming search over the relaxed class of vector fields that the proof actually covers: 1 700
  configurations with two exponentially stable zeros on a class line, every one violating global ACR;
- a scan of 47 847 dimension-two networks on three species: no ACR system with two stable states, and opposite
  determinant signs in every class with two nondegenerate states.

The corrections are needed but do not touch the main results. The mathematics of the main theorem is unchanged.
One side statement in §6 (item 3) was false as written. It claimed that no three-species network with at most 5
reactions has three positive steady states, and I found an exact two-dimensional counterexample with 4 reactions.
I restricted it to full-dimensional networks, which is what its enumeration proves and all the paper needs. The rest
is minor or editorial. Details are in §6 and §7 below.

## 1. Source and statement

- **Quote.** The §6 sentence quoted in `README.md` and `paper.tex` matches arXiv v1 word for word, apart from the
  typography of the inner quotation marks.
- **KNS bounds.** Theorem 1.2 there gives, for at-most-bimolecular networks with ACR and nondegenerate
  multistationarity: n ≥ 3, at least 3 reactant complexes, at least 5 complexes, dimension ≥ 2. All four are
  attained by Example 4.10, {X1+X2 → 2X3, X3 → X1, 2X3 → 2X2}. The paper's table and text state this correctly.
- **Definitions.**
  - ACR follows KNS Definition 2.6: a positive steady state exists, and the same x_i at all positive steady states.
  - Multistability means two exponentially stable positive steady states in one class, for one κ. "Exponentially
    stable" means all eigenvalues of Df|_S have negative real part.
  - This agrees with Jiao–Tang–Zeng ("nondegenerate, and all non-zero eigenvalues of the Jacobian have negative real
    parts"). The two notions are equivalent because Df maps R^n into S.
- **Logic.** Multistability implies nondegenerate multistationarity, so every KNS lower bound applies. The question
  is existential in κ, as in KNS Theorem 1.2.

## 2. The dimension-two theorem (Theorem 15), step by step

Notation as in the paper:
- ACR is in X1 with value c, and H = {x1 = c};
- f = ψα + (x1 − c)Pℓ in the normal form.

| Step | Check | Result |
|---|---|---|
| Lemma 5 (coefficient cones) | (y′ − y)_i = y′_i ≥ 0 off supp(y); Σ(y′ − y) = \|y′\| − 2 ≤ 0 if \|y\| = 2; both pass to nonnegative combinations. | correct |
| Lemma 6 (e1 ∉ S^⊥) | F(x) = (f(x), Wx) has invertible derivative at a nondegenerate x*, so the inverse function theorem varies the conserved quantity x1 if e1 ∈ S^⊥. Then S ∩ e1^⊥ is a line Rℓ, and each class meets H in a segment parallel to ℓ. | correct |
| Lemma 7 (f = hα on H) | f^{-1}(0) is an (n−2)-manifold near p and lies in H by ACR. Irreducible factors vanishing on open pieces are found by a finite-union (Baire) argument. Coprime real polynomials stay coprime over C, so their common zero set has dimension ≤ n − 3 and so do its real points. The degree-1 factor case uses nondegeneracy of all steady states and the isolation of p. The case n = 2 (H a line, "dimension −1" = empty) also works. | correct; one sentence added (α ≠ 0) |
| Lemma 8 (α1 = 0) | ℓ is an eigenvector with eigenvalue γh_L′(t), and h_L′(0), h_L′(τ) have opposite signs. | correct |
| Normal form | {α, ℓ} is a basis of S; φ vanishes on H and has degree ≤ 2, so φ = (x1 − c)P with P affine. | correct |
| Lemma 10 (det formula) | On H, Df = α∇ψᵀ + Pℓe1ᵀ. The matrix in the basis (α, ℓ) is [[∇ψ·α, ∇ψ·ℓ], [P, 0]]. | correct (also checked numerically: `scripts/check_dimension_two_example.py`, and my `global_step_example.py`: det 12 = −P ∇ψ·ℓ) |
| Corollary 11 | Stable in dimension 2 ⟹ det > 0; ∇ψ·ℓ = h_L′ has opposite signs at p and q. | correct |
| Lemma 12 (V⁺ connected, ∇ψ_H ≠ 0) | See below. | correct; typo fixed (j = k) |
| Lemma 13 (global ACR) | Π^± are convex and nonempty, accumulating at {c}×Z along e1 − (π1/πk)e_k. ψ has no zero on them, by ACR. If ε₊ ≠ ε₋, then P_H divides ψ_H and h_L ≡ 0. | correct |
| Lemma 14 (tangency) | z minimizes εψ_H on Z, so ∇ψ_H(z) ⊥ T; εQ ≥ 0 on T; the identity s(μ + 2B(v,ν) + sQ(ν)) = −Q(v) is exact for a quadratic ψ_H, and the sign conclusion follows. | correct |
| Theorem 15 | θ is locally constant: μ(x′) = ∇ψ_H·∇P_H/\|∇P_H\|² is continuous on V⁺ and equals the μ of Lemma 14 on V⁺ ∩ Z. V⁺ is connected, hence θ is constant. | correct |

**Lemma 12 in detail** (this is where the products are bimolecular).

The quadratic part Q of ψ_H has coefficients a_jk with v_{X_j+X_k} = a_jk α. If a_jk > 0, then:
- α_i ≥ 0 off {1, j, k};
- Σ_{i∈{j,k}} α_i ≤ −1.

So the set N = {i ≥ 2 : α_i < 0} satisfies N ⊆ {j,k} and N ∩ {j,k} ≠ ∅. Since Q(ℓ′) ≠ 0, we get 1 ≤ |N| ≤ 2. For
i ∈ N, every monomial of ψ_H without x_i has a coefficient ≤ 0. I re-derived the case analysis:
- **N = {a, b}:** ψ_H = βx_ax_b − (affine with coefficients ≥ 0). V⁺ is a graph over a convex set, or a slice of a
  hyperplane.
- **N = {a}:** ψ_H = β_a x_a² + x_a(K + R) − G, with G affine with coefficients ≥ 0. V⁺ is a graph over a convex
  set, or convex. If G > 0 and β_a > 0, the root is the unique positive one and is simple.

In every case ∇ψ_H ≠ 0 on V⁺, and no sub-case is missing: G > 0 or G ≡ 0 on the open orthant, because G is affine
with nonnegative coefficients. The paper wrote "α_j + α_k ≤ −1", which double-counts when j = k (monomial x_j²). The
conclusion is unaffected; corrected to a sum over the set {j, k}.

**Adversarial points examined.**
- **Which species has ACR.** This is only a relabelling.
- **α1 = 0.** Lemma 8.
- **e1 ∈ S^⊥.** Lemma 6.
- **f|_H ≡ 0, or α = 0.** p would not be isolated in its class; I added one sentence making this explicit.
- **h_L of degree ≤ 1.** Excluded in Lemma 7.
- **P not depending on x2, …, xn.** Then P_H is constant and cannot change sign; Lemma 13 notes this.
- **Π^± empty.** Excluded by the explicit direction e1 − (π1/πk)e_k.
- **n = 2.** All steps degenerate correctly.
- **A boundary root.** All sets used are open subsets of the positive orthant.
- **The kinetic subspace is smaller than S.** Then every steady state is degenerate, so it is irrelevant.

**Hypotheses actually used.** The proof uses deg f ≤ 2, f ∈ S with dim S = 2, global ACR, and the two cone
properties (C1) and (C2) of Lemma 5. It uses nothing else about the network. This is what makes the relaxed numerical
test below a test of the proof itself.

**Numerical tests** (floating point, searches for counterexamples, not proofs).

1. `verification/dim2_relaxed_lp.py` (run with 1500 samples per setting, seed 2026; log `dim2_relaxed_lp.log`).
   - **Setup.** The search covers vector fields f = ψα + (x1 − c)Pℓ (n = 3 and 4) satisfying only (C1)–(C2).
   - **Sampling.** α and ℓ are drawn in stoichiometric subspaces spanned by bimolecular reaction vectors, with
     special and random values of α along ℓ, and also generically. p and q are drawn on a common class line in H.
   - **Linear programs.** For fixed (α, ℓ, c, p, q), the local conditions are linear in the coefficients of
     (ψ, P): zeros at p and q, P(p)P(q) < 0, det > 0 and trace < 0 at both points, (C1) and (C2). They are solved
     by HiGHS with random objectives.
   - **Global check.** For each feasible (ψ, P), a search for zeros of ψ on Π^± follows. Such a zero means a
     positive steady state off H, so no ACR.
   - **Theorem mode** (all local conditions). 1 524 locally feasible configurations for n = 3 and 176 for n = 4.
     Every one violates global ACR.
   - **Global-lemma mode** (no stability conditions, only P(p)P(q) < 0 and deg h_L = 2). This tests Lemmas 12–14
     and the final step on their own. There are 708 + 72 feasible configurations, and every one violates global
     ACR.
   - **Control without (C2).** The same holds there: all 2 998 + 530 sampled configurations violate global ACR.
     This agrees with Remark 16 of the paper; I found no counterexample without (C2) either.
   - **Runtime.** 823 s.
2. `verification/global_step_example.py` (exact). A genuine bimolecular network of dimension 2 shows that the local
   configuration excluded by Lemmas 13–14 really occurs when ACR fails:
   - the network is {X2+X3 → 2X1, X1+X2 → 2X2, X1+X3 → 2X3, X1 → X2, X1 → X3, X2 ⇄ X3};
   - its normal form has c = 1 and P = x2 − x3;
   - in the class x1 + x2 + x3 = 5, the points (1, 1, 3) and (1, 3, 1) are both exponentially stable (trace −13,
     det 12) and P changes sign between them;
   - a third steady state, (37/11 − 4√34/11, 9/11 + 2√34/11, 9/11 + 2√34/11), has x1 ≠ 1.

   So the global ACR argument is genuinely needed, and the paper uses it correctly.
3. `verification/dim2_network_scan.py` (log `dim2_network_scan.log`).
   - **Networks.** All at-most-bimolecular networks on 3 species with dim S = 2 and 3 or 4 reactions (up to species
     permutation), plus a sample of 60 000 networks with 5 reactions, after the filters "positive flux" and "≥ 3
     reactant complexes": 47 847 networks.
   - **Method.** 3 random rate vectors each. Steady states are computed in random classes by a resultant and
     Newton polishing. ACR is detected numerically and re-checked on 60 more classes.
   - **Result** (1 136 s on 12 cores). There are 48 907 (network, rate vector) pairs with ACR. Of these, 107 pairs
     (59 networks) have a class with at least two nondegenerate positive steady states; KNS Example 4.10 is among
     them, which checks that the search works. None has a class with two exponentially stable positive steady
     states. In all 4 875 classes with exactly two nondegenerate positive steady states, the determinants of Df|_S
     have opposite signs.
   - **Limitation.** Systems that have ACR only on a measure-zero set of rate constants cannot be found this way.
4. The solver's `scripts/check_dimension_two_example.py` passes: 20 exact instances of KNS Example 4.10, opposite
   determinant signs, Lemma 10 formula.

## 3. Consequences (Corollaries 17, 19, Lemma 18, Remark 20)

**Corollary 17 (dimension ≥ 3, at least 4 reactant complexes).**
- KNS Proposition 3.10 covers any one-dimensional network with ACR: no nondegenerate multistationarity.
- Theorem 15 excludes dimension 2.
- KNS Theorem 4.1(1) gives n ≥ 3, which is needed for Theorem 5.5.
- Theorem 5.5 (n ≥ 3, k ≥ 1 conservation laws, nondegenerately multistationary with ACR) gives at least
  n − k + 1 = dim S + 1 ≥ 4 reactant complexes.
- Theorem 5.1 (full-dimensional) gives at least n + 2 ≥ 5.

All hypotheses are met. Correct.

**Lemma 18 (δ ≥ t − ℓ + 1).**
- K_κ = Y(Im A_κ) contains the image of Df, so nondegeneracy gives K_κ = S.
- dim ker A_κ = t, with kernel vectors supported on terminal complexes (Feinberg–Horn).
- Ψ(x*) > 0 on a nonterminal complex, so A_κΨ(x*) ≠ 0 lies in ker Y ∩ Im A_κ.

Hence s ≤ m − t − 1. Correct.

**Corollary 19 (at least 6 complexes).**
- δ ≥ 1, by the deficiency-zero theorem, valid for all networks (KNS Lemma 2.4, first sentence).
- So m = ℓ + s + δ ≥ 5, with equality only for ℓ = 1, s = 3, δ = 1.
- Then t = 1: by Lemma 18 if a nonterminal complex exists, and by weak reversibility otherwise.
- The deficiency-one theorem then allows at most one positive steady state per class.

Correct.

**Deficiency-one theorem: hypotheses and reference.** The uniqueness statement needs:
1. each linkage class has exactly one terminal strong linkage class;
2. δ_θ ≤ 1;
3. Σδ_θ = δ.

It does not need weak reversibility. This is how Feinberg's theorem is stated in the Joshi–Shiu survey
(arXiv:1412.5257, Theorem 3.2, attributed to Feinberg). The paper cites Feinberg, ARMA 132 (1995), which is the
standard source.

I could not open the 1995 paper to give a theorem number. The citation is acceptable but has no precise
location (minor, not changed).

**Remark 20 (the gap in KNS Theorem 4.9(2)).** KNS Lemma 2.5 is stated for weakly reversible networks. Their proof of
Theorem 4.9(2) applies it to an arbitrary network with at most 4 complexes, ℓ = 1, δ = 1, where t = 2 would not be
excluded. Lemma 18 gives t = 1 there, so the KNS bound stands. The general bound m ≥ s + 3 for two nondegenerate
states in a class is also correct.

## 4. Examples and certificates

**Solver's scripts, re-run under WSL (Fedora 44, Python 3.12, sympy 1.14, python-flint 0.9).**

| Script | Result |
|---|---|
| `scripts/check_dimension_two_example.py` | PASS |
| `complexes/scripts/certify_examples.py` | ALL PASS |
| `complexes/scripts/eigen_enclosures.py` | ALL PASS |
| `examples/scripts/certify_example.py` | all PASS |
| `examples/scripts/check_gale_formula.py` | PASS on 40 instances |
| `examples/scripts/e2_reactant_sets.py 8` | 5 PASS, 5 s; 26 sets for n = 4, 8 structures for n = 5, 1 for n = 6, none for n = 7, 8 |
| `examples/scripts/enum5_reactions.py` (10 processes, 114 s) | PASS |

The counts printed by `enum5_reactions.py` agree with `examples/report.md`:
- 3 188 646 networks in total;
- 1 251 698 of rank 3;
- 628 446 consistent;
- pole counts 36 303 / 265 885 / 253 227 / 73 031;
- sign changes of g′: 298 879 with none, 27 379 with one, none with two or more.

**My own certificate:** `verification/certify_examples_own.py` (0.4 s, ALL PASS). It is written separately from the
solver's scripts and builds f from the reaction lists. For N1 (the main example), N2 and G7 it checks:
- **counts**, computed from the reaction graph: complexes, reactant complexes, dimension, linkage classes, terminal
  strong linkage classes, deficiency. N1: 6, 6, 3, ℓ = 1, t = 1, δ = 2. N2: 6, 6, 3, ℓ = 2, t = 2, δ = 1. G7: 7, 6, 3;
- **the ACR identity**, symbolically for all rate constants: N1 gives z(2k1 − k5x), N2 gives x(k5 − k7x), G7 gives
  y(k1x − k2);
- **existence of a positive steady state for every κ** (N1: a cubic with constant term k3k6c² > 0 and leading
  coefficient −k2k4);
- **the positive steady states**, computed with a lex Gröbner basis of ⟨f, 1 − txyz⟩, exact real-root isolation and
  back substitution. There are exactly the three claimed points for each network;
- **stability**, from exact characteristic polynomials and Routh–Hurwitz, with arb enclosures (200 bits) of all
  eigenvalues certifying the sign of every real part.

For N1 the characteristic polynomials are:
- λ³ + 45λ² + 206λ + 30 at (1, 19/4, 1/2), stable;
- λ³ + 46λ² + 47λ − 50 at (1, 3, 1), with eigenvalue +0.6476…;
- λ³ + 81λ² + 632λ + 1050 at (1, 7/4, 7/2), stable.

I also checked the three Jacobian matrices printed in Theorem 21 entry by entry by hand.

**Sharpness claims.** The values in the paper's table are confirmed:
- dimension 3, species 3 and complexes 6 are attained by N1;
- the reactant-complex minimum lies between 4 (proved) and 6 (N1, G7), and the paper presents it as open.

**§6 items 3–4.** Item 3 (at least 6 reactions on three species) is correct after the scope correction below. The
enumeration covers full-dimensional networks with 5 distinct reactant complexes. Fewer reactant complexes are covered
by KNS Propositions 3.1(3) and 3.5. Full-dimensionality follows from Corollary 17. Two stable states force a third
steady state through the Gale-dual function: in dimension 3, stable means det J < 0.

Item 4 (four reactant complexes) is presented as partial and partly numerical. I re-ran its script but did not
re-prove Propositions 4.5–4.6 of `examples/report.md`. They are not needed for the conjecture.

## 5. Novelty spot-check (27 October 2026)

- **OpenAlex.** OpenAlex lists 4 works citing KNS (W4390503227), the same four as in `references.md`. None treats
  ACR with multistability.
- **arXiv, ACR with bistability or multistability.** An abstract search for "concentration robustness" together with
  bistab* or multistab* returns only Joshi–Nguyen–Johnston, arXiv:2605.11389 (12 May 2026). That paper has six
  species and dimension four, and does not discuss minimal sizes.
- **arXiv, newest ACR papers.** The newest ACR abstract is Chen, arXiv:2608.18509 (19 Aug 2026). It classifies
  degenerate two-dimensional zero-one networks, with no multistability. It is not listed in `references.md`; this is
  harmless.
- **GitHub.** `gh search repos "concentration robustness"` shows only ACR tools and `mmislan/palomar-p082-acr-zero-divisors`
  (zero-divisor criteria, no multistability). The repositories of `ipitchford` have no ACR or multistability title.
- **Zenodo.** The search for "absolute concentration robustness" gives 3 records, none on multistability.

The conclusion of `references.md` (open as of 27 October 2026) is supported.

## 6. Issues found

| # | Severity | Location | Issue |
|---|---|---|---|
| 1 | moderate (false as stated, not used) | `paper.tex` §6 item 3; `README.md`; `examples/report.md` Theorem 2.5; `examples/scripts/enum5_reactions.py` messages | The text said that no bimolecular three-species network with ≤ 5 reactions has three positive steady states in a class. This is **false**. The two-dimensional 4-reaction network {0 → Z, Y+Z → X, X+Z → Y+Z, X+Y → 2Y} with κ = (10, 10, 1, 1) has three nondegenerate positive steady states (5,1,1), (4,2,1/2), (3,3,1/3) in the class x + y = 6, two of them exponentially stable (`verification/item3_scope_example.py`, exact). I found it with an exploratory scan of dimension-two networks. The enumeration only covers full-dimensional networks with 5 distinct reactant complexes, and that is all the consequence needs. The consequence actually used, at least 6 reactions for ACR with multistability on three species, is correct, because such a network is full-dimensional by Corollary 17. |
| 2 | minor | `paper.tex` §6 item 4 | "For 4 species no class contains three positive steady states" is proved only when one steady state is nondegenerate (the coset structure needs rank N = 3). The statement needed, and proved, is that no class contains two exponentially stable ones. |
| 3 | typo | Lemma 12 | "α_j + α_k ≤ −1" double-counts when j = k; the correct inequality is Σ_{i∈{j,k}} α_i ≤ −1. The conclusion is unaffected. |
| 4 | gap of one sentence | Lemma 7 | α ≠ 0 was asserted but not argued. If f\|_H ≡ 0, p is not isolated in its class. |
| 5 | editorial | `scripts/check_dimension_two_example.py` | Its docstring referred to "Theorem 3.10" and "Lemma 3.5", which do not exist in the compiled numbering (Theorem 15, Lemma 10). |
| 6 | editorial | `paper.tex` header; `README.md` | The status read "preliminary / verification pending". |
| 7 | remark, not changed | Remark 20; references | The deficiency-one theorem is cited without a theorem number. The hypotheses used are the right ones (checked against the Joshi–Shiu survey, Theorem 3.2). |
| 8 | cosmetic, not changed | `paper.pdf` | Six overfull lines (the summary table, two displays, long file names in §5–6) were already there before my edits; I left them. |

No issue affects the dimension-two theorem, the bounds 3 / 4 / 6, or the example N1.

## 7. Corrections applied

1. **`paper.tex` §6, item 3.** Restated as "no full-dimensional bimolecular network with 3 species and at most 5
   reactions has three nondegenerate positive steady states, or two exponentially stable ones". The text now says
   how fewer reactant complexes are covered (KNS Propositions 3.1(3), 3.5) and why two stable states force a third
   one. Full-dimensionality is now derived from Corollary 17. A sentence gives the two-dimensional counterexample to
   the unrestricted statement. The same qualifier ("full-dimensional") is added to Theorem 2.5 and the summary table
   of `examples/report.md`, and to the docstring and PASS messages of `examples/scripts/enum5_reactions.py` (text
   only; the code is unchanged).
2. **`paper.tex` §6, item 4.** "no class contains three positive steady states" became "no class contains two
   exponentially stable positive steady states", with the reason.
3. **`paper.tex` Lemma 12.** Σ_{i∈{j,k}} α_i in place of α_j + α_k.
4. **`paper.tex` Lemma 7.** Added: "Likewise f|_H ≢ 0 (otherwise p would not be isolated), so the vector α obtained
   above is nonzero."
5. **`paper.tex` header.** Subtitle "verified"; credits now name the verification.
6. **`README.md`.**
   - Status line: verified, verdict, link.
   - Attribution: verification.
   - "Further results": the scope of the 5-reaction enumeration.
   - Reproduce list: the verifier's scripts.
   - Files table: `verification/`.
7. **`scripts/check_dimension_two_example.py`.** Theorem numbers in the docstring and output line now match the
   compiled paper; the code is unchanged and the script still passes.
8. **`paper.pdf`.** Recompiled with `py -3 templates/build_tex.py problems/kaihnsa-acr-multistability/paper.tex`.

## 8. Verifier's files (`verification/`)

| File | Content | Runtime |
|---|---|---|
| `certify_examples_own.py` | Exact certificate of N1, N2, G7 from the reaction lists (Gröbner basis, real roots, Routh–Hurwitz, arb) | < 1 s |
| `global_step_example.py` | Exact two-stable configuration on H without ACR, showing that the global step of the proof is needed | < 2 s |
| `item3_scope_example.py` | Exact two-dimensional 4-reaction network with three nondegenerate positive steady states in a class (scope of §6 item 3) | < 1 s |
| `dim2_relaxed_lp.py`, `dim2_relaxed_lp.log` | LP search against Lemmas 12–14 and Theorem 15 on the relaxed class (C1)–(C2), with a control without (C2) | about 14 min (one core) |
| `dim2_network_scan.py`, `dim2_network_scan.log` | Scan of dimension-two networks on 3 species: ACR detection, stability counts, determinant signs | about 19 min on 12 cores |
