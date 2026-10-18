Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Verification: Banaji–Feliu Remark 5.8 (quadratic (n, n+2, n) networks with n+1 equilibria)

## Verdict

**CONFIRMED WITH CORRECTIONS.**

The mathematics is correct. For every n ≥ 1 the network 𝒩ₙ(α, β) is a quadratic (n, n+2, n) network that satisfies
all hypotheses of Banaji–Feliu Theorem 5.5. With αⱼ = 64^(j−1), βⱼ = 1 and the stated rate constants it has exactly
n+1 positive nondegenerate equilibria, so the bound is attained. The side remark on Remark 5.6 is also correct.

**The result is not new.**
- An unrefereed manuscript posted on Zenodo and GitHub on 4 October 2026, two weeks before this work, answers Remark 5.8
  affirmatively for every n ≥ 2 (Anonymous, maintained by I. Pitchford, doi:10.5281/zenodo.23146115).
- The equilibrium equations of 𝒩ₙ are exactly the extremal fewnomial system G_ε of Phillipson–Rojas (Contemp. Math. 605,
  2013, Theorem 1.6). The R = 64 networks with the rate constants of the proof have the same Gale data as G_{1/8}.

Following PROTOCOL §1.3, the status is changed from "solved (preliminary)" to **"already solved"**. What remains of
the contribution here is an apparently independent rediscovery with a short, self-contained, verified proof.

## 1. Source and statement (read in full)

- Banaji–Feliu, arXiv:2409.06877v5, downloaded from arXiv by me (identical in size to the local copy). The published
  version (PMC13481394) was also checked: same numbering and wording for Theorem 5.5 and Remarks 5.6 and 5.8.
- The quote of Remark 5.8 in `paper.tex`/`README.md` is verbatim (p. 41).
- **Meaning of "quadratic"** (BF, Sect. 2, p. 5), verbatim: "A CRN is quadratic if all its sources have molecularity two
  or less; equivalently, each polynomial on the right hand side of (2.1) has degree two or less. It is bimolecular if
  all sources and products have molecularity two or less."
  Products are therefore unrestricted. Remark 5.8 itself treats product molecularity as a secondary measure ("no
  quadratic (2, 4, 2) network with product molecularity of less than 6 …"), so large products (64^(n−1)) are admissible.
  The solver's (3,5,3) example has products of molecularity ≤ 20.
- **Hypotheses of Theorem 5.5**, checked by hand and in code from the reaction lists:
  - k = 2, so m = n+2 reactions;
  - n+2 distinct sources of molecularity ≤ 2;
  - rank Γ = n, so the only stoichiometric class is ℝⁿ_{>0};
  - rank[A|1] = n+1, so the source rank equals n = r and the network is not source deficient (BF §4.1);
  - dynamically nontrivial (ker Γ = span(α, β) contains a positive vector);
  - each product differs from its source.
- **Nondegeneracy** in BF's sense means that J_g(x) is nonsingular on im Γ = ℝⁿ, i.e. det J ≠ 0. This is what the
  manuscript proves (Lemma 5) and what I checked.
- **Compatibility checks.** BF Remark 5.4(1) and 5.7 say the bound can only be attained if A has one or two rows
  containing a 2 and w has n entries ±2 and two entries ±1. The family meets this: one row 2X₁, and
  w = (∓1, 1, −2, 2, …). Within the family, n = 2 with three equilibria needs α₃β₂ > 9α₂β₃, so the product molecularity
  is ≥ 6. This agrees with BF's Mathematica claim for (2,4,2).

## 2. Recomputed independently (own code, different method)

All scripts are in `verification/` and were run in WSL, single-threaded, with ≤ 3 GB. They never use the Gale reduction.
Logs are in `verification/*.log`.

### `direct_check.py`

Each network is written as a list of reactions with exact rational rate constants. The script builds the mass-action
vector field g and then runs three independent checks.

- **msolve** (exact Gröbner basis plus certified real root isolation) on g = 0, t·x₁⋯xₙ = 1. It counts the solution
  boxes that are strictly positive, and checks that every other box has a negative coordinate.
- **Singular** checks that ⟨g, det Dg, t·x₁⋯xₙ − 1⟩ = ⟨1⟩, so no equilibrium in (ℂ*)ⁿ is degenerate. It also
  computes the number of complex solutions.
- **PARI/GP** computes a Sturm count of the positive roots of the eliminant in x₁. The eliminant comes from Singular
  and is squarefree and nonzero at 0.

Results:

| Case | Positive nondegenerate equilibria |
|---|---|
| n = 1: 0→X (1), 2X→3X (1), X→0 (3) | 2 |
| Small n = 3 example, κ₁ = **1/1600** (my choice; others 1) | 4 |
| Negative control, same network with κ₁ = 1/100 | 0 (the test discriminates) |
| **R = 64 family, n = 1…10**, with exactly the proof's constants κ₁ = e^c = 2^(−(6n+3)) (n odd), 2^(6n+3) (n even), others 1 | n+1 |
| Computer-found (3,5,3) example, κ(2Z→Z) = 9971888893/778 | 4 |

- **Counts are exact.** In every case all complex solutions in (ℂ*)ⁿ are real and positive, and their number equals
  the count in the table.
- **Rate constants of the proof.** The proof's constants are rational, which I derived myself: M₁ = −nL (n odd) or
  (n+1)L (n even). The solver's script used a nearby rational C instead, so the exact constants of the proof are checked
  here for the first time.

### `uniform_signs_arb.py`

Ball arithmetic (Arb, 200 bits), R = 64, n = 1…300. It checks that (−1)^(i+1)(f(Rⁱ) − c) > log 2 − 4/63 for every
i ≤ n, and the closed form of e^c.

### `prior_art_check.py`

The same pipeline is applied to three groups of networks, all of which pass:

- the earlier claim's (3,5,3) example, with product molecularity 10 and rate 2/3137: 4 positive nondegenerate
  equilibria;
- the earlier claim's all-n family for n = 2…6: n+1 equilibria;
- 𝒩ₙ realizing Phillipson–Rojas G_{1/4} (R = 16) and G_{1/8} (R = 64) for n = 2…6: n+1 equilibria. The script also
  checks, for n ≤ 100, that the Gale constant c of the G_{1/8} realization equals the constant c of the proof.

## 3. Proof re-read line by line (`paper.tex`)

Every step of Lemmas 3–6, Proposition 7 and Theorem 1 was rederived by hand.

- **Lemma 3 (structure).** Γ = [a | b | −Iₙ], so ker Γ = {(s, t, sa + tb)}. The affine dependency w was checked
  coefficient by coefficient, including n = 1.
- **Lemma 4 (bijection).** I rederived it without Gale language from the triangular structure:
  xₙ = (αₙ₊₁κ₀ + βₙ₊₁κ₁x₁²)/κₙ₊₁, then x_{n−1}, …, x₁ by back substitution. This gives the same scalar equation f(u) = c.
- **Lemma 5 (nondegeneracy).** The argument with w′ (w′·1 = 1, im S = {w·q = w′·q = 0}) is correct.
- **Lemma 6 (degree).** The leading coefficient is ∏βⱼ(1 + Σwⱼ) = ±∏βⱼ. Rolle's theorem gives ≤ n+1 zeros, and
  simplicity when there are exactly n+1. The original text left the "≤ n+1" step implicit; it is now stated.
- **Proposition 7 (uniform estimate).**
  - log(R^m + Rⁱ) = max(m, i)L + log(1 + R^(−|m−i|)).
  - |εᵢ| ≤ Σ_d 2·2·R^(−d) = 4/(R−1).
  - M_{i+1} − M_i = L(1 + 2Σ_{m≤i}(−1)^m) = ∓L.
  - The margin is L/2 − 2 log 2 = log 2 for R = 64, which exceeds 4/63.
  - The limits at 0⁺ and ∞ are correct, and the n+2 alternating signs give n+1 zeros.
  - The estimate is correct and uniform in n. It would even work for R = 32 (margin ½ log 2 > 4/31).
- **Gap check.** No gap found. The proof is complete and elementary.

## 4. Solver's scripts (re-run)

- `certify_family.py`: exit 0, **100** PASS lines, identical to `certificates/certify_family.log`. It checks the
  Gale-dual polynomial P_C for a rational C near the midpoint of the admissible window, using a Sturm count and
  gcd(P, P′) = 1. The checks behave as follows:
  - The "every local minimum below every local maximum" check uses 60-digit floats. It only serves to choose C and is
    not a certificate.
  - The "det J ≠ 0" check is a weak float test (|d| > 0). Nondegeneracy actually follows from squarefreeness via
    Lemma 5, so this weakness does not matter.
  - The PASS lines do test what they announce.
  - The script does not test the exact rate constants of the proof; `direct_check.py` now does.
- `other_example_n3.py`: exit 0. It certifies 4 simple positive roots of P_C for C = 778/9971888893. It does not test
  nondegeneracy directly (it relies on Lemma 5, whose argument applies to any full-rank (n, n+2, n) network with
  rank[A|1] = n+1). `direct_check.py` confirms nondegeneracy directly.
- The manuscript said "101 checks". The correct number is 100, and it has been corrected.

## 5. Side remark (Remark 5.6 vs BBH Lemma 20)

- **Both texts read.** BF Remark 5.6, verbatim: "by [5, Lemma 20], quadratic rank-one networks do not admit multiple
  positive nondegenerate equilibria, regardless of the number of species". [5] is Banaji–Boros–Hofbauer, Nonlinear
  Dyn. 2024; I read arXiv:2406.13451v1, whose numbering matches BF's citations, e.g. "[5, Remark 38]".
- **BBH Lemma 20 says the opposite.** A quadratic rank-one network admits multiple positive nondegenerate equilibria if
  and only if it has at most two nontrivial species and, after removal of trivial species, contains three reactions
  simply equivalent to (1) 0→X, X→0, 2X→3X, or (2a), (2b), (2c). BBH Theorem 18 and Corollary 19 ("A *bimolecular*,
  rank-one mass-action network admits at most one positive nondegenerate equilibrium") confirm this reading.
- **Counterexample.** ẋ = 1 − 3x + x² has the roots (3 ± √5)/2, so the claim of Remark 5.6 is false as stated, and the
  bound n+1 = 2 is attained for n = 1.
- **Not new either.** The same discrepancy is noted in the earlier manuscript's prior-art notes (PRIOR_ART.md,
  "Source caveat").

## 6. Novelty (the main problem found)

The solver's check covered Semantic Scholar, OpenAlex and mmislan. I confirmed it: there are two citing works (Feliu
et al., 2024), both irrelevant, and the 23 mmislan repositories contain nothing on this question. The solver missed
the following:

- **Earlier answer.** Anonymous, *Sharp positive-equilibrium counts in quadratic mass-action networks with two more
  reactions than species*, v1.0.0-candidate, unrefereed.
  - Zenodo doi:10.5281/zenodo.23146115, status "published" on 2026-10-04 according to the Zenodo API; GitHub
    ipitchford/quadratic-mass-action-sharpness, created 2026-10-04.
  - It explicitly answers Remark 5.8 for every n ≥ 2.
  - Its networks have sources 0, X₁, Xᵢ+Xᵢ₊₁, 2Xₙ, which is the same support as 𝒩ₙ after reversing the order of the
    species, but with different reaction vectors.
  - It gives a (3,5,3) example with products of molecularity ≤ 10.
  - I verified its explicit networks exactly (n = 2…6, plus the (3,5,3) example). I did not review its general proof.
- **Same polynomial system.** Phillipson–Rojas 2013, Theorem 1.6 (arXiv:1011.4128v5, p. 3, read from the rendered
  page):
  - G_ε: x₁x₂ = ε + x₁², xᵢxᵢ₊₁ = 1 + ε^(2i−3)x₁², xₙ = 1 + ε^(2n−3)x₁²;
  - it has exactly n+1 nondegenerate positive roots for ε = 1/4;
  - this is the equilibrium system of 𝒩ₙ with αⱼ = 16^(j−2), βⱼ = 1, κ = (1/4, 1, 1, 4, 4³, …), i.e. R = 16.
- **Sharp circuit bound.** Bihan 2007 (J. London Math. Soc. 75) shows that n+1 is the sharp maximum number of positive
  solutions for circuit-supported systems.
- **Banaji's later work.** Banaji, arXiv:2609.07265 (Sep 2026), still quotes only the three-equilibrium (2,4,2) example.

## 7. Corrections applied (all inside this folder)

- **`paper.tex`** (recompiled to `paper.pdf`, 6 pages):
  - subtitle and credits updated (verification done; rediscovery);
  - abstract: "the answer is not new" paragraph added;
  - new Section 5 "Prior work" (Bihan, Phillipson–Rojas with the exact dictionary to 𝒩ₙ, the 4 Oct 2026
    manuscript);
  - Remark 2: precise paraphrase of BBH Lemma 20, and the ambiguous "κ = (1,1,3)" replaced by explicit rate constants
    per reaction;
  - Lemma 6: the "at most n+1 zeros" step made explicit;
  - Theorem 1 proof: the closed form κ₁ = 2^∓(6n+3) (rational) added;
  - small example: explicit κ₁ = 1/1600 instead of "a suitable rational κ₁";
  - second example: rate constant stated;
  - mention of the earlier (3,5,3) example with product molecularity 10;
  - "101 checks" corrected to "100 checks";
  - paragraph on the independent verification added;
  - section labels added; one over-wide display tightened.
- **`references.bib`**: added PhillipsonRojas2013, Bihan2007 and Anon2026.
- **`README.md`**:
  - status changed to "already solved" (independent, verified proof);
  - verifier attribution added;
  - "Prior work" section added;
  - file table and reproduction commands extended.
- **`references.md`**: new sources added and the novelty check updated.
- **`verification/`**: new scripts and logs. The solver's `scripts/` and `certificates/` are unchanged.

## 8. Reservations

- **Earlier manuscript only partly checked.** I verified its explicit examples, not its general proof. It is
  unrefereed and anonymous, and its priority is a matter of public record (Zenodo timestamp). The algebraic core in any
  case goes back to Phillipson–Rojas (2013), whose real case I did not re-prove. I only checked the dictionary between
  G_ε and 𝒩ₙ, and the cases n = 2…6.
- **Not verified for independence.** Whether the solver had seen these works is outside my knowledge. The text says
  "apparently independent".
- **Not edited (outside this folder).** The root `README.md` results table should be updated to "already solved"
  (independent verified proof), with the verifier credit.
- **Cosmetic.** Two displays remain overfull by 6–8 pt.
