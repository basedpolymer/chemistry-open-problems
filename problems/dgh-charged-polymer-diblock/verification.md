# Independent verification: DGH charged di-block, critical temperature

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

Date: 18 October 2026. Files not opened, as instructed: `PROGRESS.md`, `BRIEF-codex.md`, `verification-request.md`, `verification-codex-internal.md`, `travail/`, `rapports/` (see the disclosure in §8).

## 1. Verdict

**CONFIRMED WITH CORRECTIONS.** The corrections are minor and editorial; none affects the mathematics.

The theorem in `paper.tex` is correct, and it answers Poisat's question (3) in the sense the source intends. In the original convention, the critical inverse temperature of the balanced di-block is the unique root of
\(\sum_{\ell\ge1}e^{-\beta\ell^2/2}=1\), namely \(\beta_*=0.6981317007991994390913413701801127970\ldots\). This is \(2\beta_0\) in Poisat's own notation (3.23). In temperature units, \(T_c=1/\beta_*=1.43239\ldots\).

- **Priority.** The number is not new: DGH (1992) state it as the di-block collapse temperature 1.43. The collapse side was argued by DGH and proved by Poisat (Prop. 3.7 and Cor. 3.8).
- **What is new.** No earlier work I could find proves the matching high-temperature side for \(\beta<\beta_*\), and Poisat still lists the question as open in October 2026. The solver supplies that side, with a short and correct argument.
- **Status.** "solved", with the priority qualification about the value that the manuscript already states.

## 2. Source and interpretation (checked myself)

- **Verbatim quote.** I downloaded arXiv:2502.18964v2 myself; its SHA-256 matches `certificates/source_manifest.json`. I rendered p. 21 as an image. Item (3) of Section 7 reads, verbatim: "What is the critical temperature for the di-block polymer from (3.4)? The best estimate so far is \(\log(\frac{1+\sqrt5}{2})\le\beta_c^{\rm di-block}\le2\log2\). Approximate numerical values for these bounds are 0.481 and 1.386." The README and paper quote it correctly. In v1, the same question cites (2.18).
- **Which convention.** The bounds are in the **original** convention (2.3)–(2.5):
  - the lower bound \(\log\varphi\) is Prop. 3.4, which is stated for the original \(p_{i,n}^{\beta,\omega}\) and holds for all charge sequences;
  - the upper bound \(2\log2\) is a crude bound on \(2\beta_0\), because \(U(\log2)=\sum_t2^{-t^2}\approx0.564<1\). Poisat himself writes "the chain collapses at inverse temperatures \(\beta>2\beta_0\)" (p. 9).
  
  In the alternative convention, the interval would be \([0.24,0.69]\), not \([0.481,1.386]\). Hence \(\beta_*=2\beta_0\approx0.698\) lies inside the quoted interval, and the factor two (\(t_c=\beta_*/2\)) is right. I checked (2.13) numerically by brute force, using the ordered double sum including \(i=j\) (verification/independent_check.py).
- **The (3.4) misprint.** Poisat's display assigns +1 to \(1\le i<n\) and −1 to \(n\le i<2n\), so \(\omega_{2n}\) is undefined. His own (3.31)–(3.32) use the balanced word \(+^n-^n\): the mixed block \((i,n+j]\) has charge \(n-i-j\), and the fully collapsed state is neutral. DGH also use \(q_i=1\) for \(i\le N/2\) and −1 otherwise. The solver's reading is therefore correct. Section 8 of the paper covers bounded imbalance, and I re-derived (25).
- **Meaning of "critical temperature".** Poisat frames the di-block inside Conjecture 3.1(3): a high-temperature regime where all \(p_i\ge p_{\min}>0\), against a low-temperature regime. His lower bound is exactly such a positivity statement, and his "low-temperature regime" for the di-block is collapse (Prop. 3.7). The paper's (6) identifies one threshold for both notions. Independently, I checked that the di-block free energy (4) has a **first-order kink** at \(\beta_*\): \(F'_{\rm db}\) jumps from −0.2982 to 1/2. All three natural notions therefore give the same \(\beta_*\): uniform positivity of cuts, collapse, and non-analyticity of the free energy. The answer does not settle for a weaker variant.

## 3. What DGH (1992) and Derrida–Higgs (1994) actually did (read in full)

- **DGH 1992**, EPL 18, 361, pp. 364–365. They give the pure-chain eigenvalue equation \(\lambda=\sum_n\lambda^{-n}e^{(n-n^2)/2T}\). They then state: "One can show that a chain with \(q_i=1\) for \(i\le N/2\) and −1 for \(i>N/2\) will collapse at a temperature \(T_c^{\rm pure}=1.43\) determined by the equation \(\lambda=\exp[1/2T]\) since below \(T_c^{\rm pure}\) the free energy of the completely collapsed state is lower than \(-T\log\lambda\)."
  - Substituting \(\lambda=e^{1/2T}\) gives exactly \(\sum_ne^{-n^2/2T}=1\).
  - This is a one-line free-energy comparison for the collapse side. No uniform high-temperature statement is made for the di-block: their only high-temperature bound is the general \(T_c\le2.08\), which corresponds to \(\beta<\log\varphi\).
  - DGH thus **predicted** the value and sketched collapse; they did not prove the transition.
- **Derrida–Higgs 1994**, J. Phys. A 27, 5485. Zero-temperature scaling for Gaussian charges with fixed total charge. Nothing on the di-block threshold.

## 4. Proof review, line by line (re-derived myself)

| Step | Paper | Check |
|---|---|---|
| Block energy | (8) | Correct: a block of length \(\ell\) and charge \(q\) contributes \((q^2-\ell)/2\) for \(\pm1\) charges, so \(Z=e^{\beta N/2}Q\). |
| Homogeneous renewal | (9)–(10) | Correct. \(u_m=A_mr^m\) is a renewal sequence with probability kernel \(p_\ell=w_\ell r^\ell\), so \(u_m\le1\). The only zero on \(\lvert z\rvert=r\) is \(r\), because \(w_1>0\); it is simple. Cauchy estimates give a geometric rate, and the infimum is positive. \(\beta=0\) is handled separately. |
| Threshold | (11) | Correct: \(W_\beta(1)\) is strictly decreasing from \(\infty\) to 0. |
| Domination and gluing | (12) | Correct: \(e^{-\beta q_B^2/2}\ge e^{-\beta\ell_B^2/2}\), and compositions of \(uv\) with a cut at the junction factorize. Poisat's Lemma 3.5 is the case of lengths 1 and 2 (I added a remark saying so). |
| No-cut decomposition | (13) | Correct: exactly one block covers both \(i\) and \(i+1\); \(k\le i-1\), \(j\ge i+1\); the endpoint cases use \(Q(\emptyset)=1\). |
| Odds bound | (14)–(15) | Correct: \(Q(0,k]/Q(0,i]\le1/Q(k,i]\le1/A_{i-k}\), and symmetrically on the right. \(\sum_a1/A_a\le r/(c(1-r))<\infty\) iff \(r<1\). The bound holds for every sequence, every \(N\) and every \(i\). The "elementary alternative" (truncated kernel, convex-combination induction) is also correct. |
| Di-block identity | (16)–(17) | Correct, with \((a-b)^2=(i-j)^2\) and \(D_n\ge1\) from the neutral collapsed state. I checked it as an exact integer-polynomial identity for \(2n\le40\) using a site-charge transfer matrix (not the composition recursion). |
| Collapse | (3) | Correct: \(\sum_{i+j\ge M}A_iA_j\le s^{-M}G(s)^2\) for \(1<s<r\), and the interface-cut term satisfies \(A_n^2\le s^{-2n}G^2\). A cut at position \(i\) forces \(2n-L_n\ge\min(i,2n-i)\). |
| Non-collapse below \(\beta_*\) | §5, end | Correct: the number of cuts is at most \(2n-L_n\), and bounded convergence finishes the argument. |
| Free energy | (18)–(21) | Correct. Using \(\sum_{i<n}r^{-i}\le r^{-n}\cdot r/(1-r)\), the at-criticality bound \(1+n\Theta\), and the \(\beta\to\beta/2\) map with \(-\beta-\log2\) shift from (2.15). |
| Critical asymptotics | (22), (5) | Correct: diagonal Cesàro averages tend to \(a^2\), dominated by \(e^{-\beta d^2/2}\); \(\Theta_{\beta_*}=1+2W_{\beta_*}(1)=3\). |
| Off-critical asymptotics | (23)–(24) | Correct: domination by \(c^{-2}r^{a+b}\), respectively monotone convergence. |
| Endpoint imbalance | (25) | Correct. |
| Numerical bounds | (26)–(28) | Correct. The ratio bound for the Gaussian tail holds. The resolvent bound \(0\le a_m-b_m\le(m+1)T\) holds because all scaled renewal coefficients are \(\le1\). The diagonal truncation bound holds. |

No gap found. The theorem as stated (items 1–4 and (6)) is proved.

## 5. What I recomputed (`verification/independent_check.py`, log `verification/independent_check.log`, all PASS, about 7 s, one thread)

1. **Two independent methods agree.** Brute-force enumeration of the original pair Hamiltonian \(\sum_{i<j}\omega_i\omega_j1\{s_i=s_j\}\) matches a site-charge transfer matrix exactly (integer polynomials), for \(N\le14\), on di-blocks and random words.
2. **Convention.** Poisat's (2.13) was checked from the alternative definition directly.
3. **Identity (16)** holds as an exact integer-polynomial identity for all \(2n\le40\).
4. **The root.**
   - mpmath (60 digits) and the Jacobi theta formulation \(\vartheta_3(0,e^{-\beta/2})=3\) give the same \(\beta_*\).
   - My own Arb bracket with my own tail bound gives \(0.698131700799199439091341370180112797008487712347<\beta_*<\ldots348\), consistent with the paper's 50-digit enclosure (7).
   - Arb proves \(\beta_*>2\pi/9\).
   - **New analytic explanation.** Poisson summation gives \(\Theta_{2\pi/9}=3(1+2e^{-9\pi}+\cdots)>3\). Hence \(\beta_*>2\pi/9\) analytically, and \(\beta_*-2\pi/9=\frac{8\pi}{9}e^{-9\pi}(1+O(e^{-9\pi}))\). The leading term \(1.46760832178\times10^{-12}\) matches the true difference \(1.46760832187\times10^{-12}\) to a relative \(6\times10^{-11}\).
5. **Uniform lower bound (2).** I used my own full-chain forward/backward DP in log space.
   - For \(\beta\in\{0.3,0.5,0.65,0.69\}\), every cut probability on homogeneous, di-block, periodic and random words (\(N=50,400\)) is \(\ge\delta_\beta\) of (15). At \(\beta=0.69\), \(\delta=2.1\times10^{-5}\) and the observed minimum is \(4.7\times10^{-3}\), so the bound is loose but valid.
   - The odds inequality (14) holds for **all 4096 words of length 12** and all cuts, at \(\beta=0.69\) and \(\beta=0.75\).
6. **Di-block.**
   - The full-chain DP agrees with \(A_n^2/D_n\) to \(10^{-9}\) at \(2n=600\), at three temperatures.
   - At \(\beta_*\): \(3n\,p_{n,2n}=1.000156,\ 1.0000156,\ 1.00000156\) for \(n=10^3,10^4,10^5\). This agrees with (5), with a \(1.56/n\) correction.
   - Below \(\beta_*\), \(p_{n,2n}\to1/K_\beta\) (agreement to \(10^{-6}\) at \(2n=3000\)).
   - Above \(\beta_*\), cut probabilities decay exponentially away from the ends.
   - The finite-volume free energies approach (4).

## 6. Solver scripts re-run (`verification/solver_scripts_rerun.log`; also `run_all.py` on a scratch copy, exit 0, 24.7 s)

All three scripts print PASS (107 PASS lines, no FAIL). What they actually test:

- **`verify_exact.py`.** Identity (16) against exhaustive enumeration, as integer polynomials, for \(N\le16\). For all words of length \(\le8\) at \(q=4/5\): recurrence = enumeration of \(\sum Q_r^2\), domination, gluing, and the exact cut decomposition (13). It does **not** test the odds bound (14) or the pair-Hamiltonian normalization.
- **`verify_independent.py`.** Despite its name, this script was written by the solver's internal-review sub-agent. It tests the normalization and (14) for lengths \(\le7\), and gives a coarse rational root bracket.
- **`verify_intervals.py`.**
  - Certified root bracket; separation from \(2\pi/9\).
  - Arb DP up to \(N=100000\), with rigorous truncation balls. Stored balls are checked to contain the recomputed balls.
  - A finite-volume check that \(3n\,p\in(0.999,1.001)\) at \(n=50000\).
  
  These are finite-volume certificates. The paper correctly does not present them as proof of the limits.
- **Displayed values.** The README table values agree with my re-run.

## 7. Novelty search (18 Oct 2026)

- **OpenAlex.** Poisat (W4414926404) has 0 citing works. DGH 1992 (W1965552398) has 26 citing works, all listed; the latest is from 2016. They cover undirected annealed or quenched charged polymers (Chen, Khoshnevisan, Hu, Wouts, Caravenna–den Hollander–Pétrélis–Poisat), directed polymers with random interaction (Mukherji–Bhattacharjee), polyampholytes, and Derrida–Higgs 1994. None treats the di-block threshold.
- **Semantic Scholar.** Poisat has 0 citations. The DGH record is barely indexed (1 citation; one HTTP 429 before success).
- **arXiv API, sorted by date.**
  - Author queries: Poisat, Pétrélis, Caravenna, den Hollander. Poisat's only charged-polymer item is 2502.18964 itself; nothing newer on this model from any of the four.
  - Keyword queries: "charged polymer", "diblock AND charged", "polyampholyte". Nothing relevant.
- **Poisat's publication page.** No other 2024–2026 charged-polymer item.
- **Poisat's seminar abstracts.** Lille, and I2M on 29 Sep 2026: a review of the freezing conjecture; no di-block result announced.
- **GitHub.**
  - `mmislan`: 23 repos (CRNT / RAF); none on charged polymers.
  - `ipitchford`: 166 repos; the keyword grep (polymer, charge, Poisat, Derrida, diblock, freezing, renewal) hits only two unrelated spin-chain and Hartree repos.
  - `gh search repos` on five key phrases: no results.
- **Zenodo.** Key phrases, "2502.18964", "weak freezing", "Poisat": nothing relevant.
- **Web searches.** Four queries: nothing beyond the primary sources.

**Conclusion.** No prior complete answer was found. The priority qualification (value and collapse side from DGH and Poisat) is stated correctly in the paper and README.

## 8. Corrections applied (diff in `verification/verifier-changes.diff`)

1. **`paper.tex`.**
   - `\PaperCredits` and subtitle updated: independent verification recorded.
   - Remark after (7): theta/Poisson-summation proof that \(\beta_*>2\pi/9\), with the asymptotic difference; \(T_c=1/\beta_*\); \(\beta_*=2\beta_0\).
   - Remark after Section 4: Lemma 3.5 is the length-1/2 case of (12). Since (2) holds for all sequences, \(\mathcal F_n^{\beta,\omega}(p)=0\) for \(p<\delta_\beta\), so any random-charge freezing point satisfies \(\beta_c\ge\beta_*\). This is explicitly *not* a solution of question (2).
   - Enclosure (7) re-typeset: it overflowed the page by about 195 pt. One section title shortened.
   - Section 9: provenance of `verify_independent.py`.
   - Section 10: last bullet updated.
   - Recompiled: `paper.pdf`, 8 pages.
2. **`README.md`.** Status changed to "solved" (independently verified), attribution updated, theta remark added. The file list now includes `verification.md` and `verification/`, and clarifies that `verify_independent.py` and `certificates/independent-verification.log` belong to the internal review.
3. **`references.md`, `scripts/verify_intervals.py`.** Dangling references to a non-shipped `notes.md` now point to `paper.tex`. The runtime in the docstring was corrected from "1–3 minutes" to the measured 15–20 s.
4. **`scripts/verify_independent.py`.** The docstring said "added by verifier"; it now names the internal-review sub-agent.

No mathematical statement of the solver needed correction.

**Disclosure.** After my analysis and computations were complete, a `grep notes.md` over the folder printed three single lines from forbidden files. One was the internal review's one-line verdict ("The theorem in `notes.md` is correct ..."). I did not open those files, and this did not influence the checks above.

## 9. Reservations

- **Interpretation.** Poisat did not define \(\beta_c^{\rm di-block}\) formally, and his upper bound \(2\log2\) instead of \(2\beta_0\) is slightly puzzling. However, every natural reading gives the same \(\beta_*\) (§2): cut positivity as in Conjecture 3.1(3), collapse, and free-energy non-analyticity. The answer is an implicit transcendental number, a root of \(\vartheta_3(0,e^{-\beta/2})=3\); no closed form is claimed.
- **Weight of the contribution.** The decisive new step, (12)+(14), is short and natural: it extends Poisat's Lemma 3.5 from lengths 1–2 to all lengths. A referee may view the result as a modest, if complete, resolution.
- **Proof status.** The proof is a human-readable argument, not formalized. The finite certificates check identities and finite volumes only.
- **Novelty search coverage.** The search cannot exclude unindexed work, theses, or talks (e.g. a remark in a lecture). Given Poisat's October 2026 listing, prior resolution is unlikely.
- **Root README.** The results table row 4 still says "solved (preliminary)". It is outside this folder, so I did not edit it. The owner should update it to "solved" and add the verifier credit.
