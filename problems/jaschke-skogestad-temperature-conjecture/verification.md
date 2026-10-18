Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Independent verification: Jäschke–Skogestad temperature conjecture

**Verdict: CONFIRMED WITH CORRECTIONS.**

The mathematical result is correct.
- The published recursion (18)/(A.15) is wrong for every \(M_j \ge 3\).
- The corrected cumulative recursion is proved for every \(M_j\).
- The \(M_j = 3\) proof of the paper contains a real algebraic error, at the step the solver identified.

The corrections concern how strongly some claims were worded (Appendix B, the \(M \ge 4\) argument), missing nuance,
and the manuscript's presentation. I did not read `PROGRESS.md`, and I launched no sub-agent.

Date: 18 October 2026.

## 1. Source, as read by the verifier

I downloaded the publisher PDF from Skogestad's site myself (`publications/2014/jaschke-hen-paper-cce-morari-issue/impress.pdf`, 15 pages, CCE 70 (2014) 35–49) and read the full text. I also checked the formulas on page images: pp. 38–39 at 130 dpi, pp. 45–48 at 130 dpi, and (A.12), (B.1), (9)–(10) cropped at 300–400 dpi.

**Transcription by the solver: correct.**
- Model (9), (10), (12): counter-current AMTD.
- Cost (14): \(J=-\sum p_{ij}Q_{ij}\).
- Theorem 2 / (16)–(18) and Theorem 3 / (19)–(20).
- The recursion \(a_{i,j}=(\theta_{i,j}-\theta_{i-1,j})(\theta_{i,j}+\theta_{i-1,j}-a_{i-1,j})/(\theta_{hi,j}-\theta_{i-1,j})\).
- The verbatim conjecture sentences on p. 39, A.3, A.4 and the Conclusions.

**The recursion subtracts only \(a_{i-1}\), unambiguously.** (18), the efficiency form (29) in §6.1, (A.13), (A.15) and the fully expanded (B.1) all write it that way.

**The printed (B.1) contains misprints.**
- \(\theta^2_{h1,1}/\theta_{h1,1}\) appears in the \(p_{1,1}\) term and inside the \(p_{4,1}\) term.
- \(\theta_{h3,1}-\theta_{2,1}\) appears instead of \(\theta_{h2,1}-\theta_{1,1}\) inside the \(p_{3,1}\) term.

**Additional misprints noticed by the verifier (harmless).**
- The printed (9) has \(w_i\) for \(w_j\).
- The printed (10) has \(T^{out}_{i-1,j}\) for \(T^{out}_{hi,j}\).
- (A.8) and (A.10), presented as the \(M_1=2\) computation, contain \(\theta_{h3,1}\), which does not exist in an \(M_1=2\) network. The factor \((\theta_{h3,1}-\theta_{2,1})\) cancels, so the printed raw expressions were apparently taken from the 3-exchanger run. This does not affect anything.

## 2. What I recomputed myself (scripts in `verification/`, no code shared with `scripts/`)

**Derivation by hand, before reading the solver's proof.**
- Under AMTD, \(Q_i = w\,e_i(\theta_{h,i}-\theta_{i-1})\) with \(e_i = 2/(1+w(2/UA_i+1/w_{h,i}))\).
- This gives \(w\,\partial e_i/\partial w=-e_i(2-e_i)/2\).
- Then \(\partial Q_i/\partial w = c_i\Delta_i/2 - c_i\,w\,\partial\theta_{i-1}/\partial w\).
- The cumulative balance \(\sum_{l<i}Q_l = w\theta_{i-1}\) then gives \(2\,\partial Q_i/\partial w = c_i(\theta_i+\theta_{i-1}-\sum_{l<i}2\,\partial Q_l/\partial w)\).
- This is the solver's corrected recursion, with the same factor 2 and the same cumulative sum.

**Line-by-line check of the solver's proof (paper.tex §3): no error found.**
- Lemma 1, including the factor 2 and the positivity of the bracket.
- Lemma 2.
- Steps 1–3 of Proposition C.
- The AMTD specialisation.
- The recursion \(d_i=c_i(A_{i-2}-d_{i-1})\).

**Scripts and results.**

| Script | Method | Result |
|---|---|---|
| `v1_symbolic.py` (log `v1_symbolic.log`, ~5 min) | Each exchanger's 3 equations solved by `sympy.solve`, derivative by `sympy.diff`, identities by `cancel`. | Corrected recursion \(2\,\partial Q_i/\partial w=a_i\) is an identity for **M = 1, 2, 3 with finite hot capacities** (all symbols generic) and **M = 1–4 with constant hot temperatures**. Published recursion exact for \(i\le2\) and wrong at \(i=3\); branch totals differ; \(a^P_3-a_3=c_3a_1\) identically. Single exchanger: \(\partial J/\partial w=-p\theta_1^2/(2\theta_{h1})\), and the true \(\theta_1=2\theta_{h1}/(1+w(1/w_h+2/UA))\), so the paper's (A.5) is wrong. ALL PASS. |
| `v2_exact.py` (log `v2_exact.log`, seconds) | Exact dual numbers over \(\mathbb{Q}\); each exchanger solved as a raw 3×3 linear system, without the affine lemma. | **R:** 320 random exact draws, M = 1–8, finite and infinite hot capacities mixed: corrected recursion exact in every draw; published branch total wrong in every draw with M ≥ 3. **X1:** `certificates/counterexample.json` re-checked from its raw parameters only. Physical point (heating, no temperature cross, \(0<c_i<1\)), paper CV exactly 0, gradient exactly \(-4972391424/1518467647\), equal to \((T^*_{J,1}-T^*_{J,2})/2\); \(Q_{tot}(99/2)>Q_{tot}(50)\); stored temperatures agree. **X2, X3:** my own exact physical counterexamples (paper CV = 0, gradient ≠ 0): \(M_1=3\) with constant hot temperatures, the case the paper reports as proved (gradient −20/11); \(M_1=4\); \(M_1=5\). ALL PASS. |
| `v3_paper_formulas.py` | My own transcription of the raw Maple coefficients printed in (A.8), (A.10), (A.12). | (A.8) gives \(a_1\); (A.10) gives \(a_2\). **(A.12) line 1 equals the corrected \(a_3\)**. (A.12) line 2 (the "simplified" form) equals the published \(a_3\). Their difference is \(c_3a_1\). Each of the 12 possible single-term sign flips of the numerator destroys the identity, so the match is not a transcription accident. My transcription coincides term by term with the solver's. ALL PASS. |
| `v4_appendixB.py` | Table B.6 data, exact rational functions of \(u\), exact root isolation. | **B1:** corrected CV \(=-(2/w_0)\,dJ/du\) identically. **B2:** unique optimum \(u_{opt}=0.7141194665\), \(J=-3938.8167\), \(T_{end}=158.77375\), which matches the paper. **B3:** paper CV root \(u_P=0.7167062190\), where \(dJ/du\neq0\) (exact root counting): an \(M_1=4\) counterexample to Theorem 3 on the paper's own data. **B4:** all prices 1: optimum 0.72267, paper CV root 0.72505 (counterexample to Theorem 2, \(M_1=4\)). ALL PASS. |

**Exploratory output of `v4_appendixB.py` (no PASS/FAIL).**
- The literal printed (B.1), with its misprints, has a single root in (0,1), at \(u=0.90394\).
- At the paper's reported \(u=0.7152\), the AMTD model gives \(T_{end}=158.7741\), the paper's value.
- But at \(u=0.7152\) neither CV vanishes: the paper's CV is 0.168, the corrected CV is −0.119.
- The LMTD optimum is 0.7014.
- So the reported 0.7152 is reproduced by no formula I tried.

## 3. Re-run of the solver's scripts

All four were re-run under WSL. Logs are in `verification/rerun_*.log`.

| Script | Result | What it really checks |
|---|---|---|
| `verify_symbolic.py` | ALL PASS, 81 PASS lines, 54.6 s, 225 MB | Cramer solution of the raw model over \(\mathbb{Z}[\ldots]\) fractions with python-flint. S1 (corrected recursion) is a true identity for M ≤ 5 with finite \(w_h\) and M ≤ 8 with \(w_h=\infty\). L1 transfers \(w_h=\infty\) to finite \(w_h\); this is correct, because \(UA'\) is independent of \(w\). S2/S3 non-vanishing is certified by exact evaluation at a single rational point, which is valid (a nonzero value implies a nonzero rational function). S4 is the (A.12) check. |
| `counterexample_exact.py` | ALL PASS | A **generator plus checker**: it rebuilds the counterexample from hard-coded parameters and **rewrites** the JSON passed as argument; it does not read it. The output is deterministic and the content was unchanged. My `v2_exact.py` supplies the missing "read and check the certificate" step. Exact dual numbers on the raw 3×3 systems; the optimum location is numerical (labelled as such). |
| `reproduce_appendixB.py` | ALL PASS | Exact root isolation; agrees with my `v4`. |
| `random_fd_tests.py` | ALL PASS | 60-digit finite differences. Numerical evidence only, correctly labelled. Part C (LMTD form of Proposition C) is numerical only, but Proposition C has a general proof that I checked (§2). |

## 4. Judgement of each claim

1. **"The published recursion is false for every \(M_j\ge3\); the correct one is \(a_i=c_i(\theta_i+\theta_{i-1}-\sum_{l<i}a_l)\); short general proof."** CONFIRMED.
   - Proof checked line by line, and re-derived independently.
   - Symbolic identities and exact random tests agree.
2. **"The published proof of \(M_j=3\) loses a term in a manual simplification of (A.12)."** CONFIRMED.
   - The raw first line of (A.12) is exactly the corrected \(a_3\). The second line, joined to it by "=", is the published \(a_3\), and the two differ by \(c_3a_1\neq0\).
   - The wording "error in the published proof of \(M_j=3\)" is accurate.
   - Precision added: the paper's \(M_j=3\) proof was only for constant hot temperatures. The error makes that "proved" case false as well. My X2 is an exact physical counterexample with constant hot temperatures.
3. **Exact rational counterexample (\(M_1=3\), \(M_2=1\)), gradient \(-4972391424/1518467647\).** CONFIRMED independently from the JSON parameters.
   - The operating point is physical.
   - Minor wording fix: the largest terminal-difference ratio is 1.64. That is outside the paper's ±40 % (\(\sqrt2\)) criterion, and AMTD exceeds LMTD there by about 2 %. "Reasonable approximation" was replaced by this figure; it is irrelevant to the claim, since the model is AMTD by assumption.
4. **Appendix B: "the discrepancy attributed to round-off is in fact the recursion error"; factor −1/2 in (16)/(19).**
   - The factor −1/2 (sign from (A.6), factor 2 from the mis-solved (A.5)) is CONFIRMED. As the solver says, it is a common normalisation that does not affect the controlled variable. I added the explicit statement that, taken literally, (16)/(19) are false even for \(M_j=1,2\).
   - The Appendix B statement was **too strong** and is **corrected**. The paper's CV, evaluated exactly on the paper's own AMTD model, vanishes at 0.71671, not at the optimum 0.71412. So the recursion error alone produces a split error of 0.0026, which cannot be blamed on round-off.
   - However, the paper's reported 0.7152 is reproduced by neither recursion. Its \(T_{end}\) is consistent with an AMTD evaluation at \(u=0.7152\), but no CV vanishes there. One therefore cannot say what the authors actually computed, nor that the reported discrepancy "is" the recursion error.
   - Also, the paper attributes the **end-temperature** discrepancy (not explicitly the split) to "flatness of the optimum and numerical round-off errors".
   - README, abstract and §3.5 now say exactly this.
5. **"The conjecture is false."** CONFIRMED, with a nuance added.
   - **What is false:** the conjectured formula, both for \(M_j\ge4\) and for \(M_j=3\) with finite hot capacities. It is false even up to any common normalisation, because the paper's CV vanishes where the gradient does not. Direct counterexamples exist for M = 3, 4, 5 (exact), and the paper's own data give one for M = 4.
   - **What is true:** the paper's qualitative principle (under AMTD, temperature measurements suffice) holds for every \(M_j\), and is now proved by Theorem A. This sentence was added to the README and to the paper.
   - **Gap in the general \(M\ge4\) argument, fixed.** Corollary B, item 3, appended small-\(\varepsilon\) exchangers but did not say that branch 2's hot inlet must be re-tuned so that the paper's CV is again exactly zero. I added the general construction (\(\theta_{h1,2}=T^P_{J,1}/c^2\)) and pointers to the direct \(M=4,5\) instances.

## 5. Novelty (quick check by the verifier)

- **Crossref** (`api.crossref.org/works/10.1016/j.compchemeng.2014.03.020`): `relation` is empty, and there is no `update-to` / `updated-by`. No erratum is registered.
- **OpenAlex:** 36 citing works. I listed all 36 titles. Those touching parallel-branch HEN operation are:
  - Jäschke & Skogestad, ESCAPE-24 2014 (≤ 2 exchangers per branch, so unaffected);
  - Trafczyński et al., CET 70 (2018). I read the full text: a real crude preheat network with 26 exchangers in four parallel branches, controlled with \(JT_A-JT_B\) and \(JT_C-JT_D\). No formula is written out, and the network also has hot-stream interactions between branches;
  - Mohanan & Jogwar (2022–2023): a different approach;
  - Bernardino, Krishnamoorthy & Skogestad (2022);
  - Ikarimoto MSc thesis (2019).
- **Web searches** ("Jäschke temperature" with erratum / correction / three exchangers) found only the original paper, its preprints, theses and the patent page.
- **Conclusion:** consistent with the solver's `references.md`. No published correction or proof for \(M_j\ge3\) was found.
- **Not done by me:** reading the theses (Aaltvedt 2013, Ikarimoto 2019) or the patent. The solver's statements about them come from its sub-agent ([sub] in `references.md`) and remain unverified by me.

## 6. Corrections applied by the verifier

**`paper.tex` / `paper.pdf`** (recompiled: 10 pages, no warnings).

*Presentation.* The previous PDF printed every "- " bullet list and every "1. 2. 3." list inline as running text, and doubled the subsection numbers ("3.1 3.1 Lemma 1"). These were rewritten as `itemize`/`enumerate`, with automatic numbering. Fixes:
- removed an internal sentence referring to "the scout's summary (rapports/scout_E_mining.md)";
- added a `references.bib` (the protocol requires one, and it was missing), with `\autocite` and a bibliography;
- removed overfull lines (`\url`, a displayed fraction);
- credits and subtitle updated ("independently verified").

*Content.*
- Abstract reworded:
  - the conjecture is false, and so is the M = 3 case reported as proved;
  - the error is located;
  - the Appendix B wording is softened.
- §1: noted the printed misprints in (9)–(10), and added (29) as a further witness of the recursion.
- After Theorem A: added the explicit statement that the qualitative claim of the paper's title holds for every \(M_j\).
- Corollary B, item 3: general construction with re-tuned \(\theta_{h1,2}\); \(M\ge4\) argument completed; direct instances cited.
- Item 4: precision on (A.13).
- §3.4: the role of `counterexample_exact.py` (it regenerates the JSON) is stated, and the AMTD/LMTD remark is replaced by the 2 % figure.
- §3.5: Appendix B rewritten as in §4 above. It now includes:
  - the price-1 instance B4;
  - the root of the literal misprinted (B.1) (0.904);
  - \(T_{end}(0.7152)\);
  - the CV values at 0.7152.
- §4: new subsections "Independent verification" and "Literature".
- §4.4: (16)/(19) are literally false for every \(M_j\).

**`README.md`.**
- Attribution and status: "solved (independently verified)".
- Result items 2–4 reworded as above.
- Reproduce section: notes that `counterexample_exact.py` rewrites the JSON, and adds the verifier's commands.
- File table updated.

**`references.bib`:** new file.

Not modified: `scripts/`, `certificates/`, `references.md`, `PROGRESS.md`. The original versions of `paper.tex`, `paper.pdf` and `README.md` were kept in the verifier's scratchpad only. **The results table of the root `README.md` is outside this folder and was not updated:** the owner should mark this problem "solved (verified)" and add the verifier credit there.

## 7. Remaining reservations

- **Exact symbolic coverage of the corrected formula.**
  - The general proof (half a page, checked) covers every \(M\).
  - The machine identities cover M ≤ 5 with finite hot capacities and M ≤ 8 with constant hot temperatures (solver), and M ≤ 3 / M ≤ 4 respectively (verifier).
  - The finite-\(w_h\) cases up to 8 follow from L1. Larger M rests on the proof alone, which is fine.
- **The origin of the paper's reported 0.7152 is unknown** (see §4, item 4).
- **The LMTD statements are side results.** They are not needed for the disproof.
  - Proposition C is proved in general.
  - The "median 9 % error" figure is numerical evidence only.
  - The question whether an exact LMTD formula without hot outlet temperatures exists is left open, as the paper says.
- **The model is the paper's, and nothing more.**
  - Ideal AMTD counter-current exchangers, constant \(c_p\), hot-side data independent of the split (no crossovers).
  - The \(M_j\ge3\) error therefore matters for plants with three or more exchangers in series per branch, such as the 26-exchanger network of Trafczyński et al.
  - The practical losses are small (flat optimum: 0.066 % in the counterexample), but the split error is first order.
