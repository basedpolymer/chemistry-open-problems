# Application of the reviews R1 and R2 to the manuscript of result 17

Applied by R1 (Claude Opus 5.5, effort xhigh), 28 October 2026, at the coordinator's request (the protocol has the
verifier apply corrections). Files changed in the problem folder: `paper.tex`, `paper.pdf` (rebuilt),
`references.bib`, `references.md`, `scripts/certify_duality.py` (one claim added), `certificates/certify_duality.log`
(rewritten by the script). Nothing else outside `gac-duality/review/R1/` was edited; nothing was committed.

Mechanics. `scripts/apply_both.py` imports R2's replacement list from `../R2/scripts/apply_corrections.py`, replaces
four R2 entries by merged versions, adds R1's replacements and proof insertions, checks that every old text occurs
exactly once in the original and that no two old texts overlap, and writes `paper.tex` and `references.bib`
(`logs/apply_both.log`: 47 R2 + 42 R1 replacements in `paper.tex`, 9 in `references.bib`, 3 entries added, 0
problems). `scripts/apply_refs_md.py` applies M01–M10 and one follow-up to `references.md` (`logs/apply_refs_md.log`,
11/11). Two presentation fixes found while building were then made by hand (§3). Line endings of each file were kept
(LF for `paper.tex`, `references.md`; CRLF for `references.bib`).

Precedence (coordinator): R1 on mathematics, R2 on citations and wording; nothing of either was dropped silently.

## 1. R1's report (`review.md`)

| id | content | status |
|---|---|---|
| R-1 | abstract: cluster claims only for equal-conflict networks | **merged** into R2's C01 (same content; R2's wording) |
| R-2 | abstract: periodic orbit "for some rate vector" | applied |
| R-3 | twin pairs ↔ doubling pairs (zero rows), twin classes of nonzero rows ↔ doubling classes; Table 6 row | applied (R1-R3a, R1-R3b) |
| R-4 | q = m needs nonzero sources | applied in a **stronger form** after R2 §9 item 2 (a species in no reaction is also a cluster): "for connected EC networks with (H1) and (H2)"; the counterexample was replaced by a connected one, 0 → X₁, 0 → X₂, X₁+X₂ → 0 (m = 2, q = 3), certified in `scripts/r1_exact.py` |
| R-5 | proof of Thm 68: factor c in (log V)˙ | applied; the last-entry argument was also written out (R1-R5b) |
| R-6 | "17–35 %" → 11–33 % | **merged** with R2's C33 (R2's wording plus R1's note on the one-network family (4,3,4)) |
| R-7 | Example 62, network Q: "κ₂ ≠ κ₃ (not analysed on the hypersurface)", also fails for WR pure networks | applied (R1-R7a, b) |
| R-8 | Thm 63: λ in the simplex | applied |
| R-9 | §6.3 standing hypothesis (H1)–(H3) | applied |
| R-10 | sketch of Thm 65: isolation sentence | applied |
| R-11 | Prop. 51: net production | applied (R1-R11a, b) |
| R-12 | Example 67: reaction order | applied |
| R-13 | Prop. 41: N1* steady state | applied |
| R-14 | NC paragraph: degree argument cited | applied (with "an exact Gröbner-basis count of reviewer R1 agrees") |
| R-15 | Conjecture F paragraph: who proved what | applied |
| R-16 | §9.2 `\zcref` slip; Thm 6.17 sentence | **merged**: R2's C40 (identical) and C41 (R2's fuller wording) applied instead |
| R-17 | Prop. 66 sketch label and pointer; §9.2 sentence on sketches | applied (proof sketch inserted; R1-R17b); the path in the pointer was replaced by a reference to the concordance table to remove an 18 pt overfull box |
| R-18 | pointers for Thm 70 and Prop. 73 | applied (R1-R18a, b) |
| R-19 | Figure 46 caption: unbalanced joins | applied |
| R-20 | proof of Prop. 48 | applied; also "the source rows of Ye" → "the rows of Ye" in the statement (R1-R20b) |
| R-21 | proof of Cor. 61 | applied |
| R-21b | sources of the two "smallest" claims | applied (R1-R21b1, b2) |
| R-21c | §9.1 certificate coverage | **preferred fix applied**: the converse example of Prop. 49 added to `certify_duality.py` (claim "S5 converse (agent E)"), 129 → 130 in the abstract (inside the merged C03) and §9.1 (R1-cert1), and "the examples of Prop. 49 and of its converse" added to the list (R1-cert2); R2's C38 wording used for the sentence; R1's minimal fix not needed |
| R-22 | proofs or pointers for Props 23, 71, 74, Cors 32, 34, 43, 59, Lemmas 57, 64, Thm 47 | all applied |

Additional insertions for R2 §6 (statements without proof that R1's report had not covered): Lemma 28 (proof of the
Perron–Frobenius lemma, as re-derived by V2), Prop. 25 and Prop. 49 (pointer to the certificate claims). Lemma 2 and
Prop. 8 already cite the archived report in their headers; left as they are.

## 2. R2's report (`../R2/review.md`)

| ids | status |
|---|---|
| C01 | applied (also covers R-1) |
| C02 | applied |
| C03 | **merged**: R2's text with "130 claims" |
| C04 | **merged**: "(e.g. for connected networks with (H1) and (H2))" instead of "(e.g. under (H1) and (H2))", because under (H1) and (H2) alone a reaction 0 → 0 has the zero source (R-4) |
| C05, C06, C07, C08, C09, C10, C11, C12, C13, C14 | applied |
| C15, C16, C17, C18, C19, C20, C21, C22, C23, C24 | applied |
| C25 | applied (checked: the hypotheses of Desel–Esparza Thm 10.16 are invariant under the reverse dual) |
| C26 | applied |
| C27 | **merged**: Prop. 3.1 verified by R1 in the arXiv v1 PDF (statement: for all c\*, c′ the set (c′+S)≥ ∩ {c\*∘e^ṽ : ṽ ∈ S̃^⊥} has at most one element iff σ(S) ∩ σ(S̃^⊥) = {0}; Prop. 3.9 is the isomorphism statement). The journal numbering is unchecked, so the paper cites "Prop. 3.1 of the arXiv version" and the bibliography entry gains `eprint = 1209.6488` (R1-B10). R2's reading of the use (§9 item 5) is the intended one: with S̃ = D and u(k) ranging over all of ℝᵐ_{>0}, Theorem 35(3) is condition 2 |
| C28, C29, C30, C31, C32 | applied |
| C33 | **merged** (also covers R-6) |
| C34 | applied (checked against D2 §2: edge flow = sum of the dual reactions whose source avoids Z; one-way/two-way/zero as in D2) |
| C35, C36, C37 | applied |
| C38 | applied (with the certificate claim of R-21c it is now accurate) |
| C39 | applied |
| C40 | applied (also covers R-16) |
| C41 | applied (also covers R-16) |
| C42, C43, C44 | applied |
| C45, C46, C47 | applied (C47: `\pdfsuppressptexinfo=-1`; see §4) |
| B01–B08 | applied |
| B09 | applied (FamiliPalsson2003, JamshidiPalsson2009, Duvall2024 appended) |
| M01 | **merged**: also records that R1 checked Prop. 3.1 in the arXiv PDF and that the paper therefore cites the arXiv version |
| M02–M10 | applied |
| (follow-up) | `references.md`, novelty list: "Blanchini–Giordano p. 6" → "p. 493" (R1-M11), for consistency with C29/M02 |

R2's §9 points for R1: (1) doubling pair/class with a zero source: settled by R-3. (2) q = m: settled by R-4 in the
stronger form above. (3) abstract outside EC: C01 applied; R1's counterexample (A+B → C, B+D → C, C → A+D, C → 2B, with
(H1)–(H3), minimal siphon {A,C,D} meeting the cluster {A,B,D} twice) shows the restriction is needed. (4) statements
without proof: all now have a proof, a sketch label with pointer, or a pointer (§1). (5) Müller–Regensburger: as C27.

## 3. Changes made while building (not in either report)

- The template prints the optional argument of `proof` in parentheses, so `\begin{proof}[Proof (sketch)]` came out as
  "Proof (Proof (sketch))." (already in the 38-page PDF). All five sketch proofs now use `\begin{proof}[sketch]`
  ("Proof (sketch).").
- The pointer inserted in Prop. 66's proof sketch named a long path and produced an 18.3 pt overfull box; it now refers
  to the concordance table.

## 4. Build, certificates, privacy

- Build from the repository root, `py -3 templates/build_tex.py problems/iosif-siphon-precluster-duality/paper.tex`:
  **41 pages**, no overfull box, no undefined reference or citation (`logs/build.log`). The PDF text contains no "??",
  no "zcref", and "independent" only in "linearly independent", "rate-independent" and "are independent and solvable"
  (equations).
- `scripts/certify_duality.py`: **130/130 PASS** on Windows (`py -3`, Python 3.14) and on WSL (Python 3.12); the
  committed log `certificates/certify_duality.log` is the WSL run ("claims 130: 130 PASS, 0 FAIL"). Copies:
  `logs/certify_duality_windows_after.log`, `logs/certify_duality_wsl_after.log`.
- `scripts/certify_examples.py`: **142/142 PASS** on Windows and WSL (`logs/certify_examples_*_after.log`).
- `scripts/r1_exact.py`: 62/62 PASS (one claim added for the connected q ≠ m example).
- **PDF paths.** `scripts/pdf_path_check.py` searches the raw bytes and every decompressed Flate stream (object streams
  included). Old `paper.pdf`: "Users", "PTEX.FileName" and the user-name fragment each found once (the figure's
  `/PTEX.FileName`). New `paper.pdf`: all three **absent** (the only raw match of "PTEX" is inside compressed binary
  data). The user-name fragment is passed on the command line and is stored in no file.
- Absolute paths with a user name: none in `paper.tex`, `references.bib`, `references.md`,
  `scripts/certify_duality.py`, `certificates/certify_duality.log`, or in any file of this folder (`logs/build.log`
  records the build output with the path shortened).

## 5. Numbering after the changes

No numbered item was added or removed (only proofs), so the numbers are unchanged (`logs/numbering_after.log`, same
as `logs/numbering.log`); pages from the rebuilt PDF.

| result | number | page |
|---|---|---|
| Iosif's dual is the weighted reverse dual | Proposition 3 (Figure 4 shows it) | 7 |
| transversality of minimal siphons | Theorem 20 | 14 |
| cluster rank equation: no critical siphon on either side, persistence | Theorem 31 (WF-EC case: Corollary 32) | 17 |
| proportional conflict: self-dual, existence, uniqueness, stability | Theorem 42 | 21 |
| NC3 (three steady states, one a saddle) | Proposition 44 | 22 |
| classes on which the duality transfers structure | Figure 46 | 24 |
| Theorem O (unique equilibrium per class, periodic orbit) | Theorem 54 | 26 |
| its numerical illustration | Figure 55 | 26 |
| rank one: exact persistence test | Theorem 58 | 27 |
| rank one: all-or-none dichotomy for WR networks | Theorem 60 | 28 |
| Theorem W (no dominated facet, any rank) | Theorem 63 | 29 |
| Theorem A (ALL dual with an OUT facet) | Theorem 68 | 32 |
| Theorem AF (the family N_{u,R}) | Theorem 70 | 32 |

Page count: **41** (was 38).
