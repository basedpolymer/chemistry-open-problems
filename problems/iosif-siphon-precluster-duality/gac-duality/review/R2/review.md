# Review R2 of the rewritten manuscript of result 17: attribution, consistency, presentation

Reviewer R2: **Claude Opus 5.5 (effort xhigh)**, fresh context, 27–28 October 2026. Brief:
`gac-duality/BRIEF-R.md` (scope R2), `PROTOCOL.md` §3–§5, `gac-duality/BRIEF-common.md`. Object: `paper.tex` /
`paper.pdf` of the problem folder as of 27 October 2026, 23:54 (38 pages). I did not take part in the work. I edited
nothing outside `gac-duality/review/R2/`. Paths below are relative to the problem folder
(`problems/iosif-siphon-precluster-duality/`). Line numbers are those of the current `paper.tex`.

**Method.** I read the paper in full, `references.bib`, `references.md`, `gac-duality/writer/CHANGES.md`,
`gac-duality/literature2/novelty.md`, the novelty verdict of agent L (`gac-duality/literature/notes.md` §10, §18,
§19) and the verdict and attribution sections of V1–V3. I compared every cited theorem/page number with the source on
disk (`certificates/sources/`, `gac-duality/literature/sources/`, `gac-duality/gac/certificates/sources/`,
`gac-duality/literature2/sources/`) or with the open-access text (Müller–Regensburger arXiv:1209.6488v1, Conradi et
al. PLoS 2017, both fetched to a scratch folder outside the repository). Every DOI and arXiv id was resolved through
Crossref, DataCite and the arXiv API (`scripts/check_ids*.py`). I compared every count of the paper with the note it
cites. I built a copy of the paper, applied **all** replacements below to a second copy
(`scripts/apply_corrections.py`: 47 + 8 replacements and 3 new entries, each old text found exactly once) and compiled it
(`logs/build_patched.log`). Scripts in `scripts/`, logs in `logs/`, exact arithmetic only.

## 0. Verdict summary

| # | Item | Verdict |
|---|---|---|
| 1 | Quotations of the source (Conj. 22, Rem. 23, Q. 24, GAC paragraph, abstract, Def. 19), page and version data | OK |
| 2 | History of the statement (HAL 2020 "Theorem", 2024 slides "Theorem (Ananiadi, I.)", 2025 poster, 2026 conjecture) | OK (poster: abstract only, not re-checked) |
| 3 | Petri-net citations (Hack; Desel–Esparza Def. 4.4, Thm 4.30(a), Prop. 5.4, Thm 6.14, Def. 6.19 and its footnote, Thm 6.21, Thms 10.16–10.17; Teruel–Silva §2, Thm 28, pp. 272/296; Jančar et al. §4; Mairesse–Nguyen Lemma 3.2, Cor. 3.2) | OK, one CORRECTION (10.16 alone is credited to Desel by the book: C24) |
| 4 | Reaction-network citations (ADLS Thm 2; Conradi–Flockerzi Fact 2, Alg. 1, Def. 3; thesis Def. 4.3, Alg. 4.4, Def. 4.6, Thm 4.17; PMDSC Cond. 3.1, Thm 3.8, Cor. 3.12; CFMW Thm 1; CDSS Thm 9; Anderson–Shiu Thm 4.6; Butler–Waltman Thm 4.1; ZGD Thm 3.10; Wiuf Thm 7.2) | OK |
| 5 | Müller–Regensburger "Prop. 3.9" for the sign-vector uniqueness criterion | ERROR: it is Prop. 3.1 in arXiv:1209.6488v1 (Prop. 3.9 is the Birch-type isomorphism under σ(S) = σ(S̃)); journal numbering not checked (C27) |
| 6 | Blanchini–Giordano "p. 6" | CORRECTION: the journal pages are 488–493; the passage is on p. 493 (C29) |
| 7 | Citations of the archived report: "Lemma 1" (transposition), "Prop. 9" (MIPS of the dual) | ERROR: they are Lemma 2 and Proposition 6 of `archive/paper-precise-forms.pdf` (Prop. 9 does not exist; Cor. 9 is another statement) (C16, C18, C20, C42, C43) |
| 8 | "Jamshidi–Palsson" named without a reference | GAP: add the reference (verified, C28, B09) |
| 9 | Prior art for networks with stoichiometric matrix ±Γᵀ (Famili–Palsson 2003, Jamshidi–Palsson 2009, Duvall 2024), listed by agent L as "must be cited, not claimed" (`literature/notes.md` §18) | GAP: not cited anywhere (C06, C15, B09) |
| 10 | Teruel–Silva's observation that EC is not self-reverse-dual, balanced-joins EC is, and their Duality Theorem (Cor. 32), listed by L, L2 (§5) and checked by V1 | GAP: not credited where the paper states the same facts (C17, C19) |
| 11 | Every DOI and arXiv id of `references.bib` (27 DOIs, 5 arXiv ids) | OK (all resolve, metadata match; ZGD lacks volume/pages; six entries lack a DOI that Crossref gives, B02, B03, B05–B07 and optionally Desel–Esparza) |
| 12 | `references.md` lists every source with what was read | OK (all 43 keys); CORRECTIONS M01–M10 to follow the changes above |
| 13 | Novelty wording ("new") | CORRECTION at three places (C21, C26, C39); everything else uses "not found" with §9.3 |
| 14 | Credits required by the brief (ZGD, Desel 1994 / DE 10.16–10.17, DE 4.30(a) and 5.4, PMDSC, Campos–Chiola–Silva, Júlvez–Recalde–Silva, BG p. 6, thesis Thm 4.17, Johnston) | OK, all present in the body; Table 1 omits Recalde/Júlvez for Thm 33 (C09); for ordinary nets the dual half of Thm 31 is also known (C25) |
| 15 | Abstract ↔ introduction ↔ body | CORRECTION: abstract states the cluster claims without "equal-conflict" (C01), the equality condition only for twins (C02), "every finite claim" certified (C03); the intro omits the hypotheses of the cluster correspondence that V1 asked for (C04) |
| 16 | Table 1 (status, verifier, novelty) | CORRECTION: sketches not marked (Thm 24, Prop. 73), Thm 70 missing, NC3 "verified" by its own finder, Thm 19's correction not marked, "Conj. F" undefined at that point (C07, C08, C11, C13, C14, C44) |
| 17 | Counts against their notes | ERROR: rank-one share "17–35 %" is 10.75–33 % in Q1's own table (C33, `logs/rank1_percent.log`); CORRECTION: "exhaustive families" for a set that includes two seeded samples, "tested sides" (C29), the unexplained "199" (C32). All other counts match (§5.3) |
| 18 | Open problems ↔ body; §9 ↔ body; concordance table | CORRECTION: open problem 6 names a notes label "Theorem C" (C37); §9.1 overclaims (C38); §9.2 lists a correction not present in the paper (C41) and has a broken `\zcref` (C40, ERROR, visible in the PDF as "zcrefthm:B"); concordance: archive numbers (C42, C43); all other notes labels exist |
| 19 | "independent" in the non-mathematical sense | CORRECTION: 6 occurrences (lines 66, 163, 1363 twice, 1438, 1620; C03, C05, C33, C36, C39). Mathematical uses (716, 831, 832, "rate-independent" 603) stay |
| 20 | Owner, personal names, absolute paths, private correspondence | `paper.tex`: OK (only the acknowledgment sentence; "in his words" quotes the source, v2 p. 7). **ERROR in `paper.pdf`**: pdfTeX writes the absolute path of the included figure, `/PTEX.FileName (C:/Users/<user name>/chemistry-open-problems/…/figures/theoremO-orbit.pdf)`, into the committed PDF. C47 (`\pdfsuppressptexinfo=-1`) removes it (patched build: no occurrence). No other PDF under `problems/` contains the user name |
| 21 | Credits (`\PaperCredits`, §9.4) | OK; match the notes' headers (P, E, E2, Q1–Q3, D2, L: Claude Opus 5.5 xhigh; L2 medium; G, A2 and the archived verification: GPT-6-Astra xhigh) |
| 22 | Status "partial result"; no numerics as proof; sketches flagged | OK for the subtitle, Fig. 55 and the simulations; GAP: 17 statements have no proof environment (list in §6), two of them marked "P." in Table 1 (C44 makes the status definition cover them; proofs or pointers are R1's call) |
| 23 | Definitions before use | GAP: deficiency, linkage class, weak reversibility, endotactic, S-/T-components, and the rank-two vocabulary (siphon edge, edge flow, one-way, two-coloured, ...) are used undefined (C34, C45, C46) |
| 24 | Notation | CORRECTION: the letter D names network D (Table 10), the space span{y_c − y_c′} (§4.1) and the dual M* (Example 50) (C30, C31); further clashes listed in §7 (no change required) |
| 25 | Figures and captions (4, 46, 55, 56) | OK (Fig. 55 captioned as illustration only; arrows and labels of Fig. 56 match Thm 58) |
| 26 | Build: overfull boxes | 2 (2.7 pt at lines 620–623, 0.7 pt at 1434–1440); both removed by C22, C35 (patched build: 39 pages, 0 overfull, no undefined reference) |
| 27 | English, length of the introduction | OK; the abstract (about 390 words) is long but every sentence carries a result |

Severity counts: **ERROR** 5 (items 5, 7, 17 (rank-one shares), 18 (`\zcref`), 20 (user path inside `paper.pdf`));
**GAP** 5 (items 8, 9, 10, 22, 23); the abstract overstatement C01 is a major CORRECTION (no counterexample outside EC
was found); the rest are CORRECTIONS or OK.

## 1. Citations against the sources

| Paper (line) | Claim | Source checked | Verdict |
|---|---|---|---|
| 75–99 | Conj. 22, Rem. 23, Q. 24 verbatim, §3 p. 7, GAC paragraph p. 7–8 | `certificates/sources/v2.txt` l. 445–466 and page markers | OK |
| 74 | "in his words … chemical species and rate constants" | v2 §3 ("This duality can be interpreted as a duality between chemical species and rate constants") | OK |
| 101–102 | source's [4] = CDSS; abstract sentence | v2 l. 565, l. 9–11 | OK |
| 103–105 | HAL 2020 theorem without proof; 2024 talk lists a theorem; 2025 poster says proved | `hal.txt` l. 154–157 ("Theorem", no proof); `beamer.txt` slide 28 "Theorem (Ananiadi, I.)"; poster: abstract not on disk | OK (poster not re-checked) |
| 106 | closing diagram writes "Refine" literally and in quotes | v2 l. 515–519 | OK |
| 233, 251, 412–415 | Defs 13, 17, Rem. 18, Def. 20, Thm 21 of the source | v2 | OK |
| 240 | dual "following §3 of the source (see also BG)" | v2 §3, BG Def. 1 | OK, but earlier Γᵀ networks are not credited (item 9, C15) |
| 247, 366, 419, 1674, 1676 | archived report Lemma 1, Prop. 3, Prop. 9, Lemma 5, Thm 7, Cor. 8–10, Thm 1, §4 | `archive/paper-precise-forms.pdf` (Lemma 2 = transposition, Prop. 3, Lemma 5, Prop. 6 = MIPS of N*, Thm 7, Cors 8–10, Thm 1, §4 = exhaustive check, 15 513 networks) | ERROR for "Lemma 1" and "Prop. 9" (C16, C18, C20, C42, C43); rest OK |
| 260–263 | Hack; DE Def. 6.19 with "sometimes called reverse dual"; TS §2 incidence Cᵀ; Jančar et al. | DE p. 130 and footnote 3; TS p. 274 ("reverse-dual, Nʳᵈ, with incidence matrix Cᵀ"); Jančar Def. 3.14 | OK |
| 321, 344–345 | DE Def. 4.4; TS Thm 28 (assumes purity); DE Rank Thm 6.14 | DE, TS p. 292–293 ("Theorem 28 … EC net … well-formed iff well-structured and rank(C) = |𝔈| − 1"), TS p. 274 (purity) | OK |
| 376–377 | Jančar et al. §4 "minimal T-traps"; Hack had no interpretation of the reverse dual of a trap | Jančar §4 "Concluding remarks"; `hack.txt` l. 2084 | OK |
| 363, 378–379, cor. 34 | EC not self-dual, balanced joins self-dual | TS p. 293 and Cor. 32 say exactly this | GAP: not credited (C17, C19) |
| 409–410 | thesis Def. 4.3, Alg. 4.4; CF Fact 2 | `thesis.txt` l. 2272; `cf2012.txt` l. 608–630 | OK |
| 412 | Burgard 2004, Pfeiffer 1999, Schuster–Hilgetag 1994 (terminology) | DOIs and titles (Crossref) | OK (not re-read, as `references.md` says) |
| 529–531 | CF Fact 2(a), 2(b), Alg. 1; thesis Def. 4.6 | `cf2012.txt` l. 608–650; `thesis.txt` l. 2292 | OK |
| 600, 838 | PMDSC Condition 3.1 (kernel basis with disjoint supports), Thm 3.8, Cor. 3.12 | `literature/sources/1102.1590.txt` l. 338, 481, 608 | OK |
| 614 | thesis Thm 4.17 (isolation ⇒ monomial parametrisation) | `thesis.txt` l. 2384 | OK; wording "new" (C21) |
| 655 | DE Thm 4.30(a), same proof | DE p. 83–84 | OK |
| 707 | DE Rank Thm 6.14, Prop. 5.4, Thm 6.21 | DE | OK |
| 798–799 | ADLS Thm 2 (conservative + every siphon contains a P-semiflow support ⇒ persistent) | `Angeli-…-MathBiosci.txt` l. 526–529 | OK |
| 807–811 | DE Thms 10.16–10.17 "due to Desel"; "their proof shows that a marking marking every S-component is live" | DE p. 222–223: 10.16 (weakly connected, positive S-/T-invariants, Rank = |C|−1 ⇒ well-formed); 10.17 (live and bounded iff marks all proper siphons), proof step (iv) uses only "marks every S-component"; bibliographic notes p. 225: "Theorem 10.16 … is due to Desel [23]" | OK except that only 10.16 is credited to Desel by the book (C24) |
| 812 | DE Prop. 5.4 | DE p. 91 | OK |
| 813 | TS pp. 272, 296 on siphons | TS p. 272 l. 79 ("no satisfactory generalisation of siphons and traps"), p. 296 l. 1197–1199 | OK |
| 840 | Campos–Chiola–Silva Lemma 5.1, Thms 5.3–5.4 | checked by V2 (author's copy); not re-opened by me | OK (V2) |
| 841 | Júlvez–Recalde–Silva eq. (18), Def. 7 | checked by V2 | OK (V2) |
| 843–846 | Johnston 2014 does not cover WF-EC; A → C, 2C → 2A | V2 §8 (ker Γ = ℝ(2,1) not a cycle space) | OK |
| 877 | Müller–Regensburger "Prop. 3.9" = sign-vector criterion | arXiv:1209.6488v1: Prop. 3.1 is "(c + S) ∩ {c* ∘ e^ṽ} has at most one element for all c, c* iff the sign-vector sets meet only in 0"; Prop. 3.9 is "F is a real analytic isomorphism if σ(im V) = σ(im Ṽ) and (+,…,+) ∈ σ(im V)^⊥" | ERROR (C27) |
| 879 | spectra of ΓG and GΓ "used by Jamshidi–Palsson and BG" | BG p. 490 ("SG and GS share the same non-zero eigenvalues"); Jamshidi–Palsson 2009 abstract ("J_x = S·G … J_v = G·S", PMC copy in `literature/sources/pmc2749767.html`) | content OK; GAP: Jamshidi–Palsson not in the bibliography (C28, B09) |
| 904–906 | Maeda–Kodama–Ohta, Moreau "not opened"; Al-Radhawi–Angeli | Crossref metadata | OK (honestly flagged); DOIs now available (B05, B06) |
| 910–911 | Mairesse–Nguyen Lemma 3.2, Cor. 3.2 for ordinary nets | `0905.3158.txt` l. 1143–1146 (non-weighted free-choice), l. 1164, l. 1260 | OK |
| 944 | CFMW Thm 1 (degree) | PLoS text, Theorem 1 "Unique and multiple equilibria" | OK |
| 1098, 1106 | Birch (Horn–Jackson, Feinberg); MR for the convex argument | standard | OK |
| 1122–1123 | BG: dual of a reversible network not WR; mass action not preserved | BG p. 493 (sixth page of the PDF, only page containing it) | OK, page CORRECTION (C29) |
| 1127, 1185, 1501–1503 | CNP 2013, GMS 2014, Anderson 2011, Pantea 2012, Anderson–Shiu 2010, CDSS | titles/abstracts; CDSS Thm 9 (codimension δ) and Thm 23 (dimension two) in `0708.3431.txt`; AS2010 Thm 4.6 in `anderson-shiu-2010.txt` | OK |
| 1366–1369 | ZGD Thm 3.10 (1-d W_I-endotactic ⇒ persistent, all rates); V3's example | V3 §7, `literature2/sources/z.txt` | OK |
| 1408–1411 | Butler–Waltman Thm 4.1 as "Butler–McGehee lemma"; acyclicity theorem; Hale–Waltman | BW p. 258 Thm 4.1 ("considerably extends Lemma A1 of [4]"), Thm 3.1 = acyclicity | OK (optional: cite the acyclicity theorem as [Thm 3.1]) |
| 1496–1500 | Craciun v3 23 Sep 2026; OpenAI 25 Sep 2026 with Lean; Wiuf 21 Sep 2026 | arXiv API (v3 updated 2026-09-23; Wiuf v1 2026-09-21); `literature/sources/openai-wr-boundedness-persistence-2026-09-25.txt` | OK |
| 1519, 1546, 1547 | CDSS Thm 9; AS2010 Thm 4.6; Wiuf Thm 7.2 (V3 checked (⋄)) | as above | OK |

## 2. Identifiers (`logs/check_ids.log`, `logs/check_ids2.log`, `logs/check_ids3.log`)

All 27 DOIs of `references.bib` resolve at Crossref (the thesis DOI 10.25673/13949 at DataCite) with matching
title, journal, volume and first page; the five arXiv ids (2603.08767, 1501.02860, 2609.24553, 2602.02352, 1010.3050)
resolve with the stated titles and version dates (2603.08767: v1 9 March, v2 15 July 2026; 1501.02860 latest
23 September 2026; 2609.24553 v1 21 September 2026; 2602.02352 v2 30 July 2026). The URLs (HAL v1, the 2024 slides,
the MIT handle of Hack's report, the Campos author copy, github.com/openai/math, the thesis landing page) answer
HTTP 200. The only "FAIL" of the script is a title-markup difference for ZGD ($\mathcal W_I$), not an error.

Missing data that Crossref gives (looked up by DOI, not guessed): Campos–Chiola–Silva 10.1109/9.106153;
Recalde–Júlvez–Silva 10.3182/20020721-6-ES-1901.00078 (IFAC Proc. Vol. 35(1) 455–460); Maeda–Kodama–Ohta
10.1109/TCS.1978.1084490 (full title "…: Nonoscillation and stability", 25(6)); Moreau 10.1109/TAC.2004.841888 (50(2));
Craciun–Nazarov–Pantea 10.1137/100812355 (73(1)); Zhang–Gao–Dochain 68(12) 7420–7435; Desel–Esparza
10.1017/CBO9780511526558 (optional). Hack's first names are "Michel Henri Théodore" (report title page and
documentation page), not "Théodéor" (B01). The bibliography style prints no volume/pages, so only the DOIs show.

## 3. `references.md`

Every one of the 43 keys of `references.bib` is listed, with what was read. To change (exact replacements M01–M10 in
§8.3): Müller–Regensburger numbering and "not re-read"; BG page; DOIs of Maeda–Kodama–Ohta, Moreau (now verified),
Campos–Chiola–Silva, Recalde–Júlvez–Silva, Craciun–Nazarov–Pantea, Desel–Esparza; ZGD volume/pages; the three new
sources (read by agent L: Famili–Palsson compound-map passages, Jamshidi–Palsson abstract and main text, Duvall §11);
the list "Prior art found and credited in the paper".

## 4. Novelty and credits

Compared with `literature2/novelty.md` (L2), L's verdict (`literature/notes.md` §18–19) and V1 §4, V2 §8, V3 §7:

- Correct and present: Zhang–Gao–Dochain Thm 3.10 (Table 1 row 58/60 and §6.1, with V3's example); Desel 1994 /
  DE 10.16–10.17 and DE 5.4 (§4.2 Credits, Table 1 row 31); DE 4.30(a) (after Thm 20, row 20); PMDSC 2012 and
  Campos–Chiola–Silva 1991, Recalde–Júlvez–Silva 2002, Júlvez–Recalde–Silva 2005 (after Thm 33); BG p. 6 (after
  Prop. 48, row 48/49); Iosif's thesis Thm 4.17 (after Prop. 18, row 18); Johnston (after Thm 33, V2's example);
  Mairesse–Nguyen; Butler–Waltman/Hale–Waltman; the max/min technique (Maeda–Kodama–Ohta, Moreau, Al-Radhawi–Angeli,
  flagged as not opened); Birch for Thm 47.
- Missing (GAP): the earlier Γᵀ networks and the ΓG/GΓ comparison (Famili–Palsson 2003, Jamshidi–Palsson 2009,
  Duvall 2024; L §10, §18 item 1–2: "must be cited, not claimed"); Teruel–Silva p. 293 and Cor. 32 (L §18 item 5, L2
  §5, V1 §1.5). Fixed by C06, C15, C17, C19, C28, B09.
- Too strong for ordinary nets: "the statement for the dual [was] not found" (line 814). For ordinary nets it
  follows from DE 10.16–10.17 applied to the reverse dual, whose hypotheses are the same (DE Prop. 6.20(3)–(5) hold for
  every net; L2 §0.1 says so). Fixed by C25.
- Wording "new" (lines 615, 842, 1619): C21, C26, C39. "Mild generalisation" (656), "recognition criterion for a
  known class, not a new case of the GAC" (1538) and "we found it stated nowhere (§9.3)" (274) are fine.
- Optional (not required by any verifier): Sackmann–Heiner–Koch's maximal common transition sets are the support
  reading of MIPS(N*) (L §19); the paper uses only the ray reading, so no change is needed.

## 5. Internal consistency

### 5.1 Abstract, introduction, body

- Abstract (1): "the seeds of preclusters are the Petri clusters, and minimal siphons cross every cluster at most
  once" is stated without "equal-conflict"; Theorem 20 and the introduction state the cluster claims for EC networks
  only (DE 4.30(a) needs free choice). C01.
- Abstract (1): "with equality exactly when twin species are coupled" after "(respectively conservative)": the
  reaction-side condition (doubling pairs fully coupled, Thm 19(1)) is missing. C02.
- Abstract: "Every example and finite claim is certified … (129 claims)" contradicts §9.1's last item (the large
  counts are cited from the notes). C03 (also adds "which we have not checked" to the GAC sentence, as in §1.2 and §7).
- Introduction §1.2: "the doubling pairs … are the reaction halves of the Petri clusters when N is equal-conflict":
  V1's correction 2 (§9.2 says it was applied) requires every species consumed and every source nonzero. C04.
- The remaining claims of §1.2 match the theorems they cite (I checked every \zcref target).

### 5.2 Table 1

Rows checked against the body, Table "Status" definitions and V1–V3 verdicts: C07 (Thm 19 has V1's correction, and
the equality condition is for doubling pairs on the reaction side), C08 and C14 (Thm 24 and Prop. 73 have proofs
marked "sketch"), C11 (NC3 was found by V2, so V2 is not its verifier; the proof in the paper is the writer's and the
certificate checks it), C13 (Thm 70, a main result proved only as a sketch, has no row; "Conj. F" is undefined at that
point), C09, C10, C12 (novelty cells aligned with the body and L2). Two rows marked "P." are statements without any
proof in the paper (Props 48, 74); C44 widens the definition of the status rather than the rows.

### 5.3 Counts

| Paper (line) | Count | Note | Verdict |
|---|---|---|---|
| 490 | 15 513 networks, five families | archive §4 (l. 550, 569) | OK |
| 1124–1126 | 82 377 WR networks "of agent E's exhaustive families"; WR dual 98; δ = 0 kept 43 of 46 709; endotactic 122 of 84 103 "tested sides"; strongly 43 of 70 185 | `computation/notes.md` l. 27–29, 231–253: the 82 377 include the two seeded samples (n = 5, 6); 84 103 = 34 927 + 35 258 + 13 918 runs, "about 1 700 networks counted twice" | CORRECTION (C29) |
| 1135–1136 | 5 252 / 924 | E l. 239 | OK |
| 1196–1202 | 141 (s ≤ 2), 103 (s = 3), 238 (s = 4), 1 735 drain sets, 1 432 without drain species, "the s = 2 count 199" | E l. 243, 257–265; E2 l. 156; Q1 table l. 254–259 (1 432 = 148 + 791 + 493); 199 is the s = 2 count over all of E's records (140 from the 82 377 + 59 from (4,3,4), (5,3,4)), not a recount of 141 | OK; wording C32 |
| 1361–1364 | eight families, up to 3 544, "ALL for 17–35 %"; 4 754 rank-one networks | Q1 table (`round2/persistence/notes.md` l. 192–201): shares 4/12, 115/460, 381/3 544 = 10.75 %, 662/2 742, 4/12, 101/426, 245/829, 1/1; E2 l. 114 (4 754) | ERROR: 11–33 % (C33; `scripts/rank1_percent.py`) |
| 1434–1438 | 18 861; 2 988; 1 048; 14 654 of 14 825; 4 873; 13 817; 171 | E2 l. 233–234; D2 l. 18–19; V3 table | OK (2 988 + 1 885 = 4 873; 1 048 + 12 769 = 13 817; sums to 18 861) |
| 1007 | 4 × 120 random PC networks | V2 §2.10 | OK |
| 1574 | 140 simulated networks | Q3 (100) + V2 (40) | OK |
| 1584 | 63 717 of 63 722 sides | Q2 l. 28, 366 | OK |
| 716 | eight relabellings for M3/M4 | V2 §2.11 | OK |
| 1018–1020, 1030–1032, 1050–1051, 1148–1149 | NC3, NC, NM, NB, M (unit rates, T = 2) | V2 §2.12, CHANGES A1–A5 (arithmetic is R1's) | OK as attributions |

### 5.4 Section 9, open problems, concordance

- §9.2: "zcref{thm:B}" lacks the backslash and prints "zcrefthm:B" in the PDF (page 32, `build/paper.txt` l. 1972):
  C40. "the citation of Desel–Esparza Thm 6.17" is listed as applied, but the paper does not state P's D(h) and cites
  no Thm 6.17; "the caveats on K8" are moot because the paper uses Q3's M3/M4: C41. "the species-side restriction of
  E's summary sentence" refers to a sentence that the paper does not reproduce; harmless.
- §9.1: "It certifies every example and finite claim of this paper" overclaims (C38); whether every example is in fact
  covered is R1's item.
- Open problem 6 is titled with a label of the notes ("Theorem C"): C37. The other eight open problems match the body
  (171 undecided, NC3/NC absolutely consistent, 63 717/63 722, Conjecture F).
- Concordance: every note label exists in the cited note (`grep` on P, Q1–Q3, D2, A2, G, E, E2, L); verifier sections
  match V1–V3; archive numbers wrong (C42, C43). Optional: add "V1 §1.8–1.9" to the row of Thm 31/Cor. 32 (V1 also
  confirmed P's Theorems C, C′).

## 6. Honesty and house rules

- "independent" (non-mathematical) at lines 66, 163, 1363 (twice), 1438, 1620: C03, C05, C33, C36, C39.
- No mention of the owner, no personal name other than cited authors, no absolute path in `paper.tex` (but one inside
  `paper.pdf`, item 20 and C47), no e-mail content (the three
  questions Q1–Q3 are framed from the source; the acknowledgment sentence is exactly the prescribed one).
- Credits: `\PaperCredits` and §9.4 agree with the notes' headers and with the brief. After both reviews are applied,
  the coordinator may want to add the reviewers (Claude Opus 5.5, effort xhigh) to "Verification".
- Status "partial result" in the subtitle. Numerical integration is presented as illustration only (Fig. 55, §5.3,
  open problem 2). Counts are labelled "checked by computation" with their notes.
- Statements without a proof environment (`scripts/list_unproved.py`, `logs/list_unproved.log`): lemma 2 and Prop. 8
  (archived report), Prop. 23, Prop. 25 (certificate), Lemma 28 (classical), Cor. 32, Cor. 34, Cor. 43, Prop. 48,
  Prop. 49 (certificate), Lemma 57, Cor. 59, Cor. 61, Lemma 64, Prop. 66 (said to be "without proofs"), Prop. 71,
  Prop. 74. GAP for R1: each needs a one-line proof or "Proof: [note, label], re-derived by Vk". C44 makes Table 1's
  status definition honest in the meantime.

## 7. Presentation and readability

- Build of the unchanged copy: 38 pages, 2 overfull boxes (lines 620–623: 2.67 pt; 1434–1440: 0.74 pt), no undefined
  reference, PDF text identical to the committed `paper.pdf`. Patched copy (all replacements): 39 pages, 0 overfull, no
  undefined reference or citation (`logs/build_patched.log`).
- Figures: Fig. 4 (reverse dual; matrices and weights match Example 50), Fig. 46 (classes; NC3 outside PC, Table 45
  networks outside the outer box, as the caption says), Fig. 55 (captioned as a numerical illustration, not part of the
  proof; panel titles match the caption: period about 2.5·10⁶), Fig. 56 (arrows match Thm 58(1)–(3)). OK.
- Definitions before use, for a Petri-net referee: deficiency, linkage class, weak reversibility, complex balance,
  (strong) endotacticity are used from §1 on and never defined (C45). For a reaction-network referee: S- and
  T-components (Prop. 12's proof, Thm 24) are undefined (C46). §6.3 uses "siphon facet", "siphon edge", "edge flow",
  "one-way", "two-way", "two-coloured", "transverse order μ", "tropically Hurwitz", "attracting sector" without
  definition (C34 defines the first group from D2's notes and points to D2 §2–3 for the rest).
- Notation clashes: D is network D (Table 10), span{y_c − y_c′} (§4.1) and M* (Example 50, fixed by C30, C31);
  S is the stoichiometric subspace and network S (Example 62); X is network X and the species X_i; H is network H
  (Prop. 49), the matrix H (Thm 54) and the functions H_c; U, W, T, R, F, Q, B, C, P name both networks and other
  objects. No change required beyond C30–C31, but the applier may rename networks S, X, H if time permits.
- Theorem 54 (§5.3) uses Theorem 68 (§6.3) in its statement and proof: a forward reference, not circular (Thm 68 does
  not use Thm 54). Acceptable.
- The shared counter makes "Figure 46", "Table 76": house template, no change.
- Table 1 starts at the foot of p. 3 with one row; acceptable for a longtable.
- English: good. Introduction (pp. 1–5 including Table 1) is of reasonable length.

## 8. Exact replacements

### 8.1 Severity of each replacement

| Severity | Ids |
|---|---|
| ERROR (wrong number, wrong count, broken output, house rule) | C16, C18, C20, C42, C43 (archive numbering); C27 (Müller–Regensburger); C33 (rank-one shares; also "independent"); C40 (`\zcref`); C47 (user path in the PDF; rebuild `paper.pdf` after applying) |
| GAP (attribution or definitions missing) | C06, C15, C28, B09 (Famili–Palsson, Jamshidi–Palsson, Duvall); C17, C19 (Teruel–Silva); C34, C45, C46 (definitions); C44 (status of unproved statements) |
| CORRECTION, major (claims stronger than the body or than the sources) | C01, C02, C03, C04 (abstract/intro); C25 (ordinary dual known); C38, C41 (§9) |
| CORRECTION, minor | C05, C36, C39 ("independent", "new"); C21, C26 ("new"); C07–C14 (Table 1); C22, C35 (overfull boxes); C23, C24, C29, C30, C31, C32, C37; B01–B08 (bibliography); M01–M10 (`references.md`) |

### 8.2 `paper.tex` and `references.bib`

Generated from `scripts/apply_corrections.py` (`--markdown`); the script applies exactly these replacements to a copy
and checks that each old text occurs once. Apply them in any order (no two overlap).

#### C01 (paper.tex, line 47)

Old:
```latex
theory (Hack 1972). Through it, siphons of \(\Net\) are the T-traps of \(\Nd\), the seeds of preclusters are the
Petri clusters, and minimal siphons cross every cluster at most once, so preclusters and siphons are related by
transversality, not refinement.
```
New:
```latex
theory (Hack 1972). Through it, siphons of \(\Net\) are the T-traps of \(\Nd\), and for equal-conflict networks the
seeds of preclusters are the Petri clusters and minimal siphons cross every cluster at most once, so preclusters and
siphons are related by transversality, not refinement.
```

#### C02 (paper.tex, line 51)

Old:
```latex
(respectively conservative) network, with equality exactly when twin species are coupled, in particular on the
```
New:
```latex
(respectively conservative) network, with equality exactly when doubling pairs (respectively twin species) are
coupled, in particular on the
```

#### C03 (paper.tex, line 65)

Old:
```latex
preprints. Every example and finite claim is certified in exact arithmetic (129 claims); the main results were
re-derived by three independent verifiers.
```
New:
```latex
preprints, which we have not checked. Every example is certified in exact arithmetic (129 claims); counts over large
families are cited from the working notes as checked by computation. The main results were re-derived by three
fresh-context verifiers.
```

#### C04 (paper.tex, line 123)

Old:
```latex
Petri clusters when \(\Net\) is equal-conflict.
```
New:
```latex
Petri clusters when \(\Net\) is equal-conflict, every species is consumed and no source is zero (e.g.\ under (H1)
and (H2)).
```

#### C05 (paper.tex, line 163)

Old:
```latex
\emph{Verification}: the independent verifier that re-derived
```
New:
```latex
\emph{Verification}: the fresh-context verifier that re-derived
```

#### C06 (paper.tex, line 182)

Old:
```latex
construction known (Hack; Desel--Esparza; Teruel--Silva); identification not found \\
```
New:
```latex
construction known (Hack; Desel--Esparza; Teruel--Silva); \(\Gamma\tp\) networks known (Famili--Palsson; Blanchini--Giordano; Duvall); identification not found \\
```

#### C07 (paper.tex, line 188)

Old:
```latex
equality iff twins coupled & P. & V1 & not found \\
```
New:
```latex
equality iff doubling (twin) pairs coupled & P. & V1 (corr.) & not found \\
```

#### C08 (paper.tex, line 191)

Old:
```latex
\ref{thm:WFFC} & structure of well-formed free-choice networks & P. & V1 &
```
New:
```latex
\ref{thm:WFFC} & structure of well-formed free-choice networks & P. (sketch) & V1 &
```

#### C09 (paper.tex, line 194)

Old:
```latex
extends Pérez Millán et al.\ 2012, Campos et al.\ 1991 \\
```
New:
```latex
extends Pérez Millán et al.\ 2012, Campos et al.\ 1991; cf.\ Recalde et al.\ 2002, Júlvez et al.\ 2005 \\
```

#### C10 (paper.tex, line 195)

Old:
```latex
Hurwitz at every point & P. & V2 & not found \\
```
New:
```latex
Hurwitz at every point & P. & V2 & not found (tools classical) \\
```

#### C11 (paper.tex, line 200)

Old:
```latex
does not give uniqueness (NC3) & P. & V2 & not found \\
```
New:
```latex
does not give uniqueness (NC3) & P. & found by V2; certificate & not found \\
```

#### C12 (paper.tex, line 205)

Old:
```latex
sufficient half weaker in Zhang--Gao--Dochain 2023 \\
```
New:
```latex
``if'' half of (1) under a stronger hypothesis: Zhang--Gao--Dochain 2023; rest not found \\
```

#### C13 (paper.tex, line 208)

Old:
```latex
\ref{thm:A} & a persistent dual with an attracting-looking facet (refutes Conj.~F) & P. & V3 & not found \\
```
New:
```latex
\ref{thm:A}, \ref{thm:AF} & a persistent dual with an OUT facet (refutes a conjecture of agents Q1 and D2); the family \(N_{u,R}\) & P.; P. (sketch) & V3 & not found \\
```

#### C14 (paper.tex, line 209)

Old:
```latex
duality keeps \(s\); prime-cycle certificate & P. & V3 &
```
New:
```latex
duality keeps \(s\); prime-cycle certificate & P.; P. (sketch) & V3 &
```

#### C15 (paper.tex, line 244)

Old:
```latex
\(\Net\) gives a loop \(y\to y\) of \(\Nd\), and a species never consumed gives the zero complex; all definitions
apply to such networks.
```
New:
```latex
\(\Net\) gives a loop \(y\to y\) of \(\Nd\), and a species never consumed gives the zero complex; all definitions
apply to such networks. Networks with stoichiometric matrix \(\pm\Gamma\tp\) appeared before: Famili and Palsson's
compound map of \(-\Gamma\tp\), on which the extreme conservation laws are drawn as extreme pathways
\autocite{FamiliPalsson2003}; the dual network of Blanchini and Giordano~\autocite{BlanchiniGiordano2022}, developed at
the same time as Iosif's (source, p.~8); and Duvall's dual of an irreversible network~\autocite{Duvall2024}. Jamshidi
and Palsson compare the Jacobians \(\Gamma G\) and \(G\Gamma\)~\autocite{JamshidiPalsson2009}. These works transfer
linear and local properties (conservation pools, spectra, structural stability, monotonicity), not siphons or
persistence.
```

#### C16 (paper.tex, line 247)

Old:
```latex
\begin{lemma}[transposition {\autocite[Lemma~1]{ArchivedReport}}]
```
New:
```latex
\begin{lemma}[transposition {\autocite[Lemma~2]{ArchivedReport}}]
```

#### C17 (paper.tex, line 363)

Old:
```latex
\(\Ye[S_c,T_c]=a_cJ\) is symmetric \\
```
New:
```latex
\(\Ye[S_c,T_c]=a_cJ\) is symmetric; Teruel--Silva Cor.~32 \\
```

#### C18 (paper.tex, line 366)

Old:
```latex
MIPS (\zcref{prop:mips}) & classes of fully coupled reactions & \autocite[Prop.~9]{ArchivedReport} \\
```
New:
```latex
MIPS (\zcref{prop:mips}) & classes of fully coupled reactions & \autocite[Prop.~6]{ArchivedReport} \\
```

#### C19 (paper.tex, line 379)

Old:
```latex
\(K_1\to K_2\), \(K_2\to K_1\) is not.
```
New:
```latex
\(K_1\to K_2\), \(K_2\to K_1\) is not. Teruel and Silva~\autocite[p.~293]{TeruelSilva1996} make the same
observation for weighted nets: equal-conflict nets with balanced joins form a self-reverse-dual class, and their Duality
Theorem~\autocite[Cor.~32]{TeruelSilva1996} states that a pure net is a well-formed balanced-joins EC net iff its
reverse dual is.
```

#### C20 (paper.tex, line 419)

Old:
```latex
\begin{proposition}[{\autocite[Prop.~3, Prop.~9]{ArchivedReport}}]\label{prop:mips}
```
New:
```latex
\begin{proposition}[{\autocite[Prop.~3, Prop.~6]{ArchivedReport}}]\label{prop:mips}
```

#### C21 (paper.tex, line 614)

Old:
```latex
the rate-free form with \(\Ph\) is
new as a statement.
```
New:
```latex
the rate-free form with \(\Ph\) was
not found in the literature as a statement (\zcref{sec:search}).
```

#### C22 (paper.tex, line 620)

Old:
```latex
\item For every consistent network, \(\Mips(\Nd)\refi\Ph(\Net)\); and \(\Ph(\Net)=\Mips(\Nd)\) iff \(\Ph(\Net)\refi
\Mips(\Nd)\) iff both have the same number of blocks iff every doubling pair is fully coupled.
```
New:
```latex
\item If \(\Net\) is consistent, then \(\Mips(\Nd)\refi\Ph(\Net)\); and \(\Ph(\Net)=\Mips(\Nd)\) iff
\(\Ph(\Net)\refi\Mips(\Nd)\) iff both have the same number of blocks iff every doubling pair is fully coupled.
```

#### C23 (paper.tex, line 725)

Old:
```latex
E.3 of agent E shows that a Galois connection
```
New:
```latex
Agent E's Proposition~E.3 shows that a Galois connection
```

#### C24 (paper.tex, line 808)

Old:
```latex
Esparza~\autocite[Thms~10.16--10.17]{DeselEsparza1995}, due to Desel~\autocite{Desel1994}, show that an ordinary net
```
New:
```latex
Esparza~\autocite[Thms~10.16--10.17]{DeselEsparza1995} (Theorem~10.16 is due to Desel~\autocite{Desel1994}) show that an ordinary net
```

#### C25 (paper.tex, line 813)

Old:
```latex
The
weighted statement, its persistence reading, \zcref{prop:K5gen} outside free choice and the statement for the dual
were not found (\zcref{sec:search}).
```
New:
```latex
For
ordinary nets the statement for the reverse dual follows in the same way, since the hypotheses of Theorem~10.16 are
invariant under it (Desel--Esparza Prop.~6.20). The weighted statement (for \(\Net\) and for \(\Nd\)), its persistence
reading and \zcref{prop:K5gen} outside free choice were not found (\zcref{sec:search}).
```

#### C26 (paper.tex, line 842)

Old:
```latex
None of these has a uniqueness or stability result under mass action. New here: weighted EC under mass action, all
real rates, the converse (2), and the dual statements below.
```
New:
```latex
None of these has a uniqueness or stability result under mass action. Not found in the literature
(\zcref{sec:search}): weighted EC under mass action, all real rates, the converse (2), and the dual statements below.
```

#### C27 (paper.tex, line 877)

Old:
```latex
criterion of Müller and Regensburger~\autocite[Prop.~3.9]{MullerRegensburger2012}
```
New:
```latex
criterion of Müller and Regensburger~\autocite[Prop.~3.1]{MullerRegensburger2012}
```

#### C28 (paper.tex, line 879)

Old:
```latex
\(G\Gamma\) used by Jamshidi--Palsson and Blanchini--Giordano~\autocite{BlanchiniGiordano2022}.
```
New:
```latex
\(G\Gamma\) used by Jamshidi and Palsson~\autocite{JamshidiPalsson2009} and by Blanchini and
Giordano~\autocite[Thm~1]{BlanchiniGiordano2022}.
```

#### C29 (paper.tex, line 1123)

Old:
```latex
mass action is not preserved~\autocite[p.~6]{BlanchiniGiordano2022}. On the 82\,377 weakly reversible networks with
(H1)--(H3) of agent E's exhaustive families, the dual is weakly reversible for 98, deficiency zero is preserved for 43
of 46\,709, the dual is endotactic for 122 of 84\,103 tested sides and strongly endotactic for 43 of 70\,185 (checked by
computation, \texttt{gac-duality/computation/notes.md} §5).
```
New:
```latex
mass action is not preserved~\autocite[p.~493]{BlanchiniGiordano2022}. On the 82\,377 weakly reversible networks with
(H1)--(H3) of agent E's exhaustive families and two seeded samples, the dual is weakly reversible for 98, deficiency
zero is preserved for 43 of 46\,709, the dual is endotactic for 122 of 84\,103 tests and strongly endotactic for 43 of
70\,185 (counted per run, about 1\,700 networks twice; checked by computation, \texttt{gac-duality/computation/notes.md}
§5).
```

#### C30 (paper.tex, line 1145)

Old:
```latex
\(D=M^*\colon 3U+2V\to2U+3V,\ V\to U\)
```
New:
```latex
\(M^*\colon 3U+2V\to2U+3V,\ V\to U\) (\(U=K_1\), \(V=K_2\) in \zcref{fig:revdual})
```

#### C31 (paper.tex, line 1146)

Old:
```latex
\(D\) is non-persistent for every rate
```
New:
```latex
\(M^*\) is non-persistent for every rate
```

#### C32 (paper.tex, line 1202)

Old:
```latex
\texttt{gac-duality/round2/persistence/notes.md} §3; the \(s=2\) drain-species count 199 was recounted by V3).
```
New:
```latex
\texttt{gac-duality/round2/persistence/notes.md} §3); over all of E's records, 199 duals with \(s=2\) have a drain
species, a count recounted by V3.
```

#### C33 (paper.tex, line 1362)

Old:
```latex
(up to 3\,544 networks per family), \(\Nd\) is persistent for every rate vector for 17--35\,\% and for none for the rest,
never SOME (Q1; recounted independently by V3 with identical numbers). E2's independent decision procedure agrees on
```
New:
```latex
(up to 3\,544 networks per family), \(\Nd\) is persistent for every rate vector for 11--33\,\% of the networks of each
family with more than one member (381 of 3\,544 to 4 of 12) and for none for the rest, never SOME (Q1; recounted by V3
with its own generator, identical numbers). E2's separate decision procedure agrees on
```

#### C34 (paper.tex, line 1390)

Old:
```latex
some class (agent D2, Lemma~2.1).
```
New:
```latex
some class (agent D2, Lemma~2.1). On an open edge whose zero set \(Z\) is a siphon of \(\Nd\) (a \emph{siphon edge};
its facet set is a \emph{siphon facet}) only the reactions of \(\Nd\) whose source avoids \(Z\) act, and they move
along the edge; their sum is the \emph{edge flow}, \emph{one-way} if all its terms have one sign, \emph{two-way} if both
signs occur, \emph{zero} if there is no such reaction. The finer notions used in \zcref{prop:none-cert}
(two-coloured edge flows, the transverse order \(\mu\), tropically Hurwitz vertices, centre manifolds, attracting
linear functions and sectors) are defined in agent D2's notes (\texttt{gac-duality/round3/dichotomy/notes.md} §2--3).
```

#### C35 (paper.tex, line 1434)

Old:
```latex
\paragraph{Counts (checked by computation).} Of the 18\,861 weakly reversible networks with (H1)--(H3) and \(s=2\) in
```
New:
```latex
\paragraph{Counts (checked by computation).} Among the 18\,861 weakly reversible networks with (H1)--(H3) and
\(s=2\) in
```

#### C36 (paper.tex, line 1438)

Old:
```latex
V3, whose independent sampled test of the ALL verdicts
```
New:
```latex
V3, whose own sampled test of the ALL verdicts
```

#### C37 (paper.tex, line 1581)

Old:
```latex
\item \emph{Weighted impure Theorem C in the Petri sense.} The siphon statement is settled by \zcref{thm:Cgen}; a liveness
theory for weighted impure EC nets (Teruel--Silva assume purity) is not.
```
New:
```latex
\item \emph{Liveness of weighted impure equal-conflict nets.} The siphon statement is settled by \zcref{thm:Cgen} and
\zcref{cor:C}; a liveness theory for weighted impure EC nets (Teruel--Silva assume purity) is not.
```

#### C38 (paper.tex, line 1598)

Old:
```latex
It certifies every example and finite claim of this paper:
```
New:
```latex
It certifies every example of this paper and every finite claim about an explicit network:
```

#### C39 (paper.tex, line 1619)

Old:
```latex
The new results were re-derived by three verifiers, Claude Opus 5.5 (effort xhigh), each with
a fresh context and an independent exact library
```
New:
```latex
The results of the working notes were re-derived by three verifiers, Claude Opus 5.5 (effort xhigh), each with
a fresh context and its own exact library
```

#### C40 (paper.tex, line 1622)

Old:
```latex
transversality, zcref{thm:B}, the WF-FC structure
```
New:
```latex
transversality, \zcref{thm:B}, the WF-FC structure
```

#### C41 (paper.tex, line 1626)

Old:
```latex
restriction of E's summary sentence; the citation of Desel--Esparza Thm~6.17; the caveats on K8.
```
New:
```latex
restriction of E's summary sentence; the caveats on K8 (Q3's pairs \(M_3,M_4\), used here, satisfy (H1)--(H3) on
both sides). V1's citation correction (Desel--Esparza Thm~6.17) concerns P's liveness statement D(h), not included here.
```

#### C42 (paper.tex, line 1674)

Old:
```latex
& archive Lemma~1; P 1.1, 1.2 &
```
New:
```latex
& archive Lemma~2; P 1.1, 1.2 &
```

#### C43 (paper.tex, line 1676)

Old:
```latex
& archive Prop.~3, 9, Thm~7, Cor.~8--10, Thm~1 &
```
New:
```latex
& archive Prop.~3, 6, Thm~7, Cor.~8--10, Thm~1 &
```

#### C44 (paper.tex, line 161)

Old:
```latex
\emph{Status}: proved here (with a complete proof or a proof sketch
marked as such), checked
```
New:
```latex
\emph{Status}: proved (here, with a complete proof or a proof sketch
marked as such; for the statements given without proof, in the working notes of \zcref{tab:concordance}, re-derived by
the verifier), checked
```

#### C47 (paper.tex, line 4)

Old:
```latex
\documentclass[11pt,a4paper]{report}

```
New:
```latex
\documentclass[11pt,a4paper]{report}
\pdfsuppressptexinfo=-1 % keep the absolute path of included figures (PTEX.FileName) out of the PDF

```

#### C45 (paper.tex, line 224)

Old:
```latex
Compatibility classes are \((x^0+S)\cap\Rnn^n\). Conjecture~22 assumes:
```
New:
```latex
Compatibility classes are \((x^0+S)\cap\Rnn^n\). The complexes and reactions form a directed graph whose
connected components are the \emph{linkage classes}; \(\Net\) is \emph{weakly reversible} if every reaction lies on a
directed cycle of this graph, and its \emph{deficiency} is the number of complexes minus the number of linkage classes
minus \(s\)~\autocite{Feinberg2019}; complex balance and (strong) endotacticity are used as
in~\autocite{HornJackson1972,CraciunNazarovPantea2013,GopalkrishnanMillerShiu2014}. Conjecture~22 assumes:
```

#### C46 (paper.tex, line 327)

Old:
```latex
of all nonzero entries of \(\Ye\) and \(\Yp\) is connected.
```
New:
```latex
of all nonzero entries of \(\Ye\) and \(\Yp\) is connected. S-components (strongly connected subnets generated by a
set of places with all their input and output transitions, in which every transition has one input and one output
place) and T-components (places and transitions exchanged) are as in Desel and
Esparza~\autocite[Defs~5.1, 5.11]{DeselEsparza1995}.
```

#### B01 (references.bib, line 218)

Old:
```bibtex
Hack, Michel Henri Th{\'e}od{\'e}or
```
New:
```bibtex
Hack, Michel Henri Th{\'e}odore
```

#### B02 (references.bib, line 282)

Old:
```bibtex
  volume  = {36},
  number  = {12},
  pages   = {1368--1382},
  year    = {1991},
```
New:
```bibtex
  volume  = {36},
  number  = {12},
  pages   = {1368--1382},
  year    = {1991},
  doi     = {10.1109/9.106153},
```

#### B03 (references.bib, line 292)

Old:
```bibtex
  booktitle = {Proceedings of the 15th IFAC World Congress},
  year      = {2002}
```
New:
```bibtex
  booktitle = {Proceedings of the 15th IFAC World Congress},
  series    = {IFAC Proceedings Volumes},
  volume    = {35},
  number    = {1},
  pages     = {455--460},
  year      = {2002},
  doi       = {10.3182/20020721-6-ES-1901.00078}
```

#### B04 (references.bib, line 329)

Old:
```bibtex
  title   = {Asymptotic behavior of nonlinear compartmental systems},
```
New:
```bibtex
  title   = {Asymptotic behavior of nonlinear compartmental systems: nonoscillation and stability},
```

#### B05 (references.bib, line 332)

Old:
```bibtex
  pages   = {372--378},
  year    = {1978},
  note    = {Cited after the second literature search; not opened}
```
New:
```bibtex
  number  = {6},
  pages   = {372--378},
  year    = {1978},
  doi     = {10.1109/TCS.1978.1084490},
  note    = {Cited after the second literature search; not opened}
```

#### B06 (references.bib, line 342)

Old:
```bibtex
  pages   = {169--182},
  year    = {2005},
  note    = {Cited after the second literature search; not opened}
```
New:
```bibtex
  number  = {2},
  pages   = {169--182},
  year    = {2005},
  doi     = {10.1109/TAC.2004.841888},
  note    = {Cited after the second literature search; not opened}
```

#### B07 (references.bib, line 387)

Old:
```bibtex
  pages      = {305--329},
  year       = {2013},
```
New:
```bibtex
  number     = {1},
  pages      = {305--329},
  year       = {2013},
  doi        = {10.1137/100812355},
```

#### B08 (references.bib, line 426)

Old:
```bibtex
  journal = {IEEE Transactions on Automatic Control},
  year    = {2023},
  doi     = {10.1109/TAC.2023.3275880}
```
New:
```bibtex
  journal = {IEEE Transactions on Automatic Control},
  volume  = {68},
  number  = {12},
  pages   = {7420--7435},
  year    = {2023},
  doi     = {10.1109/TAC.2023.3275880}
```

#### B09 (references.bib, append at the end)

```bibtex

@article{FamiliPalsson2003,
  author  = {Famili, Iman and Palsson, Bernhard {\O}.},
  title   = {The convex basis of the left null space of the stoichiometric matrix leads to the definition of
             metabolically meaningful pools},
  journal = {Biophysical Journal},
  volume  = {85},
  number  = {1},
  pages   = {16--26},
  year    = {2003},
  doi     = {10.1016/S0006-3495(03)74450-6}
}

@article{JamshidiPalsson2009,
  author  = {Jamshidi, Neema and Palsson, Bernhard {\O}.},
  title   = {Flux-concentration duality in dynamic nonequilibrium biological networks},
  journal = {Biophysical Journal},
  volume  = {97},
  number  = {5},
  pages   = {L11--L13},
  year    = {2009},
  doi     = {10.1016/j.bpj.2009.06.049}
}

@misc{Duvall2024,
  author     = {Duvall, Alon},
  title      = {Computationally checking if a reaction network is monotone or non-expansive},
  year       = {2024},
  eprint     = {2409.12487},
  eprinttype = {arXiv},
  note       = {Version 2, 27 January 2025}
}
```


### 8.3 `references.md`

Exact replacements (old text → new text); each old text occurs once in the current file.

#### M01 (section "Reaction-network background", Müller–Regensburger)

Old:
```
- S. Müller, G. Regensburger, SIAM J. Appl. Math. 72 (2012) 1926–1947, DOI 10.1137/110847056 (sign-vector criterion,
  Proposition 3.9). Cited after agents Q3, L2; title and pages matched through Crossref on 28 October 2026; not re-read.
```
New:
```
- S. Müller, G. Regensburger, SIAM J. Appl. Math. 72 (2012) 1926–1947, DOI 10.1137/110847056, arXiv:1209.6488
  (sign-vector criterion: Proposition 3.1 of arXiv v1, read by reviewer R2 on 28 October 2026; the numbering of the
  journal version was not checked). Cited after agents Q3, L2; title and pages matched through Crossref.
```

#### M02 (Blanchini–Giordano)

Old: `Remark 2, p. 6 (the dual of a reversible network is not`
New: `Remark 2, p. 493 (the sixth page; the dual of a reversible network is not`

#### M03 (Maeda–Kodama–Ohta, Moreau)

Old: `of Theorem 37; their DOIs were not verified and are not given.`
New:
```
of Theorem 37. DOIs 10.1109/TCS.1978.1084490 (full title: *Asymptotic behavior of nonlinear compartmental systems:
  nonoscillation and stability*) and 10.1109/TAC.2004.841888, matched through Crossref by reviewer R2 on 28 October 2026.
```

#### M04 (Zhang–Gao–Dochain)

Old: `Control (2023), DOI 10.1109/TAC.2023.3275880.`
New: `Control 68(12) (2023) 7420–7435, DOI 10.1109/TAC.2023.3275880.`

#### M05 (Campos–Chiola–Silva)

Old: `IEEE Trans. Automat. Control 36 (1991) 1368–1382. Author's copy`
New: `IEEE Trans. Automat. Control 36 (1991) 1368–1382, DOI 10.1109/9.106153 (Crossref). Author's copy`

#### M06 (Recalde–Júlvez–Silva)

Old: `Congress (2002). Abstract and §1–3 read (L2).`
New: `Congress (2002), IFAC Proc. Vol. 35(1) 455–460, DOI 10.3182/20020721-6-ES-1901.00078 (Crossref). Abstract and §1–3 read (L2).`

#### M07 (Craciun–Nazarov–Pantea)

Old: `J. Appl. Math. 73 (2013) 305–329, arXiv:1010.3050.`
New: `J. Appl. Math. 73 (2013) 305–329, DOI 10.1137/100812355 (Crossref), arXiv:1010.3050.`

#### M08 (Desel–Esparza)

Old: `0-521-46519-2 (author-hosted PDF). Read: Def. 4.4,`
New: `0-521-46519-2, DOI 10.1017/CBO9780511526558 (Crossref) (author-hosted PDF). Read: Def. 4.4, Defs 5.1, 5.11,`

#### M09 (new sources; insert after the Blanchini–Giordano item, i.e. after the line ending `weakly reversible; mass action is not preserved).`)

Insert:
```
- I. Famili, B. Ø. Palsson, *The convex basis of the left null space of the stoichiometric matrix leads to the
  definition of metabolically meaningful pools*, Biophys. J. 85 (2003) 16–26, DOI 10.1016/S0006-3495(03)74450-6
  (open access, PMC1303061). Compound-map passages read by agent L (`gac-duality/literature/notes.md` §10(a)): the
  compound map of −Sᵀ, with the extreme conservation laws drawn as extreme pathways.
- N. Jamshidi, B. Ø. Palsson, *Flux-concentration duality in dynamic nonequilibrium biological networks*, Biophys. J.
  97 (2009) L11–L13, DOI 10.1016/j.bpj.2009.06.049 (open access, PMC2749767). Abstract and main text read by agent L
  (§10(b)): the Jacobians S·G and G·S.
- A. Duvall, *Computationally checking if a reaction network is monotone or non-expansive*, arXiv:2409.12487 v2 (27
  January 2025). §11 ("Duality") read by agent L (§10(e)): the dual of an irreversible network has stoichiometric
  matrix Γᵀ; monotonicity and non-expansivity transfer.
```
(DOIs and the arXiv id checked by R2, `gac-duality/review/R2/logs/check_ids3.log`.)

#### M10 (section "Novelty checks", 27 October 2026, last bullet)

Old: `Google Scholar was not reachable. Prior art found and credited in the paper: Desel–Esparza Thm 4.30(a), Prop. 5.4,`
New: `Google Scholar was not reachable. Prior art found and credited in the paper: Famili–Palsson 2003, Jamshidi–Palsson 2009 and Duvall 2024 (networks with stoichiometric matrix ±Γᵀ); Teruel–Silva p. 293 and Cor. 32; Desel–Esparza Thm 4.30(a), Prop. 5.4,`

## 9. Points for R1 (mathematics) found while reviewing

These are outside my scope; I record them so that R1 can decide.

1. **Doubling pairs and doubling classes (lines 321–324).** "Two reactions with the same source form a doubling pair"
   includes the zero source, but "doubling classes" are defined with "equal nonzero source", while "the twin classes of
   N are the doubling classes of N*" includes zero rows of Ye (V1's correction to Theorem 21), which correspond to
   reactions of N* with the zero source. Under (H1) for N* (i.e. (H2) for N) no zero source occurs, so nothing in the
   theorems breaks, but the definitions should agree (e.g. "classes of reactions with equal source, the zero source
   included"). I did not write a replacement because it touches Theorem 21's hypotheses.
2. **Line 341–342** ("Under (H2) every species that occurs is consumed …, so for EC networks clusters correspond to
   sources and q = m"): a species that occurs in no reaction is a cluster of its own, so q = m needs every species to
   occur; harmless under the WF-EC hypotheses where it is used.
3. **Abstract (C01).** Outside EC I found no counterexample to "minimal siphons cross every Petri cluster at most once"
   in the time available, but neither the paper nor DE 4.30 (free choice) proves it; C01 restricts the abstract to
   what Theorem 20 states.
4. **Statements without proof** (§6 above): 17 environments, list in `logs/list_unproved.log`.
5. **Müller–Regensburger (C27).** With Prop. 3.1 the equivalence in the sentence after Theorem 35 reads: steady
   states of a WF-EC network are x*∘e^v with v ∈ D^⊥, so "at most one per class for all c, c*" is σ(S) ∩ σ(D^⊥) = {0};
   please confirm that this is the intended use (it is how I read it).

## 10. What I could not check

- The numbering of Müller–Regensburger's journal version (paywalled); I used arXiv:1209.6488v1, its only arXiv
  version.
- The 2025 poster abstract (Iosif–Ananiadi); not on disk.
- Campos–Chiola–Silva and Júlvez–Recalde–Silva theorem numbers: I relied on V2, who read them; I did not reopen them.
- Johnston 2014: relied on V2 (read from a copy outside the repository).
- Burgard 2004, Pfeiffer 1999, Schuster–Hilgetag 1994, Horn–Jackson 1972, Feinberg 2019, Maeda–Kodama–Ohta 1978,
  Moreau 2005, Al-Radhawi–Angeli 2016: metadata only (cited for terminology or flagged as not opened).
- The large counts themselves (E, E2, Q1, Q2, D2): I compared them with their notes and the verifiers' recounts; I did
  not re-run the enumerations.
- The mathematics of the replacements C04, C25 and C34 rests on the notes (V1 §1.3, DE Prop. 6.20, D2 §2); R1 should
  glance at them.
- No literature search of my own beyond the identifier checks: the novelty verdicts are L's, L2's and the verifiers'.

## 11. Scripts and logs (this folder)

| file | content | runtime |
|---|---|---|
| `scripts/check_ids.py` → `logs/check_ids.log` | every DOI (Crossref, DataCite fallback) and arXiv id of `references.bib`; bibliographic query for entries without DOI | about 1 min (network) |
| `scripts/check_ids2.py` → `logs/check_ids2.log` | direct lookup of the DOIs found for entries without DOI; resolution of the URLs and the MIT handle | about 30 s |
| `scripts/check_ids3.py` → `logs/check_ids3.log` | Famili–Palsson 2003, Jamshidi–Palsson 2009, Duvall arXiv:2409.12487 | about 10 s |
| `scripts/rank1_percent.py` → `logs/rank1_percent.log` | exact shares of ALL verdicts in Q1's eight rank-one families (exits 1: the paper's "17 %" lower bound fails) | < 1 s |
| `scripts/list_unproved.py` → `logs/list_unproved.log` | theorem-like statements of `paper.tex` without a proof environment | < 1 s |
| `scripts/apply_corrections.py` → `logs/apply_corrections.log`, `logs/replacements.md`, `build/paper.tex`, `build/references.bib` | applies C01–C47 and B01–B09 to a copy; `--markdown` prints §8.2 | < 1 s |
| `scripts/refresh_review.py` | rewrites §8.2 of this file from the script | < 1 s |
| `logs/build_patched.log`, `build/paper.pdf` | build of the patched copy (39 pages, no overfull box, no undefined reference) | about 1 min |
| `build/paper.txt`, `build/paper_patched.txt` | `pdftotext` of the unpatched and patched builds (the unpatched one equals the committed `paper.pdf`) | – |

Environment: Windows, Python 3.14 (`py -3`, standard library), MiKTeX (latexmk, pdflatex, biber), poppler
(`pdftotext`, `pdftoppm`). Logs contain no absolute path.
