# Independent verification: rigorous upper bounds for the SAW connective constant on Z^3 … Z^6

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

Date: 18 October 2026. Inputs read: `notes.md`, `paper.tex`, `README.md`, `references.md`, `scripts/`,
`certificates/`, `../../PROTOCOL.md`. `src/verify_cert.cpp` was opened only **after** my own checker had been
written and run on every certificate. `PROGRESS.md` and the solver's working logs were not read.

## Verdict: CONFIRMED WITH CORRECTIONS

All four claimed bounds are proved by the stored certificates. I re-checked them with a checker of my own, written
from the mathematical statement:

| claim | certificate | my checker | exact maximal ratio (reduced) |
|---|---|---|---|
| μ(Z³) ≤ 4.722048 | d3_k18.cert2, 206 118 560 states | **PASS**, 151 s, 3.3 GB | 73151633491513/15491508107156 = 4.722047265218… |
| μ(Z³) ≤ 4.726984 | d3_k16.cert, 11 641 894 states | **PASS**, 384 s, 0.6 GB | 129653733242276/27428428361269 = 4.726983680383… |
| μ(Z⁴) ≤ 6.801078 | d4_k14.cert, 4 235 309 states | **PASS**, 91 s, 0.3 GB | 433123876427845/63684594277489 = 6.801077738528… |
| μ(Z⁵) ≤ 8.850674 | d5_k14.cert, 7 787 256 states | **PASS**, 193 s, 0.3 GB | 622810764598795/70368744177664 = 8.850673290777… |
| μ(Z⁶) ≤ 10.885008 | d6_k12.cert, 192 224 states | **PASS**, 2.3 s | 4428155100988748/406812324202497 = 10.885007256526… |

Times are wall-clock with 4 threads in WSL. Each exact ratio equals, as a fraction, the one recorded by the solver;
for example 146303266983026/30983016214312 for K = 18 is my value multiplied by 2/2. The proof (Lemmas 1–2)
is correct. Each previous record is the best one I could find. The corrections are editorial and change no bound
(§5).

## 1. Proof review (re-derived)

I re-derived Lemma 1 (domination c_n ≤ a_n) and Lemma 2 (one-sided Collatz–Wielandt bound) line by line. Both
are correct as stated in `notes.md` §3. The specific points I was asked to examine:

* **Start state.** The start state s₀ needs F(s₀) = ∅, because invariant (2) at i = 0 requires h₀F(s₀) ⊆ ∅.
  It is the state at time 0, so "reachability" is not an issue. The checkers verify that the empty walk is listed.
  The empty walk is also always a legal child, since ∅ ⊆ W, which is harmless.
* **V > 0 on all states, or only on reachable ones?** Positivity on the set R of states reachable from s₀ is
  enough. R is closed under σ, so (A_R V)_s = (AV)_s for s ∈ R, and Lemma 2 applies to A_R with min over R. The
  certificates have V ≥ 1 on every listed state, which is stronger. Lemma 2 really needs min V > 0: with
  V_t = 0 the bound a_n ≤ λⁿV_{s₀}/min V is void. Both checkers reject a zero weight (tested).
* **Composition of frames.** In (1), g maps the child's frame to the parent's frame, relative to the new
  position. Hence h_{i+1} = h_i g(s_i,e_i) is the right composition. The label e_i = h_i⁻¹f_i is read in the
  parent's frame. The inclusion h_{i+1}F(s_{i+1}) ⊆ h_i((F(s_i)∪{0}) − e_i) = (h_iF(s_i)∪{0}) − f_i uses only
  linearity of h_i and the fact that h_i maps E to E. Injectivity holds because s_i, h_i and hence f_i = h_ie_i
  are functions of e₀…e_{i−1}.
* **Is (1) sufficient?**
  * *Old position.* The site 0 is in the union, so the child may remember it, and it then sits at −e.
  * *New position.* The new position can never be remembered: 0 ∈ gF(t) would need e ∈ F(s) ∪ {0}, which is
    excluded.
  * *Forgotten sites.* Forgetting a site only enlarges the set of accepted label sequences.
  * *Steps onto remembered sites.* These are exactly the excluded labels e ∈ F(s); a self-avoiding walk never
    takes them, by (2).

  Lemma 1 holds for arbitrary finite sets F(s), whether or not they come from self-avoiding walks.
* **From μ = lim c_n^{1/n} to the bound.** c_n ≤ a_n ≤ (V_{s₀}/min V)λⁿ, so c_n^{1/n} ≤ C^{1/n}λ → λ. Only
  limsup is needed; existence of the limit (Hammersley, submultiplicativity) is classical.
* **What soundness rests on.** Only (1) for each listed state and allowed step, F(s₀) = ∅, V ≥ 1 and
  Q·Σ_e V_{σ(s,e)} ≤ P·V_s. Generation (choice of states, candidate children, value or policy iteration,
  floating point) is irrelevant. I did not read or run the generators.

## 2. Independent checker (`verification/sawcheck.cpp`)

I wrote it from the definition in `notes.md` §3 and the format in §5. It shares no code with `src/`. In
`notes.md`, the format description left one thing implicit: where step t_i sits in the code. I inferred bits
[b·i, b·i+b), with b = 3 (d ≤ 4) or 4 (d = 5, 6), from the data: the codes are sorted, every walk starts with
t₀ = 0, and all codes decode to self-avoiding walks. This is now documented (§5). For each certificate, the
checker verifies:

* **(F) Format.**
  * magic, file size = 40 + 16N (+ 8dN);
  * codes strictly increasing and decodable: bit 63 set, l·b ≤ 57, no stray bits above the last step, step
    values < 2d;
  * 1 ≤ V ≤ 2⁶⁰ and 1 ≤ P, Q < 2⁶⁰.
* **(S) Start state.** The empty walk is listed.
* **(I) Condition (1)** for every listed s and every e ∉ F(s), with W = {p_i − e : 0 ≤ i ≤ l}. The required g is
  re-found by an exhaustive backtracking search over partial signed permutations. An axis is assigned when the
  child's walk first uses it, and every site is tested against W. This is a complete search of 𝒢 restricted to
  the axes used.
  * *SAWCERT2:* t is the stored child.
  * *SAWCERT1:* t is searched among **all** self-avoiding paths through W, not only the longest prefixes the
    solver keeps. The only pruning is k + |P_k|₁ ≤ max over the listed states of (l + |p_l|₁), which is safe
    because k + |p_k|₁ never decreases along a walk. Candidates are looked up through my own canonical form
    (every listed code was checked to be canonical in that form). The child of minimal weight is kept, then
    re-tested with the backtracking inclusion test.
* **(C) Inequality.** Q·Σ_e V_t ≤ P·V_s in unsigned 128-bit arithmetic. The bounds V ≤ 2⁶⁰ and P, Q < 2⁶⁰ keep
  every product below 2¹²⁴. The exact maximal ratio is reported as a reduced fraction.

Results:
* All 13 distinct certificates pass: d = 3 with K = 8, 10, 12, 12 (explicit), 14, 16, 18 (explicit); d = 4,
  K = 14; d = 5, K = 10, 12, 14; d = 6, K = 10, 12. Logs: `verification/run_large_v1.log` and
  `verification/run_all.log`.
* No listed walk is self-intersecting.
* The SAWCERT2 files carry no pointer on a forbidden direction.
* The maximal ratios of my min-weight choice are identical to the solver's on every SAWCERT1 file, including
  K = 16. My complete candidate search therefore finds nothing better than the solver's search per branch.
* For K = 8 and 10 my ratios reproduce the values quoted from the solver's Python brute-force checker
  (4.771089138921, 4.753600356).
* The large files are streamed: the K = 18 child table (4.9 GB) is read in chunks, with peak RSS 3.3 GB.

SHA-256 values were recomputed before use. All 14 entries of `certificates/SHA256SUMS.txt` match: 9
decompressed small files and the 5 large files in `~/work/cubic-saw/` (d3_k16, d4_k14, d5_k14_p4, d3_k18.cert,
d3_k18.cert2). Beyond the magic, `d3_k18.cert` and `d3_k18.cert2` have identical headers, codes and weights
(`cmp` on the first 3 297 896 992 bytes after the magic). The SAWCERT1 file therefore carries the same weights as
the checked explicit automaton.

`verification/peek.py` is the small script used to read the header and the first codes when inferring the format.
The first runs, made with the first version of `sawcheck.cpp`, are in `verification/run_large_v1.log`. The final
version only adds the optional a_n count and does not change the checking logic. It re-runs everything through
`verification/run_all.sh`, with output in `verification/run_all.log`.

## 3. Falsification tests (`verification/tamper.py`, logs `verification/tamper*.log`)

Each corrupted copy gets a **ground truth computed independently of both checkers**. On SAWCERT2 it is an
exact Python re-check: the untouched pairs are known to be valid, and the touched pairs are tested against (1) by
brute force over all 2^d d! signed permutations. Sample of 3187 states of the untouched K = 12 file: 0
violations of (1). Results on d3_k12 (SAWCERT1), d3_k12_explicit (SAWCERT2), d5_k10 and d6_k10:

* **Every invalid copy is rejected by both my checker and `src/verify_cert.cpp`.** The invalid copies were:
  * P − 1, Q + 1, P = 0;
  * one-step weight +10 % and start weight −10 % (SAWCERT1), zero weight, weight above 2⁶⁰;
  * corrupted magic, file truncated by one byte, two codes swapped, d + 1;
  * one-step state removed, start state removed;
  * on SAWCERT2:
    * a child redirected to a walk violating (1), set to none, or set out of range;
    * a child redirected to the start state, which satisfies (1) but breaks the inequality;
    * the children of two directions swapped, chosen so that (1) fails;
    * all children of the one-step state removed.
* Copies that are still **valid** are accepted by both checkers:
  * on SAWCERT2, one-step weight +10 % and start weight −10 %, because the start state has slack;
  * a pointer placed on a forbidden direction;
  * d = 6 → 5 on the K = 10 file. The states use at most 5 axes and the encoding is the same, so this is a
    correct, weak certificate of μ₅ ≤ 10.887669.

  My first version of the script had labelled these "must fail". Ground truth showed they are valid.
* **Stray payload bit.** A bit set above the last step of a code is rejected by my checker (strict format) and
  accepted by `verify_cert`. This is soundness-neutral: the bit does not change the decoded walk or F, so the
  automaton is unchanged.
* **K = 18 itself** (copy of the 8.2 GB file, patched in place):
  * P − 1 is rejected by both checkers: 206 054 038 rows violate the inequality. My checker took 227 s on 3
    threads, `verify_cert` 435 s on 3 threads with 8.1 GB.
  * A child of the tightest state (180715238) was redirected to the last listed walk, which violates (1)
    (checked by brute force over all 48 signed permutations, `verification/patch_k18.py`). It is also rejected
    by both checkers, each reporting exactly 1 violation of (1). Times: 201 s for mine, 414 s for
    `verify_cert` (`verification/tamper_k18_child.log`). The copy was deleted afterwards.

## 4. Cross-checks

* **c_n by brute force** (`verification/count_saw.cpp`, plain depth-first search): d = 3 up to n = 16 (92 s),
  d = 4 up to 12, d = 5 up to 10, d = 6 up to 9. All values equal the Clisby–Liang–Slade tables, which I
  downloaded myself from personal.math.ubc.ca/~slade/lacecounts/. The copies in `certificates/counts/` are
  identical to the public files.
* **a_n vs c_n** for the automata my checker actually checked (`verification/counts.log`): d = 3 with K = 8 to
  14 (including the explicit K = 12 automaton) up to n = 30; d = 5 and d = 6 with K = 10, 12 up to n = 24.
  * a_n ≥ c_n for every n, as Lemma 1 predicts.
  * a_n = c_n exactly for n ≤ 9, 10, 12, 12, 13 (d = 3, K = 8 to 14), n ≤ 10 and n ≤ 11 (d = 5, 6), which
    confirms "equality for n ≤ K − 1 at least".
  * At n = 30 and K = 14 the overcount ratio is 1.050.

## 5. Reproduction of the Pönitz–Tittmann table (`verification/pt_auto.cpp`)

This is my own implementation of the automaton A(d,k) from Pönitz–Tittmann §2:
* states are the walks that "could be the first part of a loop of length ≤ k", tested by breadth-first search
  for a shortest return path avoiding v₁…v_{l−1};
* the oldest steps are removed until the walk is admissible;
* states are reduced by symmetry, and the reachable ones are generated by breadth-first search from the empty walk;
* the growth rate is bracketed in floating point by Collatz–Wielandt.

The Table 2 values were read from the P–T PDF itself (Electron. J. Combin. 7, R21, p. 9). Every entry is
reproduced after rounding up to 4 decimals:

| d | k | my value | P–T Table 2 |
|---|---|---|---|
| 2 | 4, 6, 8, 10, 12, 14, 16, 18, 20, 22 | 2.8311772, 2.7755911, 2.7444582, 2.7247990, 2.7112523, 2.7013743, 2.6938489, 2.6879244, 2.6831389, 2.6791926 | 2.8312, 2.7756, 2.7445, 2.7248, 2.7113, 2.7014, 2.6939, 2.6880, 2.6832, 2.6792 |
| 3 | 4, 6, 8, 10, 12, **14** | 4.8645365, 4.8074109, 4.7779699, 4.7598376, 4.7475251, **4.7386308** | 4.8646, 4.8075, 4.7780, 4.7599, 4.7476, **4.7387** |
| 4 | 4, 6, 8, 10, 12, **14** | 6.8916709, 6.8512155, 6.8302905, 6.8178559, 6.8096952, **6.8039597** | 6.8917, 6.8513, 6.8303, 6.8179, 6.8097, **6.8040** |
| 5 | 4, 6, 8, **10** | 8.9104204, 8.8815540, 8.8678881, **8.8601588** | 8.9105, 8.8816, 8.8679, **8.8602** |
| 6 | 10 | 10.8885689 | (Noonan 1998: 10.8886, via CLS 2007) |

* The reachable state counts equal the solver's state counts: 2723, 41 425, 678 730 and 11 641 894 (d = 3, K =
  10 to 16); 139 267 and 4 235 309 (d = 4, K = 12, 14); 5732 (d = 5, 6, K = 10); 8 880 776 (d = 2, K = 22).
* All values agree with the solver's P–T-rule values to every printed digit. The d = 4, k = 14 entry (6.8040),
  which the solver had not recomputed, is now reproduced.
* For d = 2, k = 22 I obtain 2.6791926494. This confirms that P–T's abstract value 2.679192495 is too small,
  while their Table 2 entry 2.6792 is right, as reported by the solver and by vitplanocka (2026).
* A first run of my program gave a wrong d = 2, k = 22 value because of an encoding overflow on my side (3 bits
  per step for d = 2). It was fixed (2 bits) and the value above comes from the corrected program.

These are floating-point brackets for comparison only, not certificates.

## 6. Novelty (searched 18 Oct 2026)

* **Pönitz–Tittmann 2000** (full PDF read). Table 2: d = 3, k = 14: **4.7387**; d = 4, k = 14: **6.8040**; d = 5,
  k = 10: **8.8602**; no d = 6 column. "The values shown are true upper bounds." My reproduction confirms that
  4.7387 and 6.8040 are valid roundings of the automaton constants.
* **Clisby–Liang–Slade 2007 §7.2** (preprint se.pdf read). "These are not as good as the upper bounds 4.7387,
  6.8179, 8.8602, 10.8886 for dimensions d = 3, 4, 5, 6 [53, 56]. It is erroneously claimed in [47] that c26
  gives µ ≤ 4.7114 for d = 3; in fact c26 gives the weaker bound µ ≤ 4.7626 using (77)." [53] = Noonan 1998,
  [56] = P–T. The 4.7114 of MacDonald et al. (2000) is therefore not a valid prior bound. For d = 4, CLS quote
  P–T's k = 10 value; P–T's own k = 14 value 6.8040 is stronger and is the correct comparison, as the solver uses.
* **Slade 1995** (Kahane volume, read), Table 1: the best upper bounds at the time are Alm's 4.756, 6.832, 8.881,
  10.903 (d = 3..6). The lace expansion (Hara–Slade 1992) supplies lower bounds and 1/d expansions, not better
  numerical upper bounds. CLS 2007 (Slade a co-author) likewise list no lace-expansion upper bound. I did not
  read Hara–Slade 1992 itself.
* **Couronné 2022** (arXiv:2211.16146v2, read), d = 2 only. Final remark: in d ≥ 3 only the deletion of old
  vertices and small bridges carry over. No d ≥ 3 computation.
* **Kim–Pinna 2025** (arXiv:2509.04568v1, only version). It quotes "4.57213 ≤ µSAW ≤ 4.73871 for the cubic
  lattice (Hara et al., 1993, Pönitz and Tittmann, 2000)" as the current state of the art. Quote checked
  verbatim.
* **Alm 1993 / 2005.** Not read in full. 1993: values via Slade 1995 Table 1 (all weaker). 2005: title and
  abstract, Archimedean and Laves lattices (2D).
* **Citations of P–T.** OpenAlex: 43 citing works (titles of all 43 scanned, 4 since 2024). Semantic Scholar: 45
  citing works. None improves an upper bound for d ≥ 3. Bauerschmidt–Duminil-Copin–Goodman–Slade 2012 cite P–T
  only for d = 2.
* **arXiv API, 2024-01-01 to 2026-12-31.** Queries: `abs:"connective constant"` (17 hits), `"self-avoiding" AND
  "upper bound"` (5), `"self-avoiding walk" AND cubic` (8), `"self-avoiding walks" AND hypercubic` (1). Nothing
  on Z^d upper bounds; Kim–Pinna (surfaces) and He 2025 (weighted, Alm method) are not comparable. Abstracts of
  2601.01128, 2508.01993 and 2609.21724 read.
* **GitHub.**
  * `ipitchford` (Evidence Press): 166 repositories, no SAW or connective-constant repository. `mmislan`: 23
    repositories, only binary-polymer RAF work.
  * `gh search repos` ("connective constant", "self-avoiding walk upper bound", "self-avoiding cubic", "Ponitz
    Tittmann"): only `vitplanocka/saw-connective-constant` (Z² lower bounds, plus the P–T d = 2 erratum),
    He's weighted code, Union-Jack Alm code and Monte Carlo repositories.
  * `gh search code`: `ewang26/HorizonMath` (numerical μ₃ estimate only), `openai/math` (honeycomb preprints;
    749 preprint folders, none on Z^d SAW upper bounds), `neelsomani/gpt-erdos` problem 528 (estimates only),
    `elder-plinius/V3R1T4S` (lists only the 2D bound as open), `teorth/optimizationproblems` (only entry 38a,
    Z²).
* **Zenodo** (`"connective constant"`, `"connective constants"`, `"self-avoiding" AND "upper bound"`,
  `"self-avoiding walk" AND cubic`): only the Z² lower-bound archive 10.5281/zenodo.21546041 and a hard-square
  record.

Conclusion: no rigorous upper bound below 4.7387 (d = 3), 6.8040 (d = 4), 8.8602 (d = 5) or 10.8886 (d = 6)
was found. The new values are records by margins of 0.0167, 0.0029, 0.0095 and 0.0036.

## 7. What I only read (not recomputed)

* The generators `src/saw_prod.cpp`, `src/saw_big.cpp`, `src/saw_automaton.cpp`, `src/count_states.cpp` and the
  regeneration scripts `scripts/reproduce_*.sh` were neither read nor run. Soundness does not depend on them. The
  certificates were checked as stored, with SHA-256 matching.
* `src/verify_cert.cpp` was read after my runs. Its (C0)–(C2) are correct:
  * complete enumeration of 𝒢 for d ≤ 3 and complete backtracking for d ≥ 4;
  * 128-bit products without overflow;
  * the empty walk is always a candidate;
  * pointers on forbidden directions are ignored.

  It ignores stray payload bits and does not check the file size (trailing bytes are ignored). Both points are
  soundness-neutral.
* `scripts/verify_small_bruteforce.py`, `scripts/check_counts.py` and `scripts/check_small.sh` were read but not
  run; my own tools cover the same checks.
* Noonan 1998 (the d = 6 previous record 10.8886) was not read; the attribution rests on CLS 2007 [53, 56] and on
  P–T having no d = 6 column. Hara–Slade 1992 and Alm 1993 were not read (see §6).

## 8. Corrections made to the deliverables

1. `certificates/SHA256SUMS.txt`: size of d3_k16.cert corrected from "186 271 344" to **186 270 344** bytes
   (= 40 + 16 × 11 641 894; the hash was right).
2. `notes.md` §5 and `paper.tex` §5, format: added the bit position of step t_i (bits [b·i, b·i+b), b = 3 for
   d ≤ 4, 4 for d = 5, 6, zero above the last step), the 40-byte header, and the layout of the SAWCERT2 child
   table: entry 2d·s + j is the child of state s for direction j in the state's own frame (2a ↦ +e_a,
   2a+1 ↦ −e_a); entries for forbidden directions are 2³²−1 and are ignored.
3. `notes.md` §5, `paper.tex` §5, `README.md`: the P–T reproduction now records the verifier's independent
   values, including d = 4, k = 14 (6.8039597 ↦ 6.8040) and all d = 2 entries. The README sentence "reproduces
   every value of the Pönitz–Tittmann table" was only partly supported by the solver's runs, which lacked d = 4,
   k = 14 and d = 2, k ≤ 8, 18, 20. It is now supported.
4. `notes.md`, `paper.tex`, `README.md`: status changed from "awaiting independent verification" to verified;
   verifier attribution added (README and `\PaperCredits`); pointer to this file. §7's remark on Hara–Slade 1992
   was completed with what the verifier read.

No mathematical statement, bound or proof step needed correction.

## 9. Remaining reservations

* The proof is computer-assisted. Trust rests on two independently written checkers (mine and the solver's),
  agreeing to the exact fraction, and on the compiler and hardware. There is no formal (Lean) proof of the
  lemmas or of the checker.
* My checker and the solver's interpret the format the same way. This cannot create a false PASS: any consistent
  decoding of the file defines some forbidden-set automaton, and a PASS proves the bound for that automaton.
* The novelty search covered journals indexed by OpenAlex and Semantic Scholar, arXiv, GitHub and Zenodo. It cannot
  exclude unindexed or unpublished work. Noonan 1998, Hara–Slade 1992 and Alm 1993 were assessed through CLS 2007
  and Slade 1995, not read directly.
* The margins over the previous records are small for d = 4 (0.0029) and d = 6 (0.0036). They are exact
  consequences of the certificates, but the comparison is with values printed to 4 decimals.

## 10. Reproduction

Run `bash verification/run_all.sh` under WSL from this directory (about 35 min, at most 4 threads, peak 3.4 GB). The
script:
* rebuilds the verifier's programs and the solver's `verify_cert` into `$WORK` (default `~/verif_saw`);
* checks every SHA-256;
* runs the independent checker on all 13 certificates (large ones read in place from `~/work/cubic-saw`, read only);
* compares c_n (brute force) and a_n with the CLS tables;
* reproduces the P–T table;
* runs the falsification tests on both checkers.

Its output is `verification/run_all.log`. The two K = 18 corruption tests need an 8.2 GB scratch copy of
`d3_k18.cert2` and were run by hand:

```
python3 verification/patch_k18.py copy.cert2 P_minus_1                # then: sawcheck copy.cert2 3; verify_cert copy.cert2 3
python3 verification/patch_k18.py copy.cert2 P_restore 4722048
python3 verification/patch_k18.py copy.cert2 child 180715238 206118559
```
