# Literature and novelty review — the Halvorsen–Skogestad conjecture (lower bound and finite-stage attainment)

Reviewer: Claude Opus 5.5 (effort xhigh), manuscript and literature review. Date: 2026-10-20. This review extends the
targeted search of the stage-level argument (`total/literature/novelty.md`, same day, stopped early) and the novelty
checks recorded in `references.md`; it covers both halves of the result and states what the manuscript must cite.

## 0. Verdict

**(i) Lower bound — new.** No source proves that every arrangement of a broad class (recycles, splitters, several
feed stages, reverse flows, finite stages, impure products, non-adiabatic sections) needs at least max_j H_j, nor
uses the argument of the paper (the Underwood flux of a stream, conserved at every unit in Underwood's weights at the
feed roots θ_j; sign uniformity of the outflows of an equilibrium stage; summation over the units of one sign). Its
ingredients are classical and must be cited as such:
- Underwood's stage identity for a free parameter (Underwood 1948; Halvorsen 2001, eqs. 2.45–2.49);
- the sign of Underwood coordinates inside one section (Franklin and Forsyth 1953);
- stage-level sign comparisons of an Underwood-type function inside one section of a multi-feed multi-product
  column (Jiang 2019, Prop. 5.4.1; Jiang, Tawarmalani and Agrawal 2022, Prop. 1) — the closest technique;
- cut/vapour-balance arguments at the common roots for one arrangement or one junction (HS 2003 Part 3 §3, equality
  case for the Petlyuk arrangement; Halvorsen 2001 §6.2.1; T. J. Mathew 2021 thesis §6.3, one four-component family).
For a single two-product column the bound reduces to Underwood's classical minimum-vapour bound at the feed roots,
which is classical in substance; the network statement is not.

**(ii) Finite-stage attainment for N ≥ 3 — new.** No source proves that finite-stage steady states of the
generalized Petlyuk arrangement (N ≥ 3) reach pure products at a total vaporization tending to max_j H_j. What exists:
- HS's infinite-stage construction (HS 2003 Parts 2–3; root carry-over, Carlberg and Westerberg 1989);
- finite-stage simulations that confirm it numerically: Halvorsen 2001 §5.2.7 (four components, 30 stages per
  section, flows taken from the V_min diagram); Halvorsen and Skogestad 1999 (steady-state behaviour of a ternary
  Petlyuk column); Malinen and Tanskanen 2007 (rigorous minimum-energy computation for the fully thermally coupled
  system, simulation with a large number of stages);
- uniqueness and multiplicity of steady states of single columns (binary uniqueness going back to Rosenbrock, see
  Doherty and Perkins 1982; multiplicity in ideal two-product columns, Jacobsen and Skogestad 1991). No existence
  theorem for closed thermally coupled networks at prescribed flows, and no convergence theorem to the Underwood
  configuration, was found.
The N = 2 case (a single column) is classical in substance (McCabe–Thiele pinch argument); the paper gives a
self-contained proof.

**Closest prior art overall.** HS 2003 Part 3 itself (the conjecture, the formula max_j H_j, proved there for the
extreme peaks, and the preferred-split construction), Halvorsen's thesis (Underwood identities, junction lemma,
numerical confirmation), and, for technique, Jiang 2019 / Jiang et al. 2022 Prop. 1. For the Underwood-model
optimization literature (Nallasivam et al. 2016; Tumbalam Gooty et al. 2024): these compute optima numerically in
Underwood's infinite-stage model and do not prove the FTC optimal; the permissive form of that model admits a
five-component point below the FTC value (certified in this project), which no physical family can approach.

**No AI-generated or other competing candidate** was found (GitHub ipitchford: 176 repositories, none on
distillation; mmislan: 25, none; `gh search repos` on Petlyuk, "Halvorsen Skogestad", "minimum vapor distillation",
"Underwood distillation", "thermally coupled distillation": nothing relevant (one 2014 JavaScript repository named
Petlyuk, unrelated); Zenodo and arXiv: nothing relevant).

## 1. Searches done in this review (2026-10-20, 19:20–19:45)

| Source | Query | Outcome |
|---|---|---|
| OpenAlex | works citing HS Parts 1–3 (W2084600958, W2037012710, W2086353211): 114 + 143 + 119 citations, **235 unique works**, 114 with abstracts | keyword scan (prove/proof/conjecture/lower bound/lowest/any configuration/all configurations/finite stages/rigorous); sentences read; see §2 |
| OpenAlex | abstracts by DOI: Jacobsen–Skogestad 1991, Doherty–Perkins 1982, Agrawal–Tumbalam Gooty 2020, Malinen–Tanskanen 2007, Navarro-Amorós et al. 2013, Jiang et al. 2022, Jiang et al. 2018 | see §2–3 |
| Crossref | bibliographic queries for every reference added to the manuscript (DOIs verified) | §4 |
| Semantic Scholar | DOI lookups | rate-limited (HTTP 429); one record obtained (Jiang et al. 2018, open-access copy at OSTI) |
| Web search | Malinen–Tanskanen; "rigorously proved" + HS; existence/uniqueness of multicomponent column steady states; finite-stage convergence of Petlyuk columns to the Underwood minimum; 2025–2026 proofs of FTC optimality; Agrawal–Tumbalam Gooty 2020 | §2–3 |
| Halvorsen thesis (on disk) | §2.4.3, §3.4.7–3.4.8, §3.5, §5.2.7 | finite stages treated by examples and simulation; no convergence proof |
| GitHub, Zenodo, arXiv | as in §0 | nothing relevant |

Earlier searches of the same day (`total/literature/novelty.md`, `references.md`): OpenAlex citations of HS 1–3, FK
1986/1987, Agrawal–Fidkowski 1998, FA 2001 and the Halvorsen thesis (594 works in total), Semantic Scholar, arXiv (11
queries), Zenodo (15), GitHub, about 20 web searches in English, Russian and Chinese, open theses (NTNU, Purdue, Ulm),
and full-text reading of the sources in `certificates/sources/`. Their conclusions are confirmed here.

## 2. Findings in the citing literature

**Statements of the conjecture as open (primary sources).** HS 2002 (Distillation & Absorption, paper 6-27) and HS
2003 Part 3 (abstract; §5.3 p. 624; §6.3 p. 626; §6 p. 624 citing Petlyuk, private communication 2000, "no proof has
been found in the literature, except for the ternary case"); Fidkowski–Agrawal 2001 (p. 2723, "not yet developed a
rigorous proof"); Fidkowski 2006 (p. 2103–2104); T. J. Mathew, Purdue thesis 2021, p. 125 ("For four and more
components, the conjecture is yet to be proved or disproved") and p. 134.

**Statements of FTC optimality as known or proved (secondary; to be cited as misattributions).**
- Halvorsen–Skogestad 2011, J. Nat. Gas Sci. Eng. 3, 571–580, p. 576 ("It is also shown that ..."), refers back to
  the 2003 papers, which call it a conjecture.
- Caballero–Grossmann 2014, Comput. Chem. Eng. 61, 118–135 (accepted manuscript p. 4: "Halvorsen & Skogestad (2003c)
  proved ...").
- **New in this review:** Navarro-Amorós, Ruiz-Femenia and Caballero 2013, AIChE J. 59, 3767–3781,
  doi:10.1002/aic.14137, abstract: "Even though it has been proved that a fully thermally coupled distillation (TCD)
  system minimizes the energy used by a sequence of columns, ...".
- Nallasivam et al. 2016, AIChE J. 62, p. 2079 ("a mathematical lower bound"); Shah–Agrawal 2010 p. 1760; Z. Jiang,
  Purdue thesis 2019, p. 27 ("later demonstrated", i.e. by numerical enumeration); Kadu, Nogaja, Tawarmalani and
  Agrawal 2026 (abstract: "While the FTC design minimizes vapor duty"); A. S. Nogaja, Purdue thesis 2025.
- A web search snippet attributes "rigorously proved" to HS on an Elsevier topic page; the page could not be opened
  (HTTP 403) and the source was not identified: not cited.

**Ties and heat integration (consistent with the paper).** Jiang, Ramapriya, Tawarmalani and Agrawal 2018 (AIChE J.
64, 3410–3418): heat-and-mass-integrated "HMP" configurations reach "very similar and sometimes even the same minimum
total vapor duty" as the FTC; the lower bound shows they cannot go below it without heat exchange between process
streams. Mathew et al. 2019 (IECR 58, 21671–21678): configurations below the FTC heat duty exist only with heat
transfer between process streams, which HS exclude.

**Finite stages.** Malinen and Tanskanen 2007 (Chem. Eng. Res. Des. 85, 502–509, doi:10.1205/cherd06086; conference
paper D&A 2006): rigorous minimum-energy calculation for the FTC system by simulating columns with a large number
of stages (bounded Newton homotopy) — numerical. Halvorsen–Skogestad 1999 (J. Process Control 9, 407–424): finite-stage
steady-state behaviour of a ternary Petlyuk column — numerical. Halvorsen 2001 §5.2.7: a four-component Petlyuk
simulation with 30 stages per section at the V_min-diagram flows, "a practical confirmation" — numerical; §3.4.7: with
finitely many stages the real minimum vapour "has to be slightly higher". Nothing proves convergence.

**Existence and uniqueness of finite-stage steady states.** For binary ideal columns with molar flow specifications
uniqueness is classical (Rosenbrock, extended to non-ideal binaries by Doherty and Perkins 1982, Chem. Eng. Sci. 37,
381–392); Jacobsen and Skogestad 1991 (AIChE J. 37, 499–511) show multiple steady states in ideal two-product columns
with mass or volume inputs, or with an energy balance. A general multicomponent uniqueness result was not found. The
paper claims only existence and local uniqueness near the limit configuration, and notes (as the project's computations
show) that the exact equations have non-physical zeros far from it.

**Spectral/finite-stage tools.** The projective linear form of the constant-α, CMO section map and its eigenvectors
(pinch compositions) at the Underwood roots are classical in substance (Underwood 1948 as presented in Halvorsen 2001
§2.5; Franklin–Forsyth 1953; the exact ideal-cascade analyses of Acrivos and Amundson 1955). The M-matrix structure of
the component balances at fixed K-values is the standard structure of stage-to-stage component balances; the paper
proves what it uses (weak chained diagonal dominance, Shivakumar and Chew 1974).

## 3. Agrawal–Tawarmalani group (global optimization of configurations)

Nallasivam et al. 2013 (AIChE J. 59, 971–981) and 2016 (AIChE J. 62, 2071–2086), Tumbalam Gooty et al. 2019 (Comput.
Chem. Eng. 125, 13–30), 2023 (Comput. Chem. Eng. 178, 108369) and 2024 (Oper. Res. 72, 639–659), Jiang et al. 2018,
2019 (total cost, exergy), Jiang et al. 2022 (MFMP minimum reflux), Mathew et al. 2019, 2021, 2023 (CMO relaxation):
global optimization in Underwood's infinite-stage model, numerical; the FTC value is used as a reference or a feasible
point, never proved minimal; no network-level stage bound. Jiang 2019 Prop. 5.4.1 / Jiang et al. 2022 Prop. 1 is the
closest technique (stage-level sign comparisons of f_n(x) = Σ v_{i,n}/(α_i − x) inside one section, at the section's own
stage roots, to prove convergence to a pinch); it has no feed root, no network, no energy balance, no bound. Mathew's
2021 thesis §6.3 proves, by a vapour balance on a boundary, that one symmetric-split four-component family cannot beat
the FTC; it states the general conjecture as open. Agrawal and Tumbalam Gooty 2020 (AIChE J. 66, e16294) is a
perspective; its abstract is not available in OpenAlex and the text was not read (no claim taken from it).

## 4. What the manuscript cites (all DOIs checked on Crossref or the publisher)

Primary: HS 2003 Parts 1–3; HS 2002 (D&A paper 6-27); Halvorsen 2001 thesis. Configurations and Underwood theory:
Petlyuk, Platonov and Slavinskii 1965; Petlyuk 2004 (book); Kaibel 1987; Underwood 1948; Franklin–Forsyth 1953;
Fidkowski–Królikowski 1987; Carlberg–Westerberg 1989; Christiansen, Skogestad and Lien 1997; Agrawal 2003; Fidkowski–
Agrawal 2001; Fidkowski 2006; Shah–Agrawal 2010. Optimization: Nallasivam et al. 2016; Tumbalam Gooty et al. 2024;
Jiang et al. 2018; Jiang et al. 2022; Jiang 2019 thesis; Mathew 2021 thesis; Mathew et al. 2019; Mathew, Tawarmalani
and Agrawal 2023. Misattributions: HS 2011; Caballero–Grossmann 2014; Navarro-Amorós et al. 2013. Finite stages and
steady states: Halvorsen–Skogestad 1999; Malinen–Tanskanen 2007; Doherty–Perkins 1982; Jacobsen–Skogestad 1991; Acrivos
and Amundson 1955. Mathematics: Shivakumar–Chew 1974; Berman–Plemmons 1994.

## 5. Not checked

Paywalled or print-only items listed in `total/literature/novelty.md` §5 (Underwood 1945–1946, Petlyuk–Serafimov 1983,
Glinos–Malone, Koehler et al., IDEAS papers, reviews behind paywalls), the full text of Agrawal–Tumbalam Gooty 2020,
Malinen–Tanskanen 2007 (abstract and conference summary only), Russian and Chinese databases, Google Scholar (not
accessible). Semantic Scholar was rate-limited. None of these is likely to contain a network-level proof, given that
the same research groups describe the conjecture as open as late as 2021, but this cannot be excluded.
