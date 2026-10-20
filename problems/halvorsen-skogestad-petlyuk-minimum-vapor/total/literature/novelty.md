# Novelty check — stage-level proof of the Halvorsen–Skogestad conjecture

Agent: Claude Opus 5.5 (effort xhigh), role: research librarian (targeted novelty search). Date: 2026-10-20,
14:22–15:10. Brief: `BRIEF-literature.md`. The search was stopped early on the orchestrator's instruction (owner's
quota), so some planned checks were not done. They are listed in §5.

## 0. Verdict

**(A) Neither the argument nor the full theorem appears in the literature I could check.** Within that limit, no
published or posted work does either of the following:
- proves that every adiabatic arrangement (or every arrangement of any broad class with N ≥ 4) needs at least
  max_j H_j;
- uses the Underwood flux φ(θ) = V − Σ α_i n_i/(α_i − θ) as a quantity conserved at every unit of a network,
  evaluated at the feed roots θ_j, together with sign uniformity at equilibrium stages and a summation over the
  units of one sign.

The same authors and groups still describe the conjecture as open:
- Halvorsen–Skogestad 2002 and 2003;
- Petlyuk himself, through a private communication cited by Halvorsen;
- Fidkowski–Agrawal 2001 and Fidkowski 2006;
- the Agrawal group as late as the 2021 Purdue thesis of T. J. Mathew: "For four and more components, the
  conjecture is yet to be proved or disproved."

Three secondary sources state the result as proved or shown, without proof:
- Halvorsen–Skogestad 2011;
- Caballero–Grossmann 2014;
- Nallasivam et al. 2016, which calls the FTC vapor "a mathematical lower bound".

These sources should be cited as misattributions, not as prior proofs.

The individual ingredients are older (§3):
- **(a) Underwood's stage identity at a free parameter.** Underwood 1948; Halvorsen thesis eq. 2.45.
- **(b) The equilibrium relation behind the sign rule.** It is implicit in the derivation of Halvorsen's eq. 2.49.
  Franklin–Forsyth 1953 use the sign of Underwood coordinates relative to the section roots.
- **(c) Stage-level sign comparisons of an Underwood-type function.** Jiang 2019 thesis Prop. 5.4.1, published as
  Jiang–Tawarmalani–Agrawal, AIChE J. 2022, Prop. 1. These are used only inside one column section, at the
  section's own stage roots, to prove convergence to a pinch.
- **(d) Cut / vapor-balance arguments with Underwood section bounds.**
  - HS 2003 Part 3 §3: equality case for the Petlyuk.
  - Halvorsen thesis §6.2.1: one junction.
  - Mathew 2021 thesis §6.3: one four-component family.

I found no work that combines (a) and (b) at the feed roots over a whole network. This is the step that gives the
bound for arbitrary arrangements (recycles, splitters, multi-feed columns, reverse flows, finite stages, impure
products).

If ingredient-level overlap counted, the verdict would be "B for the ingredients, A for the argument and the
theorem". I recommend that the paper cite (a)–(d) explicitly (§3).

## 1. Sources on disk (`certificates/sources/`, read with pdftotext)

| Source | What was checked | Outcome |
|---|---|---|
| HS 2003 Part 3 (ie0108651) | abstract, §1, §3, §5.3, §6 | Conjecture stated (abstract p. 616; §6.3 p. 626). §3 (pp. 617–618) evaluates Underwood's equation at the common root θ_j on each intersection line I_j of the **Petlyuk only** (equality case). §5.3 p. 624: "We have not carried out a detailed proof for the general N-component M-product case for other than the far left and right peaks." §6, p. 624: "no proof has been found in the literature,16 except for the ternary case"; ref. 16 is "Petlyuk, F. B. Private communication, 2000". |
| HS 2003 Parts 1–2 | §3 (Underwood definitions), junction analysis | Section-level material balances at cross-sections only. No bound for general arrangements. |
| Halvorsen thesis 2001 (NTNU 2001:43) | §2.5.1–2.5.2, §3.2, §6.1–6.2.3 | Eq. 2.45 (p. 52) is Underwood's stage identity for a free parameter φ: "The parameter φ is free to choose". It is then used only at the section roots (eq. 2.46–2.50, pp. 52–53). §6.2.1 (pp. 176–179) is the junction lemma: direct coupling minimizes max(V_I1, V_I2) through **one** junction, with the upstream column fixed. p. 180: "We have not presented a complete proof, so the above conclusion is a conjecture." p. 176: "no proof has been found in the literature (Petlyuk 2000), except for the ternary case". |
| Fidkowski–Królikowski 1987 (aic.690330412) | full text | Closed-form Underwood optimization of four ternary thermally coupled systems plus DS/IS, using Vieta's formulas. p. 651: "thermally coupled systems are ordered as follows: TCS = TCS-I < TCS-R = TCS-S". No general bound. |
| Fidkowski–Agrawal 2001 (aic.690471211) | intro, pp. 2720–2723 | p. 2713 overstates FK 1987 ("minimum vapor flow in the FC configuration is smaller than in any other configuration of columns for ternary distillation"). p. 2720: "no proof can easily be given that the transition split operation is optimal for quaternary". p. 2723: "have not yet developed a rigorous proof". |
| Fidkowski 2006 (aic.10803) | pp. 2103–2105 | p. 2103: transition splits optimal for the FTC; "For other configurations it has not been proven yet". p. 2104: "(but not proven)". |
| Agrawal–Fidkowski 1998 (ie980062m) | abstract, intro | p. 3444: FC "known to have the lowest heat demand for ternary distillation". "Several studies have proven" it is lowest among the **five** ternary configurations. |
| Christiansen–Skogestad–Lien 1997 | full text | Nothing relevant ("we conjecture" refers to flow directions). |
| Carlberg–Westerberg 1989 Parts 2–3 | intro, results | Underwood equations for side strippers/enrichers and the Petlyuk; root carry-over. No bound for general arrangements. |
| Underwood 1949 (ie50480a044) | full text | Transfer-unit integration with E(x, φ). No bound. |
| Franklin–Forsyth 1953 (reprint 1997) | pp. S56–S62 | Underwood coordinates η_i = Σ α x/(α − φ_i). Sign/feasibility arguments on flats η_i = 0 within **one** section, relative to section roots (p. S62). Closest classical "sign" use; not at feed roots, not a network bound. |
| Shiras–Hanson–Gibson 1950; Nikolaides–Malone 1987 | grep + skim | Single (multi-feed) column minimum reflux only. |
| Petlyuk 2004 ch. 5 | grep + §5.2 | Underwood system; "the proof of equality of parameters" (common roots, p. 113). No network bound. |
| Bausa–v. Watzdorf–Marquardt 1998 | grep | Rectification body method, single columns. |
| Agrawal 1996 / 2003; Shah–Agrawal 2010 / 2011 | grep | Shah–Agrawal 2010 p. 1760 (FTC): "it always has the lowest minimum total vapor duty", citing FK 1986 and HS 2003 Part 3. Neither proves it for N ≥ 4. |
| Nallasivam et al. 2016 (aic.15204) | §"Ensuring global optimality" | p. 2079: FTC "has the lowest vapor duty in the search space.13,19,20 ... its vapor flow requirement is a mathematical lower bound for the total minimum vapor duty requirement for any configurations in the search space". Refs 13 = FA 2001, 19 = FK 1986, 20 = HS Part 3: none proves it. A misattribution. |
| Tumbalam Gooty et al. (arXiv 2010.12113) | grep | The FTC objective is used as a feasible point / upper reference (eq. 40), not as a proved lower bound. |

## 2. Queries and outcomes (2026-10-20)

**Citation databases**
- **OpenAlex, works citing HS Parts 1–3** (W2084600958, W2037012710, W2086353211): 235 unique works. Abstracts
  reconstructed and scanned for proof / prove / lower bound / any configuration / arbitrary / conjecture /
  minimum vapor / lowest / Underwood / general. Outcome: no proof, no lower bound for general arrangements, no
  counterexample without heat integration.
  - Abstract statements of FTC optimality as known: Kadu et al. 2026 CHERD doi:10.1016/j.cherd.2026.09.036
    ("While the FTC design minimizes vapor duty").
  - Beating the FTC only with heat exchange (TCH), excluded by HS: Mathew et al. 2019 IECR
    doi:10.1021/acs.iecr.9b04689.
- **OpenAlex, works citing FK 1986, FK 1987, Agrawal–Fidkowski 1998, FA 2001 and the Halvorsen thesis**: 359
  works, keyword scan. New relevant items:
  - Wachter–Ko–Andres 1988 (continuum MFMP minimum reflux, single column);
  - Koehler–Poellmann–Blass 1995 review (not accessed);
  - Lucia–McCallum 2010 and Lucia–Amale–Taylor 2008 (shortest-stripping-line heuristics, no proof);
  - Jiang–Tawarmalani–Agrawal 2022 (below).
- **OpenAlex free-text searches**: "Underwood equations lower bound minimum vapor distillation configuration";
  "minimum vapor flow any distillation arrangement proof"; "minimum energy distillation network arbitrary
  configuration lower bound ideal mixture"; "Petlyuk thermodynamically reversible distillation minimum energy";
  "Timoshenko Serafimov thermally coupled distillation complexes energy"; "Platonov Petlyuk rectification optimal
  cascade"; "IDEAS infinite dimensional state space distillation network minimum utility cost"; year-restricted
  2023–2026 searches ("Petlyuk minimum vapor", "Underwood minimum reflux", "fully thermally coupled distillation
  minimum energy"). Outcome: nothing containing the argument or the theorem. Numerical global-synthesis work on
  ternary stage networks (IDEAS, Takase–Hasebe 2018) finds Petlyuk-like optima; it is numerical, N = 3, not
  accessed in full.
- **Semantic Scholar**: citations of Parts 3 and 1, 182 works; the 24 not in OpenAlex were scanned by title.
  Outcome: nothing new. Keyword search: rate-limited (HTTP 429). One of four queries returned results, none
  relevant.
- **Crossref**: four bibliographic queries plus DOI lookups. Outcome: noise; the useful hits are listed in §5.

**Preprints, data and code**
- **arXiv API**, 11 queries (Underwood AND distillation; Petlyuk; "thermally coupled" AND distillation; "dividing
  wall"; "minimum reflux"; "minimum vapor" AND distillation; "distillation configurations"; "distillation
  sequences"; ...). Only relevant items: arXiv:2010.12113 (TG) and arXiv:2503.20917 (Jiang et al., MFMP
  algorithms). The latter has no network bound and no Petlyuk/FTC optimality claim (checked on the HTML
  version).
- **Zenodo API**, 15 queries (Petlyuk; "Halvorsen-Skogestad"; "generalized Petlyuk"; "Underwood roots";
  "Vmin diagram"; "fully thermally coupled"; "minimum vapor" AND Underwood; ...). Nothing relevant.
- **GitHub**:
  - `gh api users/ipitchford/repos --paginate`: 171 repositories, none on distillation. The only matches on
    "column"/"separation" are mathematics repositories.
  - `gh api users/mmislan/repos`: 24 repositories, none relevant.
  - `gh search repos` (Petlyuk; Underwood distillation; Halvorsen Skogestad; minimum vapor distillation;
    thermally coupled distillation; dividing wall column; Vmin diagram): nothing relevant.
  - `gh search code` (Petlyuk / Underwood / distillation in user:ipitchford; Petlyuk in user:mmislan; "Halvorsen
    Skogestad"): only this project's own public `references.md`.

**Web searches** (about 20, standard and extended, in English, Russian and Chinese): proof of Petlyuk minimality;
Underwood as a conserved quantity; Russian school (Petlyuk–Platonov–Slavinskii 1965, Petlyuk–Serafimov 1983,
Serafimov, Timoshenko); Koehler–Aguirre–Blass; Glinos–Malone; Jiang–Tawarmalani–Agrawal MFMP; HS 2011; DWC
reviews (Yildirim–Kiss–Kenig, Dejanović–Olujić); IDEAS. Outcome: no proof or refutation found. The Russian
monograph and 1960s papers could not be accessed.

**Theses** (open access):
- NTNU: Halvorsen 2001 (on disk); Christiansen 1998, Dwivedi 2013, Ghadrdan 2014 (grep: nothing relevant).
- Ulm: Ränger 2021 (§2.1, page not recorded: "It is often stated that this type of column has the lowest energy
  requirement of all known distillation variants").
- Purdue: Jiang 2019, Mathew 2021, Tumbalam Gooty 2020, Nogaja 2025 (below).

**Other open documents**: HS Distillation & Absorption 2002 paper 6-27; HS 2011 J. Nat. Gas Sci. Eng. (in-press
author copy); Caballero–Grossmann 2014 (accepted manuscript); Jiang–Tawarmalani–Agrawal 2022 (author copy).

## 3. Closest prior results (exact quotes, page numbers)

1. **Underwood's identity for a free parameter** (ingredient a). Halvorsen thesis eq. 2.45, p. 52: "The parameter
   φ is free to choose, and the Underwood roots are defined as the values of φ which make the left hand side of
   (2.45) unity". The derivation of eq. 2.49 (p. 53) uses y_i = α_i x_i/Σα x, the same algebra as our
   φ_vap/φ_liq = θ/Σα_k x_k. The identity is used only at section roots, to relate stage compositions inside one
   section (eq. 2.50). Origin: Underwood 1945–1948 (not accessed; print-only, no DOI).

2. **Stage-level sign argument with an Underwood-type function** (closest in technique). Z. Jiang, PhD thesis,
   Purdue 2019 ("Multicomponent Distillation – Mathematical Modeling, Global Optimization, and Process
   Intensification"), §5.4, Prop. 5.4.1, pp. 151–152. Same as Jiang, Tawarmalani, Agrawal, AIChE J. 68 (2022)
   e17929, doi:10.1002/aic.17929, §4.1, Prop. 1, p. 10 of 21.
   - Method: f_n(x) = Σ v_{i,n}/(α_i − x), i.e. −φ_vap,n(x)/x in our notation. The section relation
     x f_n(x) − (L/Σ_j v_{j,n+1}/α_j) f_{n+1}(x) = x f_{n−1}(x) − (L/Σ_j v_{j,n}/α_j) f_n(x) expresses
     conservation of the flux inside one section at an arbitrary x, combined with the VLE relation. It is
     evaluated at the stage root x_k(n): "this implies that f_{n+1}(x_k(n)) f_{n−1}(x_k(n)) < 0 since x_k(n) > 0"
     (thesis p. 151).
   - Purpose: monotone convergence of the stage roots to the pinch roots, inside one section of one MFMP column.
   - Difference from our argument: no feed root, no network, no energy balance, no bound.

3. **Sign of Underwood coordinates** (classical sign use). Franklin & Forsyth, Trans. IChemE 31 (1953), reprinted
   Trans IChemE 75 (1997) S56–S69, p. S62. The η-coordinates of x_D "will have the same sign as those of x_D",
   and distillation lines cannot cross flats of the wrong sign. This uses section roots in a single column.

4. **Junction and cross-section arguments at the common roots** (section level).
   - HS 2003 Part 3 §3, pp. 617–618: vapor through intersections I1–I3 of the **Petlyuk** computed with the
     common root θ_j, where "there is one common active Underwood root". This is an equality for one arrangement.
   - Halvorsen thesis §6.2.1, p. 179: "minimization of the vapour rate through any intersection (I1 or I2) is
     found when the common Underwood roots in column C21 equal the actual roots in the top section of C1". This
     covers one junction. p. 180: "This result is based on the simple argument that at any junction where we might
     consider another type of connection than the direct coupling, the required vapour flow through the junction,
     and thereby through a cross-section of the whole arrangement, will increase. We have not presented a
     complete proof, so the above conclusion is a conjecture."

5. **A section-level cut proof for one four-component family**. T. J. Mathew, PhD thesis, Purdue 2021 ("Advances
   in Modelling and Optimization of Multicomponent Distillation Configurations"), §6.3, pp. 127–134.
   - Result: the symmetric-split family of four-component configurations, including a heat exchange between
     condenser C and reboiler B, cannot beat the FTC.
   - Tools: Underwood section bounds, root ordering under thermal coupling, and a vapor balance (eq. 6.12, p. 132:
     "The same result is obtained by applying a vapor balance on a suitably constructed boundary in the
     configuration").
   - Status statements: p. 125, "For four and more components, the conjecture is yet to be proved or disproved";
     p. 134, "the conjecture is not yet proved for four-component configurations".
   - This is the nearest published relative of result 12's cut identity, much narrower in scope.

6. **Ternary results**. FK 1987 (p. 651), FK 1990 (AIChE J. 36, 1275–1278, doi:10.1002/aic.690360820, not
   accessed) and Agrawal–Fidkowski 1999 (AIChE J. 45, 485–496, doi:10.1002/aic.690450306, not accessed).
   Mathew 2021 p. 125 says these "analytically proved the conjecture when only the traditional TCM are employed"
   for three components. These are comparisons among regular-column ternary configurations, not arbitrary
   arrangements.

**Statements of the conjecture as open (primary)**
- HS D&A 2002, paper 6-27, p. 1: "yet it remains unsolved even for the case of ideal mixtures". p. 6: "this is a
  conjecture since it has not been generally shown for all possible column arrangements".
- HS 2003 Part 3, pp. 616, 624, 626 (quoted above).

**Misattributions (secondary; to be cited as such, if at all)**
- HS, "Energy efficient distillation", J. Nat. Gas Sci. Eng. 3 (2011) 571–580, p. 576: "It is also shown that an
  optimally operated generalized Petlyuk arrangement results in the lowest overall vapour flow requirement for
  any distillation configuration when we consider constant pressure and no external heat integration." No proof
  is given; it refers back to the 2003 papers, which call it a conjecture.
- Caballero & Grossmann, Comput. Chem. Eng. 61 (2014) 118–135, doi:10.1016/j.compchemeng.2013.10.015, accepted
  manuscript p. 4 (Introduction): "Halvorsen & Skogestad (2003c) proved that the minimum energy consumption for an
  ideal N-component mixture is always obtained for the FTC configuration."
- Z. Jiang, Purdue thesis 2019, p. 27: "it is conjectured [FA 2001, HS 2003] and later demonstrated [Nallasivam
  2016] that the FTC configuration indeed has the maximum vapor duty benefits compared to any other configuration
  in the search space". "Demonstrated" means numerical enumeration. Also p. 97: "known to always have the lowest
  total heat duty requirement among all configurations".
- Nallasivam et al. 2016, p. 2079 (quoted in §1); Shah–Agrawal 2010, p. 1760.
- A. S. Nogaja, Purdue thesis 2025: "It is known to have the lowest heat duty among all non-heat integrated
  configurations [HS 2003 Part 3]".

## 4. Assessment

The theorem, as stated in the brief, has not been proved in any source I checked: R ≥ max_j H_j for every
arrangement (networks, recycles, splitters, multi-feed columns, finite stages, impure products up to O(ε)), and
hence optimality of the Petlyuk value. The originators described it as a conjecture (2002, 2003). Petlyuk is
reported (private communication, 2000) as knowing no proof beyond the ternary case. The Agrawal group described it
as open for N ≥ 4 in 2021.

The specific argument has not been found either: the stream flux as a network-conserved quantity at θ_j, sign
uniformity per equilibrium stage, and summation over 𝒰₊. Its pieces are classical identities, used before only
within single sections at section roots: Underwood 1948 / Halvorsen eq. 2.45; Franklin–Forsyth 1953; Jiang 2019 /
2022.

## 5. What could not be checked

- Paywalled or print-only items not read. The DOIs are listed under "PDF requests" in PROGRESS.md:
  - Underwood 1945, 1946a, 1946b (J. Inst. Pet.) and 1948 (Chem. Eng. Prog.): no DOIs.
  - Petlyuk–Platonov–Slavinskii 1965 (Int. Chem. Eng. 5, 555).
  - Petlyuk–Serafimov 1983 monograph (Russian).
  - Petlyuk 2004 ch. 4, 6, 7.
  - Koehler–Poellmann–Blass 1995 review.
  - Glinos–Malone 1985 and 1988; Glinos–Nikolaides–Malone 1986.
  - FK 1986 and 1990; Agrawal–Fidkowski 1999.
  - Offers–Düssel–Stichlmair 1995.
  - Drake–Manousiouthakis 2002 (IDEAS); Takase–Hasebe 2018 (IDEAS, ternary).
  - Mathew et al. 2019 IECR (possible published version of the symmetric-split proof).
  - Agrawal–Tumbalam Gooty 2020 perspective.
  - Reviews: Yildirim–Kiss–Kenig 2011; Dejanović–Matijašević–Olujić 2010; the 2025 Elsevier book "Thermally
    Coupled Distillation".
  - Ryu–Maravelias 2021 (generalized network synthesis).
  - HS 2025 LOT paper.
  - Jiang–Agrawal 2019 review (open access, download failed).
- Textbooks not consulted: King 1980, Stichlmair–Fair 1998, Doherty–Malone 2001, Fonyó 1974, Grunberg 1956.
- Semantic Scholar keyword search: mostly rate-limited.
- Google Scholar: not accessible.
- The search was stopped early on the orchestrator's instruction. Not done: Russian-language databases (eLIBRARY),
  Chinese databases (CNKI), and a systematic scan of the Distillation & Absorption proceedings after 2002.
- Abstract-level only: works without OpenAlex abstracts (about 40% of the citing set) were checked by title only.
