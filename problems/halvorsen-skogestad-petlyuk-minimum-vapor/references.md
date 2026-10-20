# References and novelty check

## Primary source (open statement)

- I. J. Halvorsen, S. Skogestad, *Minimum energy consumption in multicomponent distillation. 3. More than three
  products and generalized Petlyuk arrangements*, Ind. Eng. Chem. Res. 42 (2003) 616–629,
  doi:10.1021/ie0108651. Read in full (author copy on Skogestad's NTNU page and publisher PDF). Open statements:
  abstract ("We also conjecture that the minimum energy requirement for the generalized extended Petlyuk
  arrangement is lower than the minimum-energy requirement for any distillation configuration when we consider
  conventional adiabatic sections and no internal heat exchange."), §6.3 ("We have not presented a complete proof,
  so the above conclusion is a conjecture."), and §5.3 ("We have not carried out a detailed proof for the general
  N-component M-product case for other than the far left and right peaks."), p. 624.

## Sources used

| Reference | Identifier | What was read |
|---|---|---|
| Halvorsen & Skogestad, Part 1, Vmin diagram for a two-product column, Ind. Eng. Chem. Res. 42 (2003) 596–604 (Part 3's reference list swaps the page ranges of Parts 1 and 2) | doi:10.1021/ie010863g | Sections 1–3 (Underwood equations, Vmin diagram) |
| Halvorsen & Skogestad, Part 2, three-product Petlyuk arrangements, Ind. Eng. Chem. Res. 42 (2003) 605–615 | doi:10.1021/ie0108649 | carry-over of roots, junction analysis |
| Halvorsen & Skogestad, Part 3 (primary) | doi:10.1021/ie0108651 | full text |
| I. J. Halvorsen, *Minimum energy requirements in complex distillation arrangements*, Dr. ing. thesis, NTNU 2001:43 | open access, skoge.folk.ntnu.no | Ch. 2.5, 3.2 (Underwood, active roots, computation procedure), 4.3–4.4 (reverse flows, carry-over), 6.2 (junction lemma, conjecture) |
| R. Tumbalam Gooty, R. Agrawal, M. Tawarmalani, *Advances in MINLP to identify energy-efficient distillation configurations*, Oper. Res. 72, 639–659; arXiv:2010.12113 (v1, 2020) | doi:10.1287/opre.2022.2340, arXiv:2010.12113 | §3–4 in full (configuration space, mass balances, Underwood constraints, numbered (A18)–(A23) in the arXiv v1 PDF, monotonicity cuts) — this is the model of the paper (verifier's correction: the numbering was given as (B1)–(B6)) |
| V. H. Shah, R. Agrawal, *A matrix method for multicomponent distillation sequences*, AIChE J. 56 (2010) 1759–1775 | doi:10.1002/aic.12118 | configuration counts (Table 1), read by Agent C |
| V. H. Shah, R. Agrawal, *Are all thermal coupling links between multicomponent distillation columns useful from an energy perspective?*, Ind. Eng. Chem. Res. 50 (2011) 1770–1777 | doi:10.1021/ie101768c | full text: some thermal couplings bring no first-law benefit (never a loss in their examples); no statement contradicting our result |
| Z. Fidkowski, L. Królikowski, *Minimum energy requirements of thermally coupled distillation systems*, AIChE J. 33 (1987) 643–653 | doi:10.1002/aic.690330412 | skimmed: closed-form ternary comparisons (Petlyuk vs side strippers/rectifiers and simple sequences) |
| Z. T. Fidkowski, R. Agrawal, *Multicomponent thermally coupled systems of distillation columns at minimum reflux*, AIChE J. 47 (2001) 2713–2724 | doi:10.1002/aic.690471211 | abstract and FTC minimum-reflux results (read by Agent C) |
| U. Nallasivam, V. H. Shah, A. A. Shenvi, J. Huff, M. Tawarmalani, R. Agrawal, *Global optimization of multicomponent distillation configurations: 2. Enumeration based global minimization algorithm*, AIChE J. 62 (2016) 2071–2086 | doi:10.1002/aic.15204 | model equations 10–13 (read by Agent C); verifier: Eqs. 14, 17–18b use an inequality V ≥ Underwood value at every root plus explicit recovery-ordering constraints, not TG's equalities at the interior roots, so the two formulations are related but not identical |
| T. J. Mathew, R. Tumbalam Gooty, M. Tawarmalani, R. Agrawal, *Thermal coupling via heat transfer: a potential route to simple distillation configurations with lower heat duty*, Ind. Eng. Chem. Res. 58 (2019) 21671–21678 | doi:10.1021/acs.iecr.9b04689 | abstract: beats the FTC only with heat exchange between streams (excluded by the conjecture) |
| A. J. V. Underwood, *Fractional distillation of multicomponent mixtures*, Chem. Eng. Prog. 44 (1948) 603–614 | — | not read; cited through Halvorsen's thesis |
| N. A. Carlberg, A. W. Westerberg, *Temperature–heat diagrams for complex columns* 2 and 3, Ind. Eng. Chem. Res. 28 (1989) 1379–1386 and 1386–1397 | doi:10.1021/ie00093a017, doi:10.1021/ie00093a018 | not read; cited through HS Part 2–3 (carry-over of roots) |

## Novelty check (2026-10-20, Agent T)

Question: has the Halvorsen–Skogestad conjecture (or the general-N minimum-vapor formula for the generalized
Petlyuk arrangement, §5.3 of Part 3) been proved or refuted since 2003?

1. **OpenAlex**, works citing HS Part 3 (W2086353211): 119 works (1 page, sorted by year). Titles of all and
   abstracts of the 2018–2026 ones scanned for "conjecture", "proof/prove", "lowest", "minimum vapor/energy":
   no proof. The Agrawal group states the property as known: e.g. Kadu, Nogaja, Tawarmalani, Agrawal, Chem. Eng.
   Res. Des. (2026), doi:10.1016/j.cherd.2026.09.036 ("the FTC design minimizes vapor duty"); Mathew et al.
   2019 (above). Halvorsen & Skogestad 2025, Sep. Purif. Technol., doi:10.1016/j.seppur.2025.132530 (liquid-only
   transfer, analytical minimum vapor equal to the standard DWC): no proof of the conjecture.
2. **Semantic Scholar** citations of doi:10.1021/ie0108651: 123 records, 2019–2026 titles scanned: same set,
   nothing new.
3. **Web search** (3 queries: "Halvorsen" "Skogestad" conjecture generalized Petlyuk minimum vapor proof; proof
   fully thermally coupled configuration minimum vapor duty among all distillation configurations; "fully thermally
   coupled" "lowest" minimum vapor duty proof 2024–2026): no proof. Secondary sources misstate the status
   (e.g. Caballero & Grossmann, Comput. Chem. Eng. 61 (2014) 118–135, doi:10.1016/j.compchemeng.2013.10.015, as quoted in a search snippet: "Halvorsen & Skogestad (2003c) proved ...").
4. **arXiv API** (Petlyuk; "thermally coupled" AND distillation; Underwood AND distillation; "minimum vapor" AND
   distillation; "dividing wall"): no relevant preprint.
5. **Zenodo API** (Petlyuk; Halvorsen Skogestad; fully thermally coupled distillation; minimum vapor distillation
   configuration; Underwood minimum vapor): no relevant record.
6. **GitHub**: `gh search repos` (Petlyuk; Underwood distillation; minimum vapor distillation; thermally coupled
   distillation): nothing relevant. AI campaigns: `gh api users/ipitchford/repos --paginate` (171 repositories,
   none on distillation), `gh api users/mmislan/repos` (24 repositories, none).

Conclusion: open for N >= 4 (and, as HS Part 3 §5.3 states, even the lower bound "minimum vapor of the
generalized Petlyuk arrangement = highest peak" had been proved only for the extreme peaks). No competing
AI-generated candidate found.

## Novelty re-check (2026-10-20, independent verifier)

- OpenAlex `cites:W2086353211`: 119 works; the 16 works dated 2024–2026 were listed and their abstracts scanned
  (keywords conjecture / proof / prove / lowest / Petlyuk / minimum vapor / fully thermally): no proof or refutation;
  the 2026 Kadu et al. (doi:10.1016/j.cherd.2026.09.036) and Halvorsen & Skogestad 2025
  (doi:10.1016/j.seppur.2025.132530) do not address the conjecture.
- Web search ("proof generalized Petlyuk arrangement minimum vapor lowest among all distillation configurations
  Halvorsen Skogestad conjecture"): only the original papers and NTNU pages, which call it a conjecture.
- Zenodo API (Petlyuk; thermally coupled distillation minimum vapor; Halvorsen Skogestad conjecture; Underwood
  minimum vapor proof): no relevant record.
- GitHub: `gh api users/ipitchford/repos --paginate` (171 repositories; the only names containing "column" or
  "distill" are mathematics/quantum-information repositories), `gh api users/mmislan/repos` (24, none relevant),
  `gh search repos Petlyuk` and `"Halvorsen Skogestad"`: nothing relevant.
- Additional sources read by the verifier: Halvorsen thesis §3.2.3–3.2.4 (active roots, definition of a
  distributing component), §4.2.1, §4.3.3–4.3.4 (reverse flows in directly coupled columns), Table 6.1;
  TG arXiv source `template.tex` (§4.4 Underwood constraints, test set, Table of eight optima); Shah & Agrawal 2010
  (definition of basic configurations, counts); Nallasivam et al. 2016 (Underwood constraints, Eqs. 14–23).

Conclusion unchanged: no prior proof or refutation found.

## Version 2 (wave 2, 2026-10-20): additional sources

| Reference | Identifier | What was read |
|---|---|---|
| I. J. Halvorsen, S. Skogestad, *Minimum energy and entropy requirements in multicomponent distillation*, Distillation & Absorption 2002, Baden-Baden, paper 6-27 | NTNU proceedings copy (skoge.folk.ntnu.no/prost/proceedings/distillation02/dokument/6-27.pdf) | full text: conjecture restated for adiabatic arrangements at constant pressure without internal heat integration; infinitely many stages assumed "because this provides a lower bound"; total vaporization as the energy measure; internal heat exchange case "unknown" |
| Halvorsen thesis (above), additional parts | — | §2.6.1 (CMO energy balance), §3.4.3–3.4.8 (pinch zones, finite stages), §4.3.1–4.4 (directly coupled sections, reverse flows, carry-over), §6.1–6.6 (definitions of adiabatic and non-adiabatic sections, minimum energy, internal heat exchange examples, Table 6.1, several pressures) |
| Z. T. Fidkowski, *Distillation configurations and their energy requirements*, AIChE J. 52 (2006) 2098–2106 | doi:10.1002/aic.10803 | full text: configuration generation, ten types of interconnecting streams, stacking; minimum vapor by Underwood; transition splits assumed optimal ("not proven") outside the FTC |
| R. Agrawal, *Synthesis of multicomponent distillation column configurations*, AIChE J. 49 (2003) 379–401 | doi:10.1002/aic.690490210 | configuration space (scanned for the definition of the class) |
| R. Agrawal, *Synthesis of distillation column configurations for a multicomponent separation*, Ind. Eng. Chem. Res. 35 (1996) 1059–1071 | doi:10.1021/ie950323h | abstract and §1–2 (4n − 6 sections, superstructures) |
| R. Agrawal, Z. T. Fidkowski, *Are thermally coupled distillation columns always thermodynamically more efficient for ternary distillations?*, Ind. Eng. Chem. Res. 37 (1998) 3444–3454 | doi:10.1021/ie980062m | introduction: FC has the lowest heat demand among five ternary configurations (proved earlier for those five); thermodynamic efficiency comparison |
| A. C. Christiansen, S. Skogestad, K. Lien, *Complex distillation arrangements: extending the Petlyuk ideas*, Comput. Chem. Eng. 21 (1997) S237–S242 | doi:10.1016/S0098-1354(97)00055-0 | scanned: extended Petlyuk arrangements, degrees of freedom |
| I. P. Nikolaides, M. F. Malone, *Approximate design of multiple-feed/side-stream distillation systems*, Ind. Eng. Chem. Res. 26 (1987) 1839–1845 | doi:10.1021/ie00069a020 | scanned: minimum reflux of multi-feed columns (controlling feed) |
| N. A. Carlberg, A. W. Westerberg, *Temperature–heat diagrams for complex columns* 2 and 3, Ind. Eng. Chem. Res. 28 (1989) 1379–1386 and 1386–1397 | doi:10.1021/ie00093a017, doi:10.1021/ie00093a018 | scanned (now available): Underwood's method for side strippers/enrichers and the Petlyuk configuration, root carry-over |
| N. L. Franklin, J. S. Forsyth, *The interpretation of minimum reflux conditions in multi-component distillation*, Trans. IChemE 31 (1953) 363 (reprint Trans IChemE 75 (1997) S56–S81) | — | introduction and geometric interpretation of Underwood's equations (single columns) |
| F. B. Petlyuk, *Distillation Theory and its Application to Optimal Design of Separation Units*, Cambridge University Press 2004, ch. 5 (pp. 108–169) | doi:10.1017/CBO9780511547102.007 | scanned: trajectories of infinite columns at finite reflux |
| J. Bausa, R. v. Watzdorf, W. Marquardt, *Shortcut methods for nonideal multicomponent distillation: 1. Simple columns*, AIChE J. 44 (1998) 2181– | doi:10.1002/aic.690441008 | abstract: rectification body method (single columns) |
| A. J. V. Underwood, *Multicomponent distillation — number of transfer units*, Ind. Eng. Chem. 41 (1949) | doi:10.1021/ie50480a044 | not used |
| R. Shiras, D. N. Hanson, C. H. Gibson, *Calculation of minimum reflux in distillation columns*, Ind. Eng. Chem. 42 (1950) | doi:10.1021/ie50485a034 | not used |
| Z. Jiang, M. Tawarmalani, R. Agrawal, *Minimum reflux calculation for multicomponent distillation in multi-feed, multi-product columns: mathematical model*, AIChE J. 68 (2022) e17929 | doi:10.1002/aic.17929 | abstract, introduction, structure (found by the targeted novelty search): single MFMP columns, not a network bound |

## Version 2: targeted novelty check of the stage-level argument (Agent F, 2026-10-20 14:10; superseded by the literature agent's report below)

(Queries and outcomes: see the list in this section, copied from Agent F's log.) Halvorsen thesis §2.5.2 (Underwood's
stage identity, eq. 2.49, used only at a section's own roots), §3.4, §6.2.1, §6.6.4; HS Distillation'02; Fidkowski–
Królikowski 1987, Fidkowski–Agrawal 2001, Fidkowski 2006, Agrawal–Fidkowski 1998, Christiansen–Skogestad–Lien 1997,
Carlberg–Westerberg 1989, Nikolaides–Malone 1987, Franklin–Forsyth 1953, Petlyuk 2004 ch. 5, Bausa et al. 1998; Jiang
et al. 2022; OpenAlex citing works of HS Part 3 filtered by "proof" (4), "prove" (16), "lower bound" (9), "bound" (12),
abstracts of the relevant ones read (doi:10.1002/aic.18125, 10.1002/aic.18895, 10.1021/acs.iecr.1c01036); OpenAlex
free searches (4); arXiv API (5 queries); Zenodo (4); GitHub (ipitchford, mmislan, `gh search repos Petlyuk`); web
searches (3). Outcome: neither the argument nor the theorem found.

## Version 2: targeted novelty search of the stage-level argument (literature agent, Claude Opus 5.5, effort xhigh)

Full report: `total/literature/novelty.md` (stopped early for quota reasons; unchecked items listed there). Verdict:
neither the argument (Underwood flux conserved at every unit, evaluated at the feed roots, sign uniformity at
equilibrium stages, summation over the units of one sign) nor the theorem (every arrangement of a broad class needs at
least max_j H_j, N ≥ 4) was found. The same groups still describe the conjecture as open (HS 2002/2003; Fidkowski–
Agrawal 2001; Fidkowski 2006; T. J. Mathew, Purdue thesis 2021). Three secondary sources state it as proved without
proof (Halvorsen–Skogestad 2011; Caballero–Grossmann 2014; Nallasivam et al. 2016): misattributions, not prior proofs.
Older ingredients: Underwood's stage identity (Underwood 1948; Halvorsen thesis eqs. 2.45–2.49); sign of Underwood
coordinates relative to section roots (Franklin–Forsyth 1953); stage-level sign comparisons inside one section
(Jiang–Tawarmalani–Agrawal, AIChE J. 2022, Prop. 1); cut/vapor-balance arguments at one junction (HS 2003 §3;
Halvorsen thesis §6.2.1).
