# References

## Primary source (read in full)

- **[JS14]** J. Jäschke, S. Skogestad, *Optimal operation of heat exchanger networks with stream split: Only
  temperature measurements are required*, Computers and Chemical Engineering 70 (2014) 35–49,
  doi:10.1016/j.compchemeng.2014.03.020. Manfred Morari special issue.
  - Publisher PDF, 15 pages, read in full:
    https://skoge.folk.ntnu.no/publications/2014/jaschke-hen-paper-cce-morari-issue/impress.pdf.
  - Pages 39, 47 and 48 were also checked as rendered images: Theorem 2/3, eq. (18), (A.10)–(A.16), (B.1), and
    Tables B.6–B.7.
  - Also downloaded: the in-press version `inpress0.pdf` (same content) and an older preprint `old-version/HenCCE.pdf`.
    In the preprint, Theorems 2/3 were claimed "proved for Mj = 1, 2, 3 (see Section 6), and conjectured" beyond.
  - Local copies are in WSL under `~/work/jaschke-skogestad/dl/`.

## Related works by the same authors (read)

- J. Jäschke, S. Skogestad, *A self-optimizing strategy for optimal operation of a preheating train for a crude oil
  unit*, Proc. ESCAPE-24, Budapest, 2014 (Computer Aided Chemical Engineering 33).
  - PDF: https://skoge.folk.ntnu.no/publications/2014/jaschke-escape24-crudeoil-preheat/full_paper160.pdf, full
    text skimmed.
  - The Mongstad application uses branches with at most 2 exchangers in series, so it is not affected by the error.
- J. Jäschke, S. Skogestad, *Optimal controlled variables for polynomial systems*, J. Process Control 22 (2012)
  167–179. This is the elimination method used in [JS14]. Not re-read.

## Literature / novelty check

The search was run on 2026-10-18 by a sub-agent. Items marked **[sub]** were read by the sub-agent and not re-read by
me. Items marked **[me]** I checked myself.

### Errata and corrections
- Crossref record `api.crossref.org/works/10.1016/j.compchemeng.2014.03.020`: empty `relation`, no `update-to` or
  `updated-by`, so no erratum or corrigendum is registered **[sub]**.
- Web searches ("Jäschke temperature", conjecture, proof, corrigendum) found nothing **[sub]**.

### Citing works
- Counts: OpenAlex 36, Semantic Scholar 38, Crossref 31. Almost all are off-topic: HEN synthesis, flexibility,
  district heating, extremum-seeking control, heat-current method. None proves, generalises or corrects the formula for
  \(M_j \ge 3\) **[sub]**.
- Relevant ones:
  - **J. Jäschke, S. Skogestad, ESCAPE-24 (2014)**, doi:10.1016/B978-0-444-63456-6.50102-2. It restates the same
    recursion, but branches have at most 2 exchangers **[me, sub]**.
  - **M. Trafczyński et al., Chem. Eng. Trans. 70 (2018)**, doi:10.3303/CET1870138. Jäschke-temperature controlled
    variables CV = JT_A − JT_B, JT_C − JT_D are used on a real crude preheat network with 26 exchangers in 4 branches,
    i.e. 6–7 per branch. The paper cites [JS14] without writing out the formula **[me: text checked; per-branch count
    from sub]**. This is a real application where the \(M_j \ge 3\) formula matters. Note that this network also has
    hot-stream interactions between branches, which are outside [JS14]'s assumptions anyway.
  - **Ikarimoto, MSc thesis, Univ. São Paulo (2019)**, doi:10.11606/D.3.2019.tde-19112019-102150. It reproduces the
    paper's \(M_1 = 4\) case with the paper's recursion (its eq. 48) and does not question it **[sub]**.
  - **Bernardino, Krishnamoorthy, Skogestad (PSE 2021+/ESCAPE-32, DYCOPS 2022)** and **Karlsen, NTNU project (2020)**:
    one exchanger per branch only **[sub]**.
  - **Jäschke, Cao, Kariwala, survey on self-optimizing control (Annu. Rev. Control 2017)**: lists HENs as an application
    only **[sub]**.
- Abstract only, all different approaches **[sub]**: Coşgun et al. ESCAPE-31 (2021); Sun, Zha, Luo, Energy (2018);
  Mohanan & Jogwar (2022–2023); Trafczyński et al., RSER (2019).

### Theses and patent not indexed as citing works
- **S. Aaltvedt, MSc thesis, NTNU (2013)**, skoge.folk.ntnu.no/diplom/diplom13/aaltvedt/ **[sub]**.
  - Uses the paper's previous-term recursion for 4 and 6 exchangers in series.
  - Reports an outlet temperature of 207.87 vs 207.84 °C, but a split of 64.15 % vs 70.66 %, explained as a "flat
    optimum". The comparison is confounded, because that plant model is not pure AMTD.
- **Aaltvedt, specialisation project (2012)**: uses the raw Maple expression for 3 exchangers (eq. 4.15), with no
  correction **[sub]**.
- **D. G. Edvardsen (2011) and A. Leruth (2012), NTNU MSc theses**: at most 2 exchangers in series **[sub]**.
- **Patent application WO2013164472A1** (PCT/EP2013/059304, GB1207770.7, "Parallel heat exchanger control"): formulas
  only for 1 exchanger and for 2 in series; no general formula **[sub]**.

### Conclusion of the novelty check
- No published proof of the \(M_j \ge 3\) case was found.
- No published notice of the error was found.
- No correct general formula was found in any notation.

Not checked:
- paywalled full texts: Coşgun 2021, Sun 2018, Trafczyński Appl. Sci. 2019 (HTTP 403);
- Google Scholar full-text search;
- a 2013 AIChE talk;
- the Palomar (`mmislan`) and formal-conjectures repositories named in PROTOCOL §1. These cover CRNT and graph theory,
  not heat-exchanger networks, so they are irrelevant here.

## Mathematical background used

None beyond the chain rule and linear algebra. The scripts use python-flint 0.9 (FLINT multivariate polynomials), sympy
and mpmath.
