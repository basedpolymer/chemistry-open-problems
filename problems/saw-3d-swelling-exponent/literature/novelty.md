# Novelty sweep 2024-2026 (and AI-driven campaigns)

Literature agent, searches run 2026-10-24. **An absence of results is not a proof that a problem is open.** The
openness of L3, U2, E1 and (T) rests on the verbatim statements in `open-statements.md` (latest by a leading author:
Slade 2019; latest survey: Guttmann-Jensen 2022, Hutchcroft ICM 2026) plus the absence of any later result in the
searches below.

## 1. What was searched

| Channel | Query / method | Result |
|---|---|---|
| arXiv API | all abstracts with "self-avoiding" in cat:math.PR since 2024-01-01 (full listing, ~40 entries) and in cat:math-ph not math.PR | No paper on displacement / nu of strict SAW on Z^d, d = 2, 3, 4. Relevant-adjacent: Lehner-Lindorfer-Panagiotis 2403.13121 (ballistic on graphs with > 1 end), Grimmett 2601.01128 (polygon counting on crystals), Liu-Slade 2605.15545 and Duminil-Copin-Markar-Panis-Slade 2605.21438 (SAW only for d >= 5), Duminil-Copin-Panis 2410.03649 (weakly SAW, d > 4), Lazar 2504.13093 (formal Fourier formulas for c_n, no exponent bounds), Khettabi-Velenik 2606.24352 (2D polymers in a force field), Qiu-Wang 2608.03497 (fractal graphs). |
| arXiv API | phrase searches since 2023: "self-avoiding walk" AND ("end-to-end" / displacement / ballistic / diffusive / "Flory" / "mean-square" / "swelling") | Only Krachun-Panagiotis 2310.17299 (hexagonal) and LLP 2403.13121 are mathematical results on displacement/ballisticity; the rest are physics/simulation papers. |
| arXiv API | "self-avoiding walk" AND "correlation length" (since 1990); "self-avoiding" AND bridge AND critical | No rigorous low-dimensional correlation-length exponent bound; only d >= 5 / weakly SAW d = 4 results and the 2D hexagonal works. |
| Semantic Scholar | works citing Madras 2014 (14 records), DCH 2013 (45, via arXiv id), KP (2), DGHM 2016 (18), Hammond TAMS 2019 (5), Hutchcroft 2018 (6), DGHM 2020 (9) | No new bound on E_n\|omega_n\|^2 for Z^d. Most recent citing items: Beffara (X-UPS survey), LLP 2024, Lazar 2025, Grimmett 2026, Qiu-Wang 2026. |
| OpenAlex | intended `cites:` queries | **Not done**: the shared free daily budget of this IP was exhausted (HTTP 429, "Insufficient budget ... resets at midnight UTC"). I did not create an API key (it would require sending the user's e-mail to a third party). Semantic Scholar was used instead. |
| Web search | "self-avoiding walk Z^3 mean-square end-to-end distance rigorous lower bound improvement Madras 4/(3d) 2025 2026"; "quantitative sub-ballisticity Z^d" | Nothing new; surfaced Hueter 2001 (see Sec. 3). |
| Zenodo API | quoted phrases: "self-avoiding walk(s)", "self-avoiding" AND ("Flory", "end-to-end", "diffusive", "mean-square displacement", "Z^3", "cubic lattice"), "sub-ballistic", "swelling exponent", "Duminil-Copin" AND "Hammond" | Nothing on nu. Recent SAW-related records concern connective-constant bounds on Z^2 (Tomar, Oct 2026; "Numaro", Jul 2026), enumeration of cubic-lattice contacts (Corti), and non-mathematical "framework" uploads. None claims a bound on E_n\|omega_n\|^2. |
| Evidence Press (ipitchford) | GitHub user listing is blocked for this session by the proxy (`users/...` endpoints, HTTP 403); used the GitHub MCP repository search `user:ipitchford`: 191 repositories listed, names and descriptions filtered for walk/avoid/polymer/Flory/swell/lattice/percol/connective/polygon; code search `"self-avoiding" user:ipitchford`; Zenodo "Evidence Press" (26 records) | Only `small-exponent-percolation` (long-range percolation) matched; nothing on SAW. Code search: 0 hits. |
| google-deepmind/formal-conjectures | GitHub code search `"self-avoiding" repo:google-deepmind/formal-conjectures`; file names of a shallow clone (2042 files) grepped for walk/avoid/polymer/flory/saw/percol/connective (file contents searched only through GitHub code search) | 0 hits; no SAW conjecture formalised there. |
| OpenAI math collection (github.com/openai/math) | `overview.tex` (catalogue dated 2026-10-06, "372 result families in 719 manuscripts"), grep for self-avoiding/honeycomb/walk/polymer/Z^3; code search `"self-avoiding" repo:openai/math` (60 files) | See Sec. 2: one family (entry 237), honeycomb lattice only. Nothing concerning Z^3. |

## 2. OpenAI catalogue, entry 237 (what exactly is claimed)

Source: `overview.tex` line 557 (fetched from raw.githubusercontent.com, saved in `certificates/sources/openai-math/`),
entry title "The three-quarter diameter exponent for honeycomb walks":

> "Proves the diameter form of Nienhuis's predicted three-quarter exponent: a uniformly chosen n-step self-avoiding
> walk on the honeycomb lattice has diameter n^{3/4+o(1)}. Its local mass and covering numbers have exponent 4/3.
> These estimates hold at every sufficiently large fixed length, simultaneously across scales, with arbitrarily high
> polynomial probability."

It links 13 preprints dated 2026-09-26, author "OpenAI". I read the abstract, Section 1 and Section 2 of "Mass and
covering exponents for fixed-length honeycomb walks". Its Theorem 1.1: for every delta, k > 0 and n >= n_0(delta, k),
with P_n-probability >= 1 - C n^{-k}: n^{3/4-delta} <= diam <= n^{3/4+delta}, plus local mass and covering bounds with
exponent 4/3; Corollary 1.4: E_n[D^p] = n^{3p/4+o(1)} and the same for the radius of gyration. The paper itself says:
"The assertions do not include convergence to a continuum curve, [...] or a fixed-length lower bound on |gamma_n - o|.
A lower bound for the diameter does not imply that the two endpoints are far apart." So, **if correct**, it would give
on H the upper bound E_n|gamma_n|^2 <= n^{3/2+o(1)} (endpoint inside the diameter) but no lower bound on the endpoint.
Its analytic inputs (Section 2, proved in companion papers) are integrable-lattice facts: a boundary-flux identity from
the parafermionic observable, B_h ≍ (1+h)^{-1/4} for the critical strip-crossing mass, polygon nesting masses with
exponent 1/12, a polygon length-square mass bound, and a one-arc bound on a cylinder with rhombus (Yang-Baxter)
parameters. `lean/docs/237.md` says the Lean formalisation covers finiteness of critical bridge masses, existence of a
free energy, and **B_h ≍ (1+h)^{-1/4}**; "The paper's h^{4/3+o(1)} mean-length law and its central-visit and
hexagonal-chord exponents are outside them", and the 3/4 spatial laws are not formalised.

Status: unrefereed, AI-generated; I did not verify any proof or run the Lean files. **Nothing in the catalogue concerns
Z^3 or Z^d**; the method depends on the honeycomb lattice's integrable structure. The owner's report is accurate in
saying the catalogue "revendique des résultats pour le réseau hexagonal" and that they do not resolve the 3D problem.
Other entries checked by keyword (215: planar O(n) models with n >= 3; 213, 509: percolation incl. Z^3) are not about SAW.

## 3. Old claims that solvers must not use

- I. Hueter, arXiv:math/0108120 (2001, v1 only, unpublished): claims nu(d) = max(1/4 + 1/d, 1/2) for SAW in Z^d, i.e.
  nu(3) = 7/12 ~ 0.5833, and arXiv:math/0108077 (2001) claims nu = 3/4 for the planar SAW (v2 comment: "The lower bound
  for the normalized first two moments of the distance of the weakly SAW from its starting point is not uniform in beta
  as beta tends to infinity"). Never published as far as I found; the later rigorous literature (BDGS 2012, Madras 2014
  "no non-trivial rigorous bounds in low dimensions", Slade 2019, Guttmann-Jensen 2022) treats the problem as open, and
  7/12 is excluded by the numerical value nu(3) = 0.58759700(40). Read: arXiv abstracts only. Treat as not established.

## 4. Conclusion of the sweep

As of 2026-10-24 I found no rigorous result improving, for strict SAW on Z^3: the lower bound E_N\|omega_N\|^2 >=
N^{4/9}/36 (Madras 2014), the qualitative upper bound o(N^2) (DCH 2013), the existence question for nu, or Kesten's
c_n <= mu^n exp(C n^{2/5} log n). The only quantitative sub-ballisticity theorem anywhere is on the hexagonal lattice
(Krachun-Panagiotis, AoP 2026); a stronger honeycomb claim exists only in the unrefereed OpenAI catalogue. Each rung
L1-L4, U1-U2, E1 and (T) is therefore, to the best of these searches, open in d = 3.
