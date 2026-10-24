# References and novelty check

## Sources read

| # | Reference | Identifier | What was read |
|---|---|---|---|
| 1 | K. Hake, *Knotting probability of equilateral hexagons*, J. Knot Theory Ramifications 28(7) (2019) 1950040 | doi:10.1142/S0218216519500408, arXiv:1810.11882 | arXiv v1 (28 Oct 2018), full text. Theorem 4.1: "The probability that an equilateral hexagon is knotted is at most (14−3π)/192." Measure: T135 action-angle coordinates, Lebesgue × Haar (Duistermaat–Heckman, her Theorem 2.22). Monte Carlo (§4): fraction of Joint Chirality-Curl (1,1) trefoils 3.426005·10⁻⁵ (s.d. 2.24·10⁻⁶, runs of 10⁷), times four types = 1.370402·10⁻⁴. Remark: the introduction's "(14−3π)/192 < 1/42" is false ((14−3π)/192 = 0.0238293 > 1/42 = 0.0238095). Published version not seen (NEEDS.md). |
| 2 | J. Cantarella, C. Shonkwiler, *The symplectic geometry of closed equilateral random walks in 3-space*, Ann. Appl. Probab. 26(1) (2016) 549–596 | doi:10.1214/15-AAP1100, arXiv:1310.5924 | arXiv v2: abstract, §2 (Theorem 1, Duistermaat–Heckman), §3–4 (standard measure, Theorem 13 Kapovich–Millson / Howard–Manon–Millson / Hitchin), §4.1 (Proposition 16: at least 1/2 of hexagons unknotted; PTSMCMC estimate 1.3·10⁻⁴ ± 0.2·10⁻⁴ from 1 318 001 samples, 173 knots), Lemma 31 (S_n invariance). **Correction to the brief:** this Monte Carlo value is Cantarella–Shonkwiler's, not Xiong et al.'s. |
| 3 | J. A. Calvo, *Geometric knot spaces and polygonal isotopy*, J. Knot Theory Ramifications 10(2) (2001) 245–267 | doi:10.1142/S0218216501000834, arXiv:math/9904037 | arXiv v2, §§1–2: Theorem 1 (Geo(3..5) connected, unknots only; Geo(6) has five components), statement of the Δ₂Δ₄Δ₆ criterion "(i) right-handed trefoil iff Δ2 = Δ4 = Δ6 = 1, (ii) left-handed iff all −1, (iii) unknot iff Δi = 0 for some i", with the remark that if both candidate edges pierce a triangle their contributions cancel; Theorem 2 (joint chirality-curl). |
| 4 | J. Geisz, C. Peterson, C. Shonkwiler, *Trefoil probabilities and polyhedra*, arXiv:2609.24860 (21 Sep 2026) | arXiv:2609.24860 | Full text (pp. 1–10). Exact trefoil probability 1/(16π²) for six uniform points **on the sphere** (a different, random-jump model). Introduction: "we are not aware of any exact results about knot probabilities for any finite n in such models (though some bounds for related models are known; e.g., [16, 35])", [16] = Hake. §5: "Returning to the original question about knot probabilities, it would be interesting to determine the exact probability of knotting, even for hexagons, in any other random knot model with a continuous state space." Lemma 2.2 restates Calvo's criterion; Lemma 2.3 (Huh–Jeon) the "only e_{i+2}, e_{i+3}" observation. |
| 5 | A. Xiong, A. J. Taylor, M. R. Dennis, S. G. Whittington, *Knot probabilities in equilateral random polygons*, J. Phys. A 54 (2021) 405001 | doi:10.1088/1751-8121/ac1fc2, arXiv:2108.07197 | arXiv v1: introduction and §2 (action-angle sampler of Cantarella et al., Alexander-polynomial invariants, N from 6 to 4000). No rigorous bound and no hexagon-specific value stated in the text read. |

## Cited, not read (classical or paywalled)

| Reference | Identifier | Used for |
|---|---|---|
| J. A. Calvo, *The embedding space of hexagonal knots*, Topology Appl. 112 (2001) 137–174 | doi:10.1016/S0166-8641(99)00229-1 | Lemma 8 (proof of the Δ criterion), only in Remark 3.4 (not used in the proof). NEEDS.md. |
| J. A. Calvo, PhD thesis, UC Santa Barbara (1998) | ProQuest 304416529 | background |
| R. H. Crowell, R. H. Fox, *Introduction to Knot Theory*, Ginn 1963 (Springer GTM 57) | — | F1 (Δ-moves, polygonal knots), F3 (colourings) |
| G. Burde, H. Zieschang, *Knots*, de Gruyter 1985 | — | F1 |
| J. H. Przytycki, *3-coloring and other elementary invariants of knots*, Banach Center Publ. 42 (1998) 275–295 | arXiv:math/0608172 | F3 |
| J. J. Duistermaat, G. J. Heckman, Invent. Math. 69 (1982) 259–268; M. Kapovich, J. J. Millson, J. Differential Geom. 44 (1996) 479–513; B. Howard, C. Manon, J. Millson, Canad. J. Math. 63 (2011) 878–937 | — | measure-preserving action-angle coordinates, via [2], Theorems 1 and 13 |

## Novelty check (24 October 2026)

* Crossref: Hake's JKTR paper has 1 citing work (is-referenced-by-count = 1), which is Geisz–Peterson–Shonkwiler [4]
  (reference [16] there). OpenAlex and Semantic Scholar APIs were rate-limited from this container (HTTP 429 /
  budget exhausted), so the citing-works check relied on Crossref and web search.
* Web searches ("knotting probability equilateral hexagons", "equilateral hexagon trefoil probability symplectic",
  with 2025/2026) return only [1], [2], [4], Hake's 2019 AMS talk slides, and stick-number papers
  (Eddy–Shonkwiler, arXiv:1909.00917, confined polygons — no bound for unconfined hexagons).
* Geisz–Peterson–Shonkwiler (Sept. 2026) explicitly list the equilateral case as known only through bounds and
  call exact knotting probabilities "even for hexagons" in other continuous models an open question.
* Zenodo search ("equilateral hexagon knot", "hexagonal trefoil probability", "knotting probability hexagon"):
  nothing relevant (one unrelated record on stick-filtered polygonal knot spaces, Ozawa, 10.5281/zenodo.23460843).
* AI campaigns: the GitHub API user listings were blocked (HTTP 403, session restricted); the web pages
  github.com/ipitchford and github.com/mmislan filtered on "knot", "hexagon", "polygon" show no repositories.
  google-deepmind/formal-conjectures and openai/math were not searchable from here (repository-scoped API only);
  their topics (graph theory, number theory) make a hit unlikely.

**Conclusion:** no improvement of Hake's bound and no positive lower bound were found; the problem is open in
the sense of the verbatim sentence of [4], §5.
