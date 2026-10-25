# Sources and novelty audit

Audit date: **25 October 2026**. Research and assessment: GPT-6-Astra (OpenAI), effort xhigh; no sub-agents.

**Conclusion:** no resolution of the full Brownian finite-mass profile uniqueness problem was located. This is a documented literature-search conclusion, not a guarantee that unpublished or unindexed work does not exist.

## Primary mathematical sources

| Source | Material inspected | Relevance and scope |
|---|---|---|
| Escobedo–Mischler–Rodríguez Ricard (2005), *On self-similarity and stationary problem for fragmentation and coagulation models*, [DOI 10.1016/j.anihpc.2004.06.001](https://doi.org/10.1016/j.anihpc.2004.06.001) | Publisher abstract and bibliographic record; full proof not read during selection | Foundational finite-mass similarity existence result, also cited by the later uniqueness papers. |
| Niethammer–Throm–Velázquez, [arXiv:1510.03361v3](https://arxiv.org/abs/1510.03361v3), revised 26 June 2017 | PDF pp. 1–4: physical kernel, open statement, profile equation, assumptions, Theorem 1.3. Printed p. 2 also rendered and visually checked. The full 99-page proof has not been audited. | Uniqueness for sufficiently small perturbations; corrects an earlier proof gap. The old 2014 result must not be used without its corrections. |
| Throm (2021), *Stability and Uniqueness of Self-similar Profiles in L1 Spaces for Perturbations of the Constant Kernel in Smoluchowski's Coagulation Equation*, [DOI 10.1007/s00220-021-03967-6](https://doi.org/10.1007/s00220-021-03967-6); [arXiv:1902.10000v1](https://arxiv.org/abs/1902.10000v1) | Preprint PDF §§1–2 and beginning of §3, printed pp. 1–7; especially Theorem 2.5 | Removes analyticity restrictions and treats larger singularity exponents. The amplitude must still be sufficiently small. Full proof not independently checked. |
| Cañizo–Throm (2021), *The scaling hypothesis for Smoluchowski's coagulation equation with bounded perturbations of the constant kernel*, [DOI 10.1016/j.jde.2020.07.036](https://doi.org/10.1016/j.jde.2020.07.036); [arXiv:1910.08015v1](https://arxiv.org/abs/1910.08015v1) | PDF introduction, kernel assumptions and Theorem 1.2, pp. 1–4 | Quantitative attraction for small bounded perturbations. Brownian ratio terms are unbounded; this theorem does not cover them. |
| Ferreira–Lukkarinen–Nota–Velázquez (2024), *Asymptotic localization in multicomponent mass conserving coagulation equations*, [DOI 10.2140/paa.2024.6.731](https://doi.org/10.2140/paa.2024.6.731); [arXiv:2203.08076v1](https://arxiv.org/abs/2203.08076v1) | Abstract and HTML §1.3, Theorem 1.6 and equation (1.25) | Its uniqueness conclusion assumes K(rθ,sθ)=Q(θ), independent of r,s. The scalar Brownian kernel depends on r/s, so it fails this hypothesis. Homogeneity under simultaneous scaling is weaker than (1.25). |

## Recent literature and discovery channels

- [Cañizo's publication list](https://canizo.org/publications/) and [Throm's university page](https://www.umu.se/en/staff/sebastian-throm/) were inspected. Their recent 2025 similarity results concern the **inelastic Boltzmann equation**, a different collision operator. Niethammer's Bonn list was found by web search, but direct retrieval timed out; it is not counted as a complete author-publication audit.
- OpenAlex work `W2219054489` (revised uniqueness preprint): the query for citing works published from 2024 returned **0**.
- OpenAlex work `W2917077656` (Throm's CMP paper): the same query returned **1**, the 2024 multicomponent localization paper checked above.
- OpenAlex work `W2980552127` (Cañizo–Throm's JDE paper): **4** citing works from 2024 onward. One is the localization paper; two concern inelastic Boltzmann profiles/relaxation, and one concerns first-passage times in a pure-birth coalescence model. These last three were screened by metadata and available author abstracts, not full-proof review. The raw response is cached locally.
- Semantic Scholar's API returned **HTTP 429**; no successful citation audit is claimed for that service.
- Web searches covered Brownian kernel, Smoluchowski, profile uniqueness, self-similarity, 2025/2026, and computer-assisted proofs. Initial-value uniqueness, linear coagulation and moment-closure approximations were distinguished from the nonlinear similarity-profile question.

## Public problem and AI-assisted research repositories

- [MColbrook/AIM problem 163](https://github.com/MColbrook/AIM/blob/8729eab98bafe409b66c5e39253bf56569a62327/problems/163-brownian-coagulation-profile.md) remains open in the inspected tree. Its own last-check date is 8 September 2026. The tree and pull-request titles were searched; no Brownian profile solution was located. Snapshot commit: `8729eab98bafe409b66c5e39253bf56569a62327`.
- GitHub repository names and descriptions for **mmislan/Palomar** and **ipitchford/Evidence Press** were inspected through the API, including pagination for Evidence Press. No matching claimed resolution was found. This was a metadata search, not a line-by-line audit of every repository.
- GitHub repository searches `Brownian coagulation` and `Smoluchowski uniqueness` returned no repositories. [Lixiviant's publication list](https://lixiviant.ai/publications/) was also consulted during candidate screening; no matching claim was located.
- Zenodo query `"Brownian" AND "coagulation"` returned **7** records, concerning application data/software and unrelated subjects; no profile uniqueness proof. Broader follow-up queries combining coagulation and uniqueness returned HTTP 400, so broader Zenodo coverage is incomplete.

## Cached source integrity

Files are in `certificates/sources/`, intentionally ignored by Git. SHA-256:

```text
78f4f783b60810ce25549f2e544db88b02d9205470708e745f28ee09bfb61d05  niethammer-throm-velazquez-1510.03361v3.pdf
df174a17aaa702e285a56a2cea6bd7ab5b0cdd7cf64e1d6fa2f7c3c43e6499d8  throm-1902.10000v1.pdf
f93a15440ead6a6bc293f14197f386bbf1e0af2a06df07c80016cb4d447a3198  canizo-throm-1910.08015v1.pdf
```
