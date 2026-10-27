# Sources and novelty disposition

Audit date: 27 October 2026. Attribution: GPT-6-Astra (effort High), as specified by the user.

## Decisive primary source

Ivar J. Halvorsen and Sigurd Skogestad, **Minimum Energy for the Four-Product Kaibel Column**,
AIChE Annual Meeting 2006, paper 216d, 8 pages.

- [Author-hosted full text](https://skoge.folk.ntnu.no/publications/2006/halvorsen_kaibel_aiche06/paper216d.pdf).
- [Conference-proceedings copy](https://skoge.folk.ntnu.no/prost/proceedings/aiche-2006/data/papers/P67832.pdf).
- [AIChE record](https://proceedings.aiche.org/conferences/aiche-annual-meeting/2006/proceeding/paper/216d-minimum-energy-four-product-kaibel-distillation-column).
- DOI: none identified for this conference contribution; no DOI is invented.

**Read:** full text of the 8-page paper and conference-record abstract. Section 2, PDF pages 2–4,
defines the model and gives the algebraic calculation. Section 4, PDF page 5, gives the strict
comparison of the peaks. Table 1, PDF page 4, supplies the numerical example reproduced here.
PDF page 3 was also rendered and inspected visually to verify the signs and the printed `min`
in equation (7). The correction to `max` follows from simultaneously imposing both section
requirements and is consistent with the source's text and table.

SHA-256 of the downloaded author-hosted PDF:

```text
a986ad8e4dad83d695484b82e6fd322b37bbac7257da18716b72d48c2a09c42b
```

Local copy: `certificates/sources/halvorsen-skogestad-2006.pdf`; rendered equation page:
`certificates/sources/equations.png`. These third-party artifacts are deliberately Git-ignored.

## Disposition of the candidate

The repository's candidate was a provisional question, explicitly marked as needing a literature
check. It was not backed by a quotation declaring it open. The 2006 primary source already
answers both the strict-comparison and analytic-expression parts in its ideal model. Under
PROTOCOL.md §1.3 this ends the novelty search: record the prior answer and stop treating the
candidate as an open research problem. A broad search for later solver campaigns is unnecessary
to establish that the question was already answered in 2006.

The present work supplies an explicit algebraic exposition and reproducibility checks. The
eliminated-root form of the formula is obtained directly from the source's two quadratic
equations; no novelty claim is made for that rewriting. Generalization outside the ideal model
would be a different question and is not inferred from this audit.

The 2001 thesis and 2003 generalized-Petlyuk article cited by the source were not freshly audited
in this task. The note attributes and uses the Underwood model as presented in the 2006 primary
source; it does not claim a new stage-by-stage proof or rely on the repository's separate result 12.
