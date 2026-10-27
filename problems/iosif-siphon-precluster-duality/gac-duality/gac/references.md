# Sources and novelty audit — agent G

GPT-6-Astra, effort xhigh. Checked 27 October 2026. This is a scoped search, not a proof of novelty. No claim here is advertised as new to the literature. All downloaded PDFs, extracted text, page images and web snapshots are in ignored `certificates/sources/`. No private correspondence is quoted.

## Primary mathematical sources actually read

| Source | Public link / identifier | Scope read and use |
|---|---|---|
| A. Iosif, *Duality in mass-action networks*, v2 (2026) | [arXiv:2603.08767v2](https://arxiv.org/abs/2603.08767v2) | Repository briefs and source material for the transpose definition and public small-codimension question; read result 17's README, manuscript definitions/proofs, extracted PDF and rendered p. 12. No audit of Iosif's entire manuscript is claimed. |
| D. F. Anderson, *A proof of the Global Attractor Conjecture in the single linkage class case* (2011) | [arXiv:1101.0761v6](https://arxiv.org/abs/1101.0761v6); DOI [10.1137/11082631X](https://doi.org/10.1137/11082631X) | Introduction, hypotheses of Theorem 4.10 and Corollary 4.11, relevant proof discussion. Corollary 4.11 is the GAC conclusion. |
| C. Pantea, *On the persistence and global stability of mass-action systems* (2012) | [arXiv:1103.0603v3](https://arxiv.org/abs/1103.0603v3) | Introduction, Theorem 5.1/Corollary 5.1 and Section 6, especially Theorem 6.3. Lemma 6.1 already uses the standard projection-rank/face-codimension identity; our new notation simply transposes it. |
| G. Craciun, F. Nazarov, C. Pantea, *Persistence and permanence of mass-action and power-law dynamical systems* (2013) | [arXiv:1010.3050v2](https://arxiv.org/abs/1010.3050v2); [author PDF](https://people.math.wisc.edu/~craciun/PAPERS/Craciun_Nazarov_Pantea_PERSISTENCE_PERMANENCE.pdf) | Introduction and hypotheses of Theorems 6.4 and 7.2. Distinguishes species dimension from stoichiometric rank. |
| M. Gopalkrishnan, E. Miller, A. Shiu, *A geometric approach to the Global Attractor Conjecture* (2014) | [arXiv:1305.5303v2](https://arxiv.org/abs/1305.5303v2); DOI [10.1137/130928170](https://doi.org/10.1137/130928170) | Introduction, strongly endotactic definitions, Theorems 1.1--1.2 and main-result scope. The permanence theorem is the external dynamical input in Proposition 4.1. Birch/interior-equilibrium discussion also consulted. |
| D. F. Anderson, A. Shiu, *The dynamics of weakly reversible population processes near facets* (2010) | [author PDF](https://web.math.wisc.edu/~dfanderson/papers/AndShiuFacetsSIAM.pdf); DOI [10.1137/090764098](https://doi.org/10.1137/090764098) | Introduction, Theorems 3.2, 3.4, 4.6 and Corollary 4.7. The GAC criterion includes vertices as well as facets and empty faces. |
| C. Xu, *First order endotactic reaction networks* (arXiv v2 2025; journal 2026) | [arXiv:2409.01598v2](https://arxiv.org/abs/2409.01598v2); DOI [10.1007/s00285-026-02378-2](https://doi.org/10.1007/s00285-026-02378-2) | Full PDF available; read introduction, source-order definition and Theorem 6.5/its statement context. Journal DOI is supplied in the brief and L's source audit; direct publisher retrieval failed, so the theorem checked here is the arXiv v2 theorem. No missing PDF is needed for this task. |
| C. Wiuf, *A proof of the Global Attractor Conjecture in a special case* (2026) | [arXiv:2609.24553v1](https://arxiv.org/html/2609.24553v1) | Sections 2--3, 6--7, condition (diamond), Theorem 7.2, examples and relevant proof. Rendered PDF pp. 4 and 14 inspected; equations checked against HTML. Our rank and Stiemke translations are proved independently. The full dynamical theorem is cited, not independently verified. |
| D. Angeli, P. De Leenheer, E. D. Sontag, *A Petri net approach to the study of persistence in chemical reaction networks* (2007) | [author PDF](https://www.sontaglab.org/FTPDIR/angeli_leenheer_sontag_math_biosciences_MBS-D-06-00188R1.pdf); DOI [10.1016/j.mbs.2007.07.003](https://doi.org/10.1016/j.mbs.2007.07.003) | Definitions, Theorem 2 and boundary-siphon proof in Section 5. The elementary persistence proof is re-derived in our Proposition 6.1. |
| F. Blanchini, G. Giordano, *Dual chemical reaction networks and implications for Lyapunov-based structural stability* (2022) | DOI [10.1109/LCSYS.2021.3081369](https://doi.org/10.1109/LCSYS.2021.3081369) | Supplied IEEE accepted version, pp. 1--6, especially assumptions, Theorems 1--3 and Remark 2. PDF pp. 2 and 5 rendered and inspected. Its Example 5 already states that WR need not be preserved: that observation is not novel here. |
| G. Craciun, *Toric Differential Inclusions and a Proof of the Global Attractor Conjecture*, v2/v3 | [v2](https://arxiv.org/abs/1501.02860v2), [v3](https://arxiv.org/html/1501.02860v3) | v3 abstract, introduction, strategy, Theorem 3.1 and relevant statements/construction in Sections 7--8; v2 Section 4.5 and start of Section 5. Compared the general-dimensional argument's organization. No full proof audit. |

## Recent general results and the precise formalization scope inspected

- [OpenAI, boundedness/persistence paper, 25 September 2026](https://github.com/openai/math/blob/main/preprints/Boundedness-and-persistence-of-weakly-reversible-mass-action-systems-September-25-2026/paper.pdf): Theorem 1.1 and Remark 3.5 read. The theorem is stronger than GAC's persistence requirement, since it applies to all WR networks with constant positive rates. No proof audit here.
- [OpenAI, uniform permanence paper, 5 October 2026](https://github.com/openai/math/blob/main/preprints/Uniform-Permanence-in-Weakly-Reversible-Mass-Action-Systems-October-5-2026/permanence.pdf): abstract, first two pages and Theorem 1.1 read. This asserts bounds eventually uniform within each positive class.
- [Formalization scope, family 149](https://github.com/openai/math/blob/main/lean/docs/149.md), [top-level implementation](https://github.com/openai/math/blob/main/lean/OAI/Analysis/MassAction/Main.lean), and [Comparator target](https://github.com/openai/math/blob/main/lean/ComparatorChallenges/MassAction.lean) read. The target ends with `sorry`, as a challenge statement; the implementation contains a proof via imported modules. The inspected scope covers boundedness/persistence, not by itself the stronger class-uniform permanence result. We did not compile Lean, inspect all dependencies or run an axiom audit.
- [Palomar endotactic permanence](https://github.com/mmislan/palomar-endotactic-permanence): README and declared scope read. It claims the broader endotactic permanence theorem with measurable bounded positive rates and reports Lean checks. No independent audit here.

These available general claims and proof artifacts preclude marketing the restricted certificate of Proposition 4.1 as a new GAC case. It is a structural route to an established strongly endotactic theorem and remains useful independently of accepting a recent general proof.

## Search log and novelty limits

Executed the required read-only repository searches:

```text
gh api "users/mmislan/repos?per_page=100"
gh api "users/ipitchford/repos?per_page=100" --paginate
gh search repos '"dual reaction network"'
gh search repos '"reverse-dual"'
gh search repos '"transpose stoichiometric matrix"'
gh search repos '"global attractor conjecture"'
gh search repos '"preclusters"'
gh search repos '"siphon trap duality"'
gh search code '"global attractor" repo:google-deepmind/formal-conjectures'
gh search code '"preclusters" repo:google-deepmind/formal-conjectures'
```

The mmislan listing found the directly relevant permanence repositories. The ipitchford metadata scan found no transpose/precluster/siphon work; reaction-related hits included reverse-reaction deletion and oscillator projects with different questions. The repository phrase searches found no directly relevant transpose-persistence theorem; the dual-reaction search returned neural-network and unrelated software projects. Both formal-conjectures code searches returned no hits. Metadata/code search is not an audit of every repository file.

Zenodo API queries used the phrases above. The combined duality/precluster/siphon query returned one unrelated engineering record; separate exact-phrase queries for `transpose stoichiometric matrix` and `siphon trap duality` returned zero hits. `global attractor conjecture` returned one record titled *Tropical Geometry and Biochemical Reaction Networks: A Mathematical Framework for Steady-State Topology*; only its metadata was inspected, and no mathematical claim here relies on it.

Web searches combined the duality phrases with persistence, GAC, Gale duals, and codimension. They recovered the expected transpose/flux-concentration and Petri-net literature, with additional prior-art context in agent L's notes. Later targeted searches for `strongly endotactic cyclic prime`, `endotactic circulant`, and `reaction network coordinate permutation duality` did not locate the exact prime-species criterion. This is insufficient to establish novelty: it is an elementary combination of row matching, cyclotomic irreducibility and the known strongly endotactic theorem. The saddle and extinction examples likewise require an independent novelty audit before any originality claim.

## Internal read-only context and credits

Read `../literature/notes.md` and `../petri/candidates.md` as they became available. Agent L supplied useful source leads and independently noted the Jacobian-factor transpose identity; it was re-derived here. Agent P's WR equal-conflict observation generalizes the disjoint-source argument behind our first-order-dual reduction. Our claims do not rely on unproved teammate candidates. Requests for E are recorded only in `requests-for-E.md`; E's independent checks, if received, are distinguished from G's own exact certificates.
