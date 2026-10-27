#!/usr/bin/env python3
"""R1: apply R2's replacements M01-M10 to references.md (review/R2/review.md §8.3), plus two follow-ups of the merge
(the Blanchini-Giordano page in the novelty list; the archived report's numbering is not in references.md).
Each old text must occur exactly once; otherwise nothing is written.  Standard library only, < 1 s.
Usage (from this folder): py -3 scripts/apply_refs_md.py [--dry-run]"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PATH = HERE.parents[2] / "references.md"
TEXT = PATH.read_text(encoding="utf-8")

M = [
    ("M01", """- S. Müller, G. Regensburger, SIAM J. Appl. Math. 72 (2012) 1926–1947, DOI 10.1137/110847056 (sign-vector criterion,
  Proposition 3.9). Cited after agents Q3, L2; title and pages matched through Crossref on 28 October 2026; not re-read.""",
     """- S. Müller, G. Regensburger, SIAM J. Appl. Math. 72 (2012) 1926–1947, DOI 10.1137/110847056, arXiv:1209.6488
  (sign-vector criterion: Proposition 3.1 of arXiv v1, read by reviewer R2 on 28 October 2026 and checked again by
  reviewer R1 from the arXiv PDF; the numbering of the journal version was not checked, so the paper cites the arXiv
  version). Cited after agents Q3, L2; title and pages matched through Crossref."""),
    ("M02", "Remark 2, p. 6 (the dual of a reversible network is not",
     "Remark 2, p. 493 (the sixth page; the dual of a reversible network is not"),
    ("M03", "of Theorem 37; their DOIs were not verified and are not given.",
     """of Theorem 37. DOIs 10.1109/TCS.1978.1084490 (full title: *Asymptotic behavior of nonlinear compartmental systems:
  nonoscillation and stability*) and 10.1109/TAC.2004.841888, matched through Crossref by reviewer R2 on 28 October 2026."""),
    ("M04", "Control (2023), DOI 10.1109/TAC.2023.3275880.",
     "Control 68(12) (2023) 7420–7435, DOI 10.1109/TAC.2023.3275880."),
    ("M05", "IEEE Trans. Automat. Control 36 (1991) 1368–1382. Author's copy",
     "IEEE Trans. Automat. Control 36 (1991) 1368–1382, DOI 10.1109/9.106153 (Crossref). Author's copy"),
    ("M06", "Congress (2002). Abstract and §1–3 read (L2).",
     "Congress (2002), IFAC Proc. Vol. 35(1) 455–460, DOI 10.3182/20020721-6-ES-1901.00078 (Crossref). Abstract and §1–3 read (L2)."),
    ("M07", "J. Appl. Math. 73 (2013) 305–329, arXiv:1010.3050.",
     "J. Appl. Math. 73 (2013) 305–329, DOI 10.1137/100812355 (Crossref), arXiv:1010.3050."),
    ("M08", "0-521-46519-2 (author-hosted PDF). Read: Def. 4.4,",
     "0-521-46519-2, DOI 10.1017/CBO9780511526558 (Crossref) (author-hosted PDF). Read: Def. 4.4, Defs 5.1, 5.11,"),
    ("M09", "  weakly reversible; mass action is not preserved).",
     """  weakly reversible; mass action is not preserved).
- I. Famili, B. Ø. Palsson, *The convex basis of the left null space of the stoichiometric matrix leads to the
  definition of metabolically meaningful pools*, Biophys. J. 85 (2003) 16–26, DOI 10.1016/S0006-3495(03)74450-6
  (open access, PMC1303061). Compound-map passages read by agent L (`gac-duality/literature/notes.md` §10(a)): the
  compound map of −Sᵀ, with the extreme conservation laws drawn as extreme pathways.
- N. Jamshidi, B. Ø. Palsson, *Flux-concentration duality in dynamic nonequilibrium biological networks*, Biophys. J.
  97 (2009) L11–L13, DOI 10.1016/j.bpj.2009.06.049 (open access, PMC2749767). Abstract and main text read by agent L
  (§10(b)): the Jacobians S·G and G·S.
- A. Duvall, *Computationally checking if a reaction network is monotone or non-expansive*, arXiv:2409.12487 v2 (27
  January 2025). §11 ("Duality") read by agent L (§10(e)): the dual of an irreversible network has stoichiometric
  matrix Γᵀ; monotonicity and non-expansivity transfer."""),
    ("M10", "Google Scholar was not reachable. Prior art found and credited in the paper: Desel–Esparza Thm 4.30(a), Prop. 5.4,",
     "Google Scholar was not reachable. Prior art found and credited in the paper: Famili–Palsson 2003, Jamshidi–Palsson 2009 and Duvall 2024 (networks with stoichiometric matrix ±Γᵀ); Teruel–Silva p. 293 and Cor. 32; Desel–Esparza Thm 4.30(a), Prop. 5.4,"),
    ("R1-M11", "Zhang–Gao–Dochain Theorem 3.10; Blanchini–Giordano p. 6;",
     "Zhang–Gao–Dochain Theorem 3.10; Blanchini–Giordano p. 493;"),
]

text, bad = TEXT, 0
for i, old, new in M:
    c = text.count(old)
    if c != 1:
        print(f"FAIL {i}: old text occurs {c} times")
        bad += 1
        continue
    text = text.replace(old, new)
    print(f"PASS {i}")
print(f"{len(M)} replacements, problems: {bad}")
if not bad and "--dry-run" not in sys.argv:
    PATH.write_bytes(text.encode("utf-8"))
    print("written: references.md")
sys.exit(1 if bad else 0)
