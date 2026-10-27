r"""R2: apply the corrections of review.md to a COPY of paper.tex / references.bib, to show that every replacement
is mechanically applicable (each old text occurs exactly once) and that the result compiles.

Environment: Windows `py -3` (standard library). Runtime: < 1 s (plus the LaTeX build, about 1 min).
Usage, from the repository root:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/apply_corrections.py
  py -3 templates/build_tex.py problems/iosif-siphon-precluster-duality/gac-duality/review/R2/build/paper.tex
It reads the problem folder's paper.tex and references.bib (never writes them) and writes build/paper.tex and
build/references.bib inside the R2 folder. The ids C01..C43 and B01..B08 are those of review.md.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent          # gac-duality/review/R2
FOLDER = HERE.parents[2]                                 # problem folder
TEX = (FOLDER / "paper.tex").read_text(encoding="utf-8")
BIB = (FOLDER / "references.bib").read_text(encoding="utf-8")

R = []  # (id, old, new)


def rep(i, old, new):
    R.append((i, old, new))


# ---------- abstract and introduction ----------
rep("C01", r"""theory (Hack 1972). Through it, siphons of \(\Net\) are the T-traps of \(\Nd\), the seeds of preclusters are the
Petri clusters, and minimal siphons cross every cluster at most once, so preclusters and siphons are related by
transversality, not refinement.""",
    r"""theory (Hack 1972). Through it, siphons of \(\Net\) are the T-traps of \(\Nd\), and for equal-conflict networks the
seeds of preclusters are the Petri clusters and minimal siphons cross every cluster at most once, so preclusters and
siphons are related by transversality, not refinement.""")
rep("C02", r"""(respectively conservative) network, with equality exactly when twin species are coupled, in particular on the""",
    r"""(respectively conservative) network, with equality exactly when doubling pairs (respectively twin species) are
coupled, in particular on the""")
rep("C03", r"""preprints. Every example and finite claim is certified in exact arithmetic (129 claims); the main results were
re-derived by three independent verifiers.""",
    r"""preprints, which we have not checked. Every example is certified in exact arithmetic (129 claims); counts over large
families are cited from the working notes as checked by computation. The main results were re-derived by three
fresh-context verifiers.""")
rep("C04", r"""Petri clusters when \(\Net\) is equal-conflict.""",
    r"""Petri clusters when \(\Net\) is equal-conflict, every species is consumed and no source is zero (e.g.\ under (H1)
and (H2)).""")
rep("C05", r"""\emph{Verification}: the independent verifier that re-derived""",
    r"""\emph{Verification}: the fresh-context verifier that re-derived""")
# ---------- Table 1 ----------
rep("C06", r"""construction known (Hack; Desel--Esparza; Teruel--Silva); identification not found \\""",
    r"""construction known (Hack; Desel--Esparza; Teruel--Silva); \(\Gamma\tp\) networks known (Famili--Palsson; Blanchini--Giordano; Duvall); identification not found \\""")
rep("C07", r"""equality iff twins coupled & P. & V1 & not found \\""",
    r"""equality iff doubling (twin) pairs coupled & P. & V1 (corr.) & not found \\""")
rep("C08", r"""\ref{thm:WFFC} & structure of well-formed free-choice networks & P. & V1 &""",
    r"""\ref{thm:WFFC} & structure of well-formed free-choice networks & P. (sketch) & V1 &""")
rep("C09", r"""extends Pérez Millán et al.\ 2012, Campos et al.\ 1991 \\""",
    r"""extends Pérez Millán et al.\ 2012, Campos et al.\ 1991; cf.\ Recalde et al.\ 2002, Júlvez et al.\ 2005 \\""")
rep("C10", r"""Hurwitz at every point & P. & V2 & not found \\""",
    r"""Hurwitz at every point & P. & V2 & not found (tools classical) \\""")
rep("C11", r"""does not give uniqueness (NC3) & P. & V2 & not found \\""",
    r"""does not give uniqueness (NC3) & P. & found by V2; certificate & not found \\""")
rep("C12", r"""sufficient half weaker in Zhang--Gao--Dochain 2023 \\""",
    r"""``if'' half of (1) under a stronger hypothesis: Zhang--Gao--Dochain 2023; rest not found \\""")
rep("C13", r"""\ref{thm:A} & a persistent dual with an attracting-looking facet (refutes Conj.~F) & P. & V3 & not found \\""",
    r"""\ref{thm:A}, \ref{thm:AF} & a persistent dual with an OUT facet (refutes a conjecture of agents Q1 and D2); the family \(N_{u,R}\) & P.; P. (sketch) & V3 & not found \\""")
rep("C14", r"""duality keeps \(s\); prime-cycle certificate & P. & V3 &""",
    r"""duality keeps \(s\); prime-cycle certificate & P.; P. (sketch) & V3 &""")
# ---------- Section 2 ----------
rep("C15", r"""\(\Net\) gives a loop \(y\to y\) of \(\Nd\), and a species never consumed gives the zero complex; all definitions
apply to such networks.""",
    r"""\(\Net\) gives a loop \(y\to y\) of \(\Nd\), and a species never consumed gives the zero complex; all definitions
apply to such networks. Networks with stoichiometric matrix \(\pm\Gamma\tp\) appeared before: Famili and Palsson's
compound map of \(-\Gamma\tp\), on which the extreme conservation laws are drawn as extreme pathways
\autocite{FamiliPalsson2003}; the dual network of Blanchini and Giordano~\autocite{BlanchiniGiordano2022}, developed at
the same time as Iosif's (source, p.~8); and Duvall's dual of an irreversible network~\autocite{Duvall2024}. Jamshidi
and Palsson compare the Jacobians \(\Gamma G\) and \(G\Gamma\)~\autocite{JamshidiPalsson2009}. These works transfer
linear and local properties (conservation pools, spectra, structural stability, monotonicity), not siphons or
persistence.""")
rep("C16", r"""\begin{lemma}[transposition {\autocite[Lemma~1]{ArchivedReport}}]""",
    r"""\begin{lemma}[transposition {\autocite[Lemma~2]{ArchivedReport}}]""")
rep("C17", r"""\(\Ye[S_c,T_c]=a_cJ\) is symmetric \\""",
    r"""\(\Ye[S_c,T_c]=a_cJ\) is symmetric; Teruel--Silva Cor.~32 \\""")
rep("C18", r"""MIPS (\zcref{prop:mips}) & classes of fully coupled reactions & \autocite[Prop.~9]{ArchivedReport} \\""",
    r"""MIPS (\zcref{prop:mips}) & classes of fully coupled reactions & \autocite[Prop.~6]{ArchivedReport} \\""")
rep("C19", r"""\(K_1\to K_2\), \(K_2\to K_1\) is not.""",
    r"""\(K_1\to K_2\), \(K_2\to K_1\) is not. Teruel and Silva~\autocite[p.~293]{TeruelSilva1996} make the same
observation for weighted nets: equal-conflict nets with balanced joins form a self-reverse-dual class, and their Duality
Theorem~\autocite[Cor.~32]{TeruelSilva1996} states that a pure net is a well-formed balanced-joins EC net iff its
reverse dual is.""")
# ---------- Section 3 ----------
rep("C20", r"""\begin{proposition}[{\autocite[Prop.~3, Prop.~9]{ArchivedReport}}]\label{prop:mips}""",
    r"""\begin{proposition}[{\autocite[Prop.~3, Prop.~6]{ArchivedReport}}]\label{prop:mips}""")
rep("C21", r"""the rate-free form with \(\Ph\) is
new as a statement.""",
    r"""the rate-free form with \(\Ph\) was
not found in the literature as a statement (\zcref{sec:search}).""")
rep("C22", r"""\item For every consistent network, \(\Mips(\Nd)\refi\Ph(\Net)\); and \(\Ph(\Net)=\Mips(\Nd)\) iff \(\Ph(\Net)\refi
\Mips(\Nd)\) iff both have the same number of blocks iff every doubling pair is fully coupled.""",
    r"""\item If \(\Net\) is consistent, then \(\Mips(\Nd)\refi\Ph(\Net)\); and \(\Ph(\Net)=\Mips(\Nd)\) iff
\(\Ph(\Net)\refi\Mips(\Nd)\) iff both have the same number of blocks iff every doubling pair is fully coupled.""")
rep("C23", r"""E.3 of agent E shows that a Galois connection""",
    r"""Agent E's Proposition~E.3 shows that a Galois connection""")
# ---------- Section 4 ----------
rep("C24", r"""Esparza~\autocite[Thms~10.16--10.17]{DeselEsparza1995}, due to Desel~\autocite{Desel1994}, show that an ordinary net""",
    r"""Esparza~\autocite[Thms~10.16--10.17]{DeselEsparza1995} (Theorem~10.16 is due to Desel~\autocite{Desel1994}) show that an ordinary net""")
rep("C25", r"""The
weighted statement, its persistence reading, \zcref{prop:K5gen} outside free choice and the statement for the dual
were not found (\zcref{sec:search}).""",
    r"""For
ordinary nets the statement for the reverse dual follows in the same way, since the hypotheses of Theorem~10.16 are
invariant under it (Desel--Esparza Prop.~6.20). The weighted statement (for \(\Net\) and for \(\Nd\)), its persistence
reading and \zcref{prop:K5gen} outside free choice were not found (\zcref{sec:search}).""")
rep("C26", r"""None of these has a uniqueness or stability result under mass action. New here: weighted EC under mass action, all
real rates, the converse (2), and the dual statements below.""",
    r"""None of these has a uniqueness or stability result under mass action. Not found in the literature
(\zcref{sec:search}): weighted EC under mass action, all real rates, the converse (2), and the dual statements below.""")
rep("C27", r"""criterion of Müller and Regensburger~\autocite[Prop.~3.9]{MullerRegensburger2012}""",
    r"""criterion of Müller and Regensburger~\autocite[Prop.~3.1]{MullerRegensburger2012}""")
rep("C28", r"""\(G\Gamma\) used by Jamshidi--Palsson and Blanchini--Giordano~\autocite{BlanchiniGiordano2022}.""",
    r"""\(G\Gamma\) used by Jamshidi and Palsson~\autocite{JamshidiPalsson2009} and by Blanchini and
Giordano~\autocite[Thm~1]{BlanchiniGiordano2022}.""")
# ---------- Section 5 ----------
rep("C29", r"""mass action is not preserved~\autocite[p.~6]{BlanchiniGiordano2022}. On the 82\,377 weakly reversible networks with
(H1)--(H3) of agent E's exhaustive families, the dual is weakly reversible for 98, deficiency zero is preserved for 43
of 46\,709, the dual is endotactic for 122 of 84\,103 tested sides and strongly endotactic for 43 of 70\,185 (checked by
computation, \texttt{gac-duality/computation/notes.md} §5).""",
    r"""mass action is not preserved~\autocite[p.~493]{BlanchiniGiordano2022}. On the 82\,377 weakly reversible networks with
(H1)--(H3) of agent E's exhaustive families and two seeded samples, the dual is weakly reversible for 98, deficiency
zero is preserved for 43 of 46\,709, the dual is endotactic for 122 of 84\,103 tests and strongly endotactic for 43 of
70\,185 (counted per run, about 1\,700 networks twice; checked by computation, \texttt{gac-duality/computation/notes.md}
§5).""")
rep("C30", r"""\(D=M^*\colon 3U+2V\to2U+3V,\ V\to U\)""",
    r"""\(M^*\colon 3U+2V\to2U+3V,\ V\to U\) (\(U=K_1\), \(V=K_2\) in \zcref{fig:revdual})""")
rep("C31", r"""\(D\) is non-persistent for every rate""", r"""\(M^*\) is non-persistent for every rate""")
rep("C32", r"""\texttt{gac-duality/round2/persistence/notes.md} §3; the \(s=2\) drain-species count 199 was recounted by V3).""",
    r"""\texttt{gac-duality/round2/persistence/notes.md} §3); over all of E's records, 199 duals with \(s=2\) have a drain
species, a count recounted by V3.""")
# ---------- Section 6 ----------
rep("C33", r"""(up to 3\,544 networks per family), \(\Nd\) is persistent for every rate vector for 17--35\,\% and for none for the rest,
never SOME (Q1; recounted independently by V3 with identical numbers). E2's independent decision procedure agrees on""",
    r"""(up to 3\,544 networks per family), \(\Nd\) is persistent for every rate vector for 11--33\,\% of the networks of each
family with more than one member (381 of 3\,544 to 4 of 12) and for none for the rest, never SOME (Q1; recounted by V3
with its own generator, identical numbers). E2's separate decision procedure agrees on""")
rep("C34", r"""some class (agent D2, Lemma~2.1).""",
    r"""some class (agent D2, Lemma~2.1). On an open edge whose zero set \(Z\) is a siphon of \(\Nd\) (a \emph{siphon edge};
its facet set is a \emph{siphon facet}) only the reactions of \(\Nd\) whose source avoids \(Z\) act, and they move
along the edge; their sum is the \emph{edge flow}, \emph{one-way} if all its terms have one sign, \emph{two-way} if both
signs occur, \emph{zero} if there is no such reaction. The finer notions used in \zcref{prop:none-cert}
(two-coloured edge flows, the transverse order \(\mu\), tropically Hurwitz vertices, centre manifolds, attracting
linear functions and sectors) are defined in agent D2's notes (\texttt{gac-duality/round3/dichotomy/notes.md} §2--3).""")
rep("C35", r"""\paragraph{Counts (checked by computation).} Of the 18\,861 weakly reversible networks with (H1)--(H3) and \(s=2\) in""",
    r"""\paragraph{Counts (checked by computation).} Among the 18\,861 weakly reversible networks with (H1)--(H3) and
\(s=2\) in""")
rep("C36", r"""V3, whose independent sampled test of the ALL verdicts""", r"""V3, whose own sampled test of the ALL verdicts""")
# ---------- Sections 8 and 9 ----------
rep("C37", r"""\item \emph{Weighted impure Theorem C in the Petri sense.} The siphon statement is settled by \zcref{thm:Cgen}; a liveness
theory for weighted impure EC nets (Teruel--Silva assume purity) is not.""",
    r"""\item \emph{Liveness of weighted impure equal-conflict nets.} The siphon statement is settled by \zcref{thm:Cgen} and
\zcref{cor:C}; a liveness theory for weighted impure EC nets (Teruel--Silva assume purity) is not.""")
rep("C38", r"""It certifies every example and finite claim of this paper:""",
    r"""It certifies every example of this paper and every finite claim about an explicit network:""")
rep("C39", r"""The new results were re-derived by three verifiers, Claude Opus 5.5 (effort xhigh), each with
a fresh context and an independent exact library""",
    r"""The results of the working notes were re-derived by three verifiers, Claude Opus 5.5 (effort xhigh), each with
a fresh context and its own exact library""")
rep("C40", r"""transversality, zcref{thm:B}, the WF-FC structure""", r"""transversality, \zcref{thm:B}, the WF-FC structure""")
rep("C41", r"""restriction of E's summary sentence; the citation of Desel--Esparza Thm~6.17; the caveats on K8.""",
    r"""restriction of E's summary sentence; the caveats on K8 (Q3's pairs \(M_3,M_4\), used here, satisfy (H1)--(H3) on
both sides). V1's citation correction (Desel--Esparza Thm~6.17) concerns P's liveness statement D(h), not included here.""")
rep("C42", r"""& archive Lemma~1; P 1.1, 1.2 &""", r"""& archive Lemma~2; P 1.1, 1.2 &""")
rep("C43", r"""& archive Prop.~3, 9, Thm~7, Cor.~8--10, Thm~1 &""", r"""& archive Prop.~3, 6, Thm~7, Cor.~8--10, Thm~1 &""")
rep("C44", r"""\emph{Status}: proved here (with a complete proof or a proof sketch
marked as such), checked""",
    r"""\emph{Status}: proved (here, with a complete proof or a proof sketch
marked as such; for the statements given without proof, in the working notes of \zcref{tab:concordance}, re-derived by
the verifier), checked""")
rep("C47", r"""\documentclass[11pt,a4paper]{report}
""",
    r"""\documentclass[11pt,a4paper]{report}
\pdfsuppressptexinfo=-1 % keep the absolute path of included figures (PTEX.FileName) out of the PDF
""")
rep("C45",r"""Compatibility classes are \((x^0+S)\cap\Rnn^n\). Conjecture~22 assumes:""",
    r"""Compatibility classes are \((x^0+S)\cap\Rnn^n\). The complexes and reactions form a directed graph whose
connected components are the \emph{linkage classes}; \(\Net\) is \emph{weakly reversible} if every reaction lies on a
directed cycle of this graph, and its \emph{deficiency} is the number of complexes minus the number of linkage classes
minus \(s\)~\autocite{Feinberg2019}; complex balance and (strong) endotacticity are used as
in~\autocite{HornJackson1972,CraciunNazarovPantea2013,GopalkrishnanMillerShiu2014}. Conjecture~22 assumes:""")
rep("C46", r"""of all nonzero entries of \(\Ye\) and \(\Yp\) is connected.""",
    r"""of all nonzero entries of \(\Ye\) and \(\Yp\) is connected. S-components (strongly connected subnets generated by a
set of places with all their input and output transitions, in which every transition has one input and one output
place) and T-components (places and transitions exchanged) are as in Desel and
Esparza~\autocite[Defs~5.1, 5.11]{DeselEsparza1995}.""")

BIB_ADD = r"""

@article{FamiliPalsson2003,
  author  = {Famili, Iman and Palsson, Bernhard {\O}.},
  title   = {The convex basis of the left null space of the stoichiometric matrix leads to the definition of
             metabolically meaningful pools},
  journal = {Biophysical Journal},
  volume  = {85},
  number  = {1},
  pages   = {16--26},
  year    = {2003},
  doi     = {10.1016/S0006-3495(03)74450-6}
}

@article{JamshidiPalsson2009,
  author  = {Jamshidi, Neema and Palsson, Bernhard {\O}.},
  title   = {Flux-concentration duality in dynamic nonequilibrium biological networks},
  journal = {Biophysical Journal},
  volume  = {97},
  number  = {5},
  pages   = {L11--L13},
  year    = {2009},
  doi     = {10.1016/j.bpj.2009.06.049}
}

@misc{Duvall2024,
  author     = {Duvall, Alon},
  title      = {Computationally checking if a reaction network is monotone or non-expansive},
  year       = {2024},
  eprint     = {2409.12487},
  eprinttype = {arXiv},
  note       = {Version 2, 27 January 2025}
}
"""
B = []


def brep(i, old, new):
    B.append((i, old, new))


brep("B01", "Hack, Michel Henri Th{\\'e}od{\\'e}or", "Hack, Michel Henri Th{\\'e}odore")
brep("B02", "  volume  = {36},\n  number  = {12},\n  pages   = {1368--1382},\n  year    = {1991},",
     "  volume  = {36},\n  number  = {12},\n  pages   = {1368--1382},\n  year    = {1991},\n  doi     = {10.1109/9.106153},")
brep("B03", "  booktitle = {Proceedings of the 15th IFAC World Congress},\n  year      = {2002}",
     "  booktitle = {Proceedings of the 15th IFAC World Congress},\n  series    = {IFAC Proceedings Volumes},\n"
     "  volume    = {35},\n  number    = {1},\n  pages     = {455--460},\n  year      = {2002},\n"
     "  doi       = {10.3182/20020721-6-ES-1901.00078}")
brep("B04", "  title   = {Asymptotic behavior of nonlinear compartmental systems},",
     "  title   = {Asymptotic behavior of nonlinear compartmental systems: nonoscillation and stability},")
brep("B05", "  pages   = {372--378},\n  year    = {1978},\n  note    = {Cited after the second literature search; not opened}",
     "  number  = {6},\n  pages   = {372--378},\n  year    = {1978},\n  doi     = {10.1109/TCS.1978.1084490},\n"
     "  note    = {Cited after the second literature search; not opened}")
brep("B06", "  pages   = {169--182},\n  year    = {2005},\n  note    = {Cited after the second literature search; not opened}",
     "  number  = {2},\n  pages   = {169--182},\n  year    = {2005},\n  doi     = {10.1109/TAC.2004.841888},\n"
     "  note    = {Cited after the second literature search; not opened}")
brep("B07", "  pages      = {305--329},\n  year       = {2013},",
     "  number     = {1},\n  pages      = {305--329},\n  year       = {2013},\n  doi        = {10.1137/100812355},")
brep("B08", "  journal = {IEEE Transactions on Automatic Control},\n  year    = {2023},\n  doi     = {10.1109/TAC.2023.3275880}",
     "  journal = {IEEE Transactions on Automatic Control},\n  volume  = {68},\n  number  = {12},\n  pages   = {7420--7435},\n"
     "  year    = {2023},\n  doi     = {10.1109/TAC.2023.3275880}")


def apply(text, reps, label):
    bad = 0
    for i, old, new in reps:
        c = text.count(old)
        if c != 1:
            print(f"FAIL {label} {i}: old text occurs {c} times")
            bad += 1
            continue
        text = text.replace(old, new)
        print(f"PASS {label} {i}")
    return text, bad


def line_of(text, old):
    i = text.find(old)
    return text.count("\n", 0, i) + 1 if i >= 0 else None


def markdown():
    """Replacement list for review.md (old text, new text, line in the current file)."""
    out = []
    for i, old, new in R:
        out.append(f"#### {i} (paper.tex, line {line_of(TEX, old)})\n\nOld:\n```latex\n{old}\n```\nNew:\n```latex\n{new}\n```\n")
    for i, old, new in B:
        out.append(f"#### {i} (references.bib, line {line_of(BIB, old)})\n\nOld:\n```bibtex\n{old}\n```\nNew:\n```bibtex\n{new}\n```\n")
    out.append(f"#### B09 (references.bib, append at the end)\n\n```bibtex{BIB_ADD}```\n")
    return "\n".join(out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--markdown":
        sys.stdout.reconfigure(encoding="utf-8")
        print(markdown())
        sys.exit(0)
    tex, b1 = apply(TEX, R, "paper.tex")
    bib, b2 = apply(BIB, B, "references.bib")
    bib += BIB_ADD
    out = HERE / "build"
    out.mkdir(exist_ok=True)
    (out / "paper.tex").write_text(tex, encoding="utf-8")
    (out / "references.bib").write_text(bib, encoding="utf-8")
    print(f"{len(R)} paper.tex replacements, {len(B)} references.bib replacements, 3 bib entries added; FAILS: {b1 + b2}")
    sys.exit(1 if b1 + b2 else 0)
