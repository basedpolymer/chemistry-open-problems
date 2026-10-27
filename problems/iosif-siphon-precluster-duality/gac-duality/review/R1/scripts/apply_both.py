#!/usr/bin/env python3
r"""R1: apply the corrections of both reviews (R1: review/R1/review.md; R2: review/R2/review.md) to paper.tex and
references.bib of the problem folder, with the merges decided by R1 (see review/R1/applied.md).

R2's replacement list is imported from review/R2/scripts/apply_corrections.py (module-level lists R, B, BIB_ADD; the
module only reads files when imported).  A few R2 entries are replaced by merged versions (MERGED below); R1's own
replacements and proof insertions follow (R1 list).  Every old text must occur exactly once in the file it is
applied to; otherwise nothing is written.

Environment: Windows `py -3` or WSL python3, standard library only.  Runtime < 1 s.
Usage (from this folder):  py -3 scripts/apply_both.py [--dry-run]
Writes ../../../paper.tex and ../../../references.bib (unless --dry-run) and prints PASS/FAIL per replacement.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent            # gac-duality/review/R1
FOLDER = HERE.parents[2]                                   # problem folder
spec = importlib.util.spec_from_file_location("r2apply", HERE.parent / "R2" / "scripts" / "apply_corrections.py")
r2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r2)

TEX_PATH, BIB_PATH = FOLDER / "paper.tex", FOLDER / "references.bib"
TEX, BIB = TEX_PATH.read_text(encoding="utf-8"), BIB_PATH.read_text(encoding="utf-8")

# ------------------------------------------------------------------------------------------------ merged R2 entries
MERGED = {}
MERGED["C03"] = r"""preprints, which we have not checked. Every example is certified in exact arithmetic (130 claims); counts over large
families are cited from the working notes as checked by computation. The main results were re-derived by three
fresh-context verifiers."""
MERGED["C04"] = r"""Petri clusters when \(\Net\) is equal-conflict, every species is consumed and no source is zero (e.g.\ for connected
networks with (H1) and (H2))."""
MERGED["C27"] = r"""criterion of Müller and Regensburger~\autocite[Prop.~3.1 of the arXiv version]{MullerRegensburger2012}"""
MERGED["C33"] = r"""(up to 3\,544 networks per family), \(\Nd\) is persistent for every rate vector for 11--33\,\% of the networks of each
family with more than one member (381 of 3\,544 to 4 of 12; the family \((4,3,4)\) has a single network, which is ALL)
and for none for the rest, never SOME (Q1; recounted by V3 with its own generator, identical numbers). E2's separate
decision procedure agrees on"""

R2_LIST = [(i, old, MERGED.get(i, new)) for i, old, new in r2.R]

# ------------------------------------------------------------------------------------------------ R1 list
R1 = []


def rep(i, old, new):
    R1.append((i, old, new))


def after(i, anchor, insertion):
    """Insert `insertion` right after the unique `anchor`."""
    R1.append((i, anchor, anchor + insertion))


rep("R1-R2", r"""class, and a periodic orbit).""", r"""class, and, for some rate vector, a periodic orbit).""")
rep("R1-R3a", r"""the twin classes of \(\Net\) are the doubling classes of \(\Nd\).""",
    r"""the twin pairs of \(\Net\) are the doubling pairs of \(\Nd\) (two zero rows give two reactions of \(\Nd\) with the
zero source), and the twin classes of nonzero rows are its doubling classes.""")
rep("R1-R3b", r"""twin classes / doubling classes & doubling classes / twin classes & definitions \\""",
    r"""twin pairs / doubling pairs & doubling pairs / twin pairs & definitions \\""")
rep("R1-R4", r"""vector), so for EC networks clusters correspond to sources and \(q=m\).""",
    r"""vector), and under (H1) a reaction with the zero source is \(0\to0\) (\(z\tp y'=0\) for a positive law \(z\)); a species
in no reaction, or a reaction \(0\to0\), is isolated in the graph of \(\Ye\) and \(\Yp\). So for connected EC networks with
(H1) and (H2), clusters correspond to sources and \(q=m\); with (H2) alone, \(0\to X_1,\ 0\to X_2,\ X_1+X_2\to0\) is
connected and EC with \(m=2\) and \(q=3\).""")
rep("R1-R5", r"""so \((\log V)^\cdot\ge\gamma T^3/8-E\varepsilon\ge\gamma T^3/16>0\) on \(\{b+c<\varepsilon\}\) with""",
    r"""so \((\log V)^\cdot=c\,[\gamma a^3-(3\delta+2\beta\gamma a^3/\alpha)c-3\beta\delta c^2/\alpha]\ge c\,(\gamma T^3/8-E\varepsilon)
\ge c\,\gamma T^3/16>0\) at the interior points of \(\{b+c<\varepsilon\}\), with""")
rep("R1-R5b", r"""A last-entry argument
excludes visits accumulating at \((T,0,0)\).""",
    r"""So \(V\) increases
along interior trajectories while they stay in this region. If \(x(t_j)\to(T,0,0)\) along an interior trajectory, the
trajectory cannot stay in the region from some time on (\(V\) would increase while \(V(x(t_j))\to0\)); the last entry
points \(x(s_j)\) before \(t_j\) satisfy \(b+c=\varepsilon\) and \(V(x(s_j))\le V(x(t_j))\to0\), so they accumulate at
\((T-\varepsilon,\varepsilon,0)\in\Omega\), a point of the open edge \(c=0\), which was excluded.""")
rep("R1-R7a", r"""ends). In rank two the corollary for pure networks fails without weak reversibility:""",
    r"""ends). In rank two the corollary for pure networks fails, already without weak reversibility:""")
rep("R1-R7b", r"""non-persistent unless \(\kappa_2=\kappa_3\).""",
    r"""non-persistent for every rate vector with \(\kappa_2\neq\kappa_3\) (the hypersurface \(\kappa_2=\kappa_3\) was not
analysed); \zcref{ex:drain}(4) shows that it also fails for weakly reversible pure networks.""")
rep("R1-R8", r"""\((y_l)_o=\sum_q\lambda_q(y_l)_q\) for every \(l\notin Z\).""",
    r"""\((y_l)_o=\sum_q\lambda_q(y_l)_q\) for every \(l\notin Z\), with \(\lambda_q\ge0\) and
\(\sum_q\lambda_q=1\).""")
rep("R1-R9", r"""Let \(s=2\). Factor""",
    r"""Let \(s=2\) and let \(\Net\) satisfy (H1)--(H3), so that the positive classes of \(\Nd\) are compact convex polygons.
Factor""")
rep("R1-R10", r"""(A2) or by the normal coordinate of an IN edge, and is isolated: an interior total orbit near \(e\) would have its
\(\alpha\)- and \(\omega\)-limits in \(\{V_e=0\}\cap\partial P\), which contains no invariant set besides \(e\) (V3's addition).""",
    r"""(A2) or by the normal coordinate of an IN edge, and is isolated: along an interior total orbit in a compact invariant
set near \(e\), \(V_e\) is strictly increasing, so its \(\alpha\)- and \(\omega\)-limit sets lie in \(\partial P\); near \(e\)
the only invariant subset of \(\partial P\) is \(\{e\}\) (boundary orbits are monotone along the edges), so both limit
sets are \(\{e\}\), which contradicts the strict increase of \(V_e\) between them since \(V_e(e)=0\) (V3's addition).""")
rep("R1-R11a", r"""consumes \(X_i\), and (b) every reaction producing \(X_i\) has a source with \(X_i\)-coefficient \(>d_0\).""",
    r"""consumes \(X_i\), and (b) every reaction with a net production of \(X_i\) (\(\Gamma_{ij}>0\)) has a source with
\(X_i\)-coefficient \(>d_0\).""")
rep("R1-R11b", r"""production terms of \(x_i\) are \(\le Cx_i^{d_0+1}\)""",
    r"""the terms of \(\dot x_i\) with \(\Gamma_{ij}>0\) are \(\le Cx_i^{d_0+1}\)""")
rep("R1-R12", r"""For \(3X_3\to2X_1+X_2\to X_2+2X_3\to3X_3\), the dual""",
    r"""For \(3X_3\to2X_1+X_2\to X_2+2X_3\to3X_3\), with the reactions numbered \(3X_3\to2X_1+X_2\), \(X_2+2X_3\to3X_3\),
\(2X_1+X_2\to X_2+2X_3\), the dual""")
rep("R1-R13", r"""with the steady state \((\kappa_B/\kappa_A,3\kappa_AK_1^2/\kappa_C,1)\).""",
    r"""with the steady state \((\kappa_B/\kappa_A,3\kappa_B^2/(\kappa_A\kappa_C),1)\).""")
rep("R1-R14", r"""since the index sum is \(+1\), the
class contains at least one further positive steady state (V2: exactly two more, by exact root isolation).""",
    r"""the class is compact without boundary steady states
(\zcref{thm:Cgen}), so by the degree argument of \zcref{thm:Ustar}(b) the indices of its steady states add up to
\((-1)^s=+1\), and the class contains at least one further positive steady state (V2: exactly two more, by exact root
isolation; an exact Gröbner-basis count of reviewer R1 agrees).""")
rep("R1-R15", r"""(D2 proved it for
two-coloured and zero edge flows with \(\mu=1\)).""",
    r"""(D2 proved it for zero
edge flows with \(\mu=1\), and for two-coloured edge flows with \(\mu=1\) with V3's repair).""")
rep("R1-R19", r"""Global stability of \(\Nd\) for \(s\ge3\), and of general PC""",
    r"""Global stability of \(\Nd\) for \(s\ge3\) with unbalanced joins, and of general PC""")
rep("R1-R21b1", r"""the smallest conservative example is \(2X_2\rightleftarrows X_1\),""",
    r"""the smallest conservative example in agent E's exhaustive families is \(2X_2\rightleftarrows X_1\),""")
rep("R1-R21b2", r"""it is the smallest weakly
reversible network with (H1)--(H3) and a NONE dual.""",
    r"""it is the smallest weakly
reversible network with (H1)--(H3) and a NONE dual (V3: one species is impossible, \(X_1\rightleftarrows X_2\) is ALL,
and V3's enumeration of the family \((2,2,4)\) finds \(X\) first).""")
rep("R1-cert1", r"""about 1~s): 129 claims, all PASS""", r"""about 1~s): 130 claims, all PASS""")
rep("R1-cert2", r"""\(\mathbb{Q}(\sqrt{17})\)), the extinction identities and drain certificates,""",
    r"""\(\mathbb{Q}(\sqrt{17})\)), the examples of \zcref{prop:ADLS} and of its converse, the extinction identities and drain
certificates,""")
rep("R1-R17b", r"""Items marked ``sketch'' (\zcref{thm:WFFC,thm:P2,thm:AF,prop:prime} and the certificates of \zcref{prop:none-cert}) have
complete proofs in the working notes, re-derived by the verifiers; V3 checked the tropical certificates (CM), (L), (S) at
sketch level and on their worked examples only.""",
    r"""Items marked ``sketch'' (\zcref{thm:WFFC,thm:P2,thm:AF,prop:prime,prop:none-cert}) have complete proofs in the
working notes, re-derived by the verifiers, except the tropical certificates (CM), (L), (S) of \zcref{prop:none-cert},
whose proofs in D2's notes are themselves sketches; V3 checked them at sketch level and on their worked examples only.

The rewritten manuscript was then reviewed by two fresh-context reviewers, Claude Opus 5.5 (effort xhigh): R1
(mathematics) and R2 (attribution, consistency and presentation); their corrections are applied here
(\texttt{gac-duality/review/R1/applied.md}).""")
rep("R1-R18a", r"""\(T=\min(1,\delta/(2(u-2)\gamma))\). Certified for \(u=3,\dots,7\), \(R=1,\dots,7\).""",
    r"""\(T=\min(1,\delta/(2(u-2)\gamma))\) (using \(T^u\le T\) for \(T\le1\)). Certified for \(u=3,\dots,7\), \(R=1,\dots,7\).
Full proof in agent A2's notes (Theorem AF), re-derived by V3 (§6).""")
rep("R1-R18b", r"""direction not orthogonal to \(S\) is nonconstant on every class, which gives strong endotacticity.""",
    r"""direction not orthogonal to \(S\) is nonconstant on every class, which gives strong endotacticity. Full proof in
agent G's notes (Proposition~4.1), re-derived by V3 (§1).""")
rep("R1-R20b", r"""the source rows of \(\Ye\) are pairwise distinct""", r"""the rows of \(\Ye\) are pairwise distinct""")

# ---- proofs and pointers inserted after statements (GAPs R-17, R-20, R-21, R-22 and R2 §6)
after("R1-R17a", r"""one of them NONE.
\end{proposition}""", r"""
\begin{proof}[Proof (sketch)]
Proofs in agent D2's notes (\texttt{gac-duality/round3/dichotomy/notes.md}, Propositions E, H, CM, L, S, B). (E) for
\(\mu=1\), as repaired by V3: on the facet there is a point \(p^*\) with zero set \(Z\), zero edge flow and normal
coefficient \(L(p^*)<0\); it is an equilibrium of its class, and the stable manifold of the eigenvalue \(L(p^*)\) is
transverse to the invariant edge and enters the interior. (H): Cauchy--Binet for a sum of rank-one \(2\times2\)
matrices; a Hurwitz vertex attracts an open set of its class. (B): the two blocks decouple. (CM), (L), (S) are standard
centre-manifold, Lyapunov and sector arguments, written as sketches in D2's notes and checked by V3 at sketch level and
on the worked examples below.
\end{proof}""")
after("R1-R20", r"""for every cycle \(C\) of \(\pi\).
\end{proposition}""", r"""
\begin{proof}
If \(\Net\) is weakly reversible, every product complex is a source complex, so equal rows of \(\Ye\) would give equal
rows of \(\Yp\) and of \(\Gamma\), excluded by (H3); the reversed network gives the same for \(\Yp\). So \(\Nd\) has \(n\)
distinct sources and \(n\) distinct products. If \(\Nd\) is weakly reversible, each product is a source, which defines
an injective map \(\pi\) with \((\Yp)_{i\cdot}=(\Ye)_{\pi(i)\cdot}\); conversely the complex graph \(i\mapsto\pi(i)\) of a
permutation is a disjoint union of cycles. It has \(n\) complexes and \(c(\pi)\) linkage classes, so the deficiency of
\(\Nd\) is \(n-c(\pi)-s\), and \(\sum_{i\in C}\Gamma_{i\cdot}=\sum_{i\in C}((\Ye)_{\pi(i)\cdot}-(\Ye)_{i\cdot})=0\) for
every cycle \(C\) of \(\pi\). The examples are certified (agent G, Prop.~2.2; V3 §1).
\end{proof}""")
after("R1-R21", r"""rate vector iff every source contains strictly more of the species its reaction consumes than of the one it produces.
\end{corollary}""", r"""
\begin{proof}
By (H1) the vector \(a\) has entries of both signs. Let \(Z\) be a face of \(\Nd\) (\zcref{cor:S1N}). A species with
\(a_i\neq0\) is a net reactant of every reaction of \(Z\) or a net product of every reaction of \(Z\). A net reactant
occurs in every source of \(Z\), so its load is at least \(|Z|\ge1\); a net product occurs in no source of \(Z\) by
hypothesis, so its load is \(0\), and net products exist by (H1). So every species of minimal load with \(a_i\neq0\) is a
net product, every face is IN, and \zcref{cor:S1N} applies. If \(n=2\), every reaction consumes one species and
produces the other, the reactions of one orientation all consume the same species, loads add over \(Z\), so the
singleton faces decide, and \(\{j\}\) is IN iff \(y_j\) contains more of the consumed than of the produced species.
\end{proof}""")
after("R1-R22-G", r"""(law \(X_1+X_2\), only minimal siphon \(\{X_1\}\)).
\end{proposition}""", r"""
\begin{proof}
Conservation-coupled species (\(w_i=\mu w_j\)) lie in the same supports of extreme laws, and by (H1) every species lies
in one. If every siphon contains the support of a nonnegative law, a minimal siphon equals such a support, which is
then the support of an extreme law, and conversely. The examples are certified (\texttt{certify\_duality.py}).
\end{proof}""")
after("R1-R22-K8", r"""siphons but different \(\Pre\) of their duals.
\end{proposition}""", r"""
\begin{proof}
Checked by computation in exact arithmetic (\texttt{certify\_duality.py}, claims ``S3 K8''); agents P and Q3, re-checked
by V1 and V2.
\end{proof}""")
after("R1-R22-M", r"""some nonzero \(\xi\ge0\) satisfies \(\xi\tp G=0\).
\end{lemma}""", r"""
\begin{proof}
For \(c>0\) large, \(G+cI\) is nonnegative, irreducible and has a positive diagonal, hence is primitive. A positive null
vector of \(G\) (or of \(G\tp\)) is a positive eigenvector of \(G+cI\) (or of its transpose) for the eigenvalue \(c\), so
by the Perron--Frobenius theorem \(c\) is the spectral radius and a simple eigenvalue, and every other eigenvalue
\(\lambda+c\) satisfies \(|\lambda+c|<c\), i.e.\ \(\operatorname{Re}\lambda<0\). In the reducible case \(G+cI\) is
nonnegative with all row sums \(c\), so its spectral radius is \(c\), and the Perron--Frobenius theorem for nonnegative
matrices gives a nonzero \(\xi\ge0\) with \(\xi\tp(G+cI)=c\,\xi\tp\).
\end{proof}""")
after("R1-R22-C", r"""of \(\Net\) are the internal cycles of \(\Nd\) and conversely, and both are persistent for every rate vector.
\end{corollary}""", r"""
\begin{proof}
\zcref{thm:Cgen} with \(q=m\) (\zcref{sec:setting}); the internal cycles of \(\Nd\) are the supports of the extreme laws
of \(\Net\) (\zcref{lem:transpose}).
\end{proof}""")
after("R1-R22-Tbal", r"""\(n-s\) and \(r-s\).
\end{corollary}""", r"""
\begin{proof}
\(\Nd\) is EC with balanced joins, connected, of rank \(s\), with the same clusters and (H1), (H2) exchanged
(\zcref{tab:dictionary}); apply \zcref{thm:T} to \(\Net\) and to \(\Nd\).
\end{proof}""")
after("R1-R22-BH", r"""exactly the duals of WF-EC
networks.
\end{corollary}""", r"""
\begin{proof}
If \(\Net\) is BH, two species of one block have equal rows of \(\Ye\), so \(\Nd\) is EC with clusters the blocks,
connected, of rank \(s\), with (H1), (H2) exchanged and \(s=m(\Nd)-1\): \(\Nd\) is WF-EC; conversely the dual of a WF-EC
network is BH. Apply \zcref{thm:Tstar,thm:Ustar,cor:C,prop:GSstar} to the WF-EC network \(\Nd\) (or \zcref{thm:PC}).
\end{proof}""")
rep("R1-R22-WR", r"""critical-point equation is the steady-state equation.""",
    r"""critical-point equation is the steady-state equation. Details: agent A2's notes (Theorems E, E+), re-derived by V3
(§6).""")
after("R1-R22-ADLS", r"""(H1)--(H3) and has the critical siphon \(\{U_4\}\).
\end{proposition}""", r"""
\begin{proof}
Checked in exact arithmetic (\texttt{certify\_duality.py}, claim ``S5 H''); agent G, re-checked by V3.
\end{proof}""")
after("R1-R22-end", r"""endpoint of any class attracts nearby points of its class.
\end{lemma}""", r"""
\begin{proof}
Agent Q1's notes, Lemmas 1.1--1.2 (along a class \(F\) is a polynomial in \(\varepsilon\), and the dynamics on a segment
is one-dimensional); re-derived by V3 (§3).
\end{proof}""")
after("R1-R22-S1N", r"""every species with \(a_i\ne0\) minimising \(\ell_Z\) is a net product of \(Z\).
\end{corollary}""", r"""
\begin{proof}
\zcref{thm:S1}(1) for \(\mathcal M=\Nd\), which is conservative by (H2), with the dictionary above.
\end{proof}""")
after("R1-R22-A", r"""stays in the open edge. Every vertex in \(\Omega\) is an equilibrium.
\end{lemma}""", r"""
\begin{proof}
Agent D2's notes, Lemma A (the siphon edge is invariant in both time directions for the polynomial field); re-derived by
V3 (§5).
\end{proof}""")
after("R1-R22-L", r"""one-way or zero edge flow, for every rate
vector.
\end{proposition}""", r"""
\begin{proof}
Agent A2's notes, Proposition L (opposite active vectors on a siphon edge would give a nonnegative conservation law of
\(\Net\) vanishing on a complex, hence on the whole linkage class); re-derived by V3 (§6).
\end{proof}""")
after("R1-R22-bd", r"""(V3 checked the reading of condition (\(\diamond\)) on the PDF).
\end{proposition}""", r"""
\begin{proof}
Agent G's notes, Propositions 6.1--6.3 ((2): the tangent space of the face is the kernel of the projection of \(S\) onto
\(\mathbb{R}^\Sigma\)); re-derived by V3 (§1).
\end{proof}""")

# ------------------------------------------------------------------------------------------------ bibliography
B_LIST = list(r2.B) + [("R1-B10", "  doi     = {10.1137/110847056}",
                        "  doi     = {10.1137/110847056},\n  eprint  = {1209.6488},\n  eprinttype = {arXiv}")]


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


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    # every old text must be unique in the ORIGINAL file and the R1 and R2 old texts must not overlap
    olds = [(i, o) for i, o, _ in R2_LIST + R1]
    bad = sum(TEX.count(o) != 1 for _, o in olds)
    for (i, o) in olds:
        for (j, p) in olds:
            if i < j and (o in p or p in o):
                print(f"OVERLAP {i} {j}")
                bad += 1
    tex, b1 = apply(TEX, R2_LIST + R1, "paper.tex")
    bib, b2 = apply(BIB, B_LIST, "references.bib")
    bib += r2.BIB_ADD
    total = bad + b1 + b2
    print(f"{len(R2_LIST)} R2 + {len(R1)} R1 paper.tex replacements, {len(B_LIST)} references.bib replacements, "
          f"3 bib entries added; problems: {total}")
    if total == 0 and "--dry-run" not in sys.argv:
        for path, text in ((TEX_PATH, tex), (BIB_PATH, bib)):
            crlf = b"\r\n" in path.read_bytes()          # keep each file's original line endings
            path.write_bytes(text.replace("\n", "\r\n" if crlf else "\n").encode("utf-8"))
        print("written: paper.tex, references.bib")
    sys.exit(1 if total else 0)
