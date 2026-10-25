# β-Kirchhoff Šoltés graphs: exact proportions below 2/3, and a 2/3 barrier for splice constructions

**Credits.** Solver: Claude Opus 5.5 (effort xhigh). Literature search helper: a Claude Opus 5.5 sub-agent (effort not
recorded). Verifier: — (pending). Editor: — (pending).

**Status: partial result** (preliminary; verification pending). Problem 4 of the source is **still open for every
rational β ∈ [2/3, 1)**. The source's own construction already settles it for β < 2/3. We prove exact-proportion
results on (0, 2/3) and a barrier theorem explaining why splice constructions stop at 2/3, and we give a 2-connected
near miss.

## Problem

K. K. Gottwald, S. Majstorović Ergotić, T. Došlić, *Šoltés problem for the Kirchhoff index of a graph*, Appl. Math.
Comput. 510 (2026) 129694, doi:10.1016/j.amc.2025.129694, arXiv:2609.13751. Journal p. 3 = arXiv p. 6 (identical
wording):

> "let 0 < β ≤ 1, β ∈ ℚ. We say that a graph G is a β−Kirchhoff Šoltés graph if |K(G)| ≥ β|V(G)|. i.e. the proportion
> of good vertices in G is at least β. We pose the following problems.
> **Problem 4.** For a fixed rational number β ∈ (0,1] construct an infinite series of β−Kirchhoff Šoltés graphs."

Here K(G) = {v : Kf(G) = Kf(G − v)} (Kf = Kirchhoff index, the sum of resistance distances). A cut vertex is never good
("x is also not a good vertex since it is a cut vertex"). Problem 5 asks for an *increasing* proportion tending to a
given γ ∈ (0, 1]. Problem 17 (p. 11 of the journal version) asks which numbers in (0, 1] are exact or limiting
proportions. β = 1 is Problem 3 (only C₅ is known), which is the subject of `../kirchhoff-soltes-problem/`.

**What the source already proves.** The definition is "at least β". Theorem 15
gives the splices k·(K₄ − e) at a degree-2 vertex, with exactly 2k good vertices among 3k + 1. Since 2k/(3k + 1)
increases to 2/3, these graphs already form an infinite series of β-Kirchhoff Šoltés graphs for every rational
β < 2/3. The source also has an infinite family with proportion exactly 1/2 (K_s + K_{s+1,1}) and limits
(s + 1)/(2s + 1) and 3 − 2√2. So the open range is β ∈ [2/3, 1).

Novelty check (25 Oct 2026, details in `references.md`):
- no work cites the source;
- no arXiv, Zenodo or GitHub item addresses β ≥ 2/3, including the AI campaigns of ipitchford/Evidence Press and mmislan,
  and google-deepmind/formal-conjectures;
- the follow-up by Majstorović Ergotić and Gottwald (Bol. Soc. Mat. Mex. 31 (2025) 88, read in full, CC-BY) only builds
  graphs with one or two good vertices.

## Result

1. **Splice lemma with an arbitrary rest** (Lemma 1, Cor. 2) [proof]. If G = B ⊕ W is a splice at x, then for
   v ∈ B: Kf(G) − Kf(G − v) = δ_B(v) − m_W ρ_B(v) + Rt_W(x). Here δ_B(v) = Kf(B) − Kf(B − v) and
   ρ_B(v) = Rt_{B−v}(x) − Rt_B(x). The rest of the graph enters only through (m_W, Rt_W(x)). This generalises
   Props. 13–14 of the source and Thm. 16 of the follow-up.
2. **Three regimes** (Prop. 3) [proof]. Let a = Rt_B(x) − m_B ρ. A class of a gadget with a = 0 or a < 0 gives families
   whose proportion of good vertices tends to g/m_B. For a < 0 the tuning is a pendant path of length O(√k), with an
   exact integrality argument. A class with a > 0 needs a linear number of extra vertices, and a class with ρ ≤ 0 dies
   for large k. Example: copies of C₅ plus a pendant path, with proportions 1/4, 4/11, 2/5, … → 1/2.
3. **Every rational β ∈ (0, 2/3) is an exact proportion** (Thm. 5, Cor. 6) [proof].
   - G(k, α, γ) = k·(K₄−e) ⊕ α·F₅ ⊕ γ·F₆, where F₅ and F₆ are explicit "neutral" fillers contributing 5 and 6 non-hub vertices that make
     the hub's transmission exactly (3/4)(n − 1).
   - G(k, α, γ) has **exactly** 2k good vertices among n = 1 + 3k + 5α + 6γ. With k = pj and 5α + 6γ = (2q − 3p)j − 1
     (Frobenius number 19), the proportion is exactly p/q.
   - Hence Problem 5 holds for every γ ∈ (0, 2/3] (strictly increasing proportions tending to γ). For Problem 17,
     every rational in (0, 2/3) is an exact value (infinitely often) and every real in (0, 2/3] is a limiting value.
4. **Small graphs** (Prop. 7) [rigorous modular computation + exact checks]. The maximum number of good vertices of a
   connected graph on n = 3, …, 10 vertices is 0, 2, 5, 3, 4, 4, 4, 6. C₅ is the only connected graph with n ≤ 10 in
   which at least 2/3 of the vertices are good (all 11 989 762 connected graphs with 3 ≤ n ≤ 10). A run for n = 11
   (about 10⁹ graphs) was started but stopped before completion, so n = 11 is not covered.
5. **The 2/3 barrier for splice constructions** (Thm. 8, Prop. 10, Cor. 11).
   - Theorem [proof]: for any finite family 𝓑 of rooted graphs, every graph with a cut vertex x whose good vertices lie
     in branches from 𝓑 has |K(G)| ≤ λ(𝓑)|V(G)| + C(𝓑). Here λ(𝓑) is an explicit linear-fractional optimum: single
     classes, or balanced pairs of classes with equal ρ.
   - Exact computation in int128 arithmetic: λ(𝓑₉) = 2/3 for all branches with ≤ 9 vertices (8 915 804 classes).
     Every class of every rooted graph with ≤ 10 vertices has efficiency ≤ 2/3.
   - Hence no infinite family of these splices can maintain a proportion β > 2/3. Within the enumerated gadgets, classes covering more than 2/3 of the non-hub vertices have
     ρ/ρ* < 0.63.
6. **A 2-connected near miss** (Prop. 13) [proof + exact computation].
   - The rook graph K_a □ K_b with one hub per row has its ab rook vertices in one orbit (proportion b/(b + 1)).
   - Kf(R) − Kf(R − v) = −P(a, b)/Den is derived symbolically (isotypic decomposition under S_a × S_b) and checked
     against brute force.
   - Its zero set has a real branch a* ≈ b³/2, on which the proportion tends to 1. P is irreducible, and there is no
     integer solution for b ≤ 400.
   - With c hubs per row (proportion b/(b + c)) there are real zeros for c < b/2, e.g. a* ≈ 61.8 for (b, c) = (10, 3),
     proportion 10/13. There is no integer solution for 3 ≤ b ≤ 200, 1 ≤ c ≤ b. A finite graph with more than 2/3 good
     vertices (other than C₅) would already be new.

**What is not proved.**
- Problem 4 for β ∈ [2/3, 1), including whether any infinite series reaches proportion ≥ 2/3.
- The barrier for branches with ≥ 10 vertices (mixtures), or ≥ 11 vertices (single gadget type).
- That the rook branch has only finitely many integer points (expected; Puiseux argument only).
- The floating-point explorations are evidence only: cographs, blow-ups, joins with vertex-transitive graphs, theta
  graphs, row–column variants (Remark 12 of the paper, `scripts/exploration/`).

## How to reproduce

From this folder, under Fedora 44 / WSL. Requirements: Python 3.12, python-flint 0.6, sympy 1.14, networkx 3, numpy,
numba (exploration only), nauty 2.8.9 `geng` on the PATH, gcc.

```bash
gcc -O3 -march=native -o scripts/goodcount scripts/goodcount.c
gcc -O3 -march=native -o scripts/gadgetbound scripts/gadgetbound.c
cd scripts
python3 check_splice_lemma.py            # Lemma 1, Cor. 2, Example 4                (3 s)
python3 check_exact_proportions.py       # Theorem 5                                 (15 s)
python3 check_rook_hubs.py               # Proposition 13 (1)-(4)                    (45 s)
python3 check_rook_hubs_c.py 200         # Proposition 13 (5), c hubs per row         (15 min; 60: 75 s)
python3 test_c_programs.py               # C programs vs exact Python, n <= 8        (30 s)
python3 check_mixing_bound.py 9          # Prop. 10(2), lambda(B_9) = 2/3            (25 min)
cd ..
scripts/run_goodcount.sh 3 10            # Prop. 7 (modular upper bounds)           (5 min)
python3 scripts/verify_small_graphs.py 10    # Prop. 7 (exact checks of the listed graphs)
for n in 3 4 5 6 7 8 9 10; do geng -c -q $n | scripts/gadgetbound 2 3 \
   > certificates/search/gadgetbound_n$n.out 2> certificates/search/gadgetbound_n$n.log; done   # Prop. 10(1) (~50 min)
# n = 11 was started and stopped (not completed); about 5 h of one core, e.g. in 4 slices:
# for r in 0 1 2 3; do geng -c -q 11 $r/4 | scripts/goodcount 3 5 > certificates/search/goodcount_n11_part$r.out \
#     2> certificates/search/goodcount_n11_part$r.log; done; python3 scripts/verify_small_graphs.py 11
py -3 ../../templates/build_tex.py paper.tex   # or: python templates/build_tex.py problems/kirchhoff-soltes-beta-families/paper.tex
```

## Files

- `paper.tex`, `references.bib`, `paper.pdf`: manuscript.
- `references.md`: sources read and novelty check.
- `scripts/kfexact.py`: exact Kf, Rt, splices (python-flint).
- `scripts/check_*.py`, `verify_small_graphs.py`, `test_c_programs.py`: PASS/FAIL certificates (exit code 1 on failure).
- `scripts/goodcount.c`: rigorous modular upper bound on the number of good vertices, for every graph read.
  `scripts/gadgetbound.c`: exact int128 computation of all classes and efficiencies of all rooted graphs.
  `scripts/run_goodcount.sh`: driver for goodcount.
- `scripts/exploration/`: floating-point exploration tools (not certificates): gadget scan, class mixing, cographs,
  blow-ups, joins, theta graphs, row–column variants, hub tables.
- `certificates/*.log`: outputs of the certificate scripts. `certificates/search/`: logs of goodcount and gadgetbound.
  `certificates/exploration/`: logs of the floating-point explorations.
- `certificates/sources/` (git-ignored): source PDFs (arXiv and journal versions, follow-up, Wiener-index papers).
