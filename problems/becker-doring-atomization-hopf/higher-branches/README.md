**Attribution.** Solver (method, scripts, branches 4–8, onset counts, Theorem C, local certificate at $M=1135$):
Claude Opus 5.5 (effort xhigh), sub-agent of the orchestrator. Editing after the solver stopped (text revisions,
Theorem C in `result.md`, completeness checks in `summary.py`): GPT-6-Astra (effort xhigh). Branches 9–13 (runs of
the solver's pipeline), final tables and manuscript: Claude Opus 5.5 (effort xhigh). Builds on result 8
(`../README.md`) and its all-M extension (`../all-M/README.md`); their solvers, editor and verifier are credited
there. Verification of this extension: pending.

# Becker–Döring with atomization: the Hopf branches j ≥ 4

## Problem

R. L. Pego and J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605. Their Theorem 3.3 gives a Hopf
bifurcation at $\kappa=\kappa_j(M)$ ($\kappa=K\sqrt M$) for every $j\le k$ and $M>M_{0,k}$; right after it, §3.3 says
(verbatim):

> "We have not managed to determine analytically whether the bifurcating solutions are stable (the supercritical
> case) or not."

Result 8 and its preliminary all-M extension address this for the branches $j=1,2,3$ and every $M$. This folder treats the
branches $j\ge4$: where each branch starts, whether the spectral hypotheses of the Hopf theorem hold on it, the sign
of the first Lyapunov coefficient $\ell_1$, and the stability of the bifurcating orbits.

**Novelty check (2026-10-25, short; result 8's full check of 2026-10-18 is not repeated).** OpenAlex lists 17
works citing the paper, the same as in result 8's `references.md`; none treats the higher branches. The GitHub
listings of mmislan (27 repositories) and ipitchford (193) contain nothing on Becker–Döring or Hopf bifurcations;
`gh search repos` and the Zenodo API ("Becker-Döring" atomization; Becker Döring Hopf) give only simulation codes
and unrelated records. Details in `references.md`.

## Status: partial result (preliminary), not yet verified (25 October 2026)

* **Branches $j=4,\dots,13$, every $M$** (table below): complete computer-assisted answer. The bifurcation is
  supercritical on the centre manifold except at the onsets $M=142$ ($j=4$) and $M=1135$ ($j=13$), where it is
  subcritical; the bifurcating orbits are unstable.
* **Every branch $j\ge1$, large $M$**: supercritical for all $M$ beyond a threshold that is not explicit.
* **Open**: branches $j\ge14$ below that threshold; no argument uniform in $j$ with explicit constants.

## Result

**1. Where branch $j$ starts.** In the limit $K\to\infty$ the matrix $B/A$ tends to an explicit matrix $C(M)$ with
$\det(\nu-C)=-G(1+\nu)/\nu$, where $G(\varphi)=-\varphi^{M-2}(\varphi^3+M\varphi^2-2M\varphi+M)+(M+1)\varphi-M$. For every
integer $3\le M\le1150$ the number of unstable eigenvalue pairs of $C(M)$ is certified: it increases by one exactly at
$$M_j=15,\ 44,\ 87,\ 142,\ 210,\ 289,\ 379,\ 480,\ 591,\ 713,\ 844,\ 985,\ 1135\qquad(j=1,\dots,13),$$
which is also the number of unstable pairs of $B(K,M)$ for all large $K$. The $j$-th pair of $C$ crosses the
imaginary axis at a real value $M^*_j\in(M_j-1,M_j)$ ($M^*_4=141.98176\ldots$). Along the $j$-th Hopf curve,
continued from Pego–Velázquez's limit root at $M=\infty$, $K\to\infty$ as $M\downarrow M^*_j$. For $4\le j\le13$ this
curve is certified from end to end and $M$ is monotone along it. So **branch $j$ has exactly one Hopf point for
every integer $M\ge M_j$ and none for $M<M_j$**. At its first integer the bifurcation happens at a very large
rate ($K\approx1.5\times10^4$ for $j=4$ at $M=142$, $K\approx4.5\times10^5$ for $j=13$ at $M=1135$).

**2. Criticality, $j=4,\dots,13$, every $M\ge M_j$.** At every $\kappa_j(M)$ the following are certified:
* $\pm i\omega$ is a simple pair and the crossing is transversal;
* no other eigenvalue lies on $i\mathbb R$ (hypothesis H2, hence no resonance);
* exactly $2j-2$ zeros of the characteristic function give eigenvalues with $\operatorname{Re}\lambda>0$;
* $\ell_1<0$, **except at $(j,M)=(4,142)$ and $(13,1135)$, where $\ell_1=+5.896\times10^{-15}$ and
  $\ell_1=+2.316\times10^{-20}$** (normalization $q^{\mathsf H}q=1$).

So on the centre manifold the bifurcation is supercritical, and subcritical at the two onsets $M=142$ (branch 4) and
$M=1135$ (branch 13). In the full mass hyperplane the bifurcating periodic orbits are **unstable**, because the
equilibrium already has unstable eigenvalues there.

**3. Every branch, large $M$.** The all-M $\varepsilon$-regular formula for $\ell_1(q_1=1)/M$ was evaluated
symbolically at $\varepsilon=0$. It equals result 8's closed form
$L=-R(R-1)(7R^2-34R+23)/[t(R+1)(5R-3)((R-1)^3+R+1)]$, $R=\sqrt{1+t^2}$, $t=t_j$, **identically** in a rational
parametrization of the Hopf curve, and $L<0$. This identity was only conjectured before (formal derivation in
result 8, 53 digits in all-M for $j\le3$). Hence for **every** $j\ge1$: $\lim_{M\to\infty}\ell_1(q_1=1)/M=L(t_j)<0$, so
$\ell_1<0$ for all $M\ge\widetilde M_j$ ($\widetilde M_j$ not explicit). With Pego–Velázquez's Theorem 3.2 (H2 and
$j-1$ unstable pairs for $M>M_{0,k}$, not explicit), the bifurcation on every branch $j\ge2$ is supercritical on the
centre manifold with unstable orbits for all large $M$.

**4. Why the onsets can be subcritical (observation, not a theorem).** Near the onset,
$\widehat L\kappa^3\approx c_j+a_j\kappa$ with $c_j\approx-5$ and a small $a_j>0$, and $M-M^*_j\approx C_j/\kappa$. So
$\ell_1$ changes sign once, at $M^*_j+w_j$, where $w_j\approx0.18/j$ ($w_4\approx0.046$, $w_7\approx0.027$). An
integer onset is subcritical exactly when $M_j-M^*_j<w_j$. This happens at $M=15$ ($j=1$, $M_1-M^*_1=0.107$),
$M=142$ ($j=4$, $0.018$) and $M=1135$ ($j=13$, $0.005$), all three certified, and at no other onset with $j\le13$.
Among $j\le40$ the rule predicts one further case, $j=40$ ($M_{40}-M^*_{40}\approx0.0004$), not computed. Since
$\sum_jw_j$ diverges, infinitely many subcritical onsets are to be expected if the fractional parts of $M^*_j$ are
spread out; this is a heuristic only.

| $j$ | $M^*_j$ (approx.) | $M_j$ | $\ell_1(M_j)$ (unit norm, approx.) | $\varepsilon_c$ | $M_c$ | boxes / H2 groups | integer points |
|---|---|---|---|---|---|---|---|
| 4 | 141.9817605587 | 142 | +5.896225e-15 | 0.06610 | 229 | 1028 / 252 | 142–228 (87) |
| 5 | 209.6289461098 | 210 | -5.627984e-12 | 0.05456 | 336 | 1032 / 268 | 210–335 (126) |
| 6 | 288.7385783342 | 289 | -4.811125e-13 | 0.04653 | 462 | 1084 / 272 | 289–461 (173) |
| 7 | 378.8658344869 | 379 | -1.821088e-14 | 0.04059 | 607 | 1180 / 296 | 379–606 (228) |
| 8 | 479.6599235850 | 480 | -1.233926e-13 | 0.03609 | 768 | 1225 / 319 | 480–767 (288) |
| 9 | 590.8342776983 | 591 | -5.413867e-15 | 0.03252 | 946 | 1277 / 341 | 591–945 (355) |
| 10 | 712.1483888900 | 713 | -3.661047e-13 | 0.02962 | 1140 | 1340 / 362 | 713–1139 (427) |
| 11 | 843.3961316756 | 844 | -6.227954e-14 | 0.02721 | 1351 | 1406 / 375 | 844–1350 (507) |
| 12 | 984.3978885172 | 985 | -3.164222e-14 | 0.02519 | 1576 | 1479 / 399 | 985–1575 (591) |
| 13 | 1134.9950223691 | 1135 | +2.315987e-20 | 0.02346 | 1817 | 1557 / 449 | 1135–1816 (682) |

The decimal onsets and coefficients in this table are approximations to certified enclosures, not interval
endpoints. Full onset enclosures: `certificates/logs/kinf.log` (K4); coefficient signs are re-proved from the
stored Hopf boxes by `pointwise.py --verify`. Here $M_c$ is the first integer covered by the $\varepsilon$-chart.

## How (in short)

For each branch $j$ the proof covers three charts:

1. **$\varepsilon$-chart.** It covers every $M\ge M_c\approx1.6M^*_j$, including $M=\infty$, with all-M's
   $\varepsilon$-boxes (`l1_sweep.py`, `verify_sweep.py`). The spectral certificate (all-M Theorem 4.5) is run on
   merged boxes (`h2cover.py`). For $j\ge5$ it needs narrower contours near $\varepsilon=0$.
2. **Onset chart.** It is regular at $K=\infty$ (`onset_sweep.py`), with unknowns $s=1/(1+K)$, $M$, $\varphi$ and
   the exact identity of $F/A^3$ as a polynomial in $s$ (`check_onset_formula.py`). It follows the curve from
   $K=\infty$ to the $\varepsilon$-chart.
3. **Integer points.** At every integer $M_j\le M<M_c$ (`pointwise.py`) the checks combine:
   * result 8's Krawczyk test and closed form of $\ell_1$;
   * all-M's Krawczyk test and Theorem 4.5 at $\varepsilon=M^{-1/2}$;
   * a link to the onset curve (`link_check.py`).

Result 8's verifier's sparse code (`xcheck_sparse.py`) reproduces the positive signs at $M=142$ and $M=1135$. The two
implementations of $\ell_1$ agree on the overlap of the charts (`check_overlap.py`).

## Reproduce

From `scripts/` (WSL Fedora, Python 3.12, python-flint 0.9, sympy, mpmath; one thread each). Copy `run_branch.sh`
before running it, since `sh` reads scripts while executing them.

| script | claims | runtime |
|---|---|---|
| `kinf.py 1150 13` | K1–K4: limit matrix $C(M)$, $\det(\nu-C)=-G(1+\nu)/\nu$, count $u(M)$ for $3\le M\le1150$, onsets $M^*_j$ | ~1 h |
| `check_onset_formula.py` | O1–O4: the onset-chart identities | 1 min |
| `limit_symbolic.py` | S1–S4: $\lim\ell_1(q_1=1)/M=L(t_j)<0$ for every $j$ | 1 min |
| `run_branch.sh J EPS_C MJ MC M_HI R0 Y cL cR cL2 cR2` | all certificates of branch $J$ (the parameters of each branch: `certificates/logs/run_branch_params.txt`) | 1–3 h |
| `l1_sweep.py`, `verify_sweep.py`, `h2cover.py` | $\varepsilon$-chart: Hopf box, nondegeneracy, $\widehat L<0$; re-proof; H2 and unstable count | |
| `onset_sweep.py J M_HI OUT L1` | onset chart and both junctions | 1–10 min |
| `pointwise.py`, `pointwise.py --verify`, `link_check.py` | integer points; re-proof; link to the curve | 3–10 s per $M$ |
| `check_overlap.py L1 MC 3MC` | result 8's $\ell_1$ inside the $\varepsilon$-chart enclosures | 0.1 s per $M$ |
| `xcheck_sparse.py PW M...` | result 8's verifier's sparse code at selected $M$ (numerical) | seconds |
| `summary.py` | the table above, from the certificate files | 10 s |

Every script prints PASS/FAIL per claim and exits non-zero on failure. Logs: `certificates/logs/`.

## Files

* `result.md`: statements, proofs, certificate description, what is not proved.
* `paper.tex`, `references.bib`, `paper.pdf`: manuscript.
* `references.md`: sources and the novelty check.
* `scripts/`: the scripts above. Copied unchanged: `ad.py`, `allm_core.py`, `l1reg.py`, `h2_lib.py`, `h2j_lib.py`,
  `l1_sweep.py`, `verify_sweep.py`, `h2j_sweep.py` (all-M), `bdlib.py` (result 8).
* `certificates/`: `l1_jJ.csv`, `h2_jJ.csv`, `onset_jJ.csv`, `pw_jJ.csv` per branch; `logs/`.
