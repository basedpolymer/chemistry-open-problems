**Attribution.** Solver: Claude Opus 5.5 (effort xhigh), sub-agent of the orchestrator. Builds on result 8
(`../README.md`; solver, editor and verifier credited there). Verification of this extension: pending.

# Becker–Döring with atomization: Hopf criticality for every M (extension of result 8)

## Problem

R. L. Pego and J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605. After their Theorem 3.3 (Hopf
bifurcation at $\kappa=\kappa_j(M)$, $\kappa=K\sqrt M$, for $j\le k$ and $M>M_{0,k}$, $M_{0,k}$ not explicit),
§3.3 says (verbatim):

> "We have not managed to determine analytically whether the bifurcating solutions are stable (the supercritical
> case) or not. Many of our numerical computations, as in Fig. 1, are consistent with the presence of stable
> periodic solutions, however."

Result 8 certified $\ell_1<0$ for $16\le M\le10^5$ ($\ell_1>0$ at $M=15$) on branches 1–3, but the spectral
hypothesis H2 only for $M\le160$, and nothing for $M>10^5$. This folder closes both gaps on branches 1, 2 and 3.

## Status

**Branch 1 (the first Hopf bifurcation): solved (preliminary)** — a complete answer for every $M$, not yet
verified. **Branches 2–3: solved (preliminary)** as well, in the sense made precise below
(supercritical on the centre manifold, unstable orbits). **Branches $j\ge4$: open** (only result 8's formal limit).
For $M\ge161$ the proofs are computer-assisted (ball arithmetic, every certificate re-proved by a PASS/FAIL
script); for $M\le160$ they are result 8's certificates and its verifier's. The analytic lemmas, with explicit
constants, are in `result.md`.

## Result

For every integer $M\ge16$ the Hopf bifurcation of the constant equilibrium at $\kappa_1(M)$ (system restricted to
the mass hyperplane) is **supercritical, and the bifurcating periodic orbits are orbitally asymptotically stable**:
at $\kappa_1(M)$ the eigenvalues $\pm i\omega$ are simple, the crossing is transversal, every other eigenvalue has
negative real part, and $\ell_1<0$. At $M=15$ it is subcritical (result 8), and for $3\le M\le14$ there is no Hopf
bifurcation (result 8's verifier). So the answer to the open question for the first branch is: *stable for every
$M\ge16$, unstable (subcritical) for $M=15$*. Moreover $\ell_1(q_1=1)/M\to-0.00931139838076\ldots$ (certified enclosure;
$\ell_1(q_1=1)/M\in[-0.0093670,-0.0092545]$ for every $M\ge10^{10}$), which agrees with result 8's formal closed form
$L_1$ to about 53 significant digits (not a proof of the identity, and not needed for any claim).

For $j=2$ (every $M\ge44$) and $j=3$ (every $M\ge87$): $\pm i\omega$ simple, transversal, $\ell_1<0$ (result 8 for
$M\le10^5$, this work for all $M\ge161$), and H2 (verifier of result 8 for $M\le160$, this work for all $M\ge161$).
So the bifurcation is supercritical on the 2-dimensional centre manifold, but the bifurcating orbits are
**unstable** in the full system, because the equilibrium already has eigenvalues with $\operatorname{Re}\lambda>0$
there (certified: $2(j-1)$ zeros of the characteristic function, counted with multiplicity).

**How (for $M\ge161$).** Everything is written in Pego–Velázquez's scaled variables
$\varepsilon=M^{-1/2}$, $z=M(\varphi-1)$, $\kappa$, as expressions that are *regular at $\varepsilon=0$*:
the characteristic function $Q^\varepsilon(z;\kappa)=K^{-3}F(1+z/M)$ (exponentially small part $A^{-M}$ in a
disk), and — the main new ingredient — the first Lyapunov coefficient, rewritten in transfer-function form
$\ell_1=\omega^{-1}\operatorname{Re}[\dots]$ with resolvents computed from exact exponential-sequence formulas in
which every leading-order cancellation (boundary layers, monomer row, geometric sums) is done analytically. Then
$\varepsilon$ is treated as an interval: 562 boxes cover $\varepsilon\in[0,161^{-1/2}]$ (the box $[0,10^{-5}]$
contains $M=\infty$). On each box: a parametric Krawczyk test encloses the Hopf point for all $\varepsilon$ in the
box; simplicity and transversality; $\widehat L=\ell_1(q_1=1)/M<0$ by a mean value form with automatic
differentiation; and H2 via an explicit localization lemma (an explicit version of PV's Props. 8.1–9.2) plus a
rigorous winding-number count of the zeros of $Q^\varepsilon$ (exactly 4: $z=0$ twice and the Hopf pair). Branches 2, 3:
the same on 589 and 834 boxes, with contours along level curves of $\operatorname{Re}\lambda$ (`h2j_lib`).

## Reproduce

From this folder (`scripts/`), Python 3.13 with python-flint 0.9, sympy, mpmath; one thread.

| script | claims | runtime |
|---|---|---|
| `check_formulas.py` | F1–F7: exact identities (sympy), entrywise dense enclosure test of every block of the regular $\ell_1$ formula at $M=40,100,300$, closed-form sums vs direct summation | 2 min |
| `l1_sweep.py 1 0 0.0788111 OUT --h2` | Theorem 7.1: all $M\ge161$, branch 1 (Hopf box, simplicity, transversality, $\widehat L<0$, H2) | 12 min |
| `l1_sweep.py J 0 0.0788111 OUT` | Theorem 7.2: branches 2, 3, all $M\ge161$ (Hopf box, simplicity, transversality, $\widehat L<0$) | 3 min, 4.5 min |
| `h2j_sweep.py J IN OUT R0 Y cL cR` | Theorem 7.3: H2 and the unstable count at $\kappa_j$, $j=2,3$, on the boxes of an `l1_sweep` file | 10 min, 16 min |
| `run_all.sh` | all of the above in sequence (parameters as used), stops at the first FAIL | about 1 h |
| `verify_sweep.py FILE EPS_LO EPS_HI [--h2] [STRIDE]` | re-proves every certificate line from the stored data, and the coverage | as the sweep |
| `check_vs_result8.py J OURS R8FILES` | result 8's Hopf points ($161\le M\le10^5$) lie in our unique-zero boxes; result 8's $\ell_1$ values inside our enclosures | 1 min |
| `asymptotics.py` | A1–A4: limit $M\to\infty$ on branches 1–3, comparison with result 8's closed form | 10 s |
| `h2_sweep.py OUT 161` | H2 only, branch 1 (tool; superseded by `l1_sweep --h2`) | 1 min |

Every script prints PASS/FAIL per claim and exits non-zero on failure. Logs: `certificates/logs/`.

## Files

- `result.md`: statements, complete proofs with explicit constants, certificate description, what is not proved.
- `references.md`, `NEEDS.md` (optional DOIs), `PROGRESS.md` (log, git-ignored).
- `scripts/`: `allm_core.py` (regular $Q^\varepsilon$, parametric Krawczyk, localization, region, winding number),
  `l1reg.py` (regular $\ell_1$), `ad.py` (automatic differentiation and series), `h2_lib.py`, `h2j_lib.py`, and the
  scripts above.
- `certificates/`: `l1h2_j1.csv` (branch 1, 562 boxes), `l1_j2.csv` (589), `l1_j3.csv` (834), `h2j_j2.csv`,
  `h2j_j3.csv`; `logs/`; `sources/` (git-ignored copy of arXiv:1905.02605).
