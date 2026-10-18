**Attribution.** Solver: Claude Opus 5.5 (effort xhigh), sub-agent of the orchestrator; one nested Claude Opus 5.5 sub-agent (effort xhigh) did the citing-works search. LaTeX edition: Claude Opus 5.5 (effort xhigh). Independent verification and corrections: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning. Verdict: CONFIRMED WITH CORRECTIONS (`verification.md`).

# Becker–Döring equations with atomization: is the Hopf bifurcation supercritical?

## Problem

R. L. Pego and J. J. L. Velázquez, *Temporal oscillations in Becker–Döring equations with atomization*,
Nonlinearity 33 (2020) 1812–1846, doi:10.1088/1361-6544/ab6815, arXiv:1905.02605.

The paper's Theorem 3.3 gives a Hopf bifurcation at \(\kappa=\kappa_j(M)\) (\(\kappa=K\sqrt M\)) for every \(j\le k\)
and \(M>M_{0,k}\). Right after it, §3.3 says (verbatim):

> "We have not managed to determine analytically whether the bifurcating solutions are stable (the supercritical
> case) or not. Many of our numerical computations, as in Fig. 1, are consistent with the presence of stable
> periodic solutions, however."

**Novelty check (2026-10-18; redone independently by the verifier the same day, see `verification.md`).** 20 citing works (OpenAlex and Semantic Scholar) were checked, along with the
authors' later papers and web and arXiv searches. None settles the question; see `references.md`.

## Status: partial result, independently verified (confirmed with corrections)

**Answer.** The bifurcation is **supercritical**: the first Lyapunov coefficient is \(\ell_1<0\). There is one
exception: at \(M=15\), the first value of \(M\) where the first branch exists, it is **subcritical** (\(\ell_1>0\)).
The scope of each claim is as follows.

1. **Theorem A (computer-assisted, ball arithmetic).** For every integer \(M\) with \(15\le M\le10^5\) (\(j=1\)),
   \(44\le M\le10^5\) (\(j=2\)) and \(87\le M\le10^5\) (\(j=3\)):
   - the Hopf point \((\kappa_j(M),\omega)\) is enclosed by the Krawczyk method;
   - \(\pm i\omega\) is certified simple and the crossing transversal;
   - the first Lyapunov coefficient of the system restricted to the mass hyperplane is certified **negative**,
     except for \((M,j)=(15,1)\), where it is certified **positive**.

   One spectral hypothesis of the Hopf theorem (H2) says no other eigenvalue lies on \(i\mathbb R\). It is
   certified with Arb's rigorous eigensolver for \(M\le120\) by the solver, and for \(M\le160\) and
   \(M=250,300\) by the verifier's independent code. Beyond that it is checked only in floating point: for
   \(M\le300\) and every 50th \(M\le2000\), and not at all for \(M>2000\). The paper proves it for
   \(M>M_{0,k}\), but \(M_{0,k}\) is not explicit. Consequently:
   - the Hopf bifurcation at \(\kappa_1(M)\) is supercritical, with orbitally stable periodic solutions (within the
     mass hyperplane), rigorously for every \(16\le M\le160\);
   - for \(160<M\le10^5\) the same holds *conditionally on H2*. The sign \(\ell_1<0\) itself is certified on the
     whole range.
   - (verifier, exact algebra) For \(3\le M\le14\) no eigenvalue ever reaches \(i\mathbb R\) for \(K>0\). For
     \(15\le M\le43\) there is exactly one crossing and for \(44\le M\le57\) exactly two. So \(\kappa_1(M)\) is
     the first loss of stability for \(15\le M\le57\), and the onsets \(M=15\) (branch 1) and \(M=44\)
     (branch 2) are exact.

2. **Limit \(M\to\infty\) (formal derivation + exact algebra).** In a Lagrangian "clock" time the system reduces
   to the scalar renewal equation
   \[G/(\kappa(\kappa+G))+\int_0^1sG(T-s)\,ds=0.\]
   Its linearization is exactly the paper's \(Q(z;\kappa)\). Its Landau coefficient is computed in closed form.
   At the \(j\)-th Hopf point it has the sign of \(-(7R^2-34R+23)\), where \(R=\sqrt{1+t_j^2}=1+(\kappa_j^0)^2\).
   This is negative **for every \(j\)**, since \(R_j\ge4.60>4.045\). The predicted limit of \(\ell_1/M\)
   (normalization \(q_1=1\)) is
   \[L_j=-\frac{R(R-1)(7R^2-34R+23)}{t(R+1)(5R-3)((R-1)^3+R+1)},\]
   with \(L_1=-0.0093114\), \(L_2=-0.0079410\), \(L_3=-0.0037777\).

3. **Convergence evidence.** \(\ell_1\) is certified negative at sample points up to \(M=10^8\). There
   \(\ell_1(q_1=1)/M\) agrees with \(L_j\) to \(0.13\%\), \(0.36\%\) and \(0.67\%\) for \(j=1,2,3\). The deviation
   decays like \(M^{-1/2}\).

**Not proved.**
- The convergence \(\ell_1(q_1=1)/M\to L_j\) for \(M\to\infty\), so \(M>10^5\) is not rigorously covered.
- The derivation of the limit equation (formal).
- The spectral hypothesis H2 for \(160<M\le10^5\), except \(M=250,300\) (floating point up to \(M=2000\), unchecked beyond), hence the orbital stability for \(M>160\).
- The branch label for \(M>300\) (continuation in \(M\)); "first" crossing for \(M>57\).
- \(j\ge4\) (floating point only, \(M\le300\); the limit coefficient is negative for all \(j\)).

The full list of gaps and the points for a verifier are in `paper.pdf` (section "What is proved, what is not, and points for the verifier").

**Validation.**
- **Planar test.** The generic Kuznetsov routine reproduces the Guckenheimer–Holmes coefficient to \(10^{-12}\).
- **Direct simulation of the full system (\(M=30\)).** A small perturbation converges to a small stable cycle. Its
  amplitude matches the normal-form prediction (ratio 1.0047, 1.0024, 1.0012 for \(\delta=2,1,0.5\%\)).
- **BifurcationKit** (independent Julia code, automatic differentiation) classifies \(M=15\) as SubCritical and
  \(M=16\), \(17\) and \(30\) as SuperCritical. Its periodic-orbit continuation confirms the direction: the branch
  leaves towards \(K<K_c\) at \(M=15\) and towards \(K>K_c\) at \(M=16\).
- **The paper's Table 1** is reproduced from certified enclosures for \(M=10^2,\dots,10^6\).

## Reproduce

From the project root (Git Bash):
`outils/wsl "cd chemistry-open-problems/problems/becker-doring-atomization-hopf/scripts && python 01_check_model.py"`.
Run the other scripts the same way. Every script prints PASS/FAIL per claim and exits non-zero on failure.

| script | claims | runtime |
|---|---|---|
| `01_check_model.py` | model identities, paper's eigenvalues and Table 1, closed form = dense formula | 1 min |
| `02_validate_l1.py` | planar test, stiff simulation, BifurcationKit cross-check (needs julia) | 10 min |
| `03_certify_range.py J MMIN MMAX OUT.csv [SEED.csv]` | generates the sweep certificates (Theorem A) | ~15 ms per \(M\) |
| `04_verify_certificates.py STRIDE FILES` | re-proves each line from stored midpoints, checks coverage and signs | ~20 ms per line |
| `05_limit_coefficient.py` | limit equation: exact identities, sign for all \(j\), \(L_j\), simulation | 1 min |
| `06_large_M_samples.py` | certified \(\ell_1\) at \(M=10^4,\dots,10^8\) and convergence to \(L_j\) | 2 min |
| `07_spectrum_H2.py MMAX FILES` | rigorous full spectra (hypothesis H2 and branch label) for \(M\le\)MMAX | 15 min (MMAX = 120) |
| `bk_check.jl M Kc ...` | BifurcationKit normal form (independent code) | 2 min |
| `bk_po_branch.jl` | periodic-orbit branch by collocation at \(M=15\) (towards \(K<K_c\)) and \(M=16\) (towards \(K>K_c\)) | 3 min |

`bdlib.py` is the shared library: model, closed-form ball evaluation of \(\ell_1\), Krawczyk operator.

**Verifier's scripts** (`verification/`, run the same way from that folder; details and logs in `verification.md`):
`v01_audit_csv.py` (audit of all certificate lines), `v02_dense_certify.py` (independent dense ball-arithmetic
certificate in reduced coordinates, \(M\le300\), with H2), `v03_sparse_float.py` / `v04_sparse_arb.py` (independent
\(O(M)\) cross-checks up to \(M=10^5\), double precision / 320 bits), `v05_recheck_sample.py` (re-proof of 1360 lines
with the solver's library plus extra checks), `v06_limit_algebra.py` (limit problem), `v07_first_crossing.py` (exact
location of all imaginary-axis crossings, PARI/GP, \(M\le57\)), `v08_simulate.py` (direct simulation, \(M=16,50\)),
`v09_po_shooting.py` (periodic orbits by shooting, \(M=15,16\)).

## Files

- `paper.tex`, `paper.pdf`: manuscript, full technical write-up (model, Lyapunov coefficient on the mass hyperplane, closed form, certificate,
  limit problem and its Landau coefficient, limitations).
- `references.md`: sources and what was read.
- `PROGRESS.md`: log of the work.
- `verification.md`: independent verification report (verdict, what was recomputed, corrections, reservations).
- `verification/`: the verifier's scripts and logs (`v01`–`v09`, see `verification.md`).
- `scripts/`: see the table above.
- `certificates/`:
  - `sweep_j{1,2,3}.csv.gz` (\(M\le10^4\)) and `sweep_j{1,2,3}_b.csv.gz` (\(10^4<M\le10^5\)): one line per
    \((M,j)\) with the certified Hopf box (midpoints, radii), the \(\omega\) and \(\ell_1\) enclosures and the
    nondegeneracy checks. In the `status` column, `PASS` means \(\ell_1<0\) together with all checks. The single
    `FAIL` line \((M,j)=(15,1)\) records the certified \(\ell_1>0\), i.e. the subcritical case;
    `04_verify_certificates.py` checks exactly this two-sided claim. `NO_CROSSING` marks the first \(M\) of each
    file (14, 43, 86). Every line of the six files was re-proved by `04_verify_certificates.py 1 ...`
    (logs `04_verify_*.log`);
  - `logs/`: outputs of the PASS/FAIL scripts.
