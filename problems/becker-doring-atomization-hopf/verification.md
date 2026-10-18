Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.

# Verification: criticality of the Hopf bifurcations in the Becker–Döring equations with atomization

Date: 2026-10-18. I did not read the solver's `PROGRESS.md`, its logs or any `travail/` or `rapports/` folder. I used
only the problem statement, the primary source (downloaded myself), the deliverables (`paper.tex`, `README.md`,
`scripts/`, `certificates/`) and my own code in `verification/`.

## Verdict: CONFIRMED WITH CORRECTIONS

**Confirmed.** Theorem A is correct as stated. For every \((M,j)\) in \(15\le M\le10^5\) (\(j=1\)),
\(44\le M\le10^5\) (\(j=2\)) and \(87\le M\le10^5\) (\(j=3\)), the certificate encloses a Hopf point with a simple pair
\(\pm i\omega\), \(\omega>0\), and a transversal crossing. The sign of the first Lyapunov coefficient of the system
restricted to the mass hyperplane is certified: \(\ell_1<0\) everywhere, except at \((15,1)\), where \(\ell_1>0\).

- **Independent re-proof, \(M\le300\).** I proved the sign again with a different code, different coordinates and a
  different Krawczyk system. It covers every integer \(M\le200\) and every 10th \(M\) up to 300, on all three
  branches. At all 487 points my \(\ell_1\) is the same double as the solver's stored value.
- **\(M=15\).** The subcritical sign was confirmed by two methods that do not use the Lyapunov formula.
- **Large \(M\).** The certificate was cross-checked by an independent \(O(M)\) code.
- **Limit problem.** The algebra of the limit problem (§7) is correct.
- **Novelty.** I found no prior answer.

**Corrected.**

1. *Stability for \(16\le M\le10^5\).* The README stated that "the first Hopf bifurcation is supercritical, with
   orbitally stable periodic solutions, for every \(16\le M\le10^5\)". This was not rigorously justified. The spectral
   hypothesis H2 (no other eigenvalue on \(i\mathbb R\)) was certified only for \(M\le120\). It was checked in
   floating point only for \(M\le300\) and every 50th \(M\le2000\), and not at all for \(2000<M\le10^5\). The paper's
   \(M_{0,k}\) is not explicit and may exceed \(10^5\).

   I extended the rigorous H2 certificate to \(M\le160\), plus \(M=250\) and \(M=300\), on all branches. The statement
   now reads:
   - *rigorous* for \(16\le M\le160\);
   - *conditional on H2* for \(160<M\le10^5\). The sign of \(\ell_1\) itself is certified on the whole range.
2. *The word "first".* The phrase "first Hopf bifurcation" and the sentence "this also proves rigorously that the point
   is the \(j\)-th crossing" overstated what the eigenvalue count shows. I proved, by exact algebra, that for
   \(15\le M\le57\) the point \(\kappa_1(M)\) is the first loss of stability for \(K>0\). I also proved the onsets of
   branches 1 and 2 (\(M=15\) and \(M=44\)) exactly. The wording was corrected accordingly.
3. *A gap in the paper's Lemma 5.4.* The first step uses a convexity argument that is valid only for real
   \(\varphi\). For the certified complex \(\varphi\) it is replaced by a direct check, \(f(\varphi)\ne0\). This check
   is implicit in the certificate and I verified it explicitly on a sample.
4. *Minor overstatements.* Two were reworded:
   - "H2 observed numerically everywhere in \(\mathcal S\)";
   - "the discrete coefficient converges to the limit value".

No error was found in any certified sign, enclosure or formula.

## 1. Primary source and the question

- I downloaded arXiv:1905.02605v2 (PDF and LaTeX source) myself.
- **The open sentence.** It sits in §3.3, right after Theorem 3.3 (source lines 658–660), and the manuscript quotes it
  verbatim: "We have not managed to determine analytically whether the bifurcating solutions are stable (the
  supercritical case) or not. Many of our numerical computations, as in Fig. 1, are consistent with the presence of
  stable periodic solutions, however."
- **The model, checked line by line against eqs. (1)–(4) of the paper.**
  - Fluxes: \(J_\ell=n_\ell n_1-n_{\ell+1}\), \(1\le\ell\le N=M-1\).
  - Equations: \(\dot n_\ell=J_{\ell-1}-J_\ell\) for \(2\le\ell\le N\); \(\dot n_M=J_{M-1}-Kn_M\);
    \(\dot n_1=-J_1-\sum_{\ell=1}^NJ_\ell+MKn_M\).
  - Parameter and equilibrium: \(\kappa=K\sqrt M\) (§3.3 of the paper), constant equilibrium \(A=1+K\).
  - Linearization: the matrix \(B\) of the paper's §3.2.
  - The solver's model (paper.tex §2, `bdlib.field`, `bdlib.Bmat`) is exactly this one. I re-derived by hand the
    exact quadratic identity \(f(A\mathbb 1+v)=Bv+v_1Sv\) and the mass conservation \(\boldsymbol\ell S=0\).
- **Does the result answer the question?** The paper asks whether the bifurcating solutions at \(\kappa_j(M)\),
  \(M>M_{0,k}\), are stable. The result answers it for the finite ranges of Theorem A, not for all \(M\). The answer
  is "supercritical, with one exception at \(M=15\)". Since \(M_{0,k}\) is not explicit, the infinite range is left
  open, as the deliverables honestly state. "Partial result" is the correct status.

## 2. Novelty (re-done independently, 2026-10-18)

- **Citing works.** I fetched Semantic Scholar (16 citing works) and OpenAlex (17). Together they give the same 20
  works as `references.md`.
  - The closest is Niethammer–Pego–Schlichting–Velázquez (SIAM J. Appl. Math. 2022). It concerns a different model
    and gives only a formal expectation of supercriticality.
  - Fortin–Choi (J. Phys. A 2023) and Budzinskiy–Matveev–Krapivsky (PRE 2021) treat addition–shattering models. I
    checked their abstracts: they do not treat this system.
- **Later work of the authors** (OpenAlex author records, 2020–2026).
  - Pego: 17 works. Velázquez: 88 works, including arXiv:2605.15926, which is unrelated.
  - Niethammer: 25. Bonacini: 13.
  - None treats the stability of these orbits.
- **Recent arXiv papers.** The most recent arXiv papers on Becker–Döring were also checked: 2508.07884, 2609.29842
  and 2609.29891. All are unrelated.
- **AI campaigns.**
  - GitHub `mmislan`: 23 repositories (Palomar). None is relevant.
  - GitHub `ipitchford`: 166 repositories (Evidence Press). I scanned the full list of names and descriptions; none
    concerns Becker–Döring, atomization or nucleation.
  - Zenodo (phrase queries "Becker-Döring" AND atomization / Hopf / oscillations): 0 hits.
  - `gh search repos`: only old simulators.
- **Web search.** Two web searches found nothing.
- **Limits.** Google Scholar was not queried. The typeset (paywalled) version of the paper was not seen.
- **Conclusion.** No prior resolution was found, so the result is not "already solved".

## 3. What I recomputed myself (scripts in `verification/`, logs next to them)

All computations ran in WSL Fedora with python-flint 0.9 (Arb), numpy/scipy, sympy, mpmath and PARI/GP. They used
at most 2 threads (`OMP_NUM_THREADS=1`).

### 3.1 Independent dense certificate, \(M\le300\) (`v02_dense_certify.py`)

This is a different method, written without `bdlib.py`.

- **Coordinates.** Reduced coordinates on the mass hyperplane: \(y=(v_2,\dots,v_M)\), \(v_1=-\sum_{\ell\ge2}\ell v_\ell\).
  The reduced Jacobian \(\hat A\) is invertible, so the solver's subtle choice "\(w_{11}\in H\)" never arises. This
  tests that choice independently.
- **Jacobian and quadratic part, from the vector field itself.**
  - Jacobian: exact integer central differences, \(\hat A(K)=A_0+KA_1\), with affinity in \(K\) checked.
  - Bilinear form: polarisation \(f(A+X+Y)-f(A+X)-f(A+Y)+f(A)\).
  - Trilinear term: computed as a third difference, and certified to be 0 for \(M\le40\).
- **Hopf point.** The Krawczyk operator is applied to the full eigen-system
  \((\hat A(K)-i\omega)q=0\), \(q_k=1\), in the \(2M\) real unknowns \((\operatorname{Re}q,\operatorname{Im}q,\omega,K)\),
  at 192 bits with relative radius \(10^{-36}\). The solver's \(\kappa\) served only as a starting guess.
- **\(\ell_1\).** Kuznetsov's formula (3.20) in ball arithmetic. The adjoint vector comes from a bordered system with
  a rigorous solve.
- **Checks made at each point.**
  - Simplicity: the bordered matrix is nonsingular and \(uq\ne0\).
  - Transversality: \(\operatorname{Re}(u\hat A_1q)>0\).
  - \(\omega>0\).
  - The sign of \(\ell_1\).
  - A floating-point count of unstable pairs at \(K_c(1\mp10^{-6})\), equal to \((j-1,j)\).
- **H2 and the unstable count, rigorously.** `acb_mat.eig` on the interval matrix \(\hat A(K\text{-box})\) must give
  exactly 2 eigenvalues meeting \(i\mathbb R\) (at \(\pm i\omega\)), \(2(j-1)\) with positive real part and the rest
  with negative real part.

Results (logs `v02_dense_j*.log`, `v02_dense_210_300.log`; CSVs `v02_dense_*.csv`):

| range | points | result |
|---|---|---|
| \(j=1\), every \(M\in[15,200]\) | 186 | all PASS; H2 certified for \(M\le160\) |
| \(j=2\), every \(M\in[44,200]\) | 157 | all PASS; H2 certified for \(M\le160\) |
| \(j=3\), every \(M\in[87,200]\) | 114 | all PASS; H2 certified for \(M\le160\) |
| \(M=210,220,\dots,300\), each branch | 30 | all PASS; H2 certified at \(M=250,300\) |

- **Agreement.** At all 487 points, my \(\ell_1\), rounded to a double, equals the solver's stored double. The
  relative difference is 0, at the resolution of the float64 storage of the CSV; the built-in check threshold was
  \(10^{-9}\). \(\kappa\) agrees to \(10^{-15}\), again the storage limit.
- **\(M=15\).** \(\ell_1=+1.1951866439\times10^{-9}\pm9\times10^{-45}\) (unit normalisation). At \(M=16\),
  \(\ell_1=-1.2963\times10^{-5}\).
- **Table of §6.** The \(M\le100\) rows of the table agree with my values.

### 3.2 Exact algebra: all imaginary-axis crossings for small \(M\) (`v07_first_crossing.py`, PARI/GP)

- **Method.**
  - \(P(\lambda,K)=\det(\lambda-\hat A(K))\in\mathbb Z[K][\lambda]\), with \(P(iy)=P_e(y^2)+iyP_o(y^2)\).
  - Every \(K>0\) with a nonzero imaginary eigenvalue is a root of \(\mathrm{Res}_u(P_e,P_o)\in\mathbb Z[K]\). The
    formal degrees are used; \(P\) is monic, so there is no degree drop.
  - Its positive real roots are isolated rigorously with `polrootsreal`.
  - \(\det\hat A(K)\) has no positive root.
  - The number of unstable eigenvalues is therefore constant between consecutive roots. It is computed exactly by a
    Routh array over \(\mathbb Q\) at a rational interior point.
- **Results** (logs `v07_first_crossing.log`, `v07_first_crossing_21_60.log`):
  - \(3\le M\le14\): the resultant has no positive root, and the equilibrium is stable on \(H\) for **every** \(K>0\).
    The solver's floating-point "no crossing for \(K\le10^5\)" becomes a theorem for all \(K>0\).
  - \(15\le M\le43\): exactly one crossing, at the certified \(K_1(M)\) (e.g. \(212.8952163315\) at \(M=15\)). The
    count is 0 on \((0,K_1)\) and 2 on \((K_1,\infty)\).
  - \(44\le M\le57\): exactly two crossings, at \(K_1\) and \(K_2\) (\(K_2=232.163\) at \(M=44\)). The counts are
    \(0,2,4\).
  - So, for \(15\le M\le57\), \(\kappa_1(M)\) is the first and, for \(M\le43\), the only loss of stability. The onsets
    \(M=15\) and \(M=44\) are exact.
  - \(M\ge58\) exceeded the 400 MB PARI stack. The run was stopped and the results for \(M=58,59\) are discarded.

### 3.3 Direction of bifurcation without the Lyapunov formula

These tests check the sign convention of Kuznetsov's formula end to end. The convention is shared by the solver's
closed form and by my dense code.

- **Shooting for the small periodic orbits** (`v09_po_shooting.py`). Newton on
  \((x_0,T,K)\) with variational equations: periodicity, mass level of \(A\mathbb 1\), phase, and amplitude \(a\) along
  \(\operatorname{Re}q\).
  - \(M=15\): \(K(a)<K_c\) for \(a=0.1,0.2,0.4\), with \((K-K_c)/a^2=-3.7209\times10^{-5}\) and
    \(-3.7195\times10^{-5}\) at \(a=0.2,0.4\). This is **subcritical**. At \(a=0.1\) the value \(-6.7\times10^{-5}\) is
    at the resolution limit of double-precision shooting; the side is the same.
  - \(M=16\): \(K(a)>K_c\) with \((K-K_c)/a^2=+4.06\times10^{-3}\), constant to 0.03 %. This is **supercritical**.
- **Direct simulation of the full ODE** (`v08_simulate.py`, LSODA with analytic Jacobian, rtol \(10^{-10}\)).
  - For \(K=K_c(1+\delta)\), \(\delta=4,2,1\%\), perturbations starting at \(0.3\times\) and \(3\times\) the predicted
    amplitude converge to the same cycle.
  - Its \(n_1\) half-amplitude over \(2|q_1|\sqrt{-\mu/(\omega\ell_1)}\) is 1.0000, 0.9998, 0.9999 at \(M=16\) and
    1.0151, 1.0076, 1.0038 at \(M=50\) (deviation \(\propto\delta\)).
  - For \(K=K_c(0.98)\) the same perturbation decays to the equilibrium.
  - (The first run's pass criterion "deviation decreases with \(\delta\)" failed at \(M=16\) only because all
    deviations are \(\le2\times10^{-4}\). The criterion was corrected and the run repeated. Both logs are kept:
    `v08_simulate_run1_oldcriterion.log` and `v08_simulate.log`.)

### 3.4 Large \(M\): an independent \(O(M)\) implementation

- **Structure.** In the reduced coordinates, \(\hat A=T-A(e_2+e_M)\boldsymbol\ell^{\mathsf T}\) with \(T\)
  tridiagonal, so every solve is a banded solve plus Sherman–Morrison.
- **Hopf point.** Bordered Newton for the eigenpair and secant in \(K\), started from \(\kappa\) perturbed by
  \(10^{-4}\).
- **Bilinear form.** Taken from \(D^2J_\ell[X,Y]=X_\ell Y_1+X_1Y_\ell\); checked against polarisation.
- `v03_sparse_float.py` (double precision): 84 points, seed 7. It covers 25 random \(M\le10^5\) per branch, plus
  \((15,1)\), \((16,1)\), \((44,2)\), \((87,3)\), \((10^3,1)\), \((10^4,2)\) and \(M=10^5\) on all three branches.
  - The sign of \(\ell_1\) agrees with the certificate at every point.
  - The relative difference is \(\le1.8\times10^{-7}\) for \(M\le10^4\) and \(\le9\times10^{-5}\) up to \(10^5\). This
    is the expected loss of double precision: \(\operatorname{Re}\lambda=O(M^{-3/2})\), and the inner products
    cancel.
  - \(\kappa\) is re-located to \(\le1.6\times10^{-9}\) relative.
- `v04_sparse_arb.py` (Arb midpoints, 320 bits, same algorithm), at 12 points:
  - \((15,1)\), \((16,1)\), \((44,2)\), \((87,3)\);
  - \((10^3,1)\), \((10^4,3)\);
  - \((31416,1)\), \((54321,2)\), \((77777,3)\);
  - \(M=10^5\) for \(j=1,2,3\).

  **Result.** \(\kappa\) agrees with the certificate to \(\le2.4\times10^{-40}\). \(\ell_1\) agrees in both
  normalisations to \(\le9\times10^{-17}\), the float64 storage of the CSV. Log: `v04_sparse_arb.log`.

  **A first run at 160 bits** (log `v04_sparse_arb_160bits.log`) agreed to \(10^{-17}\) at 9 points. At three
  \(j\ge2\) points with \(M\ge77777\), it differed by only \(4\)–\(8\times10^{-10}\). I traced this to my own code: a
  midpoint iteration with near-singular shifted solves loses digits. I re-ran the solver's closed form at 192, 256
  and 384 bits (\(j=3\), \(M=77777\) and \(10^5\)), and it gives identical enclosures. My code at 320 bits then agrees to \(10^{-17}\). So the solver's
  enclosures were right, and the 160-bit discrepancy came from my code.
- **Near-degenerate onsets.** At the onset of each branch \(\ell_1\) is tiny:
  - \(+1.2\times10^{-9}\) at \((15,1)\);
  - \(-3.1\times10^{-9}\) at \((44,2)\);
  - \(-1.2\times10^{-9}\) at \((87,3)\) (unit normalisation).

  The onset values are close to degenerate (Bautin-like) points, but all three signs are certified with radii near
  \(10^{-45}\) by two independent ball-arithmetic codes.

### 3.5 The solver's certificate files

- **`v01_audit_csv.py`** (all 300 000 lines). It checks:
  - coverage of every \(M\), and a single NO_CROSSING line at \(M=14,43,86\);
  - \(\omega_{\rm lo}>0\) everywhere. This matters: Kuznetsov's formula changes sign with \(\omega\), and `04` does
    not test it;
  - \(\operatorname{Im}z>0\);
  - the same sign in both normalisations, the sign pattern, simple = 1 and transv > 0;
  - \(\operatorname{Im}z\) inside the branch window;
  - \(\kappa_j(M)\) strictly decreasing and smooth (no jump to another branch).

  All PASS.
- **`v05_recheck_sample.py`.** It re-proves 1360 lines with the solver's library from the stored midpoints: every
  \(M\le300\), 200 random \(M\in(300,10^5]\) per branch (seed 20261018), and the boundaries 14/15/16, 43/44, 86/87 and
  \(10^5\). Beyond `04`, it also checks \(\omega>0\), \(\operatorname{Im}z>0\) on the box and \(f(\varphi)\ne0\) or
  \(f(\psi)\ne0\). All PASS, in 14 s.

### 3.6 Limit problem (`v06_limit_algebra.py` and by hand)

- **By hand.** I re-derived the formal reduction to (R):
  - the boundary-layer excess mass \(M^{3/2}a\Gamma/(\kappa\gamma)\), from \(J\) constant, \(n_\ell=a+C(1+c)^\ell\)
    and \(C(1+c)^M=a(\gamma-\kappa)/\kappa\);
  - the bulk transport \(V_\tau+\gamma V_x=0\) with \(V(0)=\Gamma\);
  - the Lagrangian clock;
  - the multiple-scales amplitude equation (\(\sigma\), \(\ell\), \(b\)).
- **By script.** `v06` checks:
  - the closed forms of \(E(it)\), \(E(2it)\), \(E'(it)\), \(D(0)\) against quadrature, to 30 digits;
  - the identity \(\operatorname{Re}(b\overline{E'})=-(7R^2-34R+23)/(R(R+1)(5R-3))\) symbolically;
  - the values \(L_1,\dots,L_5\);
  - the bound \(t_1>3.92\Rightarrow R>(17+8\sqrt2)/7\) in interval arithmetic.

  All PASS. The derivation of (R) remains formal, as the manuscript says.

## 4. The solver's scripts, re-run

All runs used WSL Fedora with at most 2 threads. The log is `verification/rerun_solver_scripts.log`.

| script | what I ran | result | time |
|---|---|---|---|
| `01_check_model.py` | full | M1–M6 PASS: model, Table 1 of the paper, closed form = dense to \(\le6\times10^{-12}\) | 58 s |
| `02_validate_l1.py` | full, with Julia/BifurcationKit | V1–V3 PASS: \(M=15\) SubCritical, \(M=16,17,30\) SuperCritical | 3.5 min |
| `03_certify_range.py` | `1 14 80` | regenerated lines identical to the stored `sweep_j1` lines (midpoints, \(\ell_1\) bounds, status) | 37 s |
| `04_verify_certificates.py` | stride 1 on the three \(M\le10^4\) files (29 857 lines); stride 450 on the three `_b` files (600 lines) | PASS | 4.3 min + 19 s |
| `05_limit_coefficient.py` | full | C1–C4 PASS | 15 s |
| `06_large_M_samples.py` | full | S1, S2 PASS; deviation at \(10^8\): 0.13 %, 0.36 %, 0.67 % | 1 s |
| `07_spectrum_H2.py` | `120` on the three small files | 217 points PASS | 5.2 min |

The solver's own logs report a stride-1 re-check of all \(3\times10^5\) lines by `04`. I did not repeat the
stride-1 pass on the `_b` files (about 55 CPU-min). Instead, `v05` re-proved 600 random `_b` lines and the
boundaries.

`06` reaches \(M=10^8\) in about a second because the closed form costs \(O(1)\) in \(M\). For \(j=2\) the approach
to \(L_2\) is not monotone: the relative deviation is \(3.8\times10^{-3}\) at \(10^4\), \(5.4\times10^{-2}\) at
\(10^5\) and \(3.6\times10^{-3}\) at \(10^8\). The fitted \(M^{-0.44}\) uses only the last five samples. This is
consistent with the manuscript's cautious wording, but the rate is not established.

What the PASS lines actually test:
- **`04`.** It re-proves Krawczyk, the sign of \(\ell_1\) in the *unit* normalisation, simplicity, transversality
  and coverage. It does **not** test \(\omega>0\), \(\operatorname{Im}z\ne0\) or the branch label. I checked the first
  two separately (all lines for \(\omega\), and a sample on the box). The label is supported only by continuation for
  \(M>2000\).
- **`07`.** It checks H2 and the unstable count. Its "\(j\)-th crossing" means "\(j-1\) pairs unstable at
  \(\kappa_j\)".
- **`02` V1.** It tests the generic Kuznetsov routine on a planar system, not the closed form. The closed form is
  tested against the dense formula (`01` M6), and both share the sign convention. That is why the formula-free tests
  of §3.3 above were useful.

## 5. What I only read (not recomputed)

- **The closed-form ansatz** of §4 (`bdlib.Problem`, `Seq`, `solve`). I checked the completeness argument: rows
  \(2..N\) form a 2-dimensional affine family, and row 1 is replaced by the mass condition when \(\mu=0\). I did not
  re-implement the closed form. It is validated by agreement with my dense code (\(M\le300\)) and my \(O(M)\) codes
  (large \(M\)).
- **The paper's proofs** of Prop. 5.3 and Lemma 5.4 (read; gap noted in correction 3) and Lemma 3.1.
- **The BifurcationKit scripts.** `bk_check.jl` was read and re-run through `02` (V3). `bk_po_branch.jl` was read
  but not re-run; my own shooting (`v09`) replaces it.
- **The Krawczyk implementation** in `bdlib` (read). It was re-run through `04`/`v05`.

## 6. Corrections applied to the deliverables

- **`paper.tex`.**
  - Subtitle and `\PaperCredits`: verifier attribution added.
  - Abstract: stability scope (rigorous \(M\le160\), conditional beyond) and the verification.
  - §3: the stability of \(j=1\) orbits is stated under H2.
  - §5: corrected "\(j\)-th crossing"; documented H2 coverage, including that nothing is checked for \(M>2000\); added
    the paragraph on Lemma 5.4.
  - §6: Consequences rewritten (rigorous range, conditional range, the word "first"); onset remark updated.
  - §8: "converges" changed to "appears to converge".
  - §10: proved and not-proved lists updated; added a note that points (a)–(f) were addressed.
  - New §11 "Independent verification".
  - Formatting only: the 62 Markdown-style run-in lists ("- ...", "1. ...") were converted to LaTeX
    `itemize`/`enumerate`, since they were printed as running text. The overlong reproduction command line was also
    fixed.
  - Recompiled: `paper.pdf`, 13 pages.
- **`README.md`.** Attribution, status, scope of the stability claim, not-proved list, files.
- **`scripts/`.** No change. No error was found in the solver's scripts.
- **Root `README.md`** (results table, outside this folder): not modified, as instructed. Its status and verifier
  attribution still need updating.
- **`verification/`.** Added: v01–v09 with their logs and CSVs.

## 7. Remaining reservations

1. **\(M>10^5\) is not covered**, and neither is \(j\ge4\). The limit coefficient is formal. The convergence
   \(\ell_1(q_1=1)/M\to L_j\) is numerical evidence only.
2. **H2 for \(160<M\le10^5\)** (except \(M=250,300\)) is not certified. The orbital stability of the first-branch
   orbits is therefore conditional there. The same holds for the word "first" for \(M>57\), and for the branch label
   for \(M>300\), which rests on continuation.
3. **The large-\(M\) ball-arithmetic certificates** come from the solver's closed-form code alone. My independent
   rigorous re-proof stops at \(M=300\). Above that, my independent checks are floating point or high-precision
   midpoint computations at sample points. The solver's code itself was re-run on 1360 lines, and its logs report
   all \(3\times10^5\) lines.
4. **Novelty** was checked without Google Scholar and without the typeset version of the paper.

## 8. Computing time (2 threads at most)

| task | CPU time |
|---|---|
| `v01` audit (all 300 000 lines) | 10 s |
| `v02` dense certificates (487 points, with `eig` for 343 of them) | 61 min |
| `v03` float \(O(M)\), 84 points | 14 min |
| `v04` 320-bit \(O(M)\), 12 points (plus a 160-bit run, 13 min) | 5 min |
| `v05` re-proof of 1360 certificate lines | 14 s |
| `v06` limit algebra | 1.4 s |
| `v07` exact crossings, \(3\le M\le59\) (the \(M=60\) run was stopped) | 22 min |
| `v08` simulations \(M=16,50\) (two runs) | 25 min |
| `v09` shooting \(M=15,16\) | 20 min |
| solver scripts 01–07 (table in §4) | 16 min |
| total | about 3 h CPU, about 2 h wall-clock on 2 threads |
