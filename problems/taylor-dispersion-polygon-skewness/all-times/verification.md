# Verification

**Verdict: CONFIRMED WITH CORRECTIONS.**

Verifier: **Claude Opus 5.5 (effort xhigh)**, a separate agent with a fresh context, 5 November 2026. I received the original
statement and source, the claimed result and the deliverables, not the solver's conversation. I read no PROGRESS.md and no
`explore/` note before forming my view. No sub-agent was used. This is an agent audit, not human peer review and not a
proof-assistant certificate.

The main theorem is confirmed. For every regular n-gon with n ≥ 4, and for the disk, `K(τ) > 0` for all `τ > 0`. So for every
cross-sectionally uniform release with a symmetric longitudinal profile, the skewness is strictly positive at all times and for
every Pe > 0. The triangle theorem is also confirmed: `K` has exactly one zero `τ*`, it is simple, and `τ*` lies in
`[1.1758, 1.2434] ⊂ [1.17, 1.25]` (time unit inradius²/κ). The finite-element cross-check below puts it at `τ* ≈ 1.209`.

No mathematical error was found that affects a claim. The corrections are listed at the end. They are one bound rounded the wrong
way in the text, a remark whose stated two-sided asymptotics had only one side proved, the missing `Pe > 0`, and an incomplete
trusted-base list. The verifier also gave a precise statement of the cited Grisvard theorem and an elementary replacement for it.

## 1. Statement against the source

* **Quote.** I fetched the arXiv abstract (arXiv API, 1708.05097, identical in v1 to v3) and the v3 PDF. The quoted sentence is
  verbatim: "Interestingly, the equilateral triangle is the only regular polygon with this property, all other polygons possess
  positive skewness at all times, although this cannot cannot be proved on finite times due to the lack of closed form flow
  solutions for such geometries." The published version (Stud. Appl. Math. 141 (2018) 399–417) is paywalled, so I did not re-read
  it. Its weakened abstract is quoted from the solver's reading.
* **Model and skewness.** arXiv v3, Section 2, gives the model `T_t + Pe u T_x = ΔT`. Here `Δu = −2`, `u` vanishes on the wall, and
  `T` has Neumann conditions. The skewness is "the centered, normalized third moment" of the longitudinal distribution. The
  simulations release walkers "uniformly distributed across the cross-section". The manuscript uses the same Aris-moment setting.
  Its `μ₃/μ₂^{3/2}` is the source's skewness, which does not depend on the frame, so working with the mean-free `u = ũ − ⟨ũ⟩` is
  legitimate. A symmetric profile has `κ₃ = 0` and covers the plane source of the source. The theorem holds for every Pe > 0.
  At Pe = 0 the skewness vanishes identically, hence the added `Pe > 0`.
* **Consistency with the source's numbers.** The short-time expansion `μ₃ ≈ Pe³⟨u³⟩τ³ + Pe³Eτ⁴` uses `E = ⟨ũ⟩² − ⟨u²⟩`. It agrees
  with arXiv v3 eq. (10), `M₃ ≈ Pe³(⟨u³⟩τ³ − (⟨ũ²⟩ − 2⟨ũ⟩²)τ⁴)`. For the triangle `E = 9/175`. The long-time formula is arXiv
  eq. (13): `3Pe³⟨ug₂⟩/(2+2Pe²⟨ug₁⟩)^{3/2}`. The triangle value `−0.07076` is arXiv eq. (14). Both are reproduced exactly, see §3.

## 2. Re-derivation of the argument

I re-derived each step by hand or with sympy/mpmath. My checking scripts are in the verifier's scratch directory; they are not
deliverables.

1. **Moments (Prop. `moments`).** Let `ρ₁ = 0`. The Aris hierarchy gives `c₁ = Pe·A`, `c₂ = ρ₂ + 2τ + Pe²S` with
   `S' = ΔS + 2uA`, and `M₃' = 3Pe³⟨uS⟩`. Duhamel and self-adjointness then give `⟨uS(σ)⟩ = 2∫₀^σ⟨A'(σ−s) u A(s)⟩ds`. Integrating
   in σ with `A(0) = 0` gives `K(τ) = 2∫₀^τ⟨uA(s)A(τ−s)⟩ds`. Inserting `A = ∫P_r u` gives the triangle formula (KG), with weight
   `τ−a−b`. Next, `μ₂ = ρ₂ + 2τ + 2Pe²∫⟨uA⟩ > 0`. The long-time `K = 2⟨ug²⟩τ + O(1)` follows from `A = g − P_s g` and the
   spectral gap. The short-time `K = ⟨u³⟩τ³/3 + O(τ^{7/2})` follows from `‖P_a u − u‖ ≤ c₀√a‖∇u‖` with
   `c₀ = sup(1−e^{−x})/√x = 0.638`. Finally `⟨ug₂⟩ = 2⟨ug²⟩` by testing the weak forms. All correct.
2. **Tools (§2).** On the half-sector `T` (diameter `R`), Payne–Weinberger gives `C_P = R/π` and `ν = π²/R²` on the invariant
   mean-free space `V`. The trace identity `‖φ‖_∂² = 2‖φ‖² + 2⟨φ x·∇φ⟩` (apothem 1) gives `τ_tr² = 2C_P² + 2RC_P`. In Lemma
   `harmonic`, `m ≤ 1+R` follows from the barrier `(1−x)(x+R)`, and Rellich gives `X ≤ (t+R)Y`. I re-derived Lemma `cell`
   (the weak form of `g − g_H`) and Lemma `modal`: the energy inequality gives `d‖y‖²/ds ≤ ℓ²/2`, which is stronger than the `ℓ²`
   used, and Minkowski gives `L₁`. In Lemma `compare`, `‖A(s)‖ ≤ min(s‖u‖, ‖g‖)`. I also checked
   `sup(3x²+x³)e^{−x} = 2.823 < 3` and `0 ≤ E_kl ≤ 1/max(μ_k,μ_l)`. All correct.
3. **Short times, uniform in n (§3).** Lemma `layer`: `Γ_s` solves the heat equation with inward Neumann flux `m ≥ 0`, so `Γ_s ≥ 0`
   and `⟨Γ_s⟩ = 2s`. I re-derived the identity `F = ⟨u²⟩ − 2cs + Λ` and the formula for `Λ'`.
   * Lemma `barrier`: `ζ = Z(d,r)` is a weak supersolution, because `∫∇d·∇ψ ≥ 0` for `ψ ≥ 0` in `H¹₀` (d is concave) and
     `Z_y(0,r) = 0`. Then `Z(0,r) = 2√(r/π)` and `∫₀^s 4m̄√(r/π)dr = 8m̄s^{3/2}/(3√π)`.
   * Prop. `G`: I expanded `⟨ũ(u−2a+Γ_a)(u−2b+Γ_b)⟩ − cF(a+b)` term by term. Its sign uses only `ũ ≥ 0` and `Γ ≥ 0`.
   * Prop. `short`: all six Dirichlet weights (`τ³/3`, `τ⁴/6`, `τ⁵/15`, `16τ^{9/2}/315`, `8τ^{9/2}/63`, `32τ^{11/2}/3465`) were
     recomputed with sympy. The coefficients `8/35` and `128/3465` were re-assembled by hand.
   * The approximant step (`ũ ≤ ũ_L − θ₋`, `ũ_L ≤ θ₊ + m̄d` with `|∇ũ_L|²` subharmonic, and the enclosures of `c` and `‖u‖`) is
     correct.
   * For n ≥ 200 (apothem-1 frame) I checked `Y ≤ x/6 + d ≤ x/6 + Rd`, `c ∈ [¼ − x/12, ¼ + 5x/12]` and `E = 2c² − ‖ũ‖²`.
4. **Third moment (§4).** I re-derived Lemma `Hsmall`: the maximum principle, the Fourier series on the inscribed circle, and
   `nt − π ≤ nt·t²/3`. I also re-derived the Green decomposition of `⟨u³⟩`, the Rellich bound `6(t+R)t⁵√(q/3)`, and the bound
   `(¾ + 5t²/4 + t²/2)` on the cubic terms. I audited `sg_certify.py` (4 ≤ n ≤ 150) and its exact polynomial library
   `polylib.py`. Its decomposition `∫u³ = ∫u_N³ + 3∫φE + 3∫u_N f² + ∫f³` and all bounds are correct.
5. **Finite n (§5, `finite_n.py`, `modal_lib.py`).** I read the code line by line.
   * **Moments.** Cone moments are integrated exactly by Gauss–Legendre: `2m−1 ≥` the degree, and the nodes and weights are arb
     enclosures. I checked `Δ(p^a q^k)`, the stiffness formula `[4ab + 2n(al+bk)]p^{A−1}q^K + n²kl p^{A+n−1}q^{K−2}` and the
     recursion for `Re(w/R)^{nk}`. The particular solution `g_H` (β, γ_k) and the flux `∂_x g_H` are correct.
   * **Wall ranges.** These come from mean-value enclosures on 300 or more pieces. `sup|u_L|` comes from a grid of `T` plus the
     Lipschitz bound; I checked the covering radius.
   * **Residuals.** The interior residual `‖R_k − ⟨R_k⟩‖` and the flux residual are exact quadratures.
   * **Error budget.** `Δ₀ = 2ω + ρ + L₁` and `ℰ(τ)` match Lemma `compare`. The Taylor step on `[τ_j, τ_{j+1}]` uses negative
     parts of `K̃'` and `K̃''` and the bound `|K̃'''| ≤ 6τ_j^{−2}Σ|W|`. It is compared with `ℰ(τ_{j+1})`, and `ℰ` is increasing.
   * **Tail.** The `τ ≥ 40` tail uses `2W₀₀ − e₁ > 0` and `(2W₀₀ − e₁)·40 > C₁`.
   * **Eigenvalues and multiplicities.** The method never needs eigenvalue enclosures. The Ritz pairs are only exact dyadic trial
     data `(ψ_k, μ_k, d_k)`. Lemma `modal` holds for any such data, and every mode left out, every error in a mode, and every
     multiplicity end up in `ρ` (projection residual) and `L₁` (residual dual norms), which are enclosed rigorously. A repeated
     μ would make `E_kl` divide by zero, which gives a FAIL, never a false PASS. So the question "are all eigenvalues below the
     cut-off accounted for?" does not arise.
   * **Floating point** is used only to choose coefficients, eigenvectors, `d_k` and `μ_k` (all then exact binary numbers) and to
     pick the smaller of two valid bounds. No floating-point comparison decides a claim.
6. **n ≥ 200 and the disk (§6, `uniform.py`).**
   * **Lemma `disk`.** All identities were checked: `⟨g₀²⟩ = 1/11520` and `⟨ūg₀²⟩ = 1/184320` exactly (sympy). The orthonormality,
     `⟨ūψ_k⟩ = −2/μ_k`, `⟨g₀ψ_k⟩ = −2/μ_k²`, `⟨ūψ_k²⟩ = 1/12`, `−2(μ_k+μ_l)/(μ_k−μ_l)²` and `(μ²+96μ−1728)/(48μ³)` were checked for
     k, l ≤ 6 by mpmath quadrature at 30 digits, with largest deviation `4·10⁻³¹`.
   * **Lemma `largeerrors`.**
     * (1) is Lemma `Hsmall` rescaled by `a²`.
     * (2): I recomputed `‖e₀‖_∂² = a x² Q/32` from `nta² = π`.
     * (3) and (4): the radial ψ_k solve the Helmholtz equation exactly, so `R_k − ⟨R_k⟩ = 0`. Also
       `|∂_nψ_k| = 2a|Ψ'(r²)| ≤ 16Ψ_k|e₀|` because `Ψ'(1) = 0`.
     * (5) uses `|𝒫∖D| = |D∖𝒫| ≤ π(1−a²)`.
     * (6): the radial laws agree on `[0,a]`, `sup|F_𝒫 − F_D| ≤ 1 − a²` and `v − a ≤ x/2`. I also checked
       `b̄ ≥ 0` (from `arctan t ≥ 3t/(3+t²)`, whose difference has derivative `4t⁴/((1+t²)(3+t²)²) ≥ 0`) and
       `b̄ ≤ 2x²/45 + x³/30 ≤ x²/20` for x ≤ 1/6.
   * **The `r_W` term.** `⟨w⟩_𝒫 = 0` exactly, so the cross term `m_a m_b⟨w⟩` vanishes, and the bound
     `w₁ Lip(ūf_af_b) + (b̄/2)F² + 2w₁LUF` is complete.
   * **Uniformity in n.** The script does not rely on monotonicity claims it cannot check. It evaluates every bound at
     `t = tan(π/200)` with `a² ∈ [1 − x/3, 1]`, `ζ ∈ [0, x/5]`, Bessel enclosures on `r ∈ [a_lo, v_up]`, and `[1, v_up]` for `ζ_M`.
     Each of these enclosures contains the corresponding quantity of every n ≥ 200. The one bound with an explicit x-dependence,
     `⟨ū²⟩_𝒫`, is increasing in x. The short-time part uses the interval `x ∈ [0, t²]`.
   * **Time rescaling.** `τ̂₀/a² ≤ 0.015` is checked.
7. **Triangle (§7).** `g` is a polynomial, and the solver's `triangle.py` does not itself check `−Δg = u` (`triangle_exact.py`
   does). I re-checked with sympy that `−Δg = u` identically and `∂_n g = 0` on all three sides. I also recomputed
   `⟨ug⟩ = 9/1540`, `⟨ug²⟩ = −7101/476476000`, `⟨u²⟩ = 27/700`, `⟨u³⟩ = 27/19250`, `E = 9/175`, the long-time coefficient
   `−789√770/309400` and `S^G = 4√21/99`. The bounds `−3/10 ≤ u ≤ 11/30` and `m̄ = max|∂_xũ| = 3/2` on the side hold.
   * I re-derived eq. (Kprime) by differentiating the identity of Lemma `compare` and integrating by parts. Note `e(0) = −Ã(0)`.
   * The logic is: `h(x)` increases and then decreases, which gives `(0, 0.025]`. Then `K > 0` on `[0.025, 1]`, `K' < 0` on
     `[1, 1.6]` and `K < 0` on `[1.6, ∞)`, so the zero is unique and simple. Bisection keeps only certified signs.

**Grisvard, Thm 3.2.1.2.**
* **Use.** It is used for `ũ ∈ H²(P_n)`. That regularity enters Lemma `harmonic` (hence Lemma `cell`, `sg_certify`, `u3_uniform`
  and the earlier stage) and the Green formula with the trace `m` in Lemma `layer`.
* **Correctness.** The theorem says that on a bounded convex open set the Dirichlet problem for the Laplacian with L² data has its
  solution in `H² ∩ H¹₀`. I could not open the book. The statement and theorem number agree with how the theorem is cited in the
  literature (e.g. arXiv:2401.02399), and it is the classical result (Kadlec 1964). `P_n` is convex and the data are constant, so
  the use is correct.
* **Replacement.** For a convex polygon there is an elementary argument, now added to the proof of Lemma `harmonic`. At a vertex
  of opening `α < π`, put `q = ½ϱ²(cos(2θ−α)/cos α − 1)`. For n = 4, use
  `q = −ϱ² sin²θ − (2/π)ϱ²(log ϱ sin 2θ + θ cos 2θ)`. Then `ũ − q` is bounded, harmonic and zero on both sides, so it equals
  `Σ c_k ϱ^{kπ/α} sin(kπθ/α)`, which is in `H²` because `π/α > 1`. Away from the vertices, Schwarz reflection does the rest. I
  checked with sympy that `−Δq = 2` and that q vanishes on both sides.

## 3. Certificates re-run and compared

All scripts were copied to a scratch directory and re-run there under WSL Fedora (Python 3.12.15, python-flint 0.9.0). The
repository's logs and certificates were not overwritten.

| script | result | comparison with the committed log |
|---|---|---|
| `run_finite.py 4 200 0.015 9` (finite part, 197 values of n) | all PASS; worst margin 0.294 at n = 21, τ = 0.015; smallest `S(0.015) = 1.333·10⁻⁴` at n = 4 | `finite_4_200.json`: identical certified quantities (see below) |
| `uniform.py 200 12`, `uniform.py 0 12` | all PASS | identical |
| `u3_uniform.py 40` | PASS (min `F ≥ 1.0629·10⁻³`) | identical |
| `triangle.py 0.025 1.0 1.6` | all PASS, zero in `[1.1758, 1.2434]` | identical except the 11th digit of two printed floating-point values (numpy eigensolver); certified claims identical |
| `../scripts/sg_certify.py 4 150` | all 147 PASS | identical up to timings |
| `../scripts/triangle_exact.py`, `lt_certify.py 4 40`, `lt_uniform.py 30`, `sg_uniform.py 140 6`, `alltime_windows.py`, `table1.py` | all PASS; Table 1 enclosures reproduced (only `S^G_100` of arXiv v3 DIFFERS, as reported) | identical to `../certificates/run_all.log` |

The finite run took 30 minutes on 9 processes, and all 197 values of n PASS both the grid test and the short-time test. I
compared the new `finite_4_200.json` with the committed one, n by n:
* `M`, `ω`, `osc`, `Lip`, `δ`, `U`, `2W₀₀` and `S(0.015)` are identical.
* `Δ₀` agrees to `2·10⁻¹²` absolute and `e₁` to `2·10⁻⁷` relative.
* The worst margin agrees to `7·10⁻⁸`.
* The only visible differences are in `L₁`: at most `5·10⁻⁶` relative for n = 4, 5, and `1.8·10⁻¹³` against `2.2·10⁻¹³` for
  n = 6, a negligible absolute amount. They come from round-off in the floating-point eigenvectors that choose the trial modes. Any such choice is rigorous, and the
  certified conclusions are unchanged.

## 4. Independent numerical cross-check (different method, not a certificate)

I wrote my own floating-point solver. It does not use the solver's code or basis. It uses P1 finite elements on the half-sector
`T` (structured mesh, Dirichlet torsion problem, Neumann elsewhere), the full dense generalized eigendecomposition of the
discrete Neumann Laplacian, the formula `A(s) = Σ(1−e^{−λs})/λ (u,φ)φ`, `K(τ) = 2∫₀^τ⟨uA(s)A(τ−s)⟩ds` on 80 Gauss nodes in
`θ` (`s = τ(1−cos θ)/2`), and Richardson extrapolation over two meshes. For the disk I used a separate radial P1 FEM.

| case | quantity | independent FEM | certified / exact value |
|---|---|---|---|
| square (n = 4) | `K(0.015)` | `8.2000·10⁻¹⁰` | `K̃(τ₀) = 8.1986·10⁻¹⁰`, `ℰ = 1.13·10⁻¹⁰` |
| | `K(0.1)`, `K(1)` | `2.6405·10⁻⁷`, `2.04869·10⁻⁵` | lattice cross-check `2.64057·10⁻⁷`, `2.04871·10⁻⁵` |
| | `⟨ug²⟩`, `μ₁` | `1.25880·10⁻⁵`, `9.8696` | `W₀₀ = 1.25880·10⁻⁵`, `π²` |
| | `c`, `E`, `S^G` | `0.281154`, `0.049134`, `0.081295` | `c ∈ [0.28106, 0.28135]`, `E ≥ 0.048985`, `S^G ∈ [0.081294, 0.081297]` |
| | `κ(R=1)` | `0.138526` | `[0.13849, 0.13856]` |
| n = 21 (worst margin) | `K(0.015)` | `3.1938·10⁻¹⁰` | `K̃(τ₀) = 3.1920·10⁻¹⁰`, `ℰ = 2.22·10⁻¹⁰` |
| | `⟨ug²⟩`, `μ₁`, `⟨u³⟩` | `5.63186·10⁻⁶`, `14.5721`, `1.019·10⁻⁶` | `W₀₀ = 5.63186·10⁻⁶`, `14.5721`, `1.021·10⁻⁶` (`sg_certify`) |
| triangle (n = 3) | `K(1.0)`, `K(1.6)` | `4.8165·10⁻⁶`, `−1.06867·10⁻⁵` | `K̃ = 4.8165·10⁻⁶`, `−1.06867·10⁻⁵` |
| | `K(1.17)`, `K(1.2)`, `K(1.25)` | `+9.70·10⁻⁷`, `+2.27·10⁻⁷`, `−1.046·10⁻⁶` | zero in `[1.1758, 1.2434]` |
| | `⟨ug²⟩`, `E`, `c` | `−1.49032·10⁻⁵`, `0.0514286`, `0.300000` | `−7101/476476000`, `9/175`, `3/10` |
| disk | `K(0.0149)` | `3.058783·10⁻¹⁰` | `K_D` (12 modes) `3.058664·10⁻¹⁰`, 200 modes `3.058783·10⁻¹⁰` |

The FEM values agree with the certified `K̃` to 3 to 5 digits, and the differences are far inside the certified error `ℰ`. The
triangle's sign change is at `τ* ≈ 1.209`. The disk check validates the closed-form `W_kl` of Lemma `disk` independently: 12 modes
differ from the converged series by `1.2·10⁻¹⁴`, well inside the truncation term `ρ ≤ 1.1·10⁻⁶` of the certificate.

## 5. Gaps looked for and not found

* Normalizations: averages over `P` against `T`, and boundary norms with `|P|⁻¹∮` and the Gauss weights. They are consistent.
* The short-time window is uniform in n. `S(0.015) > 0` is certified for every 4 ≤ n ≤ 200 and, uniformly, for all n ≥ 200 and
  the disk. It is tight at n = 4, where `S = 1.33·10⁻⁴` against `E/3 ≈ 0.0164`, but it is certified in ball arithmetic.
* The tails at `τ ≥ 40`, the time rescaling at n ≥ 200, and the overlap at n = 200 are all covered.
* `⟨u³⟩ ≥ 0`, needed in Prop. `short`, is certified: by `sg_certify` for 4 ≤ n ≤ 150 and by `u3_uniform` for n ≥ 40.

## 6. Corrections applied by the verifier

1. **paper.tex, Prop. `large`.** The text said `ω ≤ 4.2·10⁻⁶`, but the certified value is `4.251·10⁻⁶` (log), so the text rounded
   an upper bound down. Corrected to `ω ≤ 4.3·10⁻⁶`. The script always used the exact value.
2. **paper.tex, Remark `E`.** The remark claimed `K = ⅓⟨u³⟩τ³ + ⅓Eτ⁴ + O(τ^{9/2})`, but only the lower bound was proved. I added
   the upper bound, which comes from the exact expansion with `⟨ũ²Γ_b⟩ ≤ max ũ·Λ(b)` and
   `⟨ũΓ_aΓ_b⟩ ≤ (2‖u‖_∞ + 2a)Λ(b)`, and noted that the remark is not used in the proofs.
3. **paper.tex, Lemma `harmonic`.** I added the precise statement of Grisvard's Thm 3.2.1.2 and the elementary vertex-expansion
   replacement.
4. **paper.tex, trusted base.** It now names every use of `H²` regularity (including the Green formula in Lemma `layer`) and the
   bounds `|J₀| ≤ 1`, `|J₁| ≤ 1/√2`.
5. **paper.tex (abstract and Theorem `main`) and README.** "Every Péclet number" becomes "every Pe > 0", since the skewness is
   identically zero at Pe = 0.
6. **Credits and status.** These were updated in `\PaperCredits` and `\PaperSubtitle`, in the Limitations paragraph, and in the
   attribution and status of `README.md` and of `../README.md`. The status is now **solved**. `paper.pdf` was recompiled.

## 7. Remaining caveats

* The proof is computer-assisted. Its trusted base is Arb/python-flint, sympy, Payne–Weinberger (Bebendorf 2003),
  `H²` regularity on convex polygons (Grisvard, or the elementary argument above), and standard semigroup theory.
* The published version of the source and Grisvard's book were not re-read by the verifier (paywalled). The arXiv v3 text was read.
* The earlier-stage scripts `lt_certify.py`, `lt_uniform.py`, `sg_uniform.py`, `alltime_windows.py` and `table1.py` were re-run
  and reproduce their logs, but their code was not audited line by line. The main theorem does not depend on them:
  * the long-time sign `⟨ug²⟩ > 0` is re-proved by the audited all-times certificates (`2W₀₀ − e₁ > 0`);
  * the short-time sign comes from the audited `sg_certify.py` and `u3_uniform.py`;
  * the triangle values were re-derived exactly by the verifier;
  * the Table 1 enclosures are a side result. They agree with my FEM for n = 4.
