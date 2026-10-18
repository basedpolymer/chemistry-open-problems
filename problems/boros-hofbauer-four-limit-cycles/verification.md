# Independent verification — Boros–Hofbauer four-limit-cycle conjecture

Verifier: Claude Opus 5.5 (effort xhigh), fresh context, no access to the solver's reasoning.
(`PROGRESS.md` was not read. No sub-agent was used. Date: 2026-10-18.)

## Verdict: **CONFIRMED WITH CORRECTIONS**

The conjecture of Boros–Hofbauer (JDDE 36 Suppl 1 (2024) 175–197, §5.1) is proved by the deliverables, in the exact
form stated there: for suitable parameter values the unique positive equilibrium of (16) is asymptotically stable and
is surrounded by four limit cycles Γ0, Γ1, Γ2, Γ3 which are unstable, stable, unstable, stable. Both independent
proofs (Theorem A: certified triple weak focus + successive perturbations; Theorem C: validated integration at
explicit parameters) are correct; every computational ingredient that carries weight was re-derived or re-computed by
the verifier with independent code and, where possible, an independent algorithm or library. The corrections are
editorial (one wrong interval in the abstract, stale file references, a weakly supported side remark, status lines);
none affects the mathematics.

## 1. Source and statement (re-read)

arXiv:2103.00972 (LaTeX source, `sections/5_three_reactions.tex`) and the published version (PMC10901957, checked
by the verifier): network with (a1,b1)=(0,0), (a2,b2)=(0,−1), (a3,b3)=(a,b), reaction vectors (0,−1), (1,−1), (−1,d),
a>0, b>−1, d>0, 1+bd>0; positive equilibrium iff d>1; scaled ODE ẋ = 1/y − x^a y^b, ẏ = K(−(d−1) − 1/y + d x^a y^b);
det J = K(d−1)a; K = a/(1+bd) makes tr J = 0; closed form of L1. Verbatim (both versions): *"we conjecture that there
exist parameter values for which the unique positive equilibrium of the differential equation (16) is asymptotically
stable, and is surrounded by four limit cycles Γ0, Γ1, Γ2, Γ3, which are unstable, stable, unstable, stable,
respectively. […] This is why we leave the existence of four limit cycles a conjecture."* The paper does not define
the normalisation of L_k; it is fixed implicitly by the printed formulas/values (see V2 below).

## 2. What the verifier recomputed independently (scripts in `verification/`, logs `*.log`)

All runs under WSL Fedora, Python 3.12, python-flint 0.9.0, mpmath 1.3.0, sympy 1.14, msolve 0.10.1.

**V1 — exact focal values by a different algorithm (`v1_focal_real_exact.py`, 3 s).** Own derivation of the normal
form: with c=1+b, s=w²=(d−1)/(1+bd), the identity a·u̇ + c·v̇ = a s[(c−1)I + E − c] (checked with sympy) and the
change (u,v) = D(cξ, a(wη−ξ)), D = acw, τ = ωt make the field polynomial in (a,c,w) with linear part (−η, ξ) (checked).
Focal values computed with the **real** Lyapunov-function method (F = ξ²+η²+ΣF_n, homogeneous linear algebra in the
monomial basis ξ^iη^{n−i}, exact sympy inverses, circle-average-zero normalisation of F_{2k+2}, residual checked at
every degree) — a different algorithm and different code from the solver's complex recursion. Result: **V̂_k = w·G_k(a,c,w²)
exactly (polynomial identity in Q[a,c,w]) for k = 1, 2, 3, 4**, G_k read from `certificates/G_k.json`. So the G_k are
exactly the Lyapunov quantities used in the sign lemma (same normalisation; the zero set and the sign of L4 at a
triple zero are normalisation-independent anyway). Also: G_1 gives exactly the paper's closed form of L1 (sympy).

**V2 — normalisation against the paper's printed values (`v2_focal_numeric.py`, generic real-method code, 60 digits,
not a proof).** With coordinates ξ = x−1, η = −(J11u+J12v)/ω, time ωt, L_k = πV_k reproduces: quadrangle §3.2 L1
closed form (two K values, rel. error < 1e−40) and L2 ≈ 0.01293 at K0 (0.012933567); chain §4.1 L1 and L2 closed forms
(rel. < 1e−40) and L3 = −625π/110592·√(7/2) at q=1/4 (all digits); BH checkpoints (a,d)=(1,165/49): L1=L2=0, L3<0, and
((1+√3961)/60, 3): L1=L2=0, L3>0. At the certified p*: |L1|,|L2|,|L3| ≈ 1e−39, L4 = −2.58012649949442448635653476266e−6,
agreeing with the solver's enclosure to 29 digits. Hence the solver's L_k are the paper's L_k.

**V3 — independent validated integrator for Theorem C (`v3_polar_validated.py`).** Different formulation: polar
coordinates in near-normal linear coordinates ξ = u, η = κ(au+cv) (κ a binary number ≈ 1/ω), with the **angle φ as
independent variable**, dR/dφ = R(cosφ ξ̇ + sinφ η̇)/(cosφ η̇ − sinφ ξ̇). Positivity of the denominator is certified
on every a-priori enclosure, so φ is monotone along the orbit and the first return to Σ = {(1+ρ,1)} (the ray φ = φ0) is
exactly R(φ0+2π)cosφ0 — no crossing detection; scalar ODE. Each step: a-priori enclosure R_k + [0,h]F(B,Φ_k) ⊂ int B,
Taylor coefficients by Picard iteration in Arb power-series arithmetic (`arb_series`), Lagrange remainder from the
N-th coefficient computed on (B, Φ_k); M = 1024 steps, N = 32, 320 bits. (A first attempt with polar coordinates in the
(u,v) plane was correct but too ill-conditioned near the major axis; replaced.) Results at the exact parameters of
`certificates/explicit_params.json`:

| ρ | verifier's δ(ρ) (Arb) | solver's δ(ρ) |
|---|---|---|
| 0.002 | [−4.6043873099906e−25 ± 1.4e−39] | −4.604387309990613716e−25 |
| 0.004 | [+1.2735661359982e−24 ± 9.1e−39] | +1.273566135998190960e−24 |
| 0.006 | [−3.5491458181055e−24 ± 2.5e−39] | −3.549145818105502444e−24 |
| 0.012 | [+1.4112630991155e−22 ± 3.2e−36] | +1.411263099115531169e−22 |
| 0.024 | [−3.0291128765253e−20 ± 4.3e−34] | −3.029112876525257070e−20 |

(solver's column: its own enclosures re-evaluated at 25 digits, radii ≤ 5e−45; the strings stored in
`certificates/validated_cycles_output.json` are 12-digit displays). The two independent integrators agree to all
13–14 printed digits (relative differences ≤ 2.2e−14, i.e. the display rounding of the verifier's log), the
per-step remainders are ≤ 2.3e−57, and the signs are −, +, −, +, −. This reproduces Theorem C with an independent
method (angular monotonicity in the normal-coordinate angle was certified along all five orbits as well).

T2 (evidence for Theorem A, not using any focal-value code): at a rational point exactly on tr J = 0 within 1e−40 of p*,
the validated displacements are δ(0.0005) = −1.604495114e−34, δ(0.001) = −8.184962034e−32, δ(0.002) = −4.160156683e−29,
δ(0.004) = −2.099188082e−26 (radii ≤ 1e−36, all negative); δ/ρ⁹ = −8.2150e−5, −8.1850e−5, −8.1253e−5, −8.0078e−5,
whose Richardson extrapolation −8.24507e−5 matches the prediction L4(1+1/s)⁴ = −8.24522e−5 from the certified L4
(relative difference 2e−5, the expected O(ρ²) residue). So the displacement really behaves like a negative multiple
of ρ⁹ near p* — L1 = L2 = L3 = 0, L4 < 0 and the normalisation are confirmed by direct integration.

**V4 — Krawczyk re-certification with another interval library (`v4_krawczyk_mpmath_iv.py`).** mpmath.iv (outward
rounding), own Newton centre, own box radius 1e−30, nested (Horner-type) evaluation, own derivative code:
K(X) ⊂ int X (max |K_i − m_i| ≤ 3.6e−53); domain conditions on X; G4 < 0 on X with
L4 ∈ [−2.58012649949442448644739e−6, −2.58012649949442448626568e−6]; det ∂(G1,G2,G3)/∂(a,b,d) ∈ −0.0347036714263651341…
(excludes 0); the verifier's zero enclosure lies inside the solver's Krawczyk box (same zero).

**V5 — Theorem B with msolve (`v5_msolve_uniqueness.py`).** {g1=g2=g3=0} has a positive-dimensional complex component
(inside {s = −1}); after saturation by a·c·s·ρ·(1+s) (all factors > 0 on the domain) the system is zero-dimensional and
msolve isolates **5 real solutions, exactly one of which lies in the domain** (a,c,s) ≈ (1.0128219987, 1.6546256211,
0.7258957623) = μ*. The others have s<0 or a<0 (they match the solver's excluded candidates t ≈ 0.6196, 0.7405, 2.7323
and s ≈ 0.11924). (msolve's parametrisation is computed multi-modularly; I treat this as an independent cross-check of
the solver's resultant + Arb argument, which I also reviewed and found correct.)

**V6 — factual statements (`v6_misc_checks.py`).** g1…g4 integral, primitive, irreducible over Q, degrees
(2k,2k,3k−1), 9/48/136/296 terms; the appendix formulas for g1, g2 equal the certificate files; the Corollary's
mass-action form (κ1=d−1, κ2=1, κ3=K^−a, equilibrium (K,1), x=KX, y=Y, t=Kτ) maps to (BH); the explicit parameters
satisfy tr J = −1.90e−22 < 0 < det J, tr² < 4det, K − a/(1+bd) = −6.03e−23.

**Solver's pipeline re-run.** `scripts/run_all.sh`: 35 PASS, 0 FAIL, plus the 4 non-rigorous polar checks
(`verification/run_all_rerun.log`), identical certificate outputs.

## 3. What the verifier checked by reading (proof review)

* **Normal form (§3).** Re-derived by hand: u̇ = f, v̇ = as(I−1) − a(1+s)f/c on T0, P = N_f/(acw), Q = (cN_I−N_f)/(ac),
  D(P+iQ) = N_f + iw(cN_I − N_f); the parametrisation Ψ:(a,c,s) ↦ (aρ/c, a, c−1, (1+s)/ρ) is a bijection onto T0
  (d>1 ⇔ cs>0, 1+bd = c/ρ), ∂s/∂d = (1+b)/(1+bd)² > 0. Correct.
* **Krawczyk certificate (`certify.py`).** Exact rational centre and radius; Arb balls contain the exact box (checked
  in the code); G(m) evaluated on exact rationals; the interval Jacobian is evaluated on balls containing the box;
  the inclusion K(X) ⊂ int X is tested strictly against the exact box endpoints; ‖I − CJ‖∞ < 1 gives nonsingularity;
  C5 via the chain rule (det of ∂(a,c,s)/∂(a,b,d) = ∂s/∂d > 0) and positive row scaling at a common zero. Correct.
* **Sign lemma, return-map lemma, trace lemma, successive perturbations (§4.2, §6.1–6.3).** Correct. The sign order is
  right: L4<0 at p*; perturb G3>0 (δ>0 near 0), then G2<0, G1>0 (all on T0 via Φ = Ψ∘G⁻¹, inverse function theorem
  from C2–C3), then decrease K (tr J = −η(1+bd) < 0, Lemma 6.2) — signs of δ from outside in −,+,−,+,− at
  ρ1>…>ρ5, hence sign-changing zeros −→+ (repelling), +→− (attracting), −→+, +→−: equilibrium stable, Γ0 unstable,
  Γ1 stable, Γ2 unstable, Γ3 stable — exactly the conjectured configuration. Lemma 4.2 is only applied on T0 and
  with exact vanishing of the lower G_j; Lemma 6.1 provides continuity in p on a 4-dim neighbourhood, needed for the
  last (off-T0) step. The Bautin-type finite-cyclicity statement is only used in a remark (exactly four, hyperbolic),
  correctly labelled "not used, not re-proved".
* **Validated integrator of Theorem C (`validated_cycles.py`).** Read line by line: Taylor recurrences of the
  polynomial lift (x,y,E,I,J) are correct; the lift is consistent (I·y and J·x are invariant = 1, then E ≡ x^a y^b —
  a one-line justification was added to the paper); the a-priori test Z + [0,h]F(B) ⊂ int B is the standard
  Picard–Lindelöf argument (F(B) is a box, hence convex); the remainder z_N(B)τ^N is the Lagrange remainder of an
  autonomous system; the crossing bookkeeping is sound (monotone y on every step whose enclosure meets y = 1,
  left-run opened above and closed below y = 1, exactly one left run, return step starts below and ends above,
  ambiguous end points re-stepped, bisection keeps y(lo)<1<y(hi)); angular monotonicity uv̇ − vů > 0 is certified on
  every enclosure. The Jordan-curve/positive-invariance argument (polar graph, segment entered with increasing θ,
  Poincaré–Bendixson with (1,1) the only equilibrium and a sink, transversality of Σ) is correct.
* **Theorem B (`uniqueness.py`).** Parametrisation of {G1=0} by (t,s) (G1=0 ⇔ ρt(2t(1+s)+1)=1, c>0 ⇔ 2t(1+s)>1),
  factor bookkeeping, resultant logic (a common zero forces R_t(t0)=R_s(s0)=0), Arb root enclosures of each
  irreducible factor, candidate filtering and final inclusion in the Krawczyk box: correct (and confirmed by V5).
* **Novelty.** Semantic Scholar (13 citing records) and OpenAlex (12) citations of arXiv:2103.00972 through 2026,
  a web search, and the `mmislan` (Palomar) repositories: no proof of the four-cycle conjecture found.

## 4. Problems found and corrections applied by the verifier

1. **Abstract of `paper.tex`: wrong interval.** It stated L4 ∈ [−2.58012649949·10⁻⁶ ± 10⁻³⁴], i.e. a 12-digit
   midpoint with a 10⁻³⁴ radius — this interval does *not* contain the true value (off by ≈ 4·10⁻¹⁸). Replaced by
   "L4 = −2.580126499494…·10⁻⁶ < 0 (rigorous enclosure of radius 10⁻³⁴)". (The theorem statement and table were
   already correct.)
2. **Remark (a) in §6** justified the need to vary the exponents by citing arXiv:2406.13451 (read only at abstract
   level, and about *quadratic* networks) as if it covered this real-exponent network. Replaced by the direct
   argument from the paper's own scaling (for a fixed network the rate constants enter only through K), with the
   citation kept as a comparison.
3. **§6.4:** added the one-line proof that the polynomial lift projects onto (BH) (invariants Iy = Jx = 1).
4. **Stale references** to a non-existent `notes.md` in `references.md`, `scripts/focal_exact.py` and
   `scripts/validated_cycles.py` replaced by the corresponding `paper.tex` sections; wrong runtimes (“1–3 min”) and
   output file name (`Vtilde_k.txt` → `.json`) in the same docstrings corrected. (Docstrings only; the pipeline was
   re-run afterwards: unchanged, all PASS.)
5. **Status and credits** updated in `paper.tex` (`\PaperCredits`, subtitle, trust-base section, results table) and
   `README.md` (attribution, status "solved", file list); `paper.pdf` recompiled (12 pages).
   *Not done (outside this folder):* the results table of the root `README.md` should be updated by the orchestrator.

## 5. Remaining reservations

* The proofs are computer-assisted; they rely on Arb/FLINT (python-flint 0.9.0). The verifier's own checks used a
  second interval library (mpmath.iv) for the Krawczyk step and msolve for Theorem B, but the verifier's validated
  integrator again uses Arb (independent code, same library), so an Arb-level bug would be a common mode for the two
  integrators (not for Theorem A's certificate, which was redone in mpmath.iv).
* Exactly four small cycles and their hyperbolicity are not proved (only cited from Bautin theory, as the paper says);
  the conjecture only asks for (at least) four with the stated stabilities, which is proved. Nothing is claimed about
  large-amplitude cycles.
* Exponents are real (generalised mass action), as in the source's setting (a>0, b>−1, d>0 real); the explicit
  parameters of Theorem C are long decimal fractions, and κ3 = K^−a is irrational — legitimate "parameter values".
* T2 of V3 (the limit δ/ρ⁹ → L4(1+1/s)⁴) is numerical evidence built on validated δ values; it is not part of any proof.
* V5 relies on msolve's multi-modular computation (probabilistic in principle); Theorem B is in any case a bonus and
  its original proof (exact resultants + Arb) was reviewed and is correct.
