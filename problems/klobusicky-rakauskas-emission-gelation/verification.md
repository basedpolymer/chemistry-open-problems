# Verification

**Verdict: CONFIRMED WITH CORRECTIONS.**

Verifier: **GPT-6-Astra**, separate agent with a fresh context, 28 October 2026. The verifier received the original source and statement, the claimed result and deliverables, but not the solver's conversation or derivation history. No sub-agent was used by the verifier. This is an agent audit, not human peer review or a proof-assistant certificate.

The theorem and every displayed formula are confirmed for all integers `k > ell >= 1`. No mathematical formula or proof argument required correction. The corrections concern the README's unqualified uniqueness/positivity wording and the protocol's requirement to quote a complete source sentence. They were applied by this verifier, together with attribution and status updates. Details appear below.

The earlier solver-authored internal audit is preserved verbatim in `verification/internal-audit-original.md`. Its pending-verification statements are historical and are superseded by this report. That audit is not being relabelled as a verification.

## Source and target checked

I first read `PROTOCOL.md`, then the original source text, especially Definition 1, the mass balance, equation (37), the monodisperse moment formula (96), and Sections 7--8. I checked the [arXiv record](https://arxiv.org/abs/2604.14076), which displayed only v1, submitted 15 April 2026. The source was read from a local copy of the arXiv v1 PDF, which is not part of the repository.

The original deterministic equation is

\[
\dot u_n=m_1^{-2}\sum_{i=\ell+1}^{n-1}i(n+\ell-i)u_i u_{n+\ell-i}
             -2nu_n/m_1,\qquad m_1=\sum_{n>\ell}nu_n.
\]

The ordered gain sum, factor 2, reaction-count time scale and actual mass denominator all agree with the manuscript. The source proposes the threshold `k/(2k-ell)` and establishes existence only locally. Its pre-gelation class requires uniform convergence of the second-moment series on compact intervals. This is the class in which the audited theorem proves uniqueness.

## Separate derivation

Write `d=k-ell`, `n_j=ell+d*j`, `p=(2k-ell)/ell` and `a=p-1=2d/ell>0`. Since shifted sizes add, only the sizes `n_j` can be created. In the full triangular equations, induction on the integer size proves that zero initial concentrations off this support remain zero; this conclusion does not rest only on a physical reachability argument.

With the *candidate* mass `k-ell*t`, put `M=1-ell*t/k` and

\[
u_{n_j}=M^{2n_j/\ell}A_j(t/k^2)^{j-1}.
\]

Direct substitution cancels the loss term and gives

\[
A_1=1,\qquad (j-1)A_j=\sum_{i=1}^{j-1}n_i n_{j-i}A_iA_{j-i}.
\]

I rederived the solution of this recursion through a one-variable generating function, without using the manuscript's characteristic computation. Let `F(z)=sum n_j A_j z^j` and `D=z d/dz`. The recursion yields

\[
DF-F=(\ell+dD)F^2,
\quad zF'(1-2dF)=F(1+\ell F).
\]

Integration as a formal identity, fixed by `F(z)=kz+O(z^2)`, gives

\[
kz=\frac{F}{(1+\ell F)^p}.
\]

Equivalently, defining `q=ell*k*z` and `Y=1+ell*F`,

\[
Y=1+qY^p,\qquad Y(0)=1.
\]

Lagrange inversion for `Y-1=q(1+(Y-1))^p` gives

\[
[q^j](Y-1)=c_j=\frac1j\binom{pj}{j-1}
             =\frac1{aj+1}\binom{pj}{j}.
\]

The coefficient comparison is precisely

\[
A_j=\frac{\ell^{j-1}k^j}{j n_j}\binom{pj}{j-1}.
\]

This derivation holds for rational noninteger `p` as well as integer `p`. Every factor in the generalized binomial coefficient is positive because `p>1`. The transformation `H=1-1/Y` recovers the manuscript's implicit equation and confirms equivalence of the two derivations.

The real map `q=(Y-1)/Y^p` increases on `1 <= Y <= p/a`, reaching

\[
q_c=a^a/p^p.
\]

Stirling's formula applied to the gamma-function expression for `binom(pj,j)` separately gives

\[
c_j\sim \frac{\sqrt p}{\sqrt{2\pi}\,a^{3/2}}\,j^{-3/2}q_c^{-j}.
\]

Thus the radius and endpoint convergence asserted in the manuscript are correct, including for `a<1`. With `q(t)=(ell*t/k)M^a`, the actual mass generating function is

\[
\sum_{j\ge1}n_j u_{n_j}(t)x^j
=\frac{k^2M^2}{\ell t}\bigl(Y(q(t)x)-1\bigr).
\]

For `0<t<t_g=k/(2k-ell)`, the value `Y(q(t))=1/M` lies on the branch starting at 1. Hence the actual mass is `kM=k-ell*t`. This verifies the nonlinear denominator *after* construction and resolves the possible circularity. The removable limit at `t=0` is `kx`.

Implicit differentiation gives

\[
qY'(q)=\frac{Y(Y-1)}{p-aY}.
\]

Applying `ell+d*x*d/dx` to the mass generating function therefore yields

\[
m_2=(k-\ell t)\left(\ell+\frac{d}{1-(2k-\ell)t/k}\right),
\]

which simplifies to the manuscript's formula. Its pole is exactly `t_g` and has positive divergent sign because `d>0`. Also `t_g<1` and `t_g<k/ell`, so neither the cluster count nor the proposed mass is exhausted there.

At the endpoint, `M=a/p`. Multiplying the coefficient asymptotic by `k^2M^2/(ell*t_g*n_j)` gives

\[
u_{n_j}(t_g)\sim\frac{k}{\sqrt{\pi d(2k-\ell)}}j^{-5/2}.
\]

Both the exponent and the constant are confirmed. This is an asymptotic in the reachable-size index `j`; it is not a statement that all integer sizes have positive concentration. The first moment converges and the second moment diverges at the endpoint.

## Analytic and scope checks

1. **Existence, not just a formal series.** On every `[0,T]` with `T<t_g`, choose `R>1` with `Rq(T)<q_c`. The nonnegative mass series at `x=R` is bounded uniformly, including its removable limit at zero. This controls every polynomial moment tail exponentially and makes the explicit componentwise C1 formulas genuine solutions with their actual mass.
2. **Uniform endpoint mass.** Set `s=ell*t/k`. The mass coefficient is proportional to `c_j q(t)^(j-1)(1-s)^(a+2)`, uniformly bounded by `c_j q_c^(j-1)`. This is summable by the preceding `j^(-3/2)` asymptotic. Endpoint mass and cluster count therefore follow without termwise differentiation of a divergent moment. Component derivatives have finite one-sided limits because each component equation has a finite gain sum and positive mass.
3. **Uniqueness class.** Uniform second-moment tails justify the mass-balance summation for every competing pre-gelation solution. The absolute weighted gain sum is bounded by `2*m1*m2`, and its tails are uniformly controlled by the first- and second-moment tails. On compact intervals of positive mass, the denominator is bounded away from zero. Thus `m1'= -ell`; the resulting triangular equations then fix every component. This does not prove uniqueness among arbitrary componentwise solutions with infinite or nonuniformly summable second moments.
4. **Maximal pre-gelation interval.** Any pre-gelation continuation through `t_g` would coincide with the explicit solution before it, and would have a bounded second moment on a compact interval containing `t_g`. The explicit pole excludes this. There is no earlier second-moment singularity because all moments were already shown uniformly finite on subcritical compact intervals.
5. **Original-model restrictions.** Nonnegativity and `m0=1-t` imply the source's component bounds `u_n<=1`. Initial data have mass `k`, not 1; I checked that no normalization inadvertently replaces this with unit mass. No stochastic convergence, giant-cluster formation, gel fraction or post-gelation dynamics follows from this theorem.

## Exact computation performed

I inspected and reran `scripts/verify_exact.py` in the supplied WSL Fedora scientific Python environment. **85 PASS, 0 FAIL; exit code 0.** Its saved machine-readable output is `certificates/exact-checks.json`. The recurrence tests include noninteger `p`; direct dimer substitution and the rational tail bounds are valid supplementary checks. They do not establish the infinite-parameter theorem on their own.

I wrote and ran the separate script `verification/check_original_series.py`. **40 PASS, 0 FAIL; exit code 0**, approximately 6 seconds. It uses only exact fractions from the Python standard library, for 10 pairs `(ell,k)`, including `(100,101)` where `a=1/50`, through Taylor order seven. It calculates the mass denominator from all concentration coefficients at each order and never inserts `m1=k-ell*t` while solving. Only afterward does it compare the coefficients with the proposed formula and check the mass, count and off-support values. The finite size cutoff is exact at each retained order because forming a cluster of index `j` requires at least `j-1` gains. Output: `verification/original-series-results.json`.

Reproduce with:

```sh
python scripts/verify_exact.py
python verification/check_original_series.py
```

No floating-point exploration is used as a certificate. The all-parameter proof remains the analytic argument above and in `paper.tex`.

## Corrections applied by the verifier

- `README.md`: replaced the unqualified assertion of uniqueness with uniqueness in the precise pre-gelation class; restricted positivity to reachable sizes and positive times. The theorem already had the correct class and formula.
- `README.md` and `paper.tex`: replaced the short dimer quotation fragment with the complete short sentence from p. 21 proposing the monodisperse threshold; paraphrased the dimer context. `references.md` now describes the revised quotation. Mathematical notation was transcribed without changing its content.
- `README.md` and `paper.tex`: added verifier attribution, updated the status, and documented the separate checks. No theorem, coefficient, asymptotic constant or solver certificate script was changed.
- Preserved the old internal report and added this report and the separate exact Taylor script/output. The PDF was intentionally not rebuilt by this verifier; the coordinating agent will rebuild the final manuscript.

## Limits and novelty

The verdict confirms the deterministic monodisperse result and its stated pre-gelation uniqueness. It does not cover general initial distributions, arbitrary weak solution classes, a hydrodynamic limit or any post-gelation model.

The verifier also made a bounded fresh web search on the arXiv identifier, exact title, and author names with gelation/exact solution. The [arXiv record](https://arxiv.org/abs/2604.14076) and [author's publication list](https://joeklobusicky.net/publications/) showed the same submitted paper and no inspected follow-up resolving the claim. This is consistent with the solver's documented bounded novelty audit; it is not proof of priority and does not turn API access failures into successful literature checks. No external publication, message to authors or Git push was performed.

## Note on this file

The verifier's attempt to save this report failed (its patch deleted and re-created the same file in one operation and was rejected), so the solver's internal audit was still in place. The text above was restored verbatim from the verifier's session log on 28 October 2026; only the attribution wording and the description of the local source copy were edited. The PDF was rebuilt after the verifier's corrections.
