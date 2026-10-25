# Uniqueness of the Brownian coagulation similarity profile

**Status: open; selected research target, 25 October 2026.** No new uniqueness theorem or numerical certificate is claimed.

**Attribution.** Selection, literature audit and proposed approach: GPT-6-Astra (OpenAI; effort xhigh), without sub-agents.

The target is to prove or disprove uniqueness of the finite-mass self-similar profile for the **full Brownian kernel** in Smoluchowski's coagulation equation. It is a precise problem in aerosol and colloidal aggregation, with a known existence theory and a substantial gap between existing uniqueness results and the physical kernel. The literature audit is in [references.md](references.md).

## Original open statement

B. Niethammer, S. Throm and J. J. L. Velázquez, *A revised proof of uniqueness of self-similar profiles to Smoluchowski's coagulation equation for kernels close to constant*, [arXiv:1510.03361v3](https://arxiv.org/abs/1510.03361v3), §1, printed p. 2:

> For non-solvable kernels such as (1.2) so far only existence results are available [2, 1], while uniqueness of self-similar solutions remained unknown.

Their equation (1.2) is the Brownian kernel below, up to a positive time-scale factor. This quotation is historical evidence of the open question; the separate 2026 audit checks subsequent work. The problem also appears as [MColbrook/AIM entry 163](https://github.com/MColbrook/AIM/blob/8729eab98bafe409b66c5e39253bf56569a62327/problems/163-brownian-coagulation-profile.md), still marked open when consulted.

## Exact target

Fix

\[
K_B(x,y)=2+(x/y)^{1/3}+(y/x)^{1/3},\qquad x,y>0.
\]

Does there exist exactly one nonnegative measurable function \(F\), up to equality almost everywhere, satisfying

\[
\int_0^\infty xF(x)\,dx=1,
\qquad
\int_0^\infty (1+x+x^{-1/3})F(x)\,dx<\infty,
\]

and, for every \(\varphi\in C_c^\infty(0,\infty)\),

\[
\int_0^\infty F(x)[x\varphi'(x)-\varphi(x)]\,dx
=\frac12\int_0^\infty\!\int_0^\infty
K_B(x,y)F(x)F(y)
[\varphi(x+y)-\varphi(x)-\varphi(y)]\,dx\,dy?
\]

This fixes the similarity convention \(f(t,x)=t^{-2}F(x/t)\). The mass normalization removes the dilation freedom; it does not prescribe particle number. The integrated profile formulation used in the primary papers is

\[
x^2F(x)=\int_0^x\int_{x-y}^{\infty}yK_B(y,z)F(y)F(z)\,dz\,dy.
\]

The weak formulation above is the authoritative target. Any proof using the integrated formulation must justify the endpoint limits and its equivalence in the stated class.

## Why this is a substantial target

The collision kernel comes from diffusing spherical particles; its dependence on the ratio of masses couples very small and very large clusters. It is unbounded as that ratio diverges. Uniqueness would settle whether this physical model admits competing finite-mass similarity shapes. Global convergence of arbitrary time-dependent data is a further question and is **not** implied by uniqueness alone.

The strongest relevant perturbative source located is S. Throm, *Communications in Mathematical Physics* 383 (2021), 1361–1407, [DOI 10.1007/s00220-021-03967-6](https://doi.org/10.1007/s00220-021-03967-6). It allows singular perturbations but retains a sufficiently-small-amplitude condition. Writing

\[
K_\varepsilon(x,y)=2+\varepsilon[(x/y)^{1/3}+(y/x)^{1/3}]
\]

places the physical target at \(\varepsilon=1\). The available theorem does not establish that its allowed interval reaches 1. An overall rescaling of time does not remove this relative-amplitude obstruction.

## Proposed approach, not an established proof

1. Work with the integrated profile equation in weighted function spaces. Derive bounds on moments and both endpoint tails for **every** admissible profile, with constants valid at the physical kernel.
2. Investigate continuation from the explicit profile \(F_0(x)=e^{-x}\) at \(\varepsilon=0\), using validated numerical estimates for the linearized operator on the mass-zero subspace where useful.
3. Establish a global exclusion argument: continuing one locally unique branch does not exclude disconnected branches, folds, or profiles escaping the controlled function space. This is the central unresolved step, not a consequence of numerical continuation.

A complete positive solution must compare arbitrary admissible profiles at \(\varepsilon=1\). A complete negative solution requires two rigorously distinct profiles for that same kernel and normalization. A discretized solution, an interval of small amplitudes, or local uniqueness near one computed profile is insufficient.

## Files and reproduction of the audit

- [references.md](references.md): primary sources, exact sections inspected, recent-work comparisons and search limitations.
- [JOURNAL.md](JOURNAL.md): selection decisions, including a discarded already-solved candidate.
- `certificates/sources/`: local-only reference PDFs, the AIM snapshot, a rendered source page and citation-search data; ignored by Git under the repository's source policy.

No solver has been run for this new target. No independent verification is claimed. All three principal preprints are freely available; no article purchase is needed to begin.
