# Independent verification

**Verdict: CONFIRMED.**

Date: 18 October 2026. Verifier: separate agent `/root/independent_verifier`.

The theorem in `notes.md` is correct for the deterministic balanced directed diblock
with binary charges and the original pair-interaction Hamiltonian. No mathematical
or numerical correction was required. This verdict confirms the argument and
certificates, **not a claim of historical novelty**: DGH already asserted the same
diblock collapse temperature in 1992. The precise qualification is below.

## Independence and material reviewed

The verifier received the original question/source, the claimed root equation,
the deliverables, and instructions to apply the repository protocol. No solver
reasoning transcript was supplied. I read `PROTOCOL.md`, `README.md`, `notes.md`,
`references.md`, the three original scripts, their parameter/data files, and the
relevant original mathematical sources. I did not read `PROGRESS.md`, the internal
brief, the verification request, or a solver research transcript.

I independently reconstructed the argument using the source Hamiltonian and then
checked the displayed proof line by line. I re-ran every supplied certificate and
added a separately written certificate that evaluates the original pair energy
directly and brackets the root using exact rational Taylor bounds. This is an
independent agent review with fresh task context, not a human peer review or a
proof-assistant verification; the original deliverables were necessarily visible.

Source reading for this review:

- Poisat, [arXiv:2502.18964v2](https://arxiv.org/abs/2502.18964v2): Sections 2--3
  and Section 7 question (3), with direct visual checks of pages 4--6, 10, and 21.
  I checked the official version metadata. I did not independently audit the
  unrelated free-energy variational theorem or appendices.
- Derrida, Griffiths and Higgs, [Europhys. Lett. 18 (1992), 361--366](https://doi.org/10.1209/0295-5075/18/4/014):
  the model and transition discussion, especially pages 364--365 and equation
  (13). These two pages were freshly rendered from the local PDF and visually
  checked. The [author-hosted copy](https://www.phys.ens.psl.eu/~derrida/PAPIERS/1992/heteropolym-92-epl.pdf)
  was also inspected.
- The original solver's novelty audit was read, but its complete web/citation
  search was not independently repeated. A priority conclusion cannot be inferred
  from a zero citation count or from an unsuccessful keyword search.

The local source PDF SHA-256 hashes are:

| File | SHA-256 |
|---|---|
| `poisat-2502.18964.pdf` | `F32596FAC3DE6D847C26B80C8E71046952EAD5781BA9D23B97B347316388F76D` |
| `dgh-1992.pdf` | `D9F421DFDCEAEF33ED6301BF0A2764B79044E55196EC7DB69CAB0A5A0F3EAD1B` |

## Source and normalization checks

Question (3) on page 21 is accurately quoted in `README.md`:

> What is the critical temperature for the di-block polymer from (3.4)?

The source defines an unnormalized sum over directed paths, a first step fixed
to one, and the Hamiltonian with (i<j). For a consecutive occupied-site block
of length (ell) and charge (q), its pair energy is exactly

\[
 \frac{q^2-\ell}{2}.
\]

Consequently the reduced partition function is

\[
 Q_\beta(\omega)=e^{-\beta |\omega|/2}Z^{\beta,\omega}
 =\sum_{\text{compositions of }\omega}\prod_B e^{-\beta q_B^2/2}.
\]

This confirms both the sign and the factor two in the proposed equation.
The source's alternative squared-charge parameter is (t=\beta/2), and its
partition function includes (2^{1-N}). The alternative free-energy formula
in the notes, including the finite-volume convention in the numerical code,
is correct.

The endpoint problem in source (3.4) is real: the displayed ranges leave the
last charge unspecified and do not describe the balanced length-(2n) word
literally. The notes state the intended balanced word explicitly. Their bounded
imbalance argument also works: the completely collapsed reduced weight remains
at least (e^{-\beta d_0^2/2}) when the difference between block lengths is
bounded by (d_0). This is a stated repair of an ill-defined display, not a
silent substitution of a different model.

## Independent reconstruction of the key argument

Let (A_m=Q_\beta(+^m)), (A_0=1), and

\[
 W_\beta(z)=\sum_{\ell\geq1}e^{-\beta\ell^2/2}z^\ell.
\]

The decomposition by the final occupied site gives

\[
 A_m=\sum_{\ell=1}^m e^{-\beta\ell^2/2}A_{m-\ell},
 \qquad \sum_{m\geq0}A_mz^m=(1-W_\beta(z))^{-1}.
\]

The positive radius (r_\beta) solves (W_\beta(r_\beta)=1). For (\beta>0),
the series is entire and has positive coefficients at every positive degree.
There is no competing zero on the circle (|z|=r_\beta). Its simple pole gives

\[
 A_mr_\beta^m\longrightarrow
 \left(\sum_{\ell\geq1}\ell e^{-\beta\ell^2/2}r_\beta^\ell\right)^{-1}>0,
\]

with an exponentially small remainder, as justified in the notes. The induction
bound (A_mr_\beta^m\leq1) and a uniform positive lower bound are correct,
including (m=0). At (\beta=0), the explicit composition count handles the
non-entire limiting case.

For an arbitrary binary word (v), every block charge has magnitude at most
the block length. Comparing each composition therefore gives (Q_\beta(v)\geq
A_{|v|}). Restricting to a cut between two words gives

\[
 Q_\beta(uv)\geq Q_\beta(u)Q_\beta(v).
\]

For a specified cut at (i), its weight is (Q(0,i]Q(i,N]). For a missing cut,
there is a unique crossing block ((k,j]), with (k<i<j), so its weight is
the sum of (Q(0,k]e^{-\beta\Omega_{(k,j]}^2/2}Q(j,N]).
Gluing a homogeneous-dominated segment onto each exterior word gives the
inequalities in the required direction:

\[
 \frac{Q(0,k]}{Q(0,i]}\leq\frac1{A_{i-k}},\qquad
 \frac{Q(j,N]}{Q(i,N]}\leq\frac1{A_{j-i}}.
\]

Since the crossing weight is at most one, the missing-cut/cut odds are at most

\[
 \left(\sum_{a\geq1}A_a^{-1}\right)^2.
\]

This is finite whenever (r_\beta<1). Empty end segments have partition
function one, so there is no missing endpoint multiplier. This independently
confirms the principal matching hot-phase bound for **all binary words**.
The finite truncated-renewal proof given in the notes is also valid and avoids
complex analysis for this part of the result.

For the balanced diblock, at most one occupied-site block contains both signs.
Decomposing by its two exterior lengths yields exactly

\[
 D_n=e^{-\beta n}Z_{2n}^{\beta,\mathrm{db}}
 =A_n^2+\sum_{i,j=0}^{n-1}A_iA_j e^{-\beta(i-j)^2/2}.
\]

The no-mixed-block term is (A_n^2), not an additional copy of a boundary term.
The fully collapsed term has weight one, hence (D_n\geq1). When (r_\beta>1),
choose (1<s<r_\beta). Summability of (A_ms^m) bounds the probability that
at least (M) monomers lie outside the mixed block by

\[
 2(1-W_\beta(s))^{-2}s^{-M}.
\]

This also bounds a cut at (i), because that cut forces at least
(\min(i,2n-i)) monomers outside the mixed block. The treatment of a cut
at the sign interface and of (M>2n) is correct.

The equation (W_\beta(1)=1) has a unique positive solution by strict
monotonicity and its two endpoint limits. Below it, a positive density of
expected cuts rules out convergence (L_n/(2n)\to1): the number of cuts is
bounded above by (2n-L_n). Above it, the exponential bound proves that
convergence. These statements establish both thresholds in equation (6) of
the notes without a monotonicity assumption on individual cut probabilities.
The definition is appropriate to the deterministic question; it does not
assert the random-charge power-law conjecture from the surrounding source.

Finally, the exact formula for (D_n) gives exponential growth (r_\beta^{-2n})
for (r_\beta<1), a bounded partition function for (r_\beta>1), and order
(n) at (r_\beta=1). Thus the stated diblock free energy follows. At the
threshold, if (A_i\to a>0), summing along diagonals and using dominated
convergence gives

\[
 \frac{D_n}{n}\longrightarrow
 a^2\sum_{d\in\mathbb Z}e^{-\beta_*d^2/2}=3a^2.
\]

Together with (p_{n,2n}=A_n^2/D_n), this proves (p_{n,2n}\sim1/(3n)).
This critical asymptotic excludes a positive uniform cut bound at the
threshold. No assertion about complete condensation exactly at criticality is
needed for either infimum/supremum in (6).

## Numerical and code audit

All original checks passed before any file was changed: integer-polynomial
identities through (N=16), all binary words through length eight, the root
enclosure, all forty finite-volume rows through (N=100000), and the critical
finite-size check. The original run completed in 13.35 seconds.

I checked the bounds used by the interval code rather than treating a PASS
message as proof. The Gaussian tail uses the first omitted term divided by
one minus an upper bound on all subsequent ratios. The scaled renewal kernels
have mass below one, so every full and truncated coefficient is at most one.
The resolvent convolution then bounds their difference by ((m+1)T).
The omitted diagonal terms are bounded by their maximal kernel weight times
the square of the coefficient sum. Symmetric Arb error balls enclose the
one-sided errors. Uncertainty in the critical parameter and in the scaling
radius is retained. Mpmath only proposes endpoints; Arb checks both inequalities
before they are used. Decimal certificate strings are outward enclosures.

The certified root lies strictly between

\[
\begin{split}
&0.69813170079919943909134137018011279700848771234773,\\
&0.69813170079919943909134137018011279700848771234774.
\end{split}
\]

The strict inequality (\beta_*>2\pi/9) is certified, not inferred from rounded
decimal displays.

I added `scripts/verify_independent.py`. It has no imports from the solver's
scripts and no floating-point acceptance decisions. At the exact rational
(q=e^{-\beta/2}=3/4), it enumerates positions and computes the source's
(i<j) pair Hamiltonian directly for every binary word through length seven.
It checks normalization, homogeneous domination, cut factorization and the
finite cut-odds bound. Separately, even/odd Taylor polynomials for (e^{-x}),
evaluated with `Fraction`, and a geometric bound on the omitted Gaussian tail
prove the coarser independent root bracket

\[
 0.69813170079919<\beta_*<0.69813170079921.
\]

For these Taylor bounds the integral remainder fixes the sign at every
(x\geq0); monotonicity of all Taylor terms is not assumed. The new script's
root check does not depend on Arb, Mpmath or the supplied fifty-decimal result.
Its final containment check only compares the two resulting intervals.

After registering the new script in `run_all.py`, I ran the complete suite
again. Every check passed, with exit code zero. `certificates/independent-verification.log`
records this run and `/usr/bin/time -v`: elapsed time 17.76 seconds, maximum
resident set size 52,364 KiB, CPU usage 97%, no swaps. Numerical workers were
limited to one thread (`ctx.threads=1` and BLAS/OpenMP environment settings).
The resource use is well below 2 GB. The runner's generic log header still says
that a certificate run alone is not independent verification; this report
records the separate review that supplies that protocol stage.

## Historical attribution and limits of the conclusion

DGH's equation (13) is

\[
 1=\sum_{\ell\geq1}\lambda^{-\ell}
       \exp\!\left(\frac{\ell-\ell^2}{2T}\right).
\]

Their diblock paragraph states a collapse temperature about (1.43), determined
by (\lambda=\exp(1/(2T))). Substitution gives exactly
(\sum_{\ell\geq1}\exp(-\ell^2/(2T))=1), hence the present equation with
(\beta=1/T). This is an explicit previous assertion of the temperature, not
merely an unexplained numerical guess. Their discussion compares the collapsed
state with the homogeneous-chain free energy. It does not display the present
uniform hot-phase proof. Poisat already makes the cold-phase statement rigorous
above this same threshold, despite the coarser upper bound quoted in question (3).

It is therefore justified to describe this dossier as an independently verified
complete proof of the **previously asserted diblock threshold**, with a matching
uniform cut bound. It is not justified to claim discovery of an unknown number,
priority for the complete proof, or resolution of the random-charge freezing
conjecture. The tension between the 1992 assertion and the 2026 open-question
list is a historical qualification, not a mathematical contradiction in this
proof. This review does not prove that no earlier complete argument exists.

## Files changed by the verifier

- Added `verification.md`.
- Added `scripts/verify_independent.py` and registered it in `scripts/run_all.py`.
- Added `certificates/independent-verification.log`; refreshed the ordinary
  `certificates/verification.log` by running the suite.
- Added fresh source-page renders `certificates/sources/verifier-dgh-4.png`,
  `verifier-dgh-5.png`, `verifier-poisat-04.png`, and `verifier-poisat-05.png`.

No changes to `notes.md`, the numerical data, or the original two certificate
scripts were necessary. `README.md` status, `PROGRESS.md`, `FINAL-REPORT.md`,
and `verification-request.md` were not edited. No commits, pushes, publications,
messages to authors, or other external writes were made. Read-only access to
the public source pages was used.
