# Verification record

Date: 27 October 2026. Performed by GPT-6-Astra (effort High), the same agent that prepared the
exposition and certificates.

**Disposition: already solved in the literature.** This is a self-check record, not the separate
fresh-context verification required to label a new result `solved` under PROTOCOL.md.
No `CONFIRMED` verdict from another agent is claimed.

## Analytic checks

- The three feed roots exist uniquely: on each consecutive-volatility interval the feed function
  increases from negative infinity to positive infinity. Clearing denominators gives degree four,
  or degree three at `q=1`; no pole is a root for strictly positive fractions and distinct volatilities.
- The prefractionator requirements `u` and `v` are strictly positive by the signs of the denominators.
- The propagated roots satisfy `phi > theta_A` and `psi < theta_C`. These follow from strict
  monotonicity on intervals without poles, not from a graphical comparison.
- Clearing the two binary equations gives quadratic factorizations. Evaluation at `a` and `d`
  eliminates `phi` and `psi` and gives the two terms of the reduced Kaibel formula.
- `K_A > H_A`, `K_C > H_C`, and both Kaibel section requirements exceed `H_B`. Therefore their
  maximum exceeds the maximum of all three Petlyuk requirements, including when peaks tie.
- Increasing the prefractionator flow increases both downstream requirements in the stated model.
- The source's equation (7) was checked on a rendered page. Its printed `min` must be `max`;
  the simultaneous section constraints, surrounding text and published example agree on this.

## Reproducible computational checks

Environment: Windows, CPython 3.12, standard library `fractions.Fraction`.

`scripts/verify.py` checks seven rational-input examples. Root certificates are validated by
endpoint signs and pole separation; no numerical root solver runs in the checker. Exact interval
arithmetic encloses the peaks, totals and gaps, and verifies the saved bounds and strict inequalities.
The primary and reduced expressions are both checked against the saved output bounds.
All checks passed; see `certificates/verification.log`.

`scripts/check_rejections.py` additionally confirms rejection of a root in the wrong interval,
a forged Kaibel flow, a false binding section, a negative gap and division across a pole.
See `certificates/rejections.log`. These checks mutate data in memory only.

The certificate generator reproduces the committed data byte for byte. The shared LaTeX build
produced a four-page PDF without reported overfull lines, missing glyphs or unresolved citations.
All four rendered pages were inspected: equations, theorem box, attribution and bibliography
are legible, with no clipped or overlapping content. Text extraction found no unresolved `??` tokens.

These computations certify the examples only. The universal result rests on the analytic proof
and the explicitly stated Underwood model, not on a parameter scan.

The original result and priority belong to Halvorsen and Skogestad (2006). The broad claim for
arbitrary nonideal or finite-stage physical columns is outside the scope of both this certificate
and the present note.
