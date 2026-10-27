# Candidates of agent Q1 (persistence of the dual) for agent E2

Agent Q1, Claude Opus 5.5 (effort xhigh), 27 October 2026. Definitions and proofs: `notes.md`. Code: `scripts/`
(`s1lib.classify(p, rx)` classifies a rank-one conservative network; `drainsets.drain_sets(p, rx)` finds drain
sets; `s1lib.dual(n, rx)` builds N*). Notation: N* = dual, κ = rate constants of N*; "ALL / NONE / SOME" =
persistent for all / no / some-but-not-all rate vectors κ.

| id | statement | hypotheses | status | suggested test |
|---|---|---|---|---|
| Q1-1 (S1) | For a conservative rank-one network M: persistent for all κ iff every face is IN; some OUT or ESC face ⇒ non-persistent for every κ; only one-sided DOM faces besides IN ⇒ SOME | rank Γ(M) = 1, positive conservation law | proved | **independent exact decision**: for sampled rational κ and sampled classes (x₀ rational), compute F(ξ) on the segment exactly and the sign of F near both endpoints (lowest-order term); compare with the verdict (ALL ⇒ never attracting; OUT ⇒ attracting in every class whose endpoint has that zero set; ESC ⇒ for the sampled κ find the class along the separating direction) |
| Q1-2 (S1-WR) | No face of N* is DOM; hence N* is ALL or NONE | N of rank 1, every linkage class has its two extreme complexes (along the reaction direction) as sources, e.g. WR | proved | on E's s = 1 records and beyond: classification never returns DOM; with Q1-1's exact decision, no rank-one WR network whose dual is persistent for one sampled κ and not for another |
| Q1-3 (C-2) | n = 2, s = 1: N* is ALL iff every source has strictly more of its consumed species than of its produced species | (H1)–(H3), n = 2 | proved | n = 2 families |
| Q1-4 (C-pure, s = 1) | If no reaction of N has a net product in its source (e.g. N pure), N* is ALL | s = 1, (H2) | proved | rank-one families |
| Q1-5 (D) | Drain set ⇒ for every κ and every class, x_Z → 0 from points with ⟨w, x_Z⟩ small | M conservative | proved | numerical integration (any precision) on a few drain-set duals of each s, random κ, starting at x_Z = η, x_{Z^c} = 1: x_Z should decay; and the exact check that Z contains no conservation-law support (proved consequence) |
| Q1-6 | N weakly reversible and pure (no species on both sides of a reaction) ⇒ N* persistent for all κ | WR, (H1)–(H3), any s | conjectured (true for s = 1 by C-pure; false without WR at s = 2, example Q) | E's WR records: pure ones; check whether N* has no critical siphon (ADLS ⇒ ALL) or a drain set / OUT facet (⇒ refuted) |
| Q1-7 (dichotomy, s = 2) | For N weakly reversible with s = 2, N* is persistent for all κ or for none | WR, (H1)–(H3), s = 2 | conjectured (open) | numerical: for s = 2 duals without drain set and with critical siphons, integrate from near-boundary points for several κ (log-uniform over 10⁻²..10²); report any network with persistence for one κ and extinction for another |
| Q1-8 (Conjecture F) | If a facet face Z of N* (rank Γ(N*)[Z,·] = 1, realizable) is OUT, then N* is non-persistent for every κ | WR N, any s | conjectured | s = 2 duals: compute facet types with s1lib.face_analysis applied to the restriction (g = the spanning vector of Γ(N*)[Z,·]); for OUT facets without a drain set, integrate numerically near the facet |
| Q1-9 | Counts: rank-one WR duals ALL / NONE (and none SOME) as in `notes.md` §2; drain sets vs drain species as in §3 | — | checked | re-run `scripts/enum_s1.py`, `scripts/scan_records.py` independently (or reimplement) |
