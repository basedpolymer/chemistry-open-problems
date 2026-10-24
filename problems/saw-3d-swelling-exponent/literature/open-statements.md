# Verbatim statements that make L3, U2 and (T) open

Literature agent, 2026-10-24. Quotes are copied from the text I read (pdftotext of the PDF named; mathematical symbols
re-typed in plain text where pdftotext garbled them, e.g. "E_n^{(1)}|omega(n)|^2"). Page numbers are those of the
version read. Ordered by date.

## 1. Madras, Canad. Math. Bull. 57 (2014) 113-118 (publisher PDF, read in full), Section 1, p. 114

> "While equation (1.1) is undoubtedly true, it has not yet been proven in every dimension. [...] However, a proof of
> equation (1.1) remains frustratingly elusive in low dimensions, namely d = 2, 3, 4. Even worse, so far there are no
> non-trivial rigorous bounds in low dimensions, i.e., known positive numbers D_L, nu_L, D_U, nu_U with nu_U < 1 such
> that D_L N^{2 nu_L} <= E_N(||omega(N)||_2^2) <= D_U N^{2 nu_U} (1.2) for all sufficiently large N."

> "The value nu_L = 1/2 would correspond to a "mean-field" bound, based on the intuition that a self-avoiding walk should
> spread out more quickly than an ordinary classical random walk. However, such a bound has yet to be proved for general
> dimensions d. Even more simply, one could reasonably expect a lower bound with nu_L = 1/d, since the average distance of
> all points of an N-step self-avoiding walk from omega(0) must be at least of order N^{1/d}. Alas, no rigorous proof is
> known even for this bound. The purpose of this note is to prove a lower bound with nu_L = 2/(3d)."

> (p. 114, on DCH) "Informally this says that nu < 1, but it does not give an upper bound nu_U < 1 for equation (1.2)."

> (p. 113) "Our result, Proposition 1.1, is the only nontrivial lower bound that has been proven in 2 or 3 dimensions".

(Equation (1.1) is E_N(||omega(N)||_2^2) ~ D N^{2 nu}.) Covers L1 (nu_L = 1/d open), L3 (nu_L = 1/2 open), U2, E1/(T).

## 2. Bauerschmidt, Duminil-Copin, Goodman, Slade, Lectures on self-avoiding walks, arXiv:1206.2092v1 (Clay Math. Proc. 15, 2012)

Section 1.5.2, p. 6:

> "Almost nothing is known rigorously about nu in dimensions 2, 3, 4. It is an open problem to show that the
> mean-square displacement grows at least as rapidly as simple random walk, and grows more slowly than ballistically,
> i.e., it has not been proved that (1.30) cn <= E_n^{(1)}|omega(n)|^2 <= C n^{2-eps}, or even that the endpoint is
> typically as far away as the surface of a ball of volume n, i.e., c n^{2/d} <= E_n^{(1)}|omega(n)|^2. Madras
> (unpublished) has shown E_n^{(1)}|omega(n)|^2 >= c n^{4/3d}."

(Written in June 2012, before DCH and Madras were published; the upper half of (1.30) was partly answered by DCH
(o(n^2), no rate), and Madras' bound appeared in 2014. The n^{2-eps} and cn statements remain open: see items 5, 6, 8.)

Section 1.6.3, p. 9: "For d = 3, there are no rigorous results for critical exponents, and no mathematically
well-defined candidate has been proposed for the scaling limit."

Section 1.3, p. 4 (related, c_n): "Returning to Z^d, in 1963, Kesten [50] proved that lim c_{n+2}/c_n = mu^2 (1.17),
but it remains an open problem (for d = 2, 3, 4) to prove that lim c_{n+1}/c_n = mu (1.18)."

After Theorem 2.6, p. 13: "It is an open problem to improve Theorem 2.6 in d = 2, 3, 4 beyond the result of Kesten [50]
shown in (1.26)."

## 3. Duminil-Copin, Hammond, Self-avoiding walk is sub-ballistic, CMP 324 (2013), arXiv:1205.0401v1, Section 1.5 "Open problems", p. 4

> "Theorem 1.1 is a first step towards proving the conjecture on mean-squared displacement. An interesting improvement
> would be a quantitative bound:
> Question 1. Show that for some eps > 0, lim_{n->inf} n^{-(2-eps)} <||gamma_n||^2> = 0."

> "The theorem rules out the extreme of fast movement by the walk. What about the other extreme? Does self-avoiding walk
> move further than simple random walk?
> Question 3. Show that lim inf_{n->inf} n^{-1} <||gamma_n||^2> > 0.
> In fact, far slower motion by self-avoiding walk has yet to be ruled out. In principle, the walk may be space-filling:
> Question 4. Show that lim_{n->inf} n^{-2/d} <||gamma_n||^2> = inf."

Section 1.2, p. 3: "There are however very few rigorous unconditional results regarding dimensions two and three. The
question of a non-trivial upper bound on mean-squared endpoint displacement is raised in the introduction of [19]."
([19] = Madras-Slade.) Question 1 = U2; Question 3 = L3; Question 4 (with lim = inf) is between L1 and L2 in d = 3.

## 4. Duminil-Copin, Glazman, Hammond, Manolescu, AoP 44 (2016), arXiv:1305.1257v2, Section 1.2-1.3, pp. 2-3

> "Unconditional rigorous statements concerning the global geometry of the model are almost absent in the low
> dimensional cases at present."

> (on extending Thm 1.1 to fixed x) "For Z^2 this may potentially be shown by a tedious ad-hoc construction. For higher
> dimension, significant difficulties occur and it is unclear how to prove this result."

## 5. Slade, Self-avoiding walk, spin systems and renormalization, Proc. R. Soc. A 475 (2019) 20180549, arXiv:1808.04476v2, Section 2.1.3

p. 4: "For d = 2, 3, 4, none of the entries in Table 1 have been proved." (Table 1 lists gamma and nu; d = 3:
nu = 0.58759700(40).)

p. 3: "A central problem in the subject is to prove the existence of the critical exponents and to show that they have
the values listed in Table 1."

p. 5, after Theorem 2.3 ((1/6) n^{2/(3d)} <= R_n <= o(n), R_n^2 = mean of ||omega(n)||_2^2):

> "The lower bound fails to prove that on average the endpoint of a SAW is at least as far away as it is for simple
> random walk, namely n^{1/2}, even though it appears obvious that the self-avoidance constraint must push the SAW
> farther than a walk without the constraint. The upper bound states that R_n/n -> 0 but there is no bound on the rate.
> In particular, it is not proved that there is a constant C such that R_n <= C n^{0.99999}. The large gap for
> d = 2, 3, 4 between the predicted results in Table 1 and those proven in Theorems 2.2-2.3 is an invitation to look
> for more tractable models that ought to be in the same universality class as SAW."

This is the most recent explicit statement by a leading author that L3 and U2 (even U1-type with exponent 0.99999) are
open; it post-dates DCH 2013 and Madras 2014.

## 6. Guttmann, Jensen, On the existence of critical exponents for self-avoiding walks, J. Phys. A 55 (2022) 415206, arXiv:2207.08433v2, Section 1, pp. 5-6

> "However, for d = 2, 3, 4 almost nothing else has been proved. It is not even known that the mean-square displacement
> grows at least as rapidly as simple random walks and slower than ballistically. That is to say, <R^2>_n >= const.n
> has not been proved, nor that <R^2>_n < const.n^{2-eps}. [...] Informally this says that nu < 1, but proving that the
> exponent is 2-eps remains an open problem."

> "we lamented that the "obvious" existence of a critical exponent for two- and three-dimensional SAWs had still not
> been proved. While we are closer, that still remains true today, some 40 years later."

## 7. Hammond, EJP 24 (2019) 1-43, arXiv:1808.10500v3, Section 1.4, p. 3 (c_n side)

> "Hara and Slade [10, Theorem 1.1] used the lace expansion to prove that c_n ~ C mu^n when d >= 5, but the problem of
> showing that c_n mu^{-n} tends to infinity as n -> inf appears to be open when d equals two or three (in which cases,
> it is anticipated to be true)."

> "It is interesting to note that, while the problem of showing that lim sup theta_n is finite is open, it is known from
> classical results that, for any dimension d >= 2, lim inf theta_n < inf."

## 8. Hutchcroft, Dimension dependence of critical phenomena in long-range percolation, ICM 2026 Proceedings vol. 5, 462-482, doi:10.1137/25M180679X (arXiv:2510.03951v1), Section 1, p. 4 of the arXiv PDF

> "Below the upper critical dimension the self-interaction is predicted to be sufficiently strong that it completely
> changes the large-scale behaviour of the model, leading to different critical exponents and scaling limits. Very
> little is known about this rigorously in the context of SAW (although there are very precise conjectures in two
> dimensions [72, 82]), and we will defer our detailed discussion of low-dimensional critical phenomena until our
> discussion of percolation in the next section."

Section 2 (about percolation, but stated generally), p. 6: "What can be said about critical exponents in dimensions
d = 3, 4, 5? Unlike in two dimensions, there are no conjectures regarding exact values of critical exponents nor any
apparent reason that the exponents should be rational."

## 9. Beffara, La marche auto-evitante, Journees mathematiques X-UPS (2016) 103-130, doi:10.5802/xups.2016-03 (CC-BY), Section after Theorem 2 (DCH), p. 113

> "Encore une fois on sait ce qui se passe en grande dimension, et pour d >= 5 la marche auto-evitante se comporte comme
> la marche simple. En dimension 2, on n'a que des conjectures, cf. 6.2 ; en dimension 3, on ne sait absolument rien de
> plus."

## 10. Other "open" sentences recorded for related quantities

- Madras, End patterns, J. Stat. Phys. 53 (1988), abstract: "It is a difficult open problem to show that as N -> inf,
  the fraction of N-step self-avoiding walks beginning with this pattern converges."
- Kesten 1963, Section 1, p. 960: "(1.4) [lim c_{n+1}/c_n = beta] still seems to be open."
- Kesten 1964, Section 1, p. 1129: "Theorem 2 is, however, still far from the expected formulas chi_n ~ A_1 n^{alpha_2}
  beta^n, gamma_{2n-1} ~ A_2 n^{alpha_3} beta^{2n}."
- Hutchcroft, ECP 23 (2018), Section 1: "it is a major open problem to determine the true asymptotics of c_n when
  d = 2, 3, 4"; "For d = 2, 3 the conjecture is wide open."
- Krachun-Panagiotis (arXiv:2310.17299v1), Section 1, p. 2: "A major obstacle in that direction seems to be the lack of
  a quantitative bound on the partition function of bridges of a given height at the critical point x_c = 1/mu that
  goes to 0 as the height goes to infinity."

## 11. Sources requested but not read

- Madras-Slade, The Self-Avoiding Walk (Birkhauser 1993; Modern Birkhauser Classics 2013), doi:10.1007/978-1-4614-6025-1:
  its open-problem statements (introduction, Section 1.4 and the notes) are NOT quoted here because I could not read
  the book (paywalled; listed in NEEDS.md). DCH (item 3) says the question of a non-trivial upper bound "is raised in the
  introduction of [19]" (= Madras-Slade).
- Slade, The self-avoiding walk: a brief survey (Surveys in Stochastic Processes, EMS 2011): not found on arXiv, not
  read.

## 12. The two doubtful citations in the owner's report

1. "discussion des basses dimensions a l'ICM 2026", doi:10.1137/25M180679X. Crossref: T. Hutchcroft, "Dimension
   Dependence of Critical Phenomena in Long-Range Percolation", Proceedings of the ICM 2026, Volume 5: Invited Lectures
   (Sections 9-11), pp. 462-482 (July 2026). arXiv:2510.03951. It is a survey of **long-range percolation**; SAW appears
   only in the introductory Section 1 (one paragraph on SAW and d_c = 4, quoted in item 8, plus a figure of a 3D SAW).
   It contains no result and no specific open problem about nu for SAW in d = 3. The citation is real but describing it
   as a discussion of low-dimensional SAW overstates it; it supports only the general sentence "very little is known
   rigorously" for SAW below d = 4.
2. Guttmann-Jensen 2022, doi:10.1088/1751-8121/ac943a, arXiv:2207.08433. Claims: they describe three ideas of
   Hammersley (super-multiplicativity of "worms" = SAWs whose endpoints share an x-coordinate, which would give
   c_n <= C mu^n n^2; a two-layer lattice free-energy expansion kappa(p) = kappa + O(-p log p); the ratio criterion
   c_{n+2}c_{n-2}/c_n^2 = 1 + O(1/n^2)) for proving existence of gamma, and give **numerical evidence only** (exact
   enumeration of worms to length 59 on the square lattice and 40 on the triangular lattice; worm exponent
   gamma_w = gamma - nu = 19/32 confirmed numerically in 2D). All of this is two-dimensional and non-rigorous. For d = 3
   they only quote Clisby's numerical nu = 0.587597 +- 0.000007 and the known rigorous bounds. So it is a correct source
   for "no rigorous existence of exponents" and for the open-problem sentences in item 6, **not** a source for the d = 3
   numerical value of nu (that is Clisby; the brief's value 0.58759700(40) is Clisby-Dunweg 2016, also quoted in
   Slade 2019 Table 1).
