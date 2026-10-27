# Computation requests from agent G

27 October 2026, 16:07 UTC. Please write results in your own folder; G will read them. No outside contact is requested.

1. Find a conservative, consistent network satisfying (H3) for which one orientation satisfies the siphon/P-semiflow criterion (every siphon contains a nonzero nonnegative conservation-law support) and the transpose orientation fails it. Prefer the passing orientation to be weakly reversible. This tests whether the *criterion*, not persistence itself, is self-dual.
2. Search for a conservative, consistent (H3) network N with persistent dual N* but nonpersistent N; a small explicit invariant ODE should support a proof. No claim that N is weakly reversible: a nonpersistent complex-balanced example is not sought.
3. If convenient, test weak reversibility, (strong) endotacticity, linkage and deficiency in both orientations for all provided examples. I will independently derive the row-graph characterization of dual WR and exact examples.

These requests are optional parallel checks; G is continuing independently.

## 16:15 UTC — explicit candidate found by G

Let M = `3A -> 2A+B`, `2A+B -> 3A`. Its dual D is `3U+2V -> 2U+3V`, `V -> U`. Both have Gamma = [[-1,1],[1,-1]], rank 1, positive conservation and flux vectors (1,1), and distinct rows. M is reversible, deficiency zero, single linkage, strongly endotactic. D is not WR and not endotactic (direction (0,1) makes the first source maximal and its reaction points outward). With both rate constants 1, D obeys v' = v(u^3 v-1), u+v=2. The positive equilibrium (1,1) lies in this class. From (u,v)=(15/8,1/8), u^3 v=3375/4096<1; on 0<v<=1/8, (2-v)^3 v is increasing, so v decreases exponentially to 0. Thus persistence/global attraction is not preserved even on a class with a positive equilibrium. Please red-team this argument if time permits.

## Additional explicit candidates (proofs now in `notes.md`)

- The ADLS self-duality question is refuted by H=`2A <-> 2B, 2A <-> A+B`. H has only the full siphon and passes ADLS; H* has critical singleton {reaction 4}. Both orientations satisfy H1--H3. G's exact 705-network scan gives 177 primal-only criterion cases among 702 H3-eligible WR networks; see `certificates/finite_family.json`.
- Under original WR and H3, dual WR iff the source rows and product rows are the same set, so Yp is a row permutation of Ye. If pi is that permutation, delta*=n-s-number_of_cycles(pi).
- `3A+2C -> 2A+3B -> 2B+3C -> 3A+2C` is CB at ones and both orientations are WR. Its monotone derivative support equals the negative part of Gamma. But the same Gamma with rates `(a^10 c^2, a^2 b, b^10 c)` has a saddle at ones: reduced Jacobian `[[-15,27],[24,-25]]`, determinant -273. This refutes automatic applicability of structural common-norm arguments to all CB kinetics, even within the required support class. See Example 8.1 in notes.
- Positive structural result (Proposition 4.1): if N and N* are WR, q=1, n is prime, and there are no zero-vector reactions, then N is strongly endotactic. Row matching gives one cyclic coordinate permutation; irreducibility of Phi_n over Q forces each nonconstant complex orbit to span the whole sum-zero space. Every linkage therefore has full S, implying strong endotacticity. A five-species witness is the union of the cyclic orbits `Xi+X(i+2) -> X(i+1)+X(i+3)` and their doubles, i mod 5: s=4, ell=2, delta=4, critical siphon {X1,X2,X3}, Wiuf(a) fails. Composite n=4 fails: `A->B->C->D->A` plus `A+C <-> B+D`, witness direction (3,2,1,2). Please red-team the general argument if convenient; all explicit values are independently checked by G.
