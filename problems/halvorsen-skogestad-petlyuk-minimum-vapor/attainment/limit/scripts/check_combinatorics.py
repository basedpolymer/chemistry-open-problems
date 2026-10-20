#!/usr/bin/env python3
"""check_combinatorics.py -- exhaustive exact check of the case analyses of ../limit.md (Agent L).

The proofs of Proposition 3.3 (c), Lemma 3.4, Lemma 4.1 and Proposition 4.5 are case analyses over the six section
types and the forcing pattern.  This script re-does them mechanically, in exact integer arithmetic, from the limit
supports of (L2) alone, for EVERY N = 3..Nmax and EVERY pattern pi in {B,T}^{N-1}:
  K1  every trace mode is closed at a node where its component is absent (lighter traces at the bottom end,
      heavier ones at the top end), and every component present at that node lies on the correct side of k;
  K2  every section attached to a node where k is absent has zero limit flow of k (Lemma 3.4), and every such node
      has first index > k or second index < k (exactly one);
  K3  counts (Lemma 4.1): per column and component, #trace rows = #trace coordinates; #major rows = #major coordinates;
      totals N(N-1)^2;
  K4  triangularity (Proposition 4.5): for every trace section whose two ends are trace nodes, h(open) = h(closed) - 1
      for a lighter trace and + 1 for a heavier trace, both ends in the same class; every non-product trace node is the
      closed end of at least one section; hence A_k is triangular with nonzero diagonal (checked by building the
      directed graph closed -> open and verifying it is acyclic).
Prints PASS/FAIL per N and exits non-zero on failure.  This is a check of the bookkeeping, not of the analysis.

Usage: python check_combinatorics.py [Nmax]      (default 12; 2^(N-1) patterns per N)
Environment: Python 3 (standard library only).  Runtime: about 1 minute for Nmax = 12.
"""
import sys
import itertools


def supports(N, pat):
    """Limit supports (L2), 1-based components; columns (i,m); returns w-supports of sections and node supports."""
    sec = {}      # ('T'|'B', i, m) -> (support set, kind 'T'|'B', type)
    for i in range(1, N + 1):
        for m in range(i + 1, N + 1):
            if m - i >= 2:
                sec[('T', i, m)] = (set(range(i, m)), 'T-nb')
                sec[('B', i, m)] = (set(range(i + 1, m + 1)), 'B-nb')
            else:
                p = i
                if pat[p - 1] == 'B':
                    sec[('T', i, m)] = ({p}, 'T-np')
                    sec[('B', i, m)] = ({p, p + 1}, 'B-f')
                else:
                    sec[('T', i, m)] = ({p, p + 1}, 'T-f')
                    sec[('B', i, m)] = ({p + 1}, 'B-np')
    node = {}
    for i in range(1, N + 1):
        for m in range(i + 1, N + 1):
            node[(i, m)] = set(range(i, m + 1))
    for p in range(1, N + 1):
        s = {p}
        if p >= 2 and pat[p - 2] == 'B':
            s.add(p - 1)
        if p <= N - 1 and pat[p - 1] == 'T':
            s.add(p + 1)
        node[(p, p)] = s
    return sec, node


def ends(s):
    kind, i, m = s
    if kind == 'T':
        return (i, m), (i, m - 1)          # bottom, top
    return (i + 1, m), (i, m)


def check(N, pat):
    sec, node = supports(N, pat)
    errs = []
    # product supports must equal the union of the supports of the two attached sections (all limit flows >= 0)
    for p in range(1, N + 1):
        u = set()
        if p >= 2:
            u |= sec[('B', p - 1, p)][0]
        if p <= N - 1:
            u |= sec[('T', p, p + 1)][0]
        if u != node[(p, p)]:
            errs.append(f'product support {p}')
    trace_rows = {}            # k -> list of (section, closed end, open end, lighter?)
    for s, (sup, typ) in sec.items():
        a, b = min(sup), max(sup)
        assert sup == set(range(a, b + 1))
        bot, top = ends(s)
        for k in range(1, N + 1):
            if k in sup:
                continue
            lighter = k < a
            c, o = (bot, top) if lighter else (top, bot)
            # K1
            if k in node[c]:
                errs.append(f'K1 closed end has k: {s} k={k}')
            if lighter and not all(j > k for j in node[c]):
                errs.append(f'K1 side (bottom) {s} k={k}')
            if not lighter and not all(j < k for j in node[c]):
                errs.append(f'K1 side (top) {s} k={k}')
            trace_rows.setdefault(k, []).append((s, c, o, lighter))
    # K2
    attached = {}
    for s in sec:
        for e in ends(s):
            attached.setdefault(e, []).append(s)
    for k in range(1, N + 1):
        for J, comp in node.items():
            if k in comp:
                continue
            if not ((J[0] > k) ^ (J[1] < k)):
                errs.append(f'K2 dichotomy {J} k={k}')
            for s in attached[J]:
                if k in sec[s][0]:
                    errs.append(f'K2 attached section carries k: {J} {s} k={k}')
    # K3 counts per column and component
    tot_rows = tot_coords = 0
    for i in range(1, N + 1):
        for m in range(i + 1, N + 1):
            rng = set(range(i, m + 1))
            for k in range(1, N + 1):
                rows = sum(1 for s in (('T', i, m), ('B', i, m)) if k not in sec[s][0])
                if m - i >= 2:
                    coords = 2 if k not in rng else (1 if k in (i, m) else 0)
                else:
                    p = i
                    if k not in rng:
                        coords = 2
                    elif pat[p - 1] == 'B':
                        coords = 1 if k == p + 1 else 0
                    else:
                        coords = 1 if k == p else 0
                if rows != coords:
                    errs.append(f'K3 trace count {(i, m)} k={k}: rows {rows} coords {coords}')
                tot_rows += rows
                tot_coords += coords
            if m - i >= 2:
                majrows = (m - i - 1) + (m - i - 1)      # T: r + middle; B: middle + s
                majcoords = (m - i) + (m - i - 2)
            else:
                majrows, majcoords = 1, 1
            if majrows != majcoords:
                errs.append(f'K3 major count {(i, m)}')
            tot_rows += majrows
            tot_coords += majcoords
    if tot_rows != N * (N - 1) ** 2 or tot_coords != N * (N - 1) ** 2:
        errs.append(f'K3 totals {tot_rows} {tot_coords}')
    # K4 triangularity / acyclicity and nonzero diagonal
    for k, lst in trace_rows.items():
        tnodes = {J for J, comp in node.items() if k not in comp}
        edges = {}
        closed_at = {J: 0 for J in tnodes}
        for (s, c, o, lighter) in lst:
            closed_at[c] = closed_at.get(c, 0) + 1
            if o in tnodes:
                h = lambda J: J[0] + J[1]
                if lighter and not (h(o) == h(c) - 1 and c[0] > k and o[0] > k):
                    errs.append(f'K4 light move {s} k={k}')
                if not lighter and not (h(o) == h(c) + 1 and c[1] < k and o[1] < k):
                    errs.append(f'K4 heavy move {s} k={k}')
                edges.setdefault(c, set()).add(o)
        for J in tnodes:
            if J[0] != J[1] and closed_at.get(J, 0) == 0:
                errs.append(f'K4 empty diagonal at {J} k={k}')
        # acyclicity (DFS)
        state = {}

        def dfs(u):
            state[u] = 1
            for v in edges.get(u, ()):
                if state.get(v) == 1:
                    return False
                if state.get(v) is None and not dfs(v):
                    return False
            state[u] = 2
            return True
        for u in tnodes:
            if state.get(u) is None and not dfs(u):
                errs.append(f'K4 cycle k={k}')
                break
    return errs


def main():
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    sys.setrecursionlimit(10000)
    allok = True
    for N in range(3, Nmax + 1):
        bad = 0
        first = None
        for pat in itertools.product('BT', repeat=N - 1):
            e = check(N, pat)
            if e:
                bad += 1
                first = first or (''.join(pat), e[:3])
        allok &= bad == 0
        print(f'N={N:2d}: {2 ** (N - 1):5d} patterns, {N * (N - 1) ** 2:5d} rows/unknowns each: '
              f'{"PASS" if bad == 0 else "FAIL " + str(bad) + " " + str(first)}')
        sys.stdout.flush()
    print('ALL PASS' if allok else 'SOME FAIL')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
