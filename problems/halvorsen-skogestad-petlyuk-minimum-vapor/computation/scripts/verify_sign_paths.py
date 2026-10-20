"""Certificate 2 (pure combinatorics, exact): no sign-consistent light->heavy path exists in the sign graph of any
basic configuration, N = 3..6, in the STRICT (exact infinite-stage) model.

Sign graph at theta_j: nodes = streams; the top section of column c joins node c (lower) to dist(c) (upper),
the bottom section joins res(c) (lower) to c (upper).  A traversal upward needs Psi > 0, downward Psi < 0, where
A_c = Psi_top(c)(theta_j), B_c = Psi_bot(c)(theta_j).  A path is pruned as soon as its sign requirements violate
one of the necessary conditions valid at every feasible point (proved analytically, see REPORT.md):
  (S3) k_c < j  => A_c >= 0 ;   l_c > j+1 => B_c >= 0
  (S4) not (A_c < 0 and B_c < 0)   [for mixed columns; for non-mixed columns S3 applies]
  (S5) strict model, j interior for c (l_c <= j <= k_c - 1) => not (A_c > 0 and B_c > 0)
  (R)  main column: A = B >= 0, and A = B = 0 if j is interior for it (strict model)
If no path survives, the cut identity with U = {nodes reachable from the light nodes} proves
R >= P_j - (1-q) at every feasible point (this is a machine check of the general invariant argument).
The same enumeration without S5 and without the interior part of (R) (= TG literal 'permissive' model) is
reported for information (expected: N=4 -> 1 path, family K = basic configuration ABCD->ABC/BCD,
ABC->AB/C, BCD->B/CD).

Environment: Python 3 standard library + hsmodel.py (enumeration).  Runtime ~10 s.  Deterministic.
Exit code non-zero if a strict-model path is found.
"""
import sys

from hsmodel import basic_configurations


def paths(cfg, j, strict):
    N = cfg['N']
    sp = cfg['splits']
    light = lambda s: s[1] <= j
    heavy = lambda s: s[0] >= j + 1
    moves = {}
    for c, (k, l) in sp.items():
        dist, res = (c[0], k), (l, c[1])
        moves.setdefault(c, []).append((dist, (c, 'A', 1)))
        moves.setdefault(dist, []).append((c, (c, 'A', -1)))
        moves.setdefault(res, []).append((c, (c, 'B', 1)))
        moves.setdefault(c, []).append((res, (c, 'B', -1)))

    def ok(reqs):
        sg = {}
        for (c, w, s) in reqs:
            if sg.setdefault((c, w), s) != s:
                return False
        for (c, w), s in sg.items():
            k, l = sp[c]
            if w == 'A' and s < 0 and k < j:
                return False
            if w == 'B' and s < 0 and l > j + 1:
                return False
            if sg.get((c, 'A')) == -1 and sg.get((c, 'B')) == -1:
                return False
            interior = l <= j <= k - 1
            if strict and interior and sg.get((c, 'A')) == 1 and sg.get((c, 'B')) == 1:
                return False
            if c == (1, N) and (s < 0 or (strict and interior)):
                return False
        return True
    found = []

    def dfs(x, seen, reqs):
        for (y, req) in moves.get(x, []):
            if y in seen or not ok(reqs + [req]):
                continue
            if heavy(y):
                found.append(reqs + [req])
            elif not light(y):
                dfs(y, seen | {y}, reqs + [req])
    for x in list(moves):
        if light(x):
            dfs(x, {x}, [])
    return found


if __name__ == '__main__':
    fail = False
    for N in (3, 4, 5, 6):
        B = basic_configurations(N)
        ns = sum(len(paths(c, j, True)) for c in B for j in range(1, N))
        perm = [(b, j) for b, c in enumerate(B) for j in range(1, N) for _ in paths(c, j, False)]
        status = 'PASS' if ns == 0 else 'FAIL'
        fail |= ns != 0
        print(f'{status}: N={N}: {len(B)} basic configurations, strict-model consistent paths = {ns}; '
              f'(info) permissive-model paths = {len(perm)} in {len(set(b for b, _ in perm))} configurations',
              flush=True)
        if N == 4:
            fams = sorted(set(b for b, _ in perm))
            names = [' ; '.join(f"{s}->{(s[0], B[b]['splits'][s][0])}/{(B[b]['splits'][s][1], s[1])}"
                                for s in B[b]['streams']) for b in fams]
            expectK = len(fams) == 1 and set(B[fams[0]]['streams']) == {(1, 4), (1, 3), (2, 4), (1, 2), (3, 4)} \
                and B[fams[0]]['splits'][(1, 3)] == (2, 3) and B[fams[0]]['splits'][(2, 4)] == (2, 3)
            print(f"{'PASS' if expectK else 'FAIL'}: N=4 permissive-model exceptions = family K only: {names}")
            fail |= not expectK
    sys.exit(1 if fail else 0)
