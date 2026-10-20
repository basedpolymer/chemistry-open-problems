"""Classification of the mixed sub-networks M_j (Agent G).

For every basic configuration (Shah-Agrawal / TG space, enumeration of computation/scripts/hsmodel.py, used
read-only) and every split index j = 1..N-1, the mixed sub-network M_j consists of the streams [i,m] with
i <= j < m.  We record:
  - branching columns (type d: distillate AND residue mixed),
  - merge nodes (mixed streams produced by two splits = mixed side draws),
  - tree  = no merge node  (Theorem G3 of general-N.md applies to every TC/HX variant),
  - optionally, the sign-relaxation status of Agent T (theory/scripts/signcheck.py, LP per TC/HX variant).
Environment: Python 3, scipy (only with --sign).  Runtime: N=5 seconds (no --sign), N=6 ~1 min.
Usage: python classify_mixed.py N [--sign]
"""
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'computation', 'scripts'))
from hsmodel import basic_configurations, all_configurations, config_name  # noqa: E402


def mixed_info(cfg, j):
    M = [s for s in cfg['streams'] if s[0] <= j < s[1]]
    branching = [s for s in M if cfg['splits'][s][0] >= j + 1 and cfg['splits'][s][1] <= j]
    merges = [s for s in M if s != (1, cfg['N']) and len(cfg['producers'][s]) == 2]
    return M, branching, merges


def ctype(cfg, s, j):
    k, l = cfg['splits'][s]
    dm, rm = k >= j + 1, l <= j
    return {(False, False): 'a', (False, True): 'b', (True, False): 'c', (True, True): 'd'}[(dm, rm)]


def diamonds_ok(cfg, j):
    """Hypothesis of Theorem G8: every merge node s is the residue of a type-(b) column c1 and the distillate
    of a type-(c) column c2, c1 and c2 are not merge nodes, and they are the distillate and residue of one
    type-(d) column."""
    M, br, mg = mixed_info(cfg, j)
    for s in mg:
        prods = cfg['producers'][s]
        c1 = [p for (t, p) in prods if t == 'bot'][0]   # s is residue of c1
        c2 = [p for (t, p) in prods if t == 'top'][0]   # s is distillate of c2
        if ctype(cfg, c1, j) != 'b' or ctype(cfg, c2, j) != 'c' or c1 in mg or c2 in mg:
            return False
        p1 = cfg['producers'][c1]
        p2 = cfg['producers'][c2]
        if len(p1) != 1 or len(p2) != 1:
            return False
        (t1, d1), (t2, d2) = p1[0], p2[0]
        if not (t1 == 'top' and t2 == 'bot' and d1 == d2 and ctype(cfg, d1, j) == 'd'):
            return False
    return True


def main():
    N = int(sys.argv[1])
    use_sign = '--sign' in sys.argv
    B = basic_configurations(N)
    stats = Counter()
    rows = []
    for b, cfg in enumerate(B):
        for j in range(1, N):
            M, br, mg = mixed_info(cfg, j)
            tree = not mg
            stats[('tree' if tree else 'merge', 'branch' if br else 'nobranch')] += 1
            if not tree and diamonds_ok(cfg, j):
                stats[('merge', 'simple-diamonds (G8)')] += 1
            rows.append((b, j, tree, len(br), len(mg)))
    print(f'N = {N}: {len(B)} basic configurations, {len(B) * (N - 1)} (basic, j) pairs')
    for k, v in sorted(stats.items()):
        print('  ', k, v)
    if use_sign:
        sys.path.insert(0, os.path.join(HERE, '..', '..', 'theory', 'scripts'))
        from signcheck import sign_hard  # noqa: E402
        A = all_configurations(N)
        hard = Counter()
        hard_nontree = []
        for c in A:
            for j in range(1, N):
                h, _ = sign_hard(c, j, c['hx'])
                if h:
                    M, br, mg = mixed_info(c, j)
                    hard[(j, 'tree' if not mg else ('G8' if diamonds_ok(c, j) else 'merge-open'))] += 1
                    if mg and not diamonds_ok(c, j):
                        hard_nontree.append((c['id'], c['basic_id'], j))
        print('sign-hard (configuration variant, j) pairs by tree/merge:')
        for k, v in sorted(hard.items()):
            print('  ', k, v)
        basics = sorted(set((b, j) for (_, b, j) in hard_nontree))
        print('sign-hard AND not covered by G3/G8, basic (b, j):', len(basics))
        for (b, j) in basics[:40]:
            ex = [i for (i, bb, jj) in hard_nontree if bb == b and jj == j]
            print(f'   basic {b} j={j}: {len(ex)} variants, e.g. {config_name(A[ex[0]])}')


if __name__ == '__main__':
    main()
