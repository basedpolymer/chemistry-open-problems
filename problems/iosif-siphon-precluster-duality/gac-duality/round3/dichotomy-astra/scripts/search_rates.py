"""Exact rate-specific corner search. Python 3.14, standard library only.

Expected runtime: up to 5 minutes with default bounds; measured at exit.
All algebra and sign certificates use Fraction. Earlier geometry and Newton
domination code is imported without modification. Candidate output is not a
complete persistence certificate until the edge dynamics is justified.
"""
from pathlib import Path
import sys
import re
import time
from fractions import Fraction as F
from itertools import combinations, product
from collections import defaultdict, Counter

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE / 'round3' / 'dichotomy' / 'scripts'))
sys.path.insert(0, str(BASE.parent / 'scripts'))
import crn
import d2cert as D
from d2lib import Net
from explore_homogeneous import from_rows


def fixed_repeller(N, v, rates, weights=None):
    off = [j for j in range(N.r) if j not in v['Zv']]
    if weights is None:
        slopes = {F(0)}
        for kind in ('linear', 'monomial'):
            pairs = defaultdict(lambda: [F(0), F(0)])
            for i, poly in v['loc'].items():
                w = tuple(N.Y[i][j] for j in off)
                g1, g2 = N.G[i][v['e1']], N.G[i][v['e2']]
                for (p, q), c in poly.items():
                    k1 = (p, q) if kind == 'linear' else (p-1, q)
                    k2 = (p, q) if kind == 'linear' else (p, q-1)
                    pairs[k1, w][0] += rates[i]*c*g1
                    pairs[k2, w][1] += rates[i]*c*g2
            for aa, bb in pairs.values():
                if bb and -aa/bb > 0:
                    slopes.add(-aa/bb)
        ss = sorted(slopes)
        slopes.update((ss[i]+ss[i+1])/2 for i in range(len(ss)-1))
        slopes.add(ss[-1]+1)
        weights = [(F(1), t) for t in sorted(slopes)] + [(F(0), F(1))]
    for lam, mu in weights:
        for kind in ('linear', 'monomial'):
            terms = defaultdict(F)
            for i, poly in v['loc'].items():
                w = tuple(N.Y[i][j] for j in off)
                g1, g2 = N.G[i][v['e1']], N.G[i][v['e2']]
                for (p, q), c in poly.items():
                    if kind == 'linear':
                        terms[(p, q), w] += rates[i]*c*(lam*g1+mu*g2)
                    else:
                        terms[(p-1, q), w] += rates[i]*c*lam*g1
                        terms[(p, q-1), w] += rates[i]*c*mu*g2
            simple = [(e, 1 if c > 0 else -1) for (e, w), c in terms.items() if c]
            if D.dominated(simple):
                return kind, lam, mu
    return None


def parse_old():
    path = BASE / 'round3' / 'dichotomy' / 'logs' / 'undecided_s2.txt'
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        text = line.split(' | ')[0]
        n = max(map(int, re.findall(r'X(\d+)', text)))
        def cx(t):
            v = [0]*n
            for term in t.strip().split(' + '):
                m = re.fullmatch(r'(\d*)X(\d+)', term)
                if m:
                    v[int(m[2])-1] = int(m[1] or '1')
            return tuple(v)
        rx = [tuple(map(cx, s.split(' -> '))) for s in text.split(', ')]
        yield crn.Network(n, rx)


def candidates(mode):
    if mode == 'old':
        yield from parse_old()
        return
    if mode == 'random':
        import random
        rng = random.Random(381207)
        for _ in range(20000):
            d = rng.randrange(4, 15)
            n = rng.choice([4,5,6])
            vecs = [(p,q,d-p-q) for p in range(d+1) for q in range(d-p+1)
                    if not p == q == d-p-q]
            rows = rng.sample(vecs, n)
            N = from_rows(rows)
            if crn.rank(N.Gamma) != 2:
                continue
            rays = N.conservation_rays()
            if all(any(v[i] > 0 for v in rays) for i in range(n)):
                yield N
        return
    maxd = int(mode)
    for d in range(2, maxd+1):
        vecs = [(p, q, d-p-q) for p in range(d+1) for q in range(d-p+1)
                if not p == q == d-p-q]
        for rows in combinations(vecs, 3):
            N = from_rows(rows)
            if crn.rank(N.Gamma) != 2:
                continue
            h = crn.nullspace(list(zip(*N.Gamma)), 3)[0]
            if all(x > 0 for x in h) or all(x < 0 for x in h):
                yield N


def main():
    start = time.perf_counter()
    print('Environment: Python standard library; exact rational certificates.', flush=True)
    print('Runtime: measured at exit; expected up to 5 minutes.', flush=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else 'old'
    stats = Counter()
    for raw in candidates(mode):
        stats['total'] += 1
        N = Net(raw.n, raw.reactions)
        verts = [v for v in D.vertices(N) if v['eq']]
        bad = [v for v in verts if D.vertex_repelling(N, v) is None]
        if not bad:
            continue
        stats['bad_corners'] += 1
        facets = D.facets(N)
        # Fixed-rate corner certificates, no claim on edge sufficiency yet.
        if any(f['siphon'] and f['type'] != 'IN' and f['fl'] == 'zero' for f in facets):
            continue
        for rates in product((1, 4), repeat=N.n):
            certs = [fixed_repeller(N, v, rates) for v in bad]
            if all(certs):
                stats['candidate'] += 1
                print('CANDIDATE', raw.reactions, 'rates', rates,
                      'vertices', [(v['Zv'], c) for v,c in zip(bad,certs)],
                      'facets', [(f['Z'], f['type'], f['fl'], f['mu']) for f in facets if f['siphon']], flush=True)
                break
        if stats['total'] % 100 == 0:
            print('COUNTS', dict(stats), flush=True)
    print('PASS: SEARCH completed exact finite candidate search.', dict(stats))
    print('Runtime seconds:', round(time.perf_counter()-start, 3))


if __name__ == '__main__':
    main()
