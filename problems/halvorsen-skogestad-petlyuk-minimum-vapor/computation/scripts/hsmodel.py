"""Core definitions for the Halvorsen-Skogestad conjecture computations (Agent C).

Environment: Python 3 (numpy, mpmath, python-flint); no solver needed for this module.

Conventions (Tumbalam Gooty-Agrawal-Tawarmalani arXiv:2010.12113 sec. 4 = Nallasivam et al. 2016):
  components 1..N, alpha_1 > ... > alpha_N > 0, feed F = 1 (or flows F_p), liquid fraction q.
  A stream [i,j] contains components i..j.  A configuration (regular-column, "basic" in Shah-Agrawal's sense)
  is a set S of mixture streams containing [1,N]; the split of [i,j] in S is
      distillate [i,k], k = max{k' < j : [i,k'] in S or k' = i},
      residue    [l,j], l = min{l' > i : [l',j] in S or l' = j},
  and must satisfy l <= k+1 (no component lost).  Every stream other than [1,N] (including pure products)
  must be produced by at least one present stream.  A submixture produced by exactly one split carries a
  binary choice: thermal coupling ('TC') or heat exchanger ('HX': a condenser if it is a distillate, a reboiler
  if it is a residue).  A stream produced by two splits is a side draw (no heat exchanger).
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction

import mpmath as mp


# --------------------------------------------------------------------------------------------------
# Feed roots and the Petlyuk (FTC) value
# --------------------------------------------------------------------------------------------------

def underwood_lhs(theta, alpha, f):
    return sum(a * x / (a - theta) for a, x in zip(alpha, f))


def feed_roots(alpha, z, q, dps=50):
    """Common Underwood roots theta_1 > ... > theta_{N-1}, alpha_{j+1} < theta_j < alpha_j, of
    sum_i alpha_i z_i / (alpha_i - theta) = (1 - q) * sum(z).  Bisection in mpmath (dps digits)."""
    with mp.workdps(dps):
        al = [mp.mpf(a) for a in alpha]
        zz = [mp.mpf(x) for x in z]
        rhs = (1 - mp.mpf(q)) * sum(zz)
        roots = []
        for j in range(len(al) - 1):
            lo, hi = al[j + 1], al[j]
            # lhs is increasing on (lo, hi) from -inf to +inf
            a, b = lo, hi
            for _ in range(dps * 4):
                m = (a + b) / 2
                if underwood_lhs(m, al, zz) - rhs > 0:
                    b = m
                else:
                    a = m
            roots.append((a + b) / 2)
        return roots


def petlyuk_peaks(alpha, z, q, dps=50):
    """Peaks P_j = sum_{i<=j} alpha_i z_i/(alpha_i - theta_j), j = 1..N-1."""
    th = feed_roots(alpha, z, q, dps)
    with mp.workdps(dps):
        return [sum(mp.mpf(alpha[i]) * mp.mpf(z[i]) / (mp.mpf(alpha[i]) - th[j]) for i in range(j + 1))
                for j in range(len(alpha) - 1)]


def petlyuk_value(alpha, z, q, dps=50):
    """(V_top, total reboiler vapor) of the generalized Petlyuk / FTC arrangement with liquid products:
    V_top = max_j P_j, reboiler vapor = V_top - (1-q)*sum(z)."""
    P = petlyuk_peaks(alpha, z, q, dps)
    vt = max(P)
    return vt, vt - (1 - mp.mpf(q)) * sum(mp.mpf(x) for x in z)


# --------------------------------------------------------------------------------------------------
# Enumeration of configurations
# --------------------------------------------------------------------------------------------------

def all_submixtures(N):
    return [(i, j) for i in range(1, N + 1) for j in range(i + 1, N + 1) if (i, j) != (1, N)]


def split_of(S, i, j):
    """Split (k, l) of mixture [i,j] given the set S of present mixture streams (pure streams implicit)."""
    k = max(kk for kk in range(i, j) if kk == i or (i, kk) in S)
    l = min(ll for ll in range(i + 1, j + 1) if ll == j or (ll, j) in S)
    return k, l


def basic_configurations(N):
    """All valid stream sets (Shah-Agrawal basic configurations).  Returns list of dicts."""
    subs = all_submixtures(N)
    out = []
    for mask in range(1 << len(subs)):
        S = {(1, N)} | {subs[t] for t in range(len(subs)) if mask >> t & 1}
        ok = True
        splits = {}
        for (i, j) in S:
            k, l = split_of(S, i, j)
            if l > k + 1:
                ok = False
                break
            splits[(i, j)] = (k, l)
        if not ok:
            continue
        # producers of each stream
        prod = {}
        for (i, j), (k, l) in splits.items():
            prod.setdefault((i, k), []).append(('top', (i, j)))
            prod.setdefault((l, j), []).append(('bot', (i, j)))
        streams = list(S) + [(p, p) for p in range(1, N + 1)]
        if any(s != (1, N) and s not in prod for s in streams):
            continue
        if any(len([x for x in prod.get(s, []) if x[0] == 'top']) > 1 or
               len([x for x in prod.get(s, []) if x[0] == 'bot']) > 1 for s in streams):
            raise RuntimeError('non-basic producer structure; should not happen')
        out.append({'N': N, 'streams': sorted(S, key=lambda s: (-(s[1] - s[0]), s)),
                    'splits': splits, 'producers': prod})
    return out


def candidates(cfg):
    """Submixtures produced by exactly one split (candidates for TC/HX choice)."""
    N = cfg['N']
    return [s for s in cfg['streams'] if s != (1, N) and len(cfg['producers'][s]) == 1]


def all_configurations(N):
    """All (basic configuration, TC/HX assignment) pairs.  'hx' maps candidate stream -> 'TC' or 'HX'."""
    out = []
    for b, cfg in enumerate(basic_configurations(N)):
        cand = candidates(cfg)
        for bits in itertools.product(['HX', 'TC'], repeat=len(cand)):
            c = dict(cfg)
            c['basic_id'] = b
            c['hx'] = dict(zip(cand, bits))
            out.append(c)
    for t, c in enumerate(out):
        c['id'] = t
    return out


def config_name(c):
    """Human-readable name: splits, with transfer-stream types."""
    N = c['N']
    lab = lambda s: ''.join(chr(64 + p) for p in range(s[0], s[1] + 1))
    parts = []
    for s in c['streams']:
        k, l = c['splits'][s]
        parts.append(f"{lab(s)}->{lab((s[0], k))}/{lab((l, s[1]))}")
    tr = []
    for s in c['streams']:
        if s == (1, N):
            continue
        if s in c.get('hx', {}):
            tr.append(f"{lab(s)}:{c['hx'][s]}")
        else:
            tr.append(f"{lab(s)}:SD")
    return '; '.join(parts) + ' | ' + ', '.join(tr)


def is_ftc(c):
    N = c['N']
    return len(c['streams']) == N * (N - 1) // 2 and all(v == 'TC' for v in c.get('hx', {}).values())


def to_json(c):
    key = lambda s: f"{s[0]},{s[1]}"
    return {'id': c['id'], 'basic_id': c['basic_id'], 'N': c['N'], 'name': config_name(c),
            'streams': [list(s) for s in c['streams']],
            'splits': {key(s): list(v) for s, v in c['splits'].items()},
            'producers': {key(s): [[t, list(p)] for t, p in v] for s, v in c['producers'].items()},
            'transfer': {key(s): v for s, v in c['hx'].items()},
            'is_ftc': is_ftc(c)}


def from_json(d):
    unkey = lambda k: tuple(int(x) for x in k.split(','))
    c = {'id': d['id'], 'basic_id': d['basic_id'], 'N': d['N'],
         'streams': [tuple(s) for s in d['streams']],
         'splits': {unkey(k): tuple(v) for k, v in d['splits'].items()},
         'producers': {unkey(k): [(t, tuple(p)) for t, p in v] for k, v in d['producers'].items()},
         'hx': {unkey(k): v for k, v in d['transfer'].items()}}
    return c


if __name__ == '__main__':
    for N in (3, 4, 5):
        B = basic_configurations(N)
        A = all_configurations(N)
        print(N, len(B), len(A) - len(B), len(A))
