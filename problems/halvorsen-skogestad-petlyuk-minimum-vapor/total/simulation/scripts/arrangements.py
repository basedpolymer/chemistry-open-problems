"""
arrangements.py -- stage-by-stage networks of classical arrangements, built with stagenet.Network.

Agent P (Claude Opus 5.5, effort xhigh).  Components are 0-based (0 = lightest).  A column is given by its
label (i, m) and split (k, l): distillate [i, k], residue [l, m] (l <= k + 1).  Conventions (result 12,
Definition 1 and TG): a mixture stream produced by one column is transferred by a thermal coupling ('TC') or a
heat exchanger ('HX': partial condenser -> saturated vapor transfer; partial reboiler -> saturated liquid
transfer); a stream produced by two columns is received at the receiving column's feed stage through two thermal
couplings; a pure product produced once is withdrawn through a total condenser (distillate) or a partial reboiler
(residue), produced twice as a liquid side draw.  Every column has one feed stage (an equilibrium stage) between
a top and a bottom section of n stages.
"""
import numpy as np
from stagenet import Network

LET = 'ABCDEFGH'


def lab(i, m):
    return LET[i:m + 1]


def regular_network(alpha, z, q, columns, transfer=None, n=40, nT=None):
    """columns: dict {(i,m): (k,l)}; transfer: dict {label string: 'TC'|'HX'} for streams produced once.
    nT: optional dict {section name: stages} overriding n."""
    N = len(alpha)
    transfer = transfer or {}
    nT = nT or {}
    net = Network(alpha)
    # who produces what
    prod = {}
    for (i, m), (k, l) in columns.items():
        prod.setdefault((i, k), []).append(((i, m), 'D'))
        prod.setdefault((l, m), []).append(((i, m), 'B'))
    for (i, m) in columns:
        net.node('C' + lab(i, m))
    for p in range(N):
        pr = prod.get((p, p), [])
        kinds = {t for _, t in pr}
        if kinds == {'D'}:
            net.node('P' + LET[p], 'cooler')
        elif kinds == {'B'}:
            net.node('P' + LET[p], 'heater')
        else:
            net.node('P' + LET[p], 'adiabatic')
    for (i, m), (k, l) in columns.items():
        c = lab(i, m)
        # top section
        dst = (i, k)
        if k == i:
            top = 'P' + LET[i]
        else:
            once = len(prod[dst]) == 1
            if once and transfer.get(lab(*dst), 'TC') == 'HX':
                net.node('cd' + c, 'cooler')
                top = 'cd' + c
                net.stream('D' + c, 'cd' + c, 'V', 'C' + lab(*dst))
            else:
                top = 'C' + lab(*dst)
        net.section('T' + c, nT.get('T' + c, n), top, 'C' + c)
        dst = (l, m)
        if l == m:
            bot = 'P' + LET[m]
        else:
            once = len(prod[dst]) == 1
            if once and transfer.get(lab(*dst), 'TC') == 'HX':
                net.node('rb' + c, 'heater')
                bot = 'rb' + c
                net.stream('B' + c, 'rb' + c, 'L', 'C' + lab(*dst))
            else:
                bot = 'C' + lab(*dst)
        net.section('B' + c, nT.get('B' + c, n), 'C' + c, bot)
    net.feed('C' + lab(0, N - 1), 1.0, q, z)
    for p in range(N):
        net.product(LET[p], 'P' + LET[p], 'L', p)
    net.build()
    return net


def ftc_columns(N):
    return {(i, m): (m - 1, i + 1) for i in range(N) for m in range(i + 1, N)}


def recovery_targets(N, z, columns):
    """Net top flows of every column for a nominal split (components p<l up, p>k down, l..k half)."""
    feeds = {(0, N - 1): np.array(z, float)}
    order = sorted(columns, key=lambda c: -(c[1] - c[0]))
    D = {}
    for (i, m) in order:
        k, l = columns[(i, m)]
        f = feeds.get((i, m), np.zeros(N))
        top = np.zeros(N); bot = np.zeros(N)
        for p in range(i, m + 1):
            if p < l:
                top[p] = f[p]
            elif p > k:
                bot[p] = f[p]
            else:
                top[p] = 0.5 * f[p]; bot[p] = 0.5 * f[p]
        D[(i, m)] = top.sum()
        for dst, v in (((i, k), top), ((l, m), bot)):
            feeds[dst] = feeds.get(dst, np.zeros(N)) + v
    return D


def init_extra(net, N, z, columns):
    D = recovery_targets(N, z, columns)
    extra = []
    for (i, m), d in D.items():
        c = lab(i, m)
        extra.append(({'V:T' + c: 1.0, 'L:T' + c: -1.0}, float(d)))
    for p in range(N):
        extra.append(({'P:' + LET[p]: 1.0}, float(z[p])))
    return extra
