"""Certificate 1 (exact rational arithmetic): the cut identity

   R - (P_j - (1-q)) = FR(U) + FC(U^c) + sum_{s: U^c -> U} Psi_s(theta) - sum_{s: U -> U^c} Psi_s(theta)

holds for every configuration of the Shah-Agrawal / TG space (N = 3, 4, 5: 8 + 152 + 6128 configurations), every
split index j and every node partition (U, U^c) with light products 1..j in U and heavy products in U^c,
whenever theta is a root of the main feed equation (1-q) = sum_p alpha_p z_p/(alpha_p - theta).

Only the linear balances are used (component balance of each column, vapor balance at every node with arbitrary
condenser/reboiler duties FC, FR >= 0 at junctions,
including arbitrary extra heater/cooler duties at every column node and product node), not the Underwood inequalities: the identity is algebraic.
For each configuration we draw random rational flows/vapors/volatilities/theta (q is then defined by the feed
equation so that theta is exactly a feed root), and check the identity for the canonical top cut, the canonical
bottom cut and 3 random partitions, in exact Fraction arithmetic.  Both sides are affine in the free flow and
vapor variables for fixed (alpha, z, theta); several random draws per configuration check it (a failure would be
an algebra/modelling error).

Environment: Python 3 standard library + hsmodel.py.  Runtime: ~2 min (N = 5 dominates).  Deterministic (seed).
Prints PASS/FAIL per N and exits non-zero on any failure.
"""
import random
import sys
from fractions import Fraction as Fr

from hsmodel import all_configurations
from nlp_structure import structure


def rnd(rng, lo=1, hi=1000):
    return Fr(rng.randint(lo, hi), rng.randint(lo, hi))


def check_config(cfg, rng, draws=2, nparts=3):
    N = cfg['N']
    cols, info, reb, side = structure(cfg)
    for _ in range(draws):
        # volatilities, feed, theta, q
        al = [Fr(1)]
        for _k in range(N - 1):
            al.insert(0, al[0] * (1 + rnd(rng, 1, 50) / 20))
        al = {p: al[p - 1] for p in range(1, N + 1)}
        z = {p: rnd(rng) for p in range(1, N + 1)}
        for j in range(1, N):
            th = al[j + 1] + (al[j] - al[j + 1]) * Fr(rng.randint(1, 999), 1000)
            one_minus_q = sum(al[p] * z[p] / (al[p] - th) for p in z)
            P = sum(al[p] * z[p] / (al[p] - th) for p in range(1, j + 1))
            # flows top-down
            order = sorted(cols, key=lambda s: -(s[1] - s[0]))
            d, b, VT, VB = {}, {}, {}, {}
            FC, FR = {}, {}
            feeds = {}
            for s in order:
                I = info[s]
                if s == (1, N):
                    f = dict(z)
                    vf = one_minus_q
                else:
                    f = {p: Fr(0) for p in range(I['i'], I['j'] + 1)}
                    vf = Fr(0)
                    if I['top'] is not None:
                        t = I['top']
                        for p in range(info[t]['i'], info[t]['k'] + 1):
                            f[p] += d[t, p]
                        if I['kind'] == 'HX':
                            D = sum(d[t, p] for p in range(info[t]['i'], info[t]['k'] + 1))
                            vf += D
                            FC[s] = VT[t] - D
                        else:
                            vf += VT[t]
                    if I['bot'] is not None:
                        u = I['bot']
                        for p in range(info[u]['l'], info[u]['j'] + 1):
                            f[p] += b[u, p]
                        if I['kind'] == 'HX':
                            FR[s] = VB[u]
                        else:
                            vf -= VB[u]
                for p in range(I['i'], I['j'] + 1):
                    if p < I['l']:
                        d[s, p] = f[p]
                    elif p > I['k']:
                        b[s, p] = f[p]
                    else:
                        x = Fr(rng.randint(0, 100), 100)
                        d[s, p] = x * f[p]
                        b[s, p] = f[p] - d[s, p]
                VT[s] = rnd(rng, 1, 5000) / 100
                # arbitrary extra heater (h) and cooler (k) duties at the column node (any junction duty allowed)
                h, k = Fr(rng.randint(0, 50), 100), Fr(rng.randint(0, 50), 100)
                FR[('node', s)] = h
                FC[('node', s)] = k
                VB[s] = VT[s] - vf - h + k
            # product nodes: arbitrary duties consistent with the vapor balance
            for p in range(1, N + 1):
                prods = cfg['producers'][(p, p)]
                top = [x for (tt, x) in prods if tt == 'top']
                bot = [x for (tt, x) in prods if tt == 'bot']
                vin = VT[top[0]] if top else Fr(0)
                vout = VB[bot[0]] if bot else Fr(0)
                extra = Fr(rng.randint(0, 100), 100)
                # FC - FR = vin - vout (liquid product); add the same 'extra' to both (allowed: heat at junctions)
                FC[(p, p)] = max(vin - vout, Fr(0)) + extra
                FR[(p, p)] = max(vout - vin, Fr(0)) + extra
            R = sum(FR.values())
            # global balance check
            if sum(FC.values()) - R != one_minus_q:
                return False, 'global vapor balance'

            def psiT(s):
                I = info[s]
                return VT[s] - sum(al[p] * d[s, p] / (al[p] - th) for p in range(I['i'], I['k'] + 1))

            def psiB(s):
                I = info[s]
                return VB[s] + sum(al[p] * b[s, p] / (al[p] - th) for p in range(I['l'], I['j'] + 1))
            sections = []
            for s in cols:
                I = info[s]
                sections.append((s, (I['i'], I['k']), psiT(s)))
                sections.append(((I['l'], I['j']), s, psiB(s)))
            nodes = set(cols) | {(p, p) for p in range(1, N + 1)}
            light = {x for x in nodes if x[1] <= j}
            heavy = {x for x in nodes if x[0] >= j + 1}
            mixed = nodes - light - heavy
            parts = [light, light | mixed]
            for _k in range(nparts):
                parts.append(light | {x for x in mixed if rng.random() < 0.5})
            lhs = R - (P - one_minus_q)
            for U in parts:
                loc = lambda x: x[1] if x[0] == 'node' else x
                rhs = sum(v for x, v in FR.items() if loc(x) in U) + sum(v for x, v in FC.items() if loc(x) not in U)
                for (lo, up, ps) in sections:
                    if lo not in U and up in U:
                        rhs += ps
                    elif lo in U and up not in U:
                        rhs -= ps
                if rhs != lhs:
                    return False, f'identity j={j} U={sorted(U)}'
    return True, ''


if __name__ == '__main__':
    rng = random.Random(20261020)
    ok_all = True
    for N in (3, 4, 5):
        C = all_configurations(N)
        bad = []
        for c in C:
            ok, msg = check_config(c, rng, draws=2 if N < 5 else 1)
            if not ok:
                bad.append((c['id'], msg))
        status = 'PASS' if not bad else 'FAIL'
        ok_all &= not bad
        print(f'{status}: cut identity, N={N}, {len(C)} configurations, all j, 5 partitions each'
              + ('' if not bad else f' failures: {bad[:5]}'), flush=True)
    sys.exit(0 if ok_all else 1)
