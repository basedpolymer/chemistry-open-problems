#!/usr/bin/env python3
"""Falsification tests (independent verifier).

Corrupt a valid certificate in many ways and run every checker given on the command line on each copy.  The
expected verdict of every copy is decided here, independently of the checkers:
  * INVALID  -- the copy violates a hypothesis of Lemmas 1-2 (proved below case by case, or by the brute-force
                re-check of this script: all 2^d d! signed permutations for condition (1), exact integers for the
                inequality, which is complete for SAWCERT2 since only the touched pairs can change);
  * VALID    -- the copy is still a correct certificate (e.g. a weight raised where the inequality has slack);
                a checker must accept it (rejecting would only be over-strict, but we report it);
  * NEUTRAL  -- a format irregularity that does not change the automaton (stray payload bits); either verdict
                is sound.
A checker passes iff it rejects every INVALID copy.  Accepting a VALID copy is also required of a correct checker.

Usage: tamper.py <valid SAWCERT1 or SAWCERT2 file> <tmpdir> <exact max ratio "num/den" of the min-weight
       automaton (sawcheck output), used for the SAWCERT1 header cases> <checker command> [...]
       a checker command is a string such as "~/verif_saw/sawcheck {} 4"  ({} = file name).
Python 3 standard library only.  Runtime: about a minute for K = 12 certificates.
"""
import itertools, os, struct, subprocess, sys
from fractions import Fraction

src, tmp, maxratio, checkers = sys.argv[1], sys.argv[2], Fraction(sys.argv[3]), sys.argv[4:]
os.makedirs(tmp, exist_ok=True)
raw = open(src, 'rb').read()
magic = raw[:8]
d, K = struct.unpack_from('<II', raw, 8)
N, P, Q = struct.unpack_from('<QQQ', raw, 16)
codes = list(struct.unpack_from('<%dQ' % N, raw, 40))
V = list(struct.unpack_from('<%dQ' % N, raw, 40 + 8 * N))
explicit = magic == b'SAWCERT2'
ch = list(struct.unpack_from('<%dI' % (2 * d * N), raw, 40 + 16 * N)) if explicit else None
B = 3 if d <= 4 else 4
NONE = 0xFFFFFFFF
GROUP = [(perm, signs) for perm in itertools.permutations(range(d)) for signs in itertools.product((1, -1), repeat=d)]


def build(magic=magic, d=d, K=K, P=P, Q=Q, codes=codes, V=V, ch=ch):
    out = magic + struct.pack('<II', d, K) + struct.pack('<QQQ', len(codes), P, Q)
    out += struct.pack('<%dQ' % len(codes), *codes) + struct.pack('<%dQ' % len(V), *V)
    if ch is not None:
        out += struct.pack('<%dI' % len(ch), *ch)
    return out


def length(c):
    return (c >> 57) & 63


def positions(c):
    p = [(0,) * d]
    for i in range(length(c)):
        t = (c >> (B * i)) & ((1 << B) - 1)
        a, s = t >> 1, (1 if t % 2 == 0 else -1)
        q = list(p[-1]); q[a] -= s; p.append(tuple(q))
    return p


def unit(j):
    u = [0] * d; u[j >> 1] = 1 if j % 2 == 0 else -1; return tuple(u)


POS = [positions(c) for c in codes]


def allowed(s):
    F = set(POS[s][1:])
    return [j for j in range(2 * d) if unit(j) not in F]


def incl_ok(s, j, t):
    e = unit(j)
    W = {tuple(p[i] - e[i] for i in range(d)) for p in POS[s]}
    Ft = POS[t][1:]
    return any(all(tuple(signs[i] * q[perm[i]] for i in range(d)) in W for q in Ft) for perm, signs in GROUP)


ALLOWED = [allowed(s) for s in range(N)]


def cert2_truth(Vx, chx, Px=P, Qx=Q, touched=()):
    """exact re-check of an explicit automaton whose untouched pairs are known to satisfy (1)"""
    if min(Vx) < 1 or max(Vx) > 1 << 60 or Px == 0:
        return 'INVALID'
    for (s, j) in touched:
        t = chx[2 * d * s + j]
        if j in ALLOWED[s] and (t >= N or not incl_ok(s, j, t)):
            return 'INVALID'
    for s in range(N):
        tot = 0
        for j in ALLOWED[s]:
            t = chx[2 * d * s + j]
            if t >= N:
                return 'INVALID'
            tot += Vx[t]
        if Qx * tot > Px * Vx[s]:
            return 'INVALID'
    return 'VALID'


cases = [('untouched', build(), 'VALID')]
lam = Fraction(P, Q)
# header: the min-weight child choice over ALL self-avoiding paths through W minimises every row sum, so
# maxratio is the smallest ratio any valid choice of children can reach with these weights.
cases.append(('P_minus_1', build(P=P - 1), 'INVALID' if Fraction(P - 1, Q) < maxratio else 'VALID'))
cases.append(('Q_plus_1', build(Q=Q + 1), 'INVALID' if Fraction(P, Q + 1) < maxratio else 'VALID'))
cases.append(('P_equals_zero', build(P=0), 'INVALID'))
idx1 = codes.index((1 << 63) | (1 << 57))          # one-step walk
V2 = V[:]; V2[idx1] = V[idx1] * 11 // 10
if explicit:
    cases.append(('one_step_weight_plus_10pct', build(V=V2), cert2_truth(V2, ch)))
else:   # the start state's only possible children are the start state and the one-step walk
    r = Fraction(2 * d * min(V2[0], V2[idx1]), V2[0])
    cases.append(('one_step_weight_plus_10pct', build(V=V2), 'INVALID' if r > lam else 'VALID?'))
V2 = V[:]; V2[0] = V[0] * 9 // 10
if explicit:
    cases.append(('start_weight_minus_10pct', build(V=V2), cert2_truth(V2, ch)))
else:
    r = Fraction(2 * d * min(V2[0], V2[idx1]), V2[0])
    cases.append(('start_weight_minus_10pct', build(V=V2), 'INVALID' if r > lam else 'VALID?'))
V2 = V[:]; V2[N // 2] = 0
cases.append(('zero_weight', build(V=V2), 'INVALID'))
V2 = V[:]; V2[N // 3] = 1 << 62
cases.append(('weight_above_2^60', build(V=V2), 'INVALID'))
cases.append(('magic_corrupted', build(magic=b'SAWCERT9'), 'INVALID'))
cases.append(('truncated_by_one_byte', build()[:-1], 'INVALID'))
c2 = codes[:]; c2[5], c2[6] = c2[6], c2[5]
cases.append(('two_codes_swapped(order)', build(codes=c2), 'INVALID'))   # violates the format (C0)
c2 = codes[:]; L = length(c2[N - 1])
c2[N - 1] |= 1 << (B * L)   # a bit above the payload: same walk, irregular code
cases.append(('stray_payload_bit_last_code', build(codes=c2), 'NEUTRAL'))
if d == 6:
    # same 4-bit encoding for d = 5; the K = 10 states use at most 5 axes, so this is a valid d = 5 certificate
    # with 10 instead of 12 directions (the bound mu_5 <= P/Q is true and weak)
    cases.append(('dimension_6_to_5', build(d=5), 'VALID' if K <= 10 else 'VALID?'))
else:
    cases.append(('dimension_plus_1', build(d=d + 1), 'INVALID'))
if not explicit:
    cases.append(('one_step_state_removed', build(codes=codes[:idx1] + codes[idx1 + 1:], V=V[:idx1] + V[idx1 + 1:]),
                  'INVALID'))   # start state: only itself as child, ratio 2d > P/Q
    cases.append(('start_state_removed', build(codes=codes[1:], V=V[1:]), 'INVALID'))
else:
    s_long = max(range(N), key=lambda s: length(codes[s]))

    def find(pred):
        for s in range(N // 2, N):
            for j in ALLOWED[s]:
                if pred(s, j, ch[2 * d * s + j]):
                    return s, j
        raise RuntimeError('no suitable pair')

    s, j = find(lambda s, j, t: t != s_long and not incl_ok(s, j, s_long))
    ch2 = ch[:]; ch2[2 * d * s + j] = s_long
    cases.append(('child_redirected_to_longest_walk', build(ch=ch2), cert2_truth(V, ch2, touched=[(s, j)])))
    s, j = find(lambda s, j, t: True)
    ch2 = ch[:]; ch2[2 * d * s + j] = NONE
    cases.append(('child_set_to_none', build(ch=ch2), cert2_truth(V, ch2, touched=[(s, j)])))
    ch2 = ch[:]; ch2[2 * d * s + j] = N
    cases.append(('child_out_of_range', build(ch=ch2), cert2_truth(V, ch2, touched=[(s, j)])))
    s, j = find(lambda s, j, t: t != 0 and Q * (sum(V[ch[2 * d * s + i]] for i in ALLOWED[s]) - V[t] + V[0]) > P * V[s])
    ch2 = ch[:]; ch2[2 * d * s + j] = 0     # (1) holds (F(start) is empty): only the inequality can catch it
    cases.append(('child_to_start_state(incl_ok)', build(ch=ch2), cert2_truth(V, ch2, touched=[(s, j)])))
    done = False
    for s in range(N // 2, N):
        for a, b in itertools.combinations(ALLOWED[s], 2):
            ta, tb = ch[2 * d * s + a], ch[2 * d * s + b]
            if ta != tb and not (incl_ok(s, a, tb) and incl_ok(s, b, ta)):
                ch2 = ch[:]; ch2[2 * d * s + a], ch2[2 * d * s + b] = tb, ta; done = True; break
        if done:
            break
    cases.append(('children_of_two_directions_swapped', build(ch=ch2), cert2_truth(V, ch2, touched=[(s, a), (s, b)])))
    ch2 = ch[:]
    for j in range(2 * d):
        ch2[2 * d * idx1 + j] = NONE
    cases.append(('one_step_state_children_removed', build(ch=ch2), cert2_truth(V, ch2, touched=[(idx1, j) for j in range(2 * d)])))
    # a pointer placed on a FORBIDDEN direction is irrelevant (that step is never taken): must be accepted
    s = next(s for s in range(N) if len(ALLOWED[s]) < 2 * d)
    jf = next(j for j in range(2 * d) if j not in ALLOWED[s])
    ch2 = ch[:]; ch2[2 * d * s + jf] = s_long
    cases.append(('pointer_on_forbidden_direction', build(ch=ch2), 'VALID'))
    bad = sum(1 for s in range(0, N, max(1, N // 3000)) for j in ALLOWED[s] if not incl_ok(s, j, ch[2 * d * s + j]))
    print('brute-force check of (1) on %d sampled states of the untouched file: %d violations'
          % (len(range(0, N, max(1, N // 3000))), bad))

ok_all = True
print('certificate %s: %s d=%d K=%d N=%d P/Q=%d/%d' % (os.path.basename(src), magic.decode(), d, K, N, P, Q))
for name, data, truth in cases:
    path = os.path.join(tmp, 'tamper_' + name.replace('(', '_').replace(')', '') + ('.cert2' if explicit else '.cert'))
    open(path, 'wb').write(data)
    res = []
    for cmd in checkers:
        r = subprocess.run(cmd.replace('{}', path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rejected = r.returncode != 0
        if truth == 'INVALID':
            good = rejected
        elif truth.startswith('VALID'):
            good = not rejected
        else:
            good = True
        ok_all &= good or truth != 'INVALID'
        res.append(('rejected' if rejected else 'accepted') + ('' if good else ' <-- UNEXPECTED'))
    print('  %-38s truth=%-8s %s' % (name, truth, ' | '.join(res)))
    os.remove(path)
print('PASS (every INVALID copy rejected by every checker)' if ok_all else 'FAIL')
sys.exit(0 if ok_all else 1)
