#!/usr/bin/env python3
"""Targeted corruption of a COPY of the K = 18 SAWCERT2 file (in place, few bytes), with ground truth.
Usage: patch_k18.py <copy.cert2> P_minus_1 | P_restore <P> | child <state> <target> | child_restore <state> <dir> <old>
'child' redirects the first allowed direction of <state> to <target> after checking, by brute force over all 48
signed permutations, that condition (1) fails for that pair; prints the direction and the old pointer."""
import itertools, struct, sys
f = open(sys.argv[1], 'r+b'); h = f.read(40)
d, K = struct.unpack_from('<II', h, 8); N, P, Q = struct.unpack_from('<QQQ', h, 16)
op = sys.argv[2]
def code(i):
    f.seek(40 + 8 * i); return struct.unpack('<Q', f.read(8))[0]
def pos(c):
    p = [(0,) * d]
    for i in range((c >> 57) & 63):
        t = (c >> (3 * i)) & 7; q = list(p[-1]); q[t >> 1] -= 1 if t % 2 == 0 else -1; p.append(tuple(q))
    return p
def unit(j):
    u = [0] * d; u[j >> 1] = 1 if j % 2 == 0 else -1; return tuple(u)
if op == 'P_minus_1':
    f.seek(24); f.write(struct.pack('<Q', P - 1)); print('P', P, '->', P - 1)
elif op == 'P_restore':
    f.seek(24); f.write(struct.pack('<Q', int(sys.argv[3]))); print('P restored')
elif op == 'child':
    s, t = int(sys.argv[3]), int(sys.argv[4])
    ps, pt = pos(code(s)), pos(code(t))
    for j in range(2 * d):
        if unit(j) in ps[1:]:
            continue
        e = unit(j); W = {tuple(p[i] - e[i] for i in range(d)) for p in ps}
        ok = any(all(tuple(sg[i] * q[pm[i]] for i in range(d)) in W for q in pt[1:])
                 for pm in itertools.permutations(range(d)) for sg in itertools.product((1, -1), repeat=d))
        off = 40 + 16 * N + 4 * (2 * d * s + j)
        f.seek(off); old = struct.unpack('<I', f.read(4))[0]
        print('state', s, 'dir', j, 'old child', old, 'new child', t, 'condition (1) holds for new pair:', ok)
        if not ok:
            f.seek(off); f.write(struct.pack('<I', t)); print('patched'); break
elif op == 'child_restore':
    s, j, old = map(int, sys.argv[3:6])
    f.seek(40 + 16 * N + 4 * (2 * d * s + j)); f.write(struct.pack('<I', old)); print('restored')
f.close()
