#!/usr/bin/env python3
"""Certificate: the counterexample is robust -- it holds on a set of parameters with nonempty interior in R^13.

Box: k_i in k*_i (1 +- eps) (k* = (1,1,1,1,300,1,3,1,1,100)), and steady-state coordinates
y = (x1, x2, x6) in y* (1 +- eps), y* = (1/3, 59, 1/4).  Let U = {(k, Phi_k(y)) : (k, y) in the open box},
where Phi_k(y) = W phi_k(y) are the totals of the class whose steady state is phi_k(y) (source, Prop. 3.6).
  E1  det D_y Phi_k(y) != 0 on the box (ball arithmetic) -> (k, y) -> (k, Phi_k(y)) is a local diffeomorphism,
      so U is open in R^13 (and nonempty);
  E2  for every (k, y) in the box all eigenvalues of the reduced Jacobian at phi_k(y) have negative real part:
      Gershgorin discs of R^-1 A R (R an approximate eigenvector matrix, A the ball enclosure of the Jacobian)
      lie in the open left half-plane -> the steady state is locally asymptotically stable;
  E3  enclosure of the totals Phi_k(y) over the box (ball arithmetic) -> interval parameters for CAPD;
  E4  po_proof_box (CAPD; parameters carried as constant state variables) with all 13 parameters as intervals:
      for every (k, T) in the box hull, the class contains an orbitally asymptotically stable periodic orbit.
Hence for every (k, T) in the open set U: the steady state is LAS, and it is not globally attracting.

Environment: Python 3.12, python-flint 0.9, CAPD 6.1 (see certify_orbits.py). Runtime: about 1 minute.
Usage: python certify_open_set.py [eps]   (default 1/10^12)
"""
import json
import math
import os
import subprocess
import sys
from fractions import Fraction

from flint import arb, arb_mat, fmpq, ctx

ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
eps = Fraction(sys.argv[1]) if len(sys.argv) > 1 else Fraction(1, 10 ** 12)
RADII = {'stable': '5e-8', 'existence': '1e-8'}
FAILS = []


def claim(name, ok, detail=''):
    print('%-4s %s  %s' % ('PASS' if ok else 'FAIL', name, detail), flush=True)
    if not ok:
        FAILS.append(name)


def box(q):
    q = Fraction(q)
    lo, hi = q * (1 - eps), q * (1 + eps)
    return arb(fmpq(lo.numerator, lo.denominator)).union(arb(fmpq(hi.numerator, hi.denominator)))

kstar = [1, 1, 1, 1, 300, 1, 3, 1, 1, 100]
K = [box(v) for v in kstar]
x1, x2, x6 = box(Fraction(1, 3)), box(59), box(Fraction(1, 4))
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = K


def phi(x1, x2, x6):
    J = k1 * k3 / (k2 + k3) * x1 * x2
    return [x1, x2, k1 / (k2 + k3) * x1 * x2, J / k4, (k6 + k7) * J / (k7 * k5 * x6), x6, J / k7,
            (k9 + k10) * J / (k10 * k8 * x6), J / k10]

x = phi(x1, x2, x6)
# totals and their Jacobian with respect to y = (x1, x2, x6) (explicit formulas; source, Prop. 3.6)
c3 = k1 / (k2 + k3); c4 = k1 * k3 / ((k2 + k3) * k4); c5 = k1 * k3 * (k6 + k7) / ((k2 + k3) * k5 * k7)
c7 = k1 * k3 / ((k2 + k3) * k7); c8 = k1 * k3 * (k9 + k10) / ((k2 + k3) * k8 * k10); c9 = k1 * k3 / ((k2 + k3) * k10)
# Ktot = x2 + (c3+c4) x1 x2 ; Ptot = x6 + (c7+c9) x1 x2 ; Stot = x1 + (c3+c4+c7+c9) x1 x2 + (c5+c8) x1 x2 / x6
a = c3 + c4; bP = c7 + c9; s1 = c3 + c4 + c7 + c9; s2 = c5 + c8
T = [x2 + a * x1 * x2, x6 + bP * x1 * x2, x1 + s1 * x1 * x2 + s2 * x1 * x2 / x6]
DT = arb_mat([[a * x2, 1 + a * x1, arb(0)],
              [bP * x2, bP * x1, arb(1)],
              [1 + s1 * x2 + s2 * x2 / x6, s1 * x1 + s2 * x1 / x6, -s2 * x1 * x2 / (x6 * x6)]])
dDT = DT.det()
claim('E1 det D_y Phi_k(y) != 0 on the box -> U is open', dDT > 0 or dDT < 0, 'det in %s' % dDT.str(8))

# reduced Jacobian A = P Dg(x) E (z = x3, x4, x5, x7, x8, x9)
REACTIONS = [((0, 1), (2,)), ((2,), (0, 1)), ((2,), (3,)), ((3,), (4, 1)), ((4, 5), (6,)),
             ((6,), (4, 5)), ((6,), (7, 5)), ((7, 5), (8,)), ((8,), (7, 5)), ((8,), (0, 5))]
N = [[0] * 10 for _ in range(9)]
for j, (re_, pr) in enumerate(REACTIONS):
    for i in re_:
        N[i][j] -= 1
    for i in pr:
        N[i][j] += 1
ZIDX = [2, 3, 4, 6, 7, 8]
E = [[0] * 6 for _ in range(9)]
for j, i in enumerate(ZIDX):
    E[i][j] = 1
E[0] = [-1] * 6; E[1] = [-1, -1, 0, 0, 0, 0]; E[5] = [0, 0, 0, -1, 0, -1]
D = [[arb(0)] * 9 for _ in range(10)]
for j, (re_, pr) in enumerate(REACTIONS):
    for a_ in re_:
        v = K[j]
        for b_ in re_:
            if b_ != a_:
                v = v * x[b_]
        D[j][a_] = D[j][a_] + v
Dg = [[sum((N[i][j] * D[j][l] for j in range(10)), arb(0)) for l in range(9)] for i in range(9)]
A = arb_mat([[sum((Dg[ZIDX[r]][l] * E[l][c] for l in range(9)), arb(0)) for c in range(6)] for r in range(6)])
import numpy as np
from flint import acb, acb_mat
Amid = np.array([[float(A[i, j].mid()) for j in range(6)] for i in range(6)])
w, Rn = np.linalg.eig(Amid)
R = acb_mat([[acb(float(Rn[i, j].real), float(Rn[i, j].imag)) for j in range(6)] for i in range(6)])
Rinv = R.inv()                                      # certified enclosure of the inverse of the exact double matrix R
Mg = Rinv * acb_mat([[acb(A[i, j]) for j in range(6)] for i in range(6)]) * R
ok2 = True; worst = None
for i in range(6):
    radius = sum((abs(Mg[i, j]) for j in range(6) if j != i), arb(0)) + Mg[i, i].rad()
    right_edge = Mg[i, i].real + radius
    ok2 = ok2 and (right_edge < 0)
    worst = right_edge if worst is None or right_edge.upper() > worst.upper() else worst
claim('E2 Gershgorin discs of R^-1 A R in Re < 0 for every (k, y) in the box -> steady state LAS', ok2,
      'max right edge of the discs <= %s' % worst.upper().str(5))
claim('E3 totals enclosed', True, 'Ktot in %s, Ptot in %s, Stot in %s' % tuple(t.str(12) for t in T))


orbits = json.load(open(os.path.join(ROOT, 'certificates', 'orbits.json')))['orbits']
capd_config = os.environ.get('CAPD_CONFIG', os.path.expanduser('~/opt/capd-src/build/bin/capd-config'))
exe = '/tmp/po_proof_box_cms'
flags = subprocess.run([capd_config, '--cflags', '--libs'], capture_output=True, text=True).stdout.split()
r = subprocess.run(['g++', os.path.join(HERE, 'po_proof_box.cpp')] + flags + ['-o', exe], capture_output=True, text=True)
claim('E4.0 po_proof_box.cpp compiles', r.returncode == 0, r.stderr[-300:])


def outward(v):
    """exactly representable double interval containing the ball v"""
    lo = float(v.lower()); hi = float(v.upper())
    lo = math.nextafter(lo, -math.inf); hi = math.nextafter(hi, math.inf)   # one ulp of slack covers the conversion
    assert arb(lo) <= v.lower() and arb(hi) >= v.upper()
    return lo.hex(), hi.hex()

pfile = '/tmp/cms_params_openset.txt'
with open(pfile, 'w') as fh:
    for i in range(10):
        fh.write('k%d %s %s\n' % ((i + 1,) + outward(K[i])))
    for name, t in zip(['K', 'P', 'S'], T):
        fh.write('%s %s %s\n' % ((name,) + outward(t)))
log = []
for o in orbits:
    if o['label'] != 'stable_mu59':
        continue
    args = [exe, pfile, str(o['section_index']), o['section_value'], str(o['direction'])] + o['z_guess'] \
        + [RADII[o['mode']], str(o['order']), o['mode'], o['label'] + '_box']
    r = subprocess.run(args, capture_output=True, text=True, timeout=1200)
    log.append('$ ' + ' '.join(args[1:]) + '\n' + open(pfile).read() + r.stdout + r.stderr)
    lines = r.stdout.splitlines()
    need = ['Krawczyk', 'positive'] + (['spectral radius'] if o['mode'] == 'stable' else [])
    ok = r.returncode == 0 and all(any(l.startswith('PASS') and n in l for l in lines) for n in need) \
        and not any(l.startswith('FAIL') for l in lines)
    info = [l for l in lines if 'I - C DF' in l or 'spectral' in l and 'PASS' not in l]
    claim('E4 %s: CAPD proof for the whole parameter box' % o['label'], ok, ' | '.join(s.split(': ', 1)[1] for s in info))
with open(os.path.join(ROOT, 'certificates', 'open_set_log.txt'), 'w') as fh:
    fh.write('eps = %s, radii = %s\n' % (eps, RADII))
    fh.write('\n'.join(log))
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('ALL CLAIMS PASS (eps = %s)' % eps)
