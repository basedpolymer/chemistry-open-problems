#!/usr/bin/env python3
"""Certificate: explicit rational compatibility classes with a locally asymptotically stable steady state and
periodic orbits (mixed-mechanism network, rate constants k = (1,1,1,1,300,1,3,1,1,100)).

For each class in certificates/orbits.json (mu = 59: totals (236/3, 6527/1800, 341699/5400);
mu = 60: totals (80, 221/60, 5791/90)):
  D1  the steady state x* = phi(1/3, mu, 1/4) is rational, positive, and lies in the class (exact);
  D2  its stability (exact rational Hurwitz determinants of the reduced Jacobian):
        mu = 59: det H1..H5 > 0, b6 > 0  -> locally asymptotically stable;
        mu = 60: det H5 < 0              -> unstable (Routh: two eigenvalues with positive real part);
  D3  the vector field hard-coded in po_proof.cpp equals the reduced field P g(x(z)) (symbolic identity);
  D4  po_proof (CAPD 6.1 rigorous C^1 Poincare map + Krawczyk operator) proves the periodic orbits:
        mu = 59: an orbitally asymptotically stable orbit and a second periodic orbit;
        mu = 60: an orbitally asymptotically stable orbit.
Prints PASS/FAIL per claim, exit code non-zero on failure. CAPD output is saved to certificates/po_proof_log.txt.

Environment: Python 3.12 + sympy 1.14 + python-flint 0.9; CAPD 6.1 (github.com/CAPDGroup/CAPD, built with
cmake/g++ in WSL Fedora), path to capd-config given by $CAPD_CONFIG (default ~/opt/capd-src/build/bin/capd-config).
Runtime: about 1 minute (mostly compiling po_proof.cpp).
"""
import json
import math
import os
from fractions import Fraction
import re
import subprocess
import sys

import sympy as sp
from flint import fmpq, arb, acb, acb_mat

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
data = json.load(open(os.path.join(ROOT, 'certificates', 'orbits.json')))
FAILS = []


def claim(name, ok, detail=''):
    print('%-4s %s  %s' % ('PASS' if ok else 'FAIL', name, detail), flush=True)
    if not ok:
        FAILS.append(name)


k = [sp.Integer(v) for v in data['rate_constants']]
REACTIONS = [((0, 1), (2,)), ((2,), (0, 1)), ((2,), (3,)), ((3,), (4, 1)), ((4, 5), (6,)),
             ((6,), (4, 5)), ((6,), (7, 5)), ((7, 5), (8,)), ((8,), (7, 5)), ((8,), (0, 5))]
Nm = sp.zeros(9, 10)
for j, (re_, pr) in enumerate(REACTIONS):
    for i in re_:
        Nm[i, j] -= 1
    for i in pr:
        Nm[i, j] += 1
Wm = sp.Matrix([[0, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 1], [1, 0, 1, 1, 1, 0, 1, 1, 1]])
ZIDX = [2, 3, 4, 6, 7, 8]
E = sp.zeros(9, 6)
for j, i in enumerate(ZIDX):
    E[i, j] = 1
E[0, :] = sp.Matrix([[-1] * 6]); E[1, :] = sp.Matrix([[-1, -1, 0, 0, 0, 0]]); E[5, :] = sp.Matrix([[0, 0, 0, -1, 0, -1]])
Psel = sp.zeros(6, 9)
for j, i in enumerate(ZIDX):
    Psel[j, i] = 1


def g(x):
    r = []
    for j, (re_, pr) in enumerate(REACTIONS):
        v = k[j]
        for i in re_:
            v = v * x[i]
        r.append(v)
    return Nm * sp.Matrix(r)


def gsym(x, kk):
    r = []
    for j, (re_, pr) in enumerate(REACTIONS):
        v = kk[j]
        for i in re_:
            v = v * x[i]
        r.append(v)
    return Nm * sp.Matrix(r)


def phi(x1, x2, x6):
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = k
    return sp.Matrix([x1, x2, k1 / (k2 + k3) * x1 * x2, k1 * k3 / ((k2 + k3) * k4) * x1 * x2,
                      k1 * k3 * (k6 + k7) / ((k2 + k3) * k5 * k7) * x1 * x2 / x6, x6,
                      k1 * k3 / ((k2 + k3) * k7) * x1 * x2, k1 * k3 * (k9 + k10) / ((k2 + k3) * k8 * k10) * x1 * x2 / x6,
                      k1 * k3 / ((k2 + k3) * k10) * x1 * x2])


def hurwitz(b):
    def c(m):
        return b[m] if 0 <= m <= 6 else 0
    return [sp.Matrix(i, i, lambda r, s: c(2 * (r + 1) - (s + 1))).det() for i in range(1, 6)]


lam = sp.symbols('lam')
xsym = sp.symbols('x1:10')
for mu in sorted(set(o['mu'] for o in data['orbits'])):
    mu_q = sp.Rational(mu)
    xs = phi(sp.Rational(1, 3), mu_q, sp.Rational(1, 4))
    T = [sp.Rational(t) for t in next(o['totals'] for o in data['orbits'] if o['mu'] == mu)]
    claim('D1 mu=%s: x* = phi(1/3,%s,1/4) > 0, g(x*) = 0, W x* = totals' % (mu, mu),
          all(v > 0 for v in xs) and g(xs) == sp.zeros(9, 1) and list(Wm * xs) == T,
          'x* = (%s)' % ', '.join(str(v) for v in xs))
    J = Nm * sp.Matrix(10, 9, lambda j, a: sp.diff(sp.prod([xsym[i] for i in REACTIONS[j][0]]) * k[j], xsym[a]))
    A = (Psel * J.subs(dict(zip(xsym, xs))) * E)
    b = sp.Poly((lam * sp.eye(6) - A).det(), lam).all_coeffs()
    H = hurwitz(b)
    if mu_q < 59.1676975:
        ok = all(h > 0 for h in H) and b[6] > 0
        claim('D2 mu=%s: det H1..H5 > 0 and b6 > 0 (exact) -> x* locally asymptotically stable' % mu, ok,
              'det H5 = %.6e' % float(H[4]))
    else:
        ok = all(h > 0 for h in H[:4]) and H[4] < 0 and b[6] > 0
        claim('D2 mu=%s: det H1..H4 > 0, det H5 < 0, b6 > 0 (exact) -> x* not locally asymptotically stable' % mu,
              ok, 'det H5 = %.6e' % float(H[4]))
    # certified eigenvalue enclosures of the exact rational matrix (arb, isolated simple eigenvalues)
    Aq = acb_mat([[acb(arb(fmpq(int(sp.Rational(A[i, j]).p), int(sp.Rational(A[i, j]).q)))) for j in range(6)] for i in range(6)])
    evs = Aq.eig()
    npos = sum(1 for z in evs if z.real > 0); nneg = sum(1 for z in evs if z.real < 0)
    claim('D2b mu=%s: certified eigenvalue enclosures of the reduced Jacobian: %d with Re > 0, %d with Re < 0' % (mu, npos, nneg),
          npos + nneg == 6 and npos == (0 if mu_q < 59.1676975 else 2),
          ' '.join(z.str(6) for z in evs))

# D3: the CAPD vector field
src = open(os.path.join(HERE, 'po_proof.cpp'), encoding='utf-8').read()
body = re.search(r'FIELD =\s*((?:\s*"[^"]*"(?:\s*//[^\n]*)?)+);', src).group(1)
field = ''.join(re.findall(r'"([^"]*)"', body))
funs = field.split('fun:')[1].rstrip(';').split(',')
Ks, Ps, Ss = sp.symbols('K P S')
ksy = sp.symbols('k1:11')
a, b_, c, d, e, f = sp.symbols('a b c d e f')
loc = {'K': Ks, 'P': Ps, 'S': Ss, 'a': a, 'b': b_, 'c': c, 'd': d, 'e': e, 'f': f}
loc.update({'k%d' % (i + 1): ksy[i] for i in range(10)})
capd_f = [sp.sympify(s_, locals=loc) for s_ in funs]
zv = [a, b_, c, d, e, f]
xz = [Ss - sum(zv), Ks - a - b_, a, b_, c, Ps - d - f, d, e, f]
ref = Psel * gsym(sp.Matrix(xz), ksy)
claim('D3 po_proof.cpp vector field == P g(x(z)) (symbolic k, totals) with x1 = S - sum z, x2 = K - x3 - x4, x6 = P - x7 - x9',
      len(capd_f) == 6 and all(sp.expand(capd_f[i] - ref[i]) == 0 for i in range(6)))


def bracket(q):
    """exactly representable doubles lo <= q <= hi (q rational), as hex strings"""
    q = Fraction(int(sp.Rational(q).p), int(sp.Rational(q).q))
    x = float(q)
    lo = x if Fraction(x) <= q else math.nextafter(x, -math.inf)
    hi = x if Fraction(x) >= q else math.nextafter(x, math.inf)
    assert Fraction(lo) <= q <= Fraction(hi)
    return lo.hex(), hi.hex()


def write_params(path, kvals, totals):
    with open(path, 'w') as fh:
        for name, v in list(zip(['k%d' % (i + 1) for i in range(10)], kvals)) + list(zip(['K', 'P', 'S'], totals)):
            if isinstance(v, tuple):          # already an interval of rationals (lo, hi)
                lo, hi = bracket(v[0])[0], bracket(v[1])[1]
            else:
                lo, hi = bracket(v)
            fh.write('%s %s %s\n' % (name, lo, hi))


# D4: compile and run CAPD
capd_config = os.environ.get('CAPD_CONFIG', os.path.expanduser('~/opt/capd-src/build/bin/capd-config'))
exe = os.path.join('/tmp', 'po_proof_cms')
flags = subprocess.run([capd_config, '--cflags', '--libs'], capture_output=True, text=True).stdout.split()
r = subprocess.run(['g++', os.path.join(HERE, 'po_proof.cpp')] + flags + ['-o', exe], capture_output=True, text=True)
claim('D4.0 po_proof.cpp compiles against CAPD', r.returncode == 0, r.stderr[-300:])
if r.returncode != 0:
    print('FAILED:', FAILS)
    sys.exit(1)
log = []
for o in data['orbits']:
    T = [sp.Rational(t) for t in o['totals']]
    pfile = os.path.join('/tmp', 'cms_params_%s.txt' % o['label'])
    write_params(pfile, k, T)
    args = [exe, pfile] + [str(o['section_index']), o['section_value'], str(o['direction'])]         + o['z_guess'] + [o['radius'], str(o['order']), o['mode'], o['label']]
    r = subprocess.run(args, capture_output=True, text=True, timeout=600)
    log.append('$ ' + ' '.join(args[1:]) + '\n' + r.stdout + r.stderr)
    lines = r.stdout.splitlines()
    need = ['Krawczyk', 'positive'] + (['spectral radius'] if o['mode'] == 'stable' else [])
    ok = r.returncode == 0 and all(any(l.startswith('PASS') and n in l for l in lines) for n in need) \
        and not any(l.startswith('FAIL') for l in lines)
    rt = [l for l in lines if 'return time in' in l]
    claim('D4 %s: CAPD proof (%s)' % (o['label'], 'periodic orbit, orbitally asymptotically stable' if o['mode'] == 'stable'
                                      else 'periodic orbit'), ok, rt[0].split(': ', 1)[1] if rt else '')
open(os.path.join(ROOT, 'certificates', 'po_proof_log.txt'), 'w').write('\n'.join(log))
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('ALL CLAIMS PASS')
