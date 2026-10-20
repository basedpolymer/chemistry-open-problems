#!/usr/bin/env python3
"""spurious_zeros.py -- Agent C (red team).  The exact finite-stage system in junction / net-flow variables (L's
Definition 1.3, A's G_n) has zeros that are NOT steady states: their reconstructed stage profiles have negative mole
fractions and S(x_k) < 0 at some stages although every junction and product composition may look reasonable.  So
the hypothesis "S > 0 along every reconstructed profile" (L Theorem 1.4 (b), A Lemma A.1 (b)) is essential, and a
numerical zero of G_n must always be re-checked stage by stage.

Example 1: N = 3, alpha = (4,2,1), z = 1/3, q = 1, delta = 0.05, n = 24.  The forced solution (eta = 0.003) is
continued to eta = 0 in the chart variables; the end point is a zero of G_n (|G_n| < 1e-50) with negative stage
compositions.  The physical unforced steady state at the same flows (Agent N's stagenet) is a different zero.
Example 2: N = 5 case 'n5_mixed' (cases.py): the physical solution at n = 24 (from v_inf) used as the Newton start
at n = 48 converges to a non-physical zero at distance ~0.17 from v_inf; Newton started from v_inf at n = 48 finds
the physical one.
PASS means: the spurious zeros are found AND identified as non-physical by the stage check (i.e. the check works).
Usage: python spurious_zeros.py
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env).  Runtime: about 2 minutes.
"""
import sys
import time
import mpmath as mp
from petnet import Problem, Design, limit_config, System, reconstruct, stage_check, impurities, svals
import cases

t0 = time.time()
fails = []
print("Example 1 (N = 3, unforced limit of the forced family, n = 24)")
mp.mp.dps = 60
pb = Problem([4, 2, 1], [1, 1, 1], 1)


def design(eta):
    return Design(pb, 0.05, eta_nb_default=0.025, eta_bin={(0, 1): ('B', eta), (1, 2): ('T', eta)})


lim = limit_config(design(mp.mpf('0.003')))
n = 24
s1 = System(design(mp.mpf('0.003')), lim, n)
vv, nF, _ = s1.newton(s1.start())
for eta in ['0.001', '3e-4', '1e-4', '3e-5', '1e-5', '3e-6', '1e-6', '3e-7', '1e-7', '0']:
    s2 = System(design(mp.mpf(eta)), lim, n)
    vv, nF, _ = s2.newton(vv)
chk = stage_check(s2, vv, reconstruct(s2, vv))
X, W = s2.unpack(vv)
print(f"  zero of G_n: |G_n| = {mp.nstr(nF, 3)}; junction compositions J[1,2] {[mp.nstr(x, 5) for x in X[(0, 1)]]}, "
      f"J[2,3] {[mp.nstr(x, 5) for x in X[(1, 2)]]}; product impurities {[mp.nstr(x, 3) for x in impurities(s2, vv)]}")
print(f"  stage check: min mole fraction {mp.nstr(chk['minx'], 4)}, min S {mp.nstr(chk['minS'], 4)}, "
      f"{len(chk['negs'])} negative entries")
ok = nF < mp.mpf(10) ** -45 and (chk['negs'] or chk['minS'] <= 0)
print(f"[{'PASS' if ok else 'FAIL'}] non-physical zero of G_n found and rejected by the stage check")
if not ok:
    fails.append('example 1')

print("Example 2 (N = 5, case n5_mixed, continuation n = 24 -> 48)")
case = cases.CASES['n5_mixed']


def build():
    p = Problem(case['alpha'], case['z'], case['q'])
    d = Design(p, case['delta'], eta_nb_default=case['eta_nb_default'], eta_bin=case['eta_bin'])
    return p, d, limit_config(d)


mp.mp.dps = 80
p5, d5, l5 = build()
sA = System(d5, l5, 24)
vA, nFA, _ = sA.newton(sA.start())
chkA = stage_check(sA, vA, reconstruct(sA, vA))
print(f"  n = 24 from v_inf: |G_n| {mp.nstr(nFA, 3)}, physical: {not chkA['negs'] and chkA['minS'] > 0}")
mp.mp.dps = 130
p5, d5, l5 = build()
sB = System(d5, l5, 48)
vB, nFB, _ = sB.newton(mp.matrix([mp.mpf(str(x)) for x in vA]))
chkB = stage_check(sB, vB, reconstruct(sB, vB))
vl = sB.start()
dB = max(abs(vB[j] - vl[j]) for j in range(len(vB)) if vl[j] != 0)
print(f"  n = 48 from the n = 24 solution: |G_n| {mp.nstr(nFB, 3)}, |v - v_inf| {mp.nstr(dB, 3)}, min x "
      f"{mp.nstr(chkB['minx'], 4)}, min S {mp.nstr(chkB['minS'], 4)}, sigma_min {mp.nstr(svals(sB.jacobian(vB))[-1], 3)}")
vC, nFC, _ = sB.newton(sB.start())
chkC = stage_check(sB, vC, reconstruct(sB, vC))
dC = max(abs(vC[j] - vl[j]) for j in range(len(vC)) if vl[j] != 0)
print(f"  n = 48 from v_inf: |G_n| {mp.nstr(nFC, 3)}, |v - v_inf| {mp.nstr(dC, 3)}, min x {mp.nstr(chkC['minx'], 4)}, "
      f"min S {mp.nstr(chkC['minS'], 4)}, sigma_min {mp.nstr(svals(sB.jacobian(vC))[-1], 3)}")
spurious = nFB < mp.mpf(10) ** -100 and (chkB['negs'] or chkB['minS'] <= 0)
physical = nFC < mp.mpf(10) ** -100 and not chkC['negs'] and chkC['minS'] > 0
ok = spurious and physical
print(f"[{'PASS' if ok else 'FAIL'}] a non-physical zero (rejected) and the physical zero near v_inf coexist at n = 48")
if not ok:
    fails.append('example 2')
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
