"""Strict-model minimisation (verifier's solver) for the configuration of the permissive counterexample, at the
same rational data: Theorem 'strict' predicts R >= V_FTC.  Also re-solves the permissive model there.
Usage: python strict_at_counterexample.py permissive_counterexample_N5.json nstart.  Numerical evidence only."""
import sys, json, random
from fractions import Fraction as Fr
import pyomo.environ as pe
from indep_nlp import configurations, name, build, feas_check
C = json.load(open(sys.argv[1])); nst = int(sys.argv[2])
cfg = [c for c in configurations(5) if name(c) == C['configuration']][0]
al = [float(Fr(a)) for a in C['alpha']]; z = [float(Fr(v)) for v in C['z']]; q = float(Fr(C['q']))
th = [float(Fr(t)) for t in C['feed_roots']]; V = float(Fr(C['V_FTC']))
solver = pe.SolverFactory('ipopt'); solver.options['tol'] = 1e-10; solver.options['max_iter'] = 3000; solver.options['print_level'] = 0
rng = random.Random(3)
for strict in (True, False):
    best = None
    for s in range(nst):
        M, F, U, dget, bget, reb = build(cfg, al, z, q, strict, rng, th, V)
        try: solver.solve(M, tee=False)
        except Exception: continue
        if feas_check(cfg, M, F, U, dget, bget, al, strict) > 1e-8: continue
        R = pe.value(M.obj); best = R if best is None else min(best, R)
    print(('strict' if strict else 'permissive'), 'min R/V_FTC over', nst, 'starts:', None if best is None else best / V)
