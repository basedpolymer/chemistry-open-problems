"""Re-optimizes a candidate permissive point with strict margins (flows >= df, every inequality slack >= dm), at the
rational data used by certify_candidate.py (rational alpha, rational feed roots, exact z), starting from the candidate.
Writes logs/<candidate>_refined.json in the same format.  Usage: python refine_candidate.py cand.json df dm
Environment: WSL Fedora, Python 3, pyomo, ipopt."""
import sys, json, random
from fractions import Fraction as Fr
import pyomo.environ as pe
from indep_nlp import configurations, name, build, feas_check
D = json.load(open(sys.argv[1])); df, dm = float(sys.argv[2]), float(sys.argv[3])
cfg = [c for c in configurations(5) if name(c) == D['name']][0]
N = 5
alpha = [Fr(a).limit_denominator(10 ** 6) for a in D['alpha']]; alpha[-1] = Fr(1)
q = Fr(D['q']).limit_denominator(1000)
def roots_f(al, z, q):
    out = []
    for r in range(N - 1):
        a, b = float(al[r + 1]), float(al[r])
        for _ in range(200):
            c = 0.5 * (a + b)
            v = sum(float(x) * float(zz) / (float(x) - c) for x, zz in zip(al, z)) - (1 - float(q))
            if v < 0: a = c
            else: b = c
        out.append(0.5 * (a + b))
    return out
theta = [Fr(t).limit_denominator(10 ** 9) for t in roots_f(alpha, D['z'], q)]
# exact z
rows = [[alpha[p] / (alpha[p] - theta[r]) for p in range(N)] for r in range(N - 1)] + [[Fr(1)] * N]
rhs = [1 - q] * (N - 1) + [Fr(1)]
M_ = [list(r) + [b] for r, b in zip(rows, rhs)]
for col in range(N):
    piv = next(i for i in range(col, N) if M_[i][col] != 0); M_[col], M_[piv] = M_[piv], M_[col]
    pv = M_[col][col]; M_[col] = [v / pv for v in M_[col]]
    for i in range(N):
        if i != col and M_[i][col] != 0:
            f = M_[i][col]; M_[i] = [a - f * b for a, b in zip(M_[i], M_[col])]
z = [M_[i][N] for i in range(N)]
al_f = [float(a) for a in alpha]; z_f = [float(v) for v in z]; th_f = [float(t) for t in theta]
H = [sum(al_f[p] * z_f[p] / (th_f[r] - al_f[p]) for p in range(r + 1, N)) for r in range(N - 1)]
V = max(H)
M, F, U, dget, bget, reb = build(cfg, al_f, z_f, float(q), False, random.Random(0), th_f, V, df=df, dm=dm)
for v in M.component_data_objects(pe.Var):
    key = f'{v.parent_component().name}[{v.index()}]'
    if key in D['solution'] and not v.fixed:
        v.value = max(D['solution'][key], df) if v.parent_component().name in ('d', 'b') else D['solution'][key]
solver = pe.SolverFactory('ipopt')
solver.options['tol'] = 1e-12; solver.options['max_iter'] = 5000; solver.options['print_level'] = 0
solver.options['bound_relax_factor'] = 0.0
solver.solve(M, tee=False)
R = pe.value(M.obj)
print('refined: R =', R, ' V_FTC =', V, ' R/V =', R / V, ' feas viol =', feas_check(cfg, M, F, U, dget, bget, al_f, False))
sol = {f'{v.parent_component().name}[{v.index()}]': v.value for v in M.component_data_objects(pe.Var)}
out = dict(D); out.update({'alpha': al_f, 'z': z_f, 'q': float(q), 'H': H, 'R': R, 'V_FTC': V, 'solution': sol})
json.dump(out, open(sys.argv[1].replace('.json', f'_refined.json'), 'w'), indent=1)
