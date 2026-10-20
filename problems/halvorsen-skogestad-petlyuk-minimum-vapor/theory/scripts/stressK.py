"""Exploration (floating point, NOT a certificate). Stress test of the recovery-ordering and family-K lemmas (N=4):
minimize over the model (Agent C's nlp.build_model)
  (i)  r_B - r_C of the root column, (ii) t_1 + t_2 + kappa_1 at theta_2, (iii) R - H_2.
Uses build_model(..., permissive=True) (TG model, excess vapor allowed).
Environment: WSL Fedora, pyomo + ipopt.  Runtime: about 15 minutes for 6 trials."""
import sys, random, math
sys.path.insert(0, '../../computation/scripts')
import pyomo.environ as pyo
from hsmodel import all_configurations, config_name, feed_roots, petlyuk_peaks
from nlp import build_model, random_start, solver, max_violation
A = all_configurations(4)

def lam_T(m, s, th):
    cols, info, reb, side, al, Fd, q = m._info
    d = info[s]
    return pyo.value(m.VT[s]) - sum(al[p]*pyo.value(m.d[s,p])/(al[p]-th) for p in range(d['i'], d['k']+1))

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
worst = [1e9, 1e9, 1e9]
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 10):
    r = [math.exp(rng.uniform(math.log(1.05), math.log(5))) for _ in range(3)]
    alpha = [r[0]*r[1]*r[2], r[1]*r[2], r[2], 1.0]
    z = [rng.uniform(0.01, 1) for _ in range(4)]; s_ = sum(z); z = [x/s_ for x in z]
    q = rng.uniform(-0.3, 1.3)
    th2 = float(feed_roots(alpha, z, q)[1])
    H2 = float(petlyuk_peaks(alpha, z, q)[1]) - (1-q)
    for cid in range(88, 104):
        cfg = A[cid]
        m = build_model(cfg, alpha, z, q, permissive=True)
        cols, info, reb, side, al, Fd, qq = m._info
        root, c1, c2 = (1,4), (1,3), (2,4)
        hx1 = cfg['hx'].get(c1) == 'HX'
        kappa1 = (m.VT[root] - sum(m.d[root,p] for p in range(1,4))) if hx1 else 0
        t1 = m.VT[c1] - sum(al[p]*m.d[c1,p]/(al[p]-th2) for p in (1,2))
        t2 = m.VT[c2] - al[2]*m.d[c2,2]/(al[2]-th2)
        objs = [m.d[root,2]/Fd[2] - m.d[root,3]/Fd[3], t1 + t2 + kappa1, m.obj.expr - H2]
        m.obj.deactivate()
        for k, ex in enumerate(objs):
            m.o2 = pyo.Objective(expr=ex)
            for st in range(4):
                random_start(m, rng, 1.0)
                try:
                    res = solver().solve(m, load_solutions=False)
                except Exception:
                    continue
                if str(res.solver.termination_condition) in ('optimal', 'locallyOptimal', 'acceptable'):
                    m.solutions.load_from(res)
                    if max_violation(m) < 1e-7:
                        v = pyo.value(ex)
                        if v < worst[k]:
                            worst[k] = v
                            print('new worst', k, v, cid, alpha, z, q, flush=True)
            m.del_component(m.o2)
print('min r_B - r_C:', worst[0], ' min t1+t2+kappa1:', worst[1], ' min R - H2:', worst[2])
