"""Verifier's numerical probe of the 20 open (basic configuration, j) pairs for N = 5 in the PERMISSIVE model,
with the verifier's own solver (indep_nlp.py), at random parameter points where H_j is the highest peak and
the other peaks are within 10 % of it (the hardest regime for excess vapor).  Variants: all thermal couplings,
and all heat exchangers.  Numerical evidence only (local solver, multistart).
Usage: python open_n5_probe.py part nparts seed      (runs pairs part, part+nparts, ...)
Environment: WSL Fedora, Python 3, pyomo, ipopt."""
import sys, json, random
import pyomo.environ as pe
from indep_nlp import configurations, name, peaks, build, feas_check

part, nparts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ntrials = int(sys.argv[4]) if len(sys.argv) > 4 else 3          # optional: trials per pair
tc_only = len(sys.argv) > 5 and sys.argv[5] == 'tc'                # optional: all-TC variant only
qchoices = [1.0, 0.5, 1.2, 0.0] if len(sys.argv) <= 6 else [float(v) for v in sys.argv[6].split(',')]
pairs = json.load(open('open_pairs_N5.json'))
cfgs = configurations(5)
solver = pe.SolverFactory('ipopt')
solver.options['tol'] = 1e-10; solver.options['max_iter'] = 3000; solver.options['print_level'] = 0
rng = random.Random(seed + part)

def point_with_top_peak(j):
    while True:
        al = [1.0]
        for _ in range(4):
            al.append(al[-1] * (1 + rng.uniform(0.08, 1.5)))
        al = al[::-1]
        z = [rng.uniform(0.05, 1) for _ in range(5)]; s = sum(z); z = [x / s for x in z]
        q = rng.choice(qchoices)
        th, H = peaks(al, z, q)
        if max(range(4), key=lambda i: H[i]) == j - 1 and sorted(H)[-2] > 0.9 * H[j - 1]:
            return al, z, q, th, H

worst = None
for idx in range(part, len(pairs), nparts):
    P = pairs[idx]
    S = set(tuple(s) for s in P['streams'])
    variants = [c for c in cfgs if set(c['streams']) == S]
    sel = [v for v in variants if all(t == 'TC' for t in v['hx'].values())] + \
          ([] if tc_only else [v for v in variants if all(t == 'HX' for t in v['hx'].values())])
    for trial in range(ntrials):
        al, z, q, th, H = point_with_top_peak(P['j'])
        V = max(H)
        for cfg in sel:
            best = None
            for s in range(10):
                M, F, U, dget, bget, reb = build(cfg, al, z, q, False, rng, th, V)
                try:
                    solver.solve(M, tee=False, load_solutions=True)
                except Exception:
                    continue
                if feas_check(cfg, M, F, U, dget, bget, al, False) > 1e-7:
                    continue
                R = pe.value(M.obj)
                if best is None or R < best:
                    best = R
                    if R < (1 - 1e-6) * V:     # dump candidate points below the bound
                        sol = {str(v.index()): v.value for v in M.component_data_objects(pe.Var)}
                        sol = {f'{v.parent_component().name}[{v.index()}]': v.value
                               for v in M.component_data_objects(pe.Var)}
                        json.dump({'basic': P['basic'], 'j': P['j'], 'name': name(cfg), 'alpha': al, 'z': z, 'q': q,
                                   'H': H, 'R': R, 'V_FTC': V, 'solution': sol},
                                  open(f"logs/candidate_{P['basic']}_{P['j']}_s{seed}_t{trial}.json", 'w'), indent=1)
            ratio = best / V if best is not None else None
            print(f"basic {P['basic']} j={P['j']} trial {trial} {'TC' if cfg is sel[0] else 'HX'}: "
                  f"ratio {ratio if ratio is None else round(ratio, 10)}  H={[round(h, 4) for h in H]} q={q}", flush=True)
            if ratio is not None and (worst is None or ratio < worst[0]):
                worst = (ratio, P['basic'], P['j'], name(cfg))
print('WORST', worst)
