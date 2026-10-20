"""val_ternary.py -- validation 2: the six ternary arrangements of Halvorsen's Table 6.1 (alpha = 4,2,1, z = 1/3,
q = 1) in the stage-by-stage model: minimum total vaporization at purity 1-eps with n stages per section, compared
with the infinite-stage values (TG transfer conventions).  Also checks R >= stage-level bound at every optimum.
Run: python val_ternary.py  (~5-10 min, 1 thread)."""
import json, sys
from study import run_config
from arrangements import ftc_columns

alpha = (4.0, 2.0, 1.0); z = (1 / 3, 1 / 3, 1 / 3); q = 1.0
configs = {
    'direct': ({(0, 2): (0, 1), (1, 2): (1, 2)}, {'BC': 'HX'}, 2.0718),
    'indirect': ({(0, 2): (1, 2), (0, 1): (0, 1)}, {'AB': 'HX'}, 2.0324),
    'side_rectifier': ({(0, 2): (0, 1), (1, 2): (1, 2)}, {'BC': 'TC'}, 1.8819),
    'side_stripper': ({(0, 2): (1, 2), (0, 1): (0, 1)}, {'AB': 'TC'}, 1.8819),
    'prefractionator_HX': (ftc_columns(3), {'AB': 'HX', 'BC': 'HX'}, 1.5556),
    'petlyuk': (ftc_columns(3), {}, 1.3657),
}
out = []
for name, (cols, tr, ref) in configs.items():
    for eps in (1e-3, 1e-5):
        for n in (20, 40, 80):
            r = run_config(alpha, z, q, cols, tr, n, eps)
            row = dict(config=name, n=n, eps=eps, R=r['R'], ref_inf=ref, ratio=r['R'] / ref, bound=r['bound'],
                       R_minus_bound=r['R'] - r['bound'], status=r['status'], max_imp=r['max_imp'],
                       resolve_ok=bool(r['resolve_ok']), sec=round(r['sec'], 1))
            out.append(row); print(json.dumps(row)); sys.stdout.flush()
json.dump(out, open('../runs/val_ternary.json', 'w'), indent=1)
