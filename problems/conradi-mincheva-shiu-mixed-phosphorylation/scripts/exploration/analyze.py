# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
import glob
import numpy as np
from model import *

a = np.vstack([np.load(f) for f in sorted(glob.glob('search2_[1-9].npy'))])
a = a[np.isfinite(a[:, 13])]
print('total Hopf points', len(a), ' positive', (a[:, 13] > 0).sum())
pos = a[a[:, 13] > 0]; neg = a[a[:, 13] < 0]
names = ['k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k7', 'k8', 'k9', 'k10', 'x1', 'x2', 'x6']
print('median log10 of each parameter (pos vs neg):')
for i, n in enumerate(names):
    print('  %4s  pos %+.2f   neg %+.2f' % (n, np.median(np.log10(pos[:, i])), np.median(np.log10(neg[:, i]))))
eqd = np.isclose(a[:, 1], a[:, 5]) & np.isclose(a[:, 1], a[:, 8])
print('equal dissociation: total', eqd.sum(), 'positive', (eqd & (a[:, 13] > 0)).sum())
# stiffness and features for positives
rows = []
for r in pos:
    k = r[:10]; x = steady_state(k, *r[10:13])
    ev = np.linalg.eigvals(reduced_A(k, x))
    stiff = np.abs(ev).max() / r[14]
    rows.append([stiff, r[15]])
rows = np.array(rows)
score = np.log10(rows[:, 1]) - 0.7 * np.log10(rows[:, 0])
idx = np.argsort(-score)
np.set_printoptions(precision=4, linewidth=250, suppress=False)
print('best (stiffness, L) positives:')
for i in idx[:25]:
    r = pos[i]
    print('stiff %.3g  L %.3g eqd %d  k=%s x1=%.4g x2=%.4g x6=%.4g' % (rows[i, 0], rows[i, 1],
          np.isclose(r[1], r[5]) and np.isclose(r[1], r[8]), np.array2string(r[:10], precision=4), r[10], r[11], r[12]))
np.save('pos_all.npy', pos); np.save('pos_stiffL.npy', rows)
