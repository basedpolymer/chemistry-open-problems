#!/usr/bin/env python3
"""Figure (floating point, illustration only): amplitude of the periodic orbits along mu = x2 for
k = (1,1,1,1,300,1,3,1,1,100), x1 = 1/3, x6 = 1/4, from the shooting continuation in po_branch.npy.
Environment: Python 3.12, numpy, scipy, matplotlib. Runtime: about 30 s.
Usage: python figure_bifurcation.py po_branch.npy ../../figures/bifurcation"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

k = np.array([1, 1, 1, 1, 300, 1, 3, 1, 1, 100], float)
x1, x6 = 1 / 3, 1 / 4
MU0 = 59.16769751409043          # certified Hopf point (certify_subcritical_hopf.py)
xs1 = np.array([0, 1, 1 / 6, 1 / 6, 2 / 675, 0, 1 / 18, 101 / 150, 1 / 600])
xs0 = np.array([x1, 0, 0, 0, 0, x6, 0, 0, 0])
ZIDX = [2, 3, 4, 6, 7, 8]
E = np.zeros((9, 6))
for j, i in enumerate(ZIDX):
    E[i, j] = 1
E[0, :] = -1; E[1, :] = [-1, -1, 0, 0, 0, 0]; E[5, :] = [0, 0, 0, -1, 0, -1]


def g(t, x):
    r = [k[0] * x[0] * x[1], k[1] * x[2], k[2] * x[2], k[3] * x[3], k[4] * x[4] * x[5], k[5] * x[6], k[6] * x[6],
         k[7] * x[7] * x[5], k[8] * x[8], k[9] * x[8]]
    return [-r[0] + r[1] + r[9], -r[0] + r[1] + r[3], r[0] - r[1] - r[2], r[2] - r[3], r[3] - r[4] + r[5],
            -r[4] + r[5] + r[6] - r[7] + r[8] + r[9], r[4] - r[5] - r[6], r[6] - r[7] + r[8], r[7] - r[8] - r[9]]

br = np.load(sys.argv[1])
rows = []
for r in br:
    z0, T, mu = r[:6], r[6], r[7]
    x0 = xs0 + mu * xs1 + E @ z0
    s = solve_ivp(g, [0, T], x0, method='DOP853', rtol=1e-10, atol=1e-12, dense_output=True)
    X = s.sol(np.linspace(0, T, 2000))
    rows.append((mu, X[4].max() - X[4].min(), r[8] > 1 + 1e-9))
rows = np.array(rows)
unst = rows[rows[:, 2] == 1]; st = rows[rows[:, 2] == 0]
iF = np.argmin(rows[:, 0])
blue, orange, ink, muted = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': muted, 'axes.labelcolor': ink, 'xtick.color': muted,
                     'ytick.color': muted})
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.plot([58.6, MU0], [0, 0], color=blue, lw=2, label='steady state / periodic orbit, stable')
ax.plot([MU0, 61.1], [0, 0], color=orange, lw=2, ls='--', label='steady state / periodic orbit, unstable')
u = np.vstack([[MU0, 0, 1], unst, rows[iF]])
ax.plot(u[:, 0], u[:, 1], color=orange, lw=2, ls='--')
sb = st[np.argsort(st[:, 1])]
sb = np.vstack([rows[iF], sb])
sb = sb[sb[:, 0] <= 61.1]
ax.plot(sb[:, 0], sb[:, 1], color=blue, lw=2)
ax.plot([MU0], [0], 'o', ms=6, color=ink, mec='white', mew=1.5, zorder=5)
ax.annotate('subcritical Hopf\n$\\mu_0\\approx 59.168$', (MU0, 0), (MU0 + 0.3, 0.12), color=ink, fontsize=8,
            arrowprops=dict(arrowstyle='-', color=muted, lw=0.8))
ax.plot([rows[iF, 0]], [rows[iF, 1]], 'o', ms=6, color=ink, mec='white', mew=1.5, zorder=5)
ax.annotate('fold of cycles\n$\\mu_F\\approx %.3f$' % rows[iF, 0], (rows[iF, 0], rows[iF, 1]),
            (rows[iF, 0] + 0.45, rows[iF, 1] + 0.02), color=ink, fontsize=8, ha='left',
            arrowprops=dict(arrowstyle='-', color=muted, lw=0.8))
ax.axvline(59, color=muted, lw=0.8, ls=':'); ax.text(59.02, 0.8, '$\\mu=59$', color=muted, fontsize=8)
ax.axvline(60, color=muted, lw=0.8, ls=':'); ax.text(60.02, 0.8, '$\\mu=60$', color=muted, fontsize=8)
ax.set_xlim(58.6, 61.1); ax.set_ylim(-0.03, 1.25)
ax.set_xlabel('$\\mu = x_2^*$ (free kinase at the steady state)')
ax.set_ylabel('amplitude of $x_5$ (S$_2$)')
for sp_ in ('top', 'right'):
    ax.spines[sp_].set_visible(False)
ax.grid(True, color='#e6e5e1', lw=0.6); ax.set_axisbelow(True)
ax.legend(loc='upper left', bbox_to_anchor=(0, 1.02), frameon=False, fontsize=7.5)
fig.tight_layout()
fig.savefig(sys.argv[2] + '.pdf'); fig.savefig(sys.argv[2] + '.svg')
print('fold of cycles at mu_F ~ %.6f, amplitude %.4f' % (rows[iF, 0], rows[iF, 1]))
for m in (59, 60):
    for lab, sub in (('unstable', unst), ('stable', st)):
        if len(sub) and sub[:, 0].min() <= m <= sub[:, 0].max():
            print('mu=%d %s amplitude ~ %.4f' % (m, lab, np.interp(m, *sub[np.argsort(sub[:, 0])][:, :2].T)))
