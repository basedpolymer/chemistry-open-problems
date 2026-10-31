#!/usr/bin/env python3
"""Figure for README: large-Peclet skewness Sk(tau) = 3 K / (2 int_0^tau <u A>)^(3/2) (floating point, illustration only)."""
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import jn_zeros
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from galerkin import Polygon, K_of_tau

taus = np.geomspace(1e-4, 30, 180)


def sk_modes(lam, ui, Tij):
    l = lam; a = ui
    K = K_of_tau(np.r_[0, l], np.r_[0, a], np.block([[np.zeros((1, 1)), np.zeros((1, len(l)))], [np.zeros((len(l), 1)), Tij]]), taus)
    m2 = np.array([2 * ((a ** 2 / l) * (t - (1 - np.exp(-l * t)) / l)).sum() for t in taus])
    return 3 * K / m2 ** 1.5


series = []
for n, name, col in [(3, 'triangle', '#eb6834'), (4, 'square', '#2a78d6'), (6, 'hexagon', '#1baf7a')]:
    P = Polygon(n, 30)
    series.append((name, col, sk_modes(P.lam[1:], P.ui[1:], P.Tij[1:, 1:])))
M = 700
j = jn_zeros(1, M)
lam = j ** 2
T = -2 * (lam[:, None] + lam[None, :]) / np.where(np.eye(M) > 0, 1, (lam[:, None] - lam[None, :]) ** 2)
T[np.arange(M), np.arange(M)] = 1 / 12
series.append(('disk', '#eda100', sk_modes(lam, -2 / lam, T)))

fig, ax = plt.subplots(figsize=(7.2, 3.9), dpi=100)
fig.patch.set_facecolor('#fcfcfb')
ax.set_facecolor('#fcfcfb')
ax.axhline(0, color='#52514e', lw=1)
for name, col, y in series:
    ax.plot(taus, y, color=col, lw=2, label=name, solid_capstyle='round')
ax.set_xscale('log')
ax.set_xlim(1e-4, 30)
ax.set_ylim(-0.05, 0.36)
ax.set_xlabel('time  (units of apothem² / diffusivity)', color='#52514e')
ax.set_ylabel('skewness, large Péclet number', color='#52514e')
for s in ['top', 'right']:
    ax.spines[s].set_visible(False)
for s in ['left', 'bottom']:
    ax.spines[s].set_color('#b8b6ae')
ax.tick_params(colors='#52514e')
ax.grid(axis='y', color='#e7e6e0', lw=.8)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc='upper right', labelcolor='#0b0b0b')
for name, col, y in series:
    ax.annotate(name, (taus[6], y[6]), xytext=(0, 5), textcoords='offset points', fontsize=9, color='#0b0b0b')
ax.annotate('triangle changes sign once', (1.2, 0), xytext=(0.012, -0.036), fontsize=9, color='#52514e',
            arrowprops=dict(arrowstyle='-', color='#52514e', lw=.8))
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, 'skewness-all-times.svg'))
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'skewness-all-times-preview.png'), dpi=110)
print([(n, float(y[0]), float(y.min()), float(y.max())) for n, c, y in series])
