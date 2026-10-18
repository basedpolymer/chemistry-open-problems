"""Numerical illustration only, no certificate. Runtime about 10 s.
Python 3, scipy, numpy, matplotlib; one thread, below 1 GB.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from explore import shoot

start=time.monotonic()
data=json.loads(Path('certificates/exploration.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
for e in [.235,data['cusp']['epsilon'],.249]:
    aa=np.linspace(.05,16,220)
    ll=[4*shoot(a,e)[0] for a in aa]
    axs[0].plot(aa,ll,label=f'epsilon = {e:.6f}')
axs[0].set(xlabel='Central value alpha',ylabel='lambda',ylim=(2.7,3.3),xlim=(2,14))
axs[0].legend(fontsize=8)
for k in [0,1]:
    rows=[x for x in data['folds'] if len(x['folds'])==2]
    axs[1].plot([x['epsilon'] for x in rows],[x['folds'][k]['alpha'] for x in rows],'.-',color=['#147d92','#b65c24'][k])
axs[1].plot(data['cusp']['epsilon'],data['cusp']['alpha'],'ko',markersize=4,label='Certified local cusp')
axs[1].set(xlabel='epsilon',ylabel='Fold central value alpha',yscale='log')
axs[1].legend(fontsize=8)
for ax in axs:ax.grid(alpha=.2)
fig.suptitle('Numerical exploration: curves and sampled folds are not global certificates',fontsize=11)
Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/exploration.png',dpi=180)
print('PASS numerical illustration generated (not a proof)')
print('Runtime seconds:',round(time.monotonic()-start,2))
