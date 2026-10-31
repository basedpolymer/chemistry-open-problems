"""Exploration only: seek a four-source bistable module with a unary source.
An ACR controller attached to that source would use just one additional source.
"""
import itertools as it,json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

rng=np.random.default_rng(51603)
secs=float(sys.argv[1]) if len(sys.argv)>1 else 100
xsources=[tuple(1 if i==j else 0 for i in range(3)) for j in range(3)]
xsources += [tuple((1 if j==i else 0)+(1 if j==k else 0) for j in range(3)) for i in range(3) for k in range(i,3)]
B=np.array([[1.,0.],[0.,1.],[-1.,-1.]])
models=[]
for src in it.combinations(xsources,4):
 if not any(sum(y)==1 for y in src):continue
 canon=min(tuple(sorted(tuple(y[i] for i in p) for y in src)) for p in it.permutations(range(3)))
 if tuple(sorted(src))!=canon:continue
 rx=[(a,z) for a,y in enumerate(src) for z in xsources if sum(y)==sum(z) and z!=y]
 G=np.array([np.array(z)-np.array(src[a]) for a,z in rx]).T
 res=linprog(np.zeros(len(rx)),A_eq=G,b_eq=np.zeros(3),bounds=(1,None),method='highs')
 if not res.success:continue
 Q=[]
 for y in src:
  factors=[B[i] for i in range(3) for _ in range(y[i])]
  Q.append([0.,0.,0.] if len(factors)==1 else [factors[0][0]*factors[1][0],factors[0][0]*factors[1][1]+factors[0][1]*factors[1][0],factors[0][1]*factors[1][1]])
 flux=[res.x]
 for _ in range(12):
  r=linprog(rng.uniform(-1,1,len(rx)),A_eq=G,b_eq=np.zeros(3),bounds=(1e-4,1),method='highs')
  if r.success:flux.append(r.x)
 models.append((src,rx,G,np.array(Q),np.array(flux)))
print('models',len(models),flush=True)
start=time.monotonic();tries=0;multi=0
while time.monotonic()-start<secs:
 src,rx,G,Qsrc,flux=models[rng.integers(len(models))]
 lam=np.exp(rng.uniform(-5,5,len(flux)))@flux
 xs=np.exp(rng.uniform(-4,4,3))
 mons=np.prod(xs[None,:]**np.array(src),axis=1)
 N=np.zeros((3,4))
 for j,(a,z) in enumerate(rx):N[:,a]+=G[:,j]*lam[j]/mons[a]
 J=N@np.diag(mons)@np.array(src)@np.diag(1/xs)
 L=(J@B)[:2];Q=(N@Qsrc)[:2]
 poly=np.convolve(Q[0],L[1])-np.convolve(Q[1],L[0])
 candidates=[xs]
 if np.max(np.abs(poly))<1e-12:continue
 for t in np.polynomial.polynomial.polyroots(poly):
  if abs(t.imag)>1e-8:continue
  t=t.real;lin=L@np.array([1.,t]);qua=Q@np.array([1.,t,t*t])
  i=np.argmax(np.abs(qua))
  if abs(qua[i])<1e-14:continue
  u=-lin[i]/qua[i];xx=xs+B@np.array([u,t*u])
  if np.min(xx)<=1e-8*np.max(xx):continue
  if any(np.linalg.norm(xx-z)<1e-6*np.linalg.norm(xx) for z in candidates):continue
  mm=np.prod(xx[None,:]**np.array(src),axis=1)
  if np.max(np.abs(N@mm))>1e-6*np.max(np.abs(N)@mm):continue
  candidates.append(xx)
 stable=[]
 for xx in candidates:
  mm=np.prod(xx[None,:]**np.array(src),axis=1)
  jj=(N@np.diag(mm)@np.array(src)@np.diag(1/xx)@B)[:2]
  ev=np.linalg.eigvals(jj)
  if max(ev.real)<-1e-7*max(abs(ev)):stable.append((xx,ev))
 tries+=1
 if len(candidates)>=3:multi+=1
 if len(stable)>=2:
  data={'sources':src,'reactions':[(src[a],z,float(lam[j]/mons[a])) for j,(a,z) in enumerate(rx)],'steady_states':[z.tolist() for z in candidates],'stable':[z.tolist() for z,ev in stable]}
  Path(__file__).with_name('five-candidate.json').write_text(json.dumps(data,indent=2)+'\n')
  print('FOUND',data,flush=True);break
else:print('no hit',tries,'three-equilibrium cases',multi,flush=True)
