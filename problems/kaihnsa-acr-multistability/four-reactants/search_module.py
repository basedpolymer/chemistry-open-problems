"""Numerical exploration, not a certificate: rank-two four-source modules.
Usage: python search_module.py n seconds [seed]. A unary/zero source permits
adding a one-new-source ACR controller without changing the module dynamics.
"""
import itertools as it,json,sys,time
from pathlib import Path
from collections import Counter
import numpy as np

n=int(sys.argv[1]);seconds=float(sys.argv[2]);seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
rng=np.random.default_rng(seed)
def complexes(n):
 out=[(0,)*n]
 for i in range(n):
  z=[0]*n;z[i]=1;out.append(tuple(z))
 for i in range(n):
  for j in range(i,n):
   z=[0]*n;z[i]+=1;z[j]+=1;out.append(tuple(z))
 return out
comp=np.array(complexes(n),dtype=np.int64)
rx=[(a,b) for a in range(len(comp)) for b in range(len(comp)) if a!=b]
g=np.array([comp[b]-comp[a] for a,b in rx],dtype=np.int64)
unique=np.unique(g,axis=0)
ij=np.array(list(it.combinations(range(len(unique)),2)))
trip=unique[ij]
cols=list(it.combinations(range(n),2))
wedges=np.stack([trip[:,0,i]*trip[:,1,j]-trip[:,0,j]*trip[:,1,i] for i,j in cols],axis=1)
div=np.gcd.reduce(np.abs(wedges),axis=1);ok=div!=0
trip=trip[ok];wedges=wedges[ok]//div[ok,None]
first=np.argmax(wedges!=0,axis=1);wedges*=np.sign(wedges[np.arange(len(wedges)),first])[:,None]
_,where=np.unique(wedges,axis=0,return_index=True);bases=trip[where]
print('spans',n,len(bases),flush=True)

def fullcone(c):
 dots=c[:,0,None]*c[None,:,1]-c[:,1,None]*c[None,:,0]
 return not np.any(np.all(dots>=0,axis=1)|np.all(dots<=0,axis=1))

def circuits(c):
 out=[];m=len(c)
 def cross(a,b):return a[0]*b[1]-a[1]*b[0]
 for i,j in it.combinations(range(m),2):
  if cross(c[i],c[j])==0 and c[i]@c[j]<0:
   w=np.zeros(m);w[i]=sum(abs(c[j]));w[j]=sum(abs(c[i]));out.append(w/w.sum())
 for i,j,k in it.combinations(range(m),3):
  w=np.array([cross(c[j],c[k]),cross(c[k],c[i]),cross(c[i],c[j])],float)
  if np.all(w>0) or np.all(w<0):
   v=np.zeros(m);v[[i,j,k]]=w/w.sum();out.append(v)
 return np.array(out)

models=[];start=time.monotonic()
for basis in bases:
 piv=next((i,j) for i,j in cols if basis[0,i]*basis[1,j]-basis[0,j]*basis[1,i])
 M=basis[:,piv];det=int(M[0,0]*M[1,1]-M[0,1]*M[1,0])
 adj=np.array([[M[1,1],-M[0,1]],[-M[1,0],M[0,0]]])
 mask=np.all((g[:,piv]@adj)@basis==det*g,axis=1)
 inds=np.flatnonzero(mask);G=g[mask]
 if not fullcone(G[:,piv]):continue
 available=sorted(set(rx[j][0] for j in inds))
 if len(available)<4:continue
 B=(adj@basis/det).T
 for src in it.combinations(available,4):
  if not any(sum(comp[a])<=1 for a in src):continue
  if not np.all(comp[list(src)].sum(axis=0)>0):continue
  sub=np.array([j for j in inds if rx[j][0] in src]);GG=g[sub]
  if not fullcone(GG[:,piv]):continue
  rays=circuits(GG[:,piv])
  if not len(rays):continue
  Y=comp[list(src)]
  qa=[]
  for y in Y:
   factors=[B[i] for i in range(n) for _ in range(y[i])]
   qa.append([0.,0.,0.] if len(factors)<2 else [factors[0][0]*factors[1][0],factors[0][0]*factors[1][1]+factors[0][1]*factors[1][0],factors[0][1]*factors[1][1]])
  sourceids=np.array([src.index(rx[j][0]) for j in sub])
  models.append((Y,sub,GG,rays,B,piv,np.array(qa),sourceids))
print('models',len(models),'preparation',time.monotonic()-start,flush=True)
if not models:sys.exit()
start=time.monotonic();tries=0;hist=Counter()
while time.monotonic()-start<seconds:
 Y,sub,G,rays,B,piv,Qsrc,sourceids=models[rng.integers(len(models))]
 chosen=rng.choice(len(rays),size=min(len(rays),rng.integers(2,6)),replace=False)
 lam=np.exp(rng.uniform(-5,5,len(chosen)))@rays[chosen]
 xs=np.exp(rng.uniform(-4,4,n));mons=np.prod(xs[None,:]**Y,axis=1)
 N=np.zeros((n,4))
 for j,a in enumerate(sourceids):N[:,a]+=G[j]*lam[j]/mons[a]
 if np.linalg.matrix_rank(N)<2:continue
 J=N@np.diag(mons)@Y@np.diag(1/xs)
 L=(J@B)[list(piv)];Q=(N@Qsrc)[list(piv)]
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
  mm=np.prod(xx[None,:]**Y,axis=1)
  if np.max(np.abs(N@mm))>1e-6*np.max(np.abs(N)@mm):continue
  candidates.append(xx)
 stable=[]
 for xx in candidates:
  mm=np.prod(xx[None,:]**Y,axis=1)
  jj=(N@np.diag(mm)@Y@np.diag(1/xx)@B)[list(piv)]
  ev=np.linalg.eigvals(jj)
  if max(ev.real)<-1e-7*max(abs(ev)):stable.append((xx,ev))
 tries+=1;hist[(len(candidates),len(stable))]+=1
 if len(stable)>=2:
  reactions=[(comp[rx[sub[j]][0]].tolist(),comp[rx[sub[j]][1]].tolist(),float(lam[j]/mons[sourceids[j]])) for j in range(len(sub)) if lam[j]>1e-14]
  data={'n':n,'reactions':reactions,'steady_states':[z.tolist() for z in candidates],'stable':[z.tolist() for z,ev in stable]}
  Path(__file__).with_name('five-candidate.json').write_text(json.dumps(data,indent=2)+'\n')
  print('FOUND',data,flush=True);break
else: print('no hit',tries,'hist',dict(hist),flush=True)
