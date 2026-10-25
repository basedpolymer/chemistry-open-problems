"""Positive-cone continuation with moving interior apices and a free box.

Python 3.12+, numpy, numba, scipy, gemmi; one process, approximately 25 minutes
for the eight default entries at level 0 and 5000 iterations per stage.
Restoring the original catalogue archive is required to regenerate starts.
Outputs to work/barrier/ are exploratory. Certify final saved witnesses with
check_barrier_witnesses.py --build; no stationarity or optimality is claimed.
The barrier preserves a sufficient star-shaped subclass, not all admissible cells.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import dfm
import numpy as np
import numba as nb
from dfm import registry, mesh, energy, relax
from dfm.lbfgs import lbfgs
from run_dfm import initial_mesh
import argparse, csv, json, time

DEFAULT_NAMES = ['cat:N55','cat:N59','cat:N60','cat2:Y02','cat2:Y03','cat2:Y04','cat2:Y05','cat2:Y06']

@nb.njit(cache=True)
def cone_barrier(X,C,L,tri,tsh,tf,fi,fj,fsj,floor=0.):
    gx=np.zeros_like(X); gc=np.zeros_like(C); value=0.; mind=1e300
    count=2*len(tri)
    ps=np.zeros((3,3))
    for t in range(len(tri)):
        f=tf[t]
        for k in range(3):
            ps[k]=X[tri[t,k]]+L@tsh[t,k].astype(np.float64)
        for side in range(2):
            cell=fi[f] if side==0 else fj[f]
            o=C[cell].copy()
            if side==1:
                o+=L@fsj[f].astype(np.float64)
            sign=1. if side==0 else -1.
            a,b,c=ps[0]-o,ps[1]-o,ps[2]-o
            da,db,dc=np.cross(b,c),np.cross(c,a),np.cross(a,b)
            d=sign*np.dot(a,da)
            mind=min(mind,d)
            if d<=floor:
                return np.inf,gx,gc,mind
            value-=np.log(d-floor)/count
            coeff=-sign/(count*(d-floor))
            gx[tri[t,0]]+=coeff*da
            gx[tri[t,1]]+=coeff*db
            gx[tri[t,2]]+=coeff*dc
            gc[cell]-=coeff*(da+db+dc)
    return value,gx,gc,mind

def sym(z):
    return np.array([[z[0],z[3],z[4]],[z[3],z[1],z[5]],[z[4],z[5],z[2]]])

def stage(M,C,mu,maxiter):
    x0=M.X.copy(); c0=C.copy(); l0=M.L.copy(); ref0=M.ref.copy()
    nx=x0.size; nc=c0.size
    floor=6*abs(np.linalg.det(l0))/M.ncell*1e-7
    z0=np.zeros(nx+nc+6)
    def fg(z):
        xb=x0+z[:nx].reshape(x0.shape)
        cb=c0+z[nx:nx+nc].reshape(c0.shape)
        B,bx,bc,mind=cone_barrier(xb,cb,l0,M.tri,M.tsh,M.tface,M.fi,M.fj,M.fsj,floor)
        T=np.eye(3)+sym(z[nx+nc:])
        if not np.isfinite(B) or np.linalg.eigvalsh(T)[0]<1e-5:
            return np.inf,np.zeros_like(z)
        M.L=T@l0; M.ref=ref0@T.T
        F,gx,vir=energy.dfm_F_grad_virial(M,xb@T.T)
        vir=vir@np.linalg.inv(T).T
        ge=np.array([vir[0,0],vir[1,1],vir[2,2],vir[0,1]+vir[1,0],vir[0,2]+vir[2,0],vir[1,2]+vir[2,1]])
        return F+mu*B,np.r_[(gx@T+mu*bx).ravel(),(mu*bc).ravel(),ge]
    z,f,g,info=lbfgs(fg,z0,m=30,maxiter=maxiter,gtol=1e-10,stall_window=200,stall_tol=1e-13,noise=1e-17)
    T=np.eye(3)+sym(z[nx+nc:])
    M.X=(x0+z[:nx].reshape(x0.shape))@T.T
    C=(c0+z[nx:nx+nc].reshape(c0.shape))@T.T
    M.L=T@l0; M.ref=ref0@T.T
    B,_,_,mind=cone_barrier(M.X,C,M.L,M.tri,M.tsh,M.tface,M.fi,M.fj,M.fsj)
    return C,info,mind,np.abs(g).max()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('names',nargs='*',default=DEFAULT_NAMES)
    ap.add_argument('--maxiter',type=int,default=5000)
    ap.add_argument('--level',type=int,default=0)
    args=ap.parse_args()
    outdir=ROOT/'work'/'barrier'
    outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name in args.names:
        M=initial_mesh(registry.get(name)); C=M.ref.copy()
        for lev in range(args.level+1):
            if lev: M=mesh.refine(M)
            for mu in [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]:
                t=time.time()
                C,info,mind,gmax=stage(M,C,mu,args.maxiter)
                vals=energy.dfm_values(M)
                row=dict(name=name,level=lev,mu=mu,F=float(vals['F']),min_cone=float(mind/(6*vals['Vbox']/M.ncell)),
                         centroid_cone=relax.cone_check(M),gmax=float(gmax),nit=info['nit'],stop=info['msg'],seconds=time.time()-t)
                print(json.dumps(row),flush=True)
                rows.append(row)
                with (outdir/'results.jsonl').open('a') as f: f.write(json.dumps(row)+'\n')
                np.savez_compressed(outdir/(name.replace(':','_')+f'_L{lev}.npz'),X=M.X,C=C,L=M.L,tri=M.tri,tsh=M.tsh,tface=M.tface,fi=M.fi,fj=M.fj,fsj=M.fsj,ncell=M.ncell)
    if args.level == 0 and sorted(args.names) == sorted(DEFAULT_NAMES):
        # A full run records its own diagnostics; a subset never overwrites
        # the complete catalogue supplement with partial data.
        with (ROOT/'data'/'barrier_continuation.csv').open('w', newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)

if __name__=='__main__': main()
