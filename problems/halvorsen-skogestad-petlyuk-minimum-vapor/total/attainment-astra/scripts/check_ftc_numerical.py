"""Closed ternary FTC numerical checks using Agent P's read-only MESH solver.

Python 3, numpy/scipy/mpmath; normally < 1 minute, <=2 threads, <3 GiB.
Use --extended to also attempt 80 stages (known numerical failures, up to 4 min).
This is numerical evidence, not a certificate or a proof of attainment.
All output files are under this script's parent task folder.
Attainment: GPT-6-Astra via Codex CLI (effort xhigh).
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "2"
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import time
import warnings
import numpy as np
import mpmath as mp
from check_prefractionator import preferred

TASK = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TASK.parent / "simulation" / "scripts"))
from arrangements import regular_network, ftc_columns


def private_warning(message, category, filename, lineno, file=None, line=None):
    # Third-party numerical warnings must not put local paths in public logs.
    print(category.__name__+": "+str(message),file=file or sys.stderr,flush=True)


warnings.showwarning = private_warning


def main():
    a = np.array([4.,2.,1.])
    z = np.ones(3)/3
    q = 1.
    th, vp0, vb0, d0, b0 = preferred(list(map(mp.mpf,a)), [mp.mpf(1)/3]*3, mp.mpf(q))
    H = [float(sum(mp.mpf(a[i])*mp.mpf(1)/3/(th[j]-mp.mpf(a[i])) for i in range(j+1,3)))
         for j in range(2)]
    hv = max(H)
    rows = []
    print("NUMERICAL ONLY: alpha=(4,2,1), z=(1/3,1/3,1/3), q=1",flush=True)
    print("Underwood maximum =",format(hv,".16g"),flush=True)
    for extra in (0.06,0.03):
        R = hv+extra
        vp = float(vp0)+extra/10
        vb = vp-(1-q)
        t = float(d0[1])
        for n in ((20,40,80) if "--extended" in sys.argv else (20,40)):
            start = time.monotonic()
            net = regular_network(a,z,q,ftc_columns(3),n=n)
            side = R-vb
            fdict = {
                "V:TABC":vp, "L:TABC":vp-z[0]-t,
                "V:BABC":vb, "L:BABC":vb+z[1]-t+z[2],
                "V:TAB":R+1-q, "L:TAB":R+1-q-z[0],
                "V:BAB":side, "L:BAB":side+t,
                "V:TBC":side, "L:TBC":side-(z[1]-t),
                "V:BBC":R, "L:BBC":R+z[2],
                "P:A":z[0],"P:B":z[1],"P:C":z[2],
                "Q:PA":-(R+1-q),"Q:PC":R,
            }
            f = net.flows_from(fdict)
            flow_res = net.check_flows(f)
            assert flow_res < 1e-12
            x, info = net.solve(f)
            absres = float(np.max(np.abs(net.residual(x,f))))
            normres = float(np.max(np.abs(x.sum(axis=1)-1)))
            pp = np.array([x[net.nid("P"+p)] for p in "ABC"])
            impurity = float(np.max(1-np.diag(pp)))
            bounds = [sum(max(0.,-z[k]*sum(a[i]*pp[k,i]/(a[i]-float(theta)) for i in range(3)))
                          for k in range(3)) for theta in th]
            good = (bool(info["ok"]) and absres < 1e-10 and normres < 1e-10
                    and float(x.min()) >= 0 and R+1e-10 >= max(bounds))
            row = dict(n=n,extra=extra,R=R,max_impurity=impurity,
                       absolute_component_residual=absres,normalization_residual=normres,
                       relative_solver_residual=float(info["res"]),solver_ok=bool(info["ok"]),
                       stage_bound=max(bounds),minimum_composition=float(x.min()),
                       product_compositions=pp.tolist(),seconds=round(time.monotonic()-start,3),passed=bool(good))
            rows.append(row)
            print(("PASS " if good else "FAIL ")+json.dumps(row,sort_keys=True),flush=True)
    # Read-only comparison; capture observed values, including any unsuccessful re-solve.
    reference_path = TASK.parent/"simulation"/"runs"/"attain_T61_FTC.json"
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    assert abs(reference["V_FTC"]-hv) < 1e-12
    print("PASS Agent P reference maximum agrees to 1e-12",flush=True)
    print("REFERENCE "+json.dumps(reference["runs"],sort_keys=True),flush=True)
    out = dict(status="numerical evidence only",alpha=a.tolist(),z=z.tolist(),q=q,H=H,runs=rows,
               reference_runs=reference["runs"])
    (TASK/"runs").mkdir(exist_ok=True)
    filename = "ftc_exploratory.json" if "--extended" in sys.argv else "ftc_numerical.json"
    (TASK/"runs"/filename).write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    if not all(row["passed"] for row in rows):
        raise AssertionError("a finite-stage numerical check failed")
    print("PASS finite-stage numerical balances; no proof of the infinite-stage limit",flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("FAIL "+str(exc),flush=True)
        sys.exit(1)
