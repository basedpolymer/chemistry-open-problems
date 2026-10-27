"""Run A2 certificates and record sanitized logs, entirely within this folder.

Environment: Python 3.14 / standard library. Certificates use exact arithmetic.
Runtime: about 1 second by default; about 3 minutes with --searches.
Every child script prints PASS/FAIL per claim. This runner exits nonzero if a
child fails. Runtime measurements are not mathematical certificates.
"""
from pathlib import Path
import sys
import os
import subprocess
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[4]
sys.dont_write_bytecode=True


def main():
    start=time.perf_counter()
    print('Environment: Python standard library; exact-certificate subprocesses.')
    print('Runtime: about 1 second; about 3 minutes with --searches.')
    jobs=[('certify_out_all.py',[]),('certify_secondary.py',[]),('cycle_candidate.py',[])]
    if '--searches' in sys.argv:
        jobs += [('explore_homogeneous.py',['6']),('search_rates.py',['old']),
                 ('search_linear_vertices.py',['10000']),('search_center_vertices.py',['10000'])]
    logs=HERE.parent/'logs'
    logs.mkdir(exist_ok=True)
    failures=[]
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8')
    for script,args in jobs:
        rel=(HERE/script).relative_to(REPO)
        done=subprocess.run([sys.executable,str(rel),*args],cwd=REPO,env=env,capture_output=True,text=True,encoding='utf-8')
        content=done.stdout+done.stderr
        content=content.replace(str(REPO),'.').replace(REPO.as_posix(),'.')
        target=logs/(Path(script).stem+'.log')
        target.write_text(content,encoding='utf-8')
        ok=done.returncode==0 and 'FAIL:' not in content
        print(('PASS' if ok else 'FAIL')+': RUN '+script+' -> logs/'+target.name,flush=True)
        if not ok:
            failures.append(script)
    print('Runtime seconds:',round(time.perf_counter()-start,6))
    return 1 if failures else 0


if __name__=='__main__':
    raise SystemExit(main())
