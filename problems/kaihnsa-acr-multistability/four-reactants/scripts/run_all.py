#!/usr/bin/env python3
"""Run the four exact certificates and retain their logs; Python 3.12+.
About one minute. This is a reproduction runner, not an independent verifier.
"""
import concurrent.futures as cf
import datetime,hashlib,json,platform,subprocess,sys,time
from pathlib import Path
import numpy,sympy

here=Path(__file__).resolve().parent
base=here.parent
out=base/'certificates';out.mkdir(exist_ok=True)
names=['certify_sources.py','certify_subspaces.py','certify_normal_forms.py','certify_five_example.py']
def run(name):
 start=time.monotonic()
 result=subprocess.run([sys.executable,str(here/name)],capture_output=True,text=True)
 log=result.stdout+result.stderr
 (out/(Path(name).stem+'.log')).write_text(log)
 passed=result.returncode==0 and 'FAIL' not in log and 'PASS' in log
 return {'script':name,'passed':passed,'exit_code':result.returncode,
         'pass_lines':sum(line.startswith('PASS ') for line in log.splitlines()),
         'seconds':round(time.monotonic()-start,3)}
with cf.ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(run,names))
for result in results:
 print(('PASS ' if result['passed'] else 'FAIL ')+result['script']+f": {result['pass_lines']} exact checks; {result['seconds']} s",flush=True)
files=[here/name for name in names]+[base/'paper.tex']
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'environment':{'python':platform.python_version(),'platform':platform.platform(),'numpy':numpy.__version__,'sympy':sympy.__version__},
        'results':results,'total_pass_lines':sum(r['pass_lines'] for r in results),
        'sha256':{str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        'independent_verification':False}
(out/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
sys.exit(0 if all(r['passed'] for r in results) else 1)
