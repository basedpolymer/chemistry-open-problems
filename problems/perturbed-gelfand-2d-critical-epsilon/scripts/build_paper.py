"""Build using the repository script; keep all temporary writes in this folder.
Windows: py -3 scripts/build_paper.py. Typically 10-60 s; no new packages.
"""
import os
from pathlib import Path
import subprocess
import sys
import time

start=time.monotonic()
folder=Path(__file__).resolve().parents[1]
repo=folder.parents[1]
scratch=folder/'tmp'/'build'
scratch.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,TMP=str(scratch),TEMP=str(scratch),TMPDIR=str(scratch))
cmd=[sys.executable,'templates/build_tex.py',folder.relative_to(repo).as_posix()+'/paper.tex']
run=subprocess.run(cmd,cwd=repo,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace')
output=run.stdout+run.stderr
for path,replacement in [(folder,'.'),(repo,'../..')]:
    output=output.replace(str(path),replacement).replace(path.as_posix(),replacement)
(folder/'certificates'/'build.log').write_text(output,encoding='utf-8')
print(output,end='')
print(('PASS' if run.returncode==0 else 'FAIL')+' shared-template PDF build')
print('Runtime seconds:',round(time.monotonic()-start,2))
raise SystemExit(run.returncode)
