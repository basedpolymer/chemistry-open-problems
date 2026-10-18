#!/usr/bin/env python3
"""Run all offline certificates. Python >=3.10, python-flint, mpmath.
Runtime about 20 seconds; one worker/thread, <250 MB RAM.
Writes certificates/verification.log; this is not independent verification.
Prints PASS/FAIL, exits nonzero on any child failure.
"""
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    start=time.perf_counter()
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
        env[key]='1'
    env['PYTHONDONTWRITEBYTECODE']='1'
    with (ROOT/'certificates'/'verification.log').open('w') as log:
        def say(text):
            print(text,flush=True)
            log.write(text+'\n')
            log.flush()
        say('Local deterministic certificate run; not independent verification.')
        say('Python '+platform.python_version()+'; platform '+platform.platform())
        for name in ['verify_exact.py','verify_intervals.py','verify_independent.py']:
            say('RUN '+name)
            proc=subprocess.Popen([sys.executable,str(ROOT/'scripts'/name)],
                                  cwd=ROOT,env=env,stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT,text=True)
            for line in proc.stdout:
                say(line.rstrip('\n'))
            if proc.wait()!=0:
                say('FAIL '+name)
                return 1
            say('PASS '+name)
        say(f'PASS all certificates; elapsed={time.perf_counter()-start:.2f}s')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
