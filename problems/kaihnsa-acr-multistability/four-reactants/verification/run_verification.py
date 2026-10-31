#!/usr/bin/env python3
"""Runs the verifier's reimplementation in order. About 12 minutes on 3 cores.

Run under WSL with the scientific environment (Python 3.12, sympy, numpy, scipy, msolve). Each step prints
PASS/FAIL lines; this runner prints one PASS or FAIL per step, writes the logs to verification/logs/ and exits
nonzero on failure. Step 7 re-runs solver scripts in a temporary directory; nothing in certificates/ is changed.
The floating-point steps (stress_count, stress_networks) are tests, not part of any proof. Sampling is seeded.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, 'logs')
STEPS = [
    ('enum_sources.py', []),
    ('enum_spans.py', []),
    ('enum_spans.py', ['n4']),
    ('compare_solver.py', []),
    ('check_forms.py', []),
    ('check_example.py', []),
    ('mutate_solver.py', []),
    ('flux_redundancy.py', []),
    ('kinetic_structure.py', []),
    ('stress_count.py', ['100']),
    ('stress_networks.py', ['12']),
]


def main():
    os.makedirs(LOGS, exist_ok=True)
    ok = True
    for script, args in STEPS:
        t0 = time.monotonic()
        r = subprocess.run([sys.executable, os.path.join(HERE, script)] + args, capture_output=True, text=True, cwd=HERE)
        log = r.stdout + r.stderr
        name = script.removesuffix('.py') + ('-' + args[0] if args else '')
        with open(os.path.join(LOGS, name + '.log'), 'w') as f:
            f.write(log)
        good = r.returncode == 0 and any(ln.startswith('RESULT PASS') for ln in log.splitlines())
        print('PASS' if good else 'FAIL', name, f'{time.monotonic() - t0:.0f} s', flush=True)
        ok = ok and good
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
