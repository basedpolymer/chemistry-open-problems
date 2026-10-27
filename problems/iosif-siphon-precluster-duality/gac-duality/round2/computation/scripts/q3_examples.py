"""Q3's explicit examples NM, NB, NT and the K8 pair M3/M4 (parts B and K8s of q3_tests2.py), runnable alone.

Environment: Python 3.12 (WSL Fedora 44), standard library + msolve 0.10.1; exact.  Usage (from
round2/computation/): python scripts/q3_examples.py   (see q3_tests2.py for the statements; same code).
"""
import sys
import time
import platform

import q3_tests2

if __name__ == "__main__":
    t0 = time.time()
    print(f"# q3_examples.py  Python {platform.python_version()} on {platform.system()}; exact; msolve 0.10.1")
    ok = q3_tests2.examples()
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)
