#!/usr/bin/env python3
"""Re-prove the complete branches 1--3 certificate chain (no sampling).

Python 3.12 or later, python-flint 0.9, SymPy, NumPy, mpmath. About 20--40 minutes,
three independent processes. Run from any directory; output defaults to
all-M/verification/recheck-2026-10-25. Each subprocess must exit successfully.
"""
import concurrent.futures
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", nargs="?", type=Path, default=BASE / "verification/recheck-2026-10-25")
parser.add_argument("--only", nargs="+", help="explicit job subset; summary is labelled partial")
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)
jobs = []

def add(name, script, *args):
    jobs.append((name, [sys.executable, str(script), *map(str, args)]))

for j, start in ((1, 15), (2, 44), (3, 87)):
    cert = BASE / "certificates" / ("l1h2_j1.csv" if j == 1 else f"l1_j{j}.csv")
    add(f"verify_j{j}", HERE / "verify_sweep.py", cert, 0, "0.078811099999999995", *(["--h2"] if j == 1 else []))
    if j > 1:
        par = (20, 20.4, -0.3, 0.25) if j == 2 else (26, 26.6, -0.15, 0.25)
        add(f"spectrum_j{j}", HERE / "h2j_sweep.py", j, cert, OUT / f"h2j_j{j}.csv", *par)
    add(f"dense_j{j}", ROOT / "verification/v02_dense_certify.py", OUT / f"dense_j{j}.csv",
        "--eig", "--range", j, start, 160)
    add(f"match_j{j}", HERE / "check_vs_result8.py", j, cert,
        ROOT / f"certificates/sweep_j{j}.csv.gz", ROOT / f"certificates/sweep_j{j}_b.csv.gz")
add("formulas", HERE / "check_formulas.py")
add("asymptotics", HERE / "asymptotics.py")
add("branch_links", HERE / "check_branch_links.py")
add("winding_regression", HERE / "check_winding.py")
add("small_stability", HERE / "check_small_stability.py")
add("limit_symbolic", HERE / "limit_symbolic.py")
add("onset_formulas", ROOT / "higher-branches/scripts/check_onset_formula.py")
for j in (1, 2, 3):
    add(f"onset_j{j}", HERE / "check_onsets.py", j, OUT)
all_names = [name for name, _ in jobs]
if args.only:
    if set(args.only) - set(all_names):
        parser.error("unknown job name")
    jobs = [job for job in jobs if job[0] in args.only]

def run(job):
    name, cmd = job
    t = time.monotonic()
    with (OUT / f"{name}.log").open("w") as fh:
        fh.write("COMMAND: " + " ".join(cmd) + "\n"); fh.flush()
        p = subprocess.run(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    result = dict(name=name, exit_code=p.returncode, seconds=round(time.monotonic() - t, 2))
    print(("PASS" if p.returncode == 0 else "FAIL"), name, result["seconds"], "s", flush=True)
    return result

if __name__ == "__main__":
    import flint
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run, job) for job in jobs]):
            results.append(future.result())
    ok = len(results) == len(jobs) and all(r["exit_code"] == 0 for r in results)
    report = dict(python=sys.version, python_flint=flint.__version__, jobs=results, passed=ok,
                  complete_suite=not args.only, required_jobs=all_names)
    (OUT / ("summary-partial.json" if args.only else "summary.json")).write_text(json.dumps(report, indent=2) + "\n")
    print("SUMMARY", "PASS" if ok else "FAIL", len(results), "jobs", "partial" if args.only else "complete", flush=True)
    sys.exit(0 if ok else 1)
