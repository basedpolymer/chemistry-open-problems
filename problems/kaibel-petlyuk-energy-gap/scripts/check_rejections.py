#!/usr/bin/env python3
"""Negative certificate checks; Python >= 3.10, standard library, < 1 second.

Ensure malformed root locations, forged flow bounds and a false bottleneck
cannot be accepted. No committed certificate is modified.
"""
from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
import json
import sys
from verify import ROOT, Interval, verify_case


def rejected(label, case):
    try:
        with redirect_stdout(StringIO()):
            verify_case(case)
    except ValueError:
        print("PASS reject " + label)
        return
    raise ValueError("accepted " + label)


def main():
    source = json.loads((ROOT / "certificates/examples.json").read_text())["cases"][0]
    bad = deepcopy(source)
    bad["roots"]["theta_A"] = ["3", "7/2"]
    rejected("root in the wrong volatility interval", bad)
    bad = deepcopy(source)
    bad["bounds"]["Kaibel"] = ["1", "100000000001/100000000000"]
    rejected("forged Kaibel minimum", bad)
    bad = deepcopy(source)
    bad["bottleneck"] = "K_C"
    rejected("wrong binding section", bad)
    bad = deepcopy(source)
    bad["bounds"]["gap"] = ["-1", "-99999999999/100000000000"]
    rejected("negative gap", bad)
    try:
        Interval(1) / Interval(-1, 1)
    except ValueError:
        print("PASS reject division across a pole")
    else:
        raise ValueError("accepted division across a pole")
    print("PASS all rejection checks")


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print("FAIL", exc)
        sys.exit(1)
