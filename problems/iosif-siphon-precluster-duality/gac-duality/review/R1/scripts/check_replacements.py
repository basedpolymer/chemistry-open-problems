#!/usr/bin/env python3
"""R1 helper: check that every `old` text quoted in review.md (as `- old ...: `...``) is a literal substring of
paper.tex, occurring exactly once, and that every R-22 anchor `after `...`` occurs too.
Standard library only.  Usage (from this folder): py -3 scripts/check_replacements.py review.md ../../../paper.tex"""
import re
import sys

review = open(sys.argv[1], encoding="utf-8").read()
paper = open(sys.argv[2], encoding="utf-8").read()
olds = re.findall(r"- old[^:`]*: `(.*?)`\n", review, flags=re.S)
anchors = re.findall(r"after `([^`]*)` add", review)
bad = 0
for s in olds + anchors:
    k = paper.count(s)
    status = "PASS" if k == 1 else "FAIL"
    bad += k != 1
    print(f"{status}: occurs {k} time(s): {s[:90]!r}")
print(f"\n{len(olds) + len(anchors)} quoted texts, {bad} not found exactly once")
sys.exit(1 if bad else 0)
