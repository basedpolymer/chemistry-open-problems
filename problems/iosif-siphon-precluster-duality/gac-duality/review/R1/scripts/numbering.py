#!/usr/bin/env python3
"""R1 helper: list the theorem-like environments, tables and figures of paper.tex in order with their labels and
source lines, numbered with one shared counter (equations share it: checked against the numbers quoted in writer/CHANGES.md), to map labels to printed numbers.
Standard library only.  Usage: py -3 scripts/numbering.py ../../../paper.tex"""
import re
import sys

src = open(sys.argv[1], encoding="utf-8").read()
envs = "equation|lemma|theorem|proposition|corollary|definition|example|remark|table|figure|longtable|conjecture|question"
pat = re.compile(r"\\begin\{(" + envs + r")\}")
n = 0
for m in pat.finditer(src):
    env = m.group(1)
    line = src.count("\n", 0, m.start()) + 1
    lab = re.search(r"\\label\{([^}]*)\}", src[m.end():m.end() + 900])
    n += 1
    print(n, env, lab.group(1) if lab else "-", "line", line)
