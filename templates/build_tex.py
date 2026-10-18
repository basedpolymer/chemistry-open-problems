#!/usr/bin/env python3
"""Compile a manuscript written directly in LaTeX with the shared TeXtured preamble.

Usage:  python templates/build_tex.py problems/<slug>/paper.tex [--keep]

Runs latexmk (pdflatex + biber) in a temporary directory with TEXINPUTS pointing to templates/textured-cr/,
then copies paper.pdf next to paper.tex. Requires a TeX distribution with latexmk, pdflatex, biber and the
packages listed in templates/textured-cr/paquets-texlive.txt (MiKTeX installs missing packages on the fly).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parent / "textured-cr"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tex")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    tex = Path(a.tex).resolve()
    src = tex.parent
    tmp = Path(tempfile.mkdtemp(prefix="paper_"))
    env = dict(os.environ)
    sep = os.pathsep
    env["TEXINPUTS"] = sep.join([src.as_posix(), TEMPLATE.as_posix() + "//", ""])
    env["BIBINPUTS"] = sep.join([src.as_posix(), TEMPLATE.as_posix() + "//", ""])
    cmd = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
           f"-outdir={tmp.as_posix()}", tex.name]
    r = subprocess.run(cmd, cwd=src, env=env, capture_output=True, text=True, errors="replace")
    log_path = tmp / (tex.stem + ".log")
    log = log_path.read_text(errors="replace") if log_path.exists() else ""
    pdf = tmp / (tex.stem + ".pdf")
    if r.returncode != 0 or not pdf.exists():
        errs = re.findall(r"^.*:\d+:.*(?:\n.*){0,4}", log, flags=re.M) or re.findall(r"^!.*(?:\n.*){0,4}", log, flags=re.M)
        print("\n".join(errs[:5]) or r.stdout[-3000:])
        sys.exit(f"compilation failed (build directory: {tmp})")
    out = src / (tex.stem + ".pdf")
    shutil.copy(pdf, out)
    pages = re.search(r"\((\d+) pages?,", log.replace("\n", ""))
    if shutil.which("pdfinfo"):
        info = subprocess.run(["pdfinfo", str(out)], capture_output=True, text=True, errors="replace").stdout
        pages = re.search(r"^Pages:\s+(\d+)", info, flags=re.M) or pages
    over = len(re.findall(r"^Overfull \\hbox", log, flags=re.M))
    undef = re.findall(r"(?:Citation|Reference) [`'](\S+)' .*undefined", log)
    missing = sorted(set(re.findall(r"Missing character: There is no (\S+)", log)))
    print(f"PDF written: {out} ({pages.group(1) if pages else '?'} pages)")
    if over:
        print(f"warning: {over} overfull line(s)")
    if undef:
        print("warning: undefined references/citations: " + ", ".join(sorted(set(undef))))
    if missing:
        print("warning: missing glyphs: " + " ".join(missing))
    if a.keep or over or undef or missing:
        print(f"build directory: {tmp}")
    else:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
