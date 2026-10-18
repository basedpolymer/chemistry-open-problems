#!/usr/bin/env python3
"""Build a manuscript PDF from a problem's Markdown write-up with the TeXtured template.

Usage:  python templates/build_pdf.py problems/<slug>/notes.md [--out paper.pdf] [--subtitle "..."] [--keep]

- The single "# " heading is the title; "## " headings become sections (numbered unless a heading carries
  `.unnumbered`).
- Math in \\( \\) and \\[ \\] (pandoc's tex_math_single_backslash), as in the notes.
- Output: `paper.pdf` next to the Markdown file unless --out is given.
Adapted from the M2SU `cr_pdf.py` script (TeXtured 1.5.1 template by J. Dujava, CC0), English version.
Requires pandoc 3.x and pdflatex (MiKTeX or TeX Live) with the packages listed in textured-cr/paquets-texlive.txt.
"""
import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "textured-cr"
READER = ("markdown+tex_math_single_backslash+lists_without_preceding_blankline"
          "-tex_math_dollars-raw_html-auto_identifiers")


def find_pandoc():
    for p in (shutil.which("pandoc"), os.path.expandvars(r"%LOCALAPPDATA%\Pandoc\pandoc.exe"),
              str(Path.home() / ".local/bin/pandoc")):
        if p and Path(p).exists():
            return p
    sys.exit("pandoc not found")


def build(md: Path, out: Path, subtitle: str, keep: bool):
    text = md.read_text(encoding="utf-8")
    titles = re.findall(r"^# (.+)$", text, flags=re.M)
    if len(titles) != 1:
        sys.exit(f"exactly one '# ' title expected, found {len(titles)}")
    title = titles[0].strip()
    short = re.split(r"\s[—:–]\s|: ", title)[0]
    plain = re.sub(r"[\\{}$]", "", title)
    pandoc = find_pandoc()
    tmp = Path(tempfile.mkdtemp(prefix="paper_"))
    (tmp / "source.md").write_text(text, encoding="utf-8")
    env = dict(os.environ)
    env["TEXINPUTS"] = os.pathsep.join([".", TEMPLATE.as_posix() + "//", ""])
    cmd = [pandoc, "source.md", "-o", "document.tex", "--standalone", "-f", READER,
           "--shift-heading-level-by=-1", "--template", str(TEMPLATE / "cr.tex"),
           "-V", f"figdir={md.parent.as_posix()}", "-V", f"short-title={short}",
           "-V", f"title-plain={plain}", "-V", "author-meta=", "-V", "subject=", "-V", "keywords=",
           "-V", f"subtitle={subtitle}", "-V", "numbered=true"]
    subprocess.run(cmd, cwd=tmp, env=env, check=True)
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "document.tex"],
                           cwd=tmp, env=env, capture_output=True, text=True, errors="replace")
        if r.returncode != 0:
            log = (tmp / "document.log").read_text(errors="replace")
            err = re.findall(r"^!.*(?:\n.*){0,6}", log, flags=re.M)
            sys.exit("pdflatex failed:\n" + "\n".join(err[:5]) + f"\n(log: {tmp / 'document.log'})")
    log = (tmp / "document.log").read_text(errors="replace")
    shutil.copy(tmp / "document.pdf", out)
    pages = re.search(r"Output written on document\.pdf \((\d+) pages?", log)
    over = len(re.findall(r"^Overfull \\hbox", log, flags=re.M))
    missing = sorted(set(re.findall(r"Missing character: There is no (\S+)", log)))
    print(f"PDF written: {out} ({pages.group(1) if pages else '?'} pages)")
    if over:
        print(f"warning: {over} overfull line(s)")
    if missing:
        print("warning: missing glyphs: " + " ".join(missing))
    if keep or over or missing:
        print(f"build directory: {tmp}")
    else:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md")
    ap.add_argument("--out")
    ap.add_argument("--subtitle")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    md = Path(a.md).resolve()
    out = Path(a.out).resolve() if a.out else md.with_name("paper.pdf")
    d = datetime.date.today()
    subtitle = a.subtitle or (f"chemistry-open-problems \\textperiodcentered{{}} "
                              f"{d.day} {d.strftime('%B')} {d.year}")
    build(md, out, subtitle, a.keep)


if __name__ == "__main__":
    main()
