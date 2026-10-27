#!/usr/bin/env python3
"""R1 helper: search a PDF for path fragments, in the raw bytes and inside every Flate-compressed stream (object
streams included).  Standard library only (re, zlib).  Usage: py -3 scripts/pdf_path_check.py file.pdf [pattern ...]
Default patterns: the user-name fragment given with --user (not stored here), "Users", "PTEX.FileName".
Exit 1 if any pattern is found."""
import re
import sys
import zlib

args = [a for a in sys.argv[1:] if not a.startswith("--user=")]
user = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--user=")]
path, pats = args[0], (args[1:] or ["Users", "PTEX.FileName"]) + user
raw = open(path, "rb").read()
chunks = [raw]
n_streams = n_ok = 0
for m in re.finditer(rb"stream\r?\n", raw):
    start = m.end()
    end = raw.find(b"endstream", start)
    if end < 0:
        continue
    n_streams += 1
    try:
        chunks.append(zlib.decompress(raw[start:end]))
        n_ok += 1
    except zlib.error:
        try:
            chunks.append(zlib.decompressobj().decompress(raw[start:end]))
            n_ok += 1
        except zlib.error:
            pass
found = 0
for p in pats:
    k = sum(c.count(p.encode()) for c in chunks)
    found += k
    print(f"{'FOUND' if k else 'absent'}: {p!r} ({k} occurrence(s))")
print(f"{n_streams} streams, {n_ok} decompressed")
sys.exit(1 if found else 0)
