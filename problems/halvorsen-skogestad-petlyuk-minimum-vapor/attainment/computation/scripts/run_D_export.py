#!/usr/bin/env python3
"""Wrapper (Agent C): run Agent D's design/scripts/export_config.py read-only with its own arguments.  D's
check_design.py parses sys.argv when imported, so sys.argv is reset around that import.  Output goes where the
arguments say (here: ../runs/D_exports/).  Usage: python run_D_export.py CASE out.json dps"""
import os
import sys
import runpy
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'design', 'scripts')
sys.path.insert(0, D)
args = sys.argv[1:]
sys.argv = ['check_design.py']
import check_design  # noqa: E402,F401
sys.argv = ['export_config.py'] + args
runpy.run_path(os.path.join(D, 'export_config.py'), run_name='__main__')
