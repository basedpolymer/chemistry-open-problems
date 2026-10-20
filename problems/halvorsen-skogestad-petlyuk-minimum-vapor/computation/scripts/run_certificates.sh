#!/bin/bash
# Runs the PASS/FAIL certificate scripts of Agent C (exact arithmetic, combinatorics, and one seeded numerical check).
# Environment: WSL Fedora, Python 3 with mpmath (from the project root: outils/wsl "bash <this script>").
# Runtime: about 30 s (WSL Fedora, 1 thread; pyomo + ipopt needed only for check_lemmas.py).  Exit code non-zero if any certificate fails.
set -e
cd "$(dirname "$0")"
python -u export_configs.py         # enumeration counts 8 / 152 / 6128 (Shah & Agrawal 2010)
python -u verify_identity_exact.py  # cut identity, exact rationals, all configurations N = 3, 4, 5
python -u verify_sign_paths.py      # no sign-consistent light->heavy path, strict model, N = 3..6
python -u check_lemmas.py 4 20 7 0   # numerical lemma checks (sign rule, recovery ordering, S5, sign graph)
echo "ALL CERTIFICATES PASSED"
