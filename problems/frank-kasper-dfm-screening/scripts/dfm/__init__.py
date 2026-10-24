"""Diblock foam model (DFM) minimiser for periodic cell partitions.
Single-threaded by design (shared machine): thread pools are pinned to 1."""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
