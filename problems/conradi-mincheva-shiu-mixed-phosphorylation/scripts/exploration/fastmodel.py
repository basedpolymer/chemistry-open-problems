# Exploration script (floating point; not a certificate). Environment: WSL Fedora, Python 3.12, numpy, scipy, numba.
# Run from scripts/exploration/ (it imports the sibling modules). Typical runtime: seconds to minutes.
"""Fast (numba) right-hand side, Jacobian, and variational system of the 9-species model."""
import numpy as np
from numba import njit


@njit(cache=True)
def g9(x, k):
    r1 = k[0] * x[0] * x[1]; r2 = k[1] * x[2]; r3 = k[2] * x[2]; r4 = k[3] * x[3]
    r5 = k[4] * x[4] * x[5]; r6 = k[5] * x[6]; r7 = k[6] * x[6]; r8 = k[7] * x[7] * x[5]
    r9 = k[8] * x[8]; r10 = k[9] * x[8]
    out = np.empty(9)
    out[0] = -r1 + r2 + r10
    out[1] = -r1 + r2 + r4
    out[2] = r1 - r2 - r3
    out[3] = r3 - r4
    out[4] = r4 - r5 + r6
    out[5] = -r5 + r6 + r7 - r8 + r9 + r10
    out[6] = r5 - r6 - r7
    out[7] = r7 - r8 + r9
    out[8] = r8 - r9 - r10
    return out


@njit(cache=True)
def J9(x, k):
    J = np.zeros((9, 9))
    # r1 = k1 x1 x2
    d1 = np.zeros(9); d1[0] = k[0] * x[1]; d1[1] = k[0] * x[0]
    d5 = np.zeros(9); d5[4] = k[4] * x[5]; d5[5] = k[4] * x[4]
    d8 = np.zeros(9); d8[7] = k[7] * x[5]; d8[5] = k[7] * x[7]
    for j in range(9):
        J[0, j] += -d1[j]; J[1, j] += -d1[j]; J[2, j] += d1[j]
        J[4, j] += -d5[j]; J[5, j] += -d5[j]; J[6, j] += d5[j]
        J[7, j] += -d8[j]; J[5, j] += -d8[j]; J[8, j] += d8[j]
    # linear reactions
    J[0, 2] += k[1]; J[1, 2] += k[1]; J[2, 2] += -k[1]
    J[2, 2] += -k[2]; J[3, 2] += k[2]
    J[3, 3] += -k[3]; J[4, 3] += k[3]; J[1, 3] += k[3]
    J[6, 6] += -k[5]; J[4, 6] += k[5]; J[5, 6] += k[5]
    J[6, 6] += -k[6]; J[7, 6] += k[6]; J[5, 6] += k[6]
    J[8, 8] += -k[8]; J[7, 8] += k[8]; J[5, 8] += k[8]
    J[8, 8] += -k[9]; J[0, 8] += k[9]; J[5, 8] += k[9]
    return J


@njit(cache=True)
def var_rhs(t, y, k):
    """y = (x (9), Phi (9x9 flattened)); full-space variational equation."""
    x = y[:9]
    J = J9(x, k)
    out = np.empty(90)
    out[:9] = g9(x, k)
    Phi = y[9:].reshape(9, 9)
    out[9:] = (J @ Phi).ravel()
    return out
