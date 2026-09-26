"""Exact replication of MATLAB's default generator ``rng(seed)`` (Mersenne Twister).

MATLAB seeds MT19937 with ``init_genrand(seed)`` and converts to doubles with
the 53-bit ``genrand_res53`` recipe, which is what ``numpy.random.RandomState``
uses for ``random_sample``.  numpy's own seeding differs (``init_by_array``),
so the state is set manually.  Verified against ``Code26/Shocks_Types.mat``
(LHsim_i, LWsim_i, Qtoss_i reproduce with zero error).
"""
from __future__ import annotations
import numpy as np


def matlab_twister(seed: int) -> np.random.RandomState:
    mt = np.zeros(624, dtype=np.uint32)
    mt[0] = seed & 0xFFFFFFFF
    for i in range(1, 624):
        prev = int(mt[i - 1])
        mt[i] = (1812433253 * (prev ^ (prev >> 30)) + i) & 0xFFFFFFFF
    rs = np.random.RandomState()
    rs.set_state(("MT19937", mt, 624, 0, 0.0))
    return rs


def rand(seed: int, *shape) -> np.ndarray:
    """MATLAB ``rng(seed); rand(shape...)`` (column-major fill order)."""
    rs = matlab_twister(seed)
    n = int(np.prod(shape)) if shape else 1
    x = rs.random_sample(n)
    if len(shape) <= 1:
        return x.reshape(shape) if shape else x[0]
    return x.reshape(tuple(reversed(shape))).T if len(shape) == 2 else x.reshape(shape, order="F")


def randi(seed: int, imin: int, imax: int, *shape) -> np.ndarray:
    """MATLAB ``rng(seed); randi([imin imax], shape...)`` = imin + floor(range*rand)."""
    u = rand(seed, *shape)
    return (imin + np.floor((imax - imin + 1) * u)).astype(int)
