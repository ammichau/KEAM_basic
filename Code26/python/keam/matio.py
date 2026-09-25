"""Loaders for the saved MATLAB solution and output files."""
from __future__ import annotations
import os
import numpy as np
import scipy.io as sio
from .params import Params
from .solve import Solution


def load_matlab_solution(folder: str) -> Solution:
    """Load paras/policies/Vfuns from a Code26/Solution/<name> folder."""
    p = Params.from_paras_mat(os.path.join(folder, "paras.mat"))
    pol = sio.loadmat(os.path.join(folder, "policies.mat"))
    vf = sio.loadmat(os.path.join(folder, "Vfuns.mat"))
    return Solution(VE=vf["VE"], VU=vf["VU"], V=vf["V"], gH=pol["gH"], gS=pol["gS"],
                    gQ=pol["gQ"].astype(np.int8), n_iter=np.zeros((p.nI, p.nT), int),
                    last_err=np.full((p.nI, p.nT), np.nan), params=p, options=None,
                    info={"source": folder})
