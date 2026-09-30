"""Solve the model in faithful mode and compare with the saved MATLAB solution.

usage: python scripts/verify_solution.py [Solution subfolder, default Baseline]
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, scipy.io as sio
from keam.params import Params, Options
from keam.solve import solve

sub = sys.argv[1] if len(sys.argv) > 1 else "Baseline"
root = os.path.join(os.path.dirname(__file__), "..", "..", "Solution", sub)
p = Params.from_paras_mat(os.path.join(root, "paras.mat"))
pol = sio.loadmat(os.path.join(root, "policies.mat"))
vf = sio.loadmat(os.path.join(root, "Vfuns.mat"))

t0 = time.time()
sol = solve(p, Options.faithful(), verbose=("-v" in sys.argv))
print(f"solved in {time.time()-t0:.1f}s; VFI iterations per (type, age):")
print(sol.n_iter[:, :3].T)

def cmp(name, a, b, working_only=True):
    if working_only:
        a, b = a[:, :3], b[:, :3]
    d = np.abs(a - b)
    print(f"{name:3s}: max|diff| = {d.max():.3e}  mean|diff| = {d.mean():.3e}  "
          f"share exact(<1e-9) = {np.mean(d < 1e-9):.4f}  range MATLAB [{b.min():.4g}, {b.max():.4g}]")
    return d

cmp("gH", sol.gH, pol["gH"])
cmp("gS", sol.gS, pol["gS"])
dq = cmp("gQ", sol.gQ.astype(float), pol["gQ"].astype(float))
print(f"     gQ disagreements: {int((dq > 0).sum())} of {dq.size} states")
cmp("VE", sol.VE, vf["VE"])
cmp("VU", sol.VU, vf["VU"])
cmp("V", sol.V, vf["V"])
np.savez(os.path.join(os.path.dirname(__file__), "..", f"solution_faithful_{sub}.npz"),
         VE=sol.VE, VU=sol.VU, V=sol.V, gH=sol.gH, gS=sol.gS, gQ=sol.gQ, n_iter=sol.n_iter)
