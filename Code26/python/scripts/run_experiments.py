"""Solve + simulate the baseline, the three comparative statics and the four cohorts.

usage: python scripts/run_experiments.py [--corrected] [--nilfdef 0.3] [--only NAME,...]
Writes output/experiments_<mode>.csv with the moment table.
"""
import sys, os, time, warnings, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd
from keam import Params, Options, solve, simulate, SimConfig
from keam.moments import moments
warnings.simplefilter("ignore")

# (kapscale, wagegapscale, rtoexpscale) exactly as in SimplerMod_May17_shell.m and the
# saved Solution/<name>/paras.mat files
EXPERIMENTS = {
    "Baseline":   dict(kapscale=1.0,   wagegapscale=1.0,  rtoexpscale=1.0),
    "RoE_incr":   dict(kapscale=1.0,   wagegapscale=1.0,  rtoexpscale=1.2),
    "Wgap_dcr":   dict(kapscale=1.0,   wagegapscale=0.8,  rtoexpscale=1.0),
    "Kap_dcr":    dict(kapscale=0.9,   wagegapscale=1.0,  rtoexpscale=1.0),
    "Cohort1950": dict(kapscale=1.01,  wagegapscale=1.07, rtoexpscale=1.095),
    "Cohort1960": dict(kapscale=0.84,  wagegapscale=1.19, rtoexpscale=1.15),
    "Cohort1970": dict(kapscale=1.0,   wagegapscale=1.4,  rtoexpscale=1.35),
    "Cohort1980": dict(kapscale=1.037, wagegapscale=1.41, rtoexpscale=1.37),
}

ap = argparse.ArgumentParser()
ap.add_argument("--corrected", action="store_true")
ap.add_argument("--nilfdef", type=float, default=0.3)
ap.add_argument("--only", type=str, default="")
a = ap.parse_args()
opts = Options.corrected() if a.corrected else Options.faithful()
cfg = SimConfig(NiLFdef=a.nilfdef)
names = a.only.split(",") if a.only else list(EXPERIMENTS)
HERE = os.path.dirname(os.path.abspath(__file__))
rows = {}
for name in names:
    t0 = time.time()
    p = Params(**EXPERIMENTS[name])
    sol = solve(p, opts)
    r = simulate(p, sol, opts, cfg)
    rows[name] = moments(r, sol)
    rows[name].update(dict(tau_wf=p.tau_wf_eff, gam_e=p.gam_e_eff, kapscale=p.kapscale))
    print(f"{name:12s} solved+simulated in {time.time()-t0:.0f}s")
df = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 60)
print(df.to_string(float_format=lambda x: f"{x:.4f}"))
out = os.path.join(HERE, "..", "output", f"experiments_{'corrected' if a.corrected else 'faithful'}.csv")
df.to_csv(out); print("saved", out)
