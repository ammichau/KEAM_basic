"""Effect of persistence in the cost-of-work shock (rho_kT: monthly probability that the shock keeps
its value) at the calibrated parameters. First an iid check on the 100-type grid (must reproduce the
calibrated moments up to the age-transition treatment), then a coarse-grid (27 types) scan over
rho_kT with a 3-node shock. Output: output/explore_persistence.json/.md.

usage: python scripts/explore_persistence.py [--calib output/final_calib_full.json]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_full.json")
ap.add_argument("--rhos", default="0,0.5,0.8,0.9,0.95")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
calib = json.load(open(os.path.join(ROOT, a.calib))); x = calib["x"]
cfg = SimConfigFinal(N=60, n_cohorts=90)
KEYS = list(C.TARGETS) + ["U rate", "n careers"]
rows = {}


def add(name, p):
    t0 = time.time(); m, _, _ = run(p, cfg); obj, _ = C.objective_from_moments(m)
    rows[name] = {"moments": {k: float(m[k]) for k in m}, "obj": obj, "seconds": time.time() - t0}
    print(name, f"{time.time() - t0:.0f}s obj {obj:.3f}", {k: round(float(m[k]), 4) for k in C.TARGETS}, flush=True)
    json.dump({"calib": a.calib, "rows": rows}, open(os.path.join(ROOT, "output", "explore_persistence.json"), "w"), indent=1)


add("full grid, iid 5 nodes (check)", C.apply_params(FinalParams(), x))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3)
add("coarse, iid 5 nodes", C.apply_params(base, x))
for r in [float(v) for v in a.rhos.split(",")]:
    add(f"coarse, 3 nodes, rho={r}", C.apply_params(base, x).replace(n_kT=3, rho_kT=r))
L = ["# Persistence of the cost-of-work shock at the calibrated parameters", "",
     f"Calibrated point `{a.calib}`; the first row is the 100-type iid check (calibrated moments: "
     + ", ".join(f"{k} {calib['moments'][k]:.4f}" for k in ["E/pop", "share NiLF", "quit/m exp", "quit/m rec"]) + ").", "",
     "| variant | obj | " + " | ".join(KEYS) + " |", "|---|---|" + "---|" * len(KEYS)]
for n, r in rows.items():
    L.append(f"| {n} | {r['obj']:.3f} | " + " | ".join(f"{r['moments'].get(k, float('nan')):.4f}" for k in KEYS) + " |")
open(os.path.join(ROOT, "output", "explore_persistence.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
