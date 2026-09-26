"""Evaluate one calibration file on the 100-type grid (or coarse) and print the objective and the
targeted moments; optionally append the moments to the calibration JSON under "moments_full".

usage: python scripts/eval_point.py --calib output/final_calib_rho_coarse.json [--coarse] [--save]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--coarse", action="store_true")
ap.add_argument("--save", action="store_true")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
path = a.calib if os.path.isabs(a.calib) else os.path.join(ROOT, a.calib)
calib = json.load(open(path))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p = C.params_from_calib(calib, base)
t0 = time.time(); m, _, _ = run(p, SimConfigFinal(N=60, n_cohorts=90)); obj, parts = C.objective_from_moments(m)
print(f"{a.calib}: {'coarse' if a.coarse else 'full'} grid, objective {obj:.4f} ({time.time() - t0:.0f}s)")
for k, tv in C.TARGETS.items():
    print(f"  {k:28s} model {m[k]:8.4f} target {tv:8.4f} dev {parts[k]:+.3f}")
for k in ["U rate", "wife share exp", "cons drop at H job loss exp (%)", "mean assets/monthly HH inc"]:
    print(f"  {k:28s} model {m[k]:8.4f}")
if a.save:
    calib["moments_" + ("coarse" if a.coarse else "full")] = {k: float(v) for k, v in m.items()}
    calib["obj_" + ("coarse" if a.coarse else "full")] = obj
    json.dump(calib, open(path, "w"), indent=1); print("saved")
