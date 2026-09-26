"""Coarse (27-type) SMM including the child-care home-production multiplier and the
home-production parameters, seeded from the best point of scripts/explore_childcare.py.
Bounds and optional parameters are injected at run time so keam/final/calibrate.py is untouched.

usage: python scripts/calibrate_childcare.py [--maxfev 180] [--tag childcare]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C

ap = argparse.ArgumentParser()
ap.add_argument("--maxfev", type=int, default=180)
ap.add_argument("--tag", type=str, default="childcare")
ap.add_argument("--x0", type=str, default="")
ap.add_argument("--extra", type=str, default="", help="extra calibrated FinalParams fields, comma-separated (e.g. alpha_h)")
ap.add_argument("--full", action="store_true", help="100-type grid (default: 27 types)")
ap.add_argument("--bound", action="append", default=[], help="override a bound: name:lo:hi (repeatable)")
ap.add_argument("--set", action="append", default=[], help="override a starting value: name=value (repeatable)")
ap.add_argument("--maxfev-note", type=str, default="")
ap.add_argument("--n-jobs", type=int, default=0)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
C.BOUNDS.update({"home_young_mult": (1.0, 3.0), "nu_h": (0.3, 0.8), "z_h": (0.3, 0.6), "alpha_h": (0.05, 1.5),
                 "e_max": (1.5, 4.0), "theta_e": (0.005, 0.05)})
extra = [n for n in a.extra.split(",") if n]
C.BOUNDS.update({"kappa_h_power": (0.0, 1.0)})
for b in a.bound:
    n, lo, hi = b.split(":"); C.BOUNDS[n] = (float(lo), float(hi))
for n in ["home_young_mult", "nu_h", "z_h"] + extra:
    if n not in C.OPTIONAL_PARAMS:
        C.OPTIONAL_PARAMS.append(n)
names = C.PARAM_NAMES + ["home_young_mult", "nu_h", "z_h"] + extra
if a.n_jobs:
    os.environ["KEAM_NJOBS"] = str(a.n_jobs)
base = FinalParams() if a.full else FinalParams(n_omega=3, n_kbar=3, n_km=3)
cfg = SimConfigFinal(N=60, n_cohorts=90)
x0 = json.load(open(os.path.join(HERE, "..", "output", "x0_km14.json")))["x"]
x0.update(km_max=4.0, home_young_mult=1.5, nu_h=0.65, z_h=0.45)
for n in extra:
    x0.setdefault(n, getattr(FinalParams(), n))
if a.x0:
    x0.update(json.load(open(a.x0))["x"])
for kv in a.set:
    n, v = kv.split("="); x0[n] = float(v)
log = os.path.join(HERE, "..", "output", f"final_calib_{a.tag}.log"); open(log, "w").close()
t0 = time.time()
best, hist, res = C.run_smm(base, x0, cfg, log, maxfev=a.maxfev, names=names)
out = dict(x=best["x"], obj=best["obj"], moments=best["m"], targets=C.TARGETS, n_eval=len(hist),
           seconds=time.time() - t0, coarse=not a.full, names=names)
json.dump(out, open(os.path.join(HERE, "..", "output", f"final_calib_{a.tag}.json"), "w"), indent=1)
print("best objective", round(best["obj"], 3), "after", len(hist), "evaluations")
for k, v in best["x"].items():
    print(f"  {k:16s} {v:.4f}")
for k, tv in C.TARGETS.items():
    print(f"  {k:32s} model {best['m'].get(k, float('nan')):8.4f}  target {tv:8.4f}")
