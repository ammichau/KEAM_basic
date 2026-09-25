"""Run the SMM calibration of the final model's 1940s cohort.

usage: python scripts/calibrate_final.py [--coarse] [--maxfev 150] [--x0 file.json] [--tag name]
Writes output/final_calib_<tag>.log (one JSON line per evaluation) and output/final_calib_<tag>.json (best).
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import run_smm, evaluate, TARGETS

ap = argparse.ArgumentParser()
ap.add_argument("--coarse", action="store_true")
ap.add_argument("--maxfev", type=int, default=150)
ap.add_argument("--x0", type=str, default="")
ap.add_argument("--tag", type=str, default="coarse")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
cfg = SimConfigFinal(N=60, n_cohorts=90)
x0 = dict(mu=1.0, kbar_max=0.075, km_max=2.27, tau_w=0.80, lam_f0=0.40, lam_u0=0.017, lam_u1=0.022,
          ybar_h=0.11, s_bar=0.25, sd_kT=0.30)
if a.x0:
    x0.update(json.load(open(a.x0))["x"])
log = os.path.join(HERE, "..", "output", f"final_calib_{a.tag}.log")
open(log, "w").close()
t0 = time.time()
best, hist, res = run_smm(base, x0, cfg, log, maxfev=a.maxfev)
out = dict(x=best["x"], obj=best["obj"], moments=best["m"], targets=TARGETS, n_eval=len(hist),
           seconds=time.time() - t0, coarse=a.coarse)
json.dump(out, open(os.path.join(HERE, "..", "output", f"final_calib_{a.tag}.json"), "w"), indent=1)
print("best objective", best["obj"], "after", len(hist), "evaluations")
for k, v in best["x"].items():
    print(f"  {k:10s} {v:.4f}")
for k, tv in TARGETS.items():
    print(f"  {k:32s} model {best['m'].get(k, float('nan')):8.4f}  target {tv:8.4f}")
