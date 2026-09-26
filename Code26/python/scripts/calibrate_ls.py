"""Least-squares polishing of the SMM calibration (scipy trust-region reflective with bounds,
finite-difference Jacobian, common simulation seed). Residuals are sqrt(weight) x deviation, so the
sum of squares equals the objective of keam.final.calibrate.

usage: python scripts/calibrate_ls.py --x0 output/final_calib_full.json [--extra lam_f_ratio]
                                      [--max-nfev 6] [--diff-step 0.04] [--tag ls] [--bound name:lo:hi]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import least_squares
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run

ap = argparse.ArgumentParser()
ap.add_argument("--x0", default="output/final_calib_full.json")
ap.add_argument("--extra", type=str, default="")
ap.add_argument("--max-nfev", type=int, default=6)
ap.add_argument("--diff-step", type=float, default=0.04)
ap.add_argument("--tag", type=str, default="ls")
ap.add_argument("--bound", action="append", default=[])
ap.add_argument("--coarse", action="store_true", help="27-type grid (default: 100 types)")
ap.add_argument("--fixed", action="append", default=[], help="fix a FinalParams field (not calibrated): name=value (repeatable)")
ap.add_argument("--set", action="append", default=[])
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
C.BOUNDS.update({"home_young_mult": (1.0, 3.0), "nu_h": (0.3, 0.8), "z_h": (0.3, 0.6), "alpha_h": (0.05, 1.5),
                 "e_max": (1.5, 4.0), "theta_e": (0.005, 0.05), "kappa_h_power": (0.0, 1.0)})
for b in a.bound:
    n, lo, hi = b.split(":"); C.BOUNDS[n] = (float(lo), float(hi))
x0 = json.load(open(os.path.join(ROOT, a.x0)))["x"]
for n in [e for e in a.extra.split(",") if e]:
    x0.setdefault(n, {"lam_f_ratio": 0.85}.get(n, getattr(FinalParams(), n, None)))
for kv in a.set:
    n, v = kv.split("="); x0[n] = float(v)
names = list(x0.keys())
lo = np.array([C.BOUNDS[n][0] for n in names]); hi = np.array([C.BOUNDS[n][1] for n in names])
z0 = np.clip(np.array([x0[n] for n in names]), lo + 1e-9, hi - 1e-9)
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
for kv in a.fixed:
    n, v = kv.split("="); base = base.replace(**{n: type(getattr(base, n))(float(v))})
cfg = SimConfigFinal(N=60, n_cohorts=90)
log = os.path.join(ROOT, "output", f"final_calib_{a.tag}.log"); open(log, "w").close()
hist = []; t0 = time.time()
keys = list(C.TARGETS.keys()); w = np.sqrt([C.WEIGHT[k] for k in keys])


def resid(z):
    x = dict(zip(names, [float(v) for v in z]))
    m, _, _ = run(C.apply_params(base, x), cfg)
    obj, parts = C.objective_from_moments(m)
    r = w * np.array([parts[k] for k in keys])
    hist.append(dict(n=len(hist) + 1, t=round(time.time() - t0), obj=obj, x=x, m={k: float(v) for k, v in m.items()},
                     dev={k: round(v, 3) for k, v in parts.items()}))
    with open(log, "a") as f:
        f.write(json.dumps({k: v for k, v in hist[-1].items() if k != "m"}) + "\n")
    return r


sol = least_squares(resid, z0, bounds=(lo, hi), method="trf", diff_step=a.diff_step, x_scale=hi - lo,
                    max_nfev=a.max_nfev, ftol=1e-6, xtol=1e-6, gtol=1e-6)
best = min(hist, key=lambda h: h["obj"])
out = dict(x=best["x"], obj=best["obj"], moments=best["m"], targets=C.TARGETS, n_eval=len(hist),
           seconds=time.time() - t0, coarse=a.coarse, names=names,
           fixed={kv.split('=')[0]: float(kv.split('=')[1]) for kv in a.fixed}, status=int(sol.status), message=sol.message)
json.dump(out, open(os.path.join(ROOT, "output", f"final_calib_{a.tag}.json"), "w"), indent=1)
print("best objective", round(best["obj"], 4), "after", len(hist), "evaluations;", sol.message)
for k, v in best["x"].items():
    print(f"  {k:16s} {v:.4f}")
for k, tv in C.TARGETS.items():
    print(f"  {k:32s} model {best['m'].get(k, float('nan')):8.4f}  target {tv:8.4f}")
