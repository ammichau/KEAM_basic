"""Check of the scale-invariant cost shock (FinalParams.kT_mult): with KPR preferences and gamma = 1 the
one-month proportional tax exp(-kappa_T) inside the aggregator is the additive shock of the log model, so the
two solutions must give the same simulated moments. Also reports the moments and solve time at gamma = 2.
Output: output/verify_kT_mult.json/.md

usage: python scripts/verify_kT_mult.py [--calib output/final_calib_v7b_full.json] [--coarse]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_v7b_full.json")
ap.add_argument("--coarse", action="store_true")
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
calib = json.load(open(os.path.join(ROOT, a.calib)))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p0 = C.params_from_calib(calib, base)
cfg = SimConfigFinal(N=60, n_cohorts=90)
keys = list(C.TARGETS) + ["U rate"]
out = {}
for label, kw in [("gamma 1, additive shock", dict(kpr=True, gamma=1.0, kT_mult=False)),
                  ("gamma 1, proportional shock", dict(kpr=True, gamma=1.0, kT_mult=True)),
                  ("gamma 2, additive shock", dict(kpr=True, gamma=2.0, kT_mult=False)),
                  ("gamma 2, proportional shock", dict(kpr=True, gamma=2.0, kT_mult=True))]:
    t0 = time.time(); m = run(p0.replace(**kw), cfg)[0]
    out[label] = dict(seconds=time.time() - t0, m={k: float(m[k]) for k in keys if k in m})
    print(label, f"{out[label]['seconds']:.0f}s", {k: round(v, 4) for k, v in out[label]["m"].items()}, flush=True)
d1 = max(abs(out["gamma 1, additive shock"]["m"][k] - out["gamma 1, proportional shock"]["m"][k]) for k in out["gamma 1, additive shock"]["m"])
res = dict(calib=a.calib, coarse=a.coarse, max_abs_diff_gamma1=d1, runs=out)
json.dump(res, open(os.path.join(ROOT, "output", "verify_kT_mult.json"), "w"), indent=1)
L = [f"# Proportional cost shock (kT_mult) check at `{a.calib}` ({'27' if a.coarse else '100'} types, sd_kT unchanged)", "",
     f"gamma = 1: largest absolute moment difference between the additive and the proportional shock {d1:.2e} (must be ~0).", "",
     "| moment | " + " | ".join(out) + " |", "|---|" + "---|" * len(out)]
for k in keys:
    L.append(f"| {k} | " + " | ".join(f"{r['m'].get(k, np.nan):.4f}" for r in out.values()) + " |")
L.append("| solve+simulate seconds | " + " | ".join(f"{r['seconds']:.0f}" for r in out.values()) + " |")
open(os.path.join(ROOT, "output", "verify_kT_mult.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
