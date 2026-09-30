"""Local elasticities of the targeted moments with respect to the calibrated parameters
(one-sided finite differences on the 100-type grid, common simulation seed).
Output: output/jacobian_<tag>.json (moments at every point) and output/jacobian_<tag>.md
(elasticity table: percent change of the moment per percent change of the parameter; for the
recession employment drop the entry is the change in percentage points per percent change).

usage: python scripts/jacobian_final.py [--calib output/final_calib_full.json] [--step 0.05]
                                        [--extra lam_f_ratio] [--tag final]
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
ap.add_argument("--step", type=float, default=0.05)
ap.add_argument("--extra", type=str, default="lam_f_ratio")
ap.add_argument("--tag", type=str, default="final")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
x0 = json.load(open(os.path.join(ROOT, a.calib)))["x"]
for n in [e for e in a.extra.split(",") if e]:
    x0.setdefault(n, {"lam_f_ratio": 0.85}.get(n, getattr(FinalParams(), n, None)))
names = list(x0.keys())
cfg = SimConfigFinal(N=60, n_cohorts=90)
out_json = os.path.join(ROOT, "output", f"jacobian_{a.tag}.json")
res = json.load(open(out_json)) if os.path.exists(out_json) else {"x0": x0, "step": a.step, "points": {}}


def evaluate(key, x):
    if key in res["points"]:
        return res["points"][key]["m"]
    t0 = time.time(); m, _, _ = run(C.apply_params(FinalParams(), x), cfg)
    res["points"][key] = {"x": x, "m": {k: float(v) for k, v in m.items()}, "seconds": time.time() - t0}
    json.dump(res, open(out_json, "w"), indent=1)
    print(key, f"{time.time() - t0:.0f}s", {k: round(m[k], 4) for k in C.TARGETS}, flush=True)
    return res["points"][key]["m"]


m0 = evaluate("base", x0)
for n in names:
    x = dict(x0); x[n] = x0[n] * (1 + a.step)
    evaluate(n, x)
# elasticity table
keys = list(C.TARGETS.keys())
L = [f"# Elasticities of the targeted moments at `{a.calib}` (step +{100*a.step:.0f}%)", "",
     "Entry: percent change of the moment per percent change of the parameter "
     "(for the recession employment drop: percentage points per percent).", "",
     "| parameter | " + " | ".join(keys) + " |", "|---|" + "---|" * len(keys)]
E = {}
for n in names:
    m1 = res["points"][n]["m"]; row = []
    for k in keys:
        if k == "dE/pop rec-exp (pts)":
            e = (m1[k] - m0[k]) / (100 * a.step)
        else:
            e = ((m1[k] - m0[k]) / m0[k]) / a.step
        row.append(e)
    E[n] = dict(zip(keys, row))
    L.append(f"| {n} | " + " | ".join(f"{v:+.2f}" for v in row) + " |")
res["elasticities"] = E
json.dump(res, open(out_json, "w"), indent=1)
open(os.path.join(ROOT, "output", f"jacobian_{a.tag}.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
