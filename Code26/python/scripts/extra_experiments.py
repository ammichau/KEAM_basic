"""Supplementary single-factor experiments sized to the 1970s employment rate (0.73):
  (a) child-care cost down: home_young_mult -> 1 + s (m_c - 1), s in [0, 1]
  (b) total cost down: kbar_max * s and km_max -> 1 + s (km_max - 1) jointly
usage: python scripts/extra_experiments.py --calib output/final_calib_full.json [--coarse]
Writes output/extra_experiments_<full|coarse>.md and .json
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params, params_from_calib
from keam.final.experiments import run, size_to_employment

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--coarse", action="store_true")
ap.add_argument("--n-jobs", type=int, default=0)
ap.add_argument("--tag", default="", help="output name suffix (default: full or coarse)")
a = ap.parse_args()
if a.n_jobs:
    os.environ["KEAM_NJOBS"] = str(a.n_jobs)
HERE = os.path.dirname(os.path.abspath(__file__))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p = params_from_calib(json.load(open(a.calib)), base)
tagsfx = ("_" + a.tag) if a.tag else ""
cfg = SimConfigFinal(N=60, n_cohorts=90)
E70 = 0.73
t0 = time.time()
base_m, _, _ = run(p, cfg)
ex = {}; scales = {}
def childcare(pp, s):
    return pp.replace(home_young_mult=1.0 + s * (pp.home_young_mult - 1.0))
def cost_all(pp, s):
    return pp.replace(kbar_max=pp.kbar_max * s, km_max=1.0 + s * (pp.km_max - 1.0))
s1, m1, _ = size_to_employment(childcare, p, cfg, E70, 0.0, 1.0)
ex[f"child-care cost x{s1:.2f}"] = m1; scales["childcare"] = s1
json.dump(dict(scales=scales, experiments={k: {kk: float(vv) for kk, vv in v.items() if isinstance(vv, (int, float, np.floating))} for k, v in ex.items()}),
          open(os.path.join(HERE, "..", "output", f"extra_experiments_partial{tagsfx}.json"), "w"), indent=1)
s2, m2, _ = size_to_employment(cost_all, p, cfg, E70, 0.02, 1.0)
ex[f"all costs x{s2:.2f}"] = m2; scales["cost_all"] = s2
KEYS = ["E/pop", "hours|E", "U rate", "quit/m exp", "quit/m rec", "E->nonE/m exp", "E->nonE/m rec",
        "dE/pop rec-exp (pts)", "wife share exp", "wife share rec", "wage gap (hourly ratio)",
        "share Lifecycle", "share PT", "share Career", "share NiLF", "cons drop at H job loss exp (%)",
        "cons drop at H job loss rec (%)", "mean assets/monthly HH inc"]
rows = {"baseline": base_m, **ex}
md = ["### Supplementary experiments sized to the 1970s employment rate (0.73)", "",
      "| moment | " + " | ".join(rows) + " |", "|---|" + "---|" * len(rows)]
for k in KEYS:
    md.append(f"| {k} | " + " | ".join(f"{r.get(k, np.nan):.4f}" for r in rows.values()) + " |")
md.append(f"\nElapsed {time.time()-t0:.0f}s.\n")
tag = a.tag or ("coarse" if a.coarse else "full")
open(os.path.join(HERE, "..", "output", f"extra_experiments_{tag}.md"), "w").write("\n".join(md))
fl = lambda d: {k: float(v) for k, v in d.items() if isinstance(v, (int, float, np.floating))}
json.dump(dict(baseline=fl(base_m), scales=scales, experiments={k: fl(v) for k, v in ex.items()}),
          open(os.path.join(HERE, "..", "output", f"extra_experiments_{tag}.json"), "w"), indent=1)
print("\n".join(md))
