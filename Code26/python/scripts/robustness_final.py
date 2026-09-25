"""Robustness of the final model at the calibrated parameters (no re-calibration).

Variants: no-assets limit (a_max small), finer asset grid, finer hours grid, alternative
unemployment thresholds s_bar, and no recession cut in the husband's income (phi_rec_H = 1).
For each variant it reports the baseline moments, the effect of the returns-to-experience
experiment at the scale found in the main results, and the change in the cyclical quit gap
(quit rec - quit exp) when the husband's risk is made acyclical (the precautionary channel).

usage: python scripts/robustness_final.py --calib output/final_calib_full.json
           --results output/final_results_full.json [--coarse]
Writes output/robustness_final.md and output/robustness_final.json.
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params
from keam.final.experiments import run, returns_scaled

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--results", required=True, help="JSON written by run_final.py (for the experiment scale)")
ap.add_argument("--coarse", action="store_true")
ap.add_argument("--n-jobs", type=int, default=0)
a = ap.parse_args()
if a.n_jobs:
    os.environ["KEAM_NJOBS"] = str(a.n_jobs)
HERE = os.path.dirname(os.path.abspath(__file__))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p0 = apply_params(base, json.load(open(a.calib))["x"])
s_roe = json.load(open(a.results))["scales"]["returns"]
cfg = SimConfigFinal(N=60, n_cohorts=90)

VARIANTS = {
    "baseline": {},
    "no assets (a_max 0.01, 5 points)": dict(a_max=0.01, nA=5),
    "asset grid 40 points": dict(nA=40),
    "asset grid 40 points, a_max 30": dict(nA=40, a_max=30.0),
    "hours grid 40 points": dict(nH=40),
    "hours grid 40, h_min 0.025": dict(nH=40, h_min=0.025),
    "U threshold s_bar 0.10": dict(s_bar=0.10),
    "U threshold s_bar 0.50": dict(s_bar=0.50),
    "phi_rec_H = 1 (no recession cut in husband income)": dict(phi_rec_H=1.0),
}
KEYS = ["E/pop", "hours|E", "U rate", "quit/m exp", "quit/m rec", "E->nonE/m exp", "E->nonE/m rec",
        "dE/pop rec-exp (pts)", "wage gap (hourly ratio)", "wife share exp", "share Lifecycle", "share PT",
        "share Career", "share NiLF", "cons drop at H job loss rec (%)", "mean assets/monthly HH inc"]
gap = lambda m: 100 * (m["quit/m rec"] - m["quit/m exp"])
res = {}
t0 = time.time()
for name, kw in VARIANTS.items():
    p = p0.replace(**kw)
    m = run(p, cfg)[0]
    m_roe = run(returns_scaled(p, s_roe), cfg)[0]
    m_acyc = run(p.replace(lamH_loss=(p.lamH_loss[0], p.lamH_loss[0]), lamH_find=(p.lamH_find[0], p.lamH_find[0])), cfg)[0]
    res[name] = dict(m={k: float(v) for k, v in m.items()},
                     roe=dict(E=float(m_roe["E/pop"]), quit_gap=gap(m_roe), dE_rec=float(m_roe["dE/pop rec-exp (pts)"])),
                     acyc=dict(quit_gap=gap(m_acyc), dE_rec=float(m_acyc["dE/pop rec-exp (pts)"])))
    print(f"{name}: E {m['E/pop']:.3f} quit gap {gap(m):+.2f} | RoE E {m_roe['E/pop']:.3f} gap {gap(m_roe):+.2f} "
          f"| acyclical H gap {gap(m_acyc):+.2f}  [{time.time() - t0:.0f}s]", flush=True)

names = list(res)
md = ["# Robustness of the final model (calibrated parameters held fixed)\n",
      f"Calibration: `{a.calib}`; returns-to-experience scale x{s_roe:.2f} from `{a.results}`. "
      "Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.\n",
      "## Baseline moments by variant\n", "| moment | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
for k in KEYS:
    md.append(f"| {k} | " + " | ".join(f"{res[n]['m'].get(k, np.nan):.4f}" for n in names) + " |")
md += ["", "## Mechanism and experiment by variant\n", "| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | "
       "acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |", "|---|---|---|---|---|---|---|---|"]
for n in names:
    r = res[n]
    md.append(f"| {n} | {gap(r['m']):+.3f} | {r['m']['dE/pop rec-exp (pts)']:+.3f} | {r['acyc']['quit_gap']:+.3f} | "
              f"{r['acyc']['dE_rec']:+.3f} | {r['roe']['E']:.3f} | {r['roe']['quit_gap']:+.3f} | {r['roe']['dE_rec']:+.3f} |")
md.append(f"\nElapsed {time.time() - t0:.0f}s.\n")
out = os.path.join(HERE, "..", "output", "robustness_final" + ("_coarse" if a.coarse else ""))
open(out + ".md", "w").write("\n".join(md)); json.dump(res, open(out + ".json", "w"), indent=1)
print("saved", out + ".md")
