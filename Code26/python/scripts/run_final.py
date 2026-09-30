"""Produce the final-model results: baseline moments, single-factor experiments sized to the
1970s employment rate, cohort accounting, and mechanism counterfactuals.

usage: python scripts/run_final.py --calib output/final_calib_coarse.json [--full] [--out RESULTS.md]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params, TARGETS, params_from_calib
from keam.final.experiments import run, compensated_wage_gap, cost_scaled, returns_scaled, size_to_employment, counterfactuals

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--full", action="store_true", help="100 types (default: coarse 27)")
ap.add_argument("--out", default="")
ap.add_argument("--n-jobs", type=int, default=0)
ap.add_argument("--tag", default="", help="output name suffix (default: full or coarse)")
a = ap.parse_args()
if a.n_jobs:
    os.environ["KEAM_NJOBS"] = str(a.n_jobs)
HERE = os.path.dirname(os.path.abspath(__file__))
calib = json.load(open(a.calib))
base = FinalParams() if a.full else FinalParams(n_omega=3, n_kbar=3, n_km=3)
p = params_from_calib(calib, base)
cfg = SimConfigFinal(N=60, n_cohorts=90)
KEYS = ["E/pop", "hours|E", "U rate", "quit/m exp", "quit/m rec", "E->nonE/m exp", "E->nonE/m rec",
        "dE/pop rec-exp (pts)", "wife share exp", "wife share rec", "wage gap (hourly ratio)",
        "share Lifecycle", "share PT", "share Career", "share NiLF", "HH income rec/exp - 1 (%)",
        "cons drop at H job loss exp (%)", "cons drop at H job loss rec (%)", "mean assets/monthly HH inc"]

def table(rows: dict, title: str):
    lines = [f"### {title}", "", "| moment | " + " | ".join(rows) + " |", "|---|" + "---|" * len(rows)]
    for k in KEYS:
        lines.append(f"| {k} | " + " | ".join(f"{r.get(k, float('nan')):.4f}" for r in rows.values()) + " |")
    return "\n".join(lines) + "\n"

t0 = time.time(); md = [f"# Final model results ({'100' if a.full else '27'} types)\n"]
base_m, _, _ = run(p, cfg)
rows = {"target": TARGETS, "baseline 1940s": base_m}
md.append(table(rows, "Baseline calibration")); print(md[-1])

# ---- single-factor experiments sized to the 1970s employment rate (slides p.35: 0.73)
E70 = 0.73
ex = {}
s_roe, m_roe, _ = size_to_employment(lambda pp, s: returns_scaled(pp, s), p, cfg, E70, 1.0, 2.5)
ex[f"RoE x{s_roe:.2f}"] = m_roe
s_wg, m_wg, _ = size_to_employment(lambda pp, s: compensated_wage_gap(pp, base_m, s), p, cfg, E70, 1.0, 1.6)
ex[f"comp. wage gap x{s_wg:.2f}"] = m_wg
s_k, m_k, _ = size_to_employment(lambda pp, s: cost_scaled(pp, s), p, cfg, E70, 0.1, 1.0)
ex[f"cost x{s_k:.2f}"] = m_k
md.append(table({"baseline": base_m, **ex}, "Single-factor experiments sized to the 1970s employment rate (0.73)")); print(md[-1])

# ---- cohort accounting: tau_w and gam_e from the slides (p.35), cost residual to match employment
coh = {1940: (0.71, 0.50, 0.62), 1950: (0.74, 0.55, 0.67), 1960: (0.77, 0.58, 0.71), 1970: (0.76, 0.68, 0.73), 1980: (0.77, 0.69, 0.72)}
crow = {}; cscale = {}
for c, (gap, ge, E) in coh.items():
    g = gap / coh[1940][0]
    pc = compensated_wage_gap(p, base_m, g).replace(gam_e=p.gam_e * ge / coh[1940][1])
    if c == 1940:
        crow["1940"] = base_m; cscale[c] = 1.0; continue
    s_c, m_c, _ = size_to_employment(lambda pp, s: cost_scaled(pp, s), pc, cfg, E, 0.2, 2.0)
    crow[f"{c} (cost x{s_c:.2f})"] = m_c; cscale[c] = s_c
md.append(table(crow, "Cohort accounting (tau_w and gamma_e from the data, cost residual)")); print(md[-1])

# ---- mechanism counterfactuals
cf = counterfactuals(p, cfg)
md.append(table(cf, "Mechanism counterfactuals (baseline parameters)")); print(md[-1])
md.append(f"\nElapsed {time.time()-t0:.0f}s. Calibration file: {a.calib}\n")
out = a.out or os.path.join(HERE, "..", "output", f"final_results_{a.tag or ('full' if a.full else 'coarse')}.md")
open(out, "w").write("\n".join(md)); print("saved", out)
fl = lambda d: {k: float(v) for k, v in d.items() if isinstance(v, (int, float, np.floating))}
json.dump(dict(calib=a.calib, full=a.full, baseline=fl(base_m),
               scales={"returns": s_roe, "wage_gap": s_wg, "cost": s_k}, experiments={k: fl(v) for k, v in ex.items()},
               cohort_cost_scale={str(k): v for k, v in cscale.items()}, cohorts={k: fl(v) for k, v in crow.items()},
               counterfactuals={k: fl(v) for k, v in cf.items()}),
          open(out[:-3] + ".json", "w"), indent=1)
