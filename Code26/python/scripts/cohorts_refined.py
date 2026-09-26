"""Cohort accounting with two unknowns per cohort: the cost-of-work scale and the wage-penalty
scale (tau_w, with the husband's income compensated) are solved jointly so that the cohort's
employment rate AND its measured within-couple wage gap match the data (slides p.35), given the
cohort's returns to experience gamma_e.

usage: python scripts/cohorts_refined.py --calib output/final_calib_full.json [--coarse]
Writes output/cohorts_refined_<tag>.md / .json
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import root
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params, params_from_calib
from keam.final.experiments import run, compensated_wage_gap, cost_scaled

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--coarse", action="store_true")
ap.add_argument("--n-jobs", type=int, default=0)
ap.add_argument("--max-eval", type=int, default=22)
ap.add_argument("--tag", default="", help="output name suffix (default: full or coarse)")
a = ap.parse_args()
if a.n_jobs:
    os.environ["KEAM_NJOBS"] = str(a.n_jobs)
HERE = os.path.dirname(os.path.abspath(__file__))
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p = params_from_calib(json.load(open(a.calib)), base)
cfg = SimConfigFinal(N=60, n_cohorts=90)
COH = {1940: (0.71, 0.50, 0.62), 1950: (0.74, 0.55, 0.67), 1960: (0.77, 0.58, 0.71), 1970: (0.76, 0.68, 0.73), 1980: (0.77, 0.69, 0.72)}
t0 = time.time()
base_m, _, _ = run(p, cfg)
gap0 = base_m["wage gap (hourly ratio)"]         # model's 1940 gap; cohort targets are scaled by data ratios
out = {"1940": dict(m={k: float(v) for k, v in base_m.items() if isinstance(v, (int, float, np.floating))},
                    cost_scale=1.0, tau_scale=1.0, tau_w=p.tau_w, n_eval=1)}
print(f"1940 baseline: E {base_m['E/pop']:.3f} gap {gap0:.3f}  [{time.time()-t0:.0f}s]", flush=True)

class Budget(Exception):
    pass

for c, (gap, ge, E) in COH.items():
    if c == 1940:
        continue
    gap_target = gap0 * gap / COH[1940][0]        # keep the model's 1940 level, apply the data's ratio
    pc = p.replace(gam_e=p.gam_e * ge / COH[1940][1])
    cache = {}; best = [None]
    def residual(x):
        s_cost, g_tau = float(x[0]), float(x[1])
        key = (round(s_cost, 5), round(g_tau, 5))
        if key not in cache:
            if len(cache) >= a.max_eval:
                raise Budget()
            pp = cost_scaled(compensated_wage_gap(pc, base_m, g_tau), s_cost)
            m, _, _ = run(pp, cfg)
            cache[key] = m
            r = np.array([m["E/pop"] - E, m["wage gap (hourly ratio)"] - gap_target])
            if best[0] is None or np.abs(r).sum() < best[0][0]:
                best[0] = (np.abs(r).sum(), s_cost, g_tau, m)
            print(f"  {c}: cost x{s_cost:.3f} tau x{g_tau:.3f} -> E {m['E/pop']:.4f} (target {E}) gap "
                  f"{m['wage gap (hourly ratio)']:.4f} (target {gap_target:.3f})  [{time.time()-t0:.0f}s]", flush=True)
        m = cache[key]
        return np.array([m["E/pop"] - E, m["wage gap (hourly ratio)"] - gap_target])
    x0 = np.array([1.5, gap / COH[1940][0]])
    try:
        sol = root(residual, x0, method="hybr", options=dict(xtol=2e-3, eps=0.03))
    except Budget:
        pass
    err, s_cost, g_tau, m = best[0]
    out[str(c)] = dict(m={k: float(v) for k, v in m.items() if isinstance(v, (int, float, np.floating))},
                       cost_scale=s_cost, tau_scale=g_tau, tau_w=p.tau_w * g_tau, gap_target=gap_target,
                       resid=err, n_eval=len(cache))
    print(f"{c}: cost x{s_cost:.3f}, tau_w {p.tau_w * g_tau:.3f}, |resid| {err:.4f} after {len(cache)} evaluations", flush=True)

KEYS = ["E/pop", "hours|E", "quit/m exp", "quit/m rec", "E->nonE/m exp", "E->nonE/m rec", "dE/pop rec-exp (pts)",
        "wife share exp", "wife share rec", "wage gap (hourly ratio)", "share Lifecycle", "share PT", "share Career",
        "share NiLF", "cons drop at H job loss exp (%)", "cons drop at H job loss rec (%)"]
names = list(out)
md = ["### Cohort accounting, refined: cost scale and tau_w solved jointly for employment and the wage gap", "",
      "| cohort | " + " | ".join(names) + " |", "|---|" + "---|" * len(names),
      "| cost scale | " + " | ".join(f"{out[n]['cost_scale']:.3f}" for n in names) + " |",
      "| tau_w | " + " | ".join(f"{out[n]['tau_w']:.3f}" for n in names) + " |",
      "| gamma_e | " + " | ".join(f"{p.gam_e * COH[int(n)][1] / COH[1940][1]:.3f}" for n in names) + " |"]
for k in KEYS:
    md.append(f"| {k} | " + " | ".join(f"{out[n]['m'].get(k, np.nan):.4f}" for n in names) + " |")
md.append(f"\nElapsed {time.time()-t0:.0f}s.\n")
tag = a.tag or ("coarse" if a.coarse else "full")
open(os.path.join(HERE, "..", "output", f"cohorts_refined_{tag}.md"), "w").write("\n".join(md))
json.dump(out, open(os.path.join(HERE, "..", "output", f"cohorts_refined_{tag}.json"), "w"), indent=1)
print("\n".join(md))
