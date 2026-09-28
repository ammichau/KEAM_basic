"""Distribution of liquid assets in the simulated baseline of a calibration: mean, quantiles and the share of
households below one, three and six months of household income, overall and by the wife's age group.

usage: python scripts/asset_distribution.py --calib output/final_calib_v7c_full.json [--tag v7c]
Writes output/asset_distribution_<tag>.json/.md. Assets are in the model's income units; the ratios divide by
the household's own monthly income (wife plus husband) in the same month.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from keam.final import calibrate as C
from keam.final.solve import solve_all
from keam.final.simulate import simulate_final, SimConfigFinal

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--tag", default="")
ap.add_argument("--coarse", action="store_true", help="27-type grid")
ap.add_argument("--fixed", action="append", default=[], help="override a FinalParams field: name=value (repeatable)")
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
calib = json.load(open(os.path.join(ROOT, a.calib)))
tag = a.tag or os.path.basename(a.calib).replace("final_calib_", "").replace("_full.json", "")
from keam.final.params import FinalParams
base = FinalParams(n_omega=3, n_kbar=3, n_km=3) if a.coarse else FinalParams()
p = C.params_from_calib(calib, base)
for kv in a.fixed:
    n, v = kv.split("="); p = p.replace(**{n: type(getattr(p, n))(float(v))})
t0 = time.time()
sol = solve_all(p, n_jobs=int(os.environ.get("KEAM_NJOBS", "0")) or None)
cfg = SimConfigFinal(N=60, n_cohorts=90)
sim = simulate_final(p, sol, cfg)
L = cfg.L; lo, hi = sim.cfg.window or (L, sim.T)
inwin = (sim.calendar >= lo) & (sim.calendar < hi)
inc = sim.inc_w + sim.inc_h
ratio = np.where(inc > 1e-9, sim.a / np.maximum(inc, 1e-9), np.nan)
age = np.zeros_like(sim.a, dtype=int)                      # 0: 25-39, 1: 40-54, 2: 55-64
months = np.arange(sim.a.shape[1])
b1, b2 = p.age_months[0], p.age_months[0] + p.age_months[1]
age[:, months >= b1] = 1; age[:, months >= b2] = 2


def stats(m):
    r = ratio[m]; r = r[np.isfinite(r)]
    return dict(n=int(r.size), mean=float(r.mean()), p10=float(np.percentile(r, 10)), p25=float(np.percentile(r, 25)),
                median=float(np.median(r)), p75=float(np.percentile(r, 75)), p90=float(np.percentile(r, 90)),
                share_below_1m=float((r < 1).mean()), share_below_3m=float((r < 3).mean()), share_below_6m=float((r < 6).mean()),
                share_at_zero=float((sim.a[m] < 1e-6).mean()), share_near_amax=float((sim.a[m] > 0.95 * p.a_max).mean()))


out = dict(calib=a.calib, seconds=time.time() - t0, all=stats(inwin), by_age={f"{n}": stats(inwin & (age == k))
           for k, n in enumerate(["25-39", "40-54", "55-64"])},
           by_husband_state={n: stats(inwin & (sim.hstat == k)) for k, n in enumerate(["E", "R", "U"])},
           wife_employed=stats(inwin & (sim.emp == 1)), wife_nonemployed=stats(inwin & (sim.emp == 0)),
           mean_assets_over_mean_income=float(sim.a[inwin].mean() / inc[inwin].mean()))
json.dump(out, open(os.path.join(ROOT, "output", f"asset_distribution_{tag}.json"), "w"), indent=1)
rows = [("all", out["all"])] + [(f"age {k}", v) for k, v in out["by_age"].items()] + \
       [(f"husband {k}", v) for k, v in out["by_husband_state"].items()] + \
       [("wife employed", out["wife_employed"]), ("wife non-employed", out["wife_nonemployed"])]
Lm = [f"# Liquid assets in months of the household's own income (`{a.calib}`{' coarse grid' if a.coarse else ''}, fixed {a.fixed}, a_max {p.a_max})", "",
      f"Mean assets over mean income: {out['mean_assets_over_mean_income']:.2f} months.", "",
      "| group | mean | p10 | p25 | median | p75 | p90 | share < 1 month | share < 3 months | share < 6 months | share at zero | share near a_max |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for n, s in rows:
    Lm.append(f"| {n} | {s['mean']:.2f} | {s['p10']:.2f} | {s['p25']:.2f} | {s['median']:.2f} | {s['p75']:.2f} | {s['p90']:.2f} | "
              f"{s['share_below_1m']:.0%} | {s['share_below_3m']:.0%} | {s['share_below_6m']:.0%} | {s['share_at_zero']:.0%} | {s['share_near_amax']:.1%} |")
open(os.path.join(ROOT, "output", f"asset_distribution_{tag}.md"), "w").write("\n".join(Lm) + "\n")
print("\n".join(Lm))
