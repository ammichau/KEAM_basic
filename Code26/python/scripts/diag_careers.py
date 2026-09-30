"""Diagnostics for the career taxonomy at the calibrated point: distribution of average annual
hours around the 400/1500-hour cutoffs, classification by (omega, kbar) type cell, and the
simulation-seed noise of the shares. Output: output/diag_careers.json/.md.

usage: python scripts/diag_careers.py [--calib output/final_calib_full.json]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params
from keam.final.solve import solve_all
from keam.final.simulate import simulate_final
from keam.final.moments import moments_final, HOURS_PER_YEAR
from keam.final.params import make_types

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_full.json")
ap.add_argument("--seeds", type=int, default=4)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
x = json.load(open(os.path.join(ROOT, a.calib)))["x"]
p = apply_params(FinalParams(), x)
t0 = time.time(); sol = solve_all(p); t_solve = time.time() - t0
cfg = SimConfigFinal(N=60, n_cohorts=90)
sim = simulate_final(p, sol, cfg); m = moments_final(sim)
lo, hi = sim.cfg.window or (cfg.L, sim.T)
m0, m1, _ = p.age_months
full = (sim.entry + m0 + m1 <= hi) & (sim.entry >= lo - (m0 + m1))
hrs_y = HOURS_PER_YEAR * sim.hours[full, :m0].mean(1)
hrs_m = HOURS_PER_YEAR * sim.hours[full, m0:m0 + m1].mean(1)
hrs_all = HOURS_PER_YEAR * sim.hours[full, :m0 + m1].mean(1)
lifecycle = (hrs_m >= 1500) & (hrs_y < 600)
cls = np.where(lifecycle, 0, np.where(hrs_all >= 1500, 2, np.where(hrs_all >= 400, 1, 3)))  # LC, PT, Career, NiLF
names = ["Lifecycle", "PT", "Career", "NiLF"]
edges = [0, 50, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2400, 4000]
hist, _ = np.histogram(hrs_all, bins=edges)
out = {"moments": {k: float(v) for k, v in m.items()}, "solve_seconds": t_solve, "n_full": int(full.sum()),
       "hours_hist": {f"[{edges[i]},{edges[i+1]})": int(hist[i]) for i in range(len(hist))}}
# classification by (omega, kbar) cell, averaged over km
O, K, M = make_types(p)
kt = sim.ktype[full]
cell = {}
for io, om in enumerate(np.unique(O)):
    for ik, kb in enumerate(np.unique(K)):
        sel = np.isclose(O[kt], om) & np.isclose(K[kt], kb)
        if sel.sum() == 0:
            continue
        cell[f"omega={om:.3f},kbar={kb:.4f}"] = {"n": int(sel.sum()), "mean hours": float(hrs_all[sel].mean()),
                                                 **{nm: float((cls[sel] == j).mean()) for j, nm in enumerate(names)}}
out["by_cell"] = cell
# seed noise
seeds = []
for s in range(a.seeds):
    ms = moments_final(simulate_final(p, sol, SimConfigFinal(N=60, n_cohorts=90, seed=1000 + s)))
    seeds.append({k: float(ms[k]) for k in ["share Lifecycle", "share PT", "share Career", "share NiLF", "E/pop",
                                            "quit/m exp", "quit/m rec", "dE/pop rec-exp (pts)"]})
out["seed_runs"] = seeds
json.dump(out, open(os.path.join(ROOT, "output", "diag_careers.json"), "w"), indent=1)
L = ["# Career taxonomy diagnostics at the calibrated point", "",
     f"Solve {t_solve:.0f} s; {int(full.sum())} full careers.", "",
     "| avg annual hours (ages 25-54) | women |", "|---|---|"]
L += [f"| {k} | {v} |" for k, v in out["hours_hist"].items()]
L += ["", "| type cell | n | mean hours | Lifecycle | PT | Career | NiLF |", "|---|---|---|---|---|---|---|"]
L += [f"| {k} | {v['n']} | {v['mean hours']:.0f} | {v['Lifecycle']:.2f} | {v['PT']:.2f} | {v['Career']:.2f} | {v['NiLF']:.2f} |"
      for k, v in cell.items()]
L += ["", "| seed | LC | PT | Career | NiLF | E/pop | quit exp | quit rec | dE |", "|---|---|---|---|---|---|---|---|---|"]
L += [f"| {1000+i} | {r['share Lifecycle']:.3f} | {r['share PT']:.3f} | {r['share Career']:.3f} | {r['share NiLF']:.3f} | "
      f"{r['E/pop']:.3f} | {r['quit/m exp']:.4f} | {r['quit/m rec']:.4f} | {r['dE/pop rec-exp (pts)']:.2f} |"
      for i, r in enumerate(seeds)]
open(os.path.join(ROOT, "output", "diag_careers.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
