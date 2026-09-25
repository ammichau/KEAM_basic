"""Is the experience process what keeps life-cycle women out after 40?

Starting from a calibration point with a large young-age cost multiplier, vary the experience
depreciation delta_e, the cap e_max and the life-cycle multiplier km_max. For each point it reports
the life-cycle share and, among women who are out at 25-39 (< 600 annual hours), the split of their
annual hours at 40-54 (< 400, 400-1,500, >= 1,500), their mean experience at 40 and the wage
loss implied by that experience relative to the initial experience.

usage: python scripts/explore_experience.py output/final_calib_coarse.json
Writes output/explore_experience.json.
"""
import sys, os, json, time, warnings, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal, solve_all, simulate_final, moments_final
from keam.final.calibrate import apply_params, objective_from_moments
from keam.final.moments import HOURS_PER_YEAR

calib = json.load(open(sys.argv[1]))
x = dict(calib["x"]); x["kbar_max"] = 0.06
base = FinalParams(n_omega=3, n_kbar=3, n_km=3)
cfg = SimConfigFinal(N=60, n_cohorts=90)
grid = [dict(km_max=km, delta_e=d, e_max=em)
        for km in (8.0, 10.0) for d in (0.005, 0.0025, 0.001) for em in (2.0,)]
grid += [dict(km_max=10.0, delta_e=0.005, e_max=3.0), dict(km_max=10.0, delta_e=0.005, e_max=2.0, sd_kT=0.1),
         dict(km_max=14.0, delta_e=0.005, e_max=2.0)]
print(f"{'km':>5s} {'dlt':>6s} {'emax':>4s} {'sdkT':>5s} {'obj':>6s} {'E/pop':>6s} {'LC':>6s} {'Car':>6s} {'NiLF':>6s} "
      f"{'quit':>6s} | out-young: {'n':>5s} {'<400':>5s} {'mid':>5s} {'>=1500':>6s} {'e@40':>5s} {'wloss%':>6s}")
rows = []
for g in grid:
    xx = dict(x); xx["km_max"] = g["km_max"]
    if "sd_kT" in g:
        xx["sd_kT"] = g["sd_kT"]
    p = apply_params(base, xx).replace(delta_e=g["delta_e"], e_max=g["e_max"])
    t0 = time.time()
    sol = solve_all(p); sim = simulate_final(p, sol, cfg); m = moments_final(sim)
    obj, _ = objective_from_moments(m)
    lo, hi = sim.cfg.window or (sim.cfg.L, sim.T)
    m0, m1, _ = p.age_months
    full = (sim.entry + m0 + m1 <= hi) & (sim.entry >= lo - (m0 + m1))
    hy = HOURS_PER_YEAR * sim.hours[full, :m0].mean(1)
    hm = HOURS_PER_YEAR * sim.hours[full, m0:m0 + m1].mean(1)
    out_y = hy < 600
    e0 = sim.e[full, 0][out_y]; e40 = sim.e[full, m0][out_y]
    wl = 100 * (1 - (1 + p.gam_e * e40 ** p.xi) / (1 + p.gam_e * e0 ** p.xi)).mean()
    r = dict(**g, obj=obj, m={k: float(v) for k, v in m.items()}, n_out_young=int(out_y.sum()),
             share_out_young=float(out_y.mean()),
             mid_lt400=float((hm[out_y] < 400).mean()), mid_400_1500=float(((hm[out_y] >= 400) & (hm[out_y] < 1500)).mean()),
             mid_ge1500=float((hm[out_y] >= 1500).mean()), e40=float(e40.mean()), wage_loss_pct=float(wl))
    rows.append(r)
    print(f"{g['km_max']:5.1f} {g['delta_e']:6.4f} {g['e_max']:4.1f} {xx['sd_kT']:5.2f} {obj:6.2f} {m['E/pop']:6.3f} "
          f"{m['share Lifecycle']:6.3f} {m['share Career']:6.3f} {m['share NiLF']:6.3f} {m['quit/m exp']:6.4f} | "
          f"{r['share_out_young']:11.3f} {r['mid_lt400']:5.2f} {r['mid_400_1500']:5.2f} {r['mid_ge1500']:6.2f} "
          f"{r['e40']:5.2f} {r['wage_loss_pct']:6.1f}  [{time.time() - t0:.0f}s]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "output", "explore_experience.json"), "w"), indent=1)
