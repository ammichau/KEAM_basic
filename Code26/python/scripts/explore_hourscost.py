"""Ingredient scan: fixed cost of work scaling with hours, kappa (h/0.4)^p (p = 0 is the paper's
pure fixed cost), and a forced lower wage penalty tau_w, around a calibration point.

usage: python scripts/explore_hourscost.py output/final_calib_childcare3.json
"""
import sys, os, json, time, warnings, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C

x = json.load(open(sys.argv[1]))["x"]
for n in ["home_young_mult", "nu_h", "z_h", "alpha_h"]:
    if n not in C.OPTIONAL_PARAMS:
        C.OPTIONAL_PARAMS.append(n)
cfg = SimConfigFinal(N=60, n_cohorts=90)
print(f"{'p':>4s} {'tau_w':>5s} {'obj':>6s} {'E/pop':>6s} {'hrs':>6s} {'LC':>6s} {'PT':>6s} {'Car':>6s} {'NiLF':>6s} {'quit':>6s} {'qrec':>6s} {'dE':>6s} {'gap':>6s}")
rows = []
for pw, tw in itertools.product([0.0, 0.5, 1.0], [x["tau_w"], 0.72]):
    base = FinalParams(n_omega=3, n_kbar=3, n_km=3, kappa_h_power=pw)
    xx = dict(x); xx["tau_w"] = tw
    t0 = time.time()
    obj, m, parts, p = C.evaluate(base, xx, cfg)
    rows.append(dict(kappa_h_power=pw, tau_w=tw, obj=obj, m={k: float(v) for k, v in m.items() if isinstance(v, (int, float))}))
    print(f"{pw:4.1f} {tw:5.2f} {obj:6.2f} {m['E/pop']:6.3f} {m['hours|E']:6.3f} {m['share Lifecycle']:6.3f} {m['share PT']:6.3f} "
          f"{m['share Career']:6.3f} {m['share NiLF']:6.3f} {m['quit/m exp']:6.4f} {m['quit/m rec']:6.4f} {m['dE/pop rec-exp (pts)']:6.2f} "
          f"{m['wage gap (hourly ratio)']:6.3f}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "output", "explore_hourscost.json"), "w"), indent=1)
