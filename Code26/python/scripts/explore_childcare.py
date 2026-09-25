"""Model-ingredient exploration: child care as higher home productivity at ages 25-39
(`home_young_mult`) versus the utility-cost multiplier `km_max`, around a calibration point.

usage: python scripts/explore_childcare.py output/x0_km14.json
"""
import sys, os, json, time, warnings, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import evaluate

x = json.load(open(sys.argv[1]))["x"]
cfg = SimConfigFinal(N=60, n_cohorts=90)
print(f"{'hm':>5s} {'km':>5s} {'obj':>6s} {'E/pop':>6s} {'hrs':>6s} {'LC':>6s} {'PT':>6s} {'Car':>6s} {'NiLF':>6s} {'quit':>6s} {'qrec':>6s} {'dE':>6s} {'gap':>6s}")
rows = []
for hm, km in itertools.product([1.0, 1.5, 2.0, 3.0], [1.0, 4.0, 8.0]):
    base = FinalParams(n_omega=3, n_kbar=3, n_km=3 if km > 1 else 1, home_young_mult=hm)
    xx = dict(x); xx["km_max"] = km
    t0 = time.time()
    obj, m, parts, p = evaluate(base, xx, cfg)
    rows.append(dict(home_young_mult=hm, km_max=km, obj=obj, m={k: float(v) for k, v in m.items() if isinstance(v, (int, float))}))
    print(f"{hm:5.2f} {km:5.1f} {obj:6.2f} {m['E/pop']:6.3f} {m['hours|E']:6.3f} {m['share Lifecycle']:6.3f} {m['share PT']:6.3f} "
          f"{m['share Career']:6.3f} {m['share NiLF']:6.3f} {m['quit/m exp']:6.4f} {m['quit/m rec']:6.4f} {m['dE/pop rec-exp (pts)']:6.2f} "
          f"{m['wage gap (hourly ratio)']:6.3f}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "output", "explore_childcare.json"), "w"), indent=1)
