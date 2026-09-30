"""Scan the life-cycle cost multiplier and the permanent cost level around a calibration point
to locate the region where life-cycle women (out at 25-39, full-time at 40-54) appear.

usage: python scripts/explore_lifecycle.py output/final_calib_coarse.json
"""
import sys, os, json, time, warnings, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params, evaluate, TARGETS

calib = json.load(open(sys.argv[1]))
x = dict(calib["x"])
base = FinalParams(n_omega=3, n_kbar=3, n_km=3)
cfg = SimConfigFinal(N=60, n_cohorts=90)
grid_km = [3.0, 4.5, 6.0, 8.0]
grid_kb = [0.06, 0.12, 0.24]
print(f"{'km_max':>7s} {'kb_max':>7s} {'obj':>6s} {'E/pop':>6s} {'LC':>6s} {'PT':>6s} {'Car':>6s} {'NiLF':>6s} {'quit':>6s} {'wsh':>6s}")
rows = []
for km, kb in itertools.product(grid_km, grid_kb):
    xx = dict(x); xx["km_max"] = km; xx["kbar_max"] = kb
    t0 = time.time()
    obj, m, parts, _ = evaluate(base, xx, cfg)
    rows.append(dict(km_max=km, kbar_max=kb, obj=obj, m=m))
    print(f"{km:7.2f} {kb:7.3f} {obj:6.2f} {m['E/pop']:6.3f} {m['share Lifecycle']:6.3f} {m['share PT']:6.3f} "
          f"{m['share Career']:6.3f} {m['share NiLF']:6.3f} {m['quit/m exp']:6.4f} {m['wife share exp']:6.3f}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "output", "explore_lifecycle.json"), "w"), indent=1)
