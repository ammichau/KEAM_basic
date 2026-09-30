"""Why is part-time too low and career/NiLF too high? Scan the home-production scale z_h and the
hours curvature nu_h (both fixed from the paper so far) around a calibration point.

usage: python scripts/explore_parttime.py output/x0_km14.json
"""
import sys, os, json, time, warnings, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final.calibrate import apply_params, evaluate

x = json.load(open(sys.argv[1]))["x"]
cfg = SimConfigFinal(N=60, n_cohorts=90)
print(f"{'z_h':>5s} {'nu_h':>5s} {'obj':>6s} {'E/pop':>6s} {'hrs':>6s} {'LC':>6s} {'PT':>6s} {'Car':>6s} {'NiLF':>6s} {'quit':>6s} {'gap':>6s}")
rows = []
for z_h, nu_h in itertools.product([0.45, 0.65, 0.9], [0.65, 0.45]):
    base = FinalParams(n_omega=3, n_kbar=3, n_km=3, z_h=z_h, nu_h=nu_h)
    t0 = time.time()
    obj, m, parts, p = evaluate(base, x, cfg)
    rows.append(dict(z_h=z_h, nu_h=nu_h, obj=obj, m={k: float(v) for k, v in m.items() if isinstance(v, (int, float))}))
    print(f"{z_h:5.2f} {nu_h:5.2f} {obj:6.2f} {m['E/pop']:6.3f} {m['hours|E']:6.3f} {m['share Lifecycle']:6.3f} {m['share PT']:6.3f} "
          f"{m['share Career']:6.3f} {m['share NiLF']:6.3f} {m['quit/m exp']:6.4f} {m['wage gap (hourly ratio)']:6.3f}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "output", "explore_parttime.json"), "w"), indent=1)
