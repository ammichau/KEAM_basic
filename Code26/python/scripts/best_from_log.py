"""Write a calibration JSON from the best evaluation in a calibration log (for a run that was
interrupted before writing its result file).

usage: python scripts/best_from_log.py --log output/final_calib_rho_coarse.log --out output/final_calib_rho_coarse.json [--fixed n_kT=3] [--coarse]
"""
import sys, os, json, argparse
ap = argparse.ArgumentParser()
ap.add_argument("--log", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--fixed", action="append", default=[]); ap.add_argument("--coarse", action="store_true")
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from keam.final.calibrate import TARGETS
h = [json.loads(l) for l in open(a.log) if l.startswith("{")]
b = min(h, key=lambda r: r["obj"])
out = dict(x=b["x"], obj=b["obj"], moments=b.get("m", {}), targets=TARGETS, n_eval=len(h), coarse=a.coarse,
           names=list(b["x"]), fixed={kv.split("=")[0]: float(kv.split("=")[1]) for kv in a.fixed},
           note=f"best of {len(h)} evaluations in {a.log} (evaluation {b['n']})")
json.dump(out, open(a.out, "w"), indent=1)
print("best evaluation", b["n"], "objective", round(b["obj"], 4), "->", a.out)
