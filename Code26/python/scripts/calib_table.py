"""Render a calibration JSON (from calibrate_final.py / calibrate_childcare.py) as a markdown table.

usage: python scripts/calib_table.py output/final_calib_childcare2.json [> file.md]
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from keam.final.calibrate import TARGETS
d = json.load(open(sys.argv[1]))
m = d.get("moments", {})
print(f"Calibration `{sys.argv[1]}`: objective {d['obj']:.3f} after {d.get('n_eval', '?')} evaluations\n")
print("| parameter | value |\n|---|---|")
for k, v in d["x"].items():
    print(f"| {k} | {v:.4f} |")
print("\n| target | data | model | deviation |\n|---|---|---|---|")
for k, tv in TARGETS.items():
    mv = m.get(k, float("nan")); sc = 1.0 if "pts" in k else tv
    print(f"| {k} | {tv:.4f} | {mv:.4f} | {100 * (mv - tv) / sc:+.1f}% |")
print("\n| untargeted moment | model |\n|---|---|")
for k in ["U rate", "wife share exp", "wife share rec", "HH income rec/exp - 1 (%)", "cons drop at H job loss exp (%)",
          "cons drop at H job loss rec (%)", "mean assets/monthly HH inc", "share e at cap", "wage gap (FTE earnings ratio)"]:
    if k in m:
        print(f"| {k} | {m[k]:.4f} |")
