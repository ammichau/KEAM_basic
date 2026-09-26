"""Assemble RESULTS.md (repository root) from the output files:
  output/final_calib_full.json   (calibration; falls back to the latest coarse round)
  output/final_results_full.json (run_final.py)
  output/robustness_final.json   (robustness_final.py, optional)

usage: python scripts/write_results.py [--calib f] [--results f] [--robust f] [--out ../../RESULTS.md]
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from keam.final.calibrate import TARGETS

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_full.json")
ap.add_argument("--results", default="output/final_results_full.json")
ap.add_argument("--robust", default="output/robustness_final.json")
ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "RESULTS.md"))
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); PY = os.path.join(HERE, "..")
def load(rel):
    p = rel if os.path.isabs(rel) else os.path.join(PY, rel)
    return json.load(open(p)) if os.path.exists(p) else None
calib = load(a.calib); res = load(a.results); rob = load(a.robust)
L = []
L.append("# Final model results: 1940s cohort calibration, trend experiments, mechanism\n")
L.append("All numbers are produced by scripts in `Code26/python/scripts`; the files cited are in "
         "`Code26/python/output`. Model specification: `FINAL_MODEL.md`.\n")
# ---- calibration
L.append("## 1. Calibration of the 1940s cohort\n")
if calib:
    m = calib.get("moments", {})
    L.append(f"Source: `{a.calib}` (objective {calib['obj']:.3f}, {calib.get('n_eval', '?')} evaluations, "
             f"{'100' if not calib.get('coarse', True) else '27'} types).\n")
    L.append("| parameter | value |\n|---|---|")
    for k, v in calib["x"].items():
        L.append(f"| {k} | {v:.4f} |")
    L.append("\n| target | data | model | deviation |\n|---|---|---|---|")
    for k, tv in TARGETS.items():
        mv = m.get(k, np.nan); sc = 1.0 if "pts" in k else tv
        L.append(f"| {k} | {tv:.4f} | {mv:.4f} | {100 * (mv - tv) / sc:+.1f}% |")
    L.append("\n| untargeted moment | model |\n|---|---|")
    for k in ["U rate", "wife share exp", "wife share rec", "HH income rec/exp - 1 (%)", "cons drop at H job loss exp (%)",
              "cons drop at H job loss rec (%)", "mean assets/monthly HH inc", "share e at cap"]:
        if k in m:
            L.append(f"| {k} | {m[k]:.4f} |")
    L.append("")
KEYS = ["E/pop", "hours|E", "U rate", "quit/m exp", "quit/m rec", "E->nonE/m exp", "E->nonE/m rec",
        "dE/pop rec-exp (pts)", "wife share exp", "wife share rec", "wage gap (hourly ratio)",
        "share Lifecycle", "share PT", "share Career", "share NiLF", "HH income rec/exp - 1 (%)",
        "cons drop at H job loss exp (%)", "cons drop at H job loss rec (%)", "mean assets/monthly HH inc"]
def table(rows, keys=KEYS):
    out = ["| moment | " + " | ".join(rows) + " |", "|---|" + "---|" * len(rows)]
    for k in keys:
        out.append(f"| {k} | " + " | ".join(f"{r.get(k, np.nan):.4f}" for r in rows.values()) + " |")
    return "\n".join(out) + "\n"
def gap(mm):
    return 100 * (mm["quit/m rec"] - mm["quit/m exp"])
if res:
    b = res["baseline"]
    L.append("## 2. Single-factor experiments sized to the 1970s employment rate\n")
    L.append(f"Source: `{a.results}`. Scales: returns to experience x{res['scales']['returns']:.3f}, compensated "
             f"wage gap x{res['scales']['wage_gap']:.3f} (husband income scaled to keep household income constant "
             f"at baseline behaviour), cost of work x{res['scales']['cost']:.3f}.\n")
    L.append(table({"baseline": b, **res["experiments"]}))
    L.append("Change in the cyclical quit gap (recession minus expansion monthly quit rate, percentage points) and in "
             "the recession employment drop relative to the baseline:\n")
    L.append("| experiment | quit gap | ΔE/pop rec-exp (pts) | E/pop |\n|---|---|---|---|")
    L.append(f"| baseline | {gap(b):+.3f} | {b['dE/pop rec-exp (pts)']:+.3f} | {b['E/pop']:.3f} |")
    for k, v in res["experiments"].items():
        L.append(f"| {k} | {gap(v):+.3f} | {v['dE/pop rec-exp (pts)']:+.3f} | {v['E/pop']:.3f} |")
    L.append("")
    L.append("## 3. Cohort accounting\n")
    L.append("τ_w and γ_e follow the slides (p.35) relative to 1940 (wage gap 0.71, 0.74, 0.77, 0.76, 0.77; γ_e 0.50, 0.55, "
             "0.58, 0.68, 0.69), with the husband's income compensated; the cost of work is scaled to reproduce each "
             "cohort's employment rate (0.62, 0.67, 0.71, 0.73, 0.72).\n")
    L.append(table(res["cohorts"]))
    L.append("## 4. Mechanism counterfactuals (baseline parameters)\n")
    cf = res["counterfactuals"]
    L.append("| counterfactual | quit exp | quit rec | quit gap (pts) | ΔE/pop rec-exp (pts) | E/pop |\n|---|---|---|---|---|---|")
    for k, v in cf.items():
        L.append(f"| {k} | {v['quit/m exp']:.4f} | {v['quit/m rec']:.4f} | {gap(v):+.3f} | {v['dE/pop rec-exp (pts)']:+.3f} | {v['E/pop']:.3f} |")
    L.append("")
    base_gap = gap(cf["baseline"])
    L.append("Decomposition of the baseline quit gap (share removed when each channel is switched off):\n")
    L.append("| channel | quit gap without it | share of baseline gap |\n|---|---|---|")
    for k in ["acyclical husband risk", "acyclical job finding", "no recession wage cut", "acyclical own job loss"]:
        if k in cf:
            g = gap(cf[k]); L.append(f"| {k} | {g:+.3f} | {100 * (base_gap - g) / base_gap if base_gap else np.nan:+.0f}% |")
    L.append("")
if rob:
    L.append("## 5. Robustness (calibrated parameters held fixed)\n")
    L.append(f"Source: `{a.robust}`.\n")
    names = list(rob)
    L.append("| variant | E/pop | quit gap | ΔE/pop rec-exp | acyclical husband risk: quit gap | RoE experiment: E/pop | RoE: quit gap |\n|---|---|---|---|---|---|---|")
    for n in names:
        r = rob[n]
        L.append(f"| {n} | {r['m']['E/pop']:.3f} | {gap(r['m']):+.3f} | {r['m']['dE/pop rec-exp (pts)']:+.3f} | "
                 f"{r['acyc']['quit_gap']:+.3f} | {r['roe']['E']:.3f} | {r['roe']['quit_gap']:+.3f} |")
    L.append("")
L.append("## 6. What is fragile\n")
L.append("* The never-working (NiLF) share is the least well fitted target; it depends on the home-production "
         "curvature in productivity (α_h) and the hours scaling of the fixed cost.\n"
         "* The experience cap e_max is calibrated; the wage gap among employed wives is largely the experience "
         "premium at the cap, so the returns-to-experience experiment interacts with it.\n"
         "* The transitory cost shock (sd σ_κ) drives the monthly quit rate; its distribution is not disciplined by "
         "micro data beyond the quit and exit rates.\n"
         "* Career shares are computed on annual hours over ages 25-54 from the model's 4,000-hour endowment; the "
         "data taxonomy uses reported annual hours.\n")
open(a.out, "w").write("\n".join(L))
print("wrote", a.out)
