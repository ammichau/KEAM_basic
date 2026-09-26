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
ap.add_argument("--extra", default="output/extra_experiments_full.json")
ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "RESULTS.md"))
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); PY = os.path.join(HERE, "..")
def load(rel):
    p = rel if os.path.isabs(rel) else os.path.join(PY, rel)
    return json.load(open(p)) if os.path.exists(p) else None
calib = load(a.calib); res = load(a.results); rob = load(a.robust); extra = load(a.extra)
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
    if extra:
        L.append("Supplementary experiments (`" + a.extra + "`): child-care cost scaled toward zero "
                 f"(x{extra['scales']['childcare']:.2f} of the excess home productivity at 25-39) and all cost "
                 f"components scaled jointly (x{extra['scales']['cost_all']:.2f}).\n")
        L.append(table({"baseline": extra["baseline"], **extra["experiments"]}))
        L.append("| experiment | quit gap | ΔE/pop rec-exp (pts) | E/pop |\n|---|---|---|---|")
        for k, v in extra["experiments"].items():
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
L.append("## 6. Summary of findings\n")
if calib and res:
    m = calib.get("moments", {}); b = res["baseline"]; cf = res["counterfactuals"]; bg = gap(b)
    ex = res["experiments"]; names = list(ex)
    L.append(f"* Calibration: employment {m.get('E/pop', np.nan):.3f} (target 0.62), hours {m.get('hours|E', np.nan):.3f} "
             f"(0.40), monthly quit rate {m.get('quit/m exp', np.nan):.4f} in expansions and {m.get('quit/m rec', np.nan):.4f} "
             f"in recessions (targets 0.034 / 0.028), recession employment drop {m.get('dE/pop rec-exp (pts)', np.nan):.2f} "
             f"points (-1.7), wage gap {m.get('wage gap (hourly ratio)', np.nan):.3f} (0.71); career shares life-cycle "
             f"{m.get('share Lifecycle', np.nan):.2f}, part-time {m.get('share PT', np.nan):.2f}, career "
             f"{m.get('share Career', np.nan):.2f}, NiLF {m.get('share NiLF', np.nan):.2f} (0.31 / 0.28 / 0.19 / 0.22). "
             f"Untargeted: unemployment rate {m.get('U rate', np.nan):.3f}, wife's income share {m.get('wife share exp', np.nan):.3f}, "
             f"consumption falls {abs(m.get('cons drop at H job loss exp (%)', np.nan)):.1f}% at the husband's job loss in "
             f"expansions and {abs(m.get('cons drop at H job loss rec (%)', np.nan)):.1f}% in recessions.")
    L.append(f"* Quits are pro-cyclical: the monthly quit rate falls by {100 * (1 - b['quit/m rec'] / b['quit/m exp']):.0f}% "
             f"in recessions ({bg:+.2f} points). Decomposition: making the husband's job-loss risk acyclical removes "
             f"{100 * (bg - gap(cf['acyclical husband risk'])) / bg:+.0f}% of the drop, making job finding acyclical removes "
             f"{100 * (bg - gap(cf['acyclical job finding'])) / bg:+.0f}%, removing the recession wage cut changes it by "
             f"{100 * (bg - gap(cf['no recession wage cut'])) / bg:+.0f}% (the wage cut works against the insurance motive), "
             f"and making the wife's own job loss acyclical changes it by {100 * (bg - gap(cf['acyclical own job loss'])) / bg:+.0f}%.")
    L.append(f"* Recession employment drop {b['dE/pop rec-exp (pts)']:+.2f} points in the baseline; "
             f"{cf['acyclical husband risk']['dE/pop rec-exp (pts)']:+.2f} without cyclical husband risk (precautionary labor "
             f"supply offsets {cf['acyclical husband risk']['dE/pop rec-exp (pts)'] - b['dE/pop rec-exp (pts)']:+.2f} points), "
             f"{cf['acyclical job finding']['dE/pop rec-exp (pts)']:+.2f} without the fall in job finding, "
             f"{cf['no recession wage cut']['dE/pop rec-exp (pts)']:+.2f} without the wage cut, "
             f"{cf['acyclical own job loss']['dE/pop rec-exp (pts)']:+.2f} without cyclical own job loss.")
    parts = []
    for n in names:
        parts.append(f"{n}: employment {ex[n]['E/pop']:.3f}, recession drop {ex[n]['dE/pop rec-exp (pts)']:+.2f} points "
                     f"(baseline {b['dE/pop rec-exp (pts)']:+.2f}), quit gap {gap(ex[n]):+.2f} (baseline {bg:+.2f}), "
                     f"career shares LC/PT/career/NiLF {ex[n]['share Lifecycle']:.2f}/{ex[n]['share PT']:.2f}/"
                     f"{ex[n]['share Career']:.2f}/{ex[n]['share NiLF']:.2f}")
    L.append("* Trend to cycle, each force sized to the 1970s employment rate: " + "; ".join(parts) + ".")
    coh = res["cohorts"]; ck = list(coh)
    L.append(f"* Cohort accounting with the data's wage-gap and returns-to-experience paths (household income compensated): "
             f"the residual cost scale is {', '.join(f'{k}: x{v:.2f}' for k, v in res['cohort_cost_scale'].items())}; "
             f"the recession employment drop goes from {coh[ck[0]]['dE/pop rec-exp (pts)']:+.2f} to {coh[ck[-1]]['dE/pop rec-exp (pts)']:+.2f} "
             f"points and the expansion quit rate from {coh[ck[0]]['quit/m exp']:.4f} to {coh[ck[-1]]['quit/m exp']:.4f}. "
             f"Caveat: tau_w is scaled by the raw data ratio, so the measured wage gap in the model rises to "
             f"{coh[ck[-1]]['wage gap (hourly ratio)']:.2f} by the last cohort (data 0.77); the next refinement is to solve "
             f"tau_w per cohort to hit the measured gap jointly with the cost residual.")
    L.append("")
L.append("## 7. What is fragile\n")
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
