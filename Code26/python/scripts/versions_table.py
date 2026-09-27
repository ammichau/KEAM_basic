"""All calibrated versions side by side: 13-target objective, each target's deviation, the precautionary and
hoarding shares of the recession quit drop and the recession employment drop with and without cyclical husband
risk (from scripts/channels.py). Moments are the 100-type ones (`moments_full` when a file has them).

usage: python scripts/versions_table.py [--versions "label=tag,..."] [--out output/versions_summary]
       each tag reads output/final_calib_<tag>_full.json and output/channels_<tag>.json (skipped if missing)
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from keam.final import calibrate as C

ap = argparse.ArgumentParser()
ap.add_argument("--versions", default="adopted iid=ls,UI cut=ui,7 wage types=om7,version 3=v3,version 4=v4,"
                                      "version 4c=v4c,version 5 (log utility)=v5,version 5b (log utility; 2nd polish)=v5b,"
                                      "version 6 (persistent shock)=v6")
ap.add_argument("--out", default="output/versions_summary")
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SHORT = {"E/pop": "E/pop", "hours|E": "hours", "share Lifecycle": "LC", "share PT": "PT", "share Career": "career",
         "share NiLF": "NiLF", "quit/m exp": "quit exp", "quit/m rec": "quit rec", "E->nonE/m exp": "exit exp",
         "E->nonE/m rec": "exit rec", "dE/pop rec-exp (pts)": "dE dev (pts)", "wage gap (hourly ratio)": "wage gap",
         "sd log UE (women)": "sd log UE"}
rows = []
for item in a.versions.split(","):
    label, tag = item.rsplit("=", 1)
    fc = os.path.join(ROOT, "output", f"final_calib_{tag}_full.json")
    if not os.path.exists(fc):
        continue
    d = json.load(open(fc)); m = d.get("moments_full", d["moments"])
    obj, dev = C.objective_from_moments(m)
    fch = os.path.join(ROOT, "output", f"channels_{tag}.json")
    ch = json.load(open(fch))["baseline"] if os.path.exists(fch) else None
    rows.append(dict(label=label, tag=tag, obj=float(obj), dev={k: float(v) for k, v in dev.items()},
                     fixed=d.get("fixed", {}), lam_f_ratio=d["x"].get("lam_f_ratio", 0.85),
                     precaution=ch and ch["precaution"], hoarding=ch and ch["hoarding"], both_off=ch and ch["both_off"],
                     gap=ch and ch["gap"], dE=ch["dE"] if ch else float(m["dE/pop rec-exp (pts)"]),
                     dE_acycH=ch and ch["dE_acycH"]))
json.dump(rows, open(os.path.join(ROOT, a.out + ".json"), "w"), indent=1)
pct = lambda v: "-" if v is None else f"{100 * v:.0f}%"
num = lambda v, f: "-" if v is None else format(v, f)
L = ["| version | objective | " + " | ".join(SHORT[k] for k in C.TARGETS) +
     " | λ_f rec/exp | quit gap | precaution | hoarding | both off | dE rec | dE acyc. husband |",
     "|---|---|" + "---|" * len(C.TARGETS) + "---|---|---|---|---|---|---|"]
for r in rows:
    L.append(f"| {r['label']} | {r['obj']:.3f} | " +
             " | ".join(f"{100 * r['dev'][k]:+.0f}%" if "pts" not in k else f"{r['dev'][k]:+.2f}" for k in C.TARGETS) +
             f" | {r['lam_f_ratio']:.2f} | {num(r['gap'], '+.2f')} | {pct(r['precaution'])} | {pct(r['hoarding'])} | "
             f"{pct(r['both_off'])} | {r['dE']:+.2f} | {num(r['dE_acycH'], '+.2f')} |")
note = ("Objective: weighted sum of squared deviations over the 13 targets of `keam/final/calibrate.py` (100-type "
        "moments). Target columns: relative deviation from the data, except dE dev (the recession employment drop, deviation in points). Precaution / "
        "hoarding: share of the recession fall in the monthly quit rate removed when the husband's risk / the wife's "
        "job-finding efficiency is made acyclical (`scripts/channels.py`). dE: recession minus expansion employment "
        "rate (points), baseline and with acyclical husband risk.")
# ---- the carry-forward rule of the plan (CLAUDE.md): among v4c and the log-utility version (v5b if it exists, else
# v5), the lowest objective whose precaution share is within 10 points of the hoarding share; the log-utility
# version is preferred if its objective is below 0.3 and its cyclical moments are within 15%
CYC = ["quit/m rec", "E->nonE/m rec", "dE/pop rec-exp (pts)", "sd log UE (women)"]
by = {r["tag"]: r for r in rows}
F = []
def cyc_dev(r):   # relative deviations of the cyclical moments (the employment drop relative to the 1.7-point target)
    return {k: r["dev"][k] / (abs(C.TARGETS[k]) if "pts" in k else 1.0) for k in CYC}
lu = by.get("v5b") or by.get("v5")
cands = [r for r in [by.get("v4c"), lu] if r and r["precaution"] is not None]
for r in cands:
    cd = cyc_dev(r); gapc = abs(r["precaution"] - r["hoarding"])
    F.append(f"* {r['label']}: objective {r['obj']:.3f}; precaution {100 * r['precaution']:.0f}% vs hoarding "
             f"{100 * r['hoarding']:.0f}% ({100 * gapc:.0f} points apart, rule: within 10); largest cyclical-moment "
             f"deviation {100 * max(abs(v) for v in cd.values()):.0f}% ({max(cd, key=lambda k: abs(cd[k]))}).")
ok = [r for r in cands if abs(r["precaution"] - r["hoarding"]) <= 0.10]
if lu in cands and lu["obj"] < 0.3 and max(abs(v) for v in cyc_dev(lu).values()) <= 0.15 and lu in ok:
    pick, why = lu, "the log-utility version fits acceptably and keeps the split balanced (balanced-growth preferences)"
elif ok:
    pick, why = min(ok, key=lambda r: r["obj"]), "lowest objective among the candidates that satisfy the split rule"
elif cands:
    pick = min(cands, key=lambda r: (abs(r["precaution"] - r["hoarding"]), r["obj"]))
    why = ("no candidate satisfies the 10-point rule; this is the one closest to parity (it also has the lowest "
           "objective)" if pick is min(cands, key=lambda r: r["obj"]) else "no candidate satisfies the 10-point rule; "
           "this is the one closest to parity")
else:
    pick = None
if pick:
    F.append(f"* Carried forward: **{pick['label']}** (`{pick['tag']}`): {why}.")
g1 = os.path.join(ROOT, "output", "channels_v4c_gamma1.json")
if os.path.exists(g1) and "v4c" in by and lu:
    b1 = json.load(open(g1))["baseline"]; v = by["v4c"]
    F.append(f"* Log utility and the precautionary channel: at the version-4c parameters (γ = 2) precaution is "
             f"{100 * v['precaution']:.0f}% and hoarding {100 * v['hoarding']:.0f}%; imposing γ = 1 without recalibrating "
             f"gives {100 * b1['precaution']:.0f}% / {100 * b1['hoarding']:.0f}% (`output/channels_v4c_gamma1.json`; employment "
             f"{b1['E']:.2f}, because the calibrated cost levels are in γ = 2 utility units), and the recalibrated "
             f"log-utility version gives {100 * lu['precaution']:.0f}% / {100 * lu['hoarding']:.0f}%. The fall is a "
             f"property of the preferences, not of the recalibration. Fit: objective {lu['obj']:.3f} versus "
             f"{v['obj']:.3f}; the largest deviations of the log-utility version are "
             + ", ".join(f"{k} {100 * lu['dev'][k]:+.0f}%" for k in sorted(lu['dev'], key=lambda k: -abs(lu['dev'][k]))[:4]
                         if 'pts' not in k) + ".")
json.dump(dict(rows=rows, pick=pick and pick["tag"]), open(os.path.join(ROOT, a.out + ".json"), "w"), indent=1)
open(os.path.join(ROOT, a.out + ".md"), "w").write(note + "\n\n" + "\n".join(L) + "\n\n" + "\n".join(F) + "\n")
print(note); print("\n".join(L)); print("\n".join(F))
