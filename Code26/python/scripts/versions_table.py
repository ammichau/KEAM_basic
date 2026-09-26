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
                                      "version 4c=v4c,version 5 (log utility)=v5")
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
open(os.path.join(ROOT, a.out + ".md"), "w").write(note + "\n\n" + "\n".join(L) + "\n")
print(note); print("\n".join(L))
