"""Sensitivity of the precaution / hoarding split to the cyclical parameters: for each parameter, the
baseline and the three counterfactual solves (husband's risk acyclical, wife's job finding acyclical,
both) are repeated at a +step perturbation, and the change in the precaution share, the hoarding share,
the quit gap, the UE-rate cyclicality and the recession employment drop is reported per percent change
of the parameter. Output: output/jacobian_channels_<tag>.json/.md (resumable).

usage: python scripts/jacobian_channels.py --calib output/final_calib_v4_full.json --tag v4 [--step 0.10]
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run, acyclical_husband, acyclical_finding

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--tag", default="v4")
ap.add_argument("--step", type=float, default=0.10)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
calib = json.load(open(os.path.join(ROOT, a.calib)))
p0 = C.params_from_calib(calib)
cfg = SimConfigFinal(N=60, n_cohorts=90)
r0 = p0.lam_f[1] / p0.lam_f[0]
# parameter -> function giving the perturbed FinalParams (multiplicative step on the named quantity)
PERT = {
    "job-finding fall in recessions (1 - lam_f ratio)": lambda p, s: p.replace(lam_f=(p.lam_f[0], p.lam_f[0] * (1 - (1 - r0) * (1 + s)))),
    "UI cut in recessions (1 - ui_rec_mult)": lambda p, s: p.replace(ui_rec_mult=1 - (1 - p.ui_rec_mult) * (1 + s)),
    "husband job-loss rate in recessions": lambda p, s: p.replace(lamH_loss=(p.lamH_loss[0], p.lamH_loss[1] * (1 + s))),
    "husband job-finding rate in recessions": lambda p, s: p.replace(lamH_find=(p.lamH_find[0], p.lamH_find[1] * (1 + s))),
    "husband job-loss rate (both states)": lambda p, s: p.replace(lamH_loss=(p.lamH_loss[0] * (1 + s), p.lamH_loss[1] * (1 + s))),
    "UI replacement (both states)": lambda p, s: p.replace(ym_share=(p.ym_share[0], p.ym_share[1], p.ym_share[2] * (1 + s))),
    "wife own job loss in recessions (lam_u1)": lambda p, s: p.replace(lam_u=(p.lam_u[0], p.lam_u[1] * (1 + s))),
    "job-finding efficiency level (lam_f0)": lambda p, s: p.replace(lam_f=(p.lam_f[0] * (1 + s), p.lam_f[1] * (1 + s))),
    "cost-shock sd (sd_kT)": lambda p, s: p.replace(sd_kT=p.sd_kT * (1 + s)),
    "recession wage cut (1 - phi_rec)": lambda p, s: p.replace(phi_rec=1 - (1 - p.phi_rec) * (1 + s), phi_rec_H=1 - (1 - p.phi_rec_H) * (1 + s)),
    "recession persistence (piz[1,1])": lambda p, s: p.replace(piz=np.array([[p.piz[0, 0], p.piz[0, 1]], [1 - p.piz[1, 1] * (1 + s), p.piz[1, 1] * (1 + s)]])),
    "asset limit a_max": lambda p, s: p.replace(a_max=p.a_max * (1 + s)),
    "risk aversion gamma": lambda p, s: p.replace(gamma=p.gamma * (1 + s)),
}
out_json = os.path.join(ROOT, "output", f"jacobian_channels_{a.tag}.json")
res = json.load(open(out_json)) if os.path.exists(out_json) else {"calib": a.calib, "step": a.step, "points": {}}


def gap(m):
    return 100 * (m["quit/m rec"] - m["quit/m exp"])


def decompose(p):
    mb = run(p, cfg)[0]; mh = run(acyclical_husband(p), cfg)[0]
    mf = run(acyclical_finding(p), cfg)[0]; mboth = run(acyclical_finding(acyclical_husband(p)), cfg)[0]
    g, gh, gf, gb = gap(mb), gap(mh), gap(mf), gap(mboth)
    return dict(gap=g, precaution=(g - gh) / g, hoarding=(g - gf) / g, both_off=(g - gb) / g,
                sdUE=float(mb["sd log UE (women)"]), dE=float(mb["dE/pop rec-exp (pts)"]), dE_acycH=float(mh["dE/pop rec-exp (pts)"]),
                quit_rec=float(mb["quit/m rec"]), E=float(mb["E/pop"]))


t0 = time.time()
if "base" not in res["points"]:
    res["points"]["base"] = decompose(p0); json.dump(res, open(out_json, "w"), indent=1)
    print("base", {k: round(v, 3) for k, v in res["points"]["base"].items()}, f"[{time.time() - t0:.0f}s]", flush=True)
b = res["points"]["base"]
for name, f in PERT.items():
    if name in res["points"]:
        continue
    res["points"][name] = decompose(f(p0, a.step)); json.dump(res, open(out_json, "w"), indent=1)
    r = res["points"][name]
    print(f"{name}: precaution {r['precaution']:.3f} hoarding {r['hoarding']:.3f} sdUE {r['sdUE']:.4f} dE {r['dE']:+.2f} [{time.time() - t0:.0f}s]", flush=True)
keys = ["precaution", "hoarding", "gap", "sdUE", "dE", "dE_acycH", "quit_rec"]
L = [f"# Sensitivity of the precaution / hoarding split (`{a.calib}`, +{100*a.step:.0f}% steps, parameters otherwise fixed)", "",
     f"Baseline: precaution {b['precaution']:.0%}, hoarding {b['hoarding']:.0%}, both off {b['both_off']:.0%}, quit gap {b['gap']:+.2f} points, "
     f"sd log UE {b['sdUE']:.4f}, recession employment drop {b['dE']:+.2f} ({b['dE_acycH']:+.2f} without cyclical husband risk).", "",
     "Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in "
     "percentage points of the rate; sd log UE in units).", "",
     "| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |",
     "|---|---|---|---|---|---|---|---|"]
for name, r in res["points"].items():
    if name == "base":
        continue
    d = lambda k, scale: (r[k] - b[k]) * scale / (100 * a.step)
    L.append(f"| {name} | {d('precaution', 100):+.2f} | {d('hoarding', 100):+.2f} | {d('gap', 1):+.3f} | {d('sdUE', 1):+.4f} | "
             f"{d('dE', 1):+.3f} | {d('dE_acycH', 1):+.3f} | {d('quit_rec', 100):+.3f} |")
open(os.path.join(ROOT, "output", f"jacobian_channels_{a.tag}.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
