"""Precautionary labor supply versus job hoarding: how the split of the recession quit drop responds to
the parameters that govern the two channels, at fixed calibrated parameters.

For each variant: baseline, husband's risk acyclical (precautionary channel off), wife's job finding
acyclical (hoarding off), both off. Shares: precaution = (gap - gap_acycH)/gap, hoarding =
(gap - gap_acycF)/gap, where gap = quit rate in recessions minus expansions (points); the residual
(wage cut, own job loss, interaction) is 1 - both shares + (gap - gap_both)/gap adjustments shown as
"both off". Output: output/channels_<tag>.md/.json.

usage: python scripts/channels.py --calib output/final_calib_ls_full.json --tag ls [--fixed name=value]
                                  [--variants all|quick]
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
ap.add_argument("--tag", default="ls")
ap.add_argument("--fixed", action="append", default=[])
ap.add_argument("--variants", default="all")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
calib = json.load(open(os.path.join(ROOT, a.calib)))
base = FinalParams()
for kv in a.fixed:
    n, v = kv.split("="); base = base.replace(**{n: type(getattr(base, n))(float(v))})
p0 = C.params_from_calib(calib, base)
cfg = SimConfigFinal(N=60, n_cohorts=90)
VARIANTS = {
    "baseline": {},
    "job finding falls 5% in recessions (ratio 0.95)": dict(lam_f=(p0.lam_f[0], 0.95 * p0.lam_f[0])),
    "job finding falls 30% in recessions (ratio 0.70)": dict(lam_f=(p0.lam_f[0], 0.70 * p0.lam_f[0])),
    "husband job loss x2.5 in recessions (data: x1.78)": dict(lamH_loss=(p0.lamH_loss[0], 2.5 * p0.lamH_loss[0])),
    "husband job finding 0.20 in recessions (data: 0.28)": dict(lamH_find=(p0.lamH_find[0], 0.20)),
    "UI replacement 15% in recessions (ui_rec_mult 0.5)": dict(ui_rec_mult=0.5),
    "UI replacement 15% always": dict(ym_share=(1.0, 0.85, 0.15)),
    "no assets": dict(a_max=0.01, nA=5),
    "risk aversion 3": dict(gamma=3.0),
    "longer recessions (persistence 0.95)": dict(piz=np.array([[0.985, 0.015], [0.05, 0.95]])),
}
if a.variants == "quick":
    VARIANTS = {k: VARIANTS[k] for k in list(VARIANTS)[:1]}
out_json = os.path.join(ROOT, "output", f"channels_{a.tag}.json")
res = json.load(open(out_json)) if os.path.exists(out_json) else {}


def gap(m):
    return 100 * (m["quit/m rec"] - m["quit/m exp"])


t0 = time.time()
for name, kw in VARIANTS.items():
    if name in res:
        continue
    p = p0.replace(**kw)
    mb = run(p, cfg)[0]; mh = run(acyclical_husband(p), cfg)[0]
    mf = run(acyclical_finding(p), cfg)[0]; mboth = run(acyclical_finding(acyclical_husband(p)), cfg)[0]
    g, gh, gf, gb = gap(mb), gap(mh), gap(mf), gap(mboth)
    res[name] = dict(E=float(mb["E/pop"]), quit_exp=float(mb["quit/m exp"]), quit_rec=float(mb["quit/m rec"]),
                     gap=g, gap_acycH=gh, gap_acycF=gf, gap_both=gb,
                     precaution=(g - gh) / g, hoarding=(g - gf) / g, both_off=(g - gb) / g,
                     dE=float(mb["dE/pop rec-exp (pts)"]), dE_acycH=float(mh["dE/pop rec-exp (pts)"]),
                     dE_acycF=float(mf["dE/pop rec-exp (pts)"]))
    json.dump(res, open(out_json, "w"), indent=1)
    r = res[name]
    print(f"{name}: gap {g:+.2f} precaution {r['precaution']:.0%} hoarding {r['hoarding']:.0%} both {r['both_off']:.0%} "
          f"dE {r['dE']:+.2f} [{time.time() - t0:.0f}s]", flush=True)
L = [f"# Precautionary labor supply versus job hoarding (`{a.calib}`, parameters held fixed across variants)", "",
     "Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap "
     "when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when "
     "the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage "
     "cut, the wife's own cyclical job loss and interactions.", "",
     "| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for name, r in res.items():
    L.append(f"| {name} | {r['E']:.3f} | {r['quit_exp']:.4f} | {r['quit_rec']:.4f} | {r['gap']:+.2f} | {r['precaution']:.0%} | "
             f"{r['hoarding']:.0%} | {r['both_off']:.0%} | {r['dE']:+.2f} | {r['dE_acycH']:+.2f} |")
open(os.path.join(ROOT, "output", f"channels_{a.tag}.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
