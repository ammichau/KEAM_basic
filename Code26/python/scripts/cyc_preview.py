"""Cyclicality of married women's quits and N->E entry at fixed parameters, in the convention of the UE target
(sd of the log rate under the two-state aggregate process, |log(rec/exp)| sqrt(pi_exp pi_rec)), for a set of
single-parameter variants around a calibration. Used to check which levers can bring the model's quit and N->E
cyclicality to the data (0.0262 and 0.0045, early window) before a calibration targets them.

usage: python scripts/cyc_preview.py --log output/final_calib_v7cnb_full.log --fixed kpr=1 ... [--full] --tag v7cnb
       (--log: the best evaluation in a calibration log; or --calib <json>)
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="")
ap.add_argument("--log", default="")
ap.add_argument("--fixed", action="append", default=[])
ap.add_argument("--full", action="store_true", help="100-type grid (default: coarse 27 types)")
ap.add_argument("--variants", default="base,ratio1.5,ratio2.0,nocut", help="comma list of variant keys")
ap.add_argument("--tag", default="preview")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
if a.log:
    rows = [json.loads(l) for l in open(os.path.join(ROOT, a.log)) if l.startswith("{")]
    best = min(rows, key=lambda r: r["obj"]); calib = dict(x=best["x"], fixed={}); src = f"{a.log} (evaluation {best['n']}, objective {best['obj']:.4f})"
else:
    calib = json.load(open(os.path.join(ROOT, a.calib))); src = a.calib
base = FinalParams() if a.full else FinalParams(n_omega=3, n_kbar=3, n_km=3)
for kv in a.fixed:
    n, v = kv.split("="); base = base.replace(**{n: type(getattr(base, n))(float(v))})
p0 = C.params_from_calib(calib, base)
ln = p0.lam_n
VAR = {
    "base": ("calibrated parameters", {}),
    "ratio1.5": ("non-search arrival x1.5 in recessions", dict(lam_n=(ln[0], 1.5 * ln[0]))),
    "ratio2.0": ("non-search arrival x2.0 in recessions", dict(lam_n=(ln[0], 2.0 * ln[0]))),
    "nocut": ("no recession wage cut for the wife (phi_rec 1)", dict(phi_rec=1.0)),
    "lamn2": ("non-search arrival level x2", dict(lam_n=(2 * ln[0], 2 * ln[1]))),
    "fall10": ("job-finding fall 10% (ratio 0.90)", dict(lam_f=(p0.lam_f[0], 0.90 * p0.lam_f[0]))),
}
cfg = SimConfigFinal(N=60, n_cohorts=90)
KEYS = ["E/pop", "share NiLF", "quit/m exp", "quit/m rec", "quit rec/exp", "sd log quit (women)", "N->E/m exp", "N->E/m rec",
        "N->E rec/exp", "sd log N->E (women)", "UE/m exp", "sd log UE (women)", "dE/pop rec-exp (pts)"]
res = []
def parse(k):   # combined variant, e.g. "lamnr=1.4;lamfr=0.7;phi=1": recession ratios of lam_n and lam_f, wife's wage cut
    d = dict(kv.split("=") for kv in k.split(";")); kw = {}; lab = []
    if "lamnr" in d:
        kw["lam_n"] = (ln[0], float(d["lamnr"]) * ln[0]); lab.append(f"non-search arrival x{float(d['lamnr']):.2f} in recessions")
    if "lamfr" in d:
        kw["lam_f"] = (p0.lam_f[0], float(d["lamfr"]) * p0.lam_f[0]); lab.append(f"job-finding ratio {float(d['lamfr']):.2f}")
    if "phi" in d:
        kw["phi_rec"] = float(d["phi"]); lab.append(f"phi_rec {float(d['phi']):.2f}")
    if "phiH" in d:
        kw["phi_rec_H"] = float(d["phiH"]); lab.append(f"phi_rec_H {float(d['phiH']):.2f}")
    return "; ".join(lab), kw


for k in [v for v in a.variants.split(",") if v]:
    lab, kw = VAR[k] if k in VAR else parse(k); t0 = time.time()
    m = run(p0.replace(**kw), cfg)[0]
    r = dict(variant=k, label=lab, **{q: float(m[q]) for q in KEYS}, seconds=time.time() - t0); res.append(r)
    print(f"{lab}: quits {r['quit/m exp']:.4f}/{r['quit/m rec']:.4f} sd {r['sd log quit (women)']:.4f}  N->E {r['N->E/m exp']:.4f}/"
          f"{r['N->E/m rec']:.4f} sd {r['sd log N->E (women)']:.4f}  sd UE {r['sd log UE (women)']:.4f}  dE {r['dE/pop rec-exp (pts)']:+.2f}"
          f"  E {r['E/pop']:.3f} ({r['seconds']:.0f}s)", flush=True)
out = os.path.join(ROOT, "output", f"cyc_preview_{a.tag}")
json.dump(dict(source=src, fixed=a.fixed, full=a.full, rows=res), open(out + ".json", "w"), indent=1)
md = [f"# Cyclicality of quits and N->E at fixed parameters (`{src}`, fixed {a.fixed}, {'100' if a.full else '27'} types)", "",
      "sd log = |log(rec/exp)| sqrt(pi_exp pi_rec), the convention of the UE target. Data (early window, trend-adjusted): "
      "sd log quit 0.0262 (ratio 0.928), sd log N->E 0.0045 (ratio 0.987), sd log UE 0.0686, employment drop -1.70.", "",
      "| variant | quit exp | quit rec | ratio | sd log quit | N->E exp | N->E rec | ratio | sd log N->E | sd log UE | dE (pts) | E/pop | NiLF |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
      "| data | 0.0226 | 0.0210 | 0.928 | 0.0262 | 0.0530 | 0.0524 | 0.987 | 0.0045 | 0.0686 | -1.70 | 0.620 | 0.220 |"]
for r in res:
    md.append(f"| {r['label']} | {r['quit/m exp']:.4f} | {r['quit/m rec']:.4f} | {r['quit rec/exp']:.3f} | {r['sd log quit (women)']:.4f} | "
              f"{r['N->E/m exp']:.4f} | {r['N->E/m rec']:.4f} | {r['N->E rec/exp']:.3f} | {r['sd log N->E (women)']:.4f} | "
              f"{r['sd log UE (women)']:.4f} | {r['dE/pop rec-exp (pts)']:+.2f} | {r['E/pop']:.3f} | {r['share NiLF']:.3f} |")
open(out + ".md", "w").write("\n".join(md) + "\n"); print("wrote", out + ".md")
