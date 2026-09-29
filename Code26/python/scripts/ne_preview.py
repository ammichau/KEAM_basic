"""Preview of the non-search offer arrival rate (lam_n) at fixed parameters: for each lam_n level, the coarse-grid
moments (N->E rate, quits, employment, never-working share) and the precaution / hoarding split of the recession quit
gap (husband's risk acyclical; wife's job finding acyclical, search and non-search offers).

usage: python scripts/ne_preview.py --calib output/final_calib_v7cmb_full.json --fixed beta=0.993 --fixed nAc=100
       --levels 0,0.03,0.06 --ratio 1.0 [--full] --tag v7cmb
"""
import sys, os, json, argparse, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
from keam.final import FinalParams, SimConfigFinal
from keam.final import calibrate as C
from keam.final.experiments import run, acyclical_husband, acyclical_finding

ap = argparse.ArgumentParser()
ap.add_argument("--calib", required=True)
ap.add_argument("--fixed", action="append", default=[])
ap.add_argument("--levels", default="0,0.03,0.06")
ap.add_argument("--ratio", type=float, default=1.0, help="recession / expansion ratio of the non-search arrival rate")
ap.add_argument("--full", action="store_true", help="100-type grid (default: coarse 27 types)")
ap.add_argument("--tag", default="preview")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
calib = json.load(open(os.path.join(ROOT, a.calib)))
base = FinalParams() if a.full else FinalParams(n_omega=3, n_kbar=3, n_km=3)
p0 = C.params_from_calib(calib, base)
for kv in a.fixed:
    n, v = kv.split("="); p0 = p0.replace(**{n: type(getattr(p0, n))(float(v))})
cfg = SimConfigFinal(N=60, n_cohorts=90)
KEYS = ["E/pop", "share NiLF", "share Career", "quit/m exp", "quit/m rec", "N->E/m exp", "N->E/m rec", "UE/m exp",
        "UE/m rec", "sd log UE (women)", "U rate", "dE/pop rec-exp (pts)", "mean assets/monthly HH inc"]
rows = []
for lev in [float(x) for x in a.levels.split(",")]:
    p = p0.replace(lam_n=(lev, a.ratio * lev))
    t0 = time.time()
    mb = run(p, cfg)[0]; mh = run(acyclical_husband(p), cfg)[0]; mf = run(acyclical_finding(p), cfg)[0]
    gap = lambda m: 100 * (m["quit/m rec"] - m["quit/m exp"])
    g0, gh, gf = gap(mb), gap(mh), gap(mf)
    prec = (g0 - gh) / g0 if g0 != 0 else float("nan"); hoard = (g0 - gf) / g0 if g0 != 0 else float("nan")
    row = dict(lam_n0=lev, lam_n_ratio=a.ratio, gap=g0, gap_acycH=gh, gap_acycF=gf, precaution=prec, hoarding=hoard,
               dE=mb["dE/pop rec-exp (pts)"], dE_acycH=mh["dE/pop rec-exp (pts)"], dE_acycF=mf["dE/pop rec-exp (pts)"],
               **{k: mb[k] for k in KEYS}, seconds=time.time() - t0)
    rows.append(row)
    print(f"lam_n0 {lev:.3f}: N->E {mb['N->E/m exp']:.4f}/{mb['N->E/m rec']:.4f}  quits {mb['quit/m exp']:.4f}/{mb['quit/m rec']:.4f}"
          f"  E {mb['E/pop']:.3f}  NiLF {mb['share NiLF']:.3f}  gap {g0:+.2f}  precaution {prec:.0%}  hoarding {hoard:.0%}"
          f"  dE {mb['dE/pop rec-exp (pts)']:+.2f} ({time.time() - t0:.0f}s)", flush=True)
out = os.path.join(ROOT, "output", f"ne_preview_{a.tag}")
json.dump(dict(calib=a.calib, fixed=a.fixed, full=a.full, rows=rows), open(out + ".json", "w"), indent=1)
md = [f"# Non-search offer arrival rate at fixed parameters (`{a.calib}`, fixed {a.fixed}, {'100' if a.full else '27'} types)", "",
      "Job finding = lam_n(Z) + lam_f(Z) s^nu. Quit gap = 100 x (recession - expansion quit rate). Precaution share = fall in the gap "
      "when the husband's risk is acyclical; hoarding share = fall when the wife's job finding (search and non-search offers) is acyclical.", "",
      "| lam_n0 | ratio | N->E exp | N->E rec | quit exp | quit rec | E/pop | NiLF | Career | UE exp | sd log UE | U rate | assets | gap | precaution | hoarding | dE | dE acyc. H | dE acyc. F |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    md.append(f"| {r['lam_n0']:.3f} | {r['lam_n_ratio']:.2f} | {r['N->E/m exp']:.4f} | {r['N->E/m rec']:.4f} | {r['quit/m exp']:.4f} | {r['quit/m rec']:.4f} | "
              f"{r['E/pop']:.3f} | {r['share NiLF']:.3f} | {r['share Career']:.3f} | {r['UE/m exp']:.3f} | {r['sd log UE (women)']:.4f} | {r['U rate']:.3f} | "
              f"{r['mean assets/monthly HH inc']:.1f} | {r['gap']:+.2f} | {r['precaution']:.0%} | {r['hoarding']:.0%} | {r['dE']:+.2f} | {r['dE_acycH']:+.2f} | {r['dE_acycF']:+.2f} |")
open(out + ".md", "w").write("\n".join(md) + "\n"); print("wrote", out + ".md")
