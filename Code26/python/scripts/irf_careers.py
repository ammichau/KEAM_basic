"""Impulse responses by career type (paper Fig. 8-9): employment rate and quit rate in the months
around recession starts, from the calibrated baseline simulated on the NBER recession dates.

usage: python scripts/irf_careers.py --calib output/final_calib_full.json
Writes output/figures/fig6_irf_employment.png/svg, fig7_irf_quits.png/svg and output/irf_careers.json
"""
import sys, os, json, argparse, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from keam.final import FinalParams, SimConfigFinal, solve_all, simulate_final
from keam.final.calibrate import apply_params, params_from_calib
from keam.final.moments import HOURS_PER_YEAR

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_full.json")
ap.add_argument("--tag", default="", help="figure sub-directory suffix and JSON name suffix")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); PY = os.path.join(HERE, ".."); OUT = os.path.join(PY, "output", "figures" + (("_" + a.tag) if a.tag else ""))
os.makedirs(OUT, exist_ok=True)
tagsfx = ("_" + a.tag) if a.tag else ""
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.linewidth": 1, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})

p = params_from_calib(json.load(open(os.path.join(PY, a.calib))))
sol = solve_all(p)
cfg = SimConfigFinal(N=300, n_cohorts=65, zmode="nber", seed=2024)
sim = simulate_final(p, sol, cfg)
L, T = cfg.L, sim.T
m0, m1, _ = p.age_months
# careers for women whose 25-54 span is inside the simulation
full = (sim.entry + m0 + m1 <= T)
hrs_y = HOURS_PER_YEAR * sim.hours[:, :m0].mean(axis=1); hrs_m = HOURS_PER_YEAR * sim.hours[:, m0:m0 + m1].mean(axis=1)
hrs_all = HOURS_PER_YEAR * sim.hours[:, :m0 + m1].mean(axis=1)
lifecycle = (hrs_m >= 1500) & (hrs_y < 600)
career = (hrs_all >= 1500) & ~lifecycle; pt = (hrs_all >= 400) & (hrs_all < 1500) & ~lifecycle; nilf = (hrs_all < 400) & ~lifecycle
groups = {"All": full, "Part-time": full & pt, "Life-cycle": full & lifecycle, "Mostly NiLF": full & nilf, "Career": full & career}
cal = sim.calendar
emp_prev = np.zeros_like(sim.emp); emp_prev[:, 1:] = sim.emp[:, :-1]
def monthly(X, mask):
    w = X.astype(float) * mask[:, None]
    return np.bincount(cal.ravel(), weights=w.ravel(), minlength=T)[:T]
z = sim.zpath
starts = [t for t in range(1, T) if z[t] == 1 and z[t - 1] == 0 and t >= (1973 - 1955) * 12]
pre, post = 6, 9
irf = {}
for g, mask in groups.items():
    pop = monthly(np.ones_like(sim.emp), mask); E = monthly(sim.emp, mask)
    Ep = monthly(emp_prev, mask); Q = monthly(sim.quit * emp_prev, mask)
    e_rate = E / np.maximum(pop, 1); q_rate = Q / np.maximum(Ep, 1)
    k3 = np.ones(3) / 3.0                                # 3-month centred moving average
    e_rate = np.convolve(e_rate, k3, mode='same'); q_rate = np.convolve(q_rate, k3, mode='same')
    e_dev = []; q_dev = []
    for s in starts:
        if s - pre < 0 or s + post >= T: continue
        e_base = e_rate[s - pre: s].mean(); q_base = q_rate[s - pre: s].mean()
        e_dev.append(100 * (e_rate[s - pre: s + post + 1] - e_base)); q_dev.append(100 * (q_rate[s - pre: s + post + 1] - q_base))
    irf[g] = dict(emp=np.mean(e_dev, axis=0).tolist(), quit=np.mean(q_dev, axis=0).tolist(), n_women=int(mask.sum()))
irf["_meta"] = dict(recession_starts=[int(s) for s in starts], pre=pre, post=post, note="deviation from the mean of the 6 months before the recession start, averaged over recessions 1973-2007; 3-month centred moving average; percentage points")
json.dump(irf, open(os.path.join(PY, "output", f"irf_careers{tagsfx}.json"), "w"), indent=1)
h = np.arange(-pre, post + 1)

def fig(key, ylabel, title, name):
    f, ax = plt.subplots(figsize=(6.6, 3.8))
    for i, g in enumerate(groups):
        ax.plot(h, irf[g][key], color=SERIES[i], linewidth=2 if g != "All" else 2.5, label=g)
    ax.axhline(0, color=GRID, linewidth=1); ax.axvline(0, color=GRID, linewidth=1)
    ax.set_xlabel("months since the start of the recession"); ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=11); ax.legend(fontsize=8.5, ncol=2)
    ax.grid(True, axis="y", color=GRID, linewidth=1); ax.set_axisbelow(True)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.tick_params(length=0)
    f.text(0.01, 0.005, f"NBER dates, {len(starts)} recessions; deviation from the 6 pre-recession months; 3-month moving average", color=INK2, fontsize=7.5)
    f.tight_layout(); f.savefig(os.path.join(OUT, name + ".png"), dpi=160); f.savefig(os.path.join(OUT, name + ".svg")); plt.close(f)

fig("emp", "employment rate, pp deviation", "Employment by career type after a recession starts", "fig6_irf_employment")
fig("quit", "monthly quit rate, pp deviation", "Quits by career type after a recession starts", "fig7_irf_quits")
for g in groups:
    print(f"{g:12s} n={irf[g]['n_women']:5d}  employment at +6: {irf[g]['emp'][pre + 6]:+.2f} pp   quits at +3: {irf[g]['quit'][pre + 3]:+.2f} pp")
