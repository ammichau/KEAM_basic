"""Paper figures from the calibrated final model.

  fig1_quit.png/svg     quit probability over experience by state (slides p.22)
  fig2_search.png/svg   search intensity over experience by state (slides p.23)
  fig3_hours.png/svg    hours over experience by state (slides p.45)
  fig4_cohorts.png/svg  refined cohort accounting: recession employment drop, quit rates
  fig5_mechanism.png/svg decomposition of the recession quit drop (counterfactuals)

usage: python scripts/figures.py --calib output/final_calib_full.json
Representative type: median productivity, median permanent cost, lowest life-cycle cost;
age group 40-54; assets at the population median.
"""
import sys, os, json, argparse, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from keam.final import FinalParams
from keam.final.params import make_types
from keam.final.calibrate import apply_params

ap = argparse.ArgumentParser()
ap.add_argument("--calib", default="output/final_calib_full.json")
ap.add_argument("--results", default="output/final_results_full.json")
ap.add_argument("--cohorts", default="output/cohorts_refined_full.json")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); PY = os.path.join(HERE, ".."); OUT = os.path.join(PY, "output", "figures")
os.makedirs(OUT, exist_ok=True)

# ---- palette (validated reference instance, light mode) and chrome
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.linewidth": 1, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})

def style(ax):
    ax.grid(True, axis="y", color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=160)
    fig.savefig(os.path.join(OUT, name + ".svg"))
    plt.close(fig)

# ---- solve all types; policies are averaged over types (equal weights) at age 40-54, median assets
from keam.final.solve import solve_all
p = apply_params(FinalParams(), json.load(open(os.path.join(PY, a.calib)))["x"])
sol = solve_all(p)
eg, ag, hg, sg = p.egrid, p.agrid, p.hgrid, p.sgrid
kT, wT = p.kT_nodes()
tau = 1
ia = int(np.argmin(np.abs(ag - 1.47)))
states = [((0, 0), "Normal times"), ((2, 0), "Husband unemployed"), ((0, 1), "Recession"), ((2, 1), "Husband unemployed + recession")]
sub = f"average over the 100 types; age 40-54; assets {ag[ia]:.2f} (population median)"

def panel(title, ylabel, series, name, ylim=None):
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for (y, z), lab in states:
        ax.plot(eg, series(y, z), color=SERIES[[s[1] for s in states].index(lab)], linewidth=2, solid_capstyle="round", solid_joinstyle="round", label=lab)
    ax.set_xlabel("experience e"); ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=11)
    if ylim: ax.set_ylim(*ylim)
    ax.legend(loc="best", fontsize=8.5); style(ax)
    fig.text(0.01, 0.005, sub, color=INK2, fontsize=7.5)
    save(fig, name)

def quit_prob(y, z):
    d = sol.VE[:, tau, :, ia, y, z] - sol.VN[:, tau, :, ia, y, z]          # (nK, nE)
    return np.mean(np.sum(wT[None, None, :] * ((d[:, :, None] - kT[None, None, :]) < 0), axis=2), axis=0)

panel("Quit probability of an employed woman", "monthly quit probability", quit_prob, "fig1_quit", ylim=(0, None))
panel("Search intensity of a non-employed woman", "search intensity s", lambda y, z: sol.gS[:, tau, :, ia, y, z].mean(axis=0), "fig2_search", ylim=(0, None))
panel("Hours of an employed woman", "hours (share of time endowment)", lambda y, z: sol.gH[:, tau, :, ia, y, z].mean(axis=0), "fig3_hours", ylim=(0, None))

# ---- cohort trend (refined accounting)
coh = json.load(open(os.path.join(PY, a.cohorts)))
names = list(coh); x = np.arange(len(names))
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
ax = axes[0]
drop = [-coh[n]["m"]["dE/pop rec-exp (pts)"] for n in names]
ax.bar(x, drop, width=0.5, color=SERIES[0], edgecolor=SURF, linewidth=2)
for xi, v in zip(x, drop):
    ax.text(xi, v + 0.04, f"{v:.2f}", ha="center", va="bottom", fontsize=8.5, color=INK)
ax.set_xticks(x); ax.set_xticklabels(names); ax.set_ylabel("percentage points")
ax.set_title("Employment drop in recessions, by cohort", loc="left", fontsize=11); style(ax)
ax = axes[1]
ax.plot(x, [100 * coh[n]["m"]["quit/m exp"] for n in names], color=SERIES[0], linewidth=2, marker="o", markersize=6, markeredgecolor=SURF, markeredgewidth=2, label="expansions")
ax.plot(x, [100 * coh[n]["m"]["quit/m rec"] for n in names], color=SERIES[1], linewidth=2, marker="o", markersize=6, markeredgecolor=SURF, markeredgewidth=2, label="recessions")
ax.set_xticks(x); ax.set_xticklabels(names); ax.set_ylabel("monthly quit rate, %"); ax.set_ylim(0, None)
ax.set_title("Quit rate, by cohort", loc="left", fontsize=11); ax.legend(fontsize=8.5); style(ax)
fig.text(0.01, 0.005, "cohort accounting with the cost scale and tau_w solved jointly (output/cohorts_refined_full.json)", color=INK2, fontsize=7.5)
save(fig, "fig4_cohorts")

# ---- mechanism decomposition
res = json.load(open(os.path.join(PY, a.results)))
cf = res["counterfactuals"]
gap = lambda m: 100 * (m["quit/m rec"] - m["quit/m exp"])
labels = ["baseline", "acyclical husband risk", "acyclical job finding", "no recession wage cut", "acyclical own job loss"]
vals = [gap(cf[l]) for l in labels]
fig, ax = plt.subplots(figsize=(7.6, 3.4))
yy = np.arange(len(labels))[::-1]
ax.barh(yy, vals, height=0.5, color=SERIES[0], edgecolor=SURF, linewidth=2)
for yi, v in zip(yy, vals):
    ax.text(v - 0.015, yi, f"{v:+.2f}", ha="right", va="center", fontsize=8.5, color=INK)
ax.set_xlim(min(vals) - 0.12, 0)
ax.set_yticks(yy); ax.set_yticklabels(labels); ax.set_xlabel("recession minus expansion monthly quit rate, percentage points")
fig.subplots_adjust(left=0.28)
ax.set_title("What drives the fall in quits during recessions", loc="left", fontsize=11); style(ax)
ax.grid(True, axis="x", color=GRID, linewidth=1); ax.grid(False, axis="y")
save(fig, "fig5_mechanism")
print("figures written to", OUT)
