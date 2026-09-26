"""Figures for the channel decomposition: (fig8) precaution and hoarding shares of the recession quit
drop across calibration versions; (fig9) sensitivity of the two shares to the cyclical parameters
(from scripts/jacobian_channels.py, if its output exists). Output: output/figures_channels/.

usage: python scripts/figures_channels.py [--versions "adopted iid=ls,version 3=v3,..."] [--jacobian v4]
"""
import sys, os, json, argparse, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.simplefilter("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--versions", default="adopted iid (15% fall, UI 30%)=ls;7 wage types=om7;recession UI cut=ui;"
                                      "v3: UI cut + 10% fall=v3;v4: UI cut + UE target=v4;v4c: UI cut + 20% fall=v4c;v5: log utility=v5")
ap.add_argument("--jacobian", default="v4")
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); PY = os.path.join(HERE, ".."); OUT = os.path.join(PY, "output", "figures_channels")
os.makedirs(OUT, exist_ok=True)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.linewidth": 1, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})


def style(ax):
    ax.grid(True, axis="x", color=GRID, linewidth=1); ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), dpi=160); fig.savefig(os.path.join(OUT, name + ".svg")); plt.close(fig)


# ---- fig8: shares by version
rows = []
for item in a.versions.split(";"):
    label, tag = item.split("=")
    f = os.path.join(PY, "output", f"channels_{tag}.json")
    if os.path.exists(f):
        d = json.load(open(f)).get("baseline")
        if d:
            rows.append((label, d))
if rows:
    fig, ax = plt.subplots(figsize=(7.2, 0.55 * len(rows) + 1.6))
    y = np.arange(len(rows))[::-1]; h = 0.36
    ax.barh(y + h / 2, [100 * r["precaution"] for _, r in rows], h, color=SERIES[0], label="precautionary labor supply (husband's risk acyclical)")
    ax.barh(y - h / 2, [100 * r["hoarding"] for _, r in rows], h, color=SERIES[1], label="job hoarding (wife's job finding acyclical)")
    for yi, (_, r) in zip(y, rows):
        ax.text(100 * r["precaution"] + 0.6, yi + h / 2, f"{100*r['precaution']:.0f}", va="center", fontsize=8.5, color=INK2)
        ax.text(100 * r["hoarding"] + 0.6, yi - h / 2, f"{100*r['hoarding']:.0f}", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([lb for lb, _ in rows])
    ax.set_xlabel("share of the recession fall in the quit rate removed when the channel is switched off (%)")
    ax.set_title("Precautionary labor supply versus job hoarding by calibration version")
    ax.legend(loc="lower right", fontsize=8.5); style(ax)
    fig.text(0.01, 0.005, "scripts/channels.py; quit gap = recession minus expansion monthly quit rate; parameters at each version's calibrated values", color=INK2, fontsize=7.5)
    save(fig, "fig8_channels_by_version")
# ---- fig9: sensitivity of the shares
fj = os.path.join(PY, "output", f"jacobian_channels_{a.jacobian}.json")
if os.path.exists(fj):
    J = json.load(open(fj)); b = J["points"]["base"]; step = J["step"]
    names = [n for n in J["points"] if n != "base"]
    if names:
        dp = [100 * (J["points"][n]["precaution"] - b["precaution"]) / (100 * step) for n in names]
        dh = [100 * (J["points"][n]["hoarding"] - b["hoarding"]) / (100 * step) for n in names]
        order = np.argsort(np.abs(np.array(dp)) + np.abs(np.array(dh)))
        fig, ax = plt.subplots(figsize=(7.6, 0.5 * len(names) + 1.8))
        y = np.arange(len(names)); h = 0.36
        ax.barh(y + h / 2, [dp[i] for i in order], h, color=SERIES[0], label="precaution share")
        ax.barh(y - h / 2, [dh[i] for i in order], h, color=SERIES[1], label="hoarding share")
        ax.axvline(0, color=INK2, linewidth=1)
        ax.set_yticks(y); ax.set_yticklabels([names[i] for i in order], fontsize=8.5)
        ax.set_xlabel("change in the share, percentage points, per +1% of the quantity")
        ax.set_title(f"What moves the split: local sensitivity at the {a.jacobian} calibration")
        ax.legend(loc="lower right", fontsize=8.5); style(ax)
        fig.text(0.01, 0.005, f"scripts/jacobian_channels.py (+{100*step:.0f}% steps, other parameters fixed)", color=INK2, fontsize=7.5)
        save(fig, "fig9_channel_sensitivity")
print("saved to", OUT, os.listdir(OUT))
