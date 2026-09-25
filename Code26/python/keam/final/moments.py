"""Moments of the simulated final model, computed on calendar months inside a window."""
from __future__ import annotations
import numpy as np
from .simulate import FinalSim

HOURS_PER_YEAR = 4000.0    # 16 hours/day endowment: h = 0.4 is 1,600 hours (paper p.23)


def _monthly_sums(sim: FinalSim, X, mask=None):
    """Sum X (N_ind, L) by calendar month, optionally restricted to a mask."""
    cal = sim.calendar
    w = X.astype(float) if mask is None else X.astype(float) * mask
    return np.bincount(cal.ravel(), weights=w.ravel(), minlength=sim.T)[: sim.T]


def moments_final(sim: FinalSim, window=None) -> dict:
    p = sim.params; L = sim.cfg.L; T = sim.T
    if window is None:
        window = sim.cfg.window or (L, T)
    lo, hi = window
    z = sim.zpath
    alive = np.ones_like(sim.emp)
    cal = sim.calendar
    inwin = (cal >= lo) & (cal < hi)
    pop = _monthly_sums(sim, alive)
    E = _monthly_sums(sim, sim.emp)
    U = _monthly_sums(sim, (sim.stat == 1))
    # flows: previous-month employment
    emp_prev = np.zeros_like(sim.emp); emp_prev[:, 1:] = sim.emp[:, :-1]
    E_prev = _monthly_sums(sim, emp_prev)
    Q = _monthly_sums(sim, sim.quit * emp_prev)
    EN = _monthly_sums(sim, ((sim.emp == 0) * emp_prev))
    H = _monthly_sums(sim, sim.hours * sim.emp)
    IW = _monthly_sums(sim, sim.inc_w); IH = _monthly_sums(sim, sim.inc_h)
    months = np.arange(T); sel = (months >= lo) & (months < hi)
    rec = sel & (z == 1); exp_ = sel & (z == 0)

    def ratio(num, den, m):
        return float(np.mean(num[m] / np.maximum(den[m], 1)))

    out = {}
    out["E/pop"] = ratio(E, pop, sel); out["E/pop exp"] = ratio(E, pop, exp_); out["E/pop rec"] = ratio(E, pop, rec)
    out["dE/pop rec-exp (pts)"] = 100 * (out["E/pop rec"] - out["E/pop exp"])
    out["hours|E"] = ratio(H, E, sel)
    out["U rate"] = ratio(U, E + U, sel)
    out["quit/m exp"] = ratio(Q, E_prev, exp_); out["quit/m rec"] = ratio(Q, E_prev, rec)
    out["E->nonE/m exp"] = ratio(EN, E_prev, exp_); out["E->nonE/m rec"] = ratio(EN, E_prev, rec)
    out["wife share exp"] = float(IW[exp_].sum() / (IW[exp_] + IH[exp_]).sum())
    out["wife share rec"] = float(IW[rec].sum() / (IW[rec] + IH[rec]).sum())
    out["HH income rec/exp - 1 (%)"] = 100 * (ratio(IW + IH, pop, rec) / ratio(IW + IH, pop, exp_) - 1)
    # wage gap: full-time-equivalent monthly earnings of employed wives / husband income when employed
    wE = (sim.wage * sim.emp * inwin).sum() / max((sim.emp * inwin).sum(), 1)
    hE_mask = (sim.hstat == 0) & inwin
    yHm = (sim.inc_h * hE_mask).sum() / max(hE_mask.sum(), 1)
    out["wage gap (FTE earnings ratio)"] = float(0.4 * wE / yHm)
    # careers: annual hours over ages 25-54 (life months 0..359), cohorts fully inside the window
    m0, m1, _ = p.age_months
    full = (sim.entry + m0 + m1 <= hi) & (sim.entry >= lo - (m0 + m1))
    hrs_y = HOURS_PER_YEAR * sim.hours[full, : m0].mean(axis=1)
    hrs_m = HOURS_PER_YEAR * sim.hours[full, m0: m0 + m1].mean(axis=1)
    hrs_all = HOURS_PER_YEAR * sim.hours[full, : m0 + m1].mean(axis=1)
    lifecycle = (hrs_m >= 1500) & (hrs_y < 600)
    career = (hrs_all >= 1500) & ~lifecycle
    pt = (hrs_all >= 400) & (hrs_all < 1500) & ~lifecycle
    nilf = (hrs_all < 400) & ~lifecycle
    n = full.sum()
    out["share Lifecycle"] = lifecycle.mean(); out["share Career"] = career.mean()
    out["share PT"] = pt.mean(); out["share NiLF"] = nilf.mean(); out["n careers"] = int(n)
    # consumption response to husband's job loss: 12 months after vs 12 months before, by Z at loss
    idx = np.argwhere(sim.hloss[:, 12: L - 12] == 1)
    if idx.size:
        i, j = idx[:, 0], idx[:, 1] + 12
        before = np.stack([sim.cons[i, j - k] for k in range(1, 13)], 1).mean(1)
        after = np.stack([sim.cons[i, j + k] for k in range(0, 12)], 1).mean(1)
        zl = z[np.minimum(sim.entry[i] + j, T - 1)]
        chg = after / before - 1
        out["cons drop at H job loss exp (%)"] = 100 * float(chg[zl == 0].mean()) if (zl == 0).any() else np.nan
        out["cons drop at H job loss rec (%)"] = 100 * float(chg[zl == 1].mean()) if (zl == 1).any() else np.nan
    out["mean assets/monthly HH inc"] = float((sim.a * inwin).sum() / max(inwin.sum(), 1)) / max(ratio(IW + IH, pop, sel), 1e-9)
    out["share e at cap"] = float(((sim.e >= p.e_max - 1e-6) * inwin).sum() / max(inwin.sum(), 1))
    out["share rec months"] = float(rec.sum() / sel.sum())
    return out
