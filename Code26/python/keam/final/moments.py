"""Moments of the simulated final model, computed on calendar months inside a window."""
from __future__ import annotations
import numpy as np
from .simulate import FinalSim, stationary

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
    # job finding: entries into employment per non-employed (all) and per unemployed searcher (s >= s_bar)
    stat_prev = np.full_like(sim.stat, 0); stat_prev[:, 1:] = sim.stat[:, :-1]
    U_prev = (stat_prev == 1); N_prev = (emp_prev == 0); N_prev[:, 0] = False
    entry = (sim.emp == 1) & N_prev
    UE = _monthly_sums(sim, entry & U_prev); NE = _monthly_sums(sim, entry)
    Us = _monthly_sums(sim, U_prev); Ns = _monthly_sums(sim, N_prev)
    out["UE/m exp"] = ratio(UE, Us, exp_); out["UE/m rec"] = ratio(UE, Us, rec)
    out["NE/m exp"] = ratio(NE, Ns, exp_); out["NE/m rec"] = ratio(NE, Ns, rec)
    # entry from non-participation (previous month non-employed with s < s_bar), the CPS N->E flow of married women
    NP_prev = (stat_prev == 2); NP_prev[:, 0] = False
    NPE = _monthly_sums(sim, entry & NP_prev); NPs = _monthly_sums(sim, NP_prev)
    out["N->E/m exp"] = ratio(NPE, NPs, exp_); out["N->E/m rec"] = ratio(NPE, NPs, rec)
    # cyclicality as a standard deviation of the log rate under the two-state aggregate process:
    # |log(rate_rec / rate_exp)| sqrt(pi_exp pi_rec)  (data: 0.0686 for women's UE rate, 0.0765 for men's)
    piz = stationary(p.piz); wz = float(np.sqrt(piz[0] * piz[1]))
    sdlog = lambda a, b: float(abs(np.log(max(a, 1e-12) / max(b, 1e-12))) * wz)
    out["sd log UE (women)"] = sdlog(out["UE/m rec"], out["UE/m exp"])
    out["sd log NE (women)"] = sdlog(out["NE/m rec"], out["NE/m exp"])
    hprev = np.full_like(sim.hstat, 0); hprev[:, 1:] = sim.hstat[:, :-1]
    hU_prev = (hprev == 2); hU_prev[:, 0] = False
    HUE = _monthly_sums(sim, (sim.hstat == 1) & hU_prev); HU = _monthly_sums(sim, hU_prev)
    out["sd log UE (husband)"] = sdlog(ratio(HUE, HU, rec), ratio(HUE, HU, exp_))
    out["wife share exp"] = float(IW[exp_].sum() / (IW[exp_] + IH[exp_]).sum())
    out["wife share rec"] = float(IW[rec].sum() / (IW[rec] + IH[rec]).sum())
    out["HH income rec/exp - 1 (%)"] = 100 * (ratio(IW + IH, pop, rec) / ratio(IW + IH, pop, exp_) - 1)
    # wage gap: hourly wage ratio; the husband is assumed to work 2,000 hours (0.5 of the endowment)
    wE = (sim.wage * sim.emp * inwin).sum() / max((sim.emp * inwin).sum(), 1)
    hE_mask = (sim.hstat == 0) & inwin
    yHm = (sim.inc_h * hE_mask).sum() / max(hE_mask.sum(), 1)
    out["wage gap (hourly ratio)"] = float(0.5 * wE / yHm)
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
    # added-worker responses to the husband's job loss: the wife's employment and hours 1-12 months
    # after the loss relative to the 12 months before (event study on hloss), by aggregate state at the loss
    if idx.size:
        emp_b = np.stack([sim.emp[i, j - k] for k in range(1, 13)], 1).mean(1)
        emp_a = np.stack([sim.emp[i, j + k] for k in range(0, 12)], 1).mean(1)
        hrs_b = np.stack([sim.hours[i, j - k] for k in range(1, 13)], 1).mean(1)
        hrs_a = np.stack([sim.hours[i, j + k] for k in range(0, 12)], 1).mean(1)
        out["added worker: wife E +12m after H loss (pp)"] = 100 * float((emp_a - emp_b).mean())
        out["added worker: wife E +12m, loss in rec (pp)"] = 100 * float((emp_a - emp_b)[zl == 1].mean()) if (zl == 1).any() else np.nan
        out["added worker: wife hours +12m after H loss (%)"] = 100 * float((hrs_a.mean() / max(hrs_b.mean(), 1e-9)) - 1)
    # joint monthly transitions (Guner, Kulikova and Valladares-Esteban, "Does the added worker effect
    # matter?"): the wife's probability of entering the labor force (from NiLF: non-employed and not
    # searching) or employment in the month the husband moves from E to U, relative to months in which he
    # stays employed. Data: entry into the labor force is 60% more likely when the husband loses his job.
    hprev = np.full_like(sim.hstat, 0); hprev[:, 1:] = sim.hstat[:, :-1]
    sprev = np.full_like(sim.stat, 0); sprev[:, 1:] = sim.stat[:, :-1]
    valid = inwin.copy(); valid[:, 0] = False
    nilf_prev = (sprev == 2) & valid
    h_EU = (hprev == 0) & (sim.hstat == 2); h_EE = (hprev == 0) & (sim.hstat == 0)
    inLF = (sim.stat <= 1); inE = (sim.emp == 1)
    def cond(num, den):
        return float((num & den).sum() / max(den.sum(), 1))
    out["AWE: P(wife NiLF->LF | H E->U)"] = cond(inLF, nilf_prev & h_EU)
    out["AWE: P(wife NiLF->LF | H stays E)"] = cond(inLF, nilf_prev & h_EE)
    out["AWE: P(wife NiLF->E | H E->U)"] = cond(inE, nilf_prev & h_EU)
    out["AWE: P(wife NiLF->E | H stays E)"] = cond(inE, nilf_prev & h_EE)
    r = out["AWE: P(wife NiLF->LF | H stays E)"]
    out["AWE: LF entry ratio (H E->U / stays E)"] = out["AWE: P(wife NiLF->LF | H E->U)"] / r if r > 0 else np.nan
    # over the following 12 months: entry into the labor force by 12 months after the husband's loss
    if idx.size:
        lf_a = np.stack([(sim.stat[i, j + k] <= 1) for k in range(0, 12)], 1).any(1)
        nilf_at = (sim.stat[i, np.maximum(j - 1, 0)] == 2)
        out["AWE: P(NiLF wife enters LF within 12m | H loss)"] = float(lf_a[nilf_at].mean()) if nilf_at.any() else np.nan
    # employment and hours of wives by the husband's current state (E / R / U), within the window
    for name, ys in [("husband E", 0), ("husband R", 1), ("husband U", 2)]:
        msk = (sim.hstat == ys) & inwin
        out[f"wife E | {name}"] = float((sim.emp * msk).sum() / max(msk.sum(), 1))
        me = msk & (sim.emp == 1)
        out[f"wife hours|E | {name}"] = float((sim.hours * me).sum() / max(me.sum(), 1))
    out["mean assets/monthly HH inc"] = float((sim.a * inwin).sum() / max(inwin.sum(), 1)) / max(ratio(IW + IH, pop, sel), 1e-9)
    out["share e at cap"] = float(((sim.e >= p.e_max - 1e-6) * inwin).sum() / max(inwin.sum(), 1))
    out["share rec months"] = float(rec.sum() / sel.sum())
    return out
