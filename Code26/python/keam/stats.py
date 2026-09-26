"""Aggregate statistics: translation of SimplerMod_May17_sim.m lines 559-1118.

Three blocks are reproduced: (1) cross-section by age group, (2) career types,
(3) business-cycle means/variances.  Column labels follow the MATLAB variable
names, with the label swaps of the xlswrite block documented in DEPARTURES.md.
"""
from __future__ import annotations
import numpy as np
from .simulate import SimResult


def _age_masks(r: SimResult):
    """Boolean masks over (individual, quarter) for the three age groups, restricted
    to the person-quarters the MATLAB loops visit (cohorts 0..Ngen-2, tt < Tsim)."""
    m = np.zeros((3, r.Nsim_i, r.Tsim + 1), bool)
    Nsim = r.cfg.Nsim
    y1, y2 = r.cfg.years_per_age
    t4 = r.params.tsize
    bounds = [(0, y1 * t4), (y1 * t4, (y1 + y2) * t4), ((y1 + y2) * t4, r.Tsim_i)]
    for gen in range(r.Ngen - 1):
        ii = slice(gen * Nsim, (gen + 1) * Nsim)
        for a, (lo, hi) in enumerate(bounds):
            tt = np.arange(gen * 4 + lo, gen * 4 + hi)
            tt = tt[tt < r.Tsim]
            if tt.size:
                m[a, ii, tt[0]: tt[-1] + 1] = True
    return m


def cross_section(r: SimResult) -> dict:
    """Lines 567-818.  Returns dict of arrays indexed by age group (young, middle, old)."""
    m = _age_masks(r)
    Ky = r.Ky_i[np.arange(r.Nsim_i) % r.cfg.Nsim]
    ftK = r.ftK_i[np.arange(r.Nsim_i) % r.cfg.Nsim]
    ktype = (Ky + ftK)[:, None] * np.ones((1, r.Tsim + 1))
    hE = (r.Hstat == 1)
    out = {}

    def s(x, a):
        return float(x[m[a]].sum()) if x.dtype != bool else float(x[m[a]].sum())

    for a in range(3):
        pop = m[a].sum()
        E = s(r.EmpI, a); U = s(r.UnempI, a); N = s(r.NiLF, a)
        d = dict(pop=pop, Erate=E / pop, Urate=U / pop, Nrate=N / pop,
                 mHours=s(r.Hours * r.EmpI, a) / E, mWage_E=s(r.Wage * r.EmpI, a) / E,
                 mExp_E=s(r.exp * r.EmpI, a) / E, HstatE_E=s(hE * r.EmpI, a) / E,
                 Ktype_E=s(ktype * r.EmpI, a) / E,
                 EErate=s(r.EE, a) / E, ENrate=s(r.EN, a) / E, EUrate=s(r.EU, a) / E,
                 Qrate=s(r.Quit, a) / E, jLossrate=s(r.Wloss, a) / E,
                 mSearch_U=s(r.Srch * r.UnempI, a) / U if U else np.nan,
                 mExp_U=s(r.exp * r.UnempI, a) / U if U else np.nan,
                 mWage_U=s(r.Wage * r.UnempI, a) / U if U else np.nan,
                 HstatE_U=s(hE * r.UnempI, a) / U if U else np.nan,
                 Ktype_U=s(ktype * r.UnempI, a) / U if U else np.nan,
                 UErate=s(r.UE, a) / U if U else np.nan,
                 mSearch_N=s(r.Srch * r.NiLF, a) / N, mExp_N=s(r.exp * r.NiLF, a) / N,
                 mWage_N=s(r.Wage * r.NiLF, a) / N, HstatE_N=s(hE * r.NiLF, a) / N,
                 Ktype_N=s(ktype * r.NiLF, a) / N, NErate=s(r.NE, a) / N)
        for k, v in d.items():
            out.setdefault(k, []).append(v)
    return {k: np.array(v) for k, v in out.items()}


def careers(r: SimResult, drop_undefined: bool | None = None) -> dict:
    """Lines 880-933: career taxonomy.  Returns shares and the per-individual labels
    (0 undefined, 1 part-time, 2 life-cycle, 3 career, 4 mostly NiLF)."""
    if drop_undefined is None:
        drop_undefined = r.options.drop_undefined_careers
    Nsim = r.cfg.Nsim
    n_ind = r.NgenCareer * Nsim
    ageEmp = np.zeros((n_ind, 3)); ageAlive = np.zeros((n_ind, 3)); ageFT = np.zeros((n_ind, 3))
    for gen in range(r.NgenCareer - 1):
        ii = np.arange(gen * Nsim, (gen + 1) * Nsim)
        for it in range(11, r.Tsim_i):                     # MATLAB it = 12:Tsim_i
            tt = gen * 4 + it
            if tt >= r.Tsim:
                break
            a = r.Age[ii, tt]
            ok = a > 0
            for g in range(3):
                sel = ii[ok & (a == g + 1)]
                ageEmp[sel, g] += r.EmpI[sel, tt]
                ageAlive[sel, g] += r.Alive[sel, tt]
                ageFT[sel, g] += r.FT[sel, tt]
    with np.errstate(invalid="ignore", divide="ignore"):
        ageEmpR = np.where(ageEmp > 0, ageEmp / ageAlive, 0.0)
        EmpR = ageEmpR[:, :2].mean(axis=1)
        ftshare = np.where(ageEmp.sum(1) > 0, ageFT.sum(1) / ageEmp.sum(1), 0.0)
    label = np.zeros(n_ind, int)
    n_class = (r.NgenCareer - 1) * Nsim
    for k in range(n_class):
        if EmpR[k] > 0.8:
            label[k] = 3 if ftshare[k] > 0.7 else 1
        elif ageEmpR[k, 0] < r.cfg.cyclefactor * ageEmpR[k, 1]:
            label[k] = 2
        elif EmpR[k] < r.cfg.NiLFdef:
            label[k] = 4
        else:
            label[k] = 0
    lab = label[:n_class]
    counts = np.array([(lab == c).sum() for c in range(5)])
    denom = counts[1:].sum() if drop_undefined else counts.sum()
    return dict(PT=counts[1] / denom, Cycle=counts[2] / denom, Career=counts[3] / denom,
                NiLF=counts[4] / denom, Undefined=counts[0] / counts.sum(),
                counts=counts, label=label, EmpR=EmpR, ageEmpR=ageEmpR)


def cycle(r: SimResult, expansion_includes_recessions: bool | None = None) -> dict:
    """Lines 955-1094: recession / 'expansion' means and variances of aggregate series."""
    if expansion_includes_recessions is None:
        expansion_includes_recessions = r.options.expansion_includes_recessions
    sl = slice(r.ScratchT - 1, r.Tsim)                    # MATLAB ScratchT:Tsim (1-based)
    z = r.zsim[sl]
    RecI = (z == 2).astype(float)
    ExpI = np.ones_like(RecI) if expansion_includes_recessions else (z == 1).astype(float)
    Pop_t = r.Alive[:, sl].sum(0)
    H_Emp = ((r.Hstat == 1) | (r.Hstat == 2)) & (r.Alive == 1)
    W_NonEmp = (r.NiLF == 1) | (r.UnempI == 1)
    W_NonEE = (r.NE == 1) | (r.UE == 1)
    W_ENonE = (r.EN == 1) | (r.EU == 1)
    E_t = r.EmpI[:, sl].sum(0)
    N_t = r.NiLF[:, sl].sum(0)
    U_t = r.UnempI[:, sl].sum(0)
    NE_t = W_NonEmp[:, sl].sum(0)

    def stat(num, den):
        x = num / den
        return dict(rec=float((x * RecI).sum() / RecI.sum()),
                    exp=float((x * ExpI).sum() / ExpI.sum()),
                    var=float(np.var(x, ddof=1) * 100))

    out = {}
    out["H_Emp"] = stat(H_Emp[:, sl].sum(0), Pop_t)
    out["Emp"] = stat(E_t, Pop_t)
    out["UnEmp"] = stat(U_t, Pop_t)
    out["NiLF"] = stat(N_t, Pop_t)
    out["NonEmp"] = stat(NE_t, N_t)                       # MATLAB divides by NiLF count (sic)
    out["EE"] = stat(r.EE[:, sl].sum(0), E_t)
    out["EN"] = stat(r.EN[:, sl].sum(0), E_t)
    out["EU"] = stat(r.EU[:, sl].sum(0), E_t)
    out["EnonE"] = stat(W_ENonE[:, sl].sum(0), E_t)
    out["NonEE"] = stat(W_NonEE[:, sl].sum(0), NE_t)
    out["Hours"] = stat(r.Hours[:, sl].sum(0), E_t)
    out["Quit"] = stat(r.Quit[:, sl].sum(0), E_t)
    out["Search"] = stat(r.Srch[:, sl].sum(0), U_t)
    out["Wage"] = stat((r.Wage * r.EmpI)[:, sl].sum(0), E_t)
    out["IncHH"] = stat(r.Inc_hh[:, sl].sum(0), Pop_t)
    out["IncW"] = stat(r.Inc_w[:, sl].sum(0), Pop_t)
    ws = stat(r.Inc_w[:, sl].sum(0), r.Inc_hh[:, sl].sum(0))
    out["Wshare"] = dict(rec=out["IncW"]["rec"] / out["IncHH"]["rec"],
                         exp=out["IncW"]["exp"] / out["IncHH"]["exp"], var=ws["var"])
    for k in ("UnEmp", "NiLF", "NonEmp"):
        out[k]["var"] = 0.0                                # MATLAB exports 0 for these
    return out


def summary_table(r: SimResult) -> str:
    cs = cross_section(r); ca = careers(r); cy = cycle(r)
    lines = ["Cross-section (young, middle, old):"]
    for k in ["Erate", "Urate", "Nrate", "mHours", "mWage_E", "mExp_E", "EErate", "ENrate",
              "EUrate", "Qrate", "jLossrate", "NErate"]:
        lines.append(f"  {k:10s} " + "  ".join(f"{v:9.5f}" for v in cs[k]))
    lines.append("Careers: " + ", ".join(f"{k}={ca[k]:.4f}" for k in ["PT", "Cycle", "Career", "NiLF", "Undefined"]))
    lines.append("Cycle (recession mean, 'expansion' mean, variance*100):")
    for k, v in cy.items():
        lines.append(f"  {k:8s} {v['rec']:9.5f} {v['exp']:9.5f} {v['var']:9.5f}")
    return "\n".join(lines)
