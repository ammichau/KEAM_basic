"""Life-cycle panel simulation: translation of SimplerMod_May17_sim.m (lines 48-554).

The MATLAB simulation loops over birth cohorts (`in`), individuals (`i`) and
life-cycle quarters (`it`).  Here the loop over `it` is kept and everything is
vectorised over individuals.  All random draws are reproduced exactly with
:mod:`keam.matlab_rng`, so in faithful mode the panel is identical to MATLAB's.

Index conventions: individual ``ii = in*Nsim + i`` (0-based), calendar quarter
``tt = 4*in + it`` (0-based), type index of individual ``ii`` is ``ii % Nsim``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.stats import norm

from .params import Params, Options
from .solve import Solution
from . import functions as F
from .matlab_rng import rand, randi


@dataclass
class SimConfig:
    y0: int = 1955
    yT: int = 2019
    yStart: int = 1973
    yCareerEnd: int = 1995
    Nsim: int = 100          # individuals (types) per birth cohort
    age0: int = 25
    ninter: int = 10
    NiLFdef: float = 0.21
    NiLFbar: float = 0.2
    cyclefactor: float = 0.5
    FTbar: float = 0.39
    Nw_init: float = 0.08
    years_per_age: tuple = (14, 14)   # young, middle; old until Tsim_i
    Tsim_years: int = 39
    seeds: dict = field(default_factory=lambda: dict(
        ftW=111, ftK=222, Ky=222, exp0=207, hstat0=585, emp0=612, LH=812, LW=519, Q=911))


@dataclass
class SimResult:
    cfg: SimConfig
    params: Params
    options: Options
    # per-individual types
    ftW_i: np.ndarray; ftK_i: np.ndarray; Ky_i: np.ndarray; Ky_inx: np.ndarray
    # calendar
    year: np.ndarray; zsim: np.ndarray; ageT: np.ndarray
    Ngen: int; Tsim: int; Tsim_i: int; Nsim_i: int; ScratchT: int; NgenCareer: int
    # panel arrays (Nsim_i, Tsim+1)
    Alive: np.ndarray; Age: np.ndarray; AgeY: np.ndarray; exp: np.ndarray; Wage: np.ndarray
    Hstat: np.ndarray; Hloss: np.ndarray; Wloss: np.ndarray
    EmpI: np.ndarray; UnempI: np.ndarray; NiLF: np.ndarray; NempI: np.ndarray
    Quit: np.ndarray; Hours: np.ndarray; Srch: np.ndarray; FT: np.ndarray
    Inc_w: np.ndarray; Inc_hh: np.ndarray
    EE: np.ndarray; EU: np.ndarray; UE: np.ndarray; NE: np.ndarray; EN: np.ndarray
    NU: np.ndarray; UN: np.ndarray
    BirthV: dict
    ggH: np.ndarray; ggS: np.ndarray; ggQ: np.ndarray


# ----------------------------------------------------------------------------
def markov_stationary(P: np.ndarray) -> np.ndarray:
    """Literal translation of MarkovStationary.m (GTH algorithm), P[i, j] = Pr(i -> j)."""
    P = np.array(P, float, copy=True)
    ns = P.shape[0]
    n = ns
    while n > 1:
        n1 = n - 1
        s = P[n - 1, :n1].sum()
        P[:n1, n - 1] = P[:n1, n - 1] / s
        for n2 in range(n1, 0, -1):
            P[:n1, n2 - 1] = P[:n1, n2 - 1] + P[:n1, n - 1] * P[n - 1, n2 - 1]
        n -= 1
    pbar = np.ones(ns)
    for j in range(2, ns + 1):
        pbar[j - 1] = np.sum(pbar[: j - 1] * P[: j - 1, j - 1])
    return pbar / pbar.sum()


def _truncnorm_cdf(x, mu, sig, lo, hi):
    a, b = norm.cdf((lo - mu) / sig), norm.cdf((hi - mu) / sig)
    return (norm.cdf((x - mu) / sig) - a) / (b - a)


def _draw_types(cfg: SimConfig, p: Params, opts: Options):
    """Lines 86-150: fixed wage type, fixed kappa offset and young kappa index."""
    Nsim, ninter = cfg.Nsim, cfg.ninter
    ftW = p.ftW
    # wage type: truncated normal on [ftW(1), ftW(nI)], 10-point grid
    Wgrid = np.linspace(ftW[0], ftW[-1], ninter)
    cdfW = _truncnorm_cdf(Wgrid, 0.5 * ftW[0] + 0.5 * ftW[-1], (ftW[-1] - ftW[0]) / 2, ftW[0], ftW[-1])
    drw = rand(cfg.seeds["ftW"], Nsim, 1)[:, 0]
    ftW_i = np.full(Nsim, ftW[-1])
    for j in range(ninter - 2, -1, -1):
        ftW_i[drw < cdfW[j]] = Wgrid[j]
    # fixed kappa offset: truncated normal on [kaplow, 0]
    kapbar = np.array([p.kaplow, 0.0])
    kapbarlow, kapbase = kapbar.min(), kapbar.max()      # NOTE: kapbase overwritten (=0)
    Kgrid = np.linspace(kapbarlow, kapbase, ninter)
    cdfK = _truncnorm_cdf(Kgrid, 0.5 * kapbarlow + 0.5 * kapbase, (kapbase - kapbarlow) / 2, kapbarlow, kapbase)
    drw = rand(cfg.seeds["ftK"], Nsim, 1)[:, 0]
    ftK_i = np.full(Nsim, kapbase)
    for j in range(ninter - 2, -1, -1):
        ftK_i[drw < cdfK[j]] = Kgrid[j]
    # young kappa: uniform over the 4 columns kap(1, 1:4)
    kap = p.kap_matrix(opts)
    Kygrid = kap[0, :4]
    cdfKy = np.cumsum(np.full(4, 0.25))
    seed_ky = cfg.seeds["Ky"] if opts.same_seed_kappa_draws else cfg.seeds["Ky"] + 1
    drw = rand(seed_ky, Nsim, 1)[:, 0]
    Ky_i = np.full(Nsim, Kygrid[3]); Ky_inx = np.full(Nsim, 3)
    for j in range(2, -1, -1):
        m = drw < cdfKy[j]
        Ky_i[m] = Kygrid[j]; Ky_inx[m] = j
    return ftW_i, ftK_i, Ky_i, Ky_inx, kapbarlow, kapbase


def _interp_policies(cfg, p, opts, sol: Solution, ftW_i, ftK_i, Ky_inx, kapbarlow, kapbase):
    """Lines 152-190: per-individual policies by linear interpolation over the
    two fixed-type dimensions (kappa offset, wage type)."""
    nI, nT, nE, nY, nZ = p.nI, p.nT, p.nE, p.nY, p.nZ
    ftW = p.ftW

    def split(g):   # (nI, ...) -> (iw, ik, j, ...)
        return g.reshape(2, 2, 4, nT, nE, nY, nZ)

    gH0, gS0, gQ0 = split(sol.gH), split(sol.gS), split(sol.gQ.astype(float))
    wK = (ftK_i - kapbarlow) / (kapbase - kapbarlow)       # weight on ik=1 (offset 0)
    wW = (ftW_i - ftW[0]) / (ftW[-1] - ftW[0])              # weight on iw=1 (high wage)

    def interp_ok(g0):
        # g1[iw] = (1-wK) g0[iw,0,Ky] + wK g0[iw,1,Ky]   (interp1 over [kapbarlow, kapbase])
        g1 = np.empty((2, cfg.Nsim, nT, nE, nY, nZ))
        for iw in range(2):
            g1[iw] = ((1 - wK)[:, None, None, None, None] * g0[iw, 0, Ky_inx]
                      + wK[:, None, None, None, None] * g0[iw, 1, Ky_inx])
        return (1 - wW)[:, None, None, None, None] * g1[0] + wW[:, None, None, None, None] * g1[1]

    def interp_reversed(g0):
        g1 = np.empty((2, cfg.Nsim, nT, nE, nY, nZ))
        for iw in range(2):
            g1[iw] = (wK[:, None, None, None, None] * g0[iw, 0, Ky_inx]
                      + (1 - wK)[:, None, None, None, None] * g0[iw, 1, Ky_inx])
        return wW[:, None, None, None, None] * g1[0] + (1 - wW)[:, None, None, None, None] * g1[1]

    ggH, ggS = interp_ok(gH0), interp_ok(gS0)
    ggQ = interp_reversed(gQ0) if opts.reversed_type_weights_gQ else interp_ok(gQ0)
    return ggH, ggS, ggQ


def _recession_dates(cfg: SimConfig):
    """Lines 222-268."""
    nq = (cfg.yT - cfg.y0 + 1) * 4
    year = np.repeat(np.arange(cfg.y0, cfg.yT + 1), 4)
    zsim = np.ones(nq, int)
    for strt, endd in [(1957.5, 1958.25), (1960.25, 1961), (1969.75, 1970.75), (1973.75, 1975),
                       (1980, 1980.5), (1981.5, 1982.75), (1990.5, 1991), (2001, 2001.75),
                       (2007.75, 2009.25)]:
        a = int(round((strt - cfg.y0) * 4)); b = int(round((endd - cfg.y0) * 4))
        zsim[a: b + 1] = 2
    return year, zsim


# ----------------------------------------------------------------------------
def simulate(p: Params, sol: Solution, opts: Options | None = None,
             cfg: SimConfig | None = None) -> SimResult:
    opts = opts or Options.faithful()
    cfg = cfg or SimConfig()
    nT, nE, nY, nZ = p.nT, p.nE, p.nY, p.nZ
    egrid = p.egrid
    tsize = p.tsize
    Nsim = cfg.Nsim
    Tsim_i = cfg.Tsim_years * tsize
    Ngen = cfg.yT - cfg.y0
    NgenCareer = cfg.yCareerEnd - cfg.y0
    ScratchT = (cfg.yStart - cfg.y0) * tsize
    Tsim = Ngen * tsize
    Nsim_i = Nsim * Ngen

    lamHbar = markov_stationary(p.lamH[0])
    y1, y2 = cfg.years_per_age
    ageT = np.ones(Tsim_i, int)
    ageT[y1 * tsize: (y1 + y2) * tsize] = 2
    ageT[(y1 + y2) * tsize:] = 3

    ftW_i, ftK_i, Ky_i, Ky_inx, kapbarlow, kapbase = _draw_types(cfg, p, opts)
    ggH, ggS, ggQ = _interp_policies(cfg, p, opts, sol, ftW_i, ftK_i, Ky_inx, kapbarlow, kapbase)
    year, zsim = _recession_dates(cfg)

    # shocks common to every cohort (lines 270-349)
    esim0 = randi(cfg.seeds["exp0"], int(np.floor(0.4 * nE)), int(np.floor(0.6 * nE)), Nsim, 1)[:, 0]
    e0 = egrid[esim0 - 1]
    u_h = rand(cfg.seeds["hstat0"], Nsim, 1)[:, 0]
    hstat0 = np.ones(Nsim, int)
    hstat0[u_h < lamHbar[2] + lamHbar[1]] = 2
    hstat0[u_h < lamHbar[2]] = 3
    u_e = rand(cfg.seeds["emp0"], Nsim, 1)[:, 0]
    unemp0 = (u_e < cfg.Nw_init)
    LH = rand(cfg.seeds["LH"], Nsim, Tsim_i)
    LW = rand(cfg.seeds["LW"], Nsim, Tsim_i)
    QT = rand(cfg.seeds["Q"], Nsim, Tsim_i)

    # panel arrays (one extra column: MATLAB grows arrays when writing tt+1)
    T1 = Tsim + 1
    Z = lambda dt=float: np.zeros((Nsim_i, T1), dt)
    Alive, Age, AgeY = Z(np.int8), Z(np.int8), Z(np.int16)
    expr, Wage = Z(), Z()
    Hstat, Hloss, Wloss = Z(np.int8), Z(np.int8), Z(np.int8)
    EmpI, UnempI, NiLF, NempI = Z(np.int8), Z(np.int8), Z(np.int8), Z(np.int8)
    Quit, Hours, Srch, FT = Z(np.int8), Z(), Z(), Z(np.int8)
    Inc_w, Inc_hh = Z(), Z()
    EE, EU, UE, NE, EN, NU, UN = (Z(np.int8) for _ in range(7))

    # value at birth (lines 360-380): NOTE type indices 1,2 and reversed weights (sic)
    iwl = 2 * (ftW_i - 0.5); iwh = 1 - iwl
    ee0 = np.searchsorted(egrid, e0, side="left") - 1
    BirthV = {}
    for name, (iy, iz) in dict(exp=(0, 0), rec=(0, 1), hunemp=(2, 0), hunemprec=(2, 1)).items():
        BirthV[name] = float(np.mean(iwl * sol.V[0, 0, ee0, iy, iz] + iwh * sol.V[1, 0, ee0, iy, iz]))

    # ----------------------------------------------------------------- helpers
    def e_index(e):
        """MATLAB: ee = find(egrid < e, 1, 'last') (0-based).  Empty -> use top index."""
        k = np.searchsorted(egrid, e, side="left") - 1
        return np.where(k < 0, nE - 1, k)

    def lookup(g, ii, age, ee, e, iy, iz):
        """Policy interpolated in experience (lines 421-424 etc.)."""
        top = ee >= nE - 1
        eel = np.minimum(ee, nE - 2)
        de = egrid[eel + 1] - egrid[eel]
        lo, hi = g[ii, age, eel, iy, iz], g[ii, age, eel + 1, iy, iz]
        w_hi = (e - egrid[eel]) / de
        if opts.reversed_e_weights:
            val = lo * w_hi + hi * (1 - w_hi)
        else:
            val = lo * (1 - w_hi) + hi * w_hi
        return np.where(top, g[ii, age, nE - 1, iy, iz], val)

    # cohorts: in = 0..Ngen-2 (the last generation is tossed), individuals i
    gens = np.repeat(np.arange(Ngen - 1), Nsim)          # cohort of each individual
    ids = np.arange((Ngen - 1) * Nsim)                    # ii
    typ = ids % Nsim                                      # i
    tt0 = gens * 4                                        # birth quarter (0-based)
    ftW_ii = ftW_i[typ]

    if opts.phantom_last_cohort:
        # MATLAB lines 270-311 initialise the birth-quarter entries of EVERY cohort,
        # including the last one that the simulation loop (in = 1:Ngen-1) never visits.
        ii_last = np.arange((Ngen - 1) * Nsim, Ngen * Nsim); tt_last = (Ngen - 1) * 4
        expr[ii_last, tt_last] = e0
        Hstat[ii_last, tt_last] = hstat0
        EmpI[ii_last, tt_last] = (~unemp0).astype(np.int8)
        UnempI[ii_last, tt_last] = unemp0.astype(np.int8)

    # persistent per-individual state variables of the MATLAB loop
    iy = np.zeros(ids.size, int); iz = np.zeros(ids.size, int); ee_last = np.zeros(ids.size, int)

    for it in range(Tsim_i):                              # 0-based life-cycle quarter
        tt = tt0 + it
        alive = tt < Tsim
        ii = ids[alive]; i = typ[alive]; t = tt[alive]
        Alive[ii, t] = 1
        age = ageT[it] - 1                                # 0-based age group
        Age[ii, t] = ageT[it]
        AgeY[ii, t] = cfg.age0 if it == 0 else cfg.age0 + int(np.ceil((it + 1) / 4))
        if it == 0:
            expr[ii, t] = e0[i]
            Hstat[ii, t] = hstat0[i]
            EmpI[ii, t] = 1                                 # line 417 (overrides Nw_init draw)
            UnempI[ii, t] = unemp0[i].astype(np.int8)
            if not opts.all_employed_at_birth:
                EmpI[ii, t] = (~unemp0[i]).astype(np.int8)
            iy_c = hstat0[i] - 1
            iz_c = zsim[t] - 1
            iy[alive] = iy_c; iz[alive] = iz_c
        e = expr[ii, t]
        Wage[ii, t] = F.wage(p, ftW_ii[alive], e) * (1.0 if opts.no_bcwage_in_sim else p.BCwage[zsim[t] - 1])
        if it > 0:
            # husband transition drawn with LAST period's (iy, iz)  (lines 479-488)
            iy_p, iz_p = iy[alive], iz[alive]
            z_for_h = iz_p if opts.lagged_z_in_sim else zsim[t] - 1
            u = LH[i, it]
            pR = p.lamH[z_for_h, iy_p, 1]; pU = p.lamH[z_for_h, iy_p, 2]
            hs = np.ones(ii.size, int)
            hs[u < pR + pU] = 2
            hs[u < pU] = 3
            Hstat[ii, t] = hs
            Hloss[ii, t] = ((hs == 3) & (Hstat[ii, t - 1] < 3)).astype(np.int8)
            iy_c = hs - 1; iz_c = zsim[t] - 1
            iy[alive] = iy_c; iz[alive] = iz_c
        emp_in = EmpI[ii, t] > 0
        # quit / job loss decision for those entering employed (lines 418-428, 492-506)
        ee_new = e_index(e)
        ee_use = ee_last[alive].copy()
        ee_use[emp_in] = ee_new[emp_in]                   # ee refreshed only when employed
        if it == 0:
            ee_use = ee_new                                # birth block always computes ee
        ee_last[alive] = ee_use
        quit = np.zeros(ii.size, bool); loss = np.zeros(ii.size, bool)
        Q = lookup(ggQ, i, age, ee_use, e, iy_c, iz_c)
        quit[emp_in] = Q[emp_in] > QT[i, it][emp_in]
        if it > 0:
            lw = LW[i, it] < p.lossW[iz_c]
            if opts.quit_before_loss:
                loss[emp_in & ~quit] = lw[emp_in & ~quit]
            else:
                loss[emp_in] = lw[emp_in]
                quit[emp_in & loss] = False
        Quit[ii, t] = quit
        Wloss[ii, t] = loss
        if it == 0:
            EmpI[ii[quit], t[quit]] = 0
        work = emp_in & ~quit & ~loss
        # ---- employed this period --------------------------------------------
        w = work
        if it > 0:
            EE[ii[w], t[w]] = EmpI[ii[w], t[w] - 1]
            UE[ii[w], t[w]] = UnempI[ii[w], t[w] - 1]
            NE[ii[w], t[w]] = NiLF[ii[w], t[w] - 1]
        EmpI[ii[w], t[w] + 1] = 1
        hrs = lookup(ggH, i[w], age, ee_use[w], e[w], iy_c[w], iz_c[w])
        Hours[ii[w], t[w]] = hrs
        FT[ii[w], t[w]] = (hrs > cfg.FTbar)
        Inc_w[ii[w], t[w]] = hrs * Wage[ii[w], t[w]]
        yh = p.wageH[age] * p.BCwage[iz_c] * p.ym[iy_c]
        Inc_hh[ii[w], t[w]] = Inc_w[ii[w], t[w]] + yh[w]
        expr[ii[w], t[w] + 1] = F.exp2(p, e[w], hrs)
        # ---- non-employed this period ------------------------------------------
        n = ~work
        NempI[ii[n], t[n]] = 1
        Hours[ii[n], t[n]] = 0.0
        expr[ii[n], t[n] + 1] = F.exp2(p, e[n], 0.0)
        srch = lookup(ggS, i[n], age, ee_use[n], e[n], iy_c[n], iz_c[n])   # stale ee when non-employed
        Srch[ii[n], t[n]] = srch
        isU = srch > cfg.NiLFbar
        iiN, tN = ii[n], t[n]
        UnempI[iiN[isU], tN[isU]] = 1
        NiLF[iiN[~isU], tN[~isU]] = 1
        if it > 0:
            EU[iiN[isU], tN[isU]] = EmpI[iiN[isU], tN[isU]]
            NU[iiN[isU], tN[isU]] = NiLF[iiN[isU], tN[isU] - 1]
            EN[iiN[~isU], tN[~isU]] = EmpI[iiN[~isU], tN[~isU]]
            UN[iiN[~isU], tN[~isU]] = UnempI[iiN[~isU], tN[~isU] - 1]
            EmpI[iiN, tN] = 0
        Inc_hh[iiN, tN] = yh[n]
        # Job finding (line 547).  Extrapolation with a stale experience index can make
        # the interpolated search intensity negative; MATLAB then raises a negative
        # number to the power nu, obtains a complex number, and `>` compares real parts.
        s_pow = np.where(srch >= 0, np.abs(srch) ** p.nu,
                         np.abs(srch) ** p.nu * np.cos(p.nu * np.pi))
        found = s_pow * p.findW[iz_c[n]] > LW[i[n], it]
        EmpI[iiN[found], tN[found] + 1] = 1

    return SimResult(cfg=cfg, params=p, options=opts, ftW_i=ftW_i, ftK_i=ftK_i, Ky_i=Ky_i,
                     Ky_inx=Ky_inx, year=year, zsim=zsim, ageT=ageT, Ngen=Ngen, Tsim=Tsim,
                     Tsim_i=Tsim_i, Nsim_i=Nsim_i, ScratchT=ScratchT, NgenCareer=NgenCareer,
                     Alive=Alive, Age=Age, AgeY=AgeY, exp=expr, Wage=Wage, Hstat=Hstat,
                     Hloss=Hloss, Wloss=Wloss, EmpI=EmpI, UnempI=UnempI, NiLF=NiLF, NempI=NempI,
                     Quit=Quit, Hours=Hours, Srch=Srch, FT=FT, Inc_w=Inc_w, Inc_hh=Inc_hh,
                     EE=EE, EU=EU, UE=UE, NE=NE, EN=EN, NU=NU, UN=UN, BirthV=BirthV,
                     ggH=ggH, ggS=ggS, ggQ=ggQ)
