"""Monthly life-cycle simulation of the final model.

Individuals are stored in life-cycle time (month since entry at age 25, 0..L-1) so that
arrays are (N_ind, L).  Calendar month = entry + life month.  Annual entry cohorts,
independent type and shock draws per cohort, bilinear policy interpolation in (e, a).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from .params import FinalParams
from .solve import FinalSolution


NBER = [("1957-08", "1958-04"), ("1960-04", "1961-02"), ("1969-12", "1970-11"), ("1973-11", "1975-03"),
        ("1980-01", "1980-07"), ("1981-07", "1982-11"), ("1990-07", "1991-03"), ("2001-03", "2001-11"),
        ("2007-12", "2009-06")]


def nber_path(y0=1955, y1=2019):
    """Monthly recession indicator (0 expansion, 1 recession), months after the peak up to the trough."""
    T = (y1 - y0 + 1) * 12
    z = np.zeros(T, int)
    for a, b in NBER:
        ya, ma = map(int, a.split("-")); yb, mb = map(int, b.split("-"))
        ia = (ya - y0) * 12 + (ma - 1) + 1; ib = (yb - y0) * 12 + (mb - 1)
        z[max(ia, 0): min(ib, T - 1) + 1] = 1
    return z


def markov_path(piz, T, seed=0, burn=240):
    rng = np.random.default_rng(seed)
    z = np.zeros(T + burn, int); u = rng.random(T + burn)
    for t in range(1, T + burn):
        z[t] = 1 if u[t] < piz[z[t - 1], 1] else 0
    return z[burn:]


def stationary(P):
    w, v = np.linalg.eig(P.T)
    x = np.real(v[:, np.argmin(np.abs(w - 1))]); return x / x.sum()


@dataclass
class SimConfigFinal:
    N: int = 60                 # women per annual cohort
    n_cohorts: int = 90         # annual entry cohorts
    L: int = 480                # months of working life (25-64)
    zmode: str = "markov"       # 'markov' or 'nber'
    seed: int = 12345
    z_seed: int = 7
    window: tuple = None        # calendar months used for moments (default: fully populated part)


@dataclass
class FinalSim:
    cfg: SimConfigFinal; params: FinalParams
    entry: np.ndarray; ktype: np.ndarray; zpath: np.ndarray; T: int
    emp: np.ndarray; stat: np.ndarray; quit: np.ndarray; loss: np.ndarray; declined: np.ndarray
    hours: np.ndarray; wage: np.ndarray; e: np.ndarray; a: np.ndarray; srch: np.ndarray
    hstat: np.ndarray; hloss: np.ndarray; inc_w: np.ndarray; inc_h: np.ndarray; cons: np.ndarray
    age: np.ndarray                                   # age group per life month (L,)

    @property
    def calendar(self):
        return self.entry[:, None] + np.arange(self.cfg.L)[None, :]


def _interp2(g, k, tau, ek, ew, ak, aw, y, z, j=0):
    """Bilinear interpolation of g[k, tau, e, a, y, z, j] at (e, a) brackets (j: cost-shock state)."""
    g00 = g[k, tau, ek, ak, y, z, j]; g01 = g[k, tau, ek, ak + 1, y, z, j]
    g10 = g[k, tau, ek + 1, ak, y, z, j]; g11 = g[k, tau, ek + 1, ak + 1, y, z, j]
    return ((1 - ew) * ((1 - aw) * g00 + aw * g01) + ew * ((1 - aw) * g10 + aw * g11))


def simulate_final(p: FinalParams, sol: FinalSolution, cfg: SimConfigFinal | None = None) -> FinalSim:
    cfg = cfg or SimConfigFinal()
    N, nC, L = cfg.N, cfg.n_cohorts, cfg.L
    if cfg.zmode == "nber":
        zpath = nber_path(); nC = min(nC, 65)
    else:
        zpath = markov_path(p.piz, 12 * nC + L, cfg.z_seed)
    T = zpath.size
    rng = np.random.default_rng(cfg.seed)
    Nind = N * nC
    entry = np.repeat(np.arange(nC) * 12, N)
    nK = sol.omega.size
    ktype = rng.integers(0, nK, Nind)
    eg, ag = p.egrid, p.agrid
    age = np.zeros(L, int); m0, m1, m2 = p.age_months
    age[m0: m0 + m1] = 1; age[m0 + m1:] = 2
    lamH = p.lamH(); yH = p.y_husband()
    lam_u = np.asarray(p.lam_u); lam_f = np.asarray(p.lam_f)
    phi = np.array([1.0, p.phi_rec])
    omega = sol.omega[ktype]
    f = p.ybar_h + p.z_h * omega ** p.alpha_h
    # draws
    uH = rng.random((Nind, L), dtype=np.float32); uJ = rng.random((Nind, L), dtype=np.float32)
    kT_nodes, kT_w = p.kT_nodes()
    kT_idx = rng.choice(kT_nodes.size, size=(Nind, L), p=kT_w)
    if p.rho_kT > 0:                      # persistent shock: keep last month's node with probability rho_kT
        keep = rng.random((Nind, L), dtype=np.float32) < p.rho_kT
        for it in range(1, L):
            kT_idx[:, it] = np.where(keep[:, it], kT_idx[:, it - 1], kT_idx[:, it])
    kT_draw = kT_nodes[kT_idx].astype(np.float32)
    jst = kT_idx if sol.VE.shape[-1] > 1 else np.zeros((Nind, L), int)     # shock state used in the policies
    # storage (life-cycle time)
    I8 = lambda: np.zeros((Nind, L), np.int8); F4 = lambda: np.zeros((Nind, L), np.float32)
    emp, stat, quit, loss, declined, hstat, hloss = I8(), I8(), I8(), I8(), I8(), I8(), I8()
    hours, wage, E_, A_, srch, inc_w, inc_h, cons = (F4() for _ in range(8))
    # initial state
    e = rng.uniform(p.e0_range[0], p.e0_range[1], Nind)
    a = np.zeros(Nind)
    y = rng.choice(3, size=Nind, p=stationary(lamH[0]))
    E_in = np.ones(Nind, bool)
    lost_prev = np.zeros(Nind, bool); found_prev = np.zeros(Nind, bool)
    ids = np.arange(Nind)
    for it in range(L):
        t = entry + it
        valid = t < T
        z = np.where(valid, zpath[np.minimum(t, T - 1)], 0)
        tau = age[it]
        ek = np.clip(np.searchsorted(eg, e, side="right") - 1, 0, eg.size - 2)
        ew = np.clip((e - eg[ek]) / (eg[ek + 1] - eg[ek]), 0, 1)
        ak = np.clip(np.searchsorted(ag, a, side="right") - 1, 0, ag.size - 2)
        aw = np.clip((a - ag[ak]) / (ag[ak + 1] - ag[ak]), 0, 1)
        VEi = _interp2(sol.VE, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        VNi = _interp2(sol.VN, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        q = E_in & (VNi > VEi - kT_draw[:, it])
        work = E_in & ~q
        w = phi[z] * p.tau_w * omega * (1 + p.gam_e * e ** p.xi)
        yh = yH[tau, y, z]
        # employed
        h = _interp2(sol.gH, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        aE = _interp2(sol.gAE, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        # non-employed
        s = _interp2(sol.gS, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        aN = _interp2(sol.gAN, ktype, tau, ek, ew, ak, aw, y, z, jst[:, it])
        f_t = f * (p.home_young_mult if tau == 0 else 1.0)
        incE = w * h + f_t * (1 - h) ** p.nu_h + yh
        incN = f_t * (1 - s) ** p.nu_h + yh
        aE = np.minimum(aE, a + incE - 1e-6); aN = np.minimum(aN, a + incN - 1e-6)
        a_next = np.where(work, aE, aN)
        c = np.where(work, incE + a - aE, incN + a - aN)
        e_next = np.where(work, np.minimum(p.e_max, (1 - p.delta_e) * e + p.theta_e * e * h ** p.psi_e),
                          (1 - p.delta_e) * e)
        lost = work & (uJ[:, it] < lam_u[z])
        found = ~work & (uJ[:, it] < lam_f[z] * s ** p.nu)
        # records
        emp[:, it] = work; quit[:, it] = q & ~found_prev; declined[:, it] = q & found_prev
        loss[:, it] = lost_prev
        stat[:, it] = np.where(work, 0, np.where(s >= p.s_bar, 1, 2))
        hours[:, it] = np.where(work, h, 0); wage[:, it] = w; E_[:, it] = e; A_[:, it] = a
        srch[:, it] = np.where(work, 0, s); hstat[:, it] = y
        inc_w[:, it] = np.where(work, w * h, 0); inc_h[:, it] = yh; cons[:, it] = c
        # transitions
        y_new = np.ones(Nind, int)
        pU = lamH[z, y, 2]; pR = lamH[z, y, 1]
        y_new[uH[:, it] < pR + pU] = 1
        y_new[uH[:, it] < pU] = 2
        if it + 1 < L:
            hloss[:, it + 1] = (y_new == 2) & (y < 2)
        y = y_new
        e, a = e_next, a_next
        E_in = (work & ~lost) | found
        lost_prev, found_prev = lost, found
    return FinalSim(cfg=cfg, params=p, entry=entry, ktype=ktype, zpath=zpath, T=T, emp=emp, stat=stat,
                    quit=quit, loss=loss, declined=declined, hours=hours, wage=wage, e=E_, a=A_, srch=srch,
                    hstat=hstat, hloss=hloss, inc_w=inc_w, inc_h=inc_h, cons=cons, age=age)
