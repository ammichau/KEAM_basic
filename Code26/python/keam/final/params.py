from __future__ import annotations
from dataclasses import dataclass, field, replace
import numpy as np
from scipy.stats import norm


@dataclass
class FinalParams:
    # ---------------- preferences (monthly) ----------------
    beta: float = 0.99          # monthly (0.886 per year); the calibration since 2026-09-29 uses 0.993 (0.919 per
                                # year) with r_a 0.00327 and the fine savings choice grid nAc 100 (see below)
    r_a: float = 0.0            # monthly net return on assets (0 = the original model; 4% per year = 0.00327)
    gamma: float = 2.0          # CRRA
    eta: float = 1.4            # curvature of hours disutility (slides p.29)
    kpr: bool = False           # King-Plosser-Rebelo preferences: u = U(log c - mu h^(1+eta)/(1+eta) - kappa)
                                # with U(x) = exp((1-gamma) x)/(1-gamma); balanced-growth consistent for any
                                # gamma, nests the separable log model at gamma = 1. The fixed costs kappa
                                # (kbar, km) enter the aggregator as exp(-kappa), a proportional consumption
                                # tax on employment. The transitory shock is additive to V^E at the quit
                                # decision unless kT_mult (then also a proportional tax, scale-invariant).
    kT_mult: bool = False       # with kpr: the transitory cost shock is a one-month proportional consumption
                                # tax on employment, u = U(log c - v(h) - kappa - kappa_T) (scale-invariant),
                                # instead of an additive shock to V^E at the quit decision. The values then carry
                                # the shock node (iid rows of the chain); nests the additive model at gamma = 1.
    ez_rra: float = 1.0         # Epstein-Zin risk aversion with unit intertemporal elasticity (use with
                                # gamma = 1): continuation values enter through the certainty equivalent
                                # -(1/theta) log E exp(-theta V) with theta = (ez_rra - 1)(1 - beta), the
                                # log form of V = c^(1-beta) [E V'^(1-rra)]^(beta/(1-rra)). 1 = expected utility.
    mu: float = 1.0             # weight on hours disutility (calibrated)
    # ---------------- wage: phi(Z) tau_w omega (1 + gam_e e^xi) ----------------
    tau_w: float = 0.71         # gender wage penalty (calibrated to the within-couple gap)
    gam_e: float = 0.5          # returns to experience (slides p.28)
    xi: float = 0.8
    phi_rec: float = 0.88       # wage penalty in recession (slides p.28)
    phi_rec_H: float = 0.88     # husband's income penalty in recession (assumption: same as wife's)
    # ---------------- experience ----------------
    delta_e: float = 0.005      # monthly depreciation (slides p.28)
    theta_e: float = 0.025
    psi_e: float = 0.66
    e_max: float = 2.0
    e0_range: tuple = (0.2, 1.0)  # initial experience: 1/10 to 1/2 of the maximum (paper p.27)
    # ---------------- home production ----------------
    ybar_h: float = 0.11
    z_h: float = 0.45
    alpha_h: float = 0.21
    nu_h: float = 0.65
    home_young_mult: float = 1.0        # multiplier on f(omega) at ages 25-39 (child care); 1 = off
    kappa_h_power: float = 0.0          # fixed cost scales with (h/0.4)^p; 0 = pure fixed cost (paper)
    # ---------------- search / job loss (monthly) ----------------
    nu: float = 0.5             # search efficiency curvature (slides p.29)
    lam_f: tuple = (0.40, 0.34)   # job-finding efficiency, expansion / recession (15% lower)
    lam_u: tuple = (0.017, 0.029)  # exogenous job loss, expansion / recession (1.7x)
    s_bar: float = 0.25          # search intensity above which a non-employed woman counts as unemployed
    # ---------------- husband (monthly) ----------------
    lamH_loss: tuple = (0.0135, 0.024)   # E -> U
    lamH_find: tuple = (0.35, 0.28)      # U -> R (re-employed with scar)
    scar_end: float = 1.0 / 40.0         # R -> E (mean scar 3.3 years)
    scar_loss_mult: float = 2.5          # R -> U relative to E -> U
    ym_share: tuple = (1.0, 0.85, 0.30)  # income share by husband state E / R / U (30% replacement)
    ui_rec_mult: float = 1.0             # multiplier on the U-state share in recessions (0.5: replacement 15%)
    yH_age: tuple = (0.89, 1.0, 0.94)    # husband income by wife's age group (paper Table 2)
    yH_scale: float = 1.0                # scaled in the compensated wage-gap experiment
    # ---------------- ageing / retirement / aggregate (monthly) ----------------
    age_months: tuple = (180, 180, 120)  # 25-39, 40-54, 55-64
    pension: float = 0.5                 # retirement income relative to middle-age husband income
    death: float = 1.0 / 240.0           # monthly death hazard in retirement (20 years)
    piz: np.ndarray = field(default_factory=lambda: np.array([[0.985, 0.015], [0.09, 0.91]]))
    # ---------------- types ----------------
    n_omega: int = 5
    sd_log_omega: float = 0.37           # sd of the wage fixed effect (slides p.30)
    n_kbar: int = 5
    kbar_max: float = 0.075              # permanent cost, truncated normal on [0, kbar_max]
    n_km: int = 4
    km_max: float = 2.27                 # life-cycle multiplier, uniform on [1, km_max]
    n_zh: int = 1                        # permanent home-productivity type: equiprobable multipliers on z_h
    zh_spread: float = 0.0               # (1 - s, 1 + s) for two points, (1 - s, 1, 1 + s) for three; 0 = off
    n_kT: int = 5                        # transitory cost-of-work shock (slides p.13), discrete normal
    sd_kT: float = 0.05
    rho_kT: float = 0.0                  # monthly probability that the shock keeps its value (0 = iid)
    # ---------------- grids ----------------
    nE: int = 20
    nA: int = 20
    a_max: float = 15.0
    liq: float = 1.0                     # share of the asset stock that can be drawn down within a month:
                                         # a' >= (1 - liq) a (1 = fully liquid, the baseline); retirement unrestricted
    nAc: int = 0                         # savings choice grid: a' chosen on a finer quadratic grid of nAc points on
                                         # [0, a_max] with the continuation value interpolated linearly in a'
                                         # (0 = the state grid, the original solver). The discrete choice suppresses
                                         # saving when beta (1 + r) is near one: a' can only move in grid-spacing lumps.
    nH: int = 20                         # hours grid on [h_min, 1]
    h_min: float = 0.05
    nS: int = 21                         # search grid on [0, 1]
    # ---------------- solver ----------------
    vf_tol: float = 1e-5
    howard_steps: int = 25
    max_iter: int = 400

    def __post_init__(self):
        self.piz = np.asarray(self.piz, float)

    # ------------------------------------------------------------ grids
    @property
    def egrid(self):
        return np.linspace(0.0, self.e_max, self.nE)

    @property
    def agrid(self):
        return self.a_max * np.linspace(0.0, 1.0, self.nA) ** 2

    @property
    def agrid_c(self):
        """Savings choice grid (the state grid unless nAc exceeds nA)."""
        if self.nAc > self.nA:
            return self.a_max * np.linspace(0.0, 1.0, self.nAc) ** 2
        return self.agrid

    @property
    def hgrid(self):
        return np.linspace(self.h_min, 1.0, self.nH)

    @property
    def sgrid(self):
        return np.linspace(0.0, 1.0, self.nS)

    @property
    def p_age(self):
        return np.array([1.0 / m for m in self.age_months])

    def kT_nodes(self):
        """Equiprobable discretisation of N(0, sd_kT^2)."""
        if self.n_kT <= 1 or self.sd_kT <= 0:
            return np.array([0.0]), np.array([1.0])
        q = (np.arange(self.n_kT) + 0.5) / self.n_kT
        nodes = norm.ppf(q) * self.sd_kT
        nodes = nodes - nodes.mean()
        return nodes, np.full(self.n_kT, 1.0 / self.n_kT)

    def kT_chain(self):
        """Cost-shock nodes, the number of shock states carried in the value functions (1 when the
        shock is iid) and the transition matrix from state to next period's node, (n_state, n_kT):
        with probability rho_kT the shock keeps its value, otherwise it is redrawn from the nodes."""
        nodes, w = self.kT_nodes()
        if nodes.size > 1 and self.kT_mult and self.rho_kT <= 0:    # the flow depends on today's node
            return nodes, nodes.size, np.repeat(w[None, :], nodes.size, axis=0)
        if self.rho_kT <= 0 or nodes.size == 1:
            return nodes, 1, w[None, :]
        P = self.rho_kT * np.eye(nodes.size) + (1 - self.rho_kT) * w[None, :]
        return nodes, nodes.size, P

    def lamH(self):
        """Husband transition matrices by aggregate state: (nZ, 3, 3), rows = from (E, R, U)."""
        out = np.zeros((2, 3, 3))
        for z in range(2):
            lu = self.lamH_loss[z]; lf = self.lamH_find[z]
            out[z, 0] = [1 - lu, 0.0, lu]
            out[z, 1] = [self.scar_end, 1 - self.scar_end - self.scar_loss_mult * lu, self.scar_loss_mult * lu]
            out[z, 2] = [0.0, lf, 1 - lf]
        return out

    def y_husband(self):
        """Husband income by (age group, state, Z): (3, 3, 2)."""
        phi = np.array([1.0, self.phi_rec_H])
        share = np.tile(np.asarray(self.ym_share, float)[:, None], (1, 2))          # (state, Z)
        share[2, 1] *= self.ui_rec_mult                                             # UI replacement in recessions
        return (self.yH_scale * np.asarray(self.yH_age)[:, None, None]
                * share[None, :, :] * phi[None, None, :])

    def replace(self, **kw):
        return replace(self, **kw)


def make_types(p: FinalParams):
    """Discrete type grid: omega (log-normal, mean 1), kbar (truncated normal), km (uniform).
    Returns arrays omega_k, kbar_k, km_k of length n_omega*n_kbar*n_km and equal weights."""
    q = (np.arange(p.n_omega) + 0.5) / p.n_omega
    lo = norm.ppf(q) * p.sd_log_omega
    omega = np.exp(lo); omega = omega / np.mean(omega)
    # truncated normal on [0, kbar_max], mean at the midpoint, sd = half the range (paper p.22)
    q = (np.arange(p.n_kbar) + 0.5) / p.n_kbar
    mu_k, sd_k = p.kbar_max / 2, p.kbar_max / 2
    if sd_k <= 0:
        kbar = np.zeros(p.n_kbar)
    else:
        a, b = norm.cdf((0 - mu_k) / sd_k), norm.cdf((p.kbar_max - mu_k) / sd_k)
        kbar = mu_k + sd_k * norm.ppf(a + q * (b - a))
    km = np.linspace(1.0, p.km_max, p.n_km) if p.n_km > 1 else np.array([1.0])
    zh = 1.0 + p.zh_spread * np.linspace(-1.0, 1.0, p.n_zh) if p.n_zh > 1 else np.array([1.0])
    O, K, M, Z = np.meshgrid(omega, kbar, km, zh, indexing="ij")
    return O.ravel(), K.ravel(), M.ravel()


def make_types4(p: FinalParams):
    """As make_types, plus the home-productivity multiplier of each type."""
    q = (np.arange(p.n_omega) + 0.5) / p.n_omega
    lo = norm.ppf(q) * p.sd_log_omega
    omega = np.exp(lo); omega = omega / np.mean(omega)
    q = (np.arange(p.n_kbar) + 0.5) / p.n_kbar
    mu_k, sd_k = p.kbar_max / 2, p.kbar_max / 2
    if sd_k <= 0:
        kbar = np.zeros(p.n_kbar)
    else:
        a, b = norm.cdf((0 - mu_k) / sd_k), norm.cdf((p.kbar_max - mu_k) / sd_k)
        kbar = mu_k + sd_k * norm.ppf(a + q * (b - a))
    km = np.linspace(1.0, p.km_max, p.n_km) if p.n_km > 1 else np.array([1.0])
    zh = 1.0 + p.zh_spread * np.linspace(-1.0, 1.0, p.n_zh) if p.n_zh > 1 else np.array([1.0])
    O, K, M, Z = np.meshgrid(omega, kbar, km, zh, indexing="ij")
    return O.ravel(), K.ravel(), M.ravel(), Z.ravel()
