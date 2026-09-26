"""Experiments for the final model.

* single-factor experiments sized to reproduce the 1970s cohort employment rate:
    - returns to experience (gam_e)
    - compensated wage gap: tau_w up, husband income scaled down so that total household income
      at BASELINE behaviour is unchanged
    - cost of work (kbar_max and km_max scaled together)
* cohort accounting with the observed tau_w and gam_e paths and the residual cost scale
* mechanism counterfactuals: acyclical husband job-loss risk, acyclical job finding, no wage cut
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from .params import FinalParams
from .solve import solve_all
from .simulate import simulate_final, SimConfigFinal
from .moments import moments_final


def run(p: FinalParams, cfg: SimConfigFinal, n_jobs=None):
    sol = solve_all(p, n_jobs=n_jobs)
    sim = simulate_final(p, sol, cfg)
    return moments_final(sim), sol, sim


def compensated_wage_gap(p: FinalParams, base_m: dict, scale: float) -> FinalParams:
    """Raise tau_w by `scale` and lower the husband's income so that, at baseline behaviour,
    expected household income is unchanged: yH_scale = 1 - share_w (scale - 1) / (1 - share_w)."""
    sw = base_m["wife share exp"]
    yH_scale = 1.0 - sw * (scale - 1.0) / (1.0 - sw)
    return p.replace(tau_w=p.tau_w * scale, yH_scale=p.yH_scale * yH_scale)


def cost_scaled(p: FinalParams, scale: float) -> FinalParams:
    return p.replace(kbar_max=p.kbar_max * scale)


def returns_scaled(p: FinalParams, scale: float) -> FinalParams:
    return p.replace(gam_e=p.gam_e * scale)


def size_to_employment(make, p, cfg, target_E, lo, hi, tol=0.002, maxit=12, n_jobs=None):
    """Find the scale such that E/pop hits target_E (bisection on a common simulation seed;
    12 steps resolve the scale to 2^-12 of the bracket, tolerance 0.2 pp of employment)."""
    cache = {}

    def f(s):
        if s not in cache:
            m, _, _ = run(make(p, s), cfg, n_jobs)
            cache[s] = m
        return cache[s]["E/pop"] - target_E

    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        s = lo if abs(flo) < abs(fhi) else hi
        return s, cache[s], cache
    for _ in range(maxit):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if abs(fm) < tol:
            return mid, cache[mid], cache
        if fm * flo < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    mid = 0.5 * (lo + hi); f(mid)
    return mid, cache[mid], cache


def acyclical_husband(p: FinalParams) -> FinalParams:
    """Husband's job-loss and job-finding rates and his unemployment income at their expansion values in
    both aggregate states (the precautionary channel switched off)."""
    return p.replace(lamH_loss=(p.lamH_loss[0], p.lamH_loss[0]), lamH_find=(p.lamH_find[0], p.lamH_find[0]),
                     ui_rec_mult=1.0)


def acyclical_finding(p: FinalParams) -> FinalParams:
    """The wife's job-finding efficiency at its expansion value in both states (job hoarding switched off)."""
    return p.replace(lam_f=(p.lam_f[0], p.lam_f[0]))


def counterfactuals(p: FinalParams, cfg: SimConfigFinal, n_jobs=None):
    out = {}
    out["baseline"] = run(p, cfg, n_jobs)[0]
    # acyclical husband job-loss risk (expansion values in both states)
    out["acyclical husband risk"] = run(acyclical_husband(p), cfg, n_jobs)[0]
    out["acyclical job finding"] = run(acyclical_finding(p), cfg, n_jobs)[0]
    out["no recession wage cut"] = run(p.replace(phi_rec=1.0, phi_rec_H=1.0), cfg, n_jobs)[0]
    out["acyclical own job loss"] = run(p.replace(lam_u=(p.lam_u[0], p.lam_u[0])), cfg, n_jobs)[0]
    return out
