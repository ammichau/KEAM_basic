"""Simulated method of moments for the 1940s cohort.

Targets (slides p.30, 35, 37; unemployment rate assumed) and the parameters that move them:
  employment rate, hours, career shares, monthly quit and E->nonE rates by cycle, employment
  drop in recessions, within-couple wage gap, unemployment rate, wife's income share.
Parameters: mu, kbar_max, km_max, tau_w, lam_f (expansion; recession = 0.85x),
  lam_u (expansion and recession), ybar_h, s_bar, sd_kT (transitory cost shock).
"""
from __future__ import annotations
import json, os, time
import numpy as np
from scipy.optimize import minimize
from .params import FinalParams
from .solve import solve_all
from .simulate import simulate_final, SimConfigFinal
from .moments import moments_final

TARGETS = {
    "E/pop": 0.62, "hours|E": 0.40,
    "share Lifecycle": 0.31, "share PT": 0.28, "share Career": 0.19, "share NiLF": 0.22,
    "quit/m exp": 0.034, "quit/m rec": 0.028, "E->nonE/m exp": 0.050, "E->nonE/m rec": 0.048,
    "dE/pop rec-exp (pts)": -1.7, "wage gap (FTE earnings ratio)": 0.71, "U rate": 0.05,
    "wife share exp": 0.225,
}
# scale for each target (deviation divided by this); percentage-point moments use absolute scales
SCALE = {k: v for k, v in TARGETS.items()}
SCALE["dE/pop rec-exp (pts)"] = 1.0
WEIGHT = {k: 1.0 for k in TARGETS}
WEIGHT.update({"E/pop": 3.0, "hours|E": 2.0, "quit/m exp": 2.0, "quit/m rec": 2.0, "wage gap (FTE earnings ratio)": 2.0})

PARAM_NAMES = ["mu", "kbar_max", "km_max", "tau_w", "lam_f0", "lam_u0", "lam_u1", "ybar_h", "s_bar", "sd_kT"]
BOUNDS = {"mu": (0.2, 5.0), "kbar_max": (0.005, 0.6), "km_max": (1.0, 6.0), "tau_w": (0.4, 1.2),
          "lam_f0": (0.05, 0.9), "lam_u0": (0.003, 0.05), "lam_u1": (0.003, 0.08), "ybar_h": (0.0, 0.6),
          "s_bar": (0.02, 0.9), "sd_kT": (0.001, 0.6)}


def apply_params(base: FinalParams, x: dict) -> FinalParams:
    kw = dict(mu=x["mu"], kbar_max=x["kbar_max"], km_max=x["km_max"], tau_w=x["tau_w"],
              lam_f=(x["lam_f0"], 0.85 * x["lam_f0"]), lam_u=(x["lam_u0"], x["lam_u1"]),
              ybar_h=x["ybar_h"], s_bar=x["s_bar"], sd_kT=x["sd_kT"])
    return base.replace(**kw)


def _to_unit(x):   # map params to R via logit of the bounded interval
    z = []
    for n in PARAM_NAMES:
        lo, hi = BOUNDS[n]; u = (x[n] - lo) / (hi - lo); u = min(max(u, 1e-6), 1 - 1e-6)
        z.append(np.log(u / (1 - u)))
    return np.array(z)


def _from_unit(z):
    x = {}
    for n, v in zip(PARAM_NAMES, z):
        lo, hi = BOUNDS[n]; x[n] = lo + (hi - lo) / (1 + np.exp(-v))
    return x


def objective_from_moments(m):
    tot = 0.0; parts = {}
    for k, tv in TARGETS.items():
        mv = m.get(k, np.nan)
        if not np.isfinite(mv):
            mv = 0.0
        d = (mv - tv) / SCALE[k]
        parts[k] = d
        tot += WEIGHT[k] * d * d
    return tot, parts


def evaluate(base: FinalParams, x: dict, cfg: SimConfigFinal, n_jobs=4):
    p = apply_params(base, x)
    sol = solve_all(p, n_jobs=n_jobs)
    sim = simulate_final(p, sol, cfg)
    m = moments_final(sim)
    obj, parts = objective_from_moments(m)
    return obj, m, parts, p


def run_smm(base: FinalParams, x0: dict, cfg: SimConfigFinal, log_path: str, maxfev=150, n_jobs=4):
    hist = []
    t0 = time.time()

    def fun(z):
        x = _from_unit(z)
        obj, m, parts, _ = evaluate(base, x, cfg, n_jobs)
        hist.append(dict(obj=obj, x=x, m={k: float(v) for k, v in m.items() if isinstance(v, (float, int))}))
        with open(log_path, "a") as fh:
            fh.write(json.dumps(dict(n=len(hist), t=round(time.time() - t0), obj=obj, x=x,
                                     dev={k: round(v, 3) for k, v in parts.items()})) + "\n")
        return obj

    res = minimize(fun, _to_unit(x0), method="Nelder-Mead",
                   options=dict(maxfev=maxfev, xatol=1e-3, fatol=1e-4, initial_simplex=None))
    best = min(hist, key=lambda h: h["obj"])
    return best, hist, res
