"""Simulated method of moments for the 1940s cohort.

Targets (slides p.30, 35, 37; unemployment rate assumed) and the parameters that move them:
  employment rate, hours, career shares, monthly quit and E->nonE rates by cycle, employment
  drop in recessions, within-couple wage gap, unemployment rate, wife's income share.
Parameters: mu, kbar_max, km_max, tau_w, lam_f (expansion; recession = 0.85x),
  lam_u (expansion and recession), ybar_h, sd_kT (transitory cost shock). s_bar (unemployment
  definition) is fixed at 0.25; the unemployment rate and the wife's income share are reported, not targeted.
Optionally the experience depreciation delta_e (`OPTIONAL_PARAMS`; `names=PARAM_NAMES + ["delta_e"]`).

Bounds: km_max up to 15, kbar_max up to 1.0. With km_max <= 6 the life-cycle career share cannot
exceed about 9% (scripts/explore_lifecycle.py, output/explore_lifecycle.out); it keeps rising beyond
10 (scripts/explore_experience.py, output/explore_experience.out).
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
    "dE/pop rec-exp (pts)": -1.7, "wage gap (hourly ratio)": 0.71,
}
# NOT targeted (model outputs on slides p.36-37): unemployment rate, wife's share of income.
# scale for each target (deviation divided by this); percentage-point moments use absolute scales
SCALE = {k: v for k, v in TARGETS.items()}
SCALE["dE/pop rec-exp (pts)"] = 1.0
WEIGHT = {k: 1.0 for k in TARGETS}
WEIGHT.update({"E/pop": 3.0, "hours|E": 2.0, "quit/m exp": 2.0, "quit/m rec": 2.0, "wage gap (hourly ratio)": 2.0,
               "share Lifecycle": 1.5})

PARAM_NAMES = ["mu", "kbar_max", "km_max", "tau_w", "lam_f0", "lam_u0", "lam_u1", "ybar_h", "sd_kT"]
BOUNDS = {"mu": (0.2, 5.0), "kbar_max": (0.005, 1.0), "km_max": (1.0, 15.0), "tau_w": (0.4, 1.2),
          "lam_f0": (0.05, 0.9), "lam_u0": (0.003, 0.05), "lam_u1": (0.003, 0.08), "ybar_h": (0.0, 0.6),
          "sd_kT": (0.001, 0.6), "delta_e": (0.001, 0.008), "lam_f_ratio": (0.5, 1.0)}
OPTIONAL_PARAMS = ["delta_e", "lam_f_ratio"]


def apply_params(base: FinalParams, x: dict) -> FinalParams:
    r_f = x.get("lam_f_ratio", 0.85)      # recession job-finding efficiency relative to expansion
    kw = dict(mu=x["mu"], kbar_max=x["kbar_max"], km_max=x["km_max"], tau_w=x["tau_w"],
              lam_f=(x["lam_f0"], r_f * x["lam_f0"]), lam_u=(x["lam_u0"], x["lam_u1"]),
              ybar_h=x["ybar_h"], sd_kT=x["sd_kT"])
    kw.update({n: x[n] for n in OPTIONAL_PARAMS if n in x and n in FinalParams.__dataclass_fields__})
    # carry through any other calibrated field of FinalParams (e.g. home_young_mult, nu_h, z_h, alpha_h)
    kw.update({k: v for k, v in x.items() if k in FinalParams.__dataclass_fields__ and k not in kw
               and k not in ("lam_f", "lam_u")})
    return base.replace(**kw)


def _to_unit(x, names=PARAM_NAMES):   # map params to R via logit of the bounded interval
    z = []
    for n in names:
        lo, hi = BOUNDS[n]; u = (x[n] - lo) / (hi - lo); u = min(max(u, 1e-6), 1 - 1e-6)
        z.append(np.log(u / (1 - u)))
    return np.array(z)


def _from_unit(z, names=PARAM_NAMES):
    x = {}
    for n, v in zip(names, z):
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


def evaluate(base: FinalParams, x: dict, cfg: SimConfigFinal, n_jobs=None):
    p = apply_params(base, x)
    sol = solve_all(p, n_jobs=n_jobs)
    sim = simulate_final(p, sol, cfg)
    m = moments_final(sim)
    obj, parts = objective_from_moments(m)
    return obj, m, parts, p


def run_smm(base: FinalParams, x0: dict, cfg: SimConfigFinal, log_path: str, maxfev=150, n_jobs=None,
            names=PARAM_NAMES):
    hist = []
    t0 = time.time()

    def fun(z):
        x = _from_unit(z, names)
        obj, m, parts, _ = evaluate(base, x, cfg, n_jobs)
        hist.append(dict(obj=obj, x=x, m={k: float(v) for k, v in m.items() if isinstance(v, (float, int))}))
        with open(log_path, "a") as fh:
            fh.write(json.dumps(dict(n=len(hist), t=round(time.time() - t0), obj=obj, x=x,
                                     dev={k: round(v, 3) for k, v in parts.items()})) + "\n")
        return obj

    z0 = _to_unit(x0, names)
    # explicit initial simplex: scipy's default perturbs each coordinate by 5% of its value, which is
    # ~0 for a parameter that starts at the midpoint of its bounds (logit = 0) and freezes it.
    simplex = np.vstack([z0] + [z0 + 0.6 * np.eye(len(z0))[i] for i in range(len(z0))])
    res = minimize(fun, z0, method="Nelder-Mead",
                   options=dict(maxfev=maxfev, xatol=1e-3, fatol=1e-4, initial_simplex=simplex))
    best = min(hist, key=lambda h: h["obj"])
    return best, hist, res


def global_screen(base: FinalParams, cfg: SimConfigFinal, n_points: int, log_path: str, seed=0, n_jobs=None,
                  names=PARAM_NAMES):
    """Latin-hypercube screening of the bounded parameter box (scipy.stats.qmc). Returns the
    evaluated points sorted by objective; each point is appended to the log."""
    from scipy.stats import qmc
    sampler = qmc.LatinHypercube(d=len(names), seed=seed)
    U = sampler.random(n_points)
    lo = np.array([BOUNDS[n][0] for n in names]); hi = np.array([BOUNDS[n][1] for n in names])
    pts = lo + U * (hi - lo)
    out = []; t0 = time.time()
    for i, row in enumerate(pts):
        x = dict(zip(names, row))
        obj, m, parts, _ = evaluate(base, x, cfg, n_jobs)
        out.append(dict(obj=obj, x=x, m={k: float(v) for k, v in m.items() if isinstance(v, (float, int))}))
        with open(log_path, "a") as fh:
            fh.write(json.dumps(dict(stage="global", n=i + 1, t=round(time.time() - t0), obj=obj, x=x,
                                     dev={k: round(v, 3) for k, v in parts.items()})) + "\n")
    out.sort(key=lambda r: r["obj"])
    return out
