"""Data moments of married women's monthly quit (E -> non-E, job leavers) and layoff (E -> non-E, job losers)
rates for the calibration: levels by NBER regime, trends, and cyclicality, on the full sample and on an
early window (the 1940s cohort's prime working years).

usage: python scripts/quit_layoff_moments.py --file <csv|xlsx|dta> [--quit eqmw_ma] [--layoff elmw_ma]
           [--early-end 1985] [--start 1976] [--end 2019] [--tag quit_layoff]
The file needs a monthly date (a `date`/`month`/`mdate` column, `year` + `month` columns, or a Stata monthly
date) and the two series (shares, or percent: detected from the level). "_ma" series are moving averages,
so the standard deviations below understate month-to-month noise but not the recession-expansion gaps.
Writes output/<tag>.json and output/<tag>.md.

Moments per window (full sample, early window, decades):
  mean quit and layoff rates; means by regime (NBER recession months, peak+1 to trough); ratio and gap;
  linear trend (per decade) and the trend-adjusted recession effect (rate on trend + recession dummy);
  sd of the log detrended series; the two-state cyclicality |log(rec/exp)| sqrt(pi_exp pi_rec) with the
  window's recession share and with the model's stationary share (0.143), comparable to the model's
  "sd log UE (women)" convention; correlation of quits and layoffs.
Suggested targets: quit/m exp and rec from the early window (trend-adjusted), lam_u0/lam_u1 (exogenous
job loss, expansion/recession) from the early window's layoff means.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from keam.final.params import FinalParams
from keam.final.simulate import stationary

NBER = [("1948-11", "1949-10"), ("1953-07", "1954-05"), ("1957-08", "1958-04"), ("1960-04", "1961-02"),
        ("1969-12", "1970-11"), ("1973-11", "1975-03"), ("1980-01", "1980-07"), ("1981-07", "1982-11"),
        ("1990-07", "1991-03"), ("2001-03", "2001-11"), ("2007-12", "2009-06"), ("2020-02", "2020-04")]

ap = argparse.ArgumentParser()
ap.add_argument("--file", required=True)
ap.add_argument("--quit", default="eqmw_ma")
ap.add_argument("--layoff", default="elmw_ma")
ap.add_argument("--early-end", type=int, default=1985, help="last year of the early window (inclusive)")
ap.add_argument("--start", type=int, default=0, help="first year used (0 = all)")
ap.add_argument("--end", type=int, default=2019, help="last year used (exclude the pandemic by default)")
ap.add_argument("--tag", default="quit_layoff")
ap.add_argument("--sheet", default=None)
ap.add_argument("--ma", type=int, default=12, help="length of the moving average behind the _ma series; the recession regressor is the same centred average of the NBER dummy (rows marked MA-consistent)")
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def read_any(f):
    ext = os.path.splitext(f)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(f, sheet_name=a.sheet or 0)
    if ext == ".dta":
        return pd.read_stata(f)
    return pd.read_csv(f)


df = read_any(a.file)
df.columns = [str(c).strip() for c in df.columns]
cols = {c.lower(): c for c in df.columns}
# ---- monthly date
date = None
for c in ("date", "mdate", "month_date", "period", "ym", "time"):
    if c in cols:
        s = df[cols[c]]
        if pd.api.types.is_numeric_dtype(s):            # Stata monthly (months since 1960-01) or yyyymm
            v = s.astype(float)
            if v.max() > 190000:                            # yyyymm
                date = pd.to_datetime(v.astype(int).astype(str), format="%Y%m")
            else:
                date = pd.to_datetime("1960-01-01") + pd.to_timedelta(0, "D")
                date = pd.PeriodIndex(year=1960 + (v // 12).astype(int), month=(v % 12 + 1).astype(int), freq="M").to_timestamp()
        else:
            date = pd.to_datetime(s)
        break
if date is None and "year" in cols and "month" in cols:
    date = pd.to_datetime(dict(year=df[cols["year"]].astype(int), month=df[cols["month"]].astype(int), day=1))
if date is None:
    raise SystemExit(f"no monthly date column found in {list(df.columns)}")
df = df.assign(date=pd.DatetimeIndex(date).to_period("M").to_timestamp())
qn = cols.get(a.quit.lower()); ln = cols.get(a.layoff.lower())
if qn is None or ln is None:
    raise SystemExit(f"series {a.quit!r} / {a.layoff!r} not in {list(df.columns)}")
d = df[["date", qn, ln]].rename(columns={qn: "quit", ln: "layoff"}).dropna().sort_values("date").reset_index(drop=True)
units = "share"
if d[["quit", "layoff"]].mean().max() > 0.5:      # percent
    d[["quit", "layoff"]] = d[["quit", "layoff"]] / 100.0; units = "percent (divided by 100)"
d["year"] = d.date.dt.year
if a.start:
    d = d[d.year >= a.start]
d = d[d.year <= a.end].reset_index(drop=True)
rec = np.zeros(len(d), bool)
for p, t in NBER:
    lo = (pd.Period(p, "M") + 1).to_timestamp(); hi = pd.Period(t, "M").to_timestamp()
    rec |= (d.date >= lo) & (d.date <= hi)
d["rec"] = rec
piz_model = stationary(FinalParams().piz); w_model = float(np.sqrt(piz_model[0] * piz_model[1]))
# length of the moving average behind the series: an MA(K) of noise gives monthly changes whose autocorrelation is
# about -0.5 at lag K and near zero elsewhere
dq = np.diff(d["quit"].values); dq = dq - dq.mean()
ac = {k: float(np.corrcoef(dq[:-k], dq[k:])[0, 1]) for k in range(2, 25)}
ma_inferred = min(ac, key=ac.get)


def stats(x: pd.DataFrame, name: str) -> dict:
    n = len(x); nr = int(x.rec.sum()); out = dict(window=name, months=n, rec_months=nr,
                                                  first=str(x.date.iloc[0].date()), last=str(x.date.iloc[-1].date()))
    t = (x.date - x.date.iloc[0]).dt.days.values / 365.25 / 10.0          # decades since the window start
    for v in ("quit", "layoff"):
        y = x[v].values
        out[f"{v} mean"] = float(y.mean())
        out[f"{v} exp"] = float(y[~x.rec].mean()); out[f"{v} rec"] = float(y[x.rec].mean()) if nr else float("nan")
        # trend + recession dummy
        X = np.column_stack([np.ones(n), t, x.rec.values.astype(float)]) if nr else np.column_stack([np.ones(n), t])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        out[f"{v} trend/decade"] = float(b[1])
        out[f"{v} rec effect (trend-adj.)"] = float(b[2]) if nr else float("nan")
        fit_exp_mid = float(b[0] + b[1] * t.mean())                          # expansion level at the window midpoint
        out[f"{v} exp (trend-adj., midpoint)"] = fit_exp_mid
        out[f"{v} rec (trend-adj., midpoint)"] = fit_exp_mid + (float(b[2]) if nr else float("nan"))
        resid = y - X @ b
        if nr and a.ma > 1:                      # MA-consistent: recession regressor averaged like the series
            rm = pd.Series(x.rec.values.astype(float)).rolling(a.ma, center=True, min_periods=1).mean().values
            Xm = np.column_stack([np.ones(n), t, rm]); bm, *_ = np.linalg.lstsq(Xm, y, rcond=None)
            out[f"{v} rec effect (trend-adj., MA-consistent)"] = float(bm[2])
            em = float(bm[0] + bm[1] * t.mean())
            out[f"{v} exp (MA-consistent, midpoint)"] = em
            out[f"{v} rec (MA-consistent, midpoint)"] = em + float(bm[2])
            out[f"{v} rec/exp (MA-consistent)"] = float((em + bm[2]) / em)
            out[f"{v} sd log two-state (MA-consistent, model share)"] = float(abs(np.log((em + bm[2]) / em)) * w_model)
        out[f"{v} sd detrended"] = float(resid.std())
        out[f"{v} sd log detrended"] = float((np.log(np.maximum(y, 1e-9)) - np.log(np.maximum(X @ b, 1e-9))).std())
        if nr:
            r = out[f"{v} rec"] / out[f"{v} exp"]; ra = out[f"{v} rec (trend-adj., midpoint)"] / fit_exp_mid
            out[f"{v} rec/exp"] = float(r); out[f"{v} rec/exp (trend-adj.)"] = float(ra)
            pr = nr / n; wd = float(np.sqrt(pr * (1 - pr)))
            out[f"{v} sd log two-state (data share)"] = float(abs(np.log(ra)) * wd)
            out[f"{v} sd log two-state (model share)"] = float(abs(np.log(ra)) * w_model)
    out["rec share"] = nr / n
    out["corr(quit, layoff)"] = float(np.corrcoef(x["quit"], x["layoff"])[0, 1]) if n > 2 else float("nan")
    return out


windows = [("full", d), (f"early (<= {a.early_end})", d[d.year <= a.early_end])]
for y0 in range(int(d.year.min()) // 10 * 10, int(d.year.max()) + 1, 10):
    sub = d[(d.year >= y0) & (d.year < y0 + 10)]
    if len(sub) >= 24:
        windows.append((f"{y0}s", sub))
res = [stats(x, nm) for nm, x in windows if len(x) >= 24]
early = res[1] if len(res) > 1 else res[0]
full = res[0]
sugg = {"quit/m exp": early["quit exp (trend-adj., midpoint)"], "quit/m rec": early["quit rec (trend-adj., midpoint)"],
        "quit/m rec (MA-consistent)": early.get("quit rec (MA-consistent, midpoint)"),
        "lam_u0": early["layoff exp (trend-adj., midpoint)"], "lam_u1": early["layoff rec (trend-adj., midpoint)"],
        "lam_u1 (MA-consistent)": early.get("layoff rec (MA-consistent, midpoint)"),
        "quit rec/exp (early)": early.get("quit rec/exp (trend-adj.)"), "quit rec/exp (full)": full.get("quit rec/exp (trend-adj.)"),
        "sd log quit two-state (early, model share)": early.get("quit sd log two-state (model share)"),
        "sd log quit two-state (full, model share)": full.get("quit sd log two-state (model share)")}
out = dict(file=a.file, series=dict(quit=qn, layoff=ln), units=units, ma_used=a.ma,
           ma_inferred=dict(lag=ma_inferred, autocorr_of_changes=ac[ma_inferred]), start=str(d.date.iloc[0].date()), end=str(d.date.iloc[-1].date()),
           early_end=a.early_end, model_rec_share=float(piz_model[1]), windows=res, suggested_targets=sugg)
json.dump(out, open(os.path.join(ROOT, "output", f"{a.tag}.json"), "w"), indent=1)

keys = ["months", "rec months", "quit mean", "quit exp", "quit rec", "quit rec/exp", "quit trend/decade", "quit rec effect (trend-adj.)",
        "quit exp (trend-adj., midpoint)", "quit rec (trend-adj., midpoint)", "quit rec/exp (trend-adj.)",
        "quit rec effect (trend-adj., MA-consistent)", "quit rec (MA-consistent, midpoint)", "quit rec/exp (MA-consistent)",
        "quit sd log detrended", "quit sd log two-state (data share)", "quit sd log two-state (model share)",
        "quit sd log two-state (MA-consistent, model share)",
        "layoff mean", "layoff exp", "layoff rec", "layoff rec/exp", "layoff trend/decade", "layoff rec effect (trend-adj.)",
        "layoff exp (trend-adj., midpoint)", "layoff rec (trend-adj., midpoint)",
        "layoff rec effect (trend-adj., MA-consistent)", "layoff rec (MA-consistent, midpoint)", "layoff rec/exp (MA-consistent)",
        "layoff sd log detrended", "corr(quit, layoff)"]
fmt = lambda v: ("-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else (f"{v:d}" if isinstance(v, int) else f"{v:.4f}"))
L = [f"# Married women's quit and layoff rates (`{os.path.basename(a.file)}`: `{qn}`, `{ln}`; {units}; "
     f"{out['start']} to {out['end']}; NBER recession months)", "",
     "Trend-adjusted values: rate regressed on a linear trend and a recession dummy within the window; the expansion level "
     "is the fitted value at the window midpoint and the recession level adds the dummy. Two-state cyclicality: "
     "|log(rec/exp)| sqrt(pi_exp pi_rec) with the window's recession share and with the model's stationary share "
     f"({piz_model[1]:.3f}), the convention of the model's `sd log UE (women)` moment. MA-consistent rows: the recession "
     f"regressor is the centred {a.ma}-month average of the NBER dummy, matching the smoothing of the _ma series "
     f"(inferred length: the autocorrelation of the quit series' monthly changes is {ac[ma_inferred]:+.2f} at lag {ma_inferred}, "
     f"near zero at other lags).", "",
     "| moment | " + " | ".join(r["window"] for r in res) + " |", "|---|" + "---|" * len(res)]
for k in keys:
    kk = "rec_months" if k == "rec months" else k
    L.append(f"| {k} | " + " | ".join(fmt(r.get(kk)) for r in res) + " |")
L += ["", f"## Suggested targets (early window, <= {a.early_end})", "", "| target | value |", "|---|---|"]
L += [f"| {k} | {fmt(v)} |" for k, v in sugg.items()]
L += ["", "Exogenous separation `lam_u` = layoff rate by regime (fixed, not calibrated); quit targets from the quit series; "
      "the E->nonE targets become redundant (quits + layoffs) and are dropped."]
open(os.path.join(ROOT, "output", f"{a.tag}.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
