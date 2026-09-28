# Married women's quit and layoff rates (`QuitLayoff2024.csv`: `eqmw_seats`, `elmw_seats`; percent (divided by 100); 1978-02-01 to 2016-12-01; NBER recession months)

Trend-adjusted values: rate regressed on a linear trend and a recession dummy within the window; the expansion level is the fitted value at the window midpoint and the recession level adds the dummy. Two-state cyclicality: |log(rec/exp)| sqrt(pi_exp pi_rec) with the window's recession share and with the model's stationary share (0.143), the convention of the model's `sd log UE (women)` moment. MA-consistent rows: the recession regressor is the centred 1-month average of the NBER dummy, matching the smoothing of the _ma series (inferred length: the autocorrelation of the quit series' monthly changes is -0.21 at lag 24, near zero at other lags).

| moment | full | early (<= 1985) | 1980s | 1990s | 2000s | 2010s |
|---|---|---|---|---|---|---|
| months | 467 | 95 | 120 | 120 | 120 | 84 |
| rec months | 56 | 22 | 22 | 8 | 26 | 0 |
| quit mean | 0.0162 | 0.0223 | 0.0203 | 0.0148 | 0.0139 | 0.0126 |
| quit exp | 0.0160 | 0.0225 | 0.0201 | 0.0145 | 0.0139 | 0.0126 |
| quit rec | 0.0176 | 0.0213 | 0.0213 | 0.0190 | 0.0139 | - |
| quit rec/exp | 1.1003 | 0.9458 | 1.0596 | 1.3126 | 1.0002 | - |
| quit trend/decade | -0.0030 | -0.0100 | -0.0036 | -0.0046 | -0.0016 | 0.0074 |
| quit rec effect (trend-adj.) | 0.0006 | -0.0016 | -0.0002 | 0.0025 | 0.0003 | - |
| quit exp (trend-adj., midpoint) | 0.0161 | 0.0226 | 0.0204 | 0.0146 | 0.0139 | 0.0126 |
| quit rec (trend-adj., midpoint) | 0.0167 | 0.0210 | 0.0202 | 0.0171 | 0.0142 | - |
| quit rec/exp (trend-adj.) | 1.0369 | 0.9280 | 0.9889 | 1.1730 | 1.0223 | - |
| quit rec effect (trend-adj., MA-consistent) | - | - | - | - | - | - |
| quit rec (MA-consistent, midpoint) | - | - | - | - | - | - |
| quit rec/exp (MA-consistent) | - | - | - | - | - | - |
| quit sd log detrended | 0.2203 | 0.1499 | 0.1447 | 0.1805 | 0.2021 | 0.2150 |
| quit sd log two-state (data share) | 0.0118 | 0.0315 | 0.0043 | 0.0398 | 0.0091 | - |
| quit sd log two-state (model share) | 0.0127 | 0.0262 | 0.0039 | 0.0558 | 0.0077 | - |
| quit sd log two-state (MA-consistent, model share) | - | - | - | - | - | - |
| layoff mean | 0.0098 | 0.0136 | 0.0127 | 0.0099 | 0.0079 | 0.0079 |
| layoff exp | 0.0095 | 0.0130 | 0.0121 | 0.0096 | 0.0076 | 0.0079 |
| layoff rec | 0.0120 | 0.0153 | 0.0153 | 0.0133 | 0.0087 | - |
| layoff rec/exp | 1.2570 | 1.1793 | 1.2706 | 1.3820 | 1.1489 | - |
| layoff trend/decade | -0.0018 | 0.0016 | -0.0054 | -0.0049 | -0.0007 | -0.0048 |
| layoff rec effect (trend-adj.) | 0.0018 | 0.0024 | 0.0011 | 0.0015 | 0.0013 | - |
| layoff exp (trend-adj., midpoint) | 0.0096 | 0.0130 | 0.0125 | 0.0098 | 0.0076 | 0.0079 |
| layoff rec (trend-adj., midpoint) | 0.0115 | 0.0154 | 0.0136 | 0.0113 | 0.0089 | - |
| layoff rec effect (trend-adj., MA-consistent) | - | - | - | - | - | - |
| layoff rec (MA-consistent, midpoint) | - | - | - | - | - | - |
| layoff rec/exp (MA-consistent) | - | - | - | - | - | - |
| layoff sd log detrended | 0.2613 | 0.1700 | 0.1972 | 0.2158 | 0.2854 | 0.2444 |
| corr(quit, layoff) | 0.3489 | -0.3169 | 0.0053 | 0.0898 | -0.4531 | -0.5036 |

## Suggested targets (early window, <= 1985)

| target | value |
|---|---|
| quit/m exp | 0.0226 |
| quit/m rec | 0.0210 |
| quit/m rec (MA-consistent) | - |
| lam_u0 | 0.0130 |
| lam_u1 | 0.0154 |
| lam_u1 (MA-consistent) | - |
| quit rec/exp (early) | 0.9280 |
| quit rec/exp (full) | 1.0369 |
| sd log quit two-state (early, model share) | 0.0262 |
| sd log quit two-state (full, model share) | 0.0127 |

Exogenous separation `lam_u` = layoff rate by regime (fixed, not calibrated); quit targets from the quit series; the E->nonE targets become redundant (quits + layoffs) and are dropped.
