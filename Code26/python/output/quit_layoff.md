# Married women's quit and layoff rates (`QuitLayoff2024_sa.csv`: `eqmw_ma`, `elmw_ma`; percent (divided by 100); 1978-01-01 to 2019-12-01; NBER recession months)

Trend-adjusted values: rate regressed on a linear trend and a recession dummy within the window; the expansion level is the fitted value at the window midpoint and the recession level adds the dummy. Two-state cyclicality: |log(rec/exp)| sqrt(pi_exp pi_rec) with the window's recession share and with the model's stationary share (0.143), the convention of the model's `sd log UE (women)` moment. MA-consistent rows: the recession regressor is the centred 13-month average of the NBER dummy, matching the smoothing of the _ma series (inferred length: the autocorrelation of the quit series' monthly changes is -0.41 at lag 13, near zero at other lags).

| moment | full | early (<= 1985) | 1970s | 1980s | 1990s | 2000s | 2010s |
|---|---|---|---|---|---|---|---|
| months | 504 | 96 | 24 | 120 | 120 | 120 | 120 |
| rec months | 56 | 22 | 0 | 22 | 8 | 26 | 0 |
| quit mean | 0.0165 | 0.0226 | 0.0262 | 0.0210 | 0.0157 | 0.0142 | 0.0133 |
| quit exp | 0.0164 | 0.0229 | 0.0262 | 0.0208 | 0.0154 | 0.0143 | 0.0133 |
| quit rec | 0.0177 | 0.0216 | - | 0.0216 | 0.0194 | 0.0139 | - |
| quit rec/exp | 1.0789 | 0.9399 | - | 1.0358 | 1.2580 | 0.9662 | - |
| quit trend/decade | -0.0027 | -0.0094 | -0.0051 | -0.0026 | -0.0055 | -0.0018 | 0.0042 |
| quit rec effect (trend-adj.) | -0.0001 | -0.0017 | - | -0.0003 | 0.0016 | -0.0001 | - |
| quit exp (trend-adj., midpoint) | 0.0165 | 0.0230 | 0.0262 | 0.0210 | 0.0156 | 0.0143 | 0.0133 |
| quit rec (trend-adj., midpoint) | 0.0165 | 0.0213 | - | 0.0207 | 0.0171 | 0.0141 | - |
| quit rec/exp (trend-adj.) | 0.9959 | 0.9256 | - | 0.9850 | 1.0997 | 0.9906 | - |
| quit rec effect (trend-adj., MA-consistent) | -0.0001 | -0.0031 | - | -0.0005 | 0.0050 | -0.0005 | - |
| quit rec (MA-consistent, midpoint) | 0.0164 | 0.0202 | - | 0.0206 | 0.0203 | 0.0138 | - |
| quit rec/exp (MA-consistent) | 0.9942 | 0.8676 | - | 0.9769 | 1.3230 | 0.9620 | - |
| quit sd log detrended | 0.1466 | 0.0464 | 0.0215 | 0.0616 | 0.0952 | 0.0980 | 0.0851 |
| quit sd log two-state (data share) | 0.0013 | 0.0325 | - | 0.0058 | 0.0237 | 0.0039 | - |
| quit sd log two-state (model share) | 0.0014 | 0.0271 | - | 0.0053 | 0.0333 | 0.0033 | - |
| quit sd log two-state (MA-consistent, model share) | 0.0020 | 0.0497 | - | 0.0082 | 0.0979 | 0.0136 | - |
| layoff mean | 0.0101 | 0.0137 | 0.0124 | 0.0128 | 0.0103 | 0.0086 | 0.0085 |
| layoff exp | 0.0099 | 0.0133 | 0.0124 | 0.0123 | 0.0102 | 0.0083 | 0.0085 |
| layoff rec | 0.0121 | 0.0151 | - | 0.0151 | 0.0129 | 0.0094 | - |
| layoff rec/exp | 1.2264 | 1.1317 | - | 1.2269 | 1.2674 | 1.1355 | - |
| layoff trend/decade | -0.0015 | 0.0019 | 0.0038 | -0.0053 | -0.0049 | -0.0005 | -0.0020 |
| layoff rec effect (trend-adj.) | 0.0015 | 0.0018 | - | 0.0006 | 0.0006 | 0.0012 | - |
| layoff exp (trend-adj., midpoint) | 0.0100 | 0.0133 | 0.0124 | 0.0127 | 0.0103 | 0.0083 | 0.0085 |
| layoff rec (trend-adj., midpoint) | 0.0115 | 0.0151 | - | 0.0133 | 0.0109 | 0.0095 | - |
| layoff rec effect (trend-adj., MA-consistent) | 0.0023 | 0.0030 | - | 0.0009 | 0.0004 | 0.0022 | - |
| layoff rec (MA-consistent, midpoint) | 0.0122 | 0.0160 | - | 0.0135 | 0.0107 | 0.0102 | - |
| layoff rec/exp (MA-consistent) | 1.2328 | 1.2322 | - | 1.0692 | 1.0397 | 1.2662 | - |
| layoff sd log detrended | 0.1437 | 0.0496 | 0.0287 | 0.0915 | 0.0573 | 0.1237 | 0.1482 |
| corr(quit, layoff) | 0.6284 | -0.6758 | -0.7025 | 0.1951 | 0.7070 | -0.7264 | -0.5599 |

## Suggested targets (early window, <= 1985)

| target | value |
|---|---|
| quit/m exp | 0.0230 |
| quit/m rec | 0.0213 |
| quit/m rec (MA-consistent) | 0.0202 |
| lam_u0 | 0.0133 |
| lam_u1 | 0.0151 |
| lam_u1 (MA-consistent) | 0.0160 |
| quit rec/exp (early) | 0.9256 |
| quit rec/exp (full) | 0.9959 |
| sd log quit two-state (early, model share) | 0.0271 |
| sd log quit two-state (full, model share) | 0.0014 |

Exogenous separation `lam_u` = layoff rate by regime (fixed, not calibrated); quit targets from the quit series; the E->nonE targets become redundant (quits + layoffs) and are dropped.
