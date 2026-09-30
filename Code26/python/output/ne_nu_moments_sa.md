# Married women's quit and layoff rates (`QuitLayoff2024_sa.csv`: `nemw_ma`, `numw_ma`; percent (divided by 100); 1978-01-01 to 2019-12-01; NBER recession months)

Trend-adjusted values: rate regressed on a linear trend and a recession dummy within the window; the expansion level is the fitted value at the window midpoint and the recession level adds the dummy. Two-state cyclicality: |log(rec/exp)| sqrt(pi_exp pi_rec) with the window's recession share and with the model's stationary share (0.143), the convention of the model's `sd log UE (women)` moment. MA-consistent rows: the recession regressor is the centred 12-month average of the NBER dummy, matching the smoothing of the _ma series.

| moment | full | early (<= 1985) | 1970s | 1980s | 1990s | 2000s | 2010s |
|---|---|---|---|---|---|---|---|
| months | 504 | 96 | 24 | 120 | 120 | 120 | 120 |
| rec months | 56 | 22 | 0 | 22 | 8 | 26 | 0 |
| quit mean | 0.0468 | 0.0449 | 0.0448 | 0.0472 | 0.0485 | 0.0485 | 0.0436 |
| quit exp | 0.0469 | 0.0452 | 0.0448 | 0.0479 | 0.0483 | 0.0487 | 0.0436 |
| quit rec | 0.0466 | 0.0439 | - | 0.0439 | 0.0509 | 0.0475 | - |
| quit rec/exp | 0.9936 | 0.9699 | - | 0.9157 | 1.0533 | 0.9751 | - |
| quit trend/decade | -0.0006 | 0.0021 | 0.0150 | 0.0099 | 0.0006 | -0.0044 | 0.0057 |
| quit rec effect (trend-adj.) | -0.0006 | -0.0013 | - | -0.0001 | 0.0028 | -0.0003 | - |
| quit exp (trend-adj., midpoint) | 0.0469 | 0.0452 | 0.0448 | 0.0472 | 0.0483 | 0.0486 | 0.0436 |
| quit rec (trend-adj., midpoint) | 0.0463 | 0.0439 | - | 0.0471 | 0.0511 | 0.0482 | - |
| quit rec/exp (trend-adj.) | 0.9873 | 0.9715 | - | 0.9982 | 1.0586 | 0.9930 | - |
| quit rec effect (trend-adj., MA-consistent) | -0.0012 | -0.0024 | - | -0.0002 | 0.0064 | -0.0011 | - |
| quit rec (MA-consistent, midpoint) | 0.0458 | 0.0431 | - | 0.0470 | 0.0544 | 0.0477 | - |
| quit rec/exp (MA-consistent) | 0.9740 | 0.9474 | - | 0.9965 | 1.1331 | 0.9781 | - |
| quit sd log detrended | 0.0645 | 0.0292 | 0.0035 | 0.0334 | 0.0402 | 0.0366 | 0.0270 |
| quit sd log two-state (data share) | 0.0040 | 0.0122 | - | 0.0007 | 0.0142 | 0.0029 | - |
| quit sd log two-state (model share) | 0.0045 | 0.0101 | - | 0.0006 | 0.0199 | 0.0025 | - |
| quit sd log two-state (MA-consistent, model share) | 0.0092 | 0.0189 | - | 0.0012 | 0.0437 | 0.0077 | - |
| layoff mean | 0.0198 | 0.0187 | 0.0155 | 0.0194 | 0.0207 | 0.0196 | 0.0204 |
| layoff exp | 0.0197 | 0.0185 | 0.0155 | 0.0194 | 0.0207 | 0.0192 | 0.0204 |
| layoff rec | 0.0203 | 0.0194 | - | 0.0194 | 0.0199 | 0.0211 | - |
| layoff rec/exp | 1.0263 | 1.0512 | - | 1.0023 | 0.9584 | 1.0985 | - |
| layoff trend/decade | 0.0002 | 0.0091 | 0.0090 | -0.0004 | -0.0042 | 0.0046 | -0.0193 |
| layoff rec effect (trend-adj.) | 0.0006 | 0.0013 | - | -0.0001 | -0.0027 | 0.0010 | - |
| layoff exp (trend-adj., midpoint) | 0.0197 | 0.0184 | 0.0155 | 0.0194 | 0.0209 | 0.0194 | 0.0204 |
| layoff rec (trend-adj., midpoint) | 0.0204 | 0.0197 | - | 0.0193 | 0.0182 | 0.0204 | - |
| layoff rec effect (trend-adj., MA-consistent) | 0.0012 | 0.0022 | - | -0.0003 | -0.0067 | 0.0022 | - |
| layoff rec (MA-consistent, midpoint) | 0.0209 | 0.0204 | - | 0.0191 | 0.0145 | 0.0214 | - |
| layoff rec/exp (MA-consistent) | 1.0616 | 1.1199 | - | 0.9842 | 0.6850 | 1.1170 | - |
| layoff sd log detrended | 0.1730 | 0.0361 | 0.0151 | 0.0704 | 0.0812 | 0.1094 | 0.0842 |
| corr(quit, layoff) | -0.4075 | 0.0214 | 0.8904 | -0.4192 | -0.7231 | -0.8408 | -0.8335 |

## Suggested targets (early window, <= 1985)

| target | value |
|---|---|
| quit/m exp | 0.0452 |
| quit/m rec | 0.0439 |
| quit/m rec (MA-consistent) | 0.0431 |
| lam_u0 | 0.0184 |
| lam_u1 | 0.0197 |
| lam_u1 (MA-consistent) | 0.0204 |
| quit rec/exp (early) | 0.9715 |
| quit rec/exp (full) | 0.9873 |
| sd log quit two-state (early, model share) | 0.0101 |
| sd log quit two-state (full, model share) | 0.0045 |

Exogenous separation `lam_u` = layoff rate by regime (fixed, not calibrated); quit targets from the quit series; the E->nonE targets become redundant (quits + layoffs) and are dropped.
