# Married women's quit and layoff rates (`QuitLayoff2024.csv`: `nemw_seats`, `numw_seats`; percent (divided by 100); 1978-02-01 to 2016-12-01; NBER recession months)

Trend-adjusted values: rate regressed on a linear trend and a recession dummy within the window; the expansion level is the fitted value at the window midpoint and the recession level adds the dummy. Two-state cyclicality: |log(rec/exp)| sqrt(pi_exp pi_rec) with the window's recession share and with the model's stationary share (0.143), the convention of the model's `sd log UE (women)` moment. MA-consistent rows: the recession regressor is the centred 1-month average of the NBER dummy, matching the smoothing of the _ma series (inferred length: the autocorrelation of the quit series' monthly changes is -0.06 at lag 14, near zero at other lags).

| moment | full | early (<= 1985) | 1980s | 1990s | 2000s | 2010s |
|---|---|---|---|---|---|---|
| months | 467 | 95 | 120 | 120 | 120 | 84 |
| rec months | 56 | 22 | 22 | 8 | 26 | 0 |
| quit mean | 0.0539 | 0.0529 | 0.0563 | 0.0566 | 0.0544 | 0.0467 |
| quit exp | 0.0540 | 0.0531 | 0.0572 | 0.0566 | 0.0548 | 0.0467 |
| quit rec | 0.0529 | 0.0521 | 0.0521 | 0.0558 | 0.0527 | - |
| quit rec/exp | 0.9798 | 0.9815 | 0.9111 | 0.9848 | 0.9621 | - |
| quit trend/decade | -0.0021 | 0.0074 | 0.0134 | -0.0003 | -0.0076 | 0.0018 |
| quit rec effect (trend-adj.) | -0.0018 | -0.0007 | 0.0003 | -0.0010 | -0.0006 | - |
| quit exp (trend-adj., midpoint) | 0.0541 | 0.0530 | 0.0562 | 0.0566 | 0.0545 | 0.0467 |
| quit rec (trend-adj., midpoint) | 0.0523 | 0.0524 | 0.0565 | 0.0557 | 0.0539 | - |
| quit rec/exp (trend-adj.) | 0.9670 | 0.9872 | 1.0049 | 0.9828 | 0.9894 | - |
| quit rec effect (trend-adj., MA-consistent) | - | - | - | - | - | - |
| quit rec (MA-consistent, midpoint) | - | - | - | - | - | - |
| quit rec/exp (MA-consistent) | - | - | - | - | - | - |
| quit sd log detrended | 0.1503 | 0.1069 | 0.1150 | 0.1401 | 0.1493 | 0.1436 |
| quit sd log two-state (data share) | 0.0109 | 0.0054 | 0.0019 | 0.0043 | 0.0044 | - |
| quit sd log two-state (model share) | 0.0117 | 0.0045 | 0.0017 | 0.0061 | 0.0037 | - |
| quit sd log two-state (MA-consistent, model share) | - | - | - | - | - | - |
| layoff mean | 0.0121 | 0.0108 | 0.0113 | 0.0129 | 0.0120 | 0.0129 |
| layoff exp | 0.0121 | 0.0106 | 0.0113 | 0.0130 | 0.0117 | 0.0129 |
| layoff rec | 0.0124 | 0.0114 | 0.0114 | 0.0128 | 0.0131 | - |
| layoff rec/exp | 1.0266 | 1.0684 | 1.0031 | 0.9861 | 1.1253 | - |
| layoff trend/decade | 0.0005 | 0.0028 | 0.0009 | -0.0023 | 0.0015 | -0.0103 |
| layoff rec effect (trend-adj.) | 0.0005 | 0.0008 | 0.0004 | -0.0012 | 0.0012 | - |
| layoff exp (trend-adj., midpoint) | 0.0120 | 0.0106 | 0.0113 | 0.0130 | 0.0117 | 0.0129 |
| layoff rec (trend-adj., midpoint) | 0.0125 | 0.0114 | 0.0116 | 0.0118 | 0.0129 | - |
| layoff rec effect (trend-adj., MA-consistent) | - | - | - | - | - | - |
| layoff rec (MA-consistent, midpoint) | - | - | - | - | - | - |
| layoff rec/exp (MA-consistent) | - | - | - | - | - | - |
| layoff sd log detrended | 0.3004 | 0.2053 | 0.2295 | 0.2642 | 0.2859 | 0.3608 |
| corr(quit, layoff) | -0.0733 | 0.1485 | 0.0902 | -0.0435 | -0.1596 | -0.1263 |

## Suggested targets (early window, <= 1985)

| target | value |
|---|---|
| quit/m exp | 0.0530 |
| quit/m rec | 0.0524 |
| quit/m rec (MA-consistent) | - |
| lam_u0 | 0.0106 |
| lam_u1 | 0.0114 |
| lam_u1 (MA-consistent) | - |
| quit rec/exp (early) | 0.9872 |
| quit rec/exp (full) | 0.9670 |
| sd log quit two-state (early, model share) | 0.0045 |
| sd log quit two-state (full, model share) | 0.0117 |

Exogenous separation `lam_u` = layoff rate by regime (fixed, not calibrated); quit targets from the quit series; the E->nonE targets become redundant (quits + layoffs) and are dropped.
