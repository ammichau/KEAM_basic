# Final model results: 1940s cohort calibration, trend experiments, mechanism

All numbers are produced by scripts in `Code26/python/scripts`; the files cited are in `Code26/python/output`. Model specification: `FINAL_MODEL.md`.

## 1. Calibration of the 1940s cohort

Source: `output/final_calib_full.json` (objective 0.149, 60 evaluations, 100 types).

| parameter | value |
|---|---|
| mu | 0.9086 |
| kbar_max | 0.0172 |
| km_max | 5.8765 |
| tau_w | 0.7571 |
| lam_f0 | 0.3910 |
| lam_u0 | 0.0169 |
| lam_u1 | 0.0196 |
| ybar_h | 0.0660 |
| sd_kT | 0.2897 |
| home_young_mult | 1.8119 |
| nu_h | 0.6876 |
| z_h | 0.4517 |
| alpha_h | 0.2505 |
| e_max | 1.9685 |
| kappa_h_power | 0.2866 |

| target | data | model | deviation |
|---|---|---|---|
| E/pop | 0.6200 | 0.6644 | +7.2% |
| hours|E | 0.4000 | 0.4141 | +3.5% |
| share Lifecycle | 0.3100 | 0.3038 | -2.0% |
| share PT | 0.2800 | 0.2531 | -9.6% |
| share Career | 0.1900 | 0.1833 | -3.5% |
| share NiLF | 0.2200 | 0.2598 | +18.1% |
| quit/m exp | 0.0340 | 0.0311 | -8.4% |
| quit/m rec | 0.0280 | 0.0237 | -15.3% |
| E->nonE/m exp | 0.0500 | 0.0485 | -3.1% |
| E->nonE/m rec | 0.0480 | 0.0431 | -10.1% |
| dE/pop rec-exp (pts) | -1.7000 | -1.5984 | +10.2% |
| wage gap (hourly ratio) | 0.7100 | 0.7477 | +5.3% |

| untargeted moment | model |
|---|---|
| U rate | 0.0459 |
| wife share exp | 0.3370 |
| wife share rec | 0.3441 |
| HH income rec/exp - 1 (%) | -14.9712 |
| cons drop at H job loss exp (%) | -5.0354 |
| cons drop at H job loss rec (%) | -7.1846 |
| mean assets/monthly HH inc | 1.4722 |
| share e at cap | 0.2488 |

## 2. Single-factor experiments sized to the 1970s employment rate

Source: `output/final_results_full.json`. Scales: returns to experience x1.375, compensated wage gap x1.084 (husband income scaled to keep household income constant at baseline behaviour), cost of work x0.100.

| moment | baseline | RoE x1.38 | comp. wage gap x1.08 | cost x0.10 |
|---|---|---|---|---|
| E/pop | 0.6644 | 0.7281 | 0.7336 | 0.7263 |
| hours|E | 0.4141 | 0.4437 | 0.4413 | 0.4194 |
| U rate | 0.0459 | 0.0473 | 0.0461 | 0.0460 |
| quit/m exp | 0.0311 | 0.0215 | 0.0204 | 0.0227 |
| quit/m rec | 0.0237 | 0.0163 | 0.0151 | 0.0165 |
| E->nonE/m exp | 0.0485 | 0.0388 | 0.0376 | 0.0399 |
| E->nonE/m rec | 0.0431 | 0.0357 | 0.0345 | 0.0359 |
| dE/pop rec-exp (pts) | -1.5984 | -1.5558 | -1.3690 | -1.9532 |
| wife share exp | 0.3370 | 0.4092 | 0.4021 | 0.3674 |
| wife share rec | 0.3441 | 0.4167 | 0.4091 | 0.3737 |
| wage gap (hourly ratio) | 0.7477 | 0.8730 | 0.8528 | 0.7678 |
| share Lifecycle | 0.3038 | 0.2935 | 0.3244 | 0.1544 |
| share PT | 0.2531 | 0.1812 | 0.2123 | 0.3415 |
| share Career | 0.1833 | 0.3375 | 0.3025 | 0.2927 |
| share NiLF | 0.2598 | 0.1877 | 0.1608 | 0.2115 |
| HH income rec/exp - 1 (%) | -14.9712 | -14.7435 | -14.8837 | -15.2162 |
| cons drop at H job loss exp (%) | -5.0354 | -4.4100 | -4.4597 | -4.8908 |
| cons drop at H job loss rec (%) | -7.1846 | -6.6915 | -6.6684 | -7.1765 |
| mean assets/monthly HH inc | 1.4722 | 1.7243 | 1.6511 | 1.4745 |

Change in the cyclical quit gap (recession minus expansion monthly quit rate, percentage points) and in the recession employment drop relative to the baseline:

| experiment | quit gap | ΔE/pop rec-exp (pts) | E/pop |
|---|---|---|---|
| baseline | -0.742 | -1.598 | 0.664 |
| RoE x1.38 | -0.524 | -1.556 | 0.728 |
| comp. wage gap x1.08 | -0.523 | -1.369 | 0.734 |
| cost x0.10 | -0.619 | -1.953 | 0.726 |

Supplementary experiments (`output/extra_experiments_full.json`): child-care cost scaled toward zero (x0.44 of the excess home productivity at 25-39) and all cost components scaled jointly (x0.14).

| moment | baseline | child-care cost x0.44 | all costs x0.14 |
|---|---|---|---|
| E/pop | 0.6644 | 0.7340 | 0.7283 |
| hours|E | 0.4141 | 0.4319 | 0.4195 |
| U rate | 0.0459 | 0.0494 | 0.0459 |
| quit/m exp | 0.0311 | 0.0206 | 0.0224 |
| quit/m rec | 0.0237 | 0.0155 | 0.0163 |
| E->nonE/m exp | 0.0485 | 0.0379 | 0.0396 |
| E->nonE/m rec | 0.0431 | 0.0348 | 0.0357 |
| dE/pop rec-exp (pts) | -1.5984 | -2.4170 | -2.0426 |
| wife share exp | 0.3370 | 0.3844 | 0.3686 |
| wife share rec | 0.3441 | 0.3904 | 0.3749 |
| wage gap (hourly ratio) | 0.7477 | 0.7866 | 0.7693 |
| share Lifecycle | 0.3038 | 0.0844 | 0.1412 |
| share PT | 0.2531 | 0.2156 | 0.3502 |
| share Career | 0.1833 | 0.4690 | 0.2963 |
| share NiLF | 0.2598 | 0.2310 | 0.2123 |
| HH income rec/exp - 1 (%) | -14.9712 | -15.4861 | -15.2400 |
| cons drop at H job loss exp (%) | -5.0354 | -4.5864 | -4.8865 |
| cons drop at H job loss rec (%) | -7.1846 | -6.9401 | -7.1631 |
| mean assets/monthly HH inc | 1.4722 | 1.4944 | 1.4759 |

| experiment | quit gap | ΔE/pop rec-exp (pts) | E/pop |
|---|---|---|---|
| child-care cost x0.44 | -0.511 | -2.417 | 0.734 |
| all costs x0.14 | -0.604 | -2.043 | 0.728 |

## 3. Cohort accounting

τ_w and γ_e follow the slides (p.35) relative to 1940 (wage gap 0.71, 0.74, 0.77, 0.76, 0.77; γ_e 0.50, 0.55, 0.58, 0.68, 0.69), with the husband's income compensated; the cost of work is scaled to reproduce each cohort's employment rate (0.62, 0.67, 0.71, 0.73, 0.72).

| moment | 1940 | 1950 (cost x1.77) | 1960 (cost x1.77) | 1970 (cost x1.77) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6644 | 0.6703 | 0.7121 | 0.7316 | 0.7324 |
| hours|E | 0.4141 | 0.4355 | 0.4515 | 0.4601 | 0.4645 |
| U rate | 0.0459 | 0.0470 | 0.0481 | 0.0485 | 0.0488 |
| quit/m exp | 0.0311 | 0.0278 | 0.0217 | 0.0191 | 0.0183 |
| quit/m rec | 0.0237 | 0.0206 | 0.0157 | 0.0138 | 0.0132 |
| E->nonE/m exp | 0.0485 | 0.0451 | 0.0390 | 0.0364 | 0.0356 |
| E->nonE/m rec | 0.0431 | 0.0400 | 0.0350 | 0.0332 | 0.0325 |
| dE/pop rec-exp (pts) | -1.5984 | -1.2115 | -1.0827 | -1.2304 | -1.0515 |
| wife share exp | 0.3370 | 0.3692 | 0.4117 | 0.4369 | 0.4447 |
| wife share rec | 0.3441 | 0.3763 | 0.4202 | 0.4456 | 0.4538 |
| wage gap (hourly ratio) | 0.7477 | 0.8231 | 0.8991 | 0.9523 | 0.9767 |
| share Lifecycle | 0.3038 | 0.3771 | 0.3881 | 0.3765 | 0.3908 |
| share PT | 0.2531 | 0.1973 | 0.1719 | 0.1454 | 0.1388 |
| share Career | 0.1833 | 0.1950 | 0.2675 | 0.3287 | 0.3283 |
| share NiLF | 0.2598 | 0.2306 | 0.1725 | 0.1494 | 0.1421 |
| HH income rec/exp - 1 (%) | -14.9712 | -14.8559 | -14.4955 | -14.4351 | -14.3072 |
| cons drop at H job loss exp (%) | -5.0354 | -4.6720 | -4.2855 | -4.0319 | -3.9346 |
| cons drop at H job loss rec (%) | -7.1846 | -6.8833 | -6.4779 | -6.2469 | -6.1531 |
| mean assets/monthly HH inc | 1.4722 | 1.6379 | 1.7272 | 1.8225 | 1.8305 |

### 3b. Cohort accounting, refined: cost scale and τ_w solved jointly

Source: `output/cohorts_refined_full.json`. For each cohort the cost scale and τ_w (husband's income compensated) are solved so that the cohort's employment rate and its measured within-couple wage gap (data ratio applied to the model's 1940 gap) both match, given the cohort's γ_e.

| | 1940 | 1950 | 1960 | 1970 | 1980 |
|---|---|---|---|---|---|
| cost scale | 1.000 | 1.192 | 0.754 | 0.089 | 0.334 |
| τ_w | 0.757 | 0.758 | 0.759 | 0.706 | 0.711 |
| E/pop | 0.6644 | 0.6700 | 0.7100 | 0.7299 | 0.7200 |
| wage gap (hourly ratio) | 0.7477 | 0.7794 | 0.8109 | 0.8004 | 0.8110 |
| quit/m exp | 0.0311 | 0.0299 | 0.0247 | 0.0227 | 0.0239 |
| quit/m rec | 0.0237 | 0.0225 | 0.0184 | 0.0167 | 0.0178 |
| dE/pop rec-exp (pts) | -1.5984 | -1.5061 | -1.6626 | -2.1952 | -2.2357 |
| wife share exp | 0.3370 | 0.3530 | 0.3793 | 0.3838 | 0.3842 |
| share Lifecycle | 0.3038 | 0.3185 | 0.2752 | 0.1506 | 0.1906 |
| share PT | 0.2531 | 0.2290 | 0.2290 | 0.2679 | 0.2352 |
| share Career | 0.1833 | 0.2006 | 0.2796 | 0.3554 | 0.3433 |
| share NiLF | 0.2598 | 0.2519 | 0.2162 | 0.2260 | 0.2308 |
| cons drop at H job loss rec (%) | -7.1846 | -7.0924 | -6.9749 | -7.1119 | -7.0844 |
| residual (|ΔE|+|Δgap|) | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0001 |

## 4. Mechanism counterfactuals (baseline parameters)

| counterfactual | quit exp | quit rec | quit gap (pts) | ΔE/pop rec-exp (pts) | E/pop |
|---|---|---|---|---|---|
| baseline | 0.0311 | 0.0237 | -0.742 | -1.598 | 0.664 |
| acyclical husband risk | 0.0315 | 0.0254 | -0.616 | -1.998 | 0.661 |
| acyclical job finding | 0.0318 | 0.0271 | -0.470 | -0.189 | 0.666 |
| no recession wage cut | 0.0304 | 0.0227 | -0.771 | -0.309 | 0.673 |
| acyclical own job loss | 0.0310 | 0.0234 | -0.755 | -0.911 | 0.668 |

Decomposition of the baseline quit gap (share removed when each channel is switched off):

| channel | quit gap without it | share of baseline gap |
|---|---|---|
| acyclical husband risk | -0.616 | +17% |
| acyclical job finding | -0.470 | +37% |
| no recession wage cut | -0.771 | -4% |
| acyclical own job loss | -0.755 | -2% |

## 5. Robustness (calibrated parameters held fixed)

Source: `output/robustness_final.json`.

| variant | E/pop | quit gap | ΔE/pop rec-exp | acyclical husband risk: quit gap | RoE experiment: E/pop | RoE: quit gap |
|---|---|---|---|---|---|---|
| baseline | 0.664 | -0.742 | -1.598 | -0.616 | 0.728 | -0.524 |
| no assets (a_max 0.01, 5 points) | 0.665 | -0.786 | -1.645 | -0.590 | 0.728 | -0.518 |
| asset grid 40 points | 0.670 | -0.718 | -1.558 | -0.620 | 0.734 | -0.491 |
| asset grid 40 points, a_max 30 | 0.667 | -0.752 | -1.485 | -0.633 | 0.731 | -0.509 |
| hours grid 40 points | 0.665 | -0.744 | -1.595 | -0.606 | 0.729 | -0.526 |
| hours grid 40, h_min 0.025 | 0.664 | -0.747 | -1.612 | -0.605 | 0.727 | -0.543 |
| U threshold s_bar 0.10 | 0.664 | -0.742 | -1.598 | -0.616 | 0.728 | -0.524 |
| U threshold s_bar 0.50 | 0.664 | -0.742 | -1.598 | -0.616 | 0.728 | -0.524 |
| phi_rec_H = 1 (no recession cut in husband income) | 0.658 | -0.614 | -1.783 | -0.523 | 0.721 | -0.456 |

## 6. Summary of findings

* Calibration: employment 0.664 (target 0.62), hours 0.414 (0.40), monthly quit rate 0.0311 in expansions and 0.0237 in recessions (targets 0.034 / 0.028), recession employment drop -1.60 points (-1.7), wage gap 0.748 (0.71); career shares life-cycle 0.30, part-time 0.25, career 0.18, NiLF 0.26 (0.31 / 0.28 / 0.19 / 0.22). Untargeted: unemployment rate 0.046, wife's income share 0.337, consumption falls 5.0% at the husband's job loss in expansions and 7.2% in recessions.
* Quits are pro-cyclical: the monthly quit rate falls by 24% in recessions (-0.74 points). Decomposition: making the husband's job-loss risk acyclical removes +17% of the drop, making job finding acyclical removes +37%, removing the recession wage cut changes it by -4% (the wage cut works against the insurance motive), and making the wife's own job loss acyclical changes it by -2%.
* Recession employment drop -1.60 points in the baseline; -2.00 without cyclical husband risk (precautionary labor supply offsets -0.40 points), -0.19 without the fall in job finding, -0.31 without the wage cut, -0.91 without cyclical own job loss.
* Trend to cycle, each force sized to the 1970s employment rate: RoE x1.38: employment 0.728, recession drop -1.56 points (baseline -1.60), quit gap -0.52 (baseline -0.74), career shares LC/PT/career/NiLF 0.29/0.18/0.34/0.19; comp. wage gap x1.08: employment 0.734, recession drop -1.37 points (baseline -1.60), quit gap -0.52 (baseline -0.74), career shares LC/PT/career/NiLF 0.32/0.21/0.30/0.16; cost x0.10: employment 0.726, recession drop -1.95 points (baseline -1.60), quit gap -0.62 (baseline -0.74), career shares LC/PT/career/NiLF 0.15/0.34/0.29/0.21.
* Cohort accounting with the data's wage-gap and returns-to-experience paths (household income compensated): the residual cost scale is 1940: x1.00, 1950: x1.77, 1960: x1.77, 1970: x1.77, 1980: x2.00; the recession employment drop goes from -1.60 to -1.05 points and the expansion quit rate from 0.0311 to 0.0183. Caveat: tau_w is scaled by the raw data ratio, so the measured wage gap in the model rises to 0.98 by the last cohort (data 0.77); the next refinement is to solve tau_w per cohort to hit the measured gap jointly with the cost residual.

## 7. What is fragile

* The never-working (NiLF) share is the least well fitted target; it depends on the home-production curvature in productivity (α_h) and the hours scaling of the fixed cost.
* The experience cap e_max is calibrated; the wage gap among employed wives is largely the experience premium at the cap, so the returns-to-experience experiment interacts with it.
* The transitory cost shock (sd σ_κ) drives the monthly quit rate; its distribution is not disciplined by micro data beyond the quit and exit rates.
* Career shares are computed on annual hours over ages 25-54 from the model's 4,000-hour endowment; the data taxonomy uses reported annual hours.
