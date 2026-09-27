# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v9_full.json`; returns-to-experience scale x1.59 from `output/final_results_v9.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6475 | 0.6254 | 0.6613 | 0.6569 | 0.6477 | 0.6509 | 0.6475 | 0.6475 | 0.6445 |
| hours|E | 0.4227 | 0.4166 | 0.4238 | 0.4264 | 0.4244 | 0.4200 | 0.4227 | 0.4227 | 0.4202 |
| U rate | 0.0416 | 0.0429 | 0.0417 | 0.0415 | 0.0416 | 0.0412 | 0.1220 | 0.0188 | 0.0419 |
| quit/m exp | 0.0338 | 0.0358 | 0.0322 | 0.0326 | 0.0338 | 0.0351 | 0.0338 | 0.0338 | 0.0341 |
| quit/m rec | 0.0254 | 0.0265 | 0.0243 | 0.0244 | 0.0254 | 0.0263 | 0.0254 | 0.0254 | 0.0269 |
| E->nonE/m exp | 0.0535 | 0.0555 | 0.0519 | 0.0523 | 0.0535 | 0.0548 | 0.0535 | 0.0535 | 0.0538 |
| E->nonE/m rec | 0.0439 | 0.0450 | 0.0428 | 0.0429 | 0.0439 | 0.0448 | 0.0439 | 0.0439 | 0.0454 |
| dE/pop rec-exp (pts) | -1.7174 | -1.7817 | -1.7401 | -1.6182 | -1.7405 | -1.9129 | -1.7174 | -1.7174 | -1.4956 |
| wage gap (hourly ratio) | 0.7354 | 0.7338 | 0.7327 | 0.7338 | 0.7355 | 0.7324 | 0.7354 | 0.7354 | 0.7206 |
| wife share exp | 0.3245 | 0.3121 | 0.3283 | 0.3295 | 0.3255 | 0.3245 | 0.3245 | 0.3245 | 0.3232 |
| share Lifecycle | 0.2898 | 0.2996 | 0.2850 | 0.2885 | 0.2900 | 0.2883 | 0.2898 | 0.2898 | 0.2827 |
| share PT | 0.2335 | 0.2179 | 0.2523 | 0.2460 | 0.2315 | 0.2313 | 0.2335 | 0.2335 | 0.2377 |
| share Career | 0.1862 | 0.1692 | 0.1902 | 0.1894 | 0.1879 | 0.1883 | 0.1862 | 0.1862 | 0.1848 |
| share NiLF | 0.2904 | 0.3133 | 0.2725 | 0.2760 | 0.2906 | 0.2921 | 0.2904 | 0.2904 | 0.2948 |
| cons drop at H job loss rec (%) | -9.3117 | -10.9330 | -8.4917 | -8.7339 | -9.3008 | -9.2933 | -9.3117 | -9.3117 | -9.3882 |
| mean assets/monthly HH inc | 1.8103 | 0.0058 | 3.0505 | 2.9836 | 1.8074 | 1.8114 | 1.8103 | 1.8103 | 1.8202 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.840 | -1.717 | -0.812 | -1.771 | 0.730 | -0.549 | -1.526 |
| no assets (a_max 0.01, 5 points) | -0.936 | -1.782 | -0.889 | -2.053 | 0.708 | -0.649 | -1.584 |
| asset grid 40 points | -0.794 | -1.740 | -0.760 | -1.780 | 0.741 | -0.514 | -1.570 |
| asset grid 40 points, a_max 30 | -0.819 | -1.618 | -0.783 | -1.728 | 0.738 | -0.526 | -1.548 |
| hours grid 40 points | -0.840 | -1.741 | -0.816 | -1.763 | 0.730 | -0.550 | -1.513 |
| hours grid 40, h_min 0.025 | -0.880 | -1.913 | -0.849 | -1.968 | 0.731 | -0.607 | -1.675 |
| U threshold s_bar 0.10 | -0.840 | -1.717 | -0.812 | -1.771 | 0.730 | -0.549 | -1.526 |
| U threshold s_bar 0.50 | -0.840 | -1.717 | -0.812 | -1.771 | 0.730 | -0.549 | -1.526 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.724 | -1.496 | -0.692 | -1.542 | 0.727 | -0.498 | -1.334 |

Elapsed 2389s.
