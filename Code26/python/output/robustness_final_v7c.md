# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v7c_full.json`; returns-to-experience scale x1.40 from `output/final_results_v7c.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6763 | 0.6744 | 0.6810 | 0.6794 | 0.6769 | 0.6760 | 0.6763 | 0.6763 | 0.6731 |
| hours|E | 0.4134 | 0.4074 | 0.4172 | 0.4173 | 0.4136 | 0.4104 | 0.4134 | 0.4134 | 0.4112 |
| U rate | 0.0397 | 0.0415 | 0.0405 | 0.0397 | 0.0397 | 0.0391 | 0.0729 | 0.0142 | 0.0400 |
| quit/m exp | 0.0333 | 0.0327 | 0.0329 | 0.0331 | 0.0333 | 0.0339 | 0.0333 | 0.0333 | 0.0339 |
| quit/m rec | 0.0253 | 0.0248 | 0.0252 | 0.0254 | 0.0253 | 0.0260 | 0.0253 | 0.0253 | 0.0266 |
| E->nonE/m exp | 0.0531 | 0.0525 | 0.0526 | 0.0529 | 0.0530 | 0.0536 | 0.0531 | 0.0531 | 0.0536 |
| E->nonE/m rec | 0.0457 | 0.0452 | 0.0455 | 0.0456 | 0.0456 | 0.0463 | 0.0457 | 0.0457 | 0.0469 |
| dE/pop rec-exp (pts) | -1.6658 | -1.5867 | -1.6794 | -1.6660 | -1.6302 | -1.7158 | -1.6658 | -1.6658 | -1.8506 |
| wage gap (hourly ratio) | 0.7118 | 0.7081 | 0.7128 | 0.7128 | 0.7121 | 0.7083 | 0.7118 | 0.7118 | 0.6976 |
| wife share exp | 0.3273 | 0.3215 | 0.3309 | 0.3311 | 0.3275 | 0.3258 | 0.3273 | 0.3273 | 0.3260 |
| share Lifecycle | 0.3056 | 0.3004 | 0.3108 | 0.3081 | 0.3094 | 0.2994 | 0.3056 | 0.3056 | 0.3046 |
| share PT | 0.2727 | 0.2906 | 0.2677 | 0.2696 | 0.2687 | 0.2698 | 0.2727 | 0.2727 | 0.2717 |
| share Career | 0.1748 | 0.1675 | 0.1781 | 0.1777 | 0.1758 | 0.1783 | 0.1748 | 0.1748 | 0.1719 |
| share NiLF | 0.2469 | 0.2415 | 0.2433 | 0.2446 | 0.2460 | 0.2525 | 0.2469 | 0.2469 | 0.2519 |
| cons drop at H job loss rec (%) | -10.0996 | -10.6967 | -8.9445 | -9.5544 | -10.1081 | -10.1606 | -10.0996 | -10.0996 | -10.0769 |
| mean assets/monthly HH inc | 1.6206 | 0.0072 | 2.8838 | 2.6268 | 1.6221 | 1.6161 | 1.6206 | 1.6206 | 1.6222 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.802 | -1.666 | -0.688 | -1.668 | 0.732 | -0.613 | -1.620 |
| no assets (a_max 0.01, 5 points) | -0.792 | -1.587 | -0.661 | -1.774 | 0.727 | -0.601 | -1.713 |
| asset grid 40 points | -0.763 | -1.679 | -0.678 | -1.650 | 0.736 | -0.577 | -1.694 |
| asset grid 40 points, a_max 30 | -0.773 | -1.666 | -0.684 | -1.578 | 0.734 | -0.590 | -1.612 |
| hours grid 40 points | -0.800 | -1.630 | -0.683 | -1.677 | 0.732 | -0.614 | -1.610 |
| hours grid 40, h_min 0.025 | -0.787 | -1.716 | -0.661 | -1.755 | 0.726 | -0.608 | -1.724 |
| U threshold s_bar 0.10 | -0.802 | -1.666 | -0.688 | -1.668 | 0.732 | -0.613 | -1.620 |
| U threshold s_bar 0.50 | -0.802 | -1.666 | -0.688 | -1.668 | 0.732 | -0.613 | -1.620 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.722 | -1.851 | -0.625 | -1.815 | 0.726 | -0.580 | -1.764 |

Elapsed 2498s.
