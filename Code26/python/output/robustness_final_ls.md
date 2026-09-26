# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_ls_full.json`; returns-to-experience scale x1.40 from `output/final_results_ls.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6641 | 0.6668 | 0.6695 | 0.6676 | 0.6649 | 0.6630 | 0.6641 | 0.6641 | 0.6577 |
| hours|E | 0.4097 | 0.3997 | 0.4171 | 0.4156 | 0.4095 | 0.4079 | 0.4097 | 0.4097 | 0.4051 |
| U rate | 0.0447 | 0.0480 | 0.0466 | 0.0456 | 0.0447 | 0.0443 | 0.1557 | 0.0182 | 0.0455 |
| quit/m exp | 0.0332 | 0.0316 | 0.0323 | 0.0326 | 0.0330 | 0.0336 | 0.0332 | 0.0332 | 0.0343 |
| quit/m rec | 0.0250 | 0.0230 | 0.0245 | 0.0246 | 0.0249 | 0.0253 | 0.0250 | 0.0250 | 0.0277 |
| E->nonE/m exp | 0.0505 | 0.0489 | 0.0496 | 0.0499 | 0.0503 | 0.0509 | 0.0505 | 0.0505 | 0.0516 |
| E->nonE/m rec | 0.0447 | 0.0428 | 0.0443 | 0.0444 | 0.0447 | 0.0451 | 0.0447 | 0.0447 | 0.0474 |
| dE/pop rec-exp (pts) | -1.6887 | -1.8162 | -1.6970 | -1.7242 | -1.7099 | -1.6619 | -1.6887 | -1.6887 | -1.8971 |
| wage gap (hourly ratio) | 0.7382 | 0.7315 | 0.7386 | 0.7388 | 0.7382 | 0.7372 | 0.7382 | 0.7382 | 0.7212 |
| wife share exp | 0.3319 | 0.3243 | 0.3375 | 0.3370 | 0.3316 | 0.3298 | 0.3319 | 0.3319 | 0.3293 |
| share Lifecycle | 0.3085 | 0.2960 | 0.3169 | 0.3106 | 0.3075 | 0.3031 | 0.3085 | 0.3085 | 0.2994 |
| share PT | 0.2502 | 0.2931 | 0.2465 | 0.2483 | 0.2483 | 0.2517 | 0.2502 | 0.2502 | 0.2554 |
| share Career | 0.1719 | 0.1598 | 0.1758 | 0.1769 | 0.1767 | 0.1700 | 0.1719 | 0.1719 | 0.1644 |
| share NiLF | 0.2694 | 0.2510 | 0.2608 | 0.2642 | 0.2675 | 0.2752 | 0.2694 | 0.2694 | 0.2808 |
| cons drop at H job loss rec (%) | -7.2569 | -7.3313 | -6.5907 | -6.9879 | -7.2083 | -7.2201 | -7.2569 | -7.2569 | -7.2577 |
| mean assets/monthly HH inc | 1.4074 | 0.0078 | 2.6759 | 2.3812 | 1.3987 | 1.3986 | 1.4074 | 1.4074 | 1.3970 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.823 | -1.689 | -0.669 | -2.228 | 0.731 | -0.574 | -1.680 |
| no assets (a_max 0.01, 5 points) | -0.856 | -1.816 | -0.650 | -2.503 | 0.732 | -0.550 | -2.080 |
| asset grid 40 points | -0.787 | -1.697 | -0.664 | -2.089 | 0.737 | -0.530 | -1.738 |
| asset grid 40 points, a_max 30 | -0.807 | -1.724 | -0.693 | -2.127 | 0.734 | -0.567 | -1.652 |
| hours grid 40 points | -0.811 | -1.710 | -0.673 | -2.202 | 0.731 | -0.574 | -1.689 |
| hours grid 40, h_min 0.025 | -0.831 | -1.662 | -0.671 | -2.333 | 0.729 | -0.584 | -1.738 |
| U threshold s_bar 0.10 | -0.823 | -1.689 | -0.669 | -2.228 | 0.731 | -0.574 | -1.680 |
| U threshold s_bar 0.50 | -0.823 | -1.689 | -0.669 | -2.228 | 0.731 | -0.574 | -1.680 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.662 | -1.897 | -0.564 | -2.271 | 0.723 | -0.487 | -1.835 |

Elapsed 2390s.
