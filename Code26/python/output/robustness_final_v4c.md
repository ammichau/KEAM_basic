# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v4c_full.json`; returns-to-experience scale x1.30 from `output/final_results_v4c.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6737 | 0.6776 | 0.6773 | 0.6768 | 0.6741 | 0.6728 | 0.6737 | 0.6737 | 0.6657 |
| hours|E | 0.4077 | 0.3972 | 0.4147 | 0.4131 | 0.4079 | 0.4058 | 0.4077 | 0.4077 | 0.4030 |
| U rate | 0.0441 | 0.0484 | 0.0460 | 0.0445 | 0.0441 | 0.0436 | 0.1608 | 0.0187 | 0.0448 |
| quit/m exp | 0.0343 | 0.0326 | 0.0338 | 0.0340 | 0.0343 | 0.0348 | 0.0343 | 0.0343 | 0.0356 |
| quit/m rec | 0.0227 | 0.0203 | 0.0231 | 0.0229 | 0.0227 | 0.0230 | 0.0227 | 0.0227 | 0.0265 |
| E->nonE/m exp | 0.0530 | 0.0513 | 0.0524 | 0.0526 | 0.0529 | 0.0535 | 0.0530 | 0.0530 | 0.0542 |
| E->nonE/m rec | 0.0440 | 0.0415 | 0.0444 | 0.0442 | 0.0440 | 0.0443 | 0.0440 | 0.0440 | 0.0477 |
| dE/pop rec-exp (pts) | -1.6494 | -1.4704 | -1.8728 | -1.7907 | -1.6494 | -1.6660 | -1.6494 | -1.6494 | -2.4526 |
| wage gap (hourly ratio) | 0.7308 | 0.7230 | 0.7319 | 0.7317 | 0.7308 | 0.7295 | 0.7308 | 0.7308 | 0.7141 |
| wife share exp | 0.3316 | 0.3231 | 0.3366 | 0.3363 | 0.3317 | 0.3299 | 0.3316 | 0.3316 | 0.3286 |
| share Lifecycle | 0.3085 | 0.2979 | 0.3181 | 0.3131 | 0.3104 | 0.3071 | 0.3085 | 0.3085 | 0.3063 |
| share PT | 0.2590 | 0.3069 | 0.2544 | 0.2560 | 0.2519 | 0.2548 | 0.2590 | 0.2590 | 0.2590 |
| share Career | 0.1673 | 0.1500 | 0.1690 | 0.1702 | 0.1735 | 0.1673 | 0.1673 | 0.1673 | 0.1588 |
| share NiLF | 0.2652 | 0.2452 | 0.2585 | 0.2606 | 0.2642 | 0.2708 | 0.2652 | 0.2652 | 0.2760 |
| cons drop at H job loss rec (%) | -8.5594 | -8.8833 | -7.5864 | -8.1973 | -8.5351 | -8.5622 | -8.5594 | -8.5594 | -8.4124 |
| mean assets/monthly HH inc | 1.4538 | 0.0078 | 2.7123 | 2.4177 | 1.4566 | 1.4542 | 1.4538 | 1.4538 | 1.4627 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -1.161 | -1.649 | -0.884 | -2.776 | 0.729 | -0.821 | -1.948 |
| no assets (a_max 0.01, 5 points) | -1.238 | -1.470 | -0.942 | -2.565 | 0.733 | -0.879 | -1.616 |
| asset grid 40 points | -1.076 | -1.873 | -0.876 | -2.494 | 0.733 | -0.752 | -2.010 |
| asset grid 40 points, a_max 30 | -1.114 | -1.791 | -0.902 | -2.563 | 0.731 | -0.797 | -1.888 |
| hours grid 40 points | -1.158 | -1.649 | -0.894 | -2.714 | 0.729 | -0.813 | -1.949 |
| hours grid 40, h_min 0.025 | -1.179 | -1.666 | -0.905 | -2.824 | 0.726 | -0.830 | -2.046 |
| U threshold s_bar 0.10 | -1.161 | -1.649 | -0.884 | -2.776 | 0.729 | -0.821 | -1.948 |
| U threshold s_bar 0.50 | -1.161 | -1.649 | -0.884 | -2.776 | 0.729 | -0.821 | -1.948 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.909 | -2.453 | -0.725 | -3.174 | 0.720 | -0.692 | -2.335 |

Elapsed 1374s.
