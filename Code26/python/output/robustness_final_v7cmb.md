# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v7cmb_full.json`; returns-to-experience scale x1.35 from `output/final_results_v7cmb.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6627 | 0.6687 | 0.6519 | 0.6430 | 0.6634 | 0.6672 | 0.6627 | 0.6627 | 0.6614 |
| hours|E | 0.4228 | 0.4191 | 0.4213 | 0.4136 | 0.4229 | 0.4167 | 0.4228 | 0.4228 | 0.4198 |
| U rate | 0.0354 | 0.0357 | 0.0355 | 0.0351 | 0.0353 | 0.0348 | 0.0483 | 0.0110 | 0.0352 |
| quit/m exp | 0.0235 | 0.0226 | 0.0249 | 0.0262 | 0.0235 | 0.0253 | 0.0235 | 0.0235 | 0.0241 |
| quit/m rec | 0.0173 | 0.0157 | 0.0190 | 0.0203 | 0.0172 | 0.0192 | 0.0173 | 0.0173 | 0.0189 |
| E->nonE/m exp | 0.0368 | 0.0359 | 0.0382 | 0.0395 | 0.0367 | 0.0386 | 0.0368 | 0.0368 | 0.0374 |
| E->nonE/m rec | 0.0325 | 0.0310 | 0.0342 | 0.0356 | 0.0324 | 0.0345 | 0.0325 | 0.0325 | 0.0341 |
| dE/pop rec-exp (pts) | -1.6247 | -1.8513 | -1.5836 | -1.0933 | -1.5869 | -1.5496 | -1.6247 | -1.6247 | -1.1235 |
| wage gap (hourly ratio) | 0.7052 | 0.6991 | 0.7111 | 0.7080 | 0.7051 | 0.7017 | 0.7052 | 0.7052 | 0.6909 |
| wife share exp | 0.3227 | 0.3187 | 0.3219 | 0.3132 | 0.3228 | 0.3215 | 0.3227 | 0.3227 | 0.3217 |
| share Lifecycle | 0.3017 | 0.3160 | 0.2696 | 0.2683 | 0.2994 | 0.2946 | 0.3017 | 0.3017 | 0.2965 |
| share PT | 0.2665 | 0.2873 | 0.2610 | 0.2687 | 0.2698 | 0.2721 | 0.2665 | 0.2665 | 0.2696 |
| share Career | 0.1892 | 0.1673 | 0.2025 | 0.1840 | 0.1890 | 0.1867 | 0.1892 | 0.1892 | 0.1856 |
| share NiLF | 0.2427 | 0.2294 | 0.2669 | 0.2790 | 0.2419 | 0.2467 | 0.2427 | 0.2427 | 0.2483 |
| cons drop at H job loss rec (%) | -10.1240 | -11.0233 | -9.8225 | -9.3372 | -10.0968 | -10.1034 | -10.1240 | -10.1240 | -10.7643 |
| mean assets/monthly HH inc | 5.7364 | 0.0073 | 15.1357 | 19.8992 | 5.7326 | 5.7494 | 5.7364 | 5.7364 | 5.6783 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.625 | -1.625 | -0.587 | -1.568 | 0.728 | -0.446 | -1.670 |
| no assets (a_max 0.01, 5 points) | -0.689 | -1.851 | -0.648 | -1.733 | 0.729 | -0.517 | -1.905 |
| asset grid 40 points | -0.593 | -1.584 | -0.514 | -1.343 | 0.720 | -0.404 | -1.478 |
| asset grid 40 points, a_max 30 | -0.592 | -1.093 | -0.521 | -1.153 | 0.710 | -0.387 | -1.315 |
| hours grid 40 points | -0.626 | -1.587 | -0.582 | -1.576 | 0.728 | -0.442 | -1.696 |
| hours grid 40, h_min 0.025 | -0.612 | -1.550 | -0.569 | -1.487 | 0.732 | -0.459 | -1.663 |
| U threshold s_bar 0.10 | -0.625 | -1.625 | -0.587 | -1.568 | 0.728 | -0.446 | -1.670 |
| U threshold s_bar 0.50 | -0.625 | -1.625 | -0.587 | -1.568 | 0.728 | -0.446 | -1.670 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.525 | -1.124 | -0.481 | -1.167 | 0.726 | -0.374 | -1.340 |

Elapsed 2784s.
