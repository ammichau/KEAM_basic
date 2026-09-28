# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v4nb_full.json`; returns-to-experience scale x1.28 from `output/final_results_v4nb.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6791 | 0.6854 | 0.6853 | 0.6817 | 0.6794 | 0.6789 | 0.6791 | 0.6791 | 0.6791 |
| hours|E | 0.4096 | 0.3964 | 0.4158 | 0.4175 | 0.4101 | 0.4076 | 0.4096 | 0.4096 | 0.4096 |
| U rate | 0.0482 | 0.0508 | 0.0487 | 0.0484 | 0.0482 | 0.0477 | 0.1916 | 0.0182 | 0.0482 |
| quit/m exp | 0.0353 | 0.0335 | 0.0344 | 0.0350 | 0.0353 | 0.0357 | 0.0353 | 0.0353 | 0.0353 |
| quit/m rec | 0.0242 | 0.0212 | 0.0235 | 0.0239 | 0.0241 | 0.0246 | 0.0242 | 0.0242 | 0.0242 |
| E->nonE/m exp | 0.0527 | 0.0509 | 0.0519 | 0.0524 | 0.0527 | 0.0531 | 0.0527 | 0.0527 | 0.0527 |
| E->nonE/m rec | 0.0497 | 0.0469 | 0.0492 | 0.0495 | 0.0496 | 0.0502 | 0.0497 | 0.0497 | 0.0497 |
| dE/pop rec-exp (pts) | -1.6702 | -0.9428 | -1.8858 | -1.7788 | -1.6154 | -1.7296 | -1.6702 | -1.6702 | -1.6702 |
| wage gap (hourly ratio) | 0.7371 | 0.7286 | 0.7379 | 0.7387 | 0.7369 | 0.7352 | 0.7371 | 0.7371 | 0.7371 |
| wife share exp | 0.3210 | 0.3127 | 0.3251 | 0.3265 | 0.3214 | 0.3200 | 0.3210 | 0.3210 | 0.3210 |
| share Lifecycle | 0.2952 | 0.2737 | 0.3044 | 0.3017 | 0.2908 | 0.2871 | 0.2952 | 0.2952 | 0.2952 |
| share PT | 0.2681 | 0.3187 | 0.2679 | 0.2640 | 0.2733 | 0.2677 | 0.2681 | 0.2681 | 0.2681 |
| share Career | 0.1794 | 0.1656 | 0.1798 | 0.1810 | 0.1775 | 0.1790 | 0.1794 | 0.1794 | 0.1794 |
| share NiLF | 0.2573 | 0.2419 | 0.2479 | 0.2533 | 0.2583 | 0.2662 | 0.2573 | 0.2573 | 0.2573 |
| cons drop at H job loss rec (%) | -9.0886 | -10.7009 | -8.0790 | -8.4866 | -9.1352 | -9.1220 | -9.0886 | -9.0886 | -9.0886 |
| mean assets/monthly HH inc | 1.9149 | 0.0069 | 2.9956 | 3.1982 | 1.9129 | 1.9067 | 1.9149 | 1.9149 | 1.9149 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -1.107 | -1.670 | -0.922 | -2.687 | 0.729 | -0.820 | -1.951 |
| no assets (a_max 0.01, 5 points) | -1.229 | -0.943 | -0.940 | -2.470 | 0.732 | -0.906 | -1.549 |
| asset grid 40 points | -1.095 | -1.886 | -0.920 | -2.563 | 0.734 | -0.813 | -2.001 |
| asset grid 40 points, a_max 30 | -1.105 | -1.779 | -0.926 | -2.550 | 0.732 | -0.814 | -1.979 |
| hours grid 40 points | -1.122 | -1.615 | -0.925 | -2.706 | 0.729 | -0.818 | -1.970 |
| hours grid 40, h_min 0.025 | -1.104 | -1.730 | -0.945 | -2.644 | 0.728 | -0.831 | -1.981 |
| U threshold s_bar 0.10 | -1.107 | -1.670 | -0.922 | -2.687 | 0.729 | -0.820 | -1.951 |
| U threshold s_bar 0.50 | -1.107 | -1.670 | -0.922 | -2.687 | 0.729 | -0.820 | -1.951 |
| phi_rec_H = 1 (no recession cut in husband income) | -1.107 | -1.670 | -0.922 | -2.687 | 0.729 | -0.820 | -1.951 |

Elapsed 1912s.
