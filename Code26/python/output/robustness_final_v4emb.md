# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v4emb_full.json`; returns-to-experience scale x1.07 from `output/final_results_v4emb.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.7172 | 0.7091 | 0.7261 | 0.7275 | 0.7174 | 0.7183 | 0.7172 | 0.7172 | 0.7101 |
| hours|E | 0.4035 | 0.3891 | 0.4145 | 0.4105 | 0.4040 | 0.4001 | 0.4035 | 0.4035 | 0.3987 |
| U rate | 0.0393 | 0.0400 | 0.0409 | 0.0408 | 0.0393 | 0.0390 | 0.1388 | 0.0099 | 0.0391 |
| quit/m exp | 0.0267 | 0.0271 | 0.0257 | 0.0256 | 0.0267 | 0.0278 | 0.0267 | 0.0267 | 0.0280 |
| quit/m rec | 0.0181 | 0.0176 | 0.0177 | 0.0179 | 0.0180 | 0.0190 | 0.0181 | 0.0181 | 0.0207 |
| E->nonE/m exp | 0.0400 | 0.0404 | 0.0390 | 0.0389 | 0.0400 | 0.0410 | 0.0400 | 0.0400 | 0.0413 |
| E->nonE/m rec | 0.0333 | 0.0329 | 0.0330 | 0.0332 | 0.0333 | 0.0342 | 0.0333 | 0.0333 | 0.0359 |
| dE/pop rec-exp (pts) | -1.6454 | -1.3860 | -1.8169 | -1.8421 | -1.6600 | -1.7176 | -1.6454 | -1.6454 | -1.6527 |
| wage gap (hourly ratio) | 0.7429 | 0.7293 | 0.7478 | 0.7435 | 0.7430 | 0.7409 | 0.7429 | 0.7429 | 0.7281 |
| wife share exp | 0.3333 | 0.3162 | 0.3452 | 0.3406 | 0.3339 | 0.3324 | 0.3333 | 0.3333 | 0.3307 |
| share Lifecycle | 0.2963 | 0.2804 | 0.3046 | 0.3227 | 0.2925 | 0.2894 | 0.2963 | 0.2963 | 0.2873 |
| share PT | 0.2896 | 0.3460 | 0.2554 | 0.2452 | 0.2854 | 0.2842 | 0.2896 | 0.2896 | 0.2990 |
| share Career | 0.2029 | 0.1640 | 0.2387 | 0.2367 | 0.2106 | 0.2102 | 0.2029 | 0.2029 | 0.1888 |
| share NiLF | 0.2112 | 0.2096 | 0.2013 | 0.1954 | 0.2115 | 0.2162 | 0.2112 | 0.2112 | 0.2250 |
| cons drop at H job loss rec (%) | -9.9832 | -10.4758 | -8.4403 | -7.5979 | -9.9563 | -9.9683 | -9.9832 | -9.9832 | -10.2695 |
| mean assets/monthly HH inc | 3.6205 | 0.0073 | 9.0475 | 8.6092 | 3.6257 | 3.6201 | 3.6205 | 3.6205 | 3.5447 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.863 | -1.645 | -0.836 | -2.011 | 0.730 | -0.802 | -1.669 |
| no assets (a_max 0.01, 5 points) | -0.952 | -1.386 | -0.847 | -2.050 | 0.724 | -0.872 | -1.375 |
| asset grid 40 points | -0.797 | -1.817 | -0.756 | -2.156 | 0.741 | -0.737 | -1.767 |
| asset grid 40 points, a_max 30 | -0.768 | -1.842 | -0.736 | -2.202 | 0.743 | -0.715 | -1.859 |
| hours grid 40 points | -0.865 | -1.660 | -0.825 | -2.038 | 0.731 | -0.802 | -1.664 |
| hours grid 40, h_min 0.025 | -0.878 | -1.718 | -0.833 | -2.074 | 0.732 | -0.813 | -1.729 |
| U threshold s_bar 0.10 | -0.863 | -1.645 | -0.836 | -2.011 | 0.730 | -0.802 | -1.669 |
| U threshold s_bar 0.50 | -0.863 | -1.645 | -0.836 | -2.011 | 0.730 | -0.802 | -1.669 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.734 | -1.653 | -0.648 | -2.153 | 0.725 | -0.670 | -1.722 |

Elapsed 2689s.
