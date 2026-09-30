# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v9n_full.json`; returns-to-experience scale x1.35 from `output/final_results_v9n.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6816 | 0.6616 | 0.6944 | 0.6900 | 0.6820 | 0.6824 | 0.6816 | 0.6816 | 0.6816 |
| hours|E | 0.4093 | 0.4023 | 0.4103 | 0.4130 | 0.4096 | 0.4057 | 0.4093 | 0.4093 | 0.4093 |
| U rate | 0.0437 | 0.0447 | 0.0435 | 0.0437 | 0.0437 | 0.0432 | 0.1372 | 0.0189 | 0.0437 |
| quit/m exp | 0.0357 | 0.0384 | 0.0342 | 0.0345 | 0.0356 | 0.0362 | 0.0357 | 0.0357 | 0.0357 |
| quit/m rec | 0.0270 | 0.0270 | 0.0254 | 0.0261 | 0.0269 | 0.0282 | 0.0270 | 0.0270 | 0.0270 |
| E->nonE/m exp | 0.0551 | 0.0578 | 0.0536 | 0.0539 | 0.0550 | 0.0556 | 0.0551 | 0.0551 | 0.0551 |
| E->nonE/m rec | 0.0485 | 0.0486 | 0.0471 | 0.0476 | 0.0484 | 0.0496 | 0.0485 | 0.0485 | 0.0485 |
| dE/pop rec-exp (pts) | -1.6678 | -1.3797 | -1.6791 | -1.6862 | -1.6671 | -1.6945 | -1.6678 | -1.6678 | -1.6678 |
| wage gap (hourly ratio) | 0.7124 | 0.7090 | 0.7097 | 0.7110 | 0.7124 | 0.7091 | 0.7124 | 0.7124 | 0.7124 |
| wife share exp | 0.3210 | 0.3092 | 0.3245 | 0.3254 | 0.3214 | 0.3197 | 0.3210 | 0.3210 | 0.3210 |
| share Lifecycle | 0.2798 | 0.2977 | 0.2752 | 0.2771 | 0.2806 | 0.2717 | 0.2798 | 0.2798 | 0.2798 |
| share PT | 0.2815 | 0.2587 | 0.2992 | 0.2940 | 0.2800 | 0.2842 | 0.2815 | 0.2815 | 0.2815 |
| share Career | 0.1758 | 0.1592 | 0.1798 | 0.1787 | 0.1767 | 0.1756 | 0.1758 | 0.1758 | 0.1758 |
| share NiLF | 0.2629 | 0.2844 | 0.2458 | 0.2502 | 0.2627 | 0.2685 | 0.2629 | 0.2629 | 0.2629 |
| cons drop at H job loss rec (%) | -8.7096 | -11.0620 | -7.6437 | -7.9131 | -8.7201 | -8.7415 | -8.7096 | -8.7096 | -8.7096 |
| mean assets/monthly HH inc | 1.8570 | 0.0054 | 3.0870 | 3.0723 | 1.8570 | 1.8535 | 1.8570 | 1.8570 | 1.8570 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.871 | -1.668 | -0.815 | -1.938 | 0.729 | -0.668 | -1.719 |
| no assets (a_max 0.01, 5 points) | -1.137 | -1.380 | -0.993 | -1.917 | 0.708 | -0.862 | -1.562 |
| asset grid 40 points | -0.878 | -1.679 | -0.835 | -1.808 | 0.742 | -0.631 | -1.734 |
| asset grid 40 points, a_max 30 | -0.849 | -1.686 | -0.816 | -1.829 | 0.739 | -0.635 | -1.695 |
| hours grid 40 points | -0.869 | -1.667 | -0.814 | -1.936 | 0.730 | -0.668 | -1.728 |
| hours grid 40, h_min 0.025 | -0.807 | -1.694 | -0.738 | -1.903 | 0.729 | -0.642 | -1.722 |
| U threshold s_bar 0.10 | -0.871 | -1.668 | -0.815 | -1.938 | 0.729 | -0.668 | -1.719 |
| U threshold s_bar 0.50 | -0.871 | -1.668 | -0.815 | -1.938 | 0.729 | -0.668 | -1.719 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.871 | -1.668 | -0.815 | -1.938 | 0.729 | -0.668 | -1.719 |

Elapsed 3461s.
