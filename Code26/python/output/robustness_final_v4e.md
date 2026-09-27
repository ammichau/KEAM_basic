# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_v4e_full.json`; returns-to-experience scale x1.26 from `output/final_results_v4e.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6831 | 0.6880 | 0.6886 | 0.6852 | 0.6834 | 0.6827 | 0.6831 | 0.6831 | 0.6738 |
| hours|E | 0.4116 | 0.3988 | 0.4176 | 0.4194 | 0.4113 | 0.4084 | 0.4116 | 0.4116 | 0.4062 |
| U rate | 0.0456 | 0.0485 | 0.0461 | 0.0460 | 0.0457 | 0.0451 | 0.1820 | 0.0176 | 0.0458 |
| quit/m exp | 0.0354 | 0.0336 | 0.0347 | 0.0352 | 0.0354 | 0.0359 | 0.0354 | 0.0354 | 0.0370 |
| quit/m rec | 0.0235 | 0.0210 | 0.0233 | 0.0237 | 0.0233 | 0.0239 | 0.0235 | 0.0235 | 0.0284 |
| E->nonE/m exp | 0.0533 | 0.0516 | 0.0526 | 0.0531 | 0.0533 | 0.0538 | 0.0533 | 0.0533 | 0.0550 |
| E->nonE/m rec | 0.0449 | 0.0425 | 0.0448 | 0.0452 | 0.0448 | 0.0454 | 0.0449 | 0.0449 | 0.0498 |
| dE/pop rec-exp (pts) | -1.7584 | -1.3551 | -1.9076 | -1.8604 | -1.7124 | -1.7913 | -1.7584 | -1.7584 | -2.2649 |
| wage gap (hourly ratio) | 0.7299 | 0.7221 | 0.7304 | 0.7308 | 0.7296 | 0.7275 | 0.7299 | 0.7299 | 0.7147 |
| wife share exp | 0.3230 | 0.3132 | 0.3270 | 0.3289 | 0.3230 | 0.3215 | 0.3230 | 0.3230 | 0.3197 |
| share Lifecycle | 0.3025 | 0.2890 | 0.3098 | 0.3058 | 0.3008 | 0.2988 | 0.3025 | 0.3025 | 0.2973 |
| share PT | 0.2550 | 0.2942 | 0.2533 | 0.2531 | 0.2567 | 0.2485 | 0.2550 | 0.2550 | 0.2529 |
| share Career | 0.1873 | 0.1750 | 0.1881 | 0.1875 | 0.1869 | 0.1875 | 0.1873 | 0.1873 | 0.1775 |
| share NiLF | 0.2552 | 0.2419 | 0.2487 | 0.2535 | 0.2556 | 0.2652 | 0.2552 | 0.2552 | 0.2723 |
| cons drop at H job loss rec (%) | -10.0754 | -10.5712 | -9.0950 | -9.6061 | -10.0637 | -10.0902 | -10.0754 | -10.0754 | -10.2837 |
| mean assets/monthly HH inc | 1.8603 | 0.0071 | 2.9729 | 3.1241 | 1.8568 | 1.8528 | 1.8603 | 1.8603 | 1.8678 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -1.198 | -1.758 | -1.058 | -2.619 | 0.732 | -0.915 | -1.671 |
| no assets (a_max 0.01, 5 points) | -1.265 | -1.355 | -1.033 | -2.658 | 0.733 | -0.924 | -1.480 |
| asset grid 40 points | -1.141 | -1.908 | -1.016 | -2.600 | 0.736 | -0.877 | -1.842 |
| asset grid 40 points, a_max 30 | -1.146 | -1.860 | -1.033 | -2.552 | 0.734 | -0.897 | -1.750 |
| hours grid 40 points | -1.209 | -1.712 | -1.053 | -2.628 | 0.732 | -0.906 | -1.768 |
| hours grid 40, h_min 0.025 | -1.198 | -1.791 | -1.052 | -2.688 | 0.731 | -0.927 | -1.789 |
| U threshold s_bar 0.10 | -1.198 | -1.758 | -1.058 | -2.619 | 0.732 | -0.915 | -1.671 |
| U threshold s_bar 0.50 | -1.198 | -1.758 | -1.058 | -2.619 | 0.732 | -0.915 | -1.671 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.869 | -2.265 | -0.722 | -2.903 | 0.723 | -0.682 | -2.219 |

Elapsed 1621s.
