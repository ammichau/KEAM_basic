# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6791 |
| hours|E | 0.4000 | 0.4096 |
| U rate | nan | 0.0482 |
| quit/m exp | 0.0340 | 0.0353 |
| quit/m rec | 0.0280 | 0.0242 |
| E->nonE/m exp | 0.0500 | 0.0527 |
| E->nonE/m rec | 0.0480 | 0.0497 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6702 |
| wife share exp | nan | 0.3210 |
| wife share rec | nan | 0.3395 |
| wage gap (hourly ratio) | 0.7100 | 0.7371 |
| share Lifecycle | 0.3100 | 0.2952 |
| share PT | 0.2800 | 0.2681 |
| share Career | 0.1900 | 0.1794 |
| share NiLF | 0.2200 | 0.2573 |
| HH income rec/exp - 1 (%) | nan | -2.4501 |
| cons drop at H job loss exp (%) | nan | -6.9240 |
| cons drop at H job loss rec (%) | nan | -9.0886 |
| mean assets/monthly HH inc | nan | 1.9149 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.28 | comp. wage gap x1.07 | cost x0.52 |
|---|---|---|---|---|
| E/pop | 0.6791 | 0.7289 | 0.7318 | 0.7292 |
| hours|E | 0.4096 | 0.4328 | 0.4310 | 0.4149 |
| U rate | 0.0482 | 0.0504 | 0.0502 | 0.0489 |
| quit/m exp | 0.0353 | 0.0260 | 0.0250 | 0.0278 |
| quit/m rec | 0.0242 | 0.0178 | 0.0166 | 0.0186 |
| E->nonE/m exp | 0.0527 | 0.0434 | 0.0424 | 0.0452 |
| E->nonE/m rec | 0.0497 | 0.0433 | 0.0423 | 0.0443 |
| dE/pop rec-exp (pts) | -1.6702 | -1.9507 | -1.7709 | -2.1875 |
| wife share exp | 0.3210 | 0.3760 | 0.3698 | 0.3449 |
| wife share rec | 0.3395 | 0.3948 | 0.3883 | 0.3627 |
| wage gap (hourly ratio) | 0.7371 | 0.8318 | 0.8156 | 0.7529 |
| share Lifecycle | 0.2952 | 0.2950 | 0.3150 | 0.2313 |
| share PT | 0.2681 | 0.2306 | 0.2506 | 0.2996 |
| share Career | 0.1794 | 0.2885 | 0.2646 | 0.2631 |
| share NiLF | 0.2573 | 0.1858 | 0.1698 | 0.2060 |
| HH income rec/exp - 1 (%) | -2.4501 | -2.1330 | -2.1972 | -2.5605 |
| cons drop at H job loss exp (%) | -6.9240 | -6.2589 | -6.2995 | -6.7194 |
| cons drop at H job loss rec (%) | -9.0886 | -8.2717 | -8.3701 | -8.8531 |
| mean assets/monthly HH inc | 1.9149 | 2.0960 | 1.9670 | 1.9319 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.77) | 1960 (cost x1.77) | 1970 (cost x1.77) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6791 | 0.6684 | 0.7107 | 0.7308 | 0.7266 |
| hours|E | 0.4096 | 0.4256 | 0.4447 | 0.4528 | 0.4572 |
| U rate | 0.0482 | 0.0497 | 0.0510 | 0.0516 | 0.0519 |
| quit/m exp | 0.0353 | 0.0325 | 0.0247 | 0.0216 | 0.0208 |
| quit/m rec | 0.0242 | 0.0217 | 0.0163 | 0.0144 | 0.0139 |
| E->nonE/m exp | 0.0527 | 0.0499 | 0.0421 | 0.0390 | 0.0382 |
| E->nonE/m rec | 0.0497 | 0.0473 | 0.0419 | 0.0401 | 0.0395 |
| dE/pop rec-exp (pts) | -1.6702 | -1.3213 | -1.4674 | -1.7145 | -1.5345 |
| wife share exp | 0.3210 | 0.3437 | 0.3870 | 0.4124 | 0.4188 |
| wife share rec | 0.3395 | 0.3632 | 0.4064 | 0.4313 | 0.4379 |
| wage gap (hourly ratio) | 0.7371 | 0.8085 | 0.8834 | 0.9367 | 0.9612 |
| share Lifecycle | 0.2952 | 0.3519 | 0.3850 | 0.3837 | 0.4000 |
| share PT | 0.2681 | 0.2154 | 0.1837 | 0.1638 | 0.1527 |
| share Career | 0.1794 | 0.1842 | 0.2523 | 0.3023 | 0.2985 |
| share NiLF | 0.2573 | 0.2485 | 0.1790 | 0.1502 | 0.1487 |
| HH income rec/exp - 1 (%) | -2.4501 | -2.0525 | -1.8220 | -1.7665 | -1.6347 |
| cons drop at H job loss exp (%) | -6.9240 | -6.5503 | -6.0532 | -5.7695 | -5.7141 |
| cons drop at H job loss rec (%) | -9.0886 | -8.6558 | -8.0941 | -7.7892 | -7.7269 |
| mean assets/monthly HH inc | 1.9149 | 1.9968 | 2.0751 | 2.1662 | 2.1812 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6791 | 0.6672 | 0.6800 | 0.6791 | 0.6882 |
| hours|E | 0.4096 | 0.4048 | 0.4089 | 0.4096 | 0.4092 |
| U rate | 0.0482 | 0.0485 | 0.0456 | 0.0482 | 0.0437 |
| quit/m exp | 0.0353 | 0.0371 | 0.0362 | 0.0353 | 0.0347 |
| quit/m rec | 0.0242 | 0.0284 | 0.0314 | 0.0242 | 0.0228 |
| E->nonE/m exp | 0.0527 | 0.0545 | 0.0537 | 0.0527 | 0.0520 |
| E->nonE/m rec | 0.0497 | 0.0540 | 0.0571 | 0.0497 | 0.0397 |
| dE/pop rec-exp (pts) | -1.6702 | -2.8139 | 0.1942 | -1.6702 | 0.7218 |
| wife share exp | 0.3210 | 0.3156 | 0.3205 | 0.3210 | 0.3227 |
| wife share rec | 0.3395 | 0.3136 | 0.3425 | 0.3395 | 0.3463 |
| wage gap (hourly ratio) | 0.7371 | 0.7368 | 0.7373 | 0.7371 | 0.7378 |
| share Lifecycle | 0.2952 | 0.2854 | 0.3094 | 0.2952 | 0.2979 |
| share PT | 0.2681 | 0.2712 | 0.2485 | 0.2681 | 0.2654 |
| share Career | 0.1794 | 0.1685 | 0.1787 | 0.1794 | 0.1852 |
| share NiLF | 0.2573 | 0.2748 | 0.2633 | 0.2573 | 0.2515 |
| HH income rec/exp - 1 (%) | -2.4501 | -0.1512 | -1.8968 | -2.4501 | -1.6378 |
| cons drop at H job loss exp (%) | -6.9240 | -7.2373 | -6.9347 | -6.9240 | -6.9399 |
| cons drop at H job loss rec (%) | -9.0886 | -7.4644 | -8.9068 | -9.0886 | -8.9486 |
| mean assets/monthly HH inc | 1.9149 | 1.8411 | 1.9128 | 1.9149 | 1.9060 |


Elapsed 1657s. Calibration file: output/final_calib_v4nb_full.json
