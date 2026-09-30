# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6763 |
| hours|E | 0.4000 | 0.4134 |
| U rate | nan | 0.0397 |
| quit/m exp | 0.0340 | 0.0333 |
| quit/m rec | 0.0280 | 0.0253 |
| E->nonE/m exp | 0.0500 | 0.0531 |
| E->nonE/m rec | 0.0480 | 0.0457 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6658 |
| wife share exp | nan | 0.3273 |
| wife share rec | nan | 0.3279 |
| wage gap (hourly ratio) | 0.7100 | 0.7118 |
| share Lifecycle | 0.3100 | 0.3056 |
| share PT | 0.2800 | 0.2727 |
| share Career | 0.1900 | 0.1748 |
| share NiLF | 0.2200 | 0.2469 |
| HH income rec/exp - 1 (%) | nan | -16.3030 |
| cons drop at H job loss exp (%) | nan | -7.0656 |
| cons drop at H job loss rec (%) | nan | -10.0996 |
| mean assets/monthly HH inc | nan | 1.6206 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.40 | comp. wage gap x1.10 | cost x0.10 |
|---|---|---|---|---|
| E/pop | 0.6763 | 0.7315 | 0.7292 | 0.7144 |
| hours|E | 0.4134 | 0.4605 | 0.4472 | 0.4211 |
| U rate | 0.0397 | 0.0417 | 0.0408 | 0.0391 |
| quit/m exp | 0.0333 | 0.0245 | 0.0245 | 0.0278 |
| quit/m rec | 0.0253 | 0.0184 | 0.0184 | 0.0209 |
| E->nonE/m exp | 0.0531 | 0.0443 | 0.0443 | 0.0476 |
| E->nonE/m rec | 0.0457 | 0.0387 | 0.0387 | 0.0412 |
| dE/pop rec-exp (pts) | -1.6658 | -1.6205 | -1.4531 | -1.8722 |
| wife share exp | 0.3273 | 0.4108 | 0.3955 | 0.3487 |
| wife share rec | 0.3279 | 0.4112 | 0.3966 | 0.3475 |
| wage gap (hourly ratio) | 0.7118 | 0.8452 | 0.8317 | 0.7265 |
| share Lifecycle | 0.3056 | 0.3231 | 0.3544 | 0.2292 |
| share PT | 0.2727 | 0.1938 | 0.2210 | 0.3008 |
| share Career | 0.1748 | 0.3010 | 0.2577 | 0.2569 |
| share NiLF | 0.2469 | 0.1821 | 0.1669 | 0.2131 |
| HH income rec/exp - 1 (%) | -16.3030 | -16.2625 | -16.1226 | -16.6391 |
| cons drop at H job loss exp (%) | -7.0656 | -6.4005 | -6.4452 | -6.9114 |
| cons drop at H job loss rec (%) | -10.0996 | -9.4733 | -9.4401 | -9.9961 |
| mean assets/monthly HH inc | 1.6206 | 1.7616 | 1.6387 | 1.6461 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x2.00) | 1960 (cost x1.77) | 1970 (cost x1.77) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6763 | 0.6701 | 0.7098 | 0.7292 | 0.7302 |
| hours|E | 0.4134 | 0.4381 | 0.4612 | 0.4800 | 0.4856 |
| U rate | 0.0397 | 0.0411 | 0.0424 | 0.0428 | 0.0432 |
| quit/m exp | 0.0333 | 0.0319 | 0.0262 | 0.0233 | 0.0225 |
| quit/m rec | 0.0253 | 0.0233 | 0.0189 | 0.0169 | 0.0161 |
| E->nonE/m exp | 0.0531 | 0.0517 | 0.0460 | 0.0430 | 0.0422 |
| E->nonE/m rec | 0.0457 | 0.0436 | 0.0393 | 0.0372 | 0.0364 |
| dE/pop rec-exp (pts) | -1.6658 | -1.5234 | -1.5415 | -1.5411 | -1.3906 |
| wife share exp | 0.3273 | 0.3592 | 0.4061 | 0.4377 | 0.4464 |
| wife share rec | 0.3279 | 0.3604 | 0.4067 | 0.4395 | 0.4496 |
| wage gap (hourly ratio) | 0.7118 | 0.7887 | 0.8647 | 0.9195 | 0.9447 |
| share Lifecycle | 0.3056 | 0.3744 | 0.3981 | 0.3869 | 0.4023 |
| share PT | 0.2727 | 0.2135 | 0.1919 | 0.1677 | 0.1623 |
| share Career | 0.1748 | 0.1756 | 0.2288 | 0.2871 | 0.2858 |
| share NiLF | 0.2469 | 0.2365 | 0.1812 | 0.1583 | 0.1496 |
| HH income rec/exp - 1 (%) | -16.3030 | -16.0753 | -16.0970 | -15.8408 | -15.5466 |
| cons drop at H job loss exp (%) | -7.0656 | -6.7728 | -6.3361 | -6.1340 | -5.9998 |
| cons drop at H job loss rec (%) | -10.0996 | -9.7860 | -9.3496 | -9.1270 | -8.9816 |
| mean assets/monthly HH inc | 1.6206 | 1.6419 | 1.7022 | 1.7679 | 1.8073 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6763 | 0.6738 | 0.6790 | 0.6835 | 0.6768 |
| hours|E | 0.4134 | 0.4111 | 0.4127 | 0.4182 | 0.4134 |
| U rate | 0.0397 | 0.0399 | 0.0384 | 0.0408 | 0.0395 |
| quit/m exp | 0.0333 | 0.0340 | 0.0339 | 0.0326 | 0.0333 |
| quit/m rec | 0.0253 | 0.0272 | 0.0296 | 0.0236 | 0.0253 |
| E->nonE/m exp | 0.0531 | 0.0538 | 0.0537 | 0.0524 | 0.0531 |
| E->nonE/m rec | 0.0457 | 0.0475 | 0.0498 | 0.0439 | 0.0449 |
| dE/pop rec-exp (pts) | -1.6658 | -1.6229 | 0.2519 | -0.5676 | -1.4773 |
| wife share exp | 0.3273 | 0.3244 | 0.3274 | 0.3291 | 0.3275 |
| wife share rec | 0.3279 | 0.3109 | 0.3311 | 0.3471 | 0.3285 |
| wage gap (hourly ratio) | 0.7118 | 0.7121 | 0.7111 | 0.7137 | 0.7118 |
| share Lifecycle | 0.3056 | 0.3040 | 0.3152 | 0.3125 | 0.3060 |
| share PT | 0.2727 | 0.2729 | 0.2623 | 0.2637 | 0.2719 |
| share Career | 0.1748 | 0.1708 | 0.1727 | 0.1858 | 0.1754 |
| share NiLF | 0.2469 | 0.2523 | 0.2498 | 0.2379 | 0.2467 |
| HH income rec/exp - 1 (%) | -16.3030 | -13.4819 | -15.8810 | -2.3043 | -16.2422 |
| cons drop at H job loss exp (%) | -7.0656 | -7.3982 | -7.0750 | -7.3056 | -7.0665 |
| cons drop at H job loss rec (%) | -10.0996 | -8.3932 | -10.0179 | -9.2465 | -10.0869 |
| mean assets/monthly HH inc | 1.6206 | 1.5610 | 1.6205 | 1.6397 | 1.6216 |


Elapsed 2735s. Calibration file: output/final_calib_v7c_full.json
