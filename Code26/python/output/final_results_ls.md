# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6641 |
| hours|E | 0.4000 | 0.4097 |
| U rate | nan | 0.0447 |
| quit/m exp | 0.0340 | 0.0332 |
| quit/m rec | 0.0280 | 0.0250 |
| E->nonE/m exp | 0.0500 | 0.0505 |
| E->nonE/m rec | 0.0480 | 0.0447 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6887 |
| wife share exp | nan | 0.3319 |
| wife share rec | nan | 0.3383 |
| wage gap (hourly ratio) | 0.7100 | 0.7382 |
| share Lifecycle | 0.3100 | 0.3085 |
| share PT | 0.2800 | 0.2502 |
| share Career | 0.1900 | 0.1719 |
| share NiLF | 0.2200 | 0.2694 |
| HH income rec/exp - 1 (%) | nan | -15.0947 |
| cons drop at H job loss exp (%) | nan | -5.0741 |
| cons drop at H job loss rec (%) | nan | -7.2569 |
| mean assets/monthly HH inc | nan | 1.4074 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.40 | comp. wage gap x1.08 | cost x0.10 |
|---|---|---|---|---|
| E/pop | 0.6641 | 0.7307 | 0.7282 | 0.7210 |
| hours|E | 0.4097 | 0.4401 | 0.4354 | 0.4162 |
| U rate | 0.0447 | 0.0464 | 0.0454 | 0.0450 |
| quit/m exp | 0.0332 | 0.0225 | 0.0225 | 0.0247 |
| quit/m rec | 0.0250 | 0.0167 | 0.0165 | 0.0179 |
| E->nonE/m exp | 0.0505 | 0.0397 | 0.0398 | 0.0419 |
| E->nonE/m rec | 0.0447 | 0.0365 | 0.0363 | 0.0377 |
| dE/pop rec-exp (pts) | -1.6887 | -1.6796 | -1.4169 | -1.9712 |
| wife share exp | 0.3319 | 0.4074 | 0.3926 | 0.3617 |
| wife share rec | 0.3383 | 0.4158 | 0.3999 | 0.3678 |
| wage gap (hourly ratio) | 0.7382 | 0.8699 | 0.8362 | 0.7595 |
| share Lifecycle | 0.3085 | 0.2956 | 0.3281 | 0.1669 |
| share PT | 0.2502 | 0.1808 | 0.2229 | 0.3269 |
| share Career | 0.1719 | 0.3333 | 0.2727 | 0.2821 |
| share NiLF | 0.2694 | 0.1902 | 0.1762 | 0.2242 |
| HH income rec/exp - 1 (%) | -15.0947 | -14.6345 | -14.8822 | -15.2151 |
| cons drop at H job loss exp (%) | -5.0741 | -4.4215 | -4.5507 | -4.9398 |
| cons drop at H job loss rec (%) | -7.2569 | -6.7607 | -6.7850 | -7.2041 |
| mean assets/monthly HH inc | 1.4074 | 1.7099 | 1.5547 | 1.4370 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.77) | 1960 (cost x1.83) | 1970 (cost x1.83) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6641 | 0.6720 | 0.7104 | 0.7306 | 0.7346 |
| hours|E | 0.4097 | 0.4313 | 0.4475 | 0.4560 | 0.4603 |
| U rate | 0.0447 | 0.0460 | 0.0470 | 0.0475 | 0.0479 |
| quit/m exp | 0.0332 | 0.0296 | 0.0232 | 0.0204 | 0.0195 |
| quit/m rec | 0.0250 | 0.0217 | 0.0167 | 0.0147 | 0.0140 |
| E->nonE/m exp | 0.0505 | 0.0469 | 0.0405 | 0.0377 | 0.0368 |
| E->nonE/m rec | 0.0447 | 0.0415 | 0.0364 | 0.0344 | 0.0337 |
| dE/pop rec-exp (pts) | -1.6887 | -1.3702 | -1.0582 | -1.2760 | -1.1518 |
| wife share exp | 0.3319 | 0.3650 | 0.4055 | 0.4314 | 0.4404 |
| wife share rec | 0.3383 | 0.3722 | 0.4143 | 0.4398 | 0.4493 |
| wage gap (hourly ratio) | 0.7382 | 0.8124 | 0.8867 | 0.9390 | 0.9638 |
| share Lifecycle | 0.3085 | 0.3777 | 0.3958 | 0.3823 | 0.3933 |
| share PT | 0.2502 | 0.1992 | 0.1723 | 0.1542 | 0.1388 |
| share Career | 0.1719 | 0.1865 | 0.2515 | 0.3090 | 0.3223 |
| share NiLF | 0.2694 | 0.2367 | 0.1804 | 0.1546 | 0.1456 |
| HH income rec/exp - 1 (%) | -15.0947 | -14.8363 | -14.4482 | -14.4858 | -14.3486 |
| cons drop at H job loss exp (%) | -5.0741 | -4.7088 | -4.3527 | -4.0974 | -3.9956 |
| cons drop at H job loss rec (%) | -7.2569 | -6.9451 | -6.5670 | -6.3352 | -6.2188 |
| mean assets/monthly HH inc | 1.4074 | 1.6443 | 1.7116 | 1.8297 | 1.8240 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6641 | 0.6607 | 0.6645 | 0.6718 | 0.6673 |
| hours|E | 0.4097 | 0.4072 | 0.4090 | 0.4105 | 0.4095 |
| U rate | 0.0447 | 0.0451 | 0.0436 | 0.0463 | 0.0434 |
| quit/m exp | 0.0332 | 0.0335 | 0.0340 | 0.0324 | 0.0330 |
| quit/m rec | 0.0250 | 0.0268 | 0.0290 | 0.0241 | 0.0245 |
| E->nonE/m exp | 0.0505 | 0.0508 | 0.0513 | 0.0497 | 0.0502 |
| E->nonE/m rec | 0.0447 | 0.0466 | 0.0487 | 0.0439 | 0.0414 |
| dE/pop rec-exp (pts) | -1.6887 | -2.2276 | -0.2177 | -0.3088 | -0.8384 |
| wife share exp | 0.3319 | 0.3300 | 0.3315 | 0.3330 | 0.3326 |
| wife share rec | 0.3383 | 0.3215 | 0.3408 | 0.3507 | 0.3412 |
| wage gap (hourly ratio) | 0.7382 | 0.7376 | 0.7382 | 0.7367 | 0.7384 |
| share Lifecycle | 0.3085 | 0.3038 | 0.3185 | 0.3056 | 0.3094 |
| share PT | 0.2502 | 0.2606 | 0.2331 | 0.2544 | 0.2471 |
| share Career | 0.1719 | 0.1642 | 0.1723 | 0.1785 | 0.1760 |
| share NiLF | 0.2694 | 0.2715 | 0.2760 | 0.2615 | 0.2675 |
| HH income rec/exp - 1 (%) | -15.0947 | -13.0107 | -14.6772 | -1.8125 | -14.7827 |
| cons drop at H job loss exp (%) | -5.0741 | -5.2423 | -5.0930 | -5.3398 | -5.0862 |
| cons drop at H job loss rec (%) | -7.2569 | -6.4003 | -7.1684 | -6.2958 | -7.1724 |
| mean assets/monthly HH inc | 1.4074 | 1.3615 | 1.4020 | 1.4017 | 1.4056 |


Elapsed 3052s. Calibration file: output/final_calib_ls_full.json
