# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6644 |
| hours|E | 0.4000 | 0.4141 |
| U rate | nan | 0.0459 |
| quit/m exp | 0.0340 | 0.0311 |
| quit/m rec | 0.0280 | 0.0237 |
| E->nonE/m exp | 0.0500 | 0.0485 |
| E->nonE/m rec | 0.0480 | 0.0431 |
| dE/pop rec-exp (pts) | -1.7000 | -1.5984 |
| wife share exp | nan | 0.3370 |
| wife share rec | nan | 0.3441 |
| wage gap (hourly ratio) | 0.7100 | 0.7477 |
| share Lifecycle | 0.3100 | 0.3038 |
| share PT | 0.2800 | 0.2531 |
| share Career | 0.1900 | 0.1833 |
| share NiLF | 0.2200 | 0.2598 |
| HH income rec/exp - 1 (%) | nan | -14.9712 |
| cons drop at H job loss exp (%) | nan | -5.0354 |
| cons drop at H job loss rec (%) | nan | -7.1846 |
| mean assets/monthly HH inc | nan | 1.4722 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.38 | comp. wage gap x1.08 | cost x0.10 |
|---|---|---|---|---|
| E/pop | 0.6644 | 0.7281 | 0.7336 | 0.7263 |
| hours|E | 0.4141 | 0.4437 | 0.4413 | 0.4194 |
| U rate | 0.0459 | 0.0473 | 0.0461 | 0.0460 |
| quit/m exp | 0.0311 | 0.0215 | 0.0204 | 0.0227 |
| quit/m rec | 0.0237 | 0.0163 | 0.0151 | 0.0165 |
| E->nonE/m exp | 0.0485 | 0.0388 | 0.0376 | 0.0399 |
| E->nonE/m rec | 0.0431 | 0.0357 | 0.0345 | 0.0359 |
| dE/pop rec-exp (pts) | -1.5984 | -1.5558 | -1.3690 | -1.9532 |
| wife share exp | 0.3370 | 0.4092 | 0.4021 | 0.3674 |
| wife share rec | 0.3441 | 0.4167 | 0.4091 | 0.3737 |
| wage gap (hourly ratio) | 0.7477 | 0.8730 | 0.8528 | 0.7678 |
| share Lifecycle | 0.3038 | 0.2935 | 0.3244 | 0.1544 |
| share PT | 0.2531 | 0.1812 | 0.2123 | 0.3415 |
| share Career | 0.1833 | 0.3375 | 0.3025 | 0.2927 |
| share NiLF | 0.2598 | 0.1877 | 0.1608 | 0.2115 |
| HH income rec/exp - 1 (%) | -14.9712 | -14.7435 | -14.8837 | -15.2162 |
| cons drop at H job loss exp (%) | -5.0354 | -4.4100 | -4.4597 | -4.8908 |
| cons drop at H job loss rec (%) | -7.1846 | -6.6915 | -6.6684 | -7.1765 |
| mean assets/monthly HH inc | 1.4722 | 1.7243 | 1.6511 | 1.4745 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.77) | 1960 (cost x1.77) | 1970 (cost x1.77) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6644 | 0.6703 | 0.7121 | 0.7316 | 0.7324 |
| hours|E | 0.4141 | 0.4355 | 0.4515 | 0.4601 | 0.4645 |
| U rate | 0.0459 | 0.0470 | 0.0481 | 0.0485 | 0.0488 |
| quit/m exp | 0.0311 | 0.0278 | 0.0217 | 0.0191 | 0.0183 |
| quit/m rec | 0.0237 | 0.0206 | 0.0157 | 0.0138 | 0.0132 |
| E->nonE/m exp | 0.0485 | 0.0451 | 0.0390 | 0.0364 | 0.0356 |
| E->nonE/m rec | 0.0431 | 0.0400 | 0.0350 | 0.0332 | 0.0325 |
| dE/pop rec-exp (pts) | -1.5984 | -1.2115 | -1.0827 | -1.2304 | -1.0515 |
| wife share exp | 0.3370 | 0.3692 | 0.4117 | 0.4369 | 0.4447 |
| wife share rec | 0.3441 | 0.3763 | 0.4202 | 0.4456 | 0.4538 |
| wage gap (hourly ratio) | 0.7477 | 0.8231 | 0.8991 | 0.9523 | 0.9767 |
| share Lifecycle | 0.3038 | 0.3771 | 0.3881 | 0.3765 | 0.3908 |
| share PT | 0.2531 | 0.1973 | 0.1719 | 0.1454 | 0.1388 |
| share Career | 0.1833 | 0.1950 | 0.2675 | 0.3287 | 0.3283 |
| share NiLF | 0.2598 | 0.2306 | 0.1725 | 0.1494 | 0.1421 |
| HH income rec/exp - 1 (%) | -14.9712 | -14.8559 | -14.4955 | -14.4351 | -14.3072 |
| cons drop at H job loss exp (%) | -5.0354 | -4.6720 | -4.2855 | -4.0319 | -3.9346 |
| cons drop at H job loss rec (%) | -7.1846 | -6.8833 | -6.4779 | -6.2469 | -6.1531 |
| mean assets/monthly HH inc | 1.4722 | 1.6379 | 1.7272 | 1.8225 | 1.8305 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6644 | 0.6608 | 0.6657 | 0.6727 | 0.6675 |
| hours|E | 0.4141 | 0.4119 | 0.4137 | 0.4146 | 0.4140 |
| U rate | 0.0459 | 0.0460 | 0.0449 | 0.0474 | 0.0445 |
| quit/m exp | 0.0311 | 0.0315 | 0.0318 | 0.0304 | 0.0310 |
| quit/m rec | 0.0237 | 0.0254 | 0.0271 | 0.0227 | 0.0234 |
| E->nonE/m exp | 0.0485 | 0.0489 | 0.0491 | 0.0477 | 0.0482 |
| E->nonE/m rec | 0.0431 | 0.0449 | 0.0465 | 0.0420 | 0.0404 |
| dE/pop rec-exp (pts) | -1.5984 | -1.9977 | -0.1888 | -0.3090 | -0.9109 |
| wife share exp | 0.3370 | 0.3350 | 0.3369 | 0.3382 | 0.3377 |
| wife share rec | 0.3441 | 0.3275 | 0.3471 | 0.3560 | 0.3468 |
| wage gap (hourly ratio) | 0.7477 | 0.7473 | 0.7478 | 0.7464 | 0.7479 |
| share Lifecycle | 0.3038 | 0.3006 | 0.3129 | 0.2979 | 0.3063 |
| share PT | 0.2531 | 0.2602 | 0.2387 | 0.2606 | 0.2496 |
| share Career | 0.1833 | 0.1750 | 0.1829 | 0.1894 | 0.1860 |
| share NiLF | 0.2598 | 0.2642 | 0.2654 | 0.2521 | 0.2581 |
| HH income rec/exp - 1 (%) | -14.9712 | -12.8653 | -14.5453 | -1.7687 | -14.6979 |
| cons drop at H job loss exp (%) | -5.0354 | -5.1934 | -5.0527 | -5.2900 | -5.0453 |
| cons drop at H job loss rec (%) | -7.1846 | -6.3677 | -7.1394 | -6.2450 | -7.1268 |
| mean assets/monthly HH inc | 1.4722 | 1.4406 | 1.4754 | 1.4780 | 1.4719 |


Elapsed 3903s. Calibration file: output/final_calib_full.json
