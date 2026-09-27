# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6475 |
| hours|E | 0.4000 | 0.4227 |
| U rate | nan | 0.0416 |
| quit/m exp | 0.0340 | 0.0338 |
| quit/m rec | 0.0280 | 0.0254 |
| E->nonE/m exp | 0.0500 | 0.0535 |
| E->nonE/m rec | 0.0480 | 0.0439 |
| dE/pop rec-exp (pts) | -1.7000 | -1.7174 |
| wife share exp | nan | 0.3245 |
| wife share rec | nan | 0.3266 |
| wage gap (hourly ratio) | 0.7100 | 0.7354 |
| share Lifecycle | 0.3100 | 0.2898 |
| share PT | 0.2800 | 0.2335 |
| share Career | 0.1900 | 0.1862 |
| share NiLF | 0.2200 | 0.2904 |
| HH income rec/exp - 1 (%) | nan | -16.1871 |
| cons drop at H job loss exp (%) | nan | -6.5250 |
| cons drop at H job loss rec (%) | nan | -9.3117 |
| mean assets/monthly HH inc | nan | 1.8103 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.59 | comp. wage gap x1.15 | cost x0.16 |
|---|---|---|---|---|
| E/pop | 0.6475 | 0.7299 | 0.7302 | 0.7316 |
| hours|E | 0.4227 | 0.4792 | 0.4678 | 0.4226 |
| U rate | 0.0416 | 0.0421 | 0.0427 | 0.0413 |
| quit/m exp | 0.0338 | 0.0210 | 0.0201 | 0.0280 |
| quit/m rec | 0.0254 | 0.0155 | 0.0145 | 0.0212 |
| E->nonE/m exp | 0.0535 | 0.0407 | 0.0398 | 0.0477 |
| E->nonE/m rec | 0.0439 | 0.0340 | 0.0331 | 0.0397 |
| dE/pop rec-exp (pts) | -1.7174 | -1.5260 | -1.3560 | -2.3831 |
| wife share exp | 0.3245 | 0.4400 | 0.4262 | 0.3605 |
| wife share rec | 0.3266 | 0.4432 | 0.4320 | 0.3614 |
| wage gap (hourly ratio) | 0.7354 | 0.9381 | 0.9252 | 0.7516 |
| share Lifecycle | 0.2898 | 0.2815 | 0.3335 | 0.1581 |
| share PT | 0.2335 | 0.1590 | 0.1783 | 0.2990 |
| share Career | 0.1862 | 0.3688 | 0.3200 | 0.3156 |
| share NiLF | 0.2904 | 0.1908 | 0.1681 | 0.2273 |
| HH income rec/exp - 1 (%) | -16.1871 | -15.7947 | -15.3455 | -16.4760 |
| cons drop at H job loss exp (%) | -6.5250 | -5.5156 | -5.6423 | -6.2415 |
| cons drop at H job loss rec (%) | -9.3117 | -8.3653 | -8.3774 | -9.0075 |
| mean assets/monthly HH inc | 1.8103 | 1.9845 | 1.8682 | 1.8349 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.21) | 1960 (cost x1.10) | 1970 (cost x1.13) | 1980 (cost x1.38) |
|---|---|---|---|---|---|
| E/pop | 0.6475 | 0.6715 | 0.7100 | 0.7297 | 0.7204 |
| hours|E | 0.4227 | 0.4502 | 0.4666 | 0.4807 | 0.4870 |
| U rate | 0.0416 | 0.0427 | 0.0429 | 0.0426 | 0.0432 |
| quit/m exp | 0.0338 | 0.0278 | 0.0226 | 0.0195 | 0.0191 |
| quit/m rec | 0.0254 | 0.0207 | 0.0164 | 0.0144 | 0.0142 |
| E->nonE/m exp | 0.0535 | 0.0475 | 0.0423 | 0.0392 | 0.0387 |
| E->nonE/m rec | 0.0439 | 0.0392 | 0.0350 | 0.0330 | 0.0328 |
| dE/pop rec-exp (pts) | -1.7174 | -1.3267 | -1.3252 | -1.3982 | -1.0683 |
| wife share exp | 0.3245 | 0.3703 | 0.4129 | 0.4434 | 0.4476 |
| wife share rec | 0.3266 | 0.3712 | 0.4191 | 0.4465 | 0.4512 |
| wage gap (hourly ratio) | 0.7354 | 0.8181 | 0.8984 | 0.9562 | 0.9776 |
| share Lifecycle | 0.2898 | 0.3269 | 0.3354 | 0.3150 | 0.3469 |
| share PT | 0.2335 | 0.2013 | 0.1694 | 0.1571 | 0.1498 |
| share Career | 0.1862 | 0.2260 | 0.2983 | 0.3517 | 0.3275 |
| share NiLF | 0.2904 | 0.2458 | 0.1969 | 0.1762 | 0.1758 |
| HH income rec/exp - 1 (%) | -16.1871 | -16.2207 | -15.2789 | -15.7974 | -15.6494 |
| cons drop at H job loss exp (%) | -6.5250 | -6.2143 | -5.7726 | -5.4389 | -5.3794 |
| cons drop at H job loss rec (%) | -9.3117 | -8.8954 | -8.5430 | -8.2774 | -8.1771 |
| mean assets/monthly HH inc | 1.8103 | 1.8560 | 1.9127 | 1.9667 | 1.9728 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6475 | 0.6455 | 0.6518 | 0.6555 | 0.6468 |
| hours|E | 0.4227 | 0.4205 | 0.4212 | 0.4267 | 0.4227 |
| U rate | 0.0416 | 0.0419 | 0.0401 | 0.0425 | 0.0421 |
| quit/m exp | 0.0338 | 0.0339 | 0.0343 | 0.0326 | 0.0338 |
| quit/m rec | 0.0254 | 0.0260 | 0.0308 | 0.0247 | 0.0254 |
| E->nonE/m exp | 0.0535 | 0.0536 | 0.0540 | 0.0523 | 0.0536 |
| E->nonE/m rec | 0.0439 | 0.0445 | 0.0493 | 0.0432 | 0.0448 |
| dE/pop rec-exp (pts) | -1.7174 | -1.7433 | 0.5954 | -0.4240 | -1.9306 |
| wife share exp | 0.3245 | 0.3217 | 0.3247 | 0.3264 | 0.3243 |
| wife share rec | 0.3266 | 0.3100 | 0.3299 | 0.3446 | 0.3258 |
| wage gap (hourly ratio) | 0.7354 | 0.7355 | 0.7340 | 0.7367 | 0.7354 |
| share Lifecycle | 0.2898 | 0.2840 | 0.2965 | 0.2938 | 0.2900 |
| share PT | 0.2335 | 0.2400 | 0.2242 | 0.2310 | 0.2329 |
| share Career | 0.1862 | 0.1842 | 0.1865 | 0.1958 | 0.1860 |
| share NiLF | 0.2904 | 0.2919 | 0.2929 | 0.2794 | 0.2910 |
| HH income rec/exp - 1 (%) | -16.1871 | -13.2625 | -15.7957 | -2.3360 | -16.2716 |
| cons drop at H job loss exp (%) | -6.5250 | -6.8518 | -6.5381 | -6.7589 | -6.5205 |
| cons drop at H job loss rec (%) | -9.3117 | -7.9115 | -9.1874 | -8.6915 | -9.3084 |
| mean assets/monthly HH inc | 1.8103 | 1.7647 | 1.8104 | 1.8459 | 1.8102 |


Elapsed 2392s. Calibration file: output/final_calib_v9_full.json
