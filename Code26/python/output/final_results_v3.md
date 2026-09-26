# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6668 |
| hours|E | 0.4000 | 0.4085 |
| U rate | nan | 0.0447 |
| quit/m exp | 0.0340 | 0.0329 |
| quit/m rec | 0.0280 | 0.0248 |
| E->nonE/m exp | 0.0500 | 0.0499 |
| E->nonE/m rec | 0.0480 | 0.0470 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6981 |
| wife share exp | nan | 0.3316 |
| wife share rec | nan | 0.3471 |
| wage gap (hourly ratio) | 0.7100 | 0.7363 |
| share Lifecycle | 0.3100 | 0.2969 |
| share PT | 0.2800 | 0.2587 |
| share Career | 0.1900 | 0.1756 |
| share NiLF | 0.2200 | 0.2687 |
| HH income rec/exp - 1 (%) | nan | -16.5952 |
| cons drop at H job loss exp (%) | nan | -5.0814 |
| cons drop at H job loss rec (%) | nan | -8.5421 |
| mean assets/monthly HH inc | nan | 1.4182 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.38 | comp. wage gap x1.08 | cost x0.10 |
|---|---|---|---|---|
| E/pop | 0.6668 | 0.7296 | 0.7292 | 0.7163 |
| hours|E | 0.4085 | 0.4375 | 0.4345 | 0.4143 |
| U rate | 0.0447 | 0.0462 | 0.0454 | 0.0452 |
| quit/m exp | 0.0329 | 0.0229 | 0.0226 | 0.0254 |
| quit/m rec | 0.0248 | 0.0171 | 0.0165 | 0.0186 |
| E->nonE/m exp | 0.0499 | 0.0399 | 0.0395 | 0.0423 |
| E->nonE/m rec | 0.0470 | 0.0393 | 0.0387 | 0.0409 |
| dE/pop rec-exp (pts) | -1.6981 | -1.8620 | -1.6917 | -2.1370 |
| wife share exp | 0.3316 | 0.4029 | 0.3916 | 0.3584 |
| wife share rec | 0.3471 | 0.4210 | 0.4091 | 0.3721 |
| wage gap (hourly ratio) | 0.7363 | 0.8603 | 0.8345 | 0.7560 |
| share Lifecycle | 0.2969 | 0.2819 | 0.3177 | 0.1740 |
| share PT | 0.2587 | 0.1940 | 0.2298 | 0.3185 |
| share Career | 0.1756 | 0.3283 | 0.2733 | 0.2767 |
| share NiLF | 0.2687 | 0.1958 | 0.1792 | 0.2308 |
| HH income rec/exp - 1 (%) | -16.5952 | -15.9574 | -16.1059 | -16.9089 |
| cons drop at H job loss exp (%) | -5.0814 | -4.4466 | -4.5321 | -4.9783 |
| cons drop at H job loss rec (%) | -8.5421 | -7.9334 | -7.9552 | -8.4663 |
| mean assets/monthly HH inc | 1.4182 | 1.7372 | 1.5848 | 1.4163 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.94) | 1960 (cost x1.94) | 1970 (cost x1.94) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6668 | 0.6698 | 0.7112 | 0.7318 | 0.7439 |
| hours|E | 0.4085 | 0.4302 | 0.4469 | 0.4550 | 0.4596 |
| U rate | 0.0447 | 0.0460 | 0.0468 | 0.0475 | 0.0480 |
| quit/m exp | 0.0329 | 0.0303 | 0.0235 | 0.0207 | 0.0189 |
| quit/m rec | 0.0248 | 0.0221 | 0.0171 | 0.0153 | 0.0138 |
| E->nonE/m exp | 0.0499 | 0.0473 | 0.0405 | 0.0376 | 0.0358 |
| E->nonE/m rec | 0.0470 | 0.0443 | 0.0393 | 0.0375 | 0.0361 |
| dE/pop rec-exp (pts) | -1.6981 | -1.2064 | -1.1779 | -1.4832 | -1.4005 |
| wife share exp | 0.3316 | 0.3625 | 0.4048 | 0.4306 | 0.4426 |
| wife share rec | 0.3471 | 0.3799 | 0.4228 | 0.4485 | 0.4617 |
| wage gap (hourly ratio) | 0.7363 | 0.8090 | 0.8841 | 0.9362 | 0.9622 |
| share Lifecycle | 0.2969 | 0.3729 | 0.3954 | 0.3758 | 0.3783 |
| share PT | 0.2587 | 0.1988 | 0.1675 | 0.1573 | 0.1485 |
| share Career | 0.1756 | 0.1837 | 0.2548 | 0.3104 | 0.3358 |
| share NiLF | 0.2687 | 0.2446 | 0.1823 | 0.1565 | 0.1373 |
| HH income rec/exp - 1 (%) | -16.5952 | -16.0890 | -15.7975 | -15.7125 | -15.4394 |
| cons drop at H job loss exp (%) | -5.0814 | -4.7363 | -4.3367 | -4.0975 | -3.9665 |
| cons drop at H job loss rec (%) | -8.5421 | -8.1588 | -7.7327 | -7.4613 | -7.3063 |
| mean assets/monthly HH inc | 1.4182 | 1.6555 | 1.7285 | 1.8499 | 1.8484 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6668 | 0.6594 | 0.6680 | 0.6748 | 0.6724 |
| hours|E | 0.4085 | 0.4047 | 0.4081 | 0.4094 | 0.4082 |
| U rate | 0.0447 | 0.0450 | 0.0439 | 0.0465 | 0.0425 |
| quit/m exp | 0.0329 | 0.0339 | 0.0333 | 0.0322 | 0.0325 |
| quit/m rec | 0.0248 | 0.0282 | 0.0273 | 0.0243 | 0.0239 |
| E->nonE/m exp | 0.0499 | 0.0509 | 0.0503 | 0.0492 | 0.0494 |
| E->nonE/m rec | 0.0470 | 0.0503 | 0.0495 | 0.0464 | 0.0405 |
| dE/pop rec-exp (pts) | -1.6981 | -2.5343 | -0.6339 | -0.2448 | -0.0905 |
| wife share exp | 0.3316 | 0.3283 | 0.3313 | 0.3323 | 0.3327 |
| wife share rec | 0.3471 | 0.3183 | 0.3492 | 0.3615 | 0.3526 |
| wage gap (hourly ratio) | 0.7363 | 0.7359 | 0.7359 | 0.7350 | 0.7365 |
| share Lifecycle | 0.2969 | 0.2910 | 0.3029 | 0.2942 | 0.3004 |
| share PT | 0.2587 | 0.2667 | 0.2498 | 0.2596 | 0.2515 |
| share Career | 0.1756 | 0.1635 | 0.1767 | 0.1837 | 0.1821 |
| share NiLF | 0.2687 | 0.2787 | 0.2706 | 0.2625 | 0.2660 |
| HH income rec/exp - 1 (%) | -16.5952 | -13.2397 | -16.2504 | -3.1574 | -15.9817 |
| cons drop at H job loss exp (%) | -5.0814 | -5.2495 | -5.0971 | -5.2748 | -5.0993 |
| cons drop at H job loss rec (%) | -8.5421 | -6.4684 | -8.4797 | -7.3060 | -8.4015 |
| mean assets/monthly HH inc | 1.4182 | 1.3629 | 1.4200 | 1.4563 | 1.4203 |


Elapsed 3014s. Calibration file: output/final_calib_v3_full.json
