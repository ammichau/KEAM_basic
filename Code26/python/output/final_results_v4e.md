# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6831 |
| hours|E | 0.4000 | 0.4116 |
| U rate | nan | 0.0456 |
| quit/m exp | 0.0340 | 0.0354 |
| quit/m rec | 0.0280 | 0.0235 |
| E->nonE/m exp | 0.0500 | 0.0533 |
| E->nonE/m rec | 0.0480 | 0.0449 |
| dE/pop rec-exp (pts) | -1.7000 | -1.7584 |
| wife share exp | nan | 0.3230 |
| wife share rec | nan | 0.3330 |
| wage gap (hourly ratio) | 0.7100 | 0.7299 |
| share Lifecycle | 0.3100 | 0.3025 |
| share PT | 0.2800 | 0.2550 |
| share Career | 0.1900 | 0.1873 |
| share NiLF | 0.2200 | 0.2552 |
| HH income rec/exp - 1 (%) | nan | -15.2417 |
| cons drop at H job loss exp (%) | nan | -6.7332 |
| cons drop at H job loss rec (%) | nan | -10.0754 |
| mean assets/monthly HH inc | nan | 1.8603 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.26 | comp. wage gap x1.06 | cost x0.52 |
|---|---|---|---|---|
| E/pop | 0.6831 | 0.7316 | 0.7290 | 0.7298 |
| hours|E | 0.4116 | 0.4324 | 0.4302 | 0.4164 |
| U rate | 0.0456 | 0.0474 | 0.0472 | 0.0462 |
| quit/m exp | 0.0354 | 0.0266 | 0.0267 | 0.0286 |
| quit/m rec | 0.0235 | 0.0174 | 0.0170 | 0.0185 |
| E->nonE/m exp | 0.0533 | 0.0445 | 0.0446 | 0.0465 |
| E->nonE/m rec | 0.0449 | 0.0389 | 0.0385 | 0.0401 |
| dE/pop rec-exp (pts) | -1.7584 | -1.6706 | -1.5118 | -2.1811 |
| wife share exp | 0.3230 | 0.3735 | 0.3649 | 0.3453 |
| wife share rec | 0.3330 | 0.3860 | 0.3764 | 0.3549 |
| wage gap (hourly ratio) | 0.7299 | 0.8156 | 0.7965 | 0.7446 |
| share Lifecycle | 0.3025 | 0.2958 | 0.3090 | 0.2362 |
| share PT | 0.2550 | 0.2196 | 0.2477 | 0.2842 |
| share Career | 0.1873 | 0.2954 | 0.2610 | 0.2692 |
| share NiLF | 0.2552 | 0.1892 | 0.1823 | 0.2104 |
| HH income rec/exp - 1 (%) | -15.2417 | -14.7595 | -14.9728 | -15.3476 |
| cons drop at H job loss exp (%) | -6.7332 | -6.0931 | -6.1640 | -6.5250 |
| cons drop at H job loss rec (%) | -10.0754 | -9.5619 | -9.5565 | -9.9191 |
| mean assets/monthly HH inc | 1.8603 | 2.0015 | 1.9061 | 1.8925 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.89) | 1960 (cost x1.94) | 1970 (cost x1.94) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.6831 | 0.6706 | 0.7099 | 0.7295 | 0.7378 |
| hours|E | 0.4116 | 0.4276 | 0.4462 | 0.4554 | 0.4604 |
| U rate | 0.0456 | 0.0472 | 0.0484 | 0.0489 | 0.0493 |
| quit/m exp | 0.0354 | 0.0334 | 0.0255 | 0.0225 | 0.0210 |
| quit/m rec | 0.0235 | 0.0212 | 0.0159 | 0.0141 | 0.0130 |
| E->nonE/m exp | 0.0533 | 0.0512 | 0.0434 | 0.0404 | 0.0389 |
| E->nonE/m rec | 0.0449 | 0.0427 | 0.0373 | 0.0356 | 0.0345 |
| dE/pop rec-exp (pts) | -1.7584 | -1.3050 | -1.2762 | -1.4206 | -1.3686 |
| wife share exp | 0.3230 | 0.3444 | 0.3860 | 0.4116 | 0.4233 |
| wife share rec | 0.3330 | 0.3563 | 0.3989 | 0.4246 | 0.4365 |
| wage gap (hourly ratio) | 0.7299 | 0.7986 | 0.8724 | 0.9253 | 0.9514 |
| share Lifecycle | 0.3025 | 0.3633 | 0.3942 | 0.3896 | 0.4006 |
| share PT | 0.2550 | 0.2006 | 0.1754 | 0.1544 | 0.1465 |
| share Career | 0.1873 | 0.1875 | 0.2487 | 0.3008 | 0.3154 |
| share NiLF | 0.2552 | 0.2485 | 0.1817 | 0.1552 | 0.1375 |
| HH income rec/exp - 1 (%) | -15.2417 | -14.8094 | -14.5648 | -14.4154 | -14.3034 |
| cons drop at H job loss exp (%) | -6.7332 | -6.3554 | -5.8543 | -5.5570 | -5.4238 |
| cons drop at H job loss rec (%) | -10.0754 | -9.7668 | -9.2468 | -8.9035 | -8.7645 |
| mean assets/monthly HH inc | 1.8603 | 1.9565 | 2.0092 | 2.1613 | 2.1723 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6831 | 0.6724 | 0.6847 | 0.6899 | 0.6880 |
| hours|E | 0.4116 | 0.4071 | 0.4108 | 0.4117 | 0.4114 |
| U rate | 0.0456 | 0.0458 | 0.0439 | 0.0466 | 0.0436 |
| quit/m exp | 0.0354 | 0.0372 | 0.0363 | 0.0350 | 0.0351 |
| quit/m rec | 0.0235 | 0.0271 | 0.0304 | 0.0233 | 0.0227 |
| E->nonE/m exp | 0.0533 | 0.0551 | 0.0542 | 0.0529 | 0.0529 |
| E->nonE/m rec | 0.0449 | 0.0487 | 0.0520 | 0.0447 | 0.0402 |
| dE/pop rec-exp (pts) | -1.7584 | -2.9074 | 0.0732 | -0.2042 | -0.5891 |
| wife share exp | 0.3230 | 0.3176 | 0.3227 | 0.3229 | 0.3238 |
| wife share rec | 0.3330 | 0.3096 | 0.3362 | 0.3440 | 0.3367 |
| wage gap (hourly ratio) | 0.7299 | 0.7304 | 0.7290 | 0.7303 | 0.7298 |
| share Lifecycle | 0.3025 | 0.2973 | 0.3152 | 0.2994 | 0.3042 |
| share PT | 0.2550 | 0.2544 | 0.2369 | 0.2606 | 0.2531 |
| share Career | 0.1873 | 0.1773 | 0.1871 | 0.1898 | 0.1908 |
| share NiLF | 0.2552 | 0.2710 | 0.2608 | 0.2502 | 0.2519 |
| HH income rec/exp - 1 (%) | -15.2417 | -12.9289 | -14.7499 | -2.0327 | -14.8386 |
| cons drop at H job loss exp (%) | -6.7332 | -7.0318 | -6.7472 | -6.9466 | -6.7439 |
| cons drop at H job loss rec (%) | -10.0754 | -8.2444 | -9.9382 | -9.0423 | -9.9923 |
| mean assets/monthly HH inc | 1.8603 | 1.8217 | 1.8592 | 1.8864 | 1.8566 |


Elapsed 1986s. Calibration file: output/final_calib_v4e_full.json
