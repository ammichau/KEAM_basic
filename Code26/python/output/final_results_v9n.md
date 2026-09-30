# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.6816 |
| hours|E | 0.4000 | 0.4093 |
| U rate | nan | 0.0437 |
| quit/m exp | 0.0340 | 0.0357 |
| quit/m rec | 0.0280 | 0.0270 |
| E->nonE/m exp | 0.0500 | 0.0551 |
| E->nonE/m rec | 0.0480 | 0.0485 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6678 |
| wife share exp | nan | 0.3210 |
| wife share rec | nan | 0.3376 |
| wage gap (hourly ratio) | 0.7100 | 0.7124 |
| share Lifecycle | 0.3100 | 0.2798 |
| share PT | 0.2800 | 0.2815 |
| share Career | 0.1900 | 0.1758 |
| share NiLF | 0.2200 | 0.2629 |
| HH income rec/exp - 1 (%) | nan | -2.6138 |
| cons drop at H job loss exp (%) | nan | -6.7609 |
| cons drop at H job loss rec (%) | nan | -8.7096 |
| mean assets/monthly HH inc | nan | 1.8570 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.35 | comp. wage gap x1.09 | cost x0.27 |
|---|---|---|---|---|
| E/pop | 0.6816 | 0.7293 | 0.7308 | 0.7303 |
| hours|E | 0.4093 | 0.4455 | 0.4397 | 0.4280 |
| U rate | 0.0437 | 0.0440 | 0.0446 | 0.0442 |
| quit/m exp | 0.0357 | 0.0271 | 0.0267 | 0.0278 |
| quit/m rec | 0.0270 | 0.0204 | 0.0198 | 0.0213 |
| E->nonE/m exp | 0.0551 | 0.0465 | 0.0461 | 0.0472 |
| E->nonE/m rec | 0.0485 | 0.0419 | 0.0413 | 0.0429 |
| dE/pop rec-exp (pts) | -1.6678 | -1.7188 | -1.5622 | -1.9300 |
| wife share exp | 0.3210 | 0.3910 | 0.3845 | 0.3546 |
| wife share rec | 0.3376 | 0.4084 | 0.4020 | 0.3706 |
| wage gap (hourly ratio) | 0.7124 | 0.8275 | 0.8233 | 0.7363 |
| share Lifecycle | 0.2798 | 0.2871 | 0.3119 | 0.1925 |
| share PT | 0.2815 | 0.2231 | 0.2506 | 0.2979 |
| share Career | 0.1758 | 0.2850 | 0.2502 | 0.2910 |
| share NiLF | 0.2629 | 0.2048 | 0.1873 | 0.2185 |
| HH income rec/exp - 1 (%) | -2.6138 | -2.1153 | -2.1484 | -2.7347 |
| cons drop at H job loss exp (%) | -6.7609 | -6.2660 | -6.2237 | -6.5686 |
| cons drop at H job loss rec (%) | -8.7096 | -8.1258 | -8.0733 | -8.4234 |
| mean assets/monthly HH inc | 1.8570 | 1.9754 | 1.9104 | 1.8693 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.83) | 1960 (cost x1.66) | 1970 (cost x1.61) | 1980 (cost x1.97) |
|---|---|---|---|---|---|
| E/pop | 0.6816 | 0.6694 | 0.7107 | 0.7315 | 0.7212 |
| hours|E | 0.4093 | 0.4187 | 0.4407 | 0.4570 | 0.4566 |
| U rate | 0.0437 | 0.0443 | 0.0449 | 0.0444 | 0.0447 |
| quit/m exp | 0.0357 | 0.0344 | 0.0279 | 0.0248 | 0.0251 |
| quit/m rec | 0.0270 | 0.0259 | 0.0207 | 0.0188 | 0.0189 |
| E->nonE/m exp | 0.0551 | 0.0538 | 0.0473 | 0.0442 | 0.0445 |
| E->nonE/m rec | 0.0485 | 0.0475 | 0.0422 | 0.0402 | 0.0404 |
| dE/pop rec-exp (pts) | -1.6678 | -1.4918 | -1.4666 | -1.6535 | -1.4916 |
| wife share exp | 0.3210 | 0.3424 | 0.3880 | 0.4190 | 0.4207 |
| wife share rec | 0.3376 | 0.3602 | 0.4064 | 0.4372 | 0.4393 |
| wage gap (hourly ratio) | 0.7124 | 0.7834 | 0.8574 | 0.9107 | 0.9319 |
| share Lifecycle | 0.2798 | 0.3156 | 0.3398 | 0.3387 | 0.3504 |
| share PT | 0.2815 | 0.2581 | 0.2348 | 0.2123 | 0.2127 |
| share Career | 0.1758 | 0.1625 | 0.2229 | 0.2706 | 0.2515 |
| share NiLF | 0.2629 | 0.2637 | 0.2025 | 0.1783 | 0.1854 |
| HH income rec/exp - 1 (%) | -2.6138 | -2.1870 | -1.8312 | -1.6460 | -1.5229 |
| cons drop at H job loss exp (%) | -6.7609 | -6.5268 | -6.1194 | -5.8812 | -5.8373 |
| cons drop at H job loss rec (%) | -8.7096 | -8.4900 | -7.9969 | -7.7503 | -7.7008 |
| mean assets/monthly HH inc | 1.8570 | 1.8904 | 1.9434 | 1.9929 | 2.0003 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.6816 | 0.6791 | 0.6859 | 0.6816 | 0.6834 |
| hours|E | 0.4093 | 0.4069 | 0.4080 | 0.4093 | 0.4092 |
| U rate | 0.0437 | 0.0439 | 0.0416 | 0.0437 | 0.0426 |
| quit/m exp | 0.0357 | 0.0359 | 0.0363 | 0.0357 | 0.0356 |
| quit/m rec | 0.0270 | 0.0280 | 0.0331 | 0.0270 | 0.0269 |
| E->nonE/m exp | 0.0551 | 0.0553 | 0.0557 | 0.0551 | 0.0549 |
| E->nonE/m rec | 0.0485 | 0.0495 | 0.0545 | 0.0485 | 0.0459 |
| dE/pop rec-exp (pts) | -1.6678 | -1.9186 | 0.7838 | -1.6678 | -1.1292 |
| wife share exp | 0.3210 | 0.3181 | 0.3211 | 0.3210 | 0.3214 |
| wife share rec | 0.3376 | 0.3201 | 0.3411 | 0.3376 | 0.3394 |
| wage gap (hourly ratio) | 0.7124 | 0.7114 | 0.7111 | 0.7124 | 0.7126 |
| share Lifecycle | 0.2798 | 0.2765 | 0.2881 | 0.2798 | 0.2798 |
| share PT | 0.2815 | 0.2869 | 0.2692 | 0.2815 | 0.2796 |
| share Career | 0.1758 | 0.1729 | 0.1771 | 0.1758 | 0.1781 |
| share NiLF | 0.2629 | 0.2637 | 0.2656 | 0.2629 | 0.2625 |
| HH income rec/exp - 1 (%) | -2.6138 | 0.5472 | -2.0927 | -2.6138 | -2.3900 |
| cons drop at H job loss exp (%) | -6.7609 | -7.0237 | -6.7764 | -6.7609 | -6.7607 |
| cons drop at H job loss rec (%) | -8.7096 | -7.1369 | -8.6093 | -8.7096 | -8.6706 |
| mean assets/monthly HH inc | 1.8570 | 1.8188 | 1.8583 | 1.8570 | 1.8573 |


Elapsed 3783s. Calibration file: output/final_calib_v9n_full.json
