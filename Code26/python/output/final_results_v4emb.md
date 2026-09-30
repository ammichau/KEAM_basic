# Final model results (100 types)

### Baseline calibration

| moment | target | baseline 1940s |
|---|---|---|
| E/pop | 0.6200 | 0.7172 |
| hours|E | 0.4000 | 0.4035 |
| U rate | nan | 0.0393 |
| quit/m exp | 0.0340 | 0.0267 |
| quit/m rec | 0.0280 | 0.0181 |
| E->nonE/m exp | 0.0500 | 0.0400 |
| E->nonE/m rec | 0.0480 | 0.0333 |
| dE/pop rec-exp (pts) | -1.7000 | -1.6454 |
| wife share exp | nan | 0.3333 |
| wife share rec | nan | 0.3414 |
| wage gap (hourly ratio) | 0.7100 | 0.7429 |
| share Lifecycle | 0.3100 | 0.2963 |
| share PT | 0.2800 | 0.2896 |
| share Career | 0.1900 | 0.2029 |
| share NiLF | 0.2200 | 0.2112 |
| HH income rec/exp - 1 (%) | nan | -15.5393 |
| cons drop at H job loss exp (%) | nan | -6.1281 |
| cons drop at H job loss rec (%) | nan | -9.9832 |
| mean assets/monthly HH inc | nan | 3.6205 |

### Single-factor experiments sized to the 1970s employment rate (0.73)

| moment | baseline | RoE x1.07 | comp. wage gap x1.02 | cost x0.89 |
|---|---|---|---|---|
| E/pop | 0.7172 | 0.7305 | 0.7300 | 0.7299 |
| hours|E | 0.4035 | 0.4101 | 0.4099 | 0.4055 |
| U rate | 0.0393 | 0.0395 | 0.0394 | 0.0392 |
| quit/m exp | 0.0267 | 0.0246 | 0.0246 | 0.0257 |
| quit/m rec | 0.0181 | 0.0166 | 0.0163 | 0.0172 |
| E->nonE/m exp | 0.0400 | 0.0378 | 0.0379 | 0.0389 |
| E->nonE/m rec | 0.0333 | 0.0318 | 0.0315 | 0.0325 |
| dE/pop rec-exp (pts) | -1.6454 | -1.6686 | -1.5789 | -1.7255 |
| wife share exp | 0.3333 | 0.3483 | 0.3467 | 0.3396 |
| wife share rec | 0.3414 | 0.3561 | 0.3548 | 0.3470 |
| wage gap (hourly ratio) | 0.7429 | 0.7668 | 0.7645 | 0.7455 |
| share Lifecycle | 0.2963 | 0.2973 | 0.3017 | 0.2806 |
| share PT | 0.2896 | 0.2719 | 0.2777 | 0.2898 |
| share Career | 0.2029 | 0.2387 | 0.2319 | 0.2310 |
| share NiLF | 0.2112 | 0.1921 | 0.1888 | 0.1985 |
| HH income rec/exp - 1 (%) | -15.5393 | -15.5140 | -15.4990 | -15.6384 |
| cons drop at H job loss exp (%) | -6.1281 | -5.9639 | -6.0151 | -6.1034 |
| cons drop at H job loss rec (%) | -9.9832 | -9.8311 | -9.8372 | -9.9616 |
| mean assets/monthly HH inc | 3.6205 | 3.8335 | 3.6697 | 3.6409 |

### Cohort accounting (tau_w and gamma_e from the data, cost residual)

| moment | 1940 | 1950 (cost x1.89) | 1960 (cost x1.94) | 1970 (cost x1.97) | 1980 (cost x2.00) |
|---|---|---|---|---|---|
| E/pop | 0.7172 | 0.6697 | 0.7089 | 0.7313 | 0.7408 |
| hours|E | 0.4035 | 0.4167 | 0.4340 | 0.4439 | 0.4500 |
| U rate | 0.0393 | 0.0399 | 0.0412 | 0.0416 | 0.0421 |
| quit/m exp | 0.0267 | 0.0268 | 0.0201 | 0.0171 | 0.0156 |
| quit/m rec | 0.0181 | 0.0187 | 0.0138 | 0.0115 | 0.0106 |
| E->nonE/m exp | 0.0400 | 0.0401 | 0.0333 | 0.0303 | 0.0289 |
| E->nonE/m rec | 0.0333 | 0.0338 | 0.0291 | 0.0269 | 0.0260 |
| dE/pop rec-exp (pts) | -1.6454 | -0.6978 | -0.5831 | -0.5471 | -0.5023 |
| wife share exp | 0.3333 | 0.3456 | 0.3860 | 0.4118 | 0.4246 |
| wife share rec | 0.3414 | 0.3548 | 0.3960 | 0.4227 | 0.4347 |
| wage gap (hourly ratio) | 0.7429 | 0.8188 | 0.8926 | 0.9468 | 0.9735 |
| share Lifecycle | 0.2963 | 0.3533 | 0.3965 | 0.4088 | 0.4121 |
| share PT | 0.2896 | 0.2031 | 0.1675 | 0.1360 | 0.1258 |
| share Career | 0.2029 | 0.1940 | 0.2521 | 0.3025 | 0.3252 |
| share NiLF | 0.2112 | 0.2496 | 0.1840 | 0.1527 | 0.1369 |
| HH income rec/exp - 1 (%) | -15.5393 | -15.2264 | -14.9215 | -14.6778 | -14.7216 |
| cons drop at H job loss exp (%) | -6.1281 | -5.8582 | -5.3668 | -5.0569 | -4.8231 |
| cons drop at H job loss rec (%) | -9.9832 | -9.6897 | -9.2255 | -8.8764 | -8.7331 |
| mean assets/monthly HH inc | 3.6205 | 4.0235 | 4.3185 | 4.7423 | 5.0353 |

### Mechanism counterfactuals (baseline parameters)

| moment | baseline | acyclical husband risk | acyclical job finding | no recession wage cut | acyclical own job loss |
|---|---|---|---|---|---|
| E/pop | 0.7172 | 0.7079 | 0.7191 | 0.7223 | 0.7195 |
| hours|E | 0.4035 | 0.3992 | 0.4024 | 0.4035 | 0.4034 |
| U rate | 0.0393 | 0.0389 | 0.0378 | 0.0395 | 0.0379 |
| quit/m exp | 0.0267 | 0.0283 | 0.0274 | 0.0264 | 0.0266 |
| quit/m rec | 0.0181 | 0.0201 | 0.0230 | 0.0184 | 0.0179 |
| E->nonE/m exp | 0.0400 | 0.0416 | 0.0407 | 0.0397 | 0.0398 |
| E->nonE/m rec | 0.0333 | 0.0354 | 0.0383 | 0.0336 | 0.0309 |
| dE/pop rec-exp (pts) | -1.6454 | -2.1072 | 0.3243 | -0.5539 | -0.9452 |
| wife share exp | 0.3333 | 0.3286 | 0.3329 | 0.3326 | 0.3337 |
| wife share rec | 0.3414 | 0.3193 | 0.3437 | 0.3532 | 0.3434 |
| wage gap (hourly ratio) | 0.7429 | 0.7446 | 0.7421 | 0.7433 | 0.7429 |
| share Lifecycle | 0.2963 | 0.2900 | 0.3021 | 0.2944 | 0.2977 |
| share PT | 0.2896 | 0.2875 | 0.2781 | 0.2900 | 0.2883 |
| share Career | 0.2029 | 0.1938 | 0.2035 | 0.2075 | 0.2042 |
| share NiLF | 0.2112 | 0.2288 | 0.2162 | 0.2081 | 0.2098 |
| HH income rec/exp - 1 (%) | -15.5393 | -13.1750 | -15.1258 | -2.1547 | -15.2799 |
| cons drop at H job loss exp (%) | -6.1281 | -6.4467 | -6.1517 | -6.5013 | -6.1390 |
| cons drop at H job loss rec (%) | -9.9832 | -7.9764 | -9.8362 | -8.8208 | -9.9036 |
| mean assets/monthly HH inc | 3.6205 | 3.5401 | 3.6207 | 3.5431 | 3.6207 |


Elapsed 3645s. Calibration file: output/final_calib_v4emb_full.json
