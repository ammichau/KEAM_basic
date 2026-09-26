# Robustness of the final model (calibrated parameters held fixed)

Calibration: `output/final_calib_full.json`; returns-to-experience scale x1.38 from `output/final_results_full.json`. Quit gap = 100 x (monthly quit rate in recessions - in expansions), percentage points.

## Baseline moments by variant

| moment | baseline | no assets (a_max 0.01, 5 points) | asset grid 40 points | asset grid 40 points, a_max 30 | hours grid 40 points | hours grid 40, h_min 0.025 | U threshold s_bar 0.10 | U threshold s_bar 0.50 | phi_rec_H = 1 (no recession cut in husband income) |
|---|---|---|---|---|---|---|---|---|---|
| E/pop | 0.6644 | 0.6654 | 0.6701 | 0.6673 | 0.6650 | 0.6638 | 0.6644 | 0.6644 | 0.6575 |
| hours|E | 0.4141 | 0.4033 | 0.4211 | 0.4195 | 0.4143 | 0.4128 | 0.4141 | 0.4141 | 0.4096 |
| U rate | 0.0459 | 0.0491 | 0.0475 | 0.0466 | 0.0459 | 0.0455 | 0.1360 | 0.0182 | 0.0464 |
| quit/m exp | 0.0311 | 0.0299 | 0.0305 | 0.0309 | 0.0311 | 0.0315 | 0.0311 | 0.0311 | 0.0323 |
| quit/m rec | 0.0237 | 0.0220 | 0.0234 | 0.0234 | 0.0236 | 0.0241 | 0.0237 | 0.0237 | 0.0262 |
| E->nonE/m exp | 0.0485 | 0.0472 | 0.0478 | 0.0482 | 0.0484 | 0.0489 | 0.0485 | 0.0485 | 0.0496 |
| E->nonE/m rec | 0.0431 | 0.0415 | 0.0428 | 0.0428 | 0.0431 | 0.0434 | 0.0431 | 0.0431 | 0.0456 |
| dE/pop rec-exp (pts) | -1.5984 | -1.6450 | -1.5583 | -1.4847 | -1.5945 | -1.6125 | -1.5984 | -1.5984 | -1.7828 |
| wage gap (hourly ratio) | 0.7477 | 0.7413 | 0.7484 | 0.7487 | 0.7476 | 0.7464 | 0.7477 | 0.7477 | 0.7308 |
| wife share exp | 0.3370 | 0.3283 | 0.3423 | 0.3416 | 0.3367 | 0.3353 | 0.3370 | 0.3370 | 0.3342 |
| share Lifecycle | 0.3038 | 0.2904 | 0.3115 | 0.3069 | 0.3040 | 0.3006 | 0.3038 | 0.3038 | 0.2933 |
| share PT | 0.2531 | 0.2971 | 0.2490 | 0.2498 | 0.2504 | 0.2531 | 0.2531 | 0.2531 | 0.2612 |
| share Career | 0.1833 | 0.1671 | 0.1865 | 0.1862 | 0.1867 | 0.1810 | 0.1833 | 0.1833 | 0.1727 |
| share NiLF | 0.2598 | 0.2454 | 0.2531 | 0.2571 | 0.2590 | 0.2652 | 0.2598 | 0.2598 | 0.2727 |
| cons drop at H job loss rec (%) | -7.1846 | -7.2633 | -6.5630 | -6.9266 | -7.1510 | -7.1557 | -7.1846 | -7.1846 | -7.2243 |
| mean assets/monthly HH inc | 1.4722 | 0.0077 | 2.6927 | 2.4161 | 1.4725 | 1.4691 | 1.4722 | 1.4722 | 1.4590 |

## Mechanism and experiment by variant

| variant | quit gap | dE/pop rec-exp | acyclical H: quit gap | acyclical H: dE/pop rec-exp | RoE: E/pop | RoE: quit gap | RoE: dE/pop rec-exp |
|---|---|---|---|---|---|---|---|
| baseline | -0.742 | -1.598 | -0.616 | -1.998 | 0.728 | -0.524 | -1.556 |
| no assets (a_max 0.01, 5 points) | -0.786 | -1.645 | -0.590 | -2.229 | 0.728 | -0.518 | -1.849 |
| asset grid 40 points | -0.718 | -1.558 | -0.620 | -1.858 | 0.734 | -0.491 | -1.649 |
| asset grid 40 points, a_max 30 | -0.752 | -1.485 | -0.633 | -1.904 | 0.731 | -0.509 | -1.595 |
| hours grid 40 points | -0.744 | -1.595 | -0.606 | -1.997 | 0.729 | -0.526 | -1.569 |
| hours grid 40, h_min 0.025 | -0.747 | -1.612 | -0.605 | -2.074 | 0.727 | -0.543 | -1.582 |
| U threshold s_bar 0.10 | -0.742 | -1.598 | -0.616 | -1.998 | 0.728 | -0.524 | -1.556 |
| U threshold s_bar 0.50 | -0.742 | -1.598 | -0.616 | -1.998 | 0.728 | -0.524 | -1.556 |
| phi_rec_H = 1 (no recession cut in husband income) | -0.614 | -1.783 | -0.523 | -2.121 | 0.721 | -0.456 | -1.711 |

Elapsed 3622s.
