# Sensitivity of the precaution / hoarding split (`output/final_calib_v9n_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 9%, hoarding 64%, both off 69%, quit gap -0.87 points, sd log UE 0.0630, recession employment drop -1.67 (-1.92 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.10 | +0.21 | -0.005 | +0.0009 | -0.026 | -0.029 | -0.006 |
| UI cut in recessions (1 - ui_rec_mult) | +0.05 | +0.04 | -0.001 | +0.0000 | +0.002 | +0.000 | -0.001 |
| husband job-loss rate in recessions | +0.10 | -0.01 | -0.001 | +0.0000 | -0.002 | +0.000 | -0.002 |
| husband job-finding rate in recessions | -0.18 | +0.02 | +0.002 | -0.0001 | +0.004 | +0.000 | +0.002 |
| husband job-loss rate (both states) | +0.18 | -0.06 | -0.002 | -0.0001 | -0.003 | -0.007 | -0.005 |
| UI replacement (both states) | -0.23 | -0.20 | +0.002 | -0.0001 | -0.002 | -0.001 | +0.002 |
| wife own job loss in recessions (lam_u1) | -0.02 | -0.12 | +0.000 | +0.0000 | -0.044 | -0.046 | +0.001 |
| job-finding efficiency level (lam_f0) | -0.55 | -0.28 | -0.001 | +0.0000 | -0.010 | -0.012 | +0.042 |
| cost-shock sd (sd_kT) | -0.72 | -0.13 | -0.001 | +0.0003 | -0.014 | -0.006 | +0.053 |
| recession wage cut (1 - phi_rec) | +0.00 | +0.00 | +0.000 | +0.0000 | +0.000 | +0.000 | +0.000 |
| expected recession duration (1 / exit probability) | +0.02 | +0.43 | +0.004 | +0.0002 | -0.025 | -0.023 | +0.002 |
| asset limit a_max | -0.40 | -0.17 | +0.001 | +0.0001 | -0.003 | -0.004 | +0.003 |
| risk aversion gamma | -0.13 | -0.13 | +0.004 | +0.0003 | -0.008 | -0.011 | +0.009 |
