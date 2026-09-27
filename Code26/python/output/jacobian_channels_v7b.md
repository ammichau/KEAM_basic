# Sensitivity of the precaution / hoarding split (`output/final_calib_v7b_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 28%, hoarding 47%, both off 64%, quit gap -0.88 points, sd log UE 0.0714, recession employment drop -1.69 (-1.84 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.05 | +0.34 | -0.006 | +0.0006 | -0.021 | -0.035 | -0.006 |
| UI cut in recessions (1 - ui_rec_mult) | +0.02 | -0.08 | -0.000 | -0.0001 | -0.004 | +0.000 | -0.001 |
| husband job-loss rate in recessions | +0.14 | -0.04 | -0.002 | -0.0001 | -0.002 | +0.000 | -0.003 |
| husband job-finding rate in recessions | -0.24 | +0.05 | +0.003 | -0.0002 | -0.003 | +0.000 | +0.004 |
| husband job-loss rate (both states) | +0.13 | +0.07 | -0.000 | -0.0000 | -0.005 | -0.005 | -0.005 |
| UI replacement (both states) | -0.01 | -0.02 | -0.000 | -0.0003 | -0.002 | +0.006 | +0.000 |
| wife own job loss in recessions (lam_u1) | -0.06 | -0.07 | +0.001 | -0.0001 | -0.052 | -0.048 | +0.002 |
| job-finding efficiency level (lam_f0) | -0.47 | -0.09 | -0.017 | -0.0012 | -0.025 | -0.032 | +0.041 |
| cost-shock sd (sd_kT) | -0.56 | -0.22 | -0.017 | +0.0002 | -0.029 | -0.044 | +0.038 |
| recession wage cut (1 - phi_rec) | -0.06 | -0.04 | +0.001 | -0.0005 | -0.014 | -0.009 | +0.002 |
| expected recession duration (1 / exit probability) | +0.09 | +0.28 | +0.005 | +0.0001 | -0.026 | -0.025 | +0.003 |
| asset limit a_max | -0.02 | +0.07 | +0.000 | +0.0001 | -0.001 | -0.004 | -0.002 |
| risk aversion gamma | -0.36 | -0.37 | -0.018 | -0.0000 | -0.010 | -0.012 | +0.025 |
