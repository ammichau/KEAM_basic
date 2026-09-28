# Sensitivity of the precaution / hoarding split (`output/final_calib_v4nb_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 21%, hoarding 56%, both off 78%, quit gap -1.11 points, sd log UE 0.0692, recession employment drop -1.67 (-2.81 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.25 | +0.21 | -0.005 | +0.0011 | -0.024 | -0.035 | -0.006 |
| UI cut in recessions (1 - ui_rec_mult) | +0.07 | -0.05 | -0.001 | +0.0001 | +0.004 | +0.000 | -0.001 |
| husband job-loss rate in recessions | +0.31 | -0.22 | -0.005 | -0.0001 | +0.019 | +0.000 | -0.007 |
| husband job-finding rate in recessions | -0.15 | +0.15 | +0.002 | +0.0002 | -0.018 | +0.000 | +0.003 |
| husband job-loss rate (both states) | +0.14 | -0.23 | +0.001 | -0.0001 | +0.000 | +0.003 | -0.012 |
| UI replacement (both states) | +0.11 | -0.10 | -0.001 | +0.0001 | -0.002 | -0.002 | +0.002 |
| wife own job loss in recessions (lam_u1) | +0.17 | +0.04 | +0.002 | -0.0001 | -0.052 | -0.055 | +0.004 |
| job-finding efficiency level (lam_f0) | +0.19 | -0.15 | -0.026 | -0.0012 | +0.039 | +0.030 | +0.041 |
| cost-shock sd (sd_kT) | +0.01 | -0.06 | -0.021 | -0.0003 | +0.009 | +0.009 | +0.047 |
| recession wage cut (1 - phi_rec) | +0.00 | +0.00 | +0.000 | +0.0000 | +0.000 | +0.000 | +0.000 |
| expected recession duration (1 / exit probability) | +0.18 | +0.29 | +0.003 | +0.0002 | -0.017 | -0.021 | +0.001 |
| asset limit a_max | +0.10 | -0.14 | -0.001 | -0.0001 | +0.004 | -0.003 | -0.002 |
| risk aversion gamma | +0.13 | -0.26 | -0.036 | +0.0000 | +0.033 | +0.019 | +0.063 |
