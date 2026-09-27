# Sensitivity of the precaution / hoarding split (`output/final_calib_v9_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 6%, hoarding 58%, both off 65%, quit gap -0.84 points, sd log UE 0.0631, recession employment drop -1.72 (-1.74 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.03 | +0.21 | -0.004 | +0.0008 | -0.023 | -0.027 | -0.005 |
| UI cut in recessions (1 - ui_rec_mult) | +0.02 | +0.01 | -0.000 | +0.0000 | +0.005 | +0.000 | +0.000 |
| husband job-loss rate in recessions | +0.17 | +0.05 | -0.002 | +0.0001 | +0.001 | +0.000 | -0.001 |
| husband job-finding rate in recessions | -0.08 | +0.09 | +0.001 | -0.0002 | +0.000 | +0.000 | +0.001 |
| husband job-loss rate (both states) | +0.24 | +0.04 | -0.003 | +0.0001 | +0.003 | +0.003 | -0.005 |
| UI replacement (both states) | +0.13 | -0.09 | +0.000 | -0.0000 | +0.001 | +0.005 | -0.001 |
| wife own job loss in recessions (lam_u1) | +0.08 | +0.06 | -0.001 | +0.0001 | -0.036 | -0.042 | +0.000 |
| job-finding efficiency level (lam_f0) | +0.06 | +0.55 | -0.015 | +0.0009 | -0.033 | -0.041 | +0.044 |
| cost-shock sd (sd_kT) | -0.15 | +0.71 | -0.018 | -0.0004 | -0.039 | -0.032 | +0.054 |
| recession wage cut (1 - phi_rec) | +0.07 | -0.08 | -0.001 | +0.0001 | -0.009 | -0.009 | +0.001 |
| expected recession duration (1 / exit probability) | +0.05 | +0.27 | +0.002 | +0.0002 | -0.027 | -0.033 | +0.002 |
| asset limit a_max | +0.04 | +0.09 | -0.001 | +0.0002 | -0.004 | -0.006 | +0.001 |
| risk aversion gamma | +0.13 | +0.39 | -0.002 | +0.0004 | -0.034 | -0.022 | +0.010 |
