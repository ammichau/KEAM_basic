# Sensitivity of the precaution / hoarding split (`output/final_calib_v7cmb_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 6%, hoarding 62%, both off 77%, quit gap -0.63 points, sd log UE 0.0743, recession employment drop -1.62 (-1.68 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.03 | +0.41 | -0.008 | +0.0007 | -0.030 | -0.031 | -0.008 |
| UI cut in recessions (1 - ui_rec_mult) | +0.01 | +0.15 | -0.000 | -0.0002 | +0.006 | +0.000 | -0.001 |
| husband job-loss rate in recessions | +0.19 | -0.15 | -0.001 | -0.0004 | +0.004 | +0.000 | -0.004 |
| husband job-finding rate in recessions | -0.05 | +0.21 | +0.000 | +0.0002 | +0.002 | +0.000 | +0.001 |
| husband job-loss rate (both states) | +0.38 | +0.06 | +0.000 | -0.0003 | +0.001 | +0.008 | -0.007 |
| UI replacement (both states) | +0.05 | +0.01 | -0.001 | +0.0002 | -0.002 | -0.001 | +0.002 |
| wife own job loss in recessions (lam_u1) | +0.14 | +0.10 | -0.001 | -0.0001 | -0.038 | -0.041 | +0.000 |
| job-finding efficiency level (lam_f0) | +0.28 | +0.09 | -0.004 | -0.0005 | +0.029 | +0.030 | +0.029 |
| cost-shock sd (sd_kT) | +0.43 | -0.06 | -0.002 | +0.0002 | +0.025 | +0.024 | +0.027 |
| recession wage cut (1 - phi_rec) | +0.05 | +0.25 | -0.005 | -0.0000 | -0.018 | -0.022 | -0.004 |
| expected recession duration (1 / exit probability) | +0.12 | +0.33 | +0.002 | +0.0001 | -0.027 | -0.029 | +0.001 |
| asset limit a_max | +0.06 | +0.12 | -0.000 | -0.0001 | +0.004 | +0.005 | -0.002 |
| risk aversion gamma | +0.11 | +0.29 | -0.007 | +0.0003 | +0.000 | -0.001 | +0.019 |
