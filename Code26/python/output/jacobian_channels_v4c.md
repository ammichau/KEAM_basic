# Sensitivity of the precaution / hoarding split (`output/final_calib_v4c_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 31%, hoarding 44%, both off 66%, quit gap -1.16 points, sd log UE 0.0686, recession employment drop -1.65 (-2.88 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.16 | +0.24 | -0.005 | +0.0006 | -0.016 | -0.025 | -0.006 |
| UI cut in recessions (1 - ui_rec_mult) | +0.13 | -0.10 | -0.002 | +0.0000 | +0.012 | +0.000 | -0.003 |
| husband job-loss rate in recessions | +0.22 | -0.16 | -0.004 | -0.0001 | +0.017 | +0.000 | -0.005 |
| husband job-finding rate in recessions | -0.17 | +0.19 | +0.003 | -0.0002 | -0.006 | +0.000 | +0.005 |
| husband job-loss rate (both states) | +0.08 | -0.20 | -0.001 | +0.0000 | +0.013 | -0.008 | -0.009 |
| UI replacement (both states) | +0.20 | +0.12 | -0.004 | -0.0002 | +0.005 | +0.000 | +0.003 |
| wife own job loss in recessions (lam_u1) | +0.12 | +0.01 | +0.000 | +0.0004 | -0.050 | -0.048 | +0.003 |
| job-finding efficiency level (lam_f0) | -0.20 | +0.27 | -0.025 | -0.0004 | +0.013 | -0.022 | +0.033 |
| cost-shock sd (sd_kT) | -0.21 | +0.37 | -0.021 | +0.0003 | -0.012 | -0.032 | +0.036 |
| recession wage cut (1 - phi_rec) | +0.07 | -0.17 | +0.000 | -0.0001 | -0.006 | -0.005 | +0.003 |
| asset limit a_max | -0.13 | +0.06 | -0.000 | -0.0001 | -0.001 | -0.001 | -0.002 |
| risk aversion gamma | +0.04 | +0.02 | -0.038 | +0.0007 | +0.055 | +0.004 | +0.029 |
| expected recession duration (1 / exit probability) | +0.19 | +0.14 | +0.002 | +0.0001 | -0.014 | -0.027 | +0.000 |
