# Sensitivity of the precaution / hoarding split (`output/final_calib_v4e_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 16%, hoarding 51%, both off 71%, quit gap -1.20 points, sd log UE 0.0705, recession employment drop -1.76 (-2.91 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.01 | +0.29 | -0.007 | +0.0012 | -0.022 | -0.019 | -0.008 |
| UI cut in recessions (1 - ui_rec_mult) | +0.04 | -0.13 | -0.001 | -0.0000 | +0.009 | +0.000 | -0.001 |
| husband job-loss rate in recessions | +0.24 | -0.23 | -0.004 | -0.0002 | +0.020 | +0.000 | -0.007 |
| husband job-finding rate in recessions | -0.13 | +0.13 | +0.002 | -0.0000 | -0.006 | +0.000 | +0.003 |
| husband job-loss rate (both states) | +0.06 | -0.29 | +0.001 | -0.0001 | +0.012 | +0.016 | -0.012 |
| UI replacement (both states) | +0.03 | -0.08 | -0.001 | -0.0002 | -0.003 | +0.000 | +0.004 |
| wife own job loss in recessions (lam_u1) | +0.05 | -0.07 | +0.002 | +0.0004 | -0.052 | -0.053 | +0.003 |
| job-finding efficiency level (lam_f0) | +0.26 | -0.30 | -0.029 | -0.0013 | +0.040 | +0.025 | +0.039 |
| cost-shock sd (sd_kT) | +0.29 | -0.34 | -0.026 | -0.0003 | +0.014 | +0.012 | +0.042 |
| recession wage cut (1 - phi_rec) | -0.14 | -0.08 | +0.000 | +0.0001 | -0.014 | -0.006 | +0.001 |
| expected recession duration (1 / exit probability) | +0.05 | +0.08 | +0.005 | +0.0000 | -0.020 | -0.023 | +0.002 |
| asset limit a_max | +0.09 | -0.05 | -0.001 | +0.0003 | +0.002 | +0.002 | -0.001 |
| risk aversion gamma | +0.15 | -0.46 | -0.048 | -0.0005 | +0.057 | +0.054 | +0.050 |
