# Sensitivity of the precaution / hoarding split (`output/final_calib_v4_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 33%, hoarding 34%, both off 59%, quit gap -0.93 points, sd log UE 0.0562, recession employment drop -1.88 (-2.76 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.02 | +0.23 | -0.003 | +0.0007 | -0.010 | -0.018 | -0.003 |
| UI cut in recessions (1 - ui_rec_mult) | +0.11 | +0.08 | -0.002 | -0.0004 | +0.007 | +0.000 | -0.002 |
| husband job-loss rate in recessions | +0.24 | +0.06 | -0.003 | -0.0000 | +0.014 | +0.000 | -0.004 |
| husband job-finding rate in recessions | -0.37 | +0.05 | +0.005 | -0.0003 | -0.012 | +0.000 | +0.007 |
| husband job-loss rate (both states) | +0.15 | +0.19 | -0.002 | -0.0001 | +0.018 | -0.004 | -0.009 |
| UI replacement (both states) | -0.08 | +0.12 | +0.000 | -0.0003 | -0.001 | -0.005 | +0.007 |
| wife own job loss in recessions (lam_u1) | -0.20 | -0.10 | +0.003 | +0.0001 | -0.057 | -0.048 | +0.006 |
| job-finding efficiency level (lam_f0) | -0.27 | +0.33 | -0.018 | -0.0005 | +0.000 | -0.022 | +0.039 |
| cost-shock sd (sd_kT) | -0.53 | +0.40 | -0.013 | +0.0003 | -0.014 | -0.026 | +0.040 |
| recession wage cut (1 - phi_rec) | -0.22 | -0.27 | +0.004 | -0.0005 | -0.019 | -0.008 | +0.005 |
| expected recession duration (1 / exit probability) | +0.22 | +0.24 | +0.004 | -0.0000 | -0.020 | -0.023 | +0.002 |
| asset limit a_max | -0.15 | +0.17 | +0.001 | -0.0004 | +0.007 | +0.000 | +0.000 |
| risk aversion gamma | +0.02 | +0.22 | -0.028 | -0.0015 | +0.048 | +0.000 | +0.034 |
