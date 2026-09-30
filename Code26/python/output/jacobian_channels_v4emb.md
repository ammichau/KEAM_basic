# Sensitivity of the precaution / hoarding split (`output/final_calib_v4emb_full.json`, +10% steps, parameters otherwise fixed)

Baseline: precaution 5%, hoarding 49%, both off 63%, quit gap -0.86 points, sd log UE 0.0775, recession employment drop -1.65 (-2.11 without cyclical husband risk).

Entries: change per +1% of the named quantity (shares in percentage points; quit gap and employment drop in percentage points of the rate; sd log UE in units).

| quantity perturbed | precaution share (pp) | hoarding share (pp) | quit gap (pp) | sd log UE | dE (pts) | dE acyc. husband (pts) | quit rec (pp) |
|---|---|---|---|---|---|---|---|
| job-finding fall in recessions (1 - lam_f ratio) | -0.03 | +0.26 | -0.005 | +0.0012 | -0.020 | -0.022 | -0.005 |
| UI cut in recessions (1 - ui_rec_mult) | +0.03 | +0.03 | -0.000 | +0.0001 | +0.002 | +0.000 | -0.000 |
| husband job-loss rate in recessions | +0.15 | -0.17 | -0.001 | -0.0000 | +0.012 | +0.000 | -0.004 |
| husband job-finding rate in recessions | -0.06 | +0.05 | +0.001 | +0.0001 | -0.002 | +0.000 | +0.002 |
| husband job-loss rate (both states) | +0.40 | -0.01 | -0.000 | -0.0001 | +0.013 | -0.002 | -0.011 |
| UI replacement (both states) | +0.09 | +0.19 | -0.002 | +0.0000 | +0.004 | -0.001 | +0.003 |
| wife own job loss in recessions (lam_u1) | -0.01 | +0.09 | +0.000 | +0.0000 | -0.040 | -0.044 | +0.001 |
| job-finding efficiency level (lam_f0) | +0.22 | -0.00 | -0.021 | +0.0003 | +0.008 | -0.011 | +0.033 |
| cost-shock sd (sd_kT) | +0.42 | +0.18 | -0.019 | +0.0003 | +0.007 | -0.013 | +0.034 |
| recession wage cut (1 - phi_rec) | -0.06 | -0.05 | -0.000 | +0.0001 | -0.011 | -0.006 | -0.000 |
| expected recession duration (1 / exit probability) | +0.15 | +0.25 | +0.002 | +0.0002 | -0.024 | -0.024 | +0.000 |
| asset limit a_max | +0.15 | -0.01 | -0.001 | +0.0004 | +0.006 | +0.001 | +0.001 |
| risk aversion gamma | +0.29 | -0.18 | -0.030 | +0.0007 | +0.042 | -0.000 | +0.033 |
