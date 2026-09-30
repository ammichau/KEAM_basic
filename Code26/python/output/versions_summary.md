Objective: weighted sum of squared deviations over the 13 targets of `keam/final/calibrate.py` (100-type moments). Target columns: relative deviation from the data, except dE dev (the recession employment drop, deviation in points). Precaution / hoarding: share of the recession fall in the monthly quit rate removed when the husband's risk / the wife's job-finding efficiency is made acyclical (`scripts/channels.py`). dE: recession minus expansion employment rate (points), baseline and with acyclical husband risk.

| version | objective | E/pop | hours | LC | PT | career | NiLF | quit exp | quit rec | exit exp | exit rec | dE dev (pts) | wage gap | sd log UE | λ_f rec/exp | quit gap | precaution | hoarding | both off | dE rec | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| adopted iid | 0.449 | +7% | +2% | -0% | -11% | -10% | +22% | -2% | -11% | +1% | -7% | +0.01 | +4% | -41% | 0.85 | -0.82 | 19% | 39% | 52% | -1.69 | -2.23 |
| UI cut | 0.307 | +6% | +2% | -2% | -7% | -16% | +25% | +0% | -12% | +2% | -5% | -0.02 | +4% | -29% | 0.85 | -0.95 | 29% | 34% | 58% | -1.72 | -2.61 |
| 7 wage types | 0.343 | +8% | +3% | +0% | -7% | -13% | +21% | -3% | -13% | +1% | -8% | +0.05 | +4% | -32% | 0.85 | -0.87 | 21% | 38% | 54% | -1.65 | -2.26 |
| version 3 | 0.759 | +8% | +2% | -4% | -8% | -8% | +22% | -3% | -11% | -0% | -2% | +0.00 | +4% | -57% | 0.90 | -0.81 | 29% | 25% | 50% | -1.70 | -2.53 |
| version 4 | 0.256 | +7% | +4% | +2% | -14% | -5% | +19% | -5% | -18% | -2% | -6% | -0.18 | +5% | -18% | 0.85 | -0.93 | 33% | 34% | 59% | -1.88 | -2.76 |
| version 4c | 0.803 | +4% | -4% | -21% | +18% | -46% | +46% | +26% | +4% | +23% | +5% | -0.25 | -4% | -8% | 0.80 | -1.16 | 31% | 44% | 66% | -1.65 | -2.88 |
| version 5 (log utility) | 0.396 | -2% | +8% | -9% | -2% | -23% | +34% | -7% | -17% | +3% | -9% | -0.15 | +15% | +17% | 0.80 | -0.84 | 7% | 56% | 66% | -1.85 | -2.10 |
| version 5b (log utility; 2nd polish) | 0.306 | -0% | +10% | +1% | -14% | -13% | +28% | -7% | -18% | +3% | -11% | -0.09 | +16% | +12% | 0.80 | -0.85 | 7% | 62% | 70% | -1.79 | -2.02 |
| version 6 (persistent shock) | 0.260 | +6% | +2% | +15% | -19% | -23% | +22% | +3% | -15% | +7% | -2% | +0.09 | +8% | +5% | 0.80 | -1.11 | 29% | 48% | 71% | -1.61 | -2.93 |

* version 4c: objective 0.803; precaution 31% vs hoarding 44% (13 points apart, rule: within 10); largest cyclical-moment deviation 15% (dE/pop rec-exp (pts)).
* version 5b (log utility; 2nd polish): objective 0.306; precaution 7% vs hoarding 62% (54 points apart, rule: within 10); largest cyclical-moment deviation 18% (quit/m rec).
* Carried forward: **version 4c** (`v4c`): no candidate satisfies the 10-point rule; this is the one closest to parity.
* Log utility and the precautionary channel: at the version-4c parameters (γ = 2) precaution is 31% and hoarding 44%; imposing γ = 1 without recalibrating gives 9% / 55% (`output/channels_v4c_gamma1.json`; employment 0.79, because the calibrated cost levels are in γ = 2 utility units), and the recalibrated log-utility version gives 7% / 62%. The fall is a property of the preferences, not of the recalibration. Fit: objective 0.306 versus 0.803; the largest deviations of the log-utility version are share NiLF +28%, quit/m rec -18%, wage gap (hourly ratio) +16%, share PT -14%.
