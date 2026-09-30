# Precautionary labor supply versus job hoarding (`output/final_calib_v4e_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.683 | 0.0354 | 0.0235 | -1.20 | 16% | 51% | 71% | -1.76 | -2.91 |
| job finding falls 5% in recessions (ratio 0.95) | 0.684 | 0.0361 | 0.0285 | -0.75 | 34% | 22% | 55% | -0.30 | -0.97 |
| job finding falls 30% in recessions (ratio 0.70) | 0.681 | 0.0352 | 0.0202 | -1.49 | 13% | 61% | 77% | -2.74 | -3.79 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.689 | 0.0344 | 0.0211 | -1.33 | 24% | 44% | 74% | -1.14 | -2.91 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.687 | 0.0349 | 0.0220 | -1.29 | 22% | 46% | 74% | -1.24 | -2.91 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.683 | 0.0354 | 0.0235 | -1.20 | 16% | 51% | 71% | -1.76 | -2.91 |
| UI replacement 15% always | 0.694 | 0.0334 | 0.0217 | -1.17 | 17% | 49% | 71% | -1.56 | -2.69 |
| no assets | 0.688 | 0.0336 | 0.0210 | -1.26 | 23% | 48% | 71% | -1.36 | -3.05 |
| risk aversion 3 | 0.526 | 0.1123 | 0.0585 | -5.37 | 15% | 30% | 58% | +1.53 | -1.15 |
| longer recessions (persistence 0.95) | 0.685 | 0.0351 | 0.0218 | -1.33 | 14% | 47% | 66% | -2.52 | -3.94 |
