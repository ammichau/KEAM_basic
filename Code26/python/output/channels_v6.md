# Precautionary labor supply versus job hoarding (`output/final_calib_v6_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.659 | 0.0350 | 0.0239 | -1.11 | 29% | 48% | 71% | -1.61 | -2.93 |
| job finding falls 5% in recessions (ratio 0.95) | 0.661 | 0.0354 | 0.0284 | -0.70 | 37% | 17% | 54% | -0.28 | -1.46 |
| job finding falls 30% in recessions (ratio 0.70) | 0.658 | 0.0348 | 0.0209 | -1.39 | 24% | 58% | 77% | -2.71 | -4.17 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.664 | 0.0345 | 0.0217 | -1.28 | 39% | 43% | 75% | -0.93 | -2.93 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.664 | 0.0346 | 0.0223 | -1.23 | 36% | 40% | 74% | -1.10 | -2.93 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.659 | 0.0350 | 0.0239 | -1.11 | 29% | 48% | 71% | -1.61 | -2.93 |
| UI replacement 15% always | 0.675 | 0.0324 | 0.0215 | -1.09 | 32% | 51% | 72% | -1.58 | -3.00 |
| no assets | 0.666 | 0.0331 | 0.0203 | -1.28 | 42% | 45% | 74% | -1.20 | -3.28 |
| risk aversion 3 | 0.533 | 0.0831 | 0.0422 | -4.08 | 38% | 35% | 70% | +3.34 | -1.77 |
| longer recessions (persistence 0.95) | 0.662 | 0.0347 | 0.0218 | -1.28 | 25% | 46% | 65% | -2.30 | -3.86 |
