# Precautionary labor supply versus job hoarding (`output/final_calib_v5_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.605 | 0.0316 | 0.0232 | -0.84 | 7% | 56% | 66% | -1.85 | -2.10 |
| job finding falls 5% in recessions (ratio 0.95) | 0.607 | 0.0324 | 0.0272 | -0.52 | 12% | 29% | 44% | -0.04 | -0.21 |
| job finding falls 30% in recessions (ratio 0.70) | 0.602 | 0.0312 | 0.0207 | -1.05 | 8% | 65% | 73% | -3.18 | -3.44 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.607 | 0.0313 | 0.0226 | -0.87 | 10% | 57% | 67% | -1.84 | -2.10 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.606 | 0.0314 | 0.0226 | -0.89 | 11% | 55% | 68% | -1.73 | -2.10 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.605 | 0.0316 | 0.0232 | -0.84 | 7% | 56% | 66% | -1.85 | -2.10 |
| UI replacement 15% always | 0.610 | 0.0306 | 0.0225 | -0.81 | 8% | 58% | 68% | -1.92 | -2.11 |
| no assets | 0.604 | 0.0315 | 0.0230 | -0.85 | 8% | 58% | 67% | -1.74 | -2.06 |
| risk aversion 3 | 0.340 | 0.1870 | 0.1365 | -5.05 | 53% | 47% | 88% | -2.34 | -5.52 |
| longer recessions (persistence 0.95) | 0.602 | 0.0320 | 0.0221 | -0.98 | 4% | 50% | 58% | -2.86 | -2.98 |
