# Precautionary labor supply versus job hoarding (`output/final_calib_v7cmba_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.659 | 0.0245 | 0.0179 | -0.66 | 10% | 64% | 79% | -1.69 | -1.70 |
| job finding falls 5% in recessions (ratio 0.95) | 0.662 | 0.0248 | 0.0214 | -0.34 | 26% | 30% | 60% | -0.11 | -0.11 |
| job finding falls 30% in recessions (ratio 0.70) | 0.656 | 0.0242 | 0.0152 | -0.90 | 4% | 73% | 85% | -2.95 | -3.03 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.663 | 0.0239 | 0.0171 | -0.68 | 13% | 56% | 80% | -1.57 | -1.70 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.662 | 0.0242 | 0.0175 | -0.67 | 11% | 59% | 79% | -1.51 | -1.70 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.659 | 0.0245 | 0.0179 | -0.66 | 10% | 64% | 79% | -1.69 | -1.70 |
| UI replacement 15% always | 0.666 | 0.0236 | 0.0169 | -0.67 | 10% | 68% | 78% | -1.77 | -1.67 |
| no assets | 0.664 | 0.0237 | 0.0166 | -0.71 | 6% | 60% | 73% | -1.82 | -1.77 |
| risk aversion 3 | 0.668 | 0.0398 | 0.0280 | -1.18 | 10% | 51% | 62% | -1.88 | -2.13 |
| longer recessions (persistence 0.95) | 0.658 | 0.0244 | 0.0167 | -0.77 | 11% | 60% | 75% | -2.85 | -2.71 |
