# Precautionary labor supply versus job hoarding (`output/final_calib_v4emb_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.717 | 0.0267 | 0.0181 | -0.86 | 5% | 49% | 63% | -1.65 | -2.11 |
| job finding falls 5% in recessions (ratio 0.95) | 0.719 | 0.0272 | 0.0215 | -0.57 | 19% | 24% | 44% | -0.19 | -0.74 |
| job finding falls 30% in recessions (ratio 0.70) | 0.716 | 0.0264 | 0.0159 | -1.05 | 4% | 58% | 70% | -2.60 | -3.11 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.722 | 0.0260 | 0.0165 | -0.95 | 13% | 45% | 66% | -1.30 | -2.11 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.721 | 0.0262 | 0.0170 | -0.92 | 11% | 47% | 65% | -1.38 | -2.11 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.717 | 0.0267 | 0.0181 | -0.86 | 5% | 49% | 63% | -1.65 | -2.11 |
| UI replacement 15% always | 0.726 | 0.0255 | 0.0168 | -0.86 | 11% | 51% | 67% | -1.46 | -2.14 |
| no assets | 0.709 | 0.0271 | 0.0176 | -0.95 | 12% | 48% | 65% | -1.39 | -2.28 |
| risk aversion 3 | 0.607 | 0.0733 | 0.0351 | -3.82 | 13% | 25% | 45% | +3.08 | +0.43 |
| longer recessions (persistence 0.95) | 0.720 | 0.0261 | 0.0165 | -0.96 | 4% | 50% | 58% | -2.76 | -3.28 |
