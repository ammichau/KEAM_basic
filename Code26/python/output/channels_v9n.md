# Precautionary labor supply versus job hoarding (`output/final_calib_v9n_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.682 | 0.0357 | 0.0270 | -0.87 | 9% | 64% | 69% | -1.67 | -1.92 |
| job finding falls 5% in recessions (ratio 0.95) | 0.685 | 0.0361 | 0.0318 | -0.44 | 9% | 28% | 39% | +0.21 | +0.02 |
| job finding falls 30% in recessions (ratio 0.70) | 0.679 | 0.0353 | 0.0235 | -1.19 | 8% | 73% | 77% | -3.24 | -3.57 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.684 | 0.0355 | 0.0261 | -0.94 | 16% | 60% | 71% | -1.42 | -1.92 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.683 | 0.0356 | 0.0265 | -0.91 | 13% | 60% | 71% | -1.53 | -1.92 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.682 | 0.0357 | 0.0270 | -0.87 | 9% | 64% | 69% | -1.67 | -1.92 |
| UI replacement 15% always | 0.682 | 0.0360 | 0.0269 | -0.90 | 12% | 62% | 69% | -1.52 | -2.00 |
| no assets | 0.662 | 0.0384 | 0.0270 | -1.14 | 20% | 61% | 74% | -1.38 | -2.18 |
| risk aversion 3 | 0.661 | 0.0386 | 0.0305 | -0.82 | 4% | 67% | 73% | -1.68 | -2.18 |
| longer recessions (persistence 0.95) | 0.682 | 0.0356 | 0.0262 | -0.94 | 5% | 63% | 67% | -2.68 | -2.97 |
