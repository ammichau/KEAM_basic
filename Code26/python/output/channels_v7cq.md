# Precautionary labor supply versus job hoarding (`output/final_calib_v7cq_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.660 | 0.0253 | 0.0184 | -0.69 | 13% | 50% | 66% | -1.73 | -1.59 |
| job finding falls 5% in recessions (ratio 0.95) | 0.664 | 0.0256 | 0.0214 | -0.43 | 26% | 19% | 44% | -0.21 | -0.08 |
| job finding falls 30% in recessions (ratio 0.70) | 0.657 | 0.0249 | 0.0154 | -0.95 | 9% | 64% | 75% | -3.33 | -3.19 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.662 | 0.0249 | 0.0172 | -0.77 | 21% | 44% | 69% | -1.85 | -1.59 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.661 | 0.0251 | 0.0177 | -0.75 | 19% | 47% | 68% | -1.76 | -1.59 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.660 | 0.0253 | 0.0184 | -0.69 | 13% | 50% | 66% | -1.73 | -1.59 |
| UI replacement 15% always | 0.666 | 0.0243 | 0.0176 | -0.67 | 14% | 51% | 67% | -1.73 | -1.60 |
| no assets | 0.660 | 0.0250 | 0.0184 | -0.66 | 12% | 48% | 65% | -1.53 | -1.65 |
| risk aversion 3 | 0.643 | 0.0411 | 0.0291 | -1.20 | 11% | 42% | 57% | -1.97 | -2.20 |
| longer recessions (persistence 0.95) | 0.657 | 0.0256 | 0.0177 | -0.79 | 13% | 47% | 64% | -3.03 | -2.82 |
