# Precautionary labor supply versus job hoarding (`output/final_calib_v4nb_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.679 | 0.0353 | 0.0242 | -1.11 | 21% | 56% | 78% | -1.67 | -2.81 |
| job finding falls 5% in recessions (ratio 0.95) | 0.680 | 0.0359 | 0.0295 | -0.65 | 37% | 25% | 63% | -0.25 | -1.27 |
| job finding falls 30% in recessions (ratio 0.70) | 0.677 | 0.0349 | 0.0208 | -1.41 | 14% | 65% | 83% | -2.80 | -4.14 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.685 | 0.0344 | 0.0218 | -1.25 | 30% | 47% | 81% | -1.00 | -2.81 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.683 | 0.0349 | 0.0226 | -1.22 | 28% | 49% | 80% | -1.22 | -2.81 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.679 | 0.0353 | 0.0242 | -1.11 | 21% | 56% | 78% | -1.67 | -2.81 |
| UI replacement 15% always | 0.690 | 0.0332 | 0.0223 | -1.09 | 23% | 54% | 77% | -1.52 | -2.82 |
| no assets | 0.685 | 0.0335 | 0.0212 | -1.23 | 31% | 50% | 78% | -0.94 | -2.98 |
| risk aversion 3 | 0.514 | 0.1135 | 0.0690 | -4.46 | 29% | 59% | 90% | +0.75 | -2.14 |
| longer recessions (persistence 0.95) | 0.681 | 0.0350 | 0.0226 | -1.24 | 22% | 49% | 75% | -2.36 | -4.03 |
