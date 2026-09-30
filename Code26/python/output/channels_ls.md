# Precautionary labor supply versus job hoarding (`output/final_calib_ls_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.664 | 0.0332 | 0.0250 | -0.82 | 19% | 39% | 52% | -1.69 | -2.23 |
| job finding falls 5% in recessions (ratio 0.95) | 0.664 | 0.0337 | 0.0275 | -0.62 | 21% | 19% | 36% | -0.74 | -1.16 |
| job finding falls 30% in recessions (ratio 0.70) | 0.662 | 0.0327 | 0.0208 | -1.19 | 18% | 58% | 67% | -2.87 | -3.69 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.668 | 0.0327 | 0.0232 | -0.95 | 30% | 39% | 59% | -1.24 | -2.23 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.667 | 0.0327 | 0.0232 | -0.95 | 30% | 34% | 59% | -1.20 | -2.23 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.668 | 0.0327 | 0.0227 | -1.00 | 33% | 36% | 61% | -1.11 | -2.23 |
| UI replacement 15% always | 0.679 | 0.0306 | 0.0218 | -0.88 | 23% | 34% | 57% | -1.09 | -1.98 |
| no assets | 0.667 | 0.0316 | 0.0230 | -0.86 | 24% | 42% | 55% | -1.82 | -2.50 |
| risk aversion 3 | 0.552 | 0.0755 | 0.0465 | -2.90 | 24% | 36% | 61% | +0.38 | -1.93 |
| longer recessions (persistence 0.95) | 0.664 | 0.0330 | 0.0239 | -0.91 | 14% | 38% | 46% | -2.76 | -3.25 |
