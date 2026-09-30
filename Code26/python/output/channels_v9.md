# Precautionary labor supply versus job hoarding (`output/final_calib_v9_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.647 | 0.0338 | 0.0254 | -0.84 | 6% | 58% | 65% | -1.72 | -1.74 |
| job finding falls 5% in recessions (ratio 0.95) | 0.651 | 0.0342 | 0.0294 | -0.47 | 10% | 26% | 38% | +0.02 | -0.06 |
| job finding falls 30% in recessions (ratio 0.70) | 0.645 | 0.0335 | 0.0227 | -1.07 | 4% | 67% | 73% | -3.04 | -3.07 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.649 | 0.0338 | 0.0249 | -0.89 | 11% | 58% | 67% | -1.63 | -1.74 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.649 | 0.0337 | 0.0249 | -0.88 | 10% | 58% | 67% | -1.60 | -1.74 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.647 | 0.0338 | 0.0254 | -0.84 | 6% | 58% | 65% | -1.72 | -1.74 |
| UI replacement 15% always | 0.646 | 0.0347 | 0.0259 | -0.88 | 5% | 60% | 60% | -1.82 | -1.95 |
| no assets | 0.625 | 0.0358 | 0.0265 | -0.94 | 11% | 61% | 74% | -1.78 | -2.17 |
| risk aversion 3 | 0.630 | 0.0390 | 0.0303 | -0.87 | 4% | 76% | 69% | -2.61 | -2.55 |
| longer recessions (persistence 0.95) | 0.645 | 0.0344 | 0.0249 | -0.94 | 7% | 60% | 66% | -2.76 | -2.98 |
