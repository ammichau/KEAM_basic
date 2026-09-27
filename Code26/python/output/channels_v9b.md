# Precautionary labor supply versus job hoarding (`output/final_calib_v9b_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.631 | 0.0328 | 0.0254 | -0.73 | 9% | 68% | 80% | -1.88 | -1.77 |
| job finding falls 5% in recessions (ratio 0.95) | 0.636 | 0.0332 | 0.0296 | -0.36 | 22% | 35% | 59% | +0.08 | +0.19 |
| job finding falls 30% in recessions (ratio 0.70) | 0.628 | 0.0325 | 0.0229 | -0.96 | 7% | 76% | 85% | -3.26 | -3.17 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.632 | 0.0326 | 0.0251 | -0.76 | 12% | 66% | 80% | -1.90 | -1.77 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.632 | 0.0326 | 0.0250 | -0.76 | 13% | 65% | 81% | -1.85 | -1.77 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.631 | 0.0328 | 0.0254 | -0.73 | 9% | 68% | 80% | -1.88 | -1.77 |
| UI replacement 15% always | 0.631 | 0.0331 | 0.0258 | -0.73 | 10% | 74% | 77% | -2.01 | -2.02 |
| no assets | 0.608 | 0.0344 | 0.0269 | -0.75 | 11% | 70% | 90% | -1.98 | -2.03 |
| risk aversion 3 | 0.616 | 0.0377 | 0.0301 | -0.75 | 10% | 84% | 86% | -2.65 | -2.50 |
| longer recessions (persistence 0.95) | 0.628 | 0.0333 | 0.0249 | -0.84 | 8% | 72% | 81% | -2.89 | -2.99 |
