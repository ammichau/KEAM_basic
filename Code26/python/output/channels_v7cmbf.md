# Precautionary labor supply versus job hoarding (`output/final_calib_v7cmbf_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.664 | 0.0238 | 0.0176 | -0.62 | 6% | 63% | 77% | -1.61 | -1.65 |
| job finding falls 5% in recessions (ratio 0.95) | 0.667 | 0.0241 | 0.0208 | -0.33 | 22% | 31% | 57% | -0.05 | -0.05 |
| job finding falls 30% in recessions (ratio 0.70) | 0.662 | 0.0234 | 0.0147 | -0.87 | 3% | 74% | 84% | -2.88 | -3.05 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.669 | 0.0231 | 0.0164 | -0.66 | 12% | 59% | 79% | -1.63 | -1.65 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.667 | 0.0234 | 0.0169 | -0.64 | 9% | 59% | 78% | -1.60 | -1.65 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.664 | 0.0238 | 0.0176 | -0.62 | 6% | 63% | 77% | -1.61 | -1.65 |
| UI replacement 15% always | 0.673 | 0.0227 | 0.0166 | -0.61 | 9% | 69% | 78% | -1.64 | -1.60 |
| no assets | 0.671 | 0.0228 | 0.0159 | -0.69 | 6% | 60% | 74% | -1.82 | -1.80 |
| risk aversion 3 | 0.673 | 0.0387 | 0.0274 | -1.13 | 9% | 51% | 63% | -1.94 | -2.20 |
| longer recessions (persistence 0.95) | 0.664 | 0.0237 | 0.0166 | -0.71 | 6% | 60% | 74% | -2.62 | -2.59 |
