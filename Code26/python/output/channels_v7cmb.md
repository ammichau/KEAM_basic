# Precautionary labor supply versus job hoarding (`output/final_calib_v7cmb_full.json`, parameters held fixed across variants)

Quit gap = recession minus expansion monthly quit rate, percentage points. Precaution share = fall in the gap when the husband's risk (job loss, job finding, recession UI cut) is made acyclical; hoarding share = fall when the wife's job-finding efficiency is made acyclical; both off = fall when both are; the remainder is the wage cut, the wife's own cyclical job loss and interactions.

| variant | E/pop | quit exp | quit rec | gap | precaution | hoarding | both off | dE base | dE acyc. husband |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 0.663 | 0.0235 | 0.0173 | -0.63 | 6% | 62% | 77% | -1.62 | -1.68 |
| job finding falls 5% in recessions (ratio 0.95) | 0.666 | 0.0239 | 0.0206 | -0.33 | 22% | 29% | 57% | -0.02 | -0.06 |
| job finding falls 30% in recessions (ratio 0.70) | 0.660 | 0.0233 | 0.0144 | -0.88 | 3% | 73% | 84% | -2.90 | -3.06 |
| husband job loss x2.5 in recessions (data: x1.78) | 0.667 | 0.0229 | 0.0163 | -0.66 | 11% | 58% | 79% | -1.61 | -1.68 |
| husband job finding 0.20 in recessions (data: 0.28) | 0.666 | 0.0232 | 0.0167 | -0.65 | 10% | 59% | 78% | -1.60 | -1.68 |
| UI replacement 15% in recessions (ui_rec_mult 0.5) | 0.663 | 0.0235 | 0.0173 | -0.63 | 6% | 62% | 77% | -1.62 | -1.68 |
| UI replacement 15% always | 0.671 | 0.0225 | 0.0160 | -0.65 | 10% | 69% | 78% | -1.73 | -1.69 |
| no assets | 0.669 | 0.0226 | 0.0157 | -0.69 | 5% | 61% | 73% | -1.85 | -1.80 |
| risk aversion 3 | 0.672 | 0.0386 | 0.0274 | -1.13 | 9% | 51% | 63% | -1.94 | -2.16 |
| longer recessions (persistence 0.95) | 0.662 | 0.0234 | 0.0159 | -0.75 | 10% | 62% | 73% | -2.83 | -2.59 |
