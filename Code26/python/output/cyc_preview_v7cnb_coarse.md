# Cyclicality of quits and N->E at fixed parameters (`output/final_calib_v7cnb_full.log (evaluation 65, objective 0.1819)`, fixed ['kpr=1', 'gamma=2.0', 'ui_rec_mult=0.5', 'beta=0.993', 'r_a=0.00327', 'a_max=60', 'nA=25', 'nAc=100'], 27 types)

sd log = |log(rec/exp)| sqrt(pi_exp pi_rec), the convention of the UE target. Data (early window, trend-adjusted): sd log quit 0.0262 (ratio 0.928), sd log N->E 0.0045 (ratio 0.987), sd log UE 0.0686, employment drop -1.70.

| variant | quit exp | quit rec | ratio | sd log quit | N->E exp | N->E rec | ratio | sd log N->E | sd log UE | dE (pts) | E/pop | NiLF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| data | 0.0226 | 0.0210 | 0.928 | 0.0262 | 0.0530 | 0.0524 | 0.987 | 0.0045 | 0.0686 | -1.70 | 0.620 | 0.220 |
| calibrated parameters | 0.0248 | 0.0199 | 0.799 | 0.0783 | 0.0680 | 0.0466 | 0.685 | 0.1323 | 0.0874 | -1.57 | 0.668 | 0.235 |
| non-search arrival x1.5 in recessions | 0.0251 | 0.0241 | 0.959 | 0.0146 | 0.0678 | 0.0606 | 0.894 | 0.0394 | 0.0525 | -0.83 | 0.669 | 0.235 |
| non-search arrival x2.0 in recessions | 0.0253 | 0.0284 | 1.121 | 0.0399 | 0.0682 | 0.0741 | 1.085 | 0.0287 | 0.0316 | -0.24 | 0.670 | 0.237 |
| no recession wage cut for the wife (phi_rec 1) | 0.0232 | 0.0172 | 0.742 | 0.1043 | 0.0675 | 0.0510 | 0.756 | 0.0977 | 0.0588 | +0.09 | 0.684 | 0.212 |
| job-finding fall 10% (ratio 0.90) | 0.0250 | 0.0210 | 0.843 | 0.0597 | 0.0681 | 0.0517 | 0.759 | 0.0963 | 0.0381 | -0.85 | 0.670 | 0.235 |
