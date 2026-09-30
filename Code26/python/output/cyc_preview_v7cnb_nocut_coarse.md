# Cyclicality of quits and N->E at fixed parameters (`output/final_calib_v7cnb_full.log (evaluation 65, objective 0.1819)`, fixed ['kpr=1', 'gamma=2.0', 'ui_rec_mult=0.5', 'beta=0.993', 'r_a=0.00327', 'a_max=60', 'nA=25', 'nAc=100'], 27 types)

sd log = |log(rec/exp)| sqrt(pi_exp pi_rec), the convention of the UE target. Data (early window, trend-adjusted): sd log quit 0.0262 (ratio 0.928), sd log N->E 0.0045 (ratio 0.987), sd log UE 0.0686, employment drop -1.70.

| variant | quit exp | quit rec | ratio | sd log quit | N->E exp | N->E rec | ratio | sd log N->E | sd log UE | dE (pts) | E/pop | NiLF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| data | 0.0226 | 0.0210 | 0.928 | 0.0262 | 0.0530 | 0.0524 | 0.987 | 0.0045 | 0.0686 | -1.70 | 0.620 | 0.220 |
| phi_rec 1.00; phi_rec_H 1.00 | 0.0232 | 0.0184 | 0.797 | 0.0795 | 0.0673 | 0.0545 | 0.810 | 0.0738 | 0.0578 | +0.04 | 0.683 | 0.211 |
| job-finding ratio 0.65; phi_rec 1.00; phi_rec_H 1.00 | 0.0229 | 0.0169 | 0.737 | 0.1066 | 0.0670 | 0.0480 | 0.716 | 0.1169 | 0.1250 | -0.98 | 0.681 | 0.209 |
| non-search arrival x1.40 in recessions; job-finding ratio 0.65; phi_rec 1.00; phi_rec_H 1.00 | 0.0232 | 0.0202 | 0.871 | 0.0484 | 0.0669 | 0.0579 | 0.866 | 0.0505 | 0.1073 | -0.57 | 0.681 | 0.211 |
| non-search arrival x1.60 in recessions; job-finding ratio 0.50; phi_rec 1.00; phi_rec_H 1.00 | 0.0233 | 0.0215 | 0.922 | 0.0286 | 0.0669 | 0.0627 | 0.937 | 0.0229 | 0.1696 | -1.02 | 0.679 | 0.211 |
