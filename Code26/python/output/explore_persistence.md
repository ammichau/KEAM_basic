# Persistence of the cost-of-work shock at the calibrated parameters

Calibrated point `output/final_calib_full.json`; the first row is the 100-type iid check (calibrated moments: E/pop 0.6644, share NiLF 0.2598, quit/m exp 0.0311, quit/m rec 0.0237).

| variant | obj | E/pop | hours|E | share Lifecycle | share PT | share Career | share NiLF | quit/m exp | quit/m rec | E->nonE/m exp | E->nonE/m rec | dE/pop rec-exp (pts) | wage gap (hourly ratio) | U rate | n careers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| full grid, iid 5 nodes (check) | 0.150 | 0.6644 | 0.4141 | 0.3038 | 0.2531 | 0.1833 | 0.2598 | 0.0312 | 0.0237 | 0.0485 | 0.0432 | -1.5936 | 0.7478 | 0.0459 | 4800.0000 |
| coarse, iid 5 nodes | 0.151 | 0.6593 | 0.4167 | 0.3200 | 0.2292 | 0.1744 | 0.2765 | 0.0325 | 0.0254 | 0.0498 | 0.0446 | -1.7010 | 0.7379 | 0.0443 | 4800.0000 |
| coarse, 3 nodes, rho=0.0 | 2.621 | 0.6932 | 0.4240 | 0.2181 | 0.4198 | 0.1985 | 0.1635 | 0.0198 | 0.0153 | 0.0372 | 0.0347 | -0.5962 | 0.7386 | 0.0419 | 4800.0000 |
| coarse, 3 nodes, rho=0.5 | 1.643 | 0.5959 | 0.4118 | 0.4090 | 0.1096 | 0.0931 | 0.3883 | 0.0424 | 0.0340 | 0.0597 | 0.0533 | -1.6528 | 0.7302 | 0.0816 | 4800.0000 |
| coarse, 3 nodes, rho=0.8 | 3.945 | 0.5587 | 0.3828 | 0.1252 | 0.5142 | 0.0017 | 0.3590 | 0.0448 | 0.0344 | 0.0621 | 0.0541 | -0.7447 | 0.6737 | 0.1360 | 4800.0000 |
| coarse, 3 nodes, rho=0.9 | 4.696 | 0.5434 | 0.3744 | 0.0327 | 0.6860 | 0.0004 | 0.2808 | 0.0338 | 0.0290 | 0.0511 | 0.0485 | -1.2087 | 0.6351 | 0.1144 | 4800.0000 |
| coarse, 3 nodes, rho=0.95 | 6.094 | 0.5563 | 0.3767 | 0.0473 | 0.6979 | 0.0060 | 0.2487 | 0.0177 | 0.0176 | 0.0350 | 0.0375 | -2.6478 | 0.6383 | 0.0840 | 4800.0000 |
