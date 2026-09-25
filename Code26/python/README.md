# Python translation of the Code26 MATLAB model

This package is a line-by-line translation of the "simpler model" MATLAB code in
`Code26/` (solver `SimplerMod_May17_splines.m`, simulator `SimplerMod_May17_sim.m`
and the helper functions).  It has two modes:

* **faithful** (`Options.faithful()`, the default): reproduces the MATLAB
  computation including every quirk documented in `../../DEPARTURES.md`.
  Verified: the eight saved solutions in `Code26/Solution/*` are reproduced to
  machine precision (max |difference| 2e-14 in the value functions, 0 in the
  policies), and the stored `Output/Baseline/SimulStats.xls` and the career
  shares in `Code26/CrossCohort.xlsx` are reproduced to 1e-6 (with
  `NiLFdef = 0.3`, see D-R6).
* **corrected** (`Options.corrected()`): every flag switched to the behaviour
  written in the paper/slides.  Individual flags can be toggled with
  `Options.faithful().with_(flag=False)`.

## Layout

| file | MATLAB source |
|---|---|
| `keam/params.py` | parameter block (lines 28-113) and `Options` flags |
| `keam/functions.py` | `utilC.m utilL.m du_dC.m du_dL.m wage.m hprod.m exp2.m de2_dh.m` |
| `keam/solve.py` | value function iteration (lines 133-393) |
| `keam/simulate.py` | panel simulation (lines 48-554), `MarkovStationary.m` |
| `keam/stats.py` | cross-section, careers and cycle statistics (lines 559-1118) |
| `keam/moments.py` | compact moment set used to compare scenarios |
| `keam/matlab_rng.py` | exact replication of MATLAB `rng(seed)` streams |
| `keam/matio.py` | loaders for the saved `.mat` solutions |
| `scripts/verify_solution.py` | compare the faithful solver with `Solution/<name>` |
| `scripts/verify_simulation.py` | compare the faithful simulator with the stored spreadsheets |
| `scripts/departure_impacts.py` | switch each departure off one at a time |
| `scripts/run_experiments.py` | baseline, comparative statics and cohorts |

## Usage

```bash
cd Code26/python
pip install numpy scipy pandas openpyxl xlrd
python scripts/verify_solution.py Baseline      # ~2 s
python scripts/verify_simulation.py 0.3         # all stored experiments
python scripts/run_experiments.py               # faithful
python scripts/run_experiments.py --corrected   # textbook behaviour
python scripts/departure_impacts.py             # one flag at a time
```

```python
from keam import Params, Options, solve, simulate, stats
p = Params()                          # baseline; Params(rtoexpscale=1.2) etc. for experiments
sol = solve(p, Options.faithful())    # arrays shaped (nI, nT, nE, nY, nZ) like MATLAB
r = simulate(p, sol, Options.faithful())
print(stats.summary_table(r))
```

Array conventions follow MATLAB: `gH[i, it, ie, iy, iz]`, ages `it = 0,1,2`
working and `3` retired, husband state `iy = 0` employed, `1` recently
unemployed, `2` unemployed, aggregate state `iz = 0` expansion, `1` recession.
