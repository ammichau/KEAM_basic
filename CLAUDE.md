# KEAM_basic: project context for Claude Code

Read this first. It is the hand-off from the cloud session that built `Code26/python`.

## What this project is

A quantitative paper (target: JEEA) on how the rise in married women's employment changed
U.S. business-cycle dynamics. Mechanism: married women's quits to non-employment fall in
recessions (precautionary labor supply against the husband's job-loss risk, plus job hoarding),
so their employment is less cyclical; forces that raise attachment (returns to experience,
lower cost of work) make it more cyclical, a closing wage gap less so. Reference documents:
`KEAM_Klein.pdf` (slides, April 2025, authoritative) and `EM_DemogBCtrends.pdf` (paper draft,
2023).

## State of the work

1. `Code26/python/keam/` is an exact Python translation of the MATLAB model in `Code26/`
   (reproduces every saved solution to 1e-14). `DEPARTURES.md` catalogues where that code
   departs from the model in the slides and quantifies each item. The MATLAB calibration
   rests on numerical artifacts; do not use it for results.
2. `Code26/python/keam/final/` is the final model (monthly, assets, explicit types,
   transitory cost shock). Specification: `FINAL_MODEL.md`. Decisions already taken by the
   author: monthly period, 30% replacement rate for the husband's unemployment income,
   compensated wage-gap experiment (household income held constant at baseline behaviour,
   income share moved toward the wife), assets with a borrowing constraint.
3. The 1940s-cohort calibration is NOT finished. A single Nelder-Mead run on the 27-type grid
   (`Code26/python/output/final_calib_coarse.*`) stalls at objective 3.6: quit rates, hours,
   employment and the wage gap are close, but the life-cycle career share stays near 2%
   against a 31% target (the young-age cost multiplier `km_max` never gets large enough for
   high-productivity women to stay out at 25-39 and return at 40). Career and NiLF shares
   are too high, the wife's income share too high.

## Plan (execute autonomously, commit and push after each step)

Work on branch `claude/hopeful-ride-vbnou4`. Use all cores (`KEAM_NJOBS` = number of cores).

1. `cd Code26/python && pip install numpy scipy pandas openpyxl xlrd`.
2. Run `python scripts/explore_lifecycle.py output/final_calib_coarse.json` and read where the
   life-cycle share becomes positive. If it never does within the bounds, widen `BOUNDS`
   (`km_max` up to 10, `kbar_max` up to 1.0) in `keam/final/calibrate.py` and check whether
   the experience process is the obstacle (a woman out for 15 years loses 55% of e; consider
   whether `delta_e` or `e_max` should be part of the calibration).
3. Calibrate on the full 100-type grid:
   `python scripts/calibrate_final.py --global 300 --starts 3 --maxfev 400 --tag full`.
   Report the fit of every target in `TARGETS`. If a target cannot be reached, say which and
   why, do not silently drop it.
4. Produce results: `python scripts/run_final.py --calib output/final_calib_full.json --full`.
   This writes baseline moments, the three single-factor experiments sized to the 1970s
   employment rate, the cohort accounting and the mechanism counterfactuals to
   `output/final_results_full.md`.
5. Write `RESULTS.md` at the repository root: calibration table (target vs model), the
   experiment tables, the cohort table, the counterfactual decomposition of the cyclical
   quit change, and a short list of what is fragile. Update the draft PR
   (https://github.com/ammichau/KEAM_basic/pull/1) description.
6. Robustness the referee will ask for: no-assets limit (`a_max` small), asset grid and
   hours grid refinement, alternative unemployment threshold, `phi_rec_H` = 1.

## Conventions

- Never modify the MATLAB files or the saved `.mat` solutions.
- Keep `Options.faithful()` reproducing the MATLAB output; verify with
  `python scripts/verify_solution.py Baseline` after touching `keam/solve.py`.
- Every number reported must come from a script in `Code26/python/scripts`; write the
  script, run it, cite its output file.
- Commit messages: plain description of the change. Push to the branch after each step.
