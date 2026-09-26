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
3. The 1940s-cohort calibration is NOT finished. Findings so far (27-type grid):
   - A first Nelder-Mead run (`output/final_calib_coarse.*`) stalled at objective 3.6 with the
     life-cycle career share near 2% (target 31%) because the young-age cost multiplier
     `km_max` was bounded at 6.
   - `scripts/explore_lifecycle.py` (`output/explore_lifecycle.json`) shows the life-cycle
     share rising with `km_max`: 2% at 3, 6% at 4.5, 9% at 6, 12% at 8 (with `kbar_max`
     0.06); objective 2.27 at (8, 0.06). Bounds were widened to `km_max` <= 15.
   - The unemployment rate and the wife's income share are model outputs on slides p.36-37,
     not calibration targets; they were removed from `TARGETS`. `s_bar` is fixed at 0.25.
     The wage-gap moment is the hourly wage ratio with the husband at 2,000 hours/year.
   - A second coarse run from (`km_max` 8, `kbar_max` 0.06) was started in the cloud
     (`output/final_calib_coarse2.*`); if it is present, start from its best point.
   Remaining tension: career and NiLF shares too high, part-time too low.
   `scripts/explore_parttime.py` (`output/explore_parttime.json`, cloud): raising `z_h` from
   0.45 to 0.65 collapses employment (0.59 -> 0.36); lowering `nu_h` from 0.65 to 0.45 raises
   part-time 16% -> 20% and lowers NiLF 37% -> 24% at the cost of employment (0.59 -> 0.67,
   too high) and the objective. Recommendation: add `nu_h` (bounds 0.3-0.8) and `z_h`
   (0.3-0.6) to the calibrated parameters (OPTIONAL_PARAMS) together with `delta_e`.
   Cloud partial round 2 (`output/final_calib_coarse2_partial.json`, objective 1.86, km_max 8).
   - NEW INGREDIENT (cloud): a child-care multiplier on home productivity at ages 25-39
     (`home_young_mult` in FinalParams, default 1 = off). `scripts/explore_childcare.py` shows it
     generates life-cycle women far more economically than the utility-cost multiplier alone.
     `scripts/calibrate_childcare.py` calibrates it jointly with nu_h, z_h (and `--extra alpha_h`).
     Round 2 (`output/final_calib_childcare2.json/.md`, 27 types): objective 0.21, all targets
     within 10% except NiLF 29% (target 22%) and the wage gap 0.82 (target 0.71). Round 3 with
     alpha_h calibrated: objective 0.19 (`childcare3`). Round 4 added e_max and kappa_h_power (hours-
     scaled fixed cost, `FinalParams.kappa_h_power`) with tau_w <= 0.78: objective 0.14
     (`output/final_calib_childcare4.md`), wage gap on target, NiLF +26% the only miss > 13%.
     A 100-type polish from round 4 is running in the cloud (`output/final_calib_full.*`);
     then `run_final.py --full`, `robustness_final.py`, `write_results.py` -> RESULTS.md.
   - The SMM initial simplex was fixed (parameters starting at a bound midpoint were frozen).
   - The routine/trigger channel does NOT deliver into a CLI remote-control session; the
     workstation session must be given its job by the author in the Claude Code app.

## Plan (execute autonomously, commit and push after each step)

Work on branch `claude/hopeful-ride-vbnou4`. Use all cores (`KEAM_NJOBS` = number of cores).

1. `cd Code26/python && pip install numpy scipy pandas openpyxl xlrd`.
2. Run `python scripts/explore_lifecycle.py output/final_calib_coarse.json` and read where the
   life-cycle share becomes positive. If it never does within the bounds, widen `BOUNDS`
   (`km_max` up to 10, `kbar_max` up to 1.0) in `keam/final/calibrate.py` and check whether
   the experience process is the obstacle (a woman out for 15 years loses 55% of e; consider
   whether `delta_e` or `e_max` should be part of the calibration).
3. Calibrate on the full 100-type grid, starting from the best coarse child-care point:
   copy `scripts/calibrate_childcare.py`, change `base = FinalParams()` (100 types), run it with
   `--x0 output/final_calib_childcare3.json --extra alpha_h --maxfev 300 --tag full` (or from
   `childcare2.json` if round 3 is absent), on all cores in the background. Then
   `python scripts/calib_table.py output/final_calib_full.json > output/final_calib_full.md`.
   Report the fit of every target; if one cannot be reached, say which and why.
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
