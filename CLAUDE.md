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
2. `Code26/python/keam/final/` is the final model (monthly, assets, explicit types, transitory
   cost shock, child-care home-productivity multiplier at 25-39, hours-scaled fixed cost).
   Specification: `FINAL_MODEL.md`. Author decisions: monthly period, 30% replacement rate for
   the husband's unemployment income, compensated wage-gap experiment, assets.
3. DONE (cloud, 4 cores, 2026-09-26): 1940s-cohort calibration on the 100-type grid. The adopted
   calibration is the least-squares polish `output/final_calib_ls_full.json/.md` (objective 0.120;
   `scripts/calibrate_ls.py`, polished from the Nelder-Mead point `output/final_calib_full.json`,
   objective 0.149). Quit rates, hours and the recession employment drop are on target; the
   never-working share is +22%, part-time and career shares about -10%. All results use the `ls`
   tag: `output/final_results_ls.*`, `extra_experiments_ls.*`, `cohorts_refined_ls.*`,
   `robustness_final_ls.*`, `figures_ls/`, `irf_careers_ls.json`; `RESULTS.md` is assembled with
   `python scripts/write_results.py --calib output/final_calib_ls_full.json --calib-prev
   output/final_calib_full.json --results output/final_results_ls.json --extra
   output/extra_experiments_ls.json --cohorts2 output/cohorts_refined_ls.json --robust
   output/robustness_final_ls.json --figdir output/figures_ls`. The earlier (unpolished) outputs
   without the tag are kept for reference. Calibration history: `output/final_calib_childcare{2,3,4}.md`.
4. EXPLORED, NOT ADOPTED (2026-09-26): the identification analysis (`output/jacobian_final.md`,
   `output/diag_careers.md`, RESULTS.md 1b) shows the never-working share and the employment rate
   cannot be lowered together. A persistent cost-of-work shock (`rho_kT` in `keam/final/solve.py`,
   nests the iid model at 0, verified to 1e-9 against the previous solver) was calibrated on the
   coarse grid (`output/final_calib_rho_coarse.json`, 0.125 on 27 types) but re-evaluates at 0.238
   on 100 types (RESULTS.md 1a). Candidates left: a finer wage-type grid (n_omega 7-9) or a permanent
   home-productivity type; a calibrated recession job-finding ratio (`lam_f_ratio`, bounds in
   `calibrate.py`) fixes the recession quit rate but drives the job-finding fall toward zero.
5. IN PROGRESS (2026-09-26, author request): (a) recession UI cut `ui_rec_mult=0.5`
   (`output/final_calib_ui_full.*`, objective 0.144, never-working share +25%: does not help); (b) 7
   wage-type points (`output/final_calib_om7_full.*`, 140 types, objective 0.135, never-working +20%:
   does not help); (c) channel decomposition `scripts/channels.py` (`output/channels_ls.md`,
   RESULTS.md 4b): precaution rises with the cyclicality of the husband's risk (job loss, finding,
   UI in recessions), hoarding with the size of the wife's job-finding fall and recession
   persistence; the recession employment drop target is generated mainly by the job-finding fall.
   (d) Version 3 under calibration: UI cut plus a 10% job-finding fall (`lam_f_ratio` 0.90 via a
   degenerate bound, own job loss in recessions free), `output/final_calib_v3_full.*`, then
   `channels.py --variants quick` for the ui, om7 and v3 versions (`output/channels_<tag>.*`).
   Author priority: a version in which precautionary labor supply has a role at least comparable
   to job hoarding; then identify the parameters/targets that govern the split.
6. Known limitations / next steps (see RESULTS.md section 7):
   - DONE: `scripts/cohorts_refined.py` solves the cost scale and tau_w jointly per cohort
     (`output/cohorts_refined_full.md`, RESULTS.md section 3b); the raw-ratio version in
     `run_final.py` is superseded for the cohort narrative;
   - the "cost of work" experiment on the permanent cost alone needs a 90% cut; the
     supplementary child-care and all-costs versions are the relevant ones;
   - NiLF share still 22% too high (see state item 4 for the diagnosis and candidates);
   - DONE: bisection in `size_to_employment` now 12 steps, tolerance 0.2 pp;
   - the routine/trigger channel does not deliver into a CLI remote-control session.

## Plan (execute autonomously, commit and push after each step)

Work on branch `claude/hopeful-ride-vbnou4`. Use all cores (`KEAM_NJOBS` = number of cores).

1. `cd Code26/python && pip install numpy scipy pandas openpyxl xlrd`.
2. (done) Refined cohort accounting: `scripts/cohorts_refined.py`.
3. (done) Calibration polish and full results pipeline for it (`scripts/run_pipeline.sh <calib> <tag>`).
4. (done) `scripts/figures.py` -> `Code26/python/output/figures/` (quit probability, search, hours by
   state over experience averaged over types; refined cohort trend; mechanism decomposition).
   `scripts/irf_careers.py` -> fig6/fig7: impulse responses by career type on the NBER dates.
5. Keep `RESULTS.md` and the PR description current.

## Conventions

- Never modify the MATLAB files or the saved `.mat` solutions.
- Keep `Options.faithful()` reproducing the MATLAB output; verify with
  `python scripts/verify_solution.py Baseline` after touching `keam/solve.py`.
- Every number reported must come from a script in `Code26/python/scripts`; write the
  script, run it, cite its output file.
- Commit messages: plain description of the change. Push to the branch after each step.
- The cloud container is reclaimed a few minutes after the session goes idle, which kills
  background jobs (three runs were lost this way on 2026-09-26). While a long job runs, keep
  the session busy with foreground waits (`timeout 590 bash -c 'until <done>; do sleep 20; done'`,
  repeated), and make long scripts resumable (`calibrate_ls.py --resume` replays its log).
