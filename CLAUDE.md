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
3. DONE (cloud, 4 cores, 2026-09-26): 1940s-cohort calibration on the 100-type grid
   (`output/final_calib_full.json/.md`, objective 0.15; all targets within 10% except NiLF
   +18% and the recession quit rate -15%), results (`output/final_results_full.md/.json`),
   supplementary experiments (`output/extra_experiments_full.*`), robustness
   (`output/robustness_final.*`) and the write-up `RESULTS.md` at the repository root.
   Calibration history: `output/final_calib_childcare{2,3,4}.md`; scans `output/explore_*.json`.
4. IN PROGRESS (2026-09-26): the iid-shock calibration cannot lower the employment rate (+7%) and
   the never-working share (+18%) together because its ever-working women return from
   non-employment too fast (`output/jacobian_final.md`, `output/diag_careers.md`). A persistent
   cost-of-work shock (`rho_kT`, `keam/final/solve.py`, nests the iid model at 0; verified to
   1e-9 against the previous solver) is being calibrated: coarse grid first
   (`scripts/calibrate_ls.py --coarse --fixed n_kT=3 --extra rho_kT --tag rho_coarse`), then a
   full-grid polish (`--tag rho_full`, 100 types, about 180 s per evaluation), then the results
   pipeline with `--calib output/final_calib_rho_full.json`. Keep the iid results until the
   persistent version is complete; report both.
5. Known limitations / next steps (see RESULTS.md section 7):
   - DONE: `scripts/cohorts_refined.py` solves the cost scale and tau_w jointly per cohort
     (`output/cohorts_refined_full.md`, RESULTS.md section 3b); the raw-ratio version in
     `run_final.py` is superseded for the cohort narrative;
   - the "cost of work" experiment on the permanent cost alone needs a 90% cut; the
     supplementary child-care and all-costs versions are the relevant ones;
   - NiLF share still 18% too high; candidates: alpha_h upper range, s_bar-independent
     participation definition, or a permanent home-productivity type;
   - bisection resolution in `size_to_employment` (8 steps) is coarse; raise `maxit`;
   - the routine/trigger channel does not deliver into a CLI remote-control session.

## Plan (execute autonomously, commit and push after each step)

Work on branch `claude/hopeful-ride-vbnou4`. Use all cores (`KEAM_NJOBS` = number of cores).

1. `cd Code26/python && pip install numpy scipy pandas openpyxl xlrd`.
2. (done) Refined cohort accounting: `scripts/cohorts_refined.py`.
3. Improve the NiLF fit: persistent cost shock (state item 4). After the full-grid calibration,
   regenerate results with `run_final.py --full`, `extra_experiments.py`, `robustness_final.py`,
   `cohorts_refined.py`, `figures.py`, `irf_careers.py`, `write_results.py` (all take `--calib`).
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
