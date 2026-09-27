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
   Versions so far (RESULTS.md 1a/4b; objective with 12 targets unless noted; precaution / hoarding shares
   of the recession quit drop): adopted `ls` 0.120, 19/39; `ui` (UI cut) 0.144, 29/34; `om7` (7 wage types)
   0.135, 21/38; `v3` (UI cut + 10% job-finding fall) 0.113, 29/25; `v4` (UI cut + job-finding fall
   calibrated to the author's UE-rate cyclicality target sd log UE = 0.0686, 13 targets) 0.256, the
   least-squares run left the fall at 15% (model sd 0.056); `v5` = v4 with log utility (gamma = 1,
   balanced growth; author request) queued after v4 (`output/chain_v5.out`). The UE-cyclicality data
   imply a job-finding fall of about 18-20%, so precaution above hoarding (v3) is not supported by it;
   near parity (ui/v4) is what the data allow.
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

Work on branch `claude/hopeful-ride-vbnou4`. Use all cores (`export KEAM_NJOBS=$(nproc)`). Setup:
`git pull origin claude/hopeful-ride-vbnou4; cd Code26/python; pip install numpy scipy pandas openpyxl xlrd`.
The cloud session (4 cores) stopped its compute on 2026-09-26 22:40 UTC so that the workstation can
run the remaining agenda without file conflicts; the tagged outputs ls, ui, om7, v3, v4 are final.

Goal: settle the 1940s-cohort calibration and measure the split of the recession fall in married
women's quits between precautionary labor supply (husband's cyclical risk) and job hoarding (the
wife's cyclical job finding). Author priority: a data-disciplined version in which precautionary
labor supply has a role at least comparable to hoarding. The women's UE-rate cyclicality target
(sd of the log rate 0.0686; men 0.0765, matched by the husband's rates) is in `calibrate.py` TARGETS.

1. Version 4c: recession UI cut and a 20% job-finding fall.
   `python3 -u scripts/calibrate_ls.py --x0 output/final_calib_v4_full.json --fixed ui_rec_mult=0.5
   --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --max-nfev 6 --diff-step 0.04 --tag v4c_full`
   then `python3 scripts/calib_table.py output/final_calib_v4c_full.json > output/final_calib_v4c_full.md`
   and `python3 -u scripts/channels.py --calib output/final_calib_v4c_full.json --tag v4c`.
2. Version 5: log utility (gamma = 1, balanced growth), same assumptions as 4c, from the rescaled start
   `output/x0_gamma1.json`: coarse Nelder-Mead, full-grid polish, channel decomposition.
   `python3 -u scripts/calibrate_childcare.py --fixed gamma=1 --fixed ui_rec_mult=0.5 --x0 output/x0_gamma1.json
   --extra alpha_h,e_max,kappa_h_power --maxfev 250 --tag v5_coarse`
   `python3 -u scripts/calibrate_ls.py --x0 output/final_calib_v5_coarse.json --fixed gamma=1.0 --fixed ui_rec_mult=0.5
   --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --max-nfev 6 --diff-step 0.04 --tag v5_full`
   `python3 -u scripts/channels.py --calib output/final_calib_v5_full.json --tag v5`
   If the coarse stage ends above objective 0.5, run a second Nelder-Mead round from its result
   (`--x0 output/final_calib_v5_coarse.json --tag v5_coarse2`) before the polish.
3. Carry forward the lowest-objective version among v4c and v5 whose precautionary share is within 10
   points of the hoarding share; prefer v5 if it fits acceptably (objective below 0.3, cyclical moments
   within 15%) because of the balanced-growth argument, and say so. Run
   `bash scripts/run_pipeline.sh output/final_calib_<tag>_full.json <tag>`, then regenerate RESULTS.md:
   `python3 scripts/write_results.py --calib output/final_calib_<tag>_full.json --calib-prev output/final_calib_ls_full.json
   --results output/final_results_<tag>.json --extra output/extra_experiments_<tag>.json --cohorts2 output/cohorts_refined_<tag>.json
   --robust output/robustness_final_<tag>.json --figdir output/figures_<tag>
   --calib-alt "adopted iid=output/final_calib_ls_full.json,version 3=output/final_calib_v3_full.json,version 4=output/final_calib_v4_full.json,version 4c=output/final_calib_v4c_full.json,version 5 (log utility)=output/final_calib_v5_full.json,recession UI cut=output/final_calib_ui_full.json,7 wage types=output/final_calib_om7_full.json"
   --channels "<tag>=output/channels_<tag>.json,adopted=output/channels_ls.json,version 4=output/channels_v4.json,version 4c=output/channels_v4c.json,version 5=output/channels_v5.json,UI cut=output/channels_ui.json"`
   (drop files that do not exist). Update the state section above and the PR description.
5. Version 6 (after 4 and 5): persistent cost-of-work shock under the version-4c assumptions. With the
   UE-rate cyclicality matched (20% job-finding fall) the iid-shock model overstates the recession quit
   drop (34% versus 18% in the data; RESULTS.md 4b/4c), because a quit driven by a one-month cost draw is
   very sensitive to re-entry prospects. A persistent shock weakens that link. Run on the full grid
   (3-node shock, about 3x the solve time):
   `python3 -u scripts/calibrate_ls.py --x0 output/final_calib_v4c_full.json --fixed ui_rec_mult=0.5 --fixed n_kT=3
   --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --extra rho_kT --set rho_kT=0.5 --max-nfev 6 --diff-step 0.04 --tag v6_full`
   then `scripts/channels.py --calib output/final_calib_v6_full.json --tag v6`. Judge it on the recession
   quit rate, the UE cyclicality and the unemployment rate (the earlier persistent version had 12%).
4. Report (RESULTS.md section 6 and the reply): a table of all versions (objective, the 13 targets'
   deviations, precaution and hoarding shares, recession employment drop with and without cyclical
   husband risk); whether log utility preserves the precautionary channel and what it does to the fit;
   which parameters and targets govern the split (section 4b).

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
