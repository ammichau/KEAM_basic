# Workstation job (2026-09-28, revised 14:30 UTC): quit-targeted re-estimation under the macro calibration

You are the executor on the author's workstation for the KEAM project (repository KEAM_basic, branch
`claude/hopeful-ride-vbnou4`). Read `CLAUDE.md` (state items 0, 3, 8 and 9) first. Work autonomously, never stop
to ask questions; commit and push after every completed step with a plain descriptive message; every reported
number must come from a script in `Code26/python/scripts` that writes an output file; never modify the MATLAB
files. Setup: `git pull origin claude/hopeful-ride-vbnou4; cd Code26/python; export KEAM_NJOBS=$(nproc)`.
Run long jobs with `nohup ... &` and poll their logs. The tags below (v4em, v7cm, v4nbm, v9nm, v4eq, v7cq, v4nbq,
v9nq and their suffixes) are yours; the cloud session will not write them.

## Step 0: stop the earlier version of this job

If calibrations with the tags v4eq, v7cq, v4nbq or v9nq (the 14:00 UTC version of this job, without the macro
calibration) are still running, kill them (`pkill -f "calibrate_ls.py"`), pull, and start below. Keep any
finished q-tag outputs; they are the beta = 0.99 comparison of step 4.

## Two changes at once, both by the author

1. **Macro calibration** (`beta`, `r_a`): discount factor 0.96 per year and a 4% annual real return on assets,
   i.e. `--fixed beta=0.99661 --fixed r_a=0.00327` (the original model had beta 0.99 per month = 0.886 per year
   and no return; `r_a` is new in `keam/final/params.py`, entering the budget constraints in the solver, the
   retirement problem and the simulator). At the old parameters this roughly doubles attachment (employment 0.76,
   quits halved) and raises assets to about 9 months of income, so the calibration must move a long way. The asset
   grid must be widened: `--fixed a_max=30 --fixed nA=25` (at a_max 15, 60% of household-months sat at the top of the
   grid, `output/asset_distribution_v7c_macro_coarse.md`); this makes each evaluation about 1.6 times slower.
2. **Quit and layoff targets** from `data/QuitLayoff2024.csv` (author; CPS monthly flows 1978-02 to 2016-12,
   seasonally adjusted, unsmoothed, percent): `eqmw_seats` = married women's monthly quit rate from employment to
   non-employment, `elmw_seats` = their layoff rate. Moments: `python3 scripts/quit_layoff_moments.py --file
   data/QuitLayoff2024.csv --quit eqmw_seats --layoff elmw_seats --early-end 1985 --ma 1` (already run:
   `output/quit_layoff.md/.json`; no smoothing anywhere, in the data or in the model). Targets for the 1940s cohort
   from the early window 1978-85 (22 NBER recession months: 1980, 1981-82), trend-adjusted (rate on a linear trend
   and a recession dummy):

| moment | value | use |
|---|---|---|
| quit/m exp | 0.0226 | target (weight 2) |
| quit/m rec | 0.0210 | target (weight 2); recession fall 7% |
| lam_u0 (exogenous job loss, expansion) | 0.0130 | fixed at the layoff rate, not calibrated |
| lam_u1 (exogenous job loss, recession) | 0.0154 | fixed at the layoff rate, not calibrated |

The `E->nonE/m exp` and `E->nonE/m rec` targets are dropped (quits plus layoffs); the other nine targets are
unchanged, so each objective has 11 targets. `lam_f_ratio` stays fixed at 0.80 (the women's UE-rate cyclicality),
as in all carried versions. Tension to watch: the data's recession fall in quits is 7% (the old target had 18%;
over 1978-2016 quits are 4% HIGHER in recessions after the trend) while the 20% job-finding fall alone gave the
model quit drops of 24-34%. Mean liquid assets are NOT targeted; the data figure (about 6 months of household
income on average) is a check: report the model's `mean assets/monthly HH inc` and the distribution.

Common flags (shell variable):

```
M='--fixed beta=0.99661 --fixed r_a=0.00327 --fixed a_max=30 --fixed nA=25'
Q='--set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1 --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec" --extra-target "quit/m exp=0.0226:2.0" --extra-target "quit/m rec=0.0210:2.0" --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --max-nfev 10 --diff-step 0.04'
```

(`--drop` keeps a parameter at its `--set` value outside the calibrated vector; `--drop-target` removes a target;
`--extra-target` on an existing name overrides its value and weight. The calibration JSON records `targets`,
`weights`, `dropped`, `dropped_targets` and `fixed`; `calib_table.py`, `versions_table.py` and
`write_paper_results.py` score each file on its own target set. Every later script (`channels.py`,
`run_pipeline.sh`, `cohorts_refined.py`, `jacobian_channels.py`, `asset_distribution.py`) reads the fixed fields
from the calibration file, so beta, r_a, a_max and nA carry through automatically.)

## Steps

1. Four calibrations, two at a time (each about 2-3 hours on 16 cores with the wider grid). Use
   `eval "python3 -u scripts/calibrate_ls.py ... $M $Q ..."` so that the quoted target names in `$Q` are parsed:
   - v4em (separable CRRA gamma 2, wage cut): `--x0 output/final_calib_v4e_full.json --fixed ui_rec_mult=0.5 $M $Q --tag v4em_full`
   - v7cm (KPR, additive shock, wage cut): `--x0 output/final_calib_v7c_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 $M $Q --tag v7cm_full`
   - v4nbm (separable, no wage cut): `--x0 output/final_calib_v4nb_full.json --fixed ui_rec_mult=0.5 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q --tag v4nbm_full`
   - v9nm (KPR proportional shock, no wage cut): `--x0 output/final_calib_v9n_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 --fixed kT_mult=1 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q --tag v9nm_full`
   The starting points are far from the new optimum (objective about 2 at the old parameters), so expect the
   evaluation limit to bind: when a run stops with the objective still falling, run a second polish from its result
   (`--x0 output/final_calib_<tag>_full.json`, same flags, `--tag <tag>b_full`) and carry the better one (a third
   polish if the second still improves by more than 10%). After each: `python3 scripts/calib_table.py
   output/final_calib_<tag>_full.json > output/final_calib_<tag>_full.md` and `python3 -u scripts/asset_distribution.py
   --calib output/final_calib_<tag>_full.json --tag <tag>` (if the share near a_max exceeds 2%, rerun the
   calibration with `--fixed a_max=45`); commit (`git add -f` the `.log` and `.out` files) and push.
2. Channel decomposition for each: `python3 -u scripts/channels.py --calib output/final_calib_<tag>_full.json --tag <tag>`
   (`--variants quick` for the two worst fits if time is short); commit and push.
3. Full pipeline for the best-fitting separable and the best-fitting KPR version:
   `bash scripts/run_pipeline.sh output/final_calib_<tag>_full.json <tag>`, then
   `python3 -u scripts/cohorts_refined.py --calib output/final_calib_<tag>_full.json --e-mode relative --tag <tag>_rel`
   and `python3 -u scripts/jacobian_channels.py --tag <tag>`; commit and push after each script.
4. Comparisons, for the best-fitting version (tag `<tag>`), in this order:
   - the same targets at the old discount factor (no macro calibration): the q-tag run from the earlier job if it
     finished, otherwise `... $Q` without `$M`, `--tag <tag without m>q_full` (e.g. v7cq); then `channels.py`;
   - job-finding fall free: `$M $Q` without `--set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001`, plus
     `--bound lam_f_ratio:0.5:1.0`, `--tag <tag>f_full`; then `channels.py`;
   - acyclical quits (the full-sample estimate, ratio 1.04): `--extra-target "quit/m rec=0.0234:2.0"` in place of
     0.0210, `--tag <tag>a_full`; then `channels.py`.
   Commit and push after each.
5. `python3 scripts/write_paper_results.py --specs "CRRA macro (v4em)=v4em,KPR macro (v7cm)=v7cm,CRRA no wage cut macro (v4nbm)=v4nbm,KPR no wage cut macro (v9nm)=v9nm" --out ../../PAPER_RESULTS_quit.md`
   (add the step-4 tags as further columns as they exist); write the outcome (objectives, the 11 targets' fit, mean
   assets against the 6-month check, precaution/hoarding shares, the comparisons) into `CLAUDE.md` state item 9;
   commit and push.
