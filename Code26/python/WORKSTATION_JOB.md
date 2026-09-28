# Workstation job (2026-09-28): quit-targeted re-estimation of the four carried versions

You are the executor on the author's workstation for the KEAM project (repository KEAM_basic, branch
`claude/hopeful-ride-vbnou4`). Read `CLAUDE.md` (state items 0, 3 and 8) first. Work autonomously, never stop
to ask questions; commit and push after every completed step with a plain descriptive message; every reported
number must come from a script in `Code26/python/scripts` that writes an output file; never modify the MATLAB
files. Setup: `git pull origin claude/hopeful-ride-vbnou4; cd Code26/python; export KEAM_NJOBS=$(nproc)`.
Run long jobs with `nohup ... &` and poll their logs. The tags below (v4eq, v7cq, v4nbq, v9nq and their
suffixes) are yours; the cloud session will not write them.

## Data and targets

`data/QuitLayoff2024.csv` (author; CPS monthly flows 1978-02 to 2016-12, seasonally adjusted, unsmoothed, percent; the `_sa.csv` file is the same data as 13-month centred moving averages, not used for targets):
`eqmw_seats` = married women's monthly quit rate from employment to non-employment, `elmw_seats` = their monthly
layoff rate to non-employment. Moments: `python3 scripts/quit_layoff_moments.py --file data/QuitLayoff2024.csv
--quit eqmw_seats --layoff elmw_seats --early-end 1985 --ma 1` (already run: `output/quit_layoff.md/.json`; no smoothing
anywhere, in the data or in the model).

Targets for the 1940s cohort come from the early window 1978-85 (22 NBER recession months: 1980, 1981-82),
trend-adjusted (rate on a linear trend and a recession dummy):

| moment | value | use |
|---|---|---|
| quit/m exp | 0.0226 | target (weight 2) |
| quit/m rec | 0.0210 | target (weight 2); recession fall 7% |
| lam_u0 (exogenous job loss, expansion) | 0.0130 | fixed at the layoff rate, not calibrated |
| lam_u1 (exogenous job loss, recession) | 0.0154 | fixed at the layoff rate, not calibrated |

The `E->nonE/m exp` and `E->nonE/m rec` targets are dropped (quits plus layoffs); the other nine targets are
unchanged, so each objective has 11 targets. `lam_f_ratio` stays fixed at 0.80 (the women's UE-rate cyclicality),
as in all carried versions. Note the tension to watch: the data's recession fall in quits is 7% (the old target
had 18%; over the full 1978-2016 sample quits are 4% HIGHER in recessions after the trend) while the 20% job-finding
fall alone gave the model quit drops of 24-34%.

Common flags (shell variable):

```
Q='--set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1 --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec" --extra-target "quit/m exp=0.0226:2.0" --extra-target "quit/m rec=0.0210:2.0" --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --max-nfev 8 --diff-step 0.04'
```

(`--drop` keeps a parameter at its `--set` value outside the calibrated vector; `--drop-target` removes a target;
`--extra-target` on an existing name overrides its value and weight. The calibration JSON records `targets`,
`weights`, `dropped` and `dropped_targets`; `calib_table.py`, `versions_table.py` and `write_paper_results.py`
score each file on its own target set.)

## Steps

1. Four calibrations, two at a time (each about 1-2 hours on 16 cores). Use `eval "python3 -u scripts/calibrate_ls.py ... $Q ..."`
   so that the quoted target names in `$Q` are parsed:
   - v4eq (separable CRRA gamma 2, wage cut): `--x0 output/final_calib_v4e_full.json --fixed ui_rec_mult=0.5 $Q --tag v4eq_full`
   - v7cq (KPR, additive shock, wage cut): `--x0 output/final_calib_v7c_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 $Q --tag v7cq_full`
   - v4nbq (separable, no wage cut): `--x0 output/final_calib_v4nb_full.json --fixed ui_rec_mult=0.5 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $Q --tag v4nbq_full`
   - v9nq (KPR proportional shock, no wage cut): `--x0 output/final_calib_v9n_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 --fixed kT_mult=1 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $Q --tag v9nq_full`
   After each: `python3 scripts/calib_table.py output/final_calib_<tag>_full.json > output/final_calib_<tag>_full.md`,
   commit (`git add -f` the `.log` and `.out` files) and push. If a run stops at the evaluation limit with the
   objective still falling, run a second polish from its result (`--x0 output/final_calib_<tag>_full.json`, same
   flags, `--tag <tag>b_full`) and carry the better one.
2. Channel decomposition for each: `python3 -u scripts/channels.py --calib output/final_calib_<tag>_full.json --tag <tag>`
   (`--variants quick` for the two worst fits if time is short); commit and push.
3. Full pipeline for the best-fitting separable and the best-fitting KPR version:
   `bash scripts/run_pipeline.sh output/final_calib_<tag>_full.json <tag>`, then
   `python3 -u scripts/cohorts_refined.py --calib output/final_calib_<tag>_full.json --e-mode relative --tag <tag>_rel`
   and `python3 -u scripts/jacobian_channels.py --tag <tag>`; commit and push after each script.
4. Exploration of the quit-cyclicality tension, for the best-fitting version (tag `<tag>`):
   - job-finding fall free: same flags without `--set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001`, plus
     `--bound lam_f_ratio:0.5:1.0`, `--tag <tag>f_full`;
   - acyclical quits (the full-sample estimate, ratio 1.04): `--extra-target "quit/m rec=0.0234:2.0"` in place of 0.0210, `--tag <tag>a_full`;
   then `channels.py` for both; commit and push.
5. `python3 scripts/write_paper_results.py --specs "CRRA quit-targeted (v4eq)=v4eq,KPR quit-targeted (v7cq)=v7cq,CRRA no wage cut quit-targeted (v4nbq)=v4nbq,KPR no wage cut quit-targeted (v9nq)=v9nq" --out ../../PAPER_RESULTS_quit.md`
   (add the step-4 tags as further columns if they exist); write the outcome (objectives, the 11 targets' fit,
   precaution/hoarding shares, the exploration) into `CLAUDE.md` state item 8; commit and push.
