# Workstation job (2026-09-29, 17:00 UTC): separate offer arrival from non-participation; re-estimation at beta 0.993

You are the executor on the author's workstation for the KEAM project (repository KEAM_basic, branch
`claude/hopeful-ride-vbnou4`). Read `CLAUDE.md` (state items 8, 9 and 10) first. Work autonomously, never stop to ask
questions; commit and push after every completed step with a plain descriptive message; every reported number must
come from a script in `Code26/python/scripts` that writes an output file; never modify the MATLAB files. Setup:
`git pull origin claude/hopeful-ride-vbnou4; cd Code26/python`. Run long jobs with `nohup ... &` and poll their logs.
The tags below (v7cn, v4en, v9nn, v4nbn and their suffixes) are yours; the cloud session will not write them.

## Step 0: stop the 15:00 UTC job

The calibrations v7ck and v4ek started at about 15:00 UTC (and anything queued after them: v9nk, v4nbk, channels) lack
the model change below: kill them (`pkill -f "calibrate_ls.py"`; kill the chain script too), delete nothing, pull.

## The model change (author, 2026-09-29)

Non-participants receive job offers without searching: job finding = `lam_n(Z) + lam_f(Z) s^nu` (new FinalParams field
`lam_n`, expansion / recession; 0 reproduces the old model exactly). `lam_n` is inferred from married women's N->E
flow in the CPS (`output/ne_nu_moments.md`, early window 1978-85, trend-adjusted: 0.0530 per month in expansions, 0.0524
in recessions, i.e. nearly acyclical), through two calibrated parameters `lam_n0` (level) and `lam_n_ratio`
(recession / expansion) and two new targets `N->E/m exp` and `N->E/m rec` (the model's monthly entry rate into
employment of women non-employed with search below `s_bar`; `moments.py`). The unemployed's job finding keeps its 20%
recession fall (`lam_f_ratio` 0.80, the women's UE-rate cyclicality). Counterfactual "acyclical job finding"
(`experiments.acyclical_finding`, used by `channels.py`, `jacobian_channels.py`, the pipeline) now holds both rates at
their expansion values. 13 targets.

## Settings

Discount factor 0.993 per month with the 4% annual return, fine savings choice grid, as in the 15:00 UTC job:

```
M='--fixed beta=0.993 --fixed r_a=0.00327 --fixed a_max=60 --fixed nA=25 --fixed nAc=100'
Q='--set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1 --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec" --extra-target "quit/m exp=0.0226:2.0" --extra-target "quit/m rec=0.0210:2.0" --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --diff-step 0.04'
N='--extra lam_n0,lam_n_ratio --set lam_n0=0.07 --set lam_n_ratio=1.0 --extra-target "N->E/m exp=0.0530:1.0" --extra-target "N->E/m rec=0.0524:1.0"'
```

Use `eval "python3 -u scripts/calibrate_ls.py ... $M $Q $N --max-nfev 5 --tag <tag>_full"` (the quoted target names need
`eval`). `--extra` must list `lam_n0,lam_n_ratio` together with the other extras of the starting file: check the `x`
dict of each `--x0` file and append any optional parameter it carries (`lam_f_ratio` is already handled by `--set`
and `--bound`). Cost: about 770 s per evaluation with all cores, about 15 hours per calibration at `--max-nfev 5`;
run two at a time with `export KEAM_NJOBS=$(( $(nproc) / 2 ))` each. Mean assets are a check (about 15 months of
household income expected), not a target; report them from `asset_distribution.py`.

## Steps

1. Calibrations, two at a time, priority order. First pair:
   - v7cn (KPR, additive shock, wage cut): `--x0 output/final_calib_v7cmb_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 $M $Q $N --max-nfev 5 --tag v7cn_full`
   - v4en (separable CRRA gamma 2, wage cut): `--x0 output/final_calib_v4emb_full.json --fixed ui_rec_mult=0.5 $M $Q $N --max-nfev 5 --tag v4en_full`
   Second pair:
   - v9nn (KPR proportional shock, no wage cut): `--x0 output/final_calib_v9nmb_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 --fixed kT_mult=1 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q $N --max-nfev 5 --tag v9nn_full`
   - v4nbn (separable, no wage cut): `--x0 output/final_calib_v4nbmc_full.json --fixed ui_rec_mult=0.5 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q $N --max-nfev 5 --tag v4nbn_full`
   After each: `python3 scripts/calib_table.py output/final_calib_<tag>_full.json > output/final_calib_<tag>_full.md`,
   `python3 -u scripts/asset_distribution.py --calib output/final_calib_<tag>_full.json --tag <tag>`; commit (`git add -f`
   the `.log` and `.out` files) and push. If a run stops at the evaluation limit with the objective still falling by
   more than 10% over its last iteration, a second polish (`--x0 output/final_calib_<tag>_full.json`, same flags,
   `--max-nfev 4`, `--tag <tag>b_full`); carry the better one.
2. Channel decomposition for each finished calibration as cores free up:
   `python3 -u scripts/channels.py --calib output/final_calib_<tag>_full.json --tag <tag>`; commit and push.
3. After all four: full pipeline for the best KPR and the best separable version, one at a time with all cores
   (about 20 hours each): `bash scripts/run_pipeline.sh output/final_calib_<tag>_full.json <tag>`, then
   `python3 -u scripts/cohorts_refined.py --calib output/final_calib_<tag>_full.json --e-mode relative --tag <tag>_rel`
   and `python3 -u scripts/jacobian_channels.py --tag <tag>`; commit and push after each script.
4. `python3 scripts/write_paper_results.py --specs "KPR wage cut (v7cn)=v7cn,CRRA wage cut (v4en)=v4en,KPR no wage cut (v9nn)=v9nn,CRRA no wage cut (v4nbn)=v4nbn" --out ../../PAPER_RESULTS_quit.md`
   with the columns that exist; write the outcome (objectives, the 13 targets' fit, mean assets, precaution / hoarding
   shares) into `CLAUDE.md` state item 10; commit and push.
