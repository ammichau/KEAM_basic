# Workstation job (2026-09-29, 15:00 UTC): quit-targeted re-estimation at beta 0.993 with the fine savings choice grid

You are the executor on the author's workstation for the KEAM project (repository KEAM_basic, branch
`claude/hopeful-ride-vbnou4`). Read `CLAUDE.md` (state items 8 and 9) first. Work autonomously, never stop to ask
questions; commit and push after every completed step with a plain descriptive message; every reported number must
come from a script in `Code26/python/scripts` that writes an output file; never modify the MATLAB files. Setup:
`git pull origin claude/hopeful-ride-vbnou4; cd Code26/python`. Run long jobs with `nohup ... &` and poll their logs.
The tags below (v7ck, v4ek, v9nk, v4nbk and their suffixes) are yours; the cloud session will not write them.
Nothing is running: the previous job (macro calibration, tags v4em.., v7cm.., v4nbm.., v9nm.., v7cq, v7cmbf, v7cmba)
is complete and superseded (its results were computed on the lumpy savings grid, CLAUDE.md state item 9).

## Settings (author decision, 2026-09-29)

1. **Discount factor 0.993 per month** (0.919 per year) with the 4% annual real return (`r_a` 0.00327). At the v7cmb
   parameters this gives mean assets of about 15 months of household income (`asset_distribution_v7cmb_c100_beta0.993.md`).
   Mean assets stay a check, not a target.
2. **Fine savings choice grid**: `nAc=100` (a' chosen on a 100-point quadratic grid with the continuation value
   interpolated in a'; the state grid stays `nA=25`), `a_max=60` (at beta 0.993 the v7cmb parameters put 4.8% of
   household-months near 45; if the share near 60 exceeds 2% in a finished calibration, polish it again with a_max 90).
3. Targets as in the previous job: `lam_u` fixed from the layoff data, the quit rate targeted (0.0226 / 0.0210), the
   E->nonE targets dropped (11 targets), `lam_f_ratio` fixed at 0.80.

```
M='--fixed beta=0.993 --fixed r_a=0.00327 --fixed a_max=60 --fixed nA=25 --fixed nAc=100'
Q='--set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1 --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec" --extra-target "quit/m exp=0.0226:2.0" --extra-target "quit/m rec=0.0210:2.0" --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001 --diff-step 0.04'
```

Use `eval "python3 -u scripts/calibrate_ls.py ... $M $Q --max-nfev 5 --tag <tag>_full"` so that the quoted target
names in `$Q` are parsed. `calibrate_ls.py` takes the starting parameters (`x`) from the `--x0` file and the fixed
fields ONLY from the `--fixed` flags, so the old files' beta, a_max and nA do not carry over.

**Cost.** One full-grid evaluation takes about 770 s with all cores at nAc 100 (`grid_check_v7cmb_full_c100.out`), 5.6
times the old solver; `--max-nfev 5` is about 70 evaluations, about 15 hours of the whole machine per calibration. Run
two calibrations at a time with `export KEAM_NJOBS=$(( $(nproc) / 2 ))` for each (same throughput, both results at
once). Every downstream script (`channels.py`, `run_pipeline.sh`, `cohorts_refined.py`, `jacobian_channels.py`,
`asset_distribution.py`) reads the fixed fields from the calibration file, so beta, r_a, a_max, nA and nAc carry
through automatically.

## Steps

1. Calibrations, in this priority order, two at a time. First pair:
   - v7ck (KPR, additive shock, wage cut): `--x0 output/final_calib_v7cmb_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 $M $Q --max-nfev 5 --tag v7ck_full`
   - v4ek (separable CRRA gamma 2, wage cut): `--x0 output/final_calib_v4emb_full.json --fixed ui_rec_mult=0.5 $M $Q --max-nfev 5 --tag v4ek_full`
   Second pair:
   - v9nk (KPR proportional shock, no wage cut): `--x0 output/final_calib_v9nmb_full.json --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 --fixed kT_mult=1 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q --max-nfev 5 --tag v9nk_full`
   - v4nbk (separable, no wage cut): `--x0 output/final_calib_v4nbmc_full.json --fixed ui_rec_mult=0.5 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0 $M $Q --max-nfev 5 --tag v4nbk_full`
   After each: `python3 scripts/calib_table.py output/final_calib_<tag>_full.json > output/final_calib_<tag>_full.md`
   and `python3 -u scripts/asset_distribution.py --calib output/final_calib_<tag>_full.json --tag <tag>` (report mean
   assets in months of income and the share near a_max); commit (`git add -f` the `.log` and `.out` files) and push.
   When a run stops at the evaluation limit with the objective still falling by more than 10% over its last
   iteration, run a second polish from its result (`--x0 output/final_calib_<tag>_full.json`, same flags,
   `--max-nfev 4`, `--tag <tag>b_full`) and carry the better one.
2. Channel decomposition for each finished calibration, as soon as it is done and cores are free:
   `python3 -u scripts/channels.py --calib output/final_calib_<tag>_full.json --tag <tag>` (10 variants, about 2-3 hours
   each with all cores; `--variants quick` for the two worse fits); commit and push.
3. After all four calibrations and channels: the full pipeline for the best-fitting KPR version and the best-fitting
   separable version, one at a time with all cores (about 20 hours each at nAc 100):
   `bash scripts/run_pipeline.sh output/final_calib_<tag>_full.json <tag>`, then
   `python3 -u scripts/cohorts_refined.py --calib output/final_calib_<tag>_full.json --e-mode relative --tag <tag>_rel`
   and `python3 -u scripts/jacobian_channels.py --tag <tag>`; commit and push after each script.
4. `python3 scripts/write_paper_results.py --specs "KPR wage cut (v7ck)=v7ck,CRRA wage cut (v4ek)=v4ek,KPR no wage cut (v9nk)=v9nk,CRRA no wage cut (v4nbk)=v4nbk" --out ../../PAPER_RESULTS_quit.md`
   with the columns that exist (rerun it as more arrive); write the outcome (objectives, the 11 targets' fit, mean
   assets, precaution / hoarding shares) into `CLAUDE.md` state item 9; commit and push.
