# Job for a Claude Code session on the workstation

Run in a terminal on the workstation, in the repository folder, on branch `claude/hopeful-ride-vbnou4`:

```
git pull origin claude/hopeful-ride-vbnou4
cd Code26/python && pip install numpy scipy pandas openpyxl xlrd
export KEAM_NJOBS=$(nproc)
nohup bash -c "python3 -u scripts/calibrate_childcare.py --fixed gamma=1 --fixed ui_rec_mult=0.5 --x0 output/x0_gamma1.json --extra alpha_h,e_max,kappa_h_power --maxfev 250 --tag v5_coarse > output/calib_v5_coarse.out 2>&1; python3 -u scripts/calibrate_ls.py --x0 output/final_calib_v5_coarse.json --fixed gamma=1.0 --fixed ui_rec_mult=0.5 --max-nfev 6 --diff-step 0.04 --tag v5_full > output/calib_v5_full.out 2>&1; python3 -u scripts/channels.py --calib output/final_calib_v5_full.json --tag v5 > output/channels_v5.out 2>&1; python3 scripts/calib_table.py output/final_calib_v5_full.json > output/final_calib_v5_full.md" > output/chain_v5_workstation.out 2>&1 &
```

This is the log-utility (gamma = 1) calibration: coarse Nelder-Mead from the rescaled start, full-grid
least-squares polish, then the ten-variant channel decomposition. With 16 cores each 100-type evaluation
takes about 15-20 s, so the whole chain needs roughly 1.5 hours. When it finishes, commit and push:

```
cd ../.. && git add Code26/python/output/final_calib_v5_* Code26/python/output/channels_v5.* && git add -f Code26/python/output/final_calib_v5_full.log && git commit -m "Log-utility calibration (version 5) and its channel decomposition" && git push origin claude/hopeful-ride-vbnou4
```

If a Claude Code session runs on the workstation, paste this file's content as its instruction; it can then
also run `bash scripts/run_pipeline.sh output/final_calib_v5_full.json v5` (full results, about 40 minutes
on 16 cores) and push the outputs.
