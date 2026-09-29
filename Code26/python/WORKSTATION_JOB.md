# Workstation job (2026-09-29, 16:30 UTC): the paper's model only (KPR, additive shock, wage cut) with the non-search
# offer arrival rate; estimate and produce the paper results

You are the executor on the author's workstation for the KEAM project (repository KEAM_basic, branch
`claude/hopeful-ride-vbnou4`). Read `CLAUDE.md` (state items 9 and 10) first. Work autonomously; never modify the
MATLAB files. The author's decision (16:20 UTC): work only on the v7c family; the other specifications are dropped.

## Step 0: stop everything from the earlier jobs

Kill every running `calibrate_ls.py` and any chain script from the 15:00 UTC job (tags v7ck, v4ek, v9nk, v4nbk) and
from the 17:00 UTC draft (v4en, v9nn, v4nbn); delete nothing. Then `git pull origin claude/hopeful-ride-vbnou4`
(commit with `scripts/run_v7cn.sh`, 16:30 UTC or later).

## Step 1: one command

```
cd Code26/python && nohup bash scripts/run_v7cn.sh > output/run_v7cn.out 2>&1 &
```

The script runs the whole chain with all cores and commits and pushes after every step: calibration `v7cn` (KPR, additive
cost shock, recession wage cut; from `final_calib_v7cmb_full.json`; beta 0.993 per month, r 4% per year, `nAc` 100,
`a_max` 60; 13 targets: the 11 quit-targeted ones plus `N->E/m exp` 0.0530 and `N->E/m rec` 0.0524, which identify the
non-search offer arrival `lam_n0` and its recession ratio `lam_n_ratio`; about 80 evaluations at about 770 s each), a
second polish `v7cnb` (about 55 evaluations; the better objective is carried), the channel decomposition, the full
results pipeline, the relative-target cohorts, the jacobian of the split, and `PAPER_RESULTS_quit.md` (root). About
two days of machine time in all. A lock directory `output/.lock_v7cn` stops a second copy: the author may launch the
same command by hand; if the lock exists, do nothing. Poll `output/run_v7cn.out` and the step logs; if a step fails,
write the reason to `BLOCKED.md`, commit, push and stop.
