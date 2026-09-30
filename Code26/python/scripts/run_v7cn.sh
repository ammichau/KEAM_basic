#!/bin/bash
# One-shot chain for the paper's model (KPR, additive cost shock, recession wage cut) with the non-search offer
# arrival rate lam_n, at beta 0.993 per month, r 4% per year, fine savings choice grid (nAc 100, a_max 60):
# calibration (13 targets: quit targets, N->E targets, lam_u fixed from the layoff data) with a second polish, channel
# decomposition, full results pipeline, relative-target cohorts, jacobian of the split, PAPER_RESULTS_quit.md.
# Commits and pushes after every step. Safe to launch twice: a lock directory stops the second copy.
# usage (workstation): cd Code26/python && nohup bash scripts/run_v7cn.sh > output/run_v7cn.out 2>&1 &
set -u
cd "$(dirname "$0")/.."
export KEAM_NJOBS=${KEAM_NJOBS:-$(nproc)}
mkdir output/.lock_v7cn 2>/dev/null || { echo "run_v7cn already running (output/.lock_v7cn exists)"; exit 0; }
trap 'rmdir output/.lock_v7cn 2>/dev/null' EXIT
BRANCH=claude/hopeful-ride-vbnou4
commit_push() {   # $1 message, $2.. tags whose .log/.out files are added
  local msg=$1; shift
  git add output ../../PAPER_RESULTS_quit.md 2>/dev/null
  for t in "$@"; do git add -f output/*"$t"*.log output/*"$t"*.out 2>/dev/null; done
  git commit -q -m "$msg" || return 0
  for d in 2 4 8 16 32; do git push -u origin "$BRANCH" && return 0; sleep $d; done
}
echo "run_v7cn start $(date -u)"
CAL=(python3 -u scripts/calibrate_ls.py --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5
     --fixed beta=0.993 --fixed r_a=0.00327 --fixed a_max=60 --fixed nA=25 --fixed nAc=100
     --set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1
     --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec"
     --extra-target "quit/m exp=0.0226:2.0" --extra-target "quit/m rec=0.0210:2.0"
     --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001
     --extra lam_n0,lam_n_ratio --set lam_n0=0.05 --set lam_n_ratio=1.0
     --extra-target "N->E/m exp=0.0530:1.0" --extra-target "N->E/m rec=0.0524:1.0"
     --diff-step 0.04)
# 1. calibration from v7cmb (the best previous KPR fit), about 80 evaluations
"${CAL[@]}" --x0 output/final_calib_v7cmb_full.json --max-nfev 6 --tag v7cn_full > output/calibrate_v7cn_full.out 2>&1
python3 scripts/calib_table.py output/final_calib_v7cn_full.json > output/final_calib_v7cn_full.md
python3 -u scripts/asset_distribution.py --calib output/final_calib_v7cn_full.json --tag v7cn > output/asset_distribution_v7cn.out 2>&1
commit_push "Calibration v7cn (KPR, non-search offer arrival, beta 0.993): $(head -1 output/final_calib_v7cn_full.md)" v7cn
# 2. second polish, about 55 evaluations; keep the better objective
"${CAL[@]}" --x0 output/final_calib_v7cn_full.json --max-nfev 4 --tag v7cnb_full > output/calibrate_v7cnb_full.out 2>&1
python3 scripts/calib_table.py output/final_calib_v7cnb_full.json > output/final_calib_v7cnb_full.md
python3 -u scripts/asset_distribution.py --calib output/final_calib_v7cnb_full.json --tag v7cnb > output/asset_distribution_v7cnb.out 2>&1
BEST=$(python3 -c "import json; a=json.load(open('output/final_calib_v7cn_full.json'))['obj']; b=json.load(open('output/final_calib_v7cnb_full.json'))['obj']; print('v7cnb' if b < a else 'v7cn')")
echo "best calibration: $BEST $(date -u)"
commit_push "Second polish v7cnb: $(head -1 output/final_calib_v7cnb_full.md); carried: $BEST" v7cnb
CALIB=output/final_calib_${BEST}_full.json
# 3. channel decomposition (precaution versus hoarding)
python3 -u scripts/channels.py --calib "$CALIB" --tag "$BEST" > "output/channels_${BEST}.out" 2>&1
commit_push "Channel decomposition for $BEST" "channels_${BEST}"
# 4. full results pipeline, relative-target cohorts, jacobian of the split
bash scripts/run_pipeline.sh "$CALIB" "$BEST" > "output/pipeline_${BEST}.out" 2>&1
commit_push "Full results pipeline for $BEST" "$BEST"
python3 -u scripts/cohorts_refined.py --calib "$CALIB" --e-mode relative --tag "${BEST}_rel" > "output/cohorts_refined_${BEST}_rel.out" 2>&1
commit_push "Refined cohort accounting with relative employment targets for $BEST" "${BEST}_rel"
python3 -u scripts/jacobian_channels.py --calib "$CALIB" --tag "$BEST" > "output/jacobian_channels_${BEST}.out" 2>&1
commit_push "Sensitivity of the precaution/hoarding split for $BEST" "jacobian_channels_${BEST}"
# 5. paper results
python3 scripts/write_paper_results.py --specs "KPR with offer arrival from non-participation (${BEST})=${BEST}" --out ../../PAPER_RESULTS_quit.md
commit_push "PAPER_RESULTS_quit.md for $BEST"
echo "run_v7cn end $(date -u)"
