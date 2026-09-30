#!/bin/bash
# One-shot chain for the paper's model (author, 2026-09-30): KPR preferences, additive cost shock, NO recession wage
# cut (wife and husband, phi_rec = phi_rec_H = 1), non-search offer arrival lam_n; beta 0.993 per month, r 4% per year,
# nAc 100, a_max 60; the cyclicality of quits and N->E entry targeted as standard deviations of the log rate (the
# convention of the UE target), not as recession levels; the recession job-finding fall stays fixed at 20%.
# Steps: stop the v7cn chain and commit its outputs; calibration v7ns and a second polish v7nsb (the better carried);
# baseline precaution/hoarding split; full results pipeline; relative-target cohorts; jacobian of the split; the channel
# variants; PAPER_RESULTS_quit.md. Commits and pushes after every step.
# usage (workstation): cd Code26/python && nohup bash scripts/run_v7ns.sh > output/run_v7ns.out 2>&1 &
set -u
cd "$(dirname "$0")/.."
export KEAM_NJOBS=${KEAM_NJOBS:-$(nproc)}
mkdir output/.lock_v7ns 2>/dev/null || { echo "run_v7ns already running (output/.lock_v7ns exists)"; exit 0; }
trap 'rmdir output/.lock_v7ns 2>/dev/null' EXIT
BRANCH=claude/hopeful-ride-vbnou4
commit_push() {   # $1 message, $2.. tags whose .log/.out files are added
  local msg=$1; shift
  git add output
  [ -f ../../PAPER_RESULTS_quit.md ] && git add ../../PAPER_RESULTS_quit.md
  for t in "$@"; do for f in output/*"$t"*.log output/*"$t"*.out; do [ -f "$f" ] && git add -f "$f"; done; done
  git commit -q -m "$msg" || return 0
  for d in 2 4 8 16 32; do git pull -q --no-rebase origin "$BRANCH" && git push -u origin "$BRANCH" && return 0; sleep $d; done
}
echo "run_v7ns start $(date -u)"
# 0. stop the v7cn chain (wage-cut version, superseded) and commit its outputs
if [ -d output/.lock_v7cn ]; then
  pkill -f "scripts/run_v7cn.sh"; sleep 2; pkill -f "v7cnb"; sleep 5; rmdir output/.lock_v7cn 2>/dev/null
  echo "v7cn chain stopped $(date -u)"
fi
commit_push "Outputs of the v7cn chain (calibrations v7cn and v7cnb, partial channel decomposition of v7cnb); chain stopped: no-wage-cut re-estimation" v7cn
CAL=(python3 -u scripts/calibrate_ls.py --fixed kpr=1 --fixed gamma=2.0 --fixed ui_rec_mult=0.5 --fixed phi_rec=1.0 --fixed phi_rec_H=1.0
     --fixed beta=0.993 --fixed r_a=0.00327 --fixed a_max=60 --fixed nA=25 --fixed nAc=100
     --set lam_u0=0.0130 --set lam_u1=0.0154 --drop lam_u0 --drop lam_u1
     --drop-target "E->nonE/m exp" --drop-target "E->nonE/m rec" --drop-target "quit/m rec"
     --extra-target "quit/m exp=0.0226:2.0" --extra-target "sd log quit (women)=0.0262:2.0:0.0686"
     --extra-target "N->E/m exp=0.0530:1.0" --extra-target "sd log N->E (women)=0.0045:1.0:0.0686"
     --set lam_f_ratio=0.80 --bound lam_f_ratio:0.8:0.8000001
     --extra lam_n0,lam_n_ratio
     --diff-step 0.04)
# 1. calibration from v7cnb, about 100 evaluations
"${CAL[@]}" --x0 output/final_calib_v7cnb_full.json --max-nfev 6 --tag v7ns_full > output/calibrate_v7ns_full.out 2>&1
python3 scripts/calib_table.py output/final_calib_v7ns_full.json > output/final_calib_v7ns_full.md
python3 -u scripts/asset_distribution.py --calib output/final_calib_v7ns_full.json --tag v7ns > output/asset_distribution_v7ns.out 2>&1
commit_push "Calibration v7ns (no wage cut; quit and N->E cyclicality targeted as sd of the log rate): $(head -1 output/final_calib_v7ns_full.md)" v7ns
# 2. second polish, about 70 evaluations; keep the better objective
"${CAL[@]}" --x0 output/final_calib_v7ns_full.json --max-nfev 4 --tag v7nsb_full > output/calibrate_v7nsb_full.out 2>&1
python3 scripts/calib_table.py output/final_calib_v7nsb_full.json > output/final_calib_v7nsb_full.md
python3 -u scripts/asset_distribution.py --calib output/final_calib_v7nsb_full.json --tag v7nsb > output/asset_distribution_v7nsb.out 2>&1
BEST=$(python3 -c "import json; a=json.load(open('output/final_calib_v7ns_full.json'))['obj']; b=json.load(open('output/final_calib_v7nsb_full.json'))['obj']; print('v7nsb' if b < a else 'v7ns')")
echo "best calibration: $BEST $(date -u)"
commit_push "Second polish v7nsb: $(head -1 output/final_calib_v7nsb_full.md); carried: $BEST" v7nsb
CALIB=output/final_calib_${BEST}_full.json
# 3. baseline precaution/hoarding split (the variants follow at the end; channels.py resumes from its JSON)
python3 -u scripts/channels.py --calib "$CALIB" --tag "$BEST" --variants quick > "output/channels_${BEST}_quick.out" 2>&1
commit_push "Precaution/hoarding split for $BEST (baseline)" "channels_${BEST}"
# 4. full results pipeline, relative-target cohorts, jacobian of the split
bash scripts/run_pipeline.sh "$CALIB" "$BEST" > "output/pipeline_${BEST}.out" 2>&1
commit_push "Full results pipeline for $BEST" "$BEST"
python3 -u scripts/cohorts_refined.py --calib "$CALIB" --e-mode relative --tag "${BEST}_rel" > "output/cohorts_refined_${BEST}_rel.out" 2>&1
commit_push "Refined cohort accounting with relative employment targets for $BEST" "${BEST}_rel"
python3 -u scripts/jacobian_channels.py --calib "$CALIB" --tag "$BEST" > "output/jacobian_channels_${BEST}.out" 2>&1
commit_push "Sensitivity of the precaution/hoarding split for $BEST" "jacobian_channels_${BEST}"
# 5. channel variants, paper results
python3 -u scripts/channels.py --calib "$CALIB" --tag "$BEST" > "output/channels_${BEST}.out" 2>&1
python3 scripts/write_paper_results.py --specs "KPR no wage cut with non-search offers (${BEST})=${BEST}" --out ../../PAPER_RESULTS_quit.md
commit_push "Channel variants and PAPER_RESULTS_quit.md for $BEST" "channels_${BEST}"
echo "run_v7ns end $(date -u)"
