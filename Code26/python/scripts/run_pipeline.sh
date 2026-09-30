#!/bin/bash
# Full results pipeline for one calibration file on the 100-type grid.
# usage: bash scripts/run_pipeline.sh output/final_calib_rho_full.json rho
# Writes output/final_results_<tag>.*, extra_experiments_<tag>.*, cohorts_refined_<tag>.*,
# robustness_final_<tag>.*, figures_<tag>/, irf_careers_<tag>.json, then RESULTS.md (write_results.py).
set -u
cd "$(dirname "$0")/.."
CALIB=$1; TAG=$2
echo "pipeline $CALIB tag=$TAG start $(date -u)"
python3 -u scripts/run_final.py --calib "$CALIB" --full --tag "$TAG" > "output/pipeline_${TAG}_run_final.out" 2>&1 && echo "run_final done $(date -u)"
python3 -u scripts/extra_experiments.py --calib "$CALIB" --tag "$TAG" > "output/pipeline_${TAG}_extra.out" 2>&1 && echo "extra done $(date -u)"
python3 -u scripts/cohorts_refined.py --calib "$CALIB" --tag "$TAG" > "output/pipeline_${TAG}_cohorts.out" 2>&1 && echo "cohorts done $(date -u)"
python3 -u scripts/figures.py --calib "$CALIB" --results "output/final_results_${TAG}.json" --cohorts "output/cohorts_refined_${TAG}.json" --tag "$TAG" > "output/pipeline_${TAG}_figures.out" 2>&1 && echo "figures done $(date -u)"
python3 -u scripts/irf_careers.py --calib "$CALIB" --tag "$TAG" > "output/pipeline_${TAG}_irf.out" 2>&1 && echo "irf done $(date -u)"
python3 -u scripts/robustness_final.py --calib "$CALIB" --results "output/final_results_${TAG}.json" --tag "$TAG" > "output/pipeline_${TAG}_robust.out" 2>&1 && echo "robustness done $(date -u)"
echo "pipeline end $(date -u)"
