#!/usr/bin/env bash
#BSUB -J "metag_split_summary.changj"
#BSUB -q v3-6t
#BSUB -n 2
#BSUB -W 24:00
#BSUB -R "span[hosts=1] rusage[mem=16GB/host]"
#BSUB -o /work/sph-huangj/analysis/metag/logs/split/run3/summary.out
#BSUB -e /work/sph-huangj/analysis/metag/logs/split/run3/summary.err

set -euo pipefail
export submit="${submit:-hpc}" WORKFLOW_WORKER=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- /work/sph-huangj/analysis/metag/releases/0fef49ee3e234ea0199c839dd5a60cf84fd54f4d0ff4d85f7f7f7f0fc11e9eff/metag
exec bash /work/sph-huangj/analysis/metag/releases/0fef49ee3e234ea0199c839dd5a60cf84fd54f4d0ff4d85f7f7f7f0fc11e9eff/metag/f/entry.sh parallel-summary-worker --config /work/sph-huangj/analysis/metag/jobs/split/run3/config.json --job-plan /work/sph-huangj/analysis/metag/jobs/split/run3/jobs.json
