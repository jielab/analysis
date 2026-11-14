#!/usr/bin/env bash
#BSUB -J "metag_split.england[1]%1"
#BSUB -q v3-6t
#BSUB -n 24
#BSUB -W 120:00
#BSUB -R "span[hosts=1] rusage[mem=1792GB/host]"
#BSUB -o /work/sph-huangj/analysis/metag/logs/split/run4/split.1.%I.out
#BSUB -e /work/sph-huangj/analysis/metag/logs/split/run4/split.1.%I.err

set -euo pipefail
export submit="${submit:-hpc}" WORKFLOW_WORKER=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- /work/sph-huangj/analysis/metag/releases/0fef49ee3e234ea0199c839dd5a60cf84fd54f4d0ff4d85f7f7f7f0fc11e9eff/metag
exec bash /work/sph-huangj/analysis/metag/releases/0fef49ee3e234ea0199c839dd5a60cf84fd54f4d0ff4d85f7f7f7f0fc11e9eff/metag/f/entry.sh parallel-split-worker --config /work/sph-huangj/analysis/metag/jobs/split/run4/config.json --job-plan /work/sph-huangj/analysis/metag/jobs/split/run4/jobs.json --batch-id "${LSB_JOBINDEX:?}" --offset 0
