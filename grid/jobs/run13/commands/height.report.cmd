#!/usr/bin/env bash
#BSUB -J prsformer.height.report
#BSUB -q 38
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=32GB/host]"
#BSUB -W 24:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run13/logs/height.report.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run13/logs/height.report.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 16 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: height.train
# Scientific command: bash /work/sph-huangj/scripts/grid/4.prsformer.sh all --traits height --analysis-root /work/sph-huangj/analysis/grid
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run13/plan.json --key height.report
