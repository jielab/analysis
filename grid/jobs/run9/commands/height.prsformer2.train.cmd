#!/usr/bin/env bash
#BSUB -J prsformer.height.prsformer2.train
#BSUB -q 2v100-32-sc
#BSUB -n 8
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -W 120:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run9/logs/height.prsformer2.train.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run9/logs/height.prsformer2.train.%J.err
#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 8 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: shared.cpu_prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/4.prsformer.sh prsformer2 --traits height --analysis-root /work/sph-huangj/analysis/grid --workflow-stage model --device cuda --threads 2
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run9/plan.json --key height.prsformer2.train
