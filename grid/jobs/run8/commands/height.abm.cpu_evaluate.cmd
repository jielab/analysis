#!/usr/bin/env bash
#BSUB -J grid.height
#BSUB -q 38
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=32GB/host]"
#BSUB -W 24:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run8/logs/height.abm.cpu_evaluate.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run8/logs/height.abm.cpu_evaluate.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 16 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: height.abm.cpu_radius_twins
# Scientific command: bash /work/sph-huangj/scripts/grid/grid.sh evaluate --traits height --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.yml --device cpu --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run8/plan.json --key height.abm.cpu_evaluate
