#!/usr/bin/env bash
#BSUB -J grid.architecture
#BSUB -q 38
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=16GB/host]"
#BSUB -W 24:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run10/logs/abm.cpu_architecture.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run10/logs/abm.cpu_architecture.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 16 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: height.abm.cpu_evaluate, ldl.abm.cpu_evaluate, t2dm.abm.cpu_evaluate
# Scientific command: bash /work/sph-huangj/scripts/grid/grid.sh architecture --traits height,ldl,t2dm --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.yml
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run10/plan.json --key abm.cpu_architecture
