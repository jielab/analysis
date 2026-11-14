#!/usr/bin/env bash
#BSUB -J grid.ldl.blockprepare
#BSUB -q 38
#BSUB -n 32
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -W 120:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run8/logs/ldl.abm.cpu_blockprepare.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run8/logs/ldl.abm.cpu_blockprepare.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: shared.cpu_prepare
# Scientific command: bash /work/sph-huangj/scripts/grid/2.csx.sh blockprepare --traits ldl --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.yml --device cpu --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run8/plan.json --key ldl.abm.cpu_blockprepare
