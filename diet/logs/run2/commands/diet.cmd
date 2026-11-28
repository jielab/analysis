#!/usr/bin/env bash
#BSUB -J diet
#BSUB -q 33
#BSUB -n 20
#BSUB -R "span[hosts=1] rusage[mem=80GB/host]"
#BSUB -W 120:00
#BSUB -o /work/sph-huangj/analysis/diet/logs/run2/logs/diet.%J.out
#BSUB -e /work/sph-huangj/analysis/diet/logs/run2/logs/diet.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 20 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
source /work/sph-huangj/scripts/0f/0os.sh
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/diet/diet.sh run --cores 20
exec python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/diet/logs/run2/plan.json --key diet
