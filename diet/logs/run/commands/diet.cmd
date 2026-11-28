#!/usr/bin/env bash
#BSUB -J diet.disco
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -W 120:00
#BSUB -o /work/sph-huangj/analysis/diet/logs/run/logs/diet.%J.out
#BSUB -e /work/sph-huangj/analysis/diet/logs/run/logs/diet.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 40 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
source /work/sph-huangj/scripts/0f/0os.sh
# Dependencies: (none)
# Scientific command: bash /work/sph-huangj/scripts/diet/1.disco.sh all --cores 40
exec python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/diet/logs/run/plan.json --key diet
