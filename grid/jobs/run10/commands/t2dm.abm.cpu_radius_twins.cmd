#!/usr/bin/env bash
#BSUB -J grid.t2dm.radius_twins
#BSUB -q 38
#BSUB -n 32
#BSUB -R "span[hosts=1] rusage[mem=256GB/host]"
#BSUB -W 120:00
#BSUB -o /work/sph-huangj/analysis/grid/jobs/run10/logs/t2dm.abm.cpu_radius_twins.%J.out
#BSUB -e /work/sph-huangj/analysis/grid/jobs/run10/logs/t2dm.abm.cpu_radius_twins.%J.err

set -euo pipefail
umask 077
[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} -ge 32 ]] || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
# Dependencies: t2dm.abm.cpu_geometry
# Scientific command: bash /work/sph-huangj/scripts/grid/grid.sh radius_twins --traits t2dm --analysis-root /work/sph-huangj/analysis/grid --block-config /work/sph-huangj/scripts/grid/config/abm.blocks.yml --device cpu --threads 8
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/0f/workflow_jobs.py worker --plan /work/sph-huangj/analysis/grid/jobs/run10/plan.json --key t2dm.abm.cpu_radius_twins
