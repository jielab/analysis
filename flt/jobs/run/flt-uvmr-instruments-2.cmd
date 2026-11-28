#!/usr/bin/env bash
set -euo pipefail
exec /work/sph-huangj/.conda/envs/le8/bin/python3 /work/sph-huangj/scripts/flt/f/workflow.py worker /work/sph-huangj/analysis/flt/jobs/run/plan.json flt-uvmr-instruments-2
