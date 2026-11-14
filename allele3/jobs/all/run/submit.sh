#!/usr/bin/env bash
set -euo pipefail
# Re-enters the receipt-aware dispatcher; never blindly resubmits each .cmd.
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/allele3/f/batch_runner.py dispatch --frozen-plan /work/sph-huangj/analysis/allele3/jobs/all/run/plan.json --out /work/sph-huangj/analysis/allele3/jobs/all/run "$@"
