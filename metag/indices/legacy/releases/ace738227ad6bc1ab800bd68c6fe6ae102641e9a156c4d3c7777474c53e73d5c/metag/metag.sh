#!/usr/bin/env bash
# One frozen job per batch: finish all selected splits, then evaluate the batch.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
case "${1:-}" in
    -h|--help|--help-full)
        cat <<'HELP'
Usage: ./metag.sh [all|split|split-group|split-species|evaluate|report] [OPTIONS]

  all           Split every sample in a batch, then evaluate/report [default]
  split         Group and selected-species split, then report; no evaluation
  split-group   Group split only
  split-species Refine selected coarse groups only
  evaluate      Evaluate existing splits
  report        Refresh completion/count reports

  --index-use optimized|legacy|full  Index selection [optimized]
      optimized  New group screening configuration: lib/group/optimized
      legacy     Preserved old group index: lib/group/legacy
      full       Existing complete category union; expensive random I/O
  --dat-list CSV          Dataset names
  --sample DATASET/SAMPLE One sample (repeatable)
  --batch 1-2|1,3|all     Per-dataset batch selection [all]
  --batch-size N          Samples per job [10]
  --groups human,plant   Detailed species groups [human,plant]; 'all' supported
  --evaluation basic|full [basic]: QC and reference controls; full also runs
                          full-library agreement and competitive alignment
  --include-uncertain    Refine uncertain bins with the full union (costly)
  --submit hpc|local|auto --plan-only --config FILE
  --queue NAME --cores N --mem-gib N --hours N --threads N

One job per batch: all splits before evaluation. An active matching index build
is automatically added as a job dependency. No extra split/evaluate jobs.
Results: data/metag/species/<index-use>/DATASET/SAMPLE
State/logs: analysis/metag/indices/<index-use>/
Repeated runs reuse compatible stages; --replace TRUE recomputes owned outputs.

Examples:
  ./metag.sh all --dat-list Ndesert --batch 1 --index-use optimized
  ./metag.sh all --dat-list Ndesert --batch 1 --index-use legacy
  ./metag.sh evaluate --sample Ndesert/SRR23862437 --index-use optimized --evaluation full
HELP
        exit 0 ;;
esac
if (($# == 0)) || [[ $1 == --* ]]; then set -- all "$@"; fi
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/metag_submit.py" "$@"
