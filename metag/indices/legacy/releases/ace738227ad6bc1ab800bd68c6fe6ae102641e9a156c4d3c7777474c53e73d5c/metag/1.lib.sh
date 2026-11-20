#!/usr/bin/env bash
# One reference library; build/update only its named optimized group index.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
case "${1:-}" in
    -h|--help|--help-full)
        cat <<'HELP'
Usage: ./1.lib.sh [build|check|status|kraken|human|sylph] [OPTIONS]

  build   Build/update lib/group/optimized; reuse prepared data and full indexes.
          Reuse a matching native capacity estimate, otherwise estimate once.
          An intact matching optimized index is reused without rebuilding.
  check   Quick catalog/size/header checks; --full submits full union SHA checks.
  status  Show legacy, optimized and full index availability.
  kraken|human|sylph  Explicit expert reference builders.

  --index-use optimized|legacy|full  Index selection [optimized]
  --max-gib N               Optimized hash cap [config: 384 GiB; not total RAM]
  --estimate auto|reuse|native  Capacity estimation [auto]
  --submit hpc|local|auto   Execution mode [auto; HPC on the cluster]
  --plan-only               Generate commands without submitting
  --config FILE            Site configuration [alongside this script]
  --queue NAME --cores N --threads N --mem-gib N --hours N

Only an owned optimized index is replaced when its inputs/settings change.
The legacy group index, prepared references and full category indexes are reused.
Default build: one job, queue 1t81c, 32 cores, 640 GiB, 120 hours.
Both named group indexes live in the same lib; no lib0 dependency or data copy.
legacy and optimized use the same classifier with different screening settings.

Examples:
  ./1.lib.sh build
  ./1.lib.sh status
  ./1.lib.sh check --index-use legacy
HELP
        exit 0 ;;
esac
source "$HERE/f/env.sh"
exec "$METAG_PYTHON" "$HERE/f/library_cli.py" "$@"
