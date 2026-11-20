#!/usr/bin/env bash
set -euo pipefail
source /work/sph-huangj/analysis/metag/indices/optimized/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/env.sh
cd /work/sph-huangj/scripts/metag
/work/sph-huangj/software/venvs/metag-20261115/bin/python /work/sph-huangj/analysis/metag/indices/optimized/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/qc_submit.py --verify-config /work/sph-huangj/analysis/metag/indices/optimized/jobs/qc/group_library/run2/effective-config.json c960872e7623e4cc4a44205920a334f60743f6b605df31e1331e7d1e1e7fe791 /work/sph-huangj/analysis/metag/indices/optimized/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag
exec /work/sph-huangj/software/venvs/metag-20261115/bin/python /work/sph-huangj/analysis/metag/indices/optimized/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/group_library.py --config /work/sph-huangj/analysis/metag/indices/optimized/jobs/qc/group_library/run2/effective-config.json --max-gib 384.0 --estimate auto --threads 32 --replace TRUE
