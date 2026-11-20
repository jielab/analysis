#!/usr/bin/env bash
set -euo pipefail
source /work/sph-huangj/analysis/metag/indices/legacy/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/env.sh
cd /work/sph-huangj/scripts/metag
/work/sph-huangj/software/venvs/metag-20261115/bin/python /work/sph-huangj/analysis/metag/indices/legacy/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/qc_submit.py --verify-config /work/sph-huangj/analysis/metag/indices/legacy/jobs/qc/hierarchical_split/run/effective-config.json 25ad61d42bd77de5425a2418bbfa8aff9898c57022b4db36dbc3c068c01d5bc5 /work/sph-huangj/analysis/metag/indices/legacy/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag
exec /work/sph-huangj/software/venvs/metag-20261115/bin/python /work/sph-huangj/analysis/metag/indices/legacy/releases/ace738227ad6bc1ab800bd68c6fe6ae102641e9a156c4d3c7777474c53e73d5c/metag/f/hierarchical_split.py --config /work/sph-huangj/analysis/metag/indices/legacy/jobs/qc/hierarchical_split/run/effective-config.json group,species,evaluate,report --replace FALSE --stage both --evaluation basic --routing compact --sample Ndesert/SRR23862437 --threads 22
