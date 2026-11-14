#!/bin/bash -l
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -W 120:00
#BSUB -J ukb.aod
#BSUB -o "/work/sph-huangj/analysis/ukb/aod/logs/run/worker.out"
#BSUB -e "/work/sph-huangj/analysis/ukb/aod/logs/run/worker.err"
set -euo pipefail
cd /work/sph-huangj/scripts/ukb
export SCRIPTS_ROOT=/work/sph-huangj/scripts
export DIR_WORK=/work/sph-huangj
export DIR_DAT=/work/sph-huangj
export DIR_BIG=/work/sph-huangj/data
export DIR_DATA=/work/sph-huangj/data
export DIR_TMP=/work/sph-huangj
export DIR_SCRATCH=/work/sph-huangj
export DIR_SOFTWARE=/work/sph-huangj/software
export DIR_FILES=/work/sph-huangj/files
export DIR_UKB=/work/sph-huangj/phe/ukb
export DIR_ANNOT=/work/sph-huangj/data/annot
export DIR_GEN=/work/sph-huangj/data/gen
export DIR_GWAS=/work/sph-huangj/data/gwas
export DIR_LD=/work/sph-huangj/data/LD
export DIR_METAG=/work/sph-huangj/data/metag
export DIR_RUNTIME_TMP=/tmp
export DIR_HOME=/work/sph-huangj
export ANALYSIS_ROOT=/work/sph-huangj/analysis
export submit=hpc
export SUBMIT=hpc
export jobs_N=8
export SUBMIT_WAIT=FALSE
export UKB_THREADS=4
export UKB_SLOTS=40
export UKB_MEMORY_GIB=160
export UKB_QUEUE=33
export UKB_WALLTIME=120:00
export UKB_PYTHON=/work/sph-huangj/.conda/envs/py311/bin/python3
export UKB_RSCRIPT=/work/sph-huangj/.conda/envs/le8/bin/Rscript
export WORKFLOW_WORKER=1 SCRIPT_CONSOLE_ACTIVE=1
[[ -n ${LSB_JOBID:-} ]] || { echo "Submit this .cmd with bsub" >&2; exit 2; }
(( ${LSB_DJOB_NUMPROC:-0} >= 40 )) || { echo "Insufficient allocated CPU slots" >&2; exit 2; }
exec bash /work/sph-huangj/scripts/ukb/aod.sh --execute --Y cvd_cad\,cvd_afib\,cvd_hfail --as-of 2026-11-12
