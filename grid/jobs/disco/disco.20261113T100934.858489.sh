#!/usr/bin/env bash
set -euo pipefail
source /work/sph-huangj/scripts/grid/f/0.common.sh
grid_activate_environment
export GRID_CSX_BACKEND=cpu GRID_CSX_GPU=0 GRID_CSX_MODULE=csx
python3 /work/sph-huangj/scripts/grid/f/0.common.py disco-inputs /work/sph-huangj/data/gen/ukb/pca_proj/ukb.discodivas.pca.tsv.gz /work/sph-huangj/files/DiscoDivas/med.g1000.4pop.tsv /work/sph-huangj/data/gen/ukb/pgs/ldl /work/sph-huangj/analysis/grid/tmp/disco/ldl 10 ''  > /work/sph-huangj/analysis/grid/logs/disco/ldl/prepare.log 2>&1
/work/sph-huangj/.conda/envs/grid/bin/Rscript /work/sph-huangj/scripts/grid/f/3.disco.R -m /work/sph-huangj/analysis/grid/tmp/disco/ldl/centers.tsv -p /work/sph-huangj/analysis/grid/tmp/disco/ldl/pca.tsv --prs.list /work/sph-huangj/analysis/grid/tmp/disco/ldl/AFR.tsv\,/work/sph-huangj/analysis/grid/tmp/disco/ldl/EAS.tsv\,/work/sph-huangj/analysis/grid/tmp/disco/ldl/EUR.tsv\,/work/sph-huangj/analysis/grid/tmp/disco/ldl/SAS.tsv -s IID\,PRS -A 1\,1\,1\,1 --regress.PCA TRUE --print.coef TRUE -o /work/sph-huangj/analysis/grid/tmp/disco/ldl/disco  > /work/sph-huangj/analysis/grid/logs/disco/ldl/disco.log 2>&1
python3 /work/sph-huangj/scripts/grid/f/0.common.py disco-output /work/sph-huangj/analysis/grid/tmp/disco/ldl/disco.tsv.gz /work/sph-huangj/analysis/grid/tmp/disco/ldl  > /work/sph-huangj/analysis/grid/logs/disco/ldl/validate.log 2>&1
python3 /work/sph-huangj/scripts/grid/f/0.common.py publish disco /work/sph-huangj/analysis/grid/tmp/disco/ldl/disco.tsv.gz /work/sph-huangj/data/gen/ukb/pgs/ldl/2disco.scores.rds --remove ''  > /work/sph-huangj/analysis/grid/logs/disco/step.1794535829226919996.858543.log 2>&1
/work/sph-huangj/.conda/envs/grid/bin/Rscript /work/sph-huangj/scripts/grid/../0f/results.R import-table /work/sph-huangj/analysis/grid/tmp/disco/ldl/disco.coef.tsv.gz /work/sph-huangj/data/gen/ukb/pgs/ldl/2disco.coefficients.xlsx  > /work/sph-huangj/analysis/grid/logs/disco/step.1794535838855203384.858543.log 2>&1
