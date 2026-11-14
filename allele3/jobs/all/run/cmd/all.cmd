#!/usr/bin/env bash
#BSUB -J "allele3.all"
#BSUB -q 33
#BSUB -n 40
#BSUB -W 120:00
#BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
#BSUB -o /work/sph-huangj/analysis/allele3/logs/jobs/all/run/all.out
#BSUB -e /work/sph-huangj/analysis/allele3/logs/jobs/all/run/all.err

set -euo pipefail
[[ ${LSB_JOBID:-} =~ ^[1-9][0-9]*$ ]] || { echo 'Submit this .cmd with bsub; no login-node payload' >&2; exit 2; }
[[ ${LSB_DJOB_NUMPROC:-0} =~ ^[0-9]+$ && ${LSB_DJOB_NUMPROC:-0} -ge 40 ]] || { echo 'Insufficient allocated CPU slots' >&2; exit 2; }
export ALLELE_1KG_PANEL=/work/sph-huangj/data/gen/1kg/38/samples.txt
export ALLELE_1KG_VCF_DIR=/work/sph-huangj/data/gen/1kg/38/vcf
export ALLELE_ALL_OF_US=/work/sph-huangj/data/annot/allele3/allofus/hg38_allofus.txt.gz
export ALLELE_ANALYSIS_ROOT=/work/sph-huangj/analysis/allele3
export ALLELE_ANCIENT_DIR=/work/sph-huangj/data/gen/ancient
export ALLELE_ARCHAIC_CACHE=/work/sph-huangj/analysis/allele3/archaic_overlap
export ALLELE_ARCHAIC_DATASETS=1kg
export ALLELE_ARCHAIC_LIFTOVER=/work/sph-huangj/software/bin/liftOver
export ALLELE_ARCHAIC_ROOT=/work/sph-huangj/data/gen
export ALLELE_ATLAS_DIR=/work/sph-huangj/data/annot/allele3/alphagenome/38
export ALLELE_CHAIN37TO38=/work/sph-huangj/files/liftOver/hg19ToHg38.over.chain.gz
export ALLELE_CHAIN38TO37=/work/sph-huangj/files/liftOver/hg38ToHg19.over.chain.gz
export ALLELE_COHORT_TABLE=/work/sph-huangj/data/gen/1kg/38/allele3/herv_cohort/herv_cohort_burden.tsv
export ALLELE_COORDINATES=/work/sph-huangj/analysis/allele3/pca3d/coordinates/all.sscore
export ALLELE_ERV_JOINED=/work/sph-huangj/data/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz
export ALLELE_EXECUTION_MODE=hpc
export ALLELE_FASTA=/work/sph-huangj/data/gen/1hgp/38/GRCH38.fasta.gz
export ALLELE_FASTA37=/work/sph-huangj/data/gen/1hgp/37/GRCH37.fasta.gz
export ALLELE_LD_BLOCK_DIR=/work/sph-huangj/data/LD/block
export ALLELE_LOW_MAPPABILITY_BED=/work/sph-huangj/files/allele3/umap/hg38.k100.lt1.bed.gz
export ALLELE_MEMORY_LIMIT_GB=160
export ALLELE_METADATA=/work/sph-huangj/analysis/allele3/pca3d/tables/sample_metadata.tsv
export ALLELE_PYTHON=/work/sph-huangj/.conda/envs/allele3-evolution/bin/python
export ALLELE_REPLACE=TRUE
export ALLELE_SCAN_STORE=/work/sph-huangj/analysis/allele3/gnomad_scan
export ALLELE_VCF_DIR=/work/sph-huangj/data/annot/allele3/gnomad/sites
export ANALYSIS_ROOT=/work/sph-huangj/analysis
export DIR_ANALYSIS=/work/sph-huangj/analysis
export DIR_ANNOT=/work/sph-huangj/data/annot
export DIR_BIG=/work/sph-huangj/data
export DIR_DAT=/work/sph-huangj
export DIR_DATA=/work/sph-huangj/data
export DIR_FILES=/work/sph-huangj/files
export DIR_GEN=/work/sph-huangj/data/gen
export DIR_GWAS=/work/sph-huangj/data/gwas
export DIR_HOME=/work/sph-huangj
export DIR_LD=/work/sph-huangj/data/LD
export DIR_METAG=/work/sph-huangj/data/metag
export DIR_RUNTIME_TMP=/work/sph-huangj/tmp
export DIR_SCRATCH=/work/sph-huangj
export DIR_SOFTWARE=/work/sph-huangj/software
export DIR_TMP=/work/sph-huangj
export DIR_UKB=/work/sph-huangj/phe/ukb
export DIR_WORK=/work/sph-huangj
export HERV_BUILD=38
export HERV_LIBRARY=/work/sph-huangj/data/annot/allele3/dfam_4.0_human/human.curated.fa
export HERV_PYTHON=/work/sph-huangj/.conda/envs/allele3-evolution/bin/python
export OMP_NUM_THREADS=1
export OS_AXT_ROOT=/work/sph-huangj/data/annot/axt/37
export OS_FILES_ROOT=/work/sph-huangj/files
export OS_GU_CONDA=/work/sph-huangj/software/gu/conda
export OS_SNP_LIST=/work/sph-huangj/phe/ukb/common/snp.lst
export PATH=/work/sph-huangj/.conda/envs/allele3-evolution/bin:/work/sph-huangj/software/bin:/work/sph-huangj/apps/SAIGE/.pixi/envs/default/bin:/work/sph-huangj/.pixi/bin:/work/sph-huangj/.local/bin:/work/sph-huangj/bin:/usr/share/Modules/bin:/opt/ibm/lsfsuite/ext/ppm/10.2/bin:/opt/ibm/lsfsuite/ext/ppm/10.2/linux2.6-glibc2.3-x86_64/bin:/opt/ibm/lsfsuite/lsf/10.1/linux2.6-glibc2.3-x86_64/etc:/opt/ibm/lsfsuite/lsf/10.1/linux2.6-glibc2.3-x86_64/bin:/opt/ibm/jre/bin:/share/base/hpctool:/usr/lpp/mmfs/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin:/usr/bin
export PCA3D_ANALYSIS_ROOT=/work/sph-huangj/analysis/allele3
export PCA3D_ANCIENT_DIR=/work/sph-huangj/data/gen/ancient
export PCA3D_COMMAND=run
export PCA3D_EIGENVEC=/work/sph-huangj/analysis/allele3/pca3d/coordinates/all.sscore
export PCA3D_HQ_VARIANTS=/work/sph-huangj/data/gen/ancient/HighQualityVariants.tsv
export PCA3D_METADATA=/work/sph-huangj/analysis/allele3/pca3d/tables/sample_metadata.tsv
export PCA3D_PGEN_PREFIX=/work/sph-huangj/data/gen/ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes
export PCA3D_PLINK2=plink2
export PCA3D_PROJECT_BASIS=/work/sph-huangj/data/gen/ancient/allele3/pca3d/pca_projection_unique
export PCA3D_PROJECT_INFO=/work/sph-huangj/data/gen/ukb/37/allele3/pca_projection/ukb_info.tsv.gz
export PCA3D_PROJECT_TARGET_GLOB='/work/sph-huangj/data/gen/ukb/37/imp/chr[0-9]*.pgen'
export PCA3D_PYTHON=/work/sph-huangj/.conda/envs/allele3-evolution/bin/python
export PCA3D_REUSE_BASIS=/work/sph-huangj/analysis/allele3/pca3d/basis
export PYTHONDONTWRITEBYTECODE=1
export SCRATCH_ANALYSIS_ROOT=/work/sph-huangj/analysis
export SCRIPTS_ROOT=/work/sph-huangj/scripts
export SHARED_PYTHON=/work/sph-huangj/.conda/envs/py311/bin/python3
export _OS_DEFAULT_ALLELE_ANALYSIS_ROOT=/work/sph-huangj/analysis/allele3
export _OS_DEFAULT_ALLELE_ARCHAIC_ROOT=/work/sph-huangj/data/gen
export _OS_DEFAULT_ALLELE_JOB_QUEUE=33
export _OS_DEFAULT_ANALYSIS_ROOT=/work/sph-huangj/analysis
export _OS_DEFAULT_DIR_ANALYSIS=/work/sph-huangj/analysis
export _OS_DEFAULT_DIR_ANNOT=/work/sph-huangj/data/annot
export _OS_DEFAULT_DIR_BIG=/work/sph-huangj/data
export _OS_DEFAULT_DIR_DAT=/work/sph-huangj
export _OS_DEFAULT_DIR_DATA=/work/sph-huangj/data
export _OS_DEFAULT_DIR_FILES=/work/sph-huangj/files
export _OS_DEFAULT_DIR_GEN=/work/sph-huangj/data/gen
export _OS_DEFAULT_DIR_GWAS=/work/sph-huangj/data/gwas
export _OS_DEFAULT_DIR_HOME=/work/sph-huangj
export _OS_DEFAULT_DIR_LD=/work/sph-huangj/data/LD
export _OS_DEFAULT_DIR_METAG=/work/sph-huangj/data/metag
export _OS_DEFAULT_DIR_RUNTIME_TMP=/tmp
export _OS_DEFAULT_DIR_SCRATCH=/work/sph-huangj
export _OS_DEFAULT_DIR_SOFTWARE=/work/sph-huangj/software
export _OS_DEFAULT_DIR_TMP=/work/sph-huangj
export _OS_DEFAULT_DIR_UKB=/work/sph-huangj/phe/ukb
export _OS_DEFAULT_DIR_WORK=/work/sph-huangj
export _OS_DEFAULT_GRID_ANALYSIS_ROOT=/work/sph-huangj/analysis/grid
export _OS_DEFAULT_GRID_CONDA_ENV=/work/sph-huangj/.conda/envs/grid
export _OS_DEFAULT_GRID_HPC_CONFIG=/work/sph-huangj/scripts/grid/config/hpc.yml
export _OS_DEFAULT_GRID_OUTPUT_ROOT=/work/sph-huangj/analysis/grid
export _OS_DEFAULT_GRID_PLANNER_PYTHON=/work/sph-huangj/.conda/envs/grid/bin/python
export _OS_DEFAULT_GRID_PYTHON=/work/sph-huangj/.conda/envs/grid/bin/python
export _OS_DEFAULT_GU_ANALYSIS_BASE=/work/sph-huangj/analysis/gu
export _OS_DEFAULT_GU_QUEUE=33
export _OS_DEFAULT_GWAS_ANALYSIS_ROOT=/work/sph-huangj/analysis/gwas
export _OS_DEFAULT_GWAS_QUEUE=33
export _OS_DEFAULT_GWAS_RUNTIME_ENV=/work/sph-huangj/scripts/gwas/f/runtime.env
export _OS_DEFAULT_HPC_CPU_QUEUE=33
export _OS_DEFAULT_HPC_GPU_MEMORY_GIB=160
export _OS_DEFAULT_HPC_GPU_QUEUE=4a100-40
export _OS_DEFAULT_HPC_HIGHMEM_QUEUE=1t81c
export _OS_DEFAULT_HPC_HUGEMEM_QUEUE=v3-6t
export _OS_DEFAULT_HPC_MEDIUM_QUEUE=38
export _OS_DEFAULT_HPC_WALLTIME=120:00
export _OS_DEFAULT_LD_ANALYSIS_ROOT=/work/sph-huangj/analysis/ld
export _OS_DEFAULT_LD_QUEUE=33
export _OS_DEFAULT_LE8_ANALYSIS_ROOT=/work/sph-huangj/analysis/le8
export _OS_DEFAULT_LE8_GPU_QUEUE=4a100-40
export _OS_DEFAULT_LE8_PYTHON=/work/sph-huangj/.conda/envs/le8/bin/python3
export _OS_DEFAULT_LE8_QUEUE=33
export _OS_DEFAULT_LE8_RSCRIPT=/work/sph-huangj/.conda/envs/le8/bin/Rscript
export _OS_DEFAULT_MEMORY_DEFAULT_GIB=160
export _OS_DEFAULT_METAG_READS_RAW=/work/sph-huangj/data/metag/raw
export _OS_DEFAULT_METAG_REFERENCE_RAW=/work/sph-huangj/data/annot/metag/raw
export _OS_DEFAULT_METAG_ROOT=/work/sph-huangj/data/metag
export _OS_DEFAULT_METAG_STATE_ROOT=/work/sph-huangj/analysis/metag
export _OS_DEFAULT_METAG_TOOLS_DIR=/work/sph-huangj/software/metag_parallel
export _OS_DEFAULT_OS_AXT_ROOT=/work/sph-huangj/data/annot/axt/37
export _OS_DEFAULT_OS_FILES_ROOT=/work/sph-huangj/files
export _OS_DEFAULT_OS_GU_CONDA=/work/sph-huangj/software/gu/conda
export _OS_DEFAULT_OS_SNP_LIST=/work/sph-huangj/phe/ukb/common/snp.lst
export _OS_DEFAULT_PERSISTENT_ANALYSIS_ROOT=/work/sph-huangj/analysis
export _OS_DEFAULT_SAIGE_ENV_DIR=/work/sph-huangj/.conda/envs/saige
export _OS_DEFAULT_SCRATCH_ANALYSIS_ROOT=/work/sph-huangj/analysis
export _OS_DEFAULT_SHARED_PYTHON=/work/sph-huangj/.conda/envs/py311/bin/python3
export _OS_DEFAULT_jobs_N=8
export jobs_N=8
export submit=hpc SUBMIT=hpc ALLELE_EXECUTION_MODE=hpc
export WORKFLOW_WORKER=1 ALLELE_IN_COMPUTE=1 ALLELE_BATCH_WORKER=1
export ALLELE_ALLOCATED_MEMORY_GIB=160
export ALLELE_ALLOCATED_CORES=40
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- /work/sph-huangj/scripts/allele3
# Payload (argument quoting shown for review): /work/sph-huangj/.conda/envs/allele3-evolution/bin/python /work/sph-huangj/scripts/allele3/f/all_workflow.py worker --out /work/sph-huangj/analysis/allele3 --steps-json '[{"name": "environment", "argv": ["bash", "/work/sph-huangj/scripts/allele3/install.sh", "--check", "--stage", "all"]}, {"name": "hotspots", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "hotspots", "--scan-dir", "/work/sph-huangj/analysis/allele3/gnomad_scan", "--out", "/work/sph-huangj/analysis/allele3/hotspots", "--with-atlas", "TRUE", "--atlas-resources", "/work/sph-huangj/scripts/allele3/config/alphagenome_resources.vcf.tsv"]}, {"name": "density", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "density", "--scan-dir", "/work/sph-huangj/analysis/allele3/gnomad_scan", "--out", "/work/sph-huangj/analysis/allele3/density"]}, {"name": "herv_enrichment", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "herv_enrichment", "--ibdmix", "FALSE", "--scan-dir", "/work/sph-huangj/analysis/allele3/gnomad_scan", "--out", "/work/sph-huangj/analysis/allele3/herv_enrichment", "--metadata", "/work/sph-huangj/analysis/allele3/pca3d/tables/sample_metadata.tsv", "--coordinates", "/work/sph-huangj/analysis/allele3/pca3d/coordinates/all.sscore", "--basis-dir", "/work/sph-huangj/analysis/allele3/pca3d/basis"]}, {"name": "allele_time", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "base_occupancy", "allele-time", "--catalog", "/work/sph-huangj/analysis/allele3/hotspots/runtime/characterize/selection_catalog.tsv", "--pgen-prefix", "/work/sph-huangj/data/gen/ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes", "--metadata", "/work/sph-huangj/analysis/allele3/pca3d/tables/sample_metadata.tsv", "--coordinates", "/work/sph-huangj/analysis/allele3/pca3d/coordinates/all.sscore", "--basis-dir", "/work/sph-huangj/analysis/allele3/pca3d/basis", "--hq-variants", "/work/sph-huangj/data/gen/ancient/HighQualityVariants.tsv", "--fasta38", "/work/sph-huangj/data/gen/1hgp/38/GRCH38.fasta.gz", "--fasta37", "/work/sph-huangj/data/gen/1hgp/37/GRCH37.fasta.gz", "--chain38to37", "/work/sph-huangj/files/liftOver/hg38ToHg19.over.chain.gz", "--chain37to38", "/work/sph-huangj/files/liftOver/hg19ToHg38.over.chain.gz", "--out", "/work/sph-huangj/analysis/allele3/base_occupancy/allele_time"]}, {"name": "space", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "base_occupancy", "space", "--time-dir", "/work/sph-huangj/analysis/allele3/base_occupancy/allele_time", "--pgen-prefix", "/work/sph-huangj/data/gen/ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes", "--metadata", "/work/sph-huangj/analysis/allele3/pca3d/tables/sample_metadata.tsv", "--coordinates", "/work/sph-huangj/analysis/allele3/pca3d/coordinates/all.sscore", "--basis-dir", "/work/sph-huangj/analysis/allele3/pca3d/basis", "--out", "/work/sph-huangj/analysis/allele3/base_occupancy/space"]}, {"name": "base_simulation", "argv": ["bash", "/work/sph-huangj/scripts/allele3/allele3.sh", "base_simulation", "--out", "/work/sph-huangj/analysis/allele3/base_simulation", "--plot"]}]'
exec /work/sph-huangj/.conda/envs/py311/bin/python3 /work/sph-huangj/scripts/allele3/f/batch_runner.py worker --frozen-plan /work/sph-huangj/analysis/allele3/jobs/all/run/plan.json --job-id all --out /work/sph-huangj/analysis/allele3/jobs/all/run
