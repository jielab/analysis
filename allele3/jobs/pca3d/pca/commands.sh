#!/usr/bin/env bash
set -euo pipefail
# Reviewable offline PCA workflow; no downloads.
# Incomplete runs recompute stages; only a checksum-verified full result is reused.
python3 - /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/pca_plan.json <<'ALLELE3_VERIFY_PCA_INPUTS'
import hashlib,json,pathlib
p=json.loads(pathlib.Path(__import__('sys').argv[1]).read_text())
for spec in p['inputs'].values():
 f=pathlib.Path(spec['path']); s=f.stat(); h=hashlib.sha256()
 with f.open('rb') as r:
  h.update(r.read(65536))
  if s.st_size>65536:
   r.seek(max(65536,s.st_size-65536)); h.update(r.read(65536))
 if spec.get('sha256'):
  full=hashlib.sha256()
  with f.open('rb') as source:
   for block in iter(lambda:source.read(1048576),b''): full.update(block)
  if full.hexdigest()!=spec['sha256']: raise SystemExit('PCA input full SHA256 changed: '+str(f))
 if s.st_size!=spec['size'] or s.st_mtime_ns!=spec['mtime_ns'] or h.hexdigest()!=spec['edge_sha256']:
  raise SystemExit('PCA input changed; generate a new plan: '+str(f))
for name,want in p['plan_file_sha256'].items():
 h=hashlib.sha256()
 with open(name,'rb') as r:
  for b in iter(lambda:r.read(1048576),b''): h.update(b)
 if h.hexdigest()!=want: raise SystemExit('PCA plan modified: '+name)
ALLELE3_VERIFY_PCA_INPUTS
if [[ -f /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/PCA_COMPLETE ]]; then
  [[ $(cat /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/PCA_COMPLETE) == 82e30192096606ff1f7a53a527355aceba64a12e6c06aa004374972cb5f172e8 ]] || exit 1
  sha256sum --check --status /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/completed.sha256
  exit 0
fi
plink2 --version
# Stage 1: see pca_plan.json for scientific assumptions.
plink2 --threads 4 --memory 16384 --pfile /work/sph-huangj/data/gen/ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes --keep /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/basis.keep --extract /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/background.ids --chr 1-22 --snps-only just-acgt --min-alleles 2 --max-alleles 2 --maf 0.05 --geno 0.1 --make-pgen --out /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/basis_common
# Stage 2: see pca_plan.json for scientific assumptions.
plink2 --threads 4 --memory 16384 --pfile /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/basis_common --indep-pairwise 200kb 1 0.1 --out /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ld
# Stage 3: see pca_plan.json for scientific assumptions.
plink2 --threads 4 --memory 16384 --pfile /work/sph-huangj/data/gen/ancient/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes --keep /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all.keep --extract /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ld.prune.in --make-pgen --out /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all_pruned
# Stage 4: see pca_plan.json for scientific assumptions.
plink2 --threads 4 --memory 16384 --pfile /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all_pruned --keep /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/basis.keep --freq counts --pca 10 allele-wts vcols=chrom,ref,alt --seed 20261021 --out /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref
# Stage 5: see pca_plan.json for scientific assumptions.
plink2 --threads 4 --memory 16384 --pfile /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all_pruned --read-freq /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref.acount --score /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref.eigenvec.allele 2 5 header-read variance-standardize cols=+scoresums --score-col-nums 6-15 --out /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all
sha256sum /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/all.sscore /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref.acount /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref.eigenvec.allele /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/ref.eigenval > /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/completed.sha256
printf '%s\n' 82e30192096606ff1f7a53a527355aceba64a12e6c06aa004374972cb5f172e8 > /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/PCA_COMPLETE.partial
mv /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/PCA_COMPLETE.partial /work/sph-huangj/analysis/allele3/tmp/pca3d/work/pca/PCA_COMPLETE
