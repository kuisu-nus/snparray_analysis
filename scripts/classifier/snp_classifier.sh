#!/bin/bash

set -e  # Exit on any error

cd $PROJECT_ROOT

# 基本用法
# python src/models/snp_classifier.py --vcf chr1.vcf.gz --pheno pheno.ped

# 更多选项
python src/models/snp_classifier.py \
  --vcf /home/sukui/01.data/03.raw_data/1kg/gwas/chr1.vcf.gz \
  --pheno /home/sukui/01.data/03.raw_data/1kg/gwas/pheno.All.ped \
  --max-snps 1000 \
  --n-features 5000 \
  --max-depth 8 \
  --learning-rate 0.05 \
  --n-estimators 200 \
  --save-model \
  --save-features