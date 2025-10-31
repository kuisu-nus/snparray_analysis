#!/bin/bash

# GWAS QC Pipeline Script
# Usage: ./gwas_qc_pipeline.sh [input_prefix] [output_prefix] [maf_threshold] [data_dir]
# data_dir: D:/03.projects/AI.PGT/data/SNP_results/PGT_TLS_ALL1/PLINK_311025_1032/1_QC
set -e  # Exit on any error

# Default parameters
DATA_DIR=${1:-.}
INPUT=${2:-PGT_TLS_ALL}

cd $DATA_DIR

echo "Starting GWAS QC Pipeline"
echo "Input: $INPUT"
echo "MAF threshold: $MAF"
echo "Data directory: $DATA_DIR"

# Step 1: Missingness analysis
echo "Step 1/1: Missingness analysis"
awk '$2 ~ /^cnvi/ {print $2}' $INPUT.map > cnv_snps.txt
plink --file $INPUT --exclude cnv_snps.txt --recode --out ${INPUT}_sex

# 替换XY为Y
echo "XY Y" > chr_fix.txt
plink --file ${INPUT}_sex --update-chr chr_fix.txt --make-bed --out ${INPUT}_sex

plink --bfile ${INPUT}_sex --impute-sex --recode --out ${INPUT}_sex
plink --bfile ${INPUT}_sex  --missing

echo "finished"