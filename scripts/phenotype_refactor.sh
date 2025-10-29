#!/bin/bash

# Phenotype Data Refactor Script
# Usage: Run this script from the project directory

ROOT="D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434"
ROOT="D:\03.projects\AI.PGT\snparray_analysis\src"
cd ${ROOT}

echo "Starting phenotype data refactoring..."
echo "Working directory: $(pwd)"

python data_process/phenotype_refactor.py \
    --ped-file "D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.ped" \
    --phenotype-file "D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.phenotype" \
    --full-phenotype-file "D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.csv" \
    --output-ped "D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.refactor.ped" \
    --output-phenotype "D:\03.projects\AI.PGT\data\SNP_results\PGT_TLS\PLINK_281025_0434\PGT_TLS.refactor.phenotype" \
    --output-full-phenotype "D:\03.projects\AI.PGT\snparray_analysis\work_dir\clinical_snparray_phenotype.refactor.csv" \
    --rename "男方姓名:PID" \
    --rename "女方姓名_clinical:MID" \
    --rename "诊断可移植:PHENO" \
    --log-level "INFO"

echo "Phenotype refactoring completed!"