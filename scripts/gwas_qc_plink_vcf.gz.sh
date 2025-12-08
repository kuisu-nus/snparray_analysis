#!/bin/bash

# GWAS QC Pipeline Script
# Usage: ./gwas_qc_pipeline.sh [input_prefix] [output_prefix] [maf_threshold] [data_dir]
# data_dir: D:/03.projects/AI.PGT/data/SNP_results/PGT_TLS_ALL1/PLINK_311025_1032/1_QC
set -e  # Exit on any error

# Default parameters
DATA_DIR=${1:-.}
INPUT=${2:-PGT_TLS_ALL.refactor}
OUTPUT=${3:-PGT_TLS_ALL.refactor}
MAF=${4:-0.05}

cd $DATA_DIR

echo "Starting GWAS QC Pipeline"
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo "MAF threshold: $MAF"
echo "Data directory: $DATA_DIR"


plink --vcf ${INPUT} --make-founders --make-bed --out ${INPUT}
plink --bfile ${INPUT}  --missing

# Rscript --no-save hist_miss.R

plink --bfile ${INPUT}  --geno 0.2 --make-bed --out ${INPUT}_step1
plink --bfile ${INPUT}_step1 --mind 0.2 --make-bed --out ${INPUT}_step2
plink --bfile ${INPUT}_step2 --geno 0.02 --make-bed --out ${INPUT}_step3
plink --bfile ${INPUT}_step3 --mind 0.02 --make-bed --out ${INPUT}_step4

# Step 2: Sex discrepancy check
# echo "Step 2/6: Sex discrepancy check"
# Rscript --no-save gender_check.R


# Step 3: Autosomal SNPs and MAF filtering
# echo "Step 3/6: Autosomal SNPs and MAF filtering"
# awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ${INPUT}_step5.bim > snp_1_22.txt
# plink --bfile ${INPUT}_step4 --extract snp_1_22.txt --make-bed --out ${INPUT}_step5

plink --bfile ${INPUT}_step4 --freq --out MAF_check
# Rscript --no-save MAF_check.R

plink --bfile ${INPUT}_step4 --maf $MAF --make-bed --out ${INPUT}_step6

# Step 4: Hardy-Weinberg equilibrium
echo "Step 4/6: Hardy-Weinberg equilibrium"
plink --bfile ${INPUT}_step6 --hardy
awk '{ if ($9 <0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe
# Rscript --no-save hwe.R

plink --bfile ${INPUT}_step6 --hwe 1e-6 --make-bed --out ${INPUT}_hwe_temp
plink --bfile ${INPUT}_hwe_temp --hwe 1e-10 --hwe-all --make-bed --out ${INPUT}_step7

# Step 5: Heterozygosity check
echo "Step 5/6: Heterozygosity check"
touch inversion.txt
plink --bfile ${INPUT}_step7 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile ${INPUT}_step7 --extract indepSNP.prune.in --het --out R_check

# Rscript --no-save check_heterozygosity_rate.R
# Rscript --no-save heterozygosity_outliers_list.R
# 杂合度检查 → 移除技术异常样本
#     ↓
# Founders筛选 → 移除已知亲属
#     ↓
# 隐性亲缘检测 → 识别未知的遗传关系
#     ↓
# 缺失率对比 → 智能选择移除个体

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}' > het_fail_ind.txt
plink --bfile ${INPUT}_step7 --remove het_fail_ind.txt --make-bed --out ${INPUT}_step8

# Step 6: Relatedness check
echo "Step 6/6: Relatedness check"
# plink --bfile ${INPUT}_step9 --filter-founders --make-bed --out ${INPUT}_step10
plink --bfile ${INPUT}_step8 --make-founders --make-bed --out ${INPUT}_step9
plink --bfile ${INPUT}_step9 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
plink --bfile ${INPUT}_step9 --missing

# Create related individuals removal list (modify as needed)
# cat > 0.2_low_call_rate_pihat.txt << EOF
# 13291 NA07045
# EOF

# plink --bfile ${INPUT}_step10 --remove 0.2_low_call_rate_pihat.txt --make-bed --out $OUTPUT

echo "GWAS QC Pipeline completed successfully!"
echo "Final output: $OUTPUT"
echo "Files needed for next tutorial:"
echo "- $OUTPUT (bed/bim/fam files)"
echo "- indepSNP.prune.in"