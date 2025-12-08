#!/bin/bash

# GWAS Association Analysis Script
# Usage: ./gwas_association.sh [data_directory] [input_prefix] [covar_file] [output_prefix] [permutations]

set -e  # Exit on any error

# Default parameters
DATA_DIR=${1:-.}
INPUT=${2:-PGT_TLS.refactor_11}
COVAR=${3:-covar_mds.txt}
OUTPUT=${4:-association_results}
PERMUTATIONS=${5:-1000000}

cd $DATA_DIR


echo "Starting GWAS Association Analysis"
echo "Input: $INPUT"
echo "Covariates: $COVAR"
echo "Output prefix: $OUTPUT"
echo "Permutations: $PERMUTATIONS"

# Step 1: Basic association analysis (no covariates)
echo "Step 1/6: Basic association analysis"
plink --bfile $INPUT --assoc --out assoc_results

# Step 2: Logistic regression with covariates
echo "Step 2/6: Logistic regression with covariates"
plink --bfile $INPUT --covar $COVAR --logistic hide-covar --out logistic_results

# Remove NA values from logistic results
awk '!/NA/' logistic_results.assoc.logistic > logistic_results.assoc_2.logistic

# Step 3: Multiple testing correction
echo "Step 3/6: Multiple testing correction"
plink --bfile $INPUT --assoc --adjust --out adjusted_assoc_results

# Step 4: Permutation testing (subset of chromosome 22)
plink --bfile ${INPUT} --assoc --mperm $PERMUTATIONS --out ${OUTPUT}_perm

# Sort permutation results
sort -gk 4 ${OUTPUT}_perm.assoc.mperm > ${OUTPUT}_perm_sorted.txt

# Step 5: Generate Manhattan plot
echo "Step 5/6: Generating Manhattan plot"
cp D:/03.projects/AI.PGT/snparray_analysis/src/RScript/Manhattan_plot.R ./
Rscript --no-save Manhattan_plot.R

# Step 6: Generate QQ plot
echo "Step 6/6: Generating QQ plot"
cp D:/03.projects/AI.PGT/snparray_analysis/src/RScript/QQ_plot.R ./
Rscript --no-save QQ_plot.R

echo "GWAS Association Analysis completed successfully!"
echo "Output files:"
echo "- ${OUTPUT}_basic.assoc (Basic association)"
echo "- ${OUTPUT}_logistic.assoc_2.logistic (Logistic regression)"
echo "- ${OUTPUT}_adjusted.assoc.adjust (Multiple testing corrected)"
echo "- ${OUTPUT}_perm_sorted.txt (Permutation results)"
echo "- Manhattan_plot.png and QQ_plot.png (Visualizations)"