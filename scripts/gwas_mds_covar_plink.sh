#!/bin/bash

# Population Stratification Analysis Script
# Usage: ./population_stratification.sh [data_directory] [input_prefix] [output_prefix] [mds_dimensions]

set -e  # Exit on any error

# Default parameters
DATA_DIR=${1:-.}
INPUT=${2:-PGT_TLS.refactor_step10}
OUTPUT=${3:-PGT_TLS.refactor_step10_mds}
MDS_DIM=${4:-10}

cd $DATA_DIR

echo "Starting Population Stratification Analysis"
echo "Input: $INPUT"
echo "Output prefix: $OUTPUT"
echo "MDS dimensions: $MDS_DIM"

# Step 1: Generate genome file using pruned SNPs
echo "Step 1/3: Generating genome file with pruned SNPs"
plink --bfile $INPUT --extract indepSNP.prune.in --genome --out $INPUT

# Step 2: Perform MDS analysis
echo "Step 2/3: Performing MDS analysis with $MDS_DIM dimensions"
plink --bfile $INPUT --read-genome ${INPUT}.genome --cluster --mds-plot $MDS_DIM --out $OUTPUT

# Step 3: Create covariate file
echo "Step 3/3: Creating MDS covariate file"
awk '{print $1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' ${OUTPUT}.mds > covar_mds.txt

echo "Population Stratification Analysis completed successfully!"
echo "Output files:"
echo "- ${OUTPUT}.mds (MDS coordinates)"
echo "- covar_mds.txt (Covariates for association analysis)"
echo "- ${INPUT}.genome (Identity-by-descent matrix)"