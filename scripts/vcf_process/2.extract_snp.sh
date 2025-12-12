#!/bin/bash

set -e
set -o pipefail

# This script extracts SNPs from the VCF file and creates a new VCF file with only SNPs.

# Input: VCF file with all variants
# Output: VCF file with only SNPs

# Usage: ./extract_snp.sh input.vcf output.vcf

# Check if input and output files are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input.vcf output.vcf"
  exit 1
fi

# Set input and output file names
VCF_FILE=${1:-input.vcf}
OUTPUT_VCF=${2:-snp_output.vcf}
SNP_POS_BED=${3:-snp_pos.txt}

# create index by tabix
tabix -p vcf $VCF_FILE

# Extract positions of SNPs
bcftools view -R $SNP_POS_BED $VCF_FILE -Oz -o $OUTPUT_VCF
tabix -p vcf $OUTPUT_VCF