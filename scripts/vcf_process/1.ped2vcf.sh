#!/bin/bash

set -e

PED_FILE=${1:-input.ped}
VCF_FILE=${2:-output.vcf}
# step 1: static and filter by MAF and missingness
plink --file ${PED_FILE} --make-bed --out $PED_FILE
 --allow-no-sex --missing
plink --bfile $PED_FILE --geno 0.2 --maf 0.05 --make-bed --out $PED_FILE
# step 2: convert to vcf
plink --bfile $PED_FILE --recode vcf --out $VCF_FILE

# get snp position
plink --bfile $PED_FILE --recode --out $PED_FILE.ped
awk '{print $1, $2, $4}' ${PED_FILE}.ped > ${VCF_FILE}_snp_pos.txt
