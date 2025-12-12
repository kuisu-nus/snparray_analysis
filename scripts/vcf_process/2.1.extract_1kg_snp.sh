#!/bin/bash

set -e
set -o pipefail

# This script extracts SNPs from the VCF file and creates a new VCF file with only SNPs.

# Input: VCF file with all variants
# Output: VCF file with only SNPs

# Usage: ./extract_snp.sh input.vcf output.vcf

# Check if input and output files are provided


CHR5=ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR1=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR15=ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR6=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR13=ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR22=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR18=ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR17=ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHRY=ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
CHR3=ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR2=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR16=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR19=ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR11=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR20=ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR10=ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHRMT=ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
CHR12=ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR21=ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR4=ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR8=ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR7=ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHRX=ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
CHR14=ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
CHR9=ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz


OneKG_ROOT=/eds-data/home/kuisu/01.data/01.gene/1kg
OUT_DIR=/eds-data/home/kuisu/01.data/01.gene/1kg_30w_snp

# Extract positions of SNPs
SNP_POS_BED=$OUT_DIR/query_position.bed

bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR1  -Oz -o $OUT_DIR/Filter_$CHR1
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR2  -Oz -o $OUT_DIR/Filter_$CHR2
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR3  -Oz -o $OUT_DIR/Filter_$CHR3
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR4  -Oz -o $OUT_DIR/Filter_$CHR4
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR5  -Oz -o $OUT_DIR/Filter_$CHR5
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR6  -Oz -o $OUT_DIR/Filter_$CHR6
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR7  -Oz -o $OUT_DIR/Filter_$CHR7
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR8  -Oz -o $OUT_DIR/Filter_$CHR8
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR9  -Oz -o $OUT_DIR/Filter_$CHR9
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR10  -Oz -o $OUT_DIR/Filter_$CHR10
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR11  -Oz -o $OUT_DIR/Filter_$CHR11
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR12  -Oz -o $OUT_DIR/Filter_$CHR12
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR13  -Oz -o $OUT_DIR/Filter_$CHR13
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR14  -Oz -o $OUT_DIR/Filter_$CHR14
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR15  -Oz -o $OUT_DIR/Filter_$CHR15
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR16  -Oz -o $OUT_DIR/Filter_$CHR16
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR17  -Oz -o $OUT_DIR/Filter_$CHR17
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR18  -Oz -o $OUT_DIR/Filter_$CHR18
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR19  -Oz -o $OUT_DIR/Filter_$CHR19
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR20  -Oz -o $OUT_DIR/Filter_$CHR20
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR21  -Oz -o $OUT_DIR/Filter_$CHR21
bcftools view -R $SNP_POS_BED $OneKG_ROOT/$CHR22  -Oz -o $OUT_DIR/Filter_$CHR22


bcftools concat Filter_$CHR1 \
    Filter_$CHR2 \
    Filter_$CHR3 \
    Filter_$CHR4 \
    Filter_$CHR5 \
    Filter_$CHR6 \
    Filter_$CHR7 \
    Filter_$CHR8 \
    Filter_$CHR9 \
    Filter_$CHR10 \
    Filter_$CHR11 \
    Filter_$CHR12 \
    Filter_$CHR13 \
    Filter_$CHR14 \
    Filter_$CHR15 \
    Filter_$CHR16 \
    Filter_$CHR17 \
    Filter_$CHR18 \
    Filter_$CHR19 \
    Filter_$CHR20 \
    Filter_$CHR21 \
    Filter_$CHR22 \
    -Oz -o genome.vcf.gz