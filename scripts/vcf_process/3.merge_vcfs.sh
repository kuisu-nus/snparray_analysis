#!/bin/bash

set -e  # Exit on any error

# vcf files to vcf.gz by bgzip
bcftools view myfile.vcf -Oz -o myfile.vcf.gz

# Index the compressed VCF files
bcftools index myfile.vcf.gz
tabix -p vcf myfile.vcf.gz

# merge vcf.gz files by same chromosome
bcftools merge chr1_part1.vcf.gz chr1_part2.vcf.gz -Oz -o chr1_merged.vcf.gz

# 从vcf.gz中提取坐标
bcftools query -f '%CHROM\t%POS0\t%POS\n' Filter.ALL.genome.vcf.gz > union_snp_pos.bed