# GWAS中群体分层控制的技术报告

## 1. Background
- 群体分层问题
群体分层（population Stratification) 是GWAS中一个重要的混杂因素。当研究样本包含不同遗传背景的亚群体时，如果表型分布在亚群间存在差异，可能导致假阳性关联结果。这种现象源于等位基因频率在自然人群中因祖先差异而系统性地变化。

- 解决方案
PCA和多维标度分析（MDS）是控制群体分层的常用方法。通过将样本投影到遗传背景空间中，可以识别和校正群体结构的影响。本教程使用1000 Genomes项目数据作为参考面板，通过MDS分析检测和校正HapMap数据中的群体分层。

## 2. 方法流程
- 准备数据
```bash
# 从前一教程复制必要文件
cp HOME/{user}/{path/1_QC_GWAS}/HapMap_3_r3_12.* HOME/{user}/{path/2_Population_stratification}
cp HOME/{user}/{path/1_QC_GWAS}/indepSNP.prune.in HOME/{user}/{path/2_Population_stratification}
```

- 1000 Genomes 参考数据下载与处理
```bash
# 下载1000 Genomes数据（约60GB）
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz

# 格式转换和质量控制
plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs

# VCF 格式转化为PLINK二进制格式，提高处理效率
# 为缺失rs标识符的SNP分配唯一ID，确保后续分析顺利进行
```

- 数据质量控制
```bash
# 多步QC流程
plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --make-bed --out 1kG_MDS
plink --bfile 1kG_MDS --mind 0.2 --make-bed --out 1kG_MDS2
plink --bfile 1kG_MDS2 --geno 0.02 --make-bed --out 1kG_MDS3
plink --bfile 1kG_MDS3 --mind 0.02 --make-bed --out 1kG_MDS4
plink --bfile 1kG_MDS4 --maf 0.05 --make-bed --out 1kG_MDS5

# 基因型缺失率：先宽松（0.2）后严格（0.02）
# 个体缺失率：先宽松（0.2）后严格（0.02）
# 次要等位基因频率（MAF）：0.05
```

- SNP交集提取
```bash
# 提供共有SNP
awk '{print$2}' HapMap_3_r3_12.bim > HapMap_SNPs.txt
plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6

awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
plink --bfile HapMap_3_r3_12 --extract 1kG_MDS6_SNPs.txt --make-bed --out HapMap_MDS
```

- 基因组版本统一
```bash
# 统一物理位置信息
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
```

- 链方向问题
```bash
# 参考等位基因统一
awk '{print$2, $5}' 1kG_MDS7.bim > 1kg_ref-list.txt
plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj

# 链方向翻转
awk '{print$2, $5, %6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
sort 1kGMDS7_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap

# 移除无法解决的SNP
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
```

## 2.1 数据合并于MDS分析

- 数据合并
```bash
plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --makebed --out HapMap_MDS2
plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --make-bed --out MDS_merge2
```

-- 多维标度分析
```bash
# 使用预修剪的SNP进行MDS
plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

# MDS原理
# 1. 基于个体间遗传相似性矩阵
# 2. 提取主要维度反映群体结构
# 3. 前几个维度通常对应地理祖先成分
```

## 2.2 群体异常值排除与协变量生成

- 协变量文件生成
```bash
# 在纯化后的数据上重新进行MDS
plink --bfile HapMap_3_r3_13 --extract indepSNP.prune.in --genome --out HapMap_3_r3_13
plink --bfile HapMap_3_r3_13 --read-genome HapMap_3_r3_13.genome --cluster --mds-plot 10 --out HapMap_3_r3_13_mds
```

## 3. 技术总结
1. 参考面板使用： 利用1000 Genomes项目作为遗传背景参考
2. 逐步QC策略： 从宽松到严格的多步质量控制
3. 数据协调：全面的SNP协调流程解决技术差异
4. 可视化验证：通过MDS图直观验证群体结构