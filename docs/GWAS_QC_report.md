# GWAS质量控制技术报告

## 1. 概述
详细记录基于HapMap数据进行GWAS质量控制流程。使用HapMap_3_r3_1数据集，包含经过筛选的欧洲人群（EUR）个体数据，并添加了模拟的二元表型性状

## 2. 数据准备
- 原始数据: HapMap_3_r3_1 (来自HapMap项目phase III)
- 人群筛选: 仅包含欧洲人群个体，避免群体分层问题
- 表型数据: 模拟的二元性状数据，仅包含 founders 个体

## 3. 质量控制流程

### 3.1 缺失率分析
```bash
# 生成个体和SNP的缺失率统计
plink --bfile HapMap_3_r3_1 --missing
# plink.imiss: 个体水平缺失率
# plink.lmiss: SNP水平缺失率
```
- 宽松过滤
```bash
plink -bfile HapMap_3_r3_1 --geno 0.2 --make-bed --out HapMap_3_r3_2
plink --bfile HapMap_3_r3_2 --mind 0.2 --make-bed --out HapMap_3_r3_3
```
- 严格过滤
```bash
# SNP缺失率 > 2% 过滤
plink --bfile HapMap_3_r3_3 --geno 0.02 --make-bed --out HapMap_3_r3_4
# 个体缺失率 > 2% 过滤
plink --bfile HapMap_3_r3_4 --mind 0.02 --make-bed --out HapMap_3_r3_5
```

### 3.2 性别不一致检查
```bash
plink --bfile HapMap_3_r3_5 --check-sex
# 判断标准:
# 女性： F Value < 0.2
# 男性： F Value: > 0.8
```
- 问题个体处理
```bash
# 方案1：直接移除问题个体
grep "PROBLEM" plink.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt
plink --bfile HapMap_3_r3_5 --remove sex_discrepancy.txt --make-bed --out HapMap_3_r3_6

# 方案2： 基于基因型推算性别(注释状态)
plink --bfile HapMap_3_r3_5 --impute-sex --make-bed --out HapMap_3_r3_6
```

### 3.3 常染色体SNP筛选与MAF过滤

- 常染色体SNP选择
```bash
# 选择chr1-22的SNP
awk '{if ($1 >=1 && $1 <= 22) print $1 )}' HapMap_3_r3_6.bim > snp_1_22.txt
plink --bfile HapMap_3_r3_6 --extract snp_1_22.txt --make-bed --out HapMap_3_r3_7

# SNP数量变化：1,430,443 → 1,398,544
```

- MAF分析与过滤
```bash
# MAF 分布分析
plink --bfile HapMap_3_r3_7 --freq --out MAF_check

# MAF 过滤(threshold = 0.05)
plink --bfile HapMap_3_r3_7 --maf 0.05 --make-bed --out HapMap_3_r3_8

# SNP数量变化：1,398,544 → 1,073,226
```

### 3.4 Hardy-Weinberg平衡检验
- HWE分析
```bash
plink --bfile HapMap_3_r3_8 --hardy

# 提取显著偏离HWE的SNP(P < 0.00001)
awk '{if (%9 < 0.00001) print $0}' plink.hwe > plinkzoomhwe.hwe
```

- 分层过滤策略
```bash
# step1：对照组严格过滤
plink --bfile HapMap_3_r3_8 --hwe 1e-6 --make-bed --out HapMap_hwe_filter_ste1

# step2：全部样本宽松过滤
plink --bfile HapMap_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out HapMap_3_r3_9

# SNP Nums: 1,073,226 → 1,070,827
```

### 3.5 杂合率异常检测
- 独立SNP选择
```bash
# 排除高连锁不平衡区域，选择独立SNP
plink --bfile HapMap_3_r3_9 --exclde inversion.txt --range --indep-pairwise 50 5 0.2 -- out indepSNP.prune.in
```

- 杂合率计算
```bash
plink --bfile HapMap_3_r3_9 --extract indepSNP.prune.in --het --out R_check
```

- result
```markdown
Pruning complete.  26236 of 47236 variants removed.

- 总变异位点: 47236
- 修建掉的位点： 26,236 (55%)
- 保留的位点: 21,000 (44.5%)


Make list wirtten to indepSNP.prune.in and indepSNP.prune.out
- indepSNP.prune.in: 保留的SNP列表
- indepSNP.prune.out: 修剪掉的高LD SNP列表

数据概况
- 样本：81人（45男， 34女，2个性别模糊）
- 初始变异位点：47,236
- 数据质量: 总基因分型率99.4%
- 表型: control: 52, case: 27
```
