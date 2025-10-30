# GWAS统计方法与多重检验校正

## 1. Background

- GWAS概述
GWAS是一种检测遗传变异与表型之间关联的统计方法。通过同时分析数十万~百万SNP，GWAS能够系统性识别与复杂疾病或性状相关的遗传位点。

- 挑战
1. 多重检验问题：同时检验大量SNP导致假阳性率升高
2. 群体分层： 样本群体结构可能造成徐建关联
3. 统计效能：需要足够样本量检测微小效应

## 2. 关联分析

- 基础关联分析（--assoc）
```bash
plink --bfile ${bfile} --assoc --out ${out}

# 技术原理：
# 1. 使用卡方检验或Fisher精确检验
# 2. 计算每个SNP的等位基因频率在病例组和对照组间的差异
# 3. 输出包含p值、优势比（OR）和95%置信区间（CI）的文件
# 优缺点：
# 优点： 计算快速，内存需求低
# 缺点： 无法校正协变量，容易受混杂因素影响
```

- 逻辑回归分析 （--logistic）
```bash
plink --bfile ${bfile} --covar covar_mds.txt --logistic --hide-covar --out logistic_results
```
统计模型
```text
logit(P(Y=1)) = β₀ + β₁SNP + β₂PC1 + ... + βₖPCk + ε
- Y: 二分类标签
- SNP：基因型计数（0, 1, 2）
- PC1...PCk 为前k个组成份（群体分层校正）
- `--covar`: 指定协变量文件
- `--hidecovar`: 输出中隐藏协变量结果，只显示SNP关联
- `--logistic`: 使用逻辑回归模型
```

## 3. 多重校验校正方法
- 传统基因组显著性阈值

Bonferrronni校正
```text
α_corrected = 0.05 / N_snps
```

对于典型GWAS
```text
α_corrected = 5 × 10⁻⁸

- 基于人类基因组连锁不平衡结构
- 约100万个独立检验的保守估计
```

- PLINK多重校验校正（--adjust）
```bash
plink --bfile ${bfile} --assoc --adjust --out ${out}

# 提供校正方法
# Bonferroni： 最保守方法
# Holm：逐步Bonferroni校正
# Sidak： 基于独立检验假设
# FDR (False Discovery Rate)
```

FDR原理
```text
对于排序后的p值p₁ ≤ p₂ ≤ ... ≤ p_m
拒绝所有满足pᵢ ≤ (i/m) × q的假设
其中q为期望的FDR水平
```

- 排列检验
```bash
# 生成SNP子集
awk '{ if ($4 >= 21595000 && $4 <= 21605000) print $2 }' PGT_TLS.refactor_11.bim > subset_snp_chr_22.txt

# 提取子集数据
plink --bfile PGT_TLS.refactor_11 --extract subset_snp_chr_22.txt --make-bed --out HapMap_subset_for_perm

# 执行100万次排列检验
plink --bfile HapMap_subset_for_perm --assoc --mperm 1000000 --out subset_1M_perm_result

# 结果排序
sort -gk 4 subset_1M_perm_result.assoc.mperm > sorted_subset.txt
```

排列检验原理

1. 随机打乱表型标签，破坏真实关联
2. 对每个排列重复GWAS分析
3. 计算经验P值
```text
p_empirical = (r + 1)/(n + 1)
其中r为排列检验中观察到的更极端统计量次数
n为总排列次数
```

## 4. 结果可视化

- Manhattan图
```r
library(qqman)
result <- read.table("logistic_results.assoc_2.logistic", header=TRUE)
resuts <- results[complete.cases(results), ]
png("manhattan.png", with=1000, height=400)
manhattan(results, chr="CHR", bp="BP", p="P", snp="SNP",
    main="Manhattan Plot: GWAS Result",
    suggestiveline=-log10(1e-5),
    genomewideline=-log10(5e-8)
    )
dev.off()

# X: 基因组位置（按染色体排列）
# Y: -log10(p value)
# 建议显著性线性： 通常设为1e-5
# 每个点代表一个SNP的关联强度
```

- QQ图(Quantile-Quantile Plot)
```r
library(qqman)
results <- read.table("logistic_results.assoc_2.logistic", header=TRUE)
results <- results[complete.cases(results), ]
png("qq.png", with=1000, height=400)
qq(results$P, main="QQ Plot: GWAS Result")
dev.off()

# X: 理论分位数
# Y: 实际分位数
# 对角线：无关联情况下的期望分布
# λ = 观察的中位数χ²统计量 / 期望的中位数χ²统计量
```

| 分析方法 | 适用场景 | 优点 | 局限性 |
|---------|----------|------|--------|
| `--assoc` | 初步筛查，计算资源有限 | 快速简单 | 无法校正协变量 |
| `--logistic` | 二分类性状，需要校正 | 控制混杂因素 | 计算较慢 |
| `--linear` | 连续性状 | 处理定量表型 | 假设正态分布 |

## 5. GWAS logistics结果解读

| CHR | SNP | BP | A1 | F_A | F_U | A2 | CHISQ | P | OR |
|-----|-----|----|----|-----|-----|----|-------|---|----|
| 1 | rs3094315 | 752566 | G | 0.2407 | 0.1765 | A | 0.9159 | 0.3385 | 1.48 |
| 1 | rs12092254 | 1113121 | A | 0.03704 | 0.09804 | G | 1.85 | 0.1737 | 0.3538 |

- CHR: 染色体编号
- SNP: SNP标识符
- BP: 物理位置
- A1: 优势等位基因
- F_A: 病例组等位基因频率
- F_U: 总体等位基因频率
- A2: 非优势等位基因
- CHISQ: 卡方统计量
- P: 显著性水平
- OR: 优势比, 效应等位基因与疾病风险的关联强度