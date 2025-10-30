################ 主脚本说明 ##########################

# 本教程使用自由可用的HapMap数据：hapmap3_r3_b36_fwd.consensus.qc。我们模拟了一个二元结果测量（即二元表型性状）并将其添加到数据集中。结果测量仅针对HapMap数据中的创始人进行模拟。该数据集将被称为PGT_TLS.refactor_1。
# 不包含我们模拟结果测量的HapMap数据也可以从 http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2010-05_phaseIII/plink_format/ 获取。
# 要执行本教程，至关重要的一点是，本教程的所有脚本都位于您UNIX工作站的同一目录中。
# 许多脚本包含了解释这些脚本工作原理的注释。注意，要完成本教程，必须执行本教程中的所有命令。
# 此脚本也可用于您自己的数据分析，要这样使用，请将HapMap文件的名称替换为您自己的数据文件的名称。
# 此外，此脚本基于二元结果测量，因此不适用于定量结果测量（这需要进行一些调整）。
# 注意，大多数GWAS研究是在种族同质人群中进行的，其中移除了群体异常值。用于本教程的HapMap数据包含多个不同的种族群体，这给分析带来了问题。
# 因此，我们在教程1-3中仅选择了完整HapMap样本中的EUR个体。这个选择已经在我们GitHub页面上的PGT_TLS.refactor_1文件中完成。
# 本教程中使用的Rscript都是从Unix命令行执行的。
# 因此，只需将“主脚本”中的所有命令复制并粘贴到Unix终端中，即可完成本教程以及我们GitHub页面上的其他教程。
# 关于所有QC步骤的详细理论解释，我们参考了伴随本教程的文章“进行全基因组关联研究的教程：质量控制和统计分析”(https://www.ncbi.nlm.nih.gov/pubmed/29484742)。


##############################################################
############### 开始分析 ###############################
##############################################################

# 将目录更改为您UNIX设备上包含“1_QC_GWAS.zip”中所有文件的文件夹。
# cd D:/03.projects/AI.PGT/data/SNP_results/PGT_TLS/PLINK_281025_0434/1_QC
setwd("D:/03.projects/AI.PGT/data/SNP_results/PGT_TLS/PLINK_281025_0434/1_QC")

### 步骤 1 ###

# 调查每个个体和每个SNP的缺失情况，并绘制直方图。
system("plink --file PGT_TLS.refactor --impute-sex --recode --out PGT_TLS.refactor_1")
system("plink --file PGT_TLS.refactor_1 --autosome --make-bed --out PGT_TLS.refactor_1")
system("plink --bfile PGT_TLS.refactor_1 --missing")
# 输出：plink.imiss 和 plink.lmiss，这些文件分别显示每个个体的缺失SNP比例和每个SNP的缺失个体比例。


# 生成图表以可视化缺失结果。
# Rscript --no-save hist_miss.R

# 删除缺失率高的SNP和个体，对此步骤及后续所有步骤的解释可以在本脚本注释中提到的文章的框1和表1中找到。
# 以下两个QC命令不会删除任何SNP或个体。然而，以这些非严格阈值开始QC是良好的做法。
# 删除缺失率 >0.2 的SNP。
system("plink --bfile PGT_TLS.refactor_1 --geno 0.2 --make-bed --out PGT_TLS.refactor_2")

# 删除缺失率 >0.2 的个体。
system("plink --bfile PGT_TLS.refactor_2 --mind 0.2 --make-bed --out PGT_TLS.refactor_3")

# 删除缺失率 >0.02 的SNP。
system("plink --bfile PGT_TLS.refactor_3 --geno 0.02 --make-bed --out PGT_TLS.refactor_4")

# 删除缺失率 >0.02 的个体。
system("plink --bfile PGT_TLS.refactor_4 --mind 0.02 --make-bed --out PGT_TLS.refactor_5")

###################################################################
### 步骤2 ####

# 检查性别差异。
# 事先确定为女性的受试者的F值必须<0.2，而事先确定为男性的受试者的F值必须>0.8。此F值基于X染色体近交（纯合性）估计。
# 不满足这些要求的受试者将被PLINK标记为“PROBLEM”。

# system("plink --bfile PGT_TLS.refactor_5 --check-sex")

# 生成图表以可视化性别检查结果。
Rscript --no-save gender_check.R
# 这些检查表明，有一名女性存在性别差异，F值为0.99。（当使用其他数据集时，通常会发现一些差异）。

# 以下两个脚本可用于处理存在性别差异的个体。
# 注意，请使用以下两个选项之一来生成bfile hapmap_r23a_6，我们将在本教程的下一步中使用此文件。

# 1) 删除存在性别差异的个体。
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
# 此命令生成一个状态为“PROBLEM”的个体列表。
plink --bfile PGT_TLS.refactor_5 --remove sex_discrepancy.txt --make-bed --out PGT_TLS.refactor_6
# 此命令移除状态为“PROBLEM”的个体列表。

# 2) 推算性别。
system("plink --bfile PGT_TLS.refactor_1 --impute-sex --make-bed --out PGT_TLS.refactor_1")
# 这会根据基因型信息将性别推算到您的数据集中。

###################################################
### 步骤 3 ###

# 生成一个仅包含常染色体SNP的bfile，并删除次要等位基因频率（MAF）低的SNP。

# 仅选择常染色体SNP（即染色体1到22）。 SNP数量: 1430443 -> 1398544
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' PGT_TLS.refactor_6.bim > snp_1_22.txt
system("plink --bfile PGT_TLS.refactor_5 --extract snp_1_22.txt --make-bed --out PGT_TLS.refactor_7")

# 生成MAF分布图。
system("plink --bfile PGT_TLS.refactor_7 --freq --out MAF_check")
Rscript --no-save MAF_check.R

# 移除MAF频率低的SNP。 SNP数量: 1398544 -> 1073226
system("plink --bfile PGT_TLS.refactor_7 --maf 0.05 --make-bed --out PGT_TLS.refactor_8")
# 剩下 1073226 个SNP
# 常规GWAS的传统MAF阈值在0.01或0.05之间，具体取决于样本量。


####################################################
### 步骤 4 ###

# 删除不符合哈迪-温伯格平衡（HWE）的SNP。 1073226 -> 1072827
# 检查所有SNP的HWE p值分布。

system("plink --bfile PGT_TLS.refactor_8 --hardy")
# 选择HWE p值低于0.00001的SNP，这是下一个Rscript生成的两个图之一所需要的，可以放大查看严重偏离的SNP。
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
Rscript --no-save hwe.R

# 默认情况下，plink中的 --hwe 选项仅针对对照组进行过滤。
# 因此，我们分两步进行，首先对对照组使用严格的HWE阈值，然后对病例数据使用较宽松的阈值。
system("plink --bfile PGT_TLS.refactor_8 --hwe 1e-6 --make-bed --out HapMap_hwe_filter_step1")

# 病例的HWE阈值仅过滤掉那些极度偏离HWE的SNP。
# 第二个HWE步骤仅针对病例，因为在对照组中所有HWE p值 < 1e-6 的SNP已经被移除。
system("plink --bfile HapMap_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out PGT_TLS.refactor_9")

# 此步骤的理论背景在我们随附的文章中给出：https://www.ncbi.nlm.nih.gov/pubmed/29484742。

############################################################
### 步骤 5 ###

# 生成您受试者杂合率分布的图表。
# 并移除杂合率偏离均值超过3个标准差的个体。

# 杂合性检查是在一组不高度相关的SNP上进行的。
# 因此，为了生成非（高度）相关SNP的列表，我们排除了高倒位区域（inversion.txt [高LD区域]）并使用命令 --indep-pairwise 来修剪SNP。
# 参数“50 5 0.2”分别代表：窗口大小、每一步移动窗口的SNP数量、以及一个SNP同时对所有其他SNP进行回归的多元相关系数。

system("plink --bfile PGT_TLS.refactor_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP")
# 注意，不要删除文件 indepSNP.prune.in，我们将在本教程的后续步骤中使用此文件。

system("plink --bfile PGT_TLS.refactor_9 --extract indepSNP.prune.in --het --out R_check")
# 此文件包含您修剪后的数据集。

# 绘制杂合率分布图
Rscript --no-save check_heterozygosity_rate.R

# 以下代码生成一个列表，包含那些杂合率偏离均值超过3个标准差的个体。
# 对于数据操作，我们建议使用UNIX。然而，在执行统计计算时，R可能更方便，因此此步骤使用Rscript：
Rscript --no-save heterozygosity_outliers_list.R

# 上述命令的输出：fail-het-qc.txt。
# 当使用我们的示例数据/HapMap数据时，此列表包含2个个体（即，有两个个体的杂合率偏离均值超过3个SD）。
# 通过从文件中移除所有引号并仅选择前两列，使此文件与PLINK兼容。
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# 移除杂合率异常值。
system("plink --bfile PGT_TLS.refactor_9 --remove het_fail_ind.txt --make-bed --out PGT_TLS.refactor_10")


############################################################
### 步骤 6 ###

# 检查您分析的数据集是否存在隐性亲缘关系至关重要。
# 假设是一个随机群体样本，我们将在本教程中排除所有pihat阈值高于0.2的个体。

# 检查pihat > 0.2的个体之间的关系。
system("plink --bfile PGT_TLS.refactor_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2")

# 已知HapMap数据集包含亲子关系。
# 以下命令将使用z值专门可视化这些亲子关系。
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome

# 生成图表以评估关系类型。
Rscript --no-save Relatedness.R

# 生成的图表显示Hapmap数据中存在大量相关个体（图例说明：PO = 亲子，UN = 无关个体），这是预期的，因为数据集就是这样构建的。
# 通常，基于家族的数据应使用特定的家族方法进行分析。在本教程中，出于演示目的，我们将这种相关性视为随机群体样本中的隐性相关性。
# 在本教程中，我们的目标是从数据集中移除所有“相关性”。
# 为了证明大部分相关性是由于亲子关系造成的，我们只纳入创始人（数据集中没有父母的个体）。

system("plink --bfile PGT_TLS.refactor_10 --filter-founders --make-bed --out PGT_TLS.refactor_11")

# 现在我们再次寻找pihat >0.2的个体。
system("plink --bfile PGT_TLS.refactor_11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders")
# 文件 'pihat_min0.2_in_founders.genome' 显示，在排除所有非创始人之后，HapMap数据中只剩下1对pihat大于0.2的个体对。
# 根据Z值，这很可能是一对全同胞或异卵双胞胎。值得注意的是，在HapMap数据中，他们没有获得相同的家族标识（FID）。

# 对于每对pihat > 0.2的“相关”个体，我们建议删除调用率较低的那个个体。 5900 <- 5995
system("plink --bfile PGT_TLS.refactor_11 --missing")
# 使用UNIX文本编辑器（例如，vi(m)）检查“相关对”中哪个个体的调用率最高。

# 生成一个FID和IID的列表，包含pihat高于0.2的个体，以检查该对中谁的调用率较低。
# 在我们的数据集中，个体 13291  NA07045 的调用率较低。
vi 0.2_low_call_rate_pihat.txt
i
13291  NA07045
# 在键盘上按 esc!
:x
# 在键盘上按回车
# 如果有多个“相关”对，可以使用与我们单个“相关”对相同的方法扩展上面生成的列表。

# 删除在pihat > 0.2的“相关”对中调用率较低的个体
plink --bfile PGT_TLS.refactor_11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out PGT_TLS.refactor_12

################################################################################################################################

# 恭喜！！您刚刚成功完成了第一个教程！您现在能够进行正确的遗传QC。

# 对于下一个教程，使用脚本：2_Main_script_MDS.txt，您需要以下文件：
# - bfile PGT_TLS.refactor_12 (即 PGT_TLS.refactor_12.fam, PGT_TLS.refactor_12.bed, 和 PGT_TLS.refactor_12.bim)
# - indepSNP.prune.in