# PCA demo
# Clear workspace
rm(list = ls())
# Set working directory
setwd("D:/03.projects/AI.PGT/data/SNP_results/PGT_TLS/PLINK_281025_0434/pca")
# generate bfile
system("plink --file PGT_TLS.refactor --autosome --make-bed --out PGT_TLS.refactor")
# Run PLINK QC
system("plink --file PGT_TLS.refactor --autosome --geno 0.1 --mind 0.1 --maf 0.05 --nonfounders --allow-no-sex --recode --out PGT_TLS.refactor.PCA")

######################################
# Principal Component Analysis - PCA #
######################################

## Genetic distances between individuals
system("plink --allow-no-sex --nonfounders --file PGT_TLS.refactor.PCA --distance-matrix --out PGT_TLS.refactor.dataForPCA")

## Load data
dist_populations<-read.table("PGT_TLS.refactor.dataForPCA.mdist",header=F)
### Extract breed names
fam <- data.frame(famids=read.table("PGT_TLS.refactor.dataForPCA.mdist.id")[,1])
### Extract individual names 
famInd <- data.frame(IID=read.table("PGT_TLS.refactor.dataForPCA.mdist.id")[,2])

## Perform PCA using the cmdscale function 
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dist_populations,eig=T,5)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

# Visualize PCA in tidyverse
# Load tidyverse
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

# PCA plot
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`,color = famids), show.legend = FALSE ) + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of wordwide goat populations",
       x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  theme_minimal()
