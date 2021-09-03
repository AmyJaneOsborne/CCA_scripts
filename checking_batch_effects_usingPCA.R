# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk
# checking for batch effects in output from plink: genotype call rate for each sample, and per variant

library(ggplot2)

imiss_sample <- read.table("Genotyping.missing-stat.imiss",header=T)
names(imiss_sample) <- c("FID","sampleID","miss_pheno","n_miss","n_geno","f_miss")

ckd_phen <- read.csv("CKD_data.csv",header=T)

# merge on sampleID to get BatchName
merge <- merge(imiss_sample,ckd_phen, by=c("sampleID"),all.x=T)
merge$sampleno <- seq.int(nrow(merge))
merge[merge$f_miss>=0.015,]

ggplot(data = merge) +
  geom_point(mapping = aes(x = sampleID, y = 100*(f_miss), color = BatchName))+
  labs(title="Missing SNP genotype call rate per sample \n(xx LD-pruned common variants; xx unrelated samples; missing call rate <=1.5%)",
        x ="Sample", y = "Missing SNP genotype call rate (%)")+
  geom_hline(yintercept=1.5, linetype="dashed", color = "black", size = 1.4)+
  scale_y_continuous(breaks=seq(0,12,1))+
  theme(
    panel.grid.major.y = element_line(colour = "grey", size = 1),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_blank()
  )
# save plot manually as PDF.


# do PCA and check for batch effects
pca_table <- read.table("pca_results_king.eigenvec",header=T,comment.char="")
names(pca_table) <- c("fid","sampleID","PC1","PC2","PC3","PC4","PC5")

# merge on sampleID to get BatchName
merge_pca <- merge(pca_table,ckd_phen, by=c("sampleID"),all.x=T)

plot(pca_table[, c("PC1", "PC2", "PC3", "PC4", "PC5")])
ggplot(data = merge_pca) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = BatchName))+
  labs(title="PCA-batch effects check \n (xx LD-pruned common variants; xx unrelated samples; missing call rate <=1.5%)",
        x ="Genotype-based global PC1", y = "Genotype-based global PC2")+
  theme(
    panel.grid.major.y = element_line(colour = "grey", size = 1),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
  )
# save plot manually as PDF.
