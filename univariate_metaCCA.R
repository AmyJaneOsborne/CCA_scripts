# This R script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
# This R script uses univariate-SNP, multivariate-PT (eGFR and BUN) metaCCA.
# metaCCA reference:
# Cichonska A, Rousu J, Marttinen P, Kangas AJ, Soininen P, Lehtimäki T, Raitakari OT, Järvelin MR, Salomaa V, Ala-Korpela M, Ripatti S, Pirinen M.
# metaCCA: summary statistics-based multivariate meta-analysis of genome-wide association studies using canonical correlation analysis.
# Bioinformatics. 2016 Jul 1;32(13):1981-9. doi: 10.1093/bioinformatics/btw052. Epub 2016 Feb 19. PMID: 27153689; PMCID: PMC4920109.

library("metaCCA")
library("dplyr")
library("sqldf")
library("rlist")
library("tidyr")
library(stringr)
`%notin%` <- Negate(`%in%`)

# The following input files can be created by using the script: data_preparation_summary_statistics.sh
Dataset <- read.table("bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt.txt",header=F)
Dataset_toRun <- Dataset[,c(3:9)]
names(Dataset_toRun) <- c("RSID", "allele_0", "allele_1","eGFR_b","eGFR_se","BUN_b","BUN_se")
data_CKD_wgenes <- read.table("bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes_multi.txt", header=T)
data_CKD_wgenes_merge <- merge(data_CKD_wgenes,Dataset_toRun,by=c("RSID"),all.x=T)
Dataset_toRun <- data_CKD_wgenes_merge[,c(1,7:12)]
Dataset_toRun <- Dataset_toRun %>% filter((Dataset_toRun$eGFR_b<0 & Dataset_toRun$BUN_b>0) | (Dataset_toRun$eGFR_b>0 & Dataset_toRun$BUN_b<0))

S_XY = distinct(Dataset_toRun)
S_XY$allele_0 <- toupper(S_XY$allele_0)
S_XY$allele_1 <- toupper(S_XY$allele_1)
S_XY$BUN_b <- as.numeric(S_XY$BUN_b)
S_XY$BUN_se <- as.numeric(S_XY$BUN_se)
nodups <- S_XY$RSID[!duplicated(S_XY$RSID)]
S_XY <- S_XY[!grepl("<|D|E|L|>|N|I", S_XY$allele_0),]
S_XY <- S_XY[!grepl("<|D|E|L|>|N|I", S_XY$allele_1),]

dups <- S_XY$RSID[duplicated(S_XY$RSID)]  # 0
nodups <- S_XY[!duplicated(S_XY$RSID),] #  121163
rownames(nodups) <- nodups[,1]

nodups$eGFR_b <- as.numeric(nodups$eGFR_b)
nodups$eGFR_se <- as.numeric(nodups$eGFR_se)
NAs <- nodups[is.na(nodups$eGFR_se),]
nodups_noNAs <- nodups[!is.na(nodups$eGFR_se),]

S_XY_study1 = nodups_noNAs[,c(2:7)]
S_XY_study1[, 1] <- as.factor(S_XY_study1[, 1])
S_XY_study1[, 2] <- as.factor(S_XY_study1[, 2])

# Estimate YY based on ALL SNPs in dataset:
data_for_YY <- read.table("bbj-a-60.eGFR.bbj-a-11.BUN.noheader.vcf")
data_for_YY_toRun <- data_for_YY[,c(3:9)]
names(data_for_YY_toRun) <- c("RSID", "allele_0", "allele_1","eGFR_b","eGFR_se","BUN_b","BUN_se")
data_for_YY_toRun_nodups <- data_for_YY_toRun[!duplicated(data_for_YY_toRun$RSID),]
rownames(data_for_YY_toRun_nodups) <- data_for_YY_toRun_nodups[,1]
YY_nodups_ready <- data_for_YY_toRun_nodups[,c(2:7)]
S_YY_study1 = estimateSyy(YY_nodups_ready)

# save for future use:
#write.table(S_YY_study1,"~/mrc_ieu_bbj/metaCCA/bbj-a-60.eGFR.bbj-a-11.BUN.gnomADg_AF_EAS.0.01_ld02_dirfilt_ALL_YY.txt",quote=F,row.names=T,col.names=T)
# Load for repeat if needed:
#S_YY_study1 = read.table("~/mrc_ieu_bbj/metaCCA/bbj-a-60.eGFR.bbj-a-11.BUN.gnomADg_AF_EAS.0.01_ld02_dirfilt_ALL_YY.txt",header=T)

# Filter out any SNPs with GWAS summary statistics standard error (se)=0
S_XY_study1_2 <- S_XY_study1[(S_XY_study1$eGFR_se!=0),]
S_XY_study1_3 <- S_XY_study1_2[(S_XY_study1_2$BUN_se!=0),]

# Sample size taken from: https://gwas.mrcieu.ac.uk/datasets/bbj-a-60/ (bigger than for bbj-a-11)
N1 = 143658
metaCCA_res = metaCcaGp(nr_studies =1, S_XY = list(S_XY_study1_3), std_info = c(0), S_YY = list(S_YY_study1),N = c(N1))

# Name columns of output:
names(metaCCA_res) <- c("r_1", "log10_pval")
# Derive P-value:
metaCCA_res$pval <- 10^(-metaCCA_res$log10_pval)
# Order results:
order_metaCCA_res <- metaCCA_res[order(-metaCCA_res$log10_pval),]
# Set RSID as new column (instead of row name)
order_metaCCA_res$RSID <- rownames(order_metaCCA_res)
# Apply Benjamini-Hochberg criterion to determine adjusted P-values based on ranking:
order_metaCCA_res$pBH <- p.adjust(order_metaCCA_res$pval,method="BH")
# Save results file:
write.table(order_metaCCA_res,"metaCCA_univariate_results.txt",row.names=T,quote=F,col.names=T)

# Check how many/which SNPs passed significance cut-off of P-value adjusted < 0.05:
sig_BH2 <- order_metaCCA_res[order_metaCCA_res$pBH<0.05,]

# Add gene information back in:
merged_chrpos <- merge(order_metaCCA_res,data_CKD_wgenes,by=c("RSID"),all.x=T)
sig_BH2_genes <- merged_chrpos[merged_chrpos$pBH<0.05,]
# Save results file with gene names:
write.table(merged_chrpos,"metaCCA_univariate_results_with_chrpos.txt",row.names=T,quote=F,col.names=T)

# Check how many genes were significant:
length(sort(unique(sig_BH2$gene)))
