# This R script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
# This R script uses multivariate-SNP, multivariate-PT (eGFR and BUN) metaCCA.
# metaCCA reference:
# Cichonska A, Rousu J, Marttinen P, Kangas AJ, Soininen P, Lehtimäki T, Raitakari OT, Järvelin MR, Salomaa V, Ala-Korpela M, Ripatti S, Pirinen M.
# metaCCA: summary statistics-based multivariate meta-analysis of genome-wide association studies using canonical correlation analysis.
# Bioinformatics. 2016 Jul 1;32(13):1981-9. doi: 10.1093/bioinformatics/btw052. Epub 2016 Feb 19. PMID: 27153689; PMCID: PMC4920109.

library("dplyr")
library("metaCCA")
library("e1071")
library("rlist")
library("stringr")
`%notin%` <- Negate(`%in%`)

# The following input files can be created by using the script: data_preparation_summary_statistics.sh
# Import XY table:
XY <- "bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes.txt"
# Import YY table (see univariate-SNP code for how this was derived)
YY <- "bbj-a-60.eGFR.bbj-a-11.BUN.gnomADg_AF_EAS.0.01_ld02_dirfilt_ALL_YY.txt"

data_CKD <- read.table(XY, header=F,skip=1)
Dataset_toRun <- data_CKD[,c(3:9)]
names(Dataset_toRun) <- c("RSID", "allele_0", "allele_1","eGFR_b","eGFR_se","BUN_b","BUN_se")
Dataset_toRun$eGFR_b <- as.numeric(Dataset_toRun$eGFR_b)
Dataset_toRun$BUN_b <- as.numeric(Dataset_toRun$BUN_b)

# Keep only SNPs with eGFR and BUN effect sizes in opposite directions:
Dataset_toRun_dir <- Dataset_toRun %>% filter((Dataset_toRun$eGFR_b<0 & Dataset_toRun$BUN_b>0) | (Dataset_toRun$eGFR_b>0 & Dataset_toRun$BUN_b<0))
nrow(Dataset_toRun_dir)

lookup_table <- data_CKD[,c(1:5,10)]
# Note that the 1000GP ref and alt alleles in XX are still the same as for the original 1000GP dataset, but
# for those that didn't match the study-data ref and alt alleles, the GTs have been "flipped" accordingly -
# see script: data_preparaion_summary_statistics.sh, lines 130 - 186
names(lookup_table) <- c("chr","pos","RSID","ref","alt","gene")
Dataset_toRun_dir$RSID <- as.character(Dataset_toRun_dir$RSID)

# Filter out any SNPs with GWAS summary statistics standard error (se)=0
Dataset_toRun_dir <- Dataset_toRun_dir[(Dataset_toRun_dir$eGFR_se!=0),]
Dataset_toRun_dir <- Dataset_toRun_dir[(Dataset_toRun_dir$BUN_se!=0),]

lookup_table$RSID <- as.character(lookup_table$RSID)
lookup_table$gene <- as.character(lookup_table$gene)
lookup_table <- lookup_table[lookup_table$RSID %in% Dataset_toRun_dir$RSID,]


# Run for each chromosome at a time:
chrom <- sort(unique(as.numeric(lookup_table$chr)))
for (j in chrom){
	print("Chrom: ")
	print(j)
	rsids_chrom <- as.character(lookup_table[lookup_table$chr==j,3])
	print(length(rsids_chrom))
	GT_filename <- paste0("~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr", j, ".txt")
	#GT_filename <- "~/uklhs/g1000p3_Prins_egfr_ure_GTs_sorted_LD02_chr1_glist_hg19_forXX.txt"
	print(GT_filename)
	GT_data2 <- read.table(GT_filename,header=F,skip=1)
	GT_data2$V3 <- as.character(GT_data2$V3)
	GT_data2 <- GT_data2[GT_data2$V3 %in% Dataset_toRun_dir$RSID,]
	GTs <- GT_data2[,c(3,10:513)]
	GTs[,1] <- as.character(GTs[,1])
	# Get rid of any duplicated SNPs:
	GTs <- GTs[!duplicated(GTs$V3),]
	rownames(GTs) <- GTs[,1]
	GTs <- GTs[ order(row.names(GTs)), ]
	GTs_t <- t(GTs[,2:504])
	GTs_t_df <- as.data.frame(GTs_t)
	# Get rid of any SNPs with reference genotype dataset (1000GP) standard deviation of 0
	getridofsd0 <- Filter(function(x) sd(x) != 0, GTs_t_df)  # ncol
	print(ncol(getridofsd0))
	S_XX_study1 <- cor(getridofsd0)
	copy_SXX <- cor(getridofsd0)
	# Any zeroes in XX genotype-genotype correlation matrix? Get rid of these zeroes:
	copy_SXX[lower.tri(copy_SXX)] <- NA
	zeroes <- rownames(which(copy_SXX==0, arr.ind=TRUE))
	S_XX_study1 <- S_XX_study1[rownames(S_XX_study1) %notin% zeroes,colnames(S_XX_study1) %notin% zeroes]
	snps_correct <- rownames(S_XX_study1)
	lookup_table2 <- lookup_table[lookup_table$RSID %in% snps_correct,]
	genes <- unique(lookup_table2$gene)
	l.genes <- length(genes)
	snps <- lookup_table2$RSID
	l.snps <- length(snps)
	# Save SNPs run to file:
	if(j==1){
		write.table(lookup_table2, file="bbj-a-60.eGFR.bbj-a-11.BUN.gnomADg_AF_EAS.0.01_ld02_dirfilt_lookup_table_run.txt", quote=F,col.names=T,row.names=F)
	}
	else{
		write.table(lookup_table2, file="bbj-a-60.eGFR.bbj-a-11.BUN.gnomADg_AF_EAS.0.01_ld02_dirfilt_lookup_table_run.txt", quote=F,append=T,row.names=F)
	}

	# Save information on how many SNPs, genes run per chromosome to file:
	phentype.l = 2
	params <- c("chr:", j, "snps: ", l.snps, "genes: ", l.genes)
	print(params)
	write.table(params, "bbj-a-60.eGFR.bbj-a-11.BUN.params.txt", col.names=F, row.names=F,quote=F,append=T)

	# Prepare new S_XY table - filter data_CKD based on RSID provided by snps variable (just been through all data filters)
	S_XY <- Dataset_toRun_dir %>% filter(Dataset_toRun_dir[,1] %in% snps)
	# Get rid of any duplicates:
	S_XY <- S_XY[!duplicated(S_XY$RSID),]
	S_XY$allele_0 <- toupper(S_XY$allele_0)
	S_XY$allele_1 <- toupper(S_XY$allele_1)
	S_XY$RSID <- as.character(S_XY$RSID)

	rownames(S_XY) <- S_XY[,1]
	S_XY_study1 = S_XY[,2:7]
	S_XY_study1$allele_1 <- as.factor(S_XY_study1$allele_1)
	S_XY_study1$allele_0 <- as.factor(S_XY_study1$allele_0)

	# Import the YY (PT-PT) correlation coefficients (precalculated based on all SNPs available):
	S_YY_study1 = read.table(YY,header=T)
	phenotypes=rownames(S_YY_study1)

	# Set other parameters needed:
	N1 = 143658
	std_info = c(0)

	nr_studies = 1
	results_v3 <- data.frame(no_SNPs = numeric(), SNPs = character(), no_PTs = numeric(), PTs = character(), r_1 = numeric(), pval_minuslog10=numeric(), pval=numeric(), genes=character(), chrom=numeric(), stringsAsFactors=FALSE)

	# for each gene in list, get all SNPs, and run metaCCA
	for (i in 1:l.genes){
		print(i)
		cat(genes[i])
		snps_for_gene <- lookup_table2[lookup_table2$gene %in% genes[i],]
		if(nrow(snps_for_gene)>0){
			snpsMetaCCA <- unique(as.character(snps[snps %in% snps_for_gene$RSID]))
			nsnps = length(snpsMetaCCA)
			ntraits = ncol(S_YY_study1)
			nind = N1
			string_snpsMetaCCA = paste(snpsMetaCCA, collapse="_")
			rowno <- nrow(results_v3)+1
			results_v3[rowno,1] <- nsnps
			results_v3[rowno,2] <- string_snpsMetaCCA
			results_v3[rowno,3] <- ntraits
			results_v3[rowno,4] <- "eGFR_BUN"
			if(nsnps==1){
				analysis_type=1
				}
			else{
				analysis_type = 2
			}

			metaCCA_res5 = metaCcaGp(nr_studies = nr_studies, S_XY = list(S_XY_study1), std_info = std_info,S_YY = list(S_YY_study1),
					N = c(N1), analysis_type = analysis_type, SNP_id = snpsMetaCCA, S_XX = list(S_XX_study1))

			names(metaCCA_res5) <- c("r_1","pval_minuslog10")
			# metaCCA - get r_1 (leading canonical correlation value)
			all.can = as.numeric(metaCCA_res5$r_1)
			all.can[all.can>0.99]=0.99
			results_v3[rowno,5] <- metaCCA_res5$r_1
			results_v3[rowno,6] <- metaCCA_res5$pval_minuslog10
			results_v3[rowno,7] <- 10^(-metaCCA_res5$pval_minuslog10)
			results_v3[rowno,8] <- paste(genes[i],collapse="_")
			results_v3[rowno,9] <- j
		}
	}

	results_v3 <- results_v3[ order(results_v3$pval), ]
	if(j==1){
		write.table(results_v3, "bbj-a-60.eGFR.bbj-a-11.BUN.LD_0.01_multi_results.txt", col.names=TRUE,quote=F,row.names=F)
	}
	else{
	write.table(results_v3, file="bbj-a-60.eGFR.bbj-a-11.BUN.LD_0.01_multi_results.txt", quote=F,append=T,row.names=F)
	}
}

# Read all results (saved in file) back in to R:
results_all <- read.table("bbj-a-60.eGFR.bbj-a-11.BUN.LD_0.01_multi_results.txt", header=T)

# Calculate Bonferroni-adjusted P-values:
results_all$pval_num <- as.numeric(levels(results_all$pval))[results_all$pval]
results_all$pBonf <- p.adjust(results_all$pval_num,method="bonferroni")
# Save results again, now with pBonf:
write.table(results_all,"bbj-a-60.eGFR.bbj-a-11.BUN.LD_0.01_multi_results_pBonf.txt", col.names=T,quote=F,row.names=F)

# Find and save results for significant genes only:
sigBonf <- results_all[results_all$pBonf<0.05,]
nrow(sigBonf)
write.table(sigBonf,"bbj-a-60.eGFR.bbj-a-11.BUN.LD_0.01_multi_results_pBonf_sig.txt", col.names=T,quote=F,row.names=F)
