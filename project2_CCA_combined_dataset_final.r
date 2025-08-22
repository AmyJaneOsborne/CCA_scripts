source("gene-based.cca3.2.R")

library(sqldf)
library(dplyr)
library(GWASTools)

`%notin%` <- Negate(`%in%`)

# F-test
  fapp = function(nsnps,ntraits,nind,all.can)
  {
    wilks = prod(1-all.can^2)
    w = nind - 1 - 0.5*(nsnps+ntraits+1)
    t = sqrt( (nsnps^2 * ntraits^2 - 4) / (nsnps^2 + ntraits^2 - 5) )
    if (nsnps*ntraits == 2) t=1
    df1 = nsnps * ntraits
    df2 = w*t - 0.5*nsnps*ntraits + 1
    f = ( (1-wilks^(1/t)) / (wilks^(1/t)) ) * df2/df1
    p = pf(f,df1,df2,lower.tail=F)
    p
  }

cca <- function(x) {
		value <- as.numeric(cancor(x,phenD_torun)$cor)[1]
		return(value)
	}

fapp_func <- function(x) {
		p <- fapp(nsnps,ntraits,nind,x)
		return(p)
	}
	

# load up phenD data
phenD <- read.csv("NCKD_CKD_PTs_top11_binary.csv", header=T)
phenD <- phenD[!is.na(phenD$SentrixPosition),]
row.names(phenD) <- phenD$SentrixPosition

# To make analyses the same as SKS dataset, exclude ""Microscopic polyangiitis - histologically proven" 
# and ""Primary reflux nephropathy - sporadic" from NCKD analysis
phenD <- phenD[phenD$primary_diagnosis %notin% c("Microscopic polyangiitis - histologically proven","Primary reflux nephropathy - sporadic"),]
# Try MN vs all other only:
phenD <- phenD[order(phenD$VCF_CKD_order),c(14),drop=F]
phenD$ds <- "nckd"
names(phenD) <- c("MN_status","ds")

# load up phenD data
phenD_SKS <- read.csv("SKS_CKD_PTs_binary.csv", header=T)
phenD_SKS <- phenD_SKS[!is.na(phenD_SKS$SentrixPosition),]
#phenD[duplicated(phenD$patient_id),]
row.names(phenD_SKS) <- phenD_SKS$SentrixPosition
phenD_SKS <- phenD_SKS[order(phenD_SKS$VCF_order_no),c(12),drop=F]
phenD_SKS$ds <- "sks"
names(phenD_SKS) <- c("MN_status","ds")

phenD.all <- rbind(phenD, phenD_SKS)

# add in ctrls:
ctrl.nckd <- read.table("GWAS_NURTURE_CONTROLS.txt", header=T)
ctrl.nckd <- ctrl.nckd[-1,]	
colnames(ctrl.nckd) <- c(colnames(ctrl.nckd[,c(1:9)]), gsub("^.{0,108}", "", colnames(ctrl.nckd[,c(10:ncol(ctrl.nckd))])))
ctrls_for_phenD <- colnames(ctrl.nckd[,c(10:ncol(ctrl.nckd))])
phenD.ctrls <- as.data.frame(ctrls_for_phenD)
phenD.ctrls$MN_status <- 0
phenD.ctrls$ds <- "NURTURE_ctrl"
rownames(phenD.ctrls) <- phenD.ctrls$ctrls_for_phenD
phenD.ctrls <- phenD.ctrls[,c(2:3)]
phenD.all <- rbind(phenD.all, phenD.ctrls)

mat_fun <- function(m){
			m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
			colnames(m2) <- colnames(m)
			rownames(m2) <- rownames(m)
			return(m2)
		}

##############################
### LOOP THROUGH EACH CHR FILE:
for (j in 1:22){
	print(paste0("Chromosome: ", j))
	inputData <- read.table(paste0("NURTURECKD_chr", j, ".GTs.txt"),header=T,skip=6)
	
	print(paste0("NCKD - Data loaded for MN status."))
	
	colnames(inputData) <- c(colnames(inputData[,c(1:9)]), gsub("^.{0,61}", "", colnames(inputData[,c(10:ncol(inputData))])))
	
	# ctrl input data:
	ctrl.nckd.chrom <- ctrl.nckd[ctrl.nckd$CHROM %in% j,]
	ctrl.nckd.chrom$CHROM <- as.numeric(ctrl.nckd.chrom$CHROM)
	ctrl.nckd.chrom$POS <- as.numeric(ctrl.nckd.chrom$POS)
	matching_SNPs <- 
		inputData[inputData$POS %in% ctrl.nckd.chrom$POS & inputData$REF %in% ctrl.nckd.chrom$REF 
			& inputData$ALT %in% ctrl.nckd.chrom$ALT,c(1,2,3,4,5)]
	print(paste0("No. matching_SNPs between NCKD and 19 new Ctrls for chr",j,": ", 
		nrow(matching_SNPs), " (was: ", nrow(inputData), ")"))
	matching_SNPs.forInput <- ctrl.nckd.chrom[ctrl.nckd.chrom$ID %in% matching_SNPs$ID,]
	inputData.forInput <- inputData[inputData$ID %in% matching_SNPs$ID,]
	inputData <- merge(inputData.forInput, matching_SNPs.forInput, 
		by = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
		
	# SKS input data:
	inputData.SKS <- read.table(paste0("SKS_chr", j, ".GTs.txt"),header=T,skip=6)
	print(paste0("SKS data loaded."))
	# cut first 60 characters out of rownames
	colnames(inputData.SKS) <- 
		c(colnames(inputData.SKS[,c(1:9)]), gsub("^.{0,61}", "", colnames(inputData.SKS[,c(10:ncol(inputData.SKS))])))
	
	
	inputData <- inputData[(inputData$ID %in% inputData.SKS$ID & inputData$ALT %in% inputData.SKS$ALT
		& inputData$POS %in% inputData.SKS$POS & inputData$CHROM %in% inputData.SKS$CHROM),]
	print(paste0("No. SNPs matching between NCKD, ctrls and SKS: ", nrow(inputData)))
	
	inputData <- merge(inputData, inputData.SKS, by = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
	print(paste0("No. cols: ", ncol(inputData)))
	# 4640
	print(paste0("No. SNPs (rows): ", nrow(inputData)))
	# 481897
	
	inputData <- inputData[!duplicated(inputData$ID),]
	inputData_torun <- inputData[,-c(1:2,4:9)]
	rownames(inputData_torun) <- inputData_torun$ID
	inputData_torun <- t(inputData_torun[,-c(1)])
	
	
	phenD <- phenD.all
	inputData_torun <- inputData_torun[rownames(inputData_torun) %in% rownames(phenD),]
	#phenD_torun <- phenD[phenD$SentrixPosition %in% rownames(inputData_torun),c(4:5)]

	# make sure sentrixposiitons in same order:
	orders <- as.data.frame(rownames(inputData_torun))
	orders$order <- seq.int(nrow(orders))
	names(orders) <- c("sp","order")
	phenD$sp <- rownames(phenD)
	phenD <- merge(phenD, orders, by = c("sp"))
	rownames(phenD) <- phenD$sp
	phenD <- phenD[order(phenD$order),]
	phenD_torun <- phenD[rownames(phenD) %in% rownames(inputData_torun),c(2),drop=F]
	# 1662 cases, now 1702 with controls
	
	write.csv(phenD_torun, "phenD_NCKD_SKS_Ctrls_sentrixpositions.csv", quote=F, row.names=T)
	
	# for all SNPS (not just in genes)
	snps = inputData$ID
	l.snps = length(snps)
	#geneD = data.matrix(as.numeric(inputData_torun),rownames.force=T)
	geneD = as.matrix(inputData_torun)
	
	if(length(which(apply(geneD, 2, var) == 0))!=0){
		geneD <- geneD[ ,-as.numeric(which(apply(geneD, 2, var) == 0))]
	} else{
		geneD_num <- mat_fun(geneD)
		
	}
	nind=nrow(geneD_num)
	print(paste0("nind: ", nind))
	# 1702
	nsnps=1
	ntraits=1
	
	#inputData_torun2 <- inputData_torun[,c(1:5)]
	print('Running CCA:')
	all.can <- as.matrix(apply(geneD_num, MARGIN = 2, cca))
	p <- apply(all.can, MARGIN = 1, fapp_func)
	
	results <- merge(all.can, as.matrix(p), by=0, all=TRUE)
	names(results) <- c('RSID','CCA_r','CCA_pval')
	write.table(results, 
		paste0("NCKD_SKS_combined_plusCtrls_1702_EUR_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_chr",j,".txt"),col.names=T,row.names=F,quote=F)

}
	
	
# load in all results and visualise:
# ================= on bash:
awk '{print $1, $2, $3, $4, $5}' NCKD_allChrs.egfr_adj_scaled.glm.linear \
	> NCKD_allChrs.egfr_adj_scaled.glm.linear.chrposrefalt
# =====================

# HERE - need to run the following:
for (j in 1:22){
	results <- read.table(paste0("NCKD_SKS_combined_plusCtrls_1702_EUR_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_chr",j,".txt"),header=T)
	results <- results[,c(1:3)]
	names(results) <- c("SNP","CCA_r","CCA_pval")
	if(j==1){
		all.results <- results
	} 
	else{
		all.results <- rbind(all.results, results)
	}
}

chrposrefalt <- read.table('NCKD_allChrs.egfr_adj_scaled.glm.linear.chrposrefalt', header=T)
names(chrposrefalt) <- c("chr","pos","SNP","ref","alt")

# ======= Manhattan plot:
library(ggplot2)
library(grid)
library(sqldf)
library(dplyr)
require(scales)
library(stringr)
`%notin%` <- Negate(`%in%`)

data_all <- merge(all.results, chrposrefalt, by = c("SNP"), all.x=T)
names(data_all) <- c("SNP","r","pval","chr","pos","ref","alt")
data_all$pBonf <- data_all$pval*(0.05/5e-8)
nrow(data_all[data_all$pBonf<0.05,])

write.table(data_all, "NCKD_SKS_combined_plusCtrls_1702_EUR_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_chr_dataAll.txt", col.names=T, quote=F, row.names=F)
data_all.pBonf <- data_all[data_all$pBonf<0.05,]
table(data_all.pBonf[,c(4)])
#   2    6   16
#  45 3414    1

# save for FUMA: rsid and pval:
write.table(data_all.pBonf[,c(1,3)], "NCKD_SKS_ctrl_combined_1702_EUR_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_pBonf_forFUMA.txt",
 quote=F, row.names=F)
write.table(data_all.pBonf, "NCKD_SKS_ctrl_combined_1702_EUR_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_pBonf.txt", quote=F, row.names=F)

### Compare SNPs to published datasets - 
# 1/ Xie et al
# 2. Sekula 2017
# 3. Stanescu 2011

stanescu_2011 <- read.csv("stanescu_2011_supp_snps.csv")
stanescu_2011$study <- "stanescu_2011"
stanescu_2011 <- stanescu_2011[,c(1,2,14)]
names(stanescu_2011) <- c("dataset_stanescu","SNP", "study")

sekula_2017 <- read.csv("sekula2017.csv")
sekula_2017$study <- "sekula_2017"
sekula_2017 <- sekula_2017[,c(1,11)]
names(sekula_2017) <- c("SNP", "study")

xie_2020 <- read.table("European_GWAS_meta-analysis.pval.5e_8.txt", header=F)
xie_2020$study <- "xie_2020"
xie_2020 <- xie_2020[,c(1,9)]
names(xie_2020) <- c("SNP", "study")

Xie2020_rep <- read.csv("Xie_2020_table_supp1_snps_selected_replication.csv")
Xie2020_rep$study <- "xie_2020_rep_5e_8"
Xie2020_rep <- Xie2020_rep[Xie2020_rep$P_combined < 5e-8,c(1,16)]
names(Xie2020_rep) <- c("SNP", "study_Xie_rep")

data_all.pBonf.stanescu <- merge(data_all.pBonf, stanescu_2011, by = c("SNP"), all.x=T)
data_all.pBonf.stanescu.sekula <- merge(data_all.pBonf.stanescu, sekula_2017, by = c("SNP"), all.x=T)
data_all.pBonf.stanescu.sekula.xie <- merge(data_all.pBonf.stanescu.sekula, xie_2020, by = c("SNP"), all.x=T)
data_all.pBonf.stanescu.sekula.xie_xie_rep <- merge(data_all.pBonf.stanescu.sekula.xie, Xie2020_rep, by = c("SNP"), all.x=T)
nrow(data_all.pBonf.stanescu.sekula.xie[is.na(data_all.pBonf.stanescu.sekula.xie$study.x) & is.na(data_all.pBonf.stanescu.sekula.xie$study.y) &
	is.na(data_all.pBonf.stanescu.sekula.xie$study),])


# add in 5 lead SNPs:
lead5 <- as.data.frame(c("rs66484345"))
lead5$lead5 <- "yes"
names(lead5) <- c("SNP","lead5")
data_all.pBonf.lead.stanescu.sekula.xie_xie_rep <- merge(data_all.pBonf.stanescu.sekula.xie_xie_rep, lead5, by = c("SNP"), all.x=T)
nrow(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep[is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep$study.x) & is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep$study.y) &
	is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep$study) & data_all.pBonf.lead.stanescu.sekula.xie_xie_rep$lead5 %in% "yes",])
# pla2r lead SNP discovered in Sekula 2017 and Xie 2020 but not replicated.

data_all.pBonf.lead.stanescu.sekula.xie_xie_rep[data_all.pBonf.lead.stanescu.sekula.xie_xie_rep$lead5 %in% "yes",]

# Check GWAS Catalog:
library(gwasrapidd)
results <- data.frame(SNP = character(), assoc_ID = numeric(), traits=character(), url=character(), stringsAsFactors=FALSE)
NCKD_SNPs <- data_all.pBonf$SNP

for (i in 1:length(NCKD_SNPs)){
	print(i)
	SNP = NCKD_SNPs[i]
	results[i,1] = SNP
	gwas.assoc.NCKD_SNPs = get_associations(variant_id = SNP,warnings = FALSE)
	if(nrow(gwas.assoc.NCKD_SNPs@associations)>0){
		print(paste0(SNP, " has prev. pub. association, ID:", gwas.assoc.NCKD_SNPs@associations$association_id))
		assoc_IDs <- gwas.assoc.NCKD_SNPs@associations$association_id
		results[i,2] = paste(assoc_IDs, collapse="_")
		gwas.traits.NCKD_SNPs = get_traits(association_id = assoc_IDs,warnings = FALSE)
		if(nrow(gwas.traits.NCKD_SNPs@traits)>0){
			print(paste0(" and traits: ", paste0(gwas.traits.NCKD_SNPs@traits$trait)))
			results[i,3] = paste(gwas.traits.NCKD_SNPs@traits$trait, collapse="_")
			results[i,4] = paste(gwas.traits.NCKD_SNPs@traits$uri, collapse="_")
		}
	}
	else{
		results[i,2] = "no assoc IDs found"
		results[i,3] = "NA"
		results[i,4] = "NA"
	}
}

# grepl for membranous, neph, kidney,
results$membranous <- NA
results[,5] <- grepl('membranous', results$traits)
results$neph <- NA
results[,6] <- grepl('neph', results$traits)
results$kidney <- NA
results[,7] <- grepl('kidney', results$traits)
results$glom <- NA
results[,8] <- grepl('glom', results$traits)
results$HLA <- NA
results[,9] <- grepl('HLA', results$traits)


write.csv(results, "NCKD_SKS_ctrls_combined_3640pBonfSNPs_MN_vs_oth_vs_published_ALLGWASCat.csv")

results.anyRelevantGWAS <- results[results$membranous %in% TRUE | results$neph %in% TRUE | results$kidney %in% TRUE | 
	results$glom %in% TRUE | results$HLA %in% TRUE,]

# Add to prev. table
data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat <- unique(
	merge(data_all.pBonf.stanescu.sekula.xie_xie_rep, results, by = c("SNP"), all.x=T))
data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat[data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$membranous %in% TRUE | data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$neph %in% TRUE | data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$kidney %in% TRUE | 
	data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$glom %in% TRUE | data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$HLA %in% TRUE,]

# concat Stanescu study results into one row:
data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat.grouped <- unique(data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat %>%
	group_by(SNP) %>%
	mutate(Stanescu_study_all = paste0(dataset_stanescu, collapse="_")) %>%
	select(SNP, r, pval, chr, pos, ref, alt, pBonf, Stanescu_study_all, study.x, study.y, study, study_Xie_rep, assoc_ID, traits, url, membranous, neph, kidney, glom, HLA))

write.csv(data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat.grouped, 
	"NCKD_SKS_ctrls_combined_3640pBonfSNPs_alldata_MN_vs_oth_vs_published_ALLGWASCat_v2.csv")

data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat[data_all.pBonf.stanescu.sekula.xie_xie_rep.GWASCat$chr==2,]

### followig not used:
length(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat[is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.x) & 
	is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.y) & is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study)
		& is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study_Xie_rep) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$membranous) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$neph) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$kidney) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$glom) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$HLA),1])
			
data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported <- 
ifelse((is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.x) & 
	is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.y) & is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study)
		& is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$study_Xie_rep) & 
		is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$membranous) & is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$neph) & 
			is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$kidney) & is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$glom) & 
				is.na(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$HLA)), data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$SNP, "")

length(unique(data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat[data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported != "",1]))
# 345

# how many of these are assoc with other non-MN related traits in GWAS catalog?
snps345_new <- data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat[data_all.pBonf.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported != "",]
snps345_new.GWASCat_unrel_traits <- merge(snps345_new, results, by = c("SNP"), all.x=T)
write.csv(snps345_new.GWASCat_unrel_traits, "MN_vs_oth_943SNPs_vs_published_GWASCat_345unreported_other_traits.csv")
#####################

# save VEP input format: chr pos . ref alt . . . 
snps_forVEP <- data_all.pBonf[,c(4, 5, 1, 6, 7)]
snps_forVEP$extra1 <- "."
snps_forVEP$extra2 <- "."
snps_forVEP$extra3 <- "."
snps_forVEP <- snps_forVEP[order(snps_forVEP$chr, snps_forVEP$pos),]
write.table(snps_forVEP, "NCKD_SKS_controls_combined_3640pBonfSNPs_forVEP.txt", quote=F, row.names=F, col.names=F)


# sort out VEP output table (multiple rows per SNps due to different transcripts..)
fromVEP <- read.table("VEP_out_3640SNPs_combined_NCKD_SKS_ctrls_MN.txt", header=T)
# cols thatn eed pasting nto one: Consequence, SYMBOL, Gene, Feature, BIOTYPE, DISTANCE, STRAND, FLAGS, HGNC_ID, TSL, SQISSPROT, TREMBL, UNIPAC, UNIPROT_ISOFORM
# get rid of: TSL	SWISSPROT	TREMBL	UNIPARC	UNIPROT_ISOFORM

fromVEP$Existing_variation <- gsub(",", ".", fromVEP$Existing_variation)

fromVEP_min <- unique(fromVEP[,-c(11, 12, 24, 28:31, 34:39)] %>% 
	group_by(Uploaded_variation, Location, Allele, REF_ALLELE, UPLOADED_ALLELE) %>% 
	mutate(Consequence = paste0(unique(gsub(",", ".", Consequence)), collapse = "_")) %>%
	mutate(IMPACT = paste0(unique(IMPACT), collapse = "_")) %>%
	mutate(SYMBOL = paste0(unique(SYMBOL), collapse = "_")) %>%
	mutate(Gene = paste0(unique(Gene), collapse = "_")) %>%
	mutate(Feature_type = paste0(unique(Feature_type), collapse = "_")) %>%
	mutate(Feature = paste0(unique(Feature), collapse = "_")) %>%
	mutate(BIOTYPE = paste0(unique(BIOTYPE), collapse = "_")) %>%
	mutate(DISTANCE = paste0(unique(DISTANCE), collapse = "_")) %>%
	#mutate(STRAND = paste0(unique(STRAND), collapse = "_")) %>%
	mutate(FLAGS = paste0(unique(gsub(",", ".", FLAGS)), collapse = "_")) %>%
	mutate(HGNC_ID = paste0(unique(HGNC_ID), collapse = "_")) %>%
	#mutate(INTRON = paste0(unique(as.character(INTRON)), collapse = "_")) %>%
	mutate(HGVSc = paste0(unique(HGVSc), collapse = "_")) %>%
	mutate(HGVSp = paste0(unique(HGVSp), collapse = "_")) %>%
	mutate(cDNA_position = paste0(unique(cDNA_position), collapse = "_")) %>%
	mutate(CDS_position = paste0(unique(CDS_position), collapse = "_")) %>%
	mutate(Protein_position = paste0(unique(Protein_position), collapse = "_")) %>%
	mutate(Amino_acids = paste0(unique(Amino_acids), collapse = "_")) %>%
	mutate(Codons = paste0(unique(Codons), collapse = "_")) %>%
	mutate(SYMBOL_SOURCE = paste0(unique(SYMBOL_SOURCE), collapse = "_")) %>%
	mutate(TRANSCRIPTION_FACTORS = paste0(unique(gsub(",", ".", TRANSCRIPTION_FACTORS)), collapse = "_")) %>%
	mutate(MOTIF_NAME = paste0(unique(gsub(",", ".", MOTIF_NAME)), collapse = "_")) %>%
	mutate(MOTIF_POS = paste0(unique(gsub(",", ".", MOTIF_POS)), collapse = "_")) %>%
	mutate(HIGH_INF_POS = paste0(unique(gsub(",", ".", HIGH_INF_POS)), collapse = "_")) %>%
	mutate(MOTIF_SCORE_CHANGE = paste0(unique(gsub(",", ".", MOTIF_SCORE_CHANGE)), collapse = "_")) %>%
	mutate(am_class = paste0(unique(gsub(",", ".", am_class)), collapse = "_")) %>%
	mutate(am_pathogenicity = paste0(unique(gsub(",", ".", am_pathogenicity)), collapse = "_")) %>%
	mutate(SIFT = paste0(unique(gsub(",", ".", SIFT)), collapse = "_")) %>%
	mutate(PolyPhen = paste0(unique(gsub(",", ".", PolyPhen)), collapse = "_")))


write.csv(fromVEP_min, "VEP_out_3640SNPs_combined_NCKD_SKS_ctrl_MN_formatted.csv", quote=F, row.names=F)
	

############################################
data_sig <- data_all[data_all$pval<0.05,]
data_all$chr <- as.numeric(data_all$chr)
data_sig$chr <- as.numeric(data_sig$chr)
data_nonsig <- data_all[data_all$SNP %notin% data_sig$SNP,]
data_sig <- data_sig[ order(data_sig$chr),]
data_nonsig <- data_nonsig[ order(data_nonsig$chr),]

gwas.dat2 <- rbind(data_sig,data_nonsig)
gwas.dat2 <- data_all[!is.na(data_all$chr),]
gwas.dat2$chr <- as.numeric(gwas.dat2$chr)
nCHR <- length(unique(gwas.dat2$chr))
gwas.dat2$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwas.dat2$chr)){
  nbp[i] <- max(gwas.dat2[gwas.dat2$chr == i,]$pos)
  gwas.dat2[gwas.dat2$chr == i,"BPcum"] <- gwas.dat2[gwas.dat2$chr == i,"pos"] + s
  s <- s + nbp[i]
}

# ==== slice:
#notsig.dat <- gwas.dat2 %>%
#  slice(sample(nrow(.), nrow(.) / 10))
#gwas.dat <- rbind(data_sig[,c(1:3,5:10)],data_nonsig[,c(1:3,5:10)])
#gwas.dat <- rbind(data_sig,data_nonsig)
#gwas.dat <- notsig.dat

# ==== slice - those under 0.05:
notsig.dat <- data_nonsig %>%
  slice(sample(nrow(.), nrow(.) / 10))
gwas.dat <- rbind(data_sig,notsig.dat)

# ==== no slice:
#gwas.dat <- gwas.dat2

gwas.dat$chr <- as.numeric(gwas.dat$chr)
gwas.dat <- gwas.dat[order(gwas.dat$chr),]
nCHR <- length(unique(gwas.dat$chr))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwas.dat$chr)){
  print(i);
  nbp[i] <- max(gwas.dat[gwas.dat$chr == i,]$pos)
  print(nbp[i]);
  gwas.dat[gwas.dat$chr == i,"BPcum"] <- gwas.dat[gwas.dat$chr == i,"pos"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(gwas.dat$pBonf)))) + 2

gwas.dat$pval_minuslog10 <- -log10(gwas.dat$pval)
gwas.dat$pBonf_minuslog10 <- -log10(gwas.dat$pBonf)
#gwas.dat$pBonf_manual <- gwas.dat$pval*nrow(gwas.dat)
sig <- 0.05

#box_xmax <- (min(gwas.dat[gwas.dat$chr==4,9]))+((max(gwas.dat[gwas.dat$chr==4,9])-(min(gwas.dat[gwas.dat$chr==4,9])))/3)

gwas.dat.nurture <- gwas.dat
gwas.dat.nurture$dataset <- "nurture-ckd"

# Add gene-based results to plot (as square shape):
#gene_based_snps <- read.csv("genebased_forplot.csv",header=T)
#gene_based_snps <- separate_rows(gene_based_snps,snp, sep=", ", convert = TRUE)
#names(gene_based_snps) <- c("SNP","gene","dataset","pval_adj")
#gwas.dat.gene_based <- merge(gene_based_snps[gene_based_snps$dataset=="nurture-ckd",],gwas.dat,by=c("SNP"),all.x=T)
#gene_based_snps_info <- gwas.dat.gene_based[,c(1:4,9,13)]
#gene_based_snps_info$pval_minuslog10 <- -log10(gene_based_snps_info$pval_adj)
#sql <- sqldf("SELECT MAX(BPcum), SNP, dataset, gene, pval_adj, chr, BPcum, pval_minuslog10 FROM gene_based_snps_info WHERE dataset == 'nurture-ckd' GROUP BY gene, dataset")
# check max y axis:
max(gwas.dat$pval_minuslog10)
# 22.0785

#manhplot_nurture <- ggplot(gwas.dat, aes(x = BPcum, y = pBonf_minuslog10, color = as.factor(chr), group = as.factor(chr))) +
manhplot_nurture <- ggplot(gwas.dat, aes(x = BPcum, y = pval_minuslog10, color = as.factor(chr), group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1,show.legend = F) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "dotted") +
  geom_hline(yintercept = -log10(5e-8), color = "black", linetype = "dashed") +
  #geom_hline(yintercept = sig, color = "black", linetype = "dashed") +
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==2,9]), xmax = max(gwas.dat[gwas.dat$chr==2,9]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==4,9]), xmax = max(gwas.dat[gwas.dat$chr==4,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==6,9]), xmax = max(gwas.dat[gwas.dat$chr==6,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==8,9]), xmax = max(gwas.dat[gwas.dat$chr==8,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==10,9]), xmax = max(gwas.dat[gwas.dat$chr==10,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==12,9]), xmax = max(gwas.dat[gwas.dat$chr==12,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==14,9]), xmax = max(gwas.dat[gwas.dat$chr==14,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==16,9]), xmax = max(gwas.dat[gwas.dat$chr==16,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==18,9]), xmax = max(gwas.dat[gwas.dat$chr==18,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==20,9]), xmax = max(gwas.dat[gwas.dat$chr==20,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==22,9]), xmax = max(gwas.dat[gwas.dat$chr==22,9]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  #annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==1,10]), xmax = box_xmax, ymin = 4.9,ymax =6, fill = "white", color="black")+
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
 # scale_y_continuous(limits = c(0, 1.4), breaks = c(seq(0, 1.4, by = 0.2)),expand=c(0.04, 0)) +
 # scale_y_continuous(limits = c(0, 1), breaks = c(seq(0, 1, by = 0.2)),expand=c(0.04, 0)) +
  scale_y_continuous(breaks = c(seq(0, 30, by = 2)),expand=c(0.04, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  #geom_point(data=sql, size=2, color="black", shape=5)+
  labs(x = "Chromosomal location",
	   y = expression(bold(-log[10](CCA~bolditalic(P)~value)))) +
	   #y = expression(bold(CCA~bolditalic(P)~value[adjusted]))) +
  #ggtitle("NURTuRE-CKD")+
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",margin = margin(t = 0, r = 12, b = 0, l = 0)),
    plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    plot.tag = element_text(face="bold", size = 20,color="black")
)

jpeg("manh_plot_NCKD_SKS_EUR_ctrls_combined_1702_CCA_binary_CKD_PT_only_MN_vs_all_oth_imp_R2_0.9.CCAout.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_nurture, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()

