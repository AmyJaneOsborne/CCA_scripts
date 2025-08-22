source("gene-based.cca3.2.R")

library(sqldf)
library(dplyr)
library(GWASTools)

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
phenD <- read.csv("CKD_PTs_top11.binary.csv", header=T)
phenD <- phenD[!is.na(phenD$SentrixPosition),]
row.names(phenD) <- phenD$SentrixPosition
# set col names: patient_id case_control egfr_adj_scaled urea_adj_scaled SentrixPosition
# To make analyses the same as SKS, exclude ""Microscopic polyangiitis - histologically proven" and ""Primary reflux nephropathy - sporadic" from NCKD analysis

# Try MN vs all other only:
phenD <- phenD[,c(1,2,10)]
write.csv(phenD, "phenD_MNvsAllOther.csv", quote=F, row.names=F, col.names=T)

##############################
### LOOP THROUGH EACH CHR FILE:

for (j in 1:22){
	print(paste0("Chromosome: ", j))
	inputData <- read.table(paste0(
	"data_chr", j, ".txt"),header=T,skip=6)
	print(paste0("NCKD - Data loaded for eGFR slope adj. and TAP imp. CCA with binary PT."))
	
	inputData <- inputData[!duplicated(inputData$ID),]
	inputData_torun <- inputData[,-c(1:2,4:9)]
	rownames(inputData_torun) <- inputData_torun$ID
	inputData_torun <- t(inputData_torun[,-c(1)])
	
	rownames(inputData_torun) <- gsub("^.{0,61}", "", rownames(inputData_torun))
	inputData_torun <- inputData_torun[rownames(inputData_torun) %in% phenD$SentrixPosition,]
	
	`%notin%` <- Negate(`%in%`)
	inputData_torun <- inputData_torun[rownames(inputData_torun) %notin% "203422530068_R05C01",]

	phenD_torun <- phenD[phenD$SentrixPosition %in% rownames(inputData_torun),c(3)]
	
	# for all SNPS (not just in genes)
	snps = inputData$ID
	l.snps = length(snps)
	geneD = as.matrix(inputData_torun)
		
	if(length(which(apply(geneD, 2, var) == 0))!=0){
		geneD <- geneD[ -as.numeric(which(apply(geneD, 2, var) == 0))]
	}
	nind=nrow(geneD)
	nsnps=1
	ntraits=1
	
	print('Running CCA:')
	all.can <- as.matrix(apply(geneD, MARGIN = 2, cca))
	p <- apply(all.can, MARGIN = 1, fapp_func)
	
	results <- merge(all.can, as.matrix(p), by=0, all=TRUE)
	names(results) <- c('RSID','CCA_r','CCA_pval')
	write.table(results, 
		paste0("NCKD_EUR_2181_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_chr",j,".txt"),col.names=T,row.names=F,quote=F)

}
	
	
# load in all results and visualise:
# ================= on bash:
awk '{print $1, $2, $3, $4, $5}' NCKD_allChrs.egfr_adj_scaled.glm.linear \
	> NCKD_allChrs.egfr_adj_scaled.glm.linear.chrposrefalt
# =====================

# HERE - need to run the following:
for (j in 1:22){
	results <- read.table(paste0("NCKD_EUR_2181_CCA_binary_CKD_PT_only_imp_R2_0.9.MN_vs_allOth_only_CCAout_chr",j,".txt"),header=T)
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
#nrow(data_all[data_all$pBonf<0.05,])

write.table(data_all, "NCKD_EUR_2181_CCA_binary_CKD_PT_only_MN_vs_all_other_imp_R2_0.9.CCAout_chr_dataAll.txt", col.names=T, quote=F, row.names=F)

data_all.pBonf <- data_all[data_all$pBonf<0.05,]
# 2003

# for LDlink:
write.table(data_all.pBonf[,c(4,5,3)], "MN_2003SNPs_pval_forLDlink.txt", quote=F, row.names=F, col.names=F)


data_all.pBonf.sks.pBonf <- data_all.pBonf[data_all.pBonf$SNP %in% sks.pBonf$SNP,]
# 943
data_all.pBonf.sks.pBonf <- merge(data_all.pBonf, sks.pBonf, by = c("SNP","chr","pos","ref","alt"))

names(data_all.pBonf.sks.pBonf) <- c("SNP","chr","pos","ref","alt","CCA_r.NCKD","CCA_pval.NCKD","pBonf.NCKD",
	"CCA_r.SKS","CCA_pval.SKS")
write.table(data_all.pBonf.sks.pBonf, 
	"NCKD_EUR_2181_CCA_binary_CKD_PT_only_MN_vs_all_oth_imp_R2_0.9.CCAout_chr_dataAll.OverlapSKSsame.pBonf.txt", 
		col.names=T, row.names=F, quote=F)

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
lead5 <- as.data.frame(c("rs3132450", "rs79350521", "rs1064994", "rs28630740", "rs2395231"))
lead5$lead5 <- "yes"
names(lead5) <- c("SNP","lead5")
snps943.lead.stanescu.sekula.xie_xie_rep <- merge(snps943.stanescu.sekula.xie_xie_rep, lead5, by = c("SNP"), all.x=T)
nrow(snps943.lead.stanescu.sekula.xie_xie_rep[is.na(snps943.lead.stanescu.sekula.xie_xie_rep$study.x) & is.na(snps943.lead.stanescu.sekula.xie$study.y) &
	is.na(snps943.lead.stanescu.sekula.xie_xie_rep$study) & snps943.lead.stanescu.sekula.xie_xie_rep$lead5 %in% "yes",])
# 3: rs2395231,  rs28630740, rs79350521

snps943.lead.stanescu.sekula.xie_xie_rep[snps943.lead.stanescu.sekula.xie_xie_rep$lead5 %in% "yes",]
# rs1064994: Sekula 2017 discovery  - dont know if replicated since only looked at replication of index variants.
# rs3132450: Xie 2020 - but didnt replicate (combined analysis Pvalue)

# Check GWAS Catalog:
library(gwasrapidd)
results <- data.frame(SNP = character(), assoc_ID = numeric(), traits=character(), url=character(), stringsAsFactors=FALSE)
NCKD_SNPs <- snps943.stanescu.sekula.xie_xie_rep$SNP

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

write.csv(results, "MN_vs_oth_943SNPs_vs_published_ALLGWASCat.csv")

results.anyRelevantGWAS <- results[results$membranous %in% TRUE | results$neph %in% TRUE | results$kidney %in% TRUE | 
	results$glom %in% TRUE | results$HLA %in% TRUE,]

# Add to prev. table
snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat <- unique(
	merge(snps943.lead.stanescu.sekula.xie_xie_rep, results.anyRelevantGWAS, by = c("SNP"), all.x=T))
snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat[snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$membranous %in% TRUE | snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$neph %in% TRUE | snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$kidney %in% TRUE | 
	snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$glom %in% TRUE | snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$HLA %in% TRUE,]

write.csv(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat, "MN_vs_oth_943SNPs_vs_published_GWASCat.csv")

length(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat[is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.x) & 
	is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.y) & is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study)
		& is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study_Xie_rep) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$membranous) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$neph) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$kidney) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$glom) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$HLA),1])
			
snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported <- 
ifelse((is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.x) & 
	is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study.y) & is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study)
		& is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$study_Xie_rep) & 
		is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$membranous) & is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$neph) & 
			is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$kidney) & is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$glom) & 
				is.na(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$HLA)), snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$SNP, "")

length(unique(snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat[snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported != "",1]))
# 345

# how many of these are assoc with other non-MN relaated traits in GWAS catalog?
snps345_new <- snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat[snps943.lead.stanescu.sekula.xie_xie_rep.GWASCat$not_prev_reported != "",]
snps345_new.GWASCat_unrel_traits <- merge(snps345_new, results, by = c("SNP"), all.x=T)
write.csv(snps345_new.GWASCat_unrel_traits, "MN_vs_oth_943SNPs_vs_published_GWASCat_345unreported_other_traits.csv")


# tidy up VEP output table (multiple rows per SNP due to different transcripts..)
fromVEP <- read.table("VEP_out_943SNPs.txt", header=T)
# cols that need pasting into one: Consequence, SYMBOL, Gene, Feature, BIOTYPE, DISTANCE, STRAND, FLAGS, HGNC_ID, TSL, SQISSPROT, TREMBL, UNIPAC, UNIPROT_ISOFORM
# get rid of: TSL	SWISSPROT	TREMBL	UNIPARC	UNIPROT_ISOFORM

fromVEP$Existing_variation <- gsub(",", ".", fromVEP$Existing_variation)

fromVEP_min <- unique(fromVEP[,-c(11, 12, 24, 28:32, 37:39)] %>% 
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


write.csv(fromVEP_min, "VEP_out_943SNPs_formatted.csv", quote=F, row.names=F)
	
		
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

# check max y axis:
max(gwas.dat$pval_minuslog10)
# 22.0785

manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = pval_minuslog10, color = as.factor(chr), group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6,show.legend = F) +
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
  scale_y_continuous(breaks = c(seq(0, 18, by = 2)),expand=c(0.04, 0)) +
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

jpeg("manh_plot_NCKD_EUR_2181_CCA_binary_CKD_PT_only_MN_vs_all_oth_imp_R2_0.9.CCAout.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()


#### Assess overlap with SKS:
nckd_results.sig.inSKS.sig <- 
	read.table("NCKD_EUR_2181_CCA_binary_CKD_PT_only_MN_vs_all_oth_imp_R2_0.9.CCAout_chr_dataAll.OverlapSKSsame.pBonf.txt",header=T)
# 943

# create input file for annovar and run annovar:
for_Annovar <- nckd_results.sig.inSKS.sig[, c(2,3,3,4,5,1)]
write.table(for_Annovar,"NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs.forAnnovar.txt", col.names=F,row.names=F,quote=F)

### in bash:
perl ~/annovar/annovar/annotate_variation.pl -out \
~/GWAS_N_APRIL2021/for_CCA/NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs_annovar_annot_hg38 -build hg38 \
~/GWAS_N_APRIL2021/for_CCA/NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs.forAnnovar.txt ~/annovar/annovar/humandb/

# GWAS catalog:
library(gwasrapidd)

SNPs <- nckd_results.sig.inSKS.sig$SNP
results <- data.frame(SNP = character(), assoc_ID = numeric(), traits=character(), url=character(), stringsAsFactors=FALSE)

for (i in 1:length(NCKD_SNPs)){
	print(i)
	SNP = SNPs[i]
	results[i,1] = SNP
	gwas.assoc.NCKD_SNPs = get_associations(variant_id = SNP,warnings = FALSE)
	if(nrow(gwas.assoc.NCKD_SNPs@associations)>0){
		print(paste0(SNP, "_has_prev._pub._association,_ID:", gwas.assoc.NCKD_SNPs@associations$association_id))
		assoc_IDs <- gwas.assoc.NCKD_SNPs@associations$association_id
		results[i,2] = paste(assoc_IDs, collapse="_")
		gwas.traits.NCKD_SNPs = get_traits(association_id = assoc_IDs,warnings = FALSE)
		if(nrow(gwas.traits.NCKD_SNPs@traits)>0){
			print(paste0(" and traits: ", paste0(gwas.traits.NCKD_SNPs@traits$trait)))
			results[i,3] = gsub(" ", "_", paste(gwas.traits.NCKD_SNPs@traits$trait, collapse="__"))
			results[i,4] = paste(gwas.traits.NCKD_SNPs@traits$uri, collapse="_")
		}
	}
	else{
		results[i,2] = "no_assoc_IDs_found"
		results[i,3] = "NA"
		results[i,4] = "NA"
	}
}

write.table(results, "NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs.gwascatalog.txt", col.names=T, quote=F, row.names=F)


# now Qtlizer:
library("Qtlizer")
inputData.RSIDs.Qlitzer <- get_qtls(nckd_results.sig.inSKS.sig$SNP, ref_version="hg38")
# 943 found
# 190664 data point(s) received

# ===== FDR5%, kidney, immune cells, WBCs, liver only:
cells <- c("Cells - Lymphoblastoid cell lines", "Cells - Monocytes LPS2","Cells - Monocytes LPS24","Cells - Monocytes IFN","Cells - CD4+ lymphocytes","Cells - Macrophages","Cells - Peripheral blood mononuclear cells",
	"Kidney - Tubulointerstitial (Nephrotic syndrome cases)","Kidney - Glomerulus (Nephrotic syndrome cases)","Cells - Monocytes","Whole blood","Kidney - Cortex",
	"Peripheral blood","Liver")
inputData.Qlitzer.FDRandFWER5.kidney.imm.blood.liver <- inputData.RSIDs.Qlitzer[(inputData.RSIDs.Qlitzer$sign_info %in%
	c("FDR<5%","FWER<5%","Q-value; FDR<5%") & inputData.RSIDs.Qlitzer$tissue %in% cells & inputData.RSIDs.Qlitzer$is_sign %in% c("true")),]
length(unique(inputData.Qlitzer.FDRandFWER5.kidney.imm.blood.liver$query_term))


for (i in 1:length(cells)){
	celltype = cells[i]
	data_cell <- inputData.Qlitzer.FDRandFWER5.kidney.imm.blood.liver[inputData.Qlitzer.FDRandFWER5.kidney.imm.blood.liver$tissue %in% celltype,]
	print(i)
	data_cell.agg <- data_cell %>%
		group_by(query_term) %>%
		# to get genes only, per cell type
		mutate(Qtlizer_genes = paste(gene,sep="_")) %>%
		select(query_term, Qtlizer_genes) %>%
		#mutate(celltype = paste(gene,colocalization,source,pmid,sep="_")) %>%
		#select(query_term, celltype) %>%
		summarise_all(funs(toString(na.omit(.))))
		#mutate(celltype = gsub(", ","__",celltype))
	if(nrow(data_cell.agg)>0){
		names(data_cell.agg) <- c("query_term", paste(gsub(" ", "_", celltype)))
		if(i==1){
			results <- data_cell.agg
		} else {
			results <- merge(results, data_cell.agg, by = c("query_term"), all=T)
		}
	}
}

names(results) <- c("query_term", "Lymphoblastoid_cell_lines", "Monocytes_LPS2", "Monocytes_LPS24", "Monocytes_IFN", "Macrophages", "Kidney_Tubulointerstitial", "Kidney_Glomerulus",
	"Monocytes", "Whole_blood", "Kidney_Cortex", "Peripheral_blood", "Liver")
results$Lymphoblastoid_cell_lines <- gsub(", ","_",results$Lymphoblastoid_cell_lines)
results$Monocytes_LPS2 <- gsub(", ","_",results$Monocytes_LPS2)
results$Monocytes_LPS24 <- gsub(", ","_",results$Monocytes_LPS24)
results$Monocytes_IFN <- gsub(", ","_",results$Monocytes_IFN)
results$Macrophages <- gsub(", ","_",results$Macrophages)
results$Kidney_Tubulointerstitial <- gsub(", ","_",results$Kidney_Tubulointerstitial)
results$Kidney_Glomerulus <- gsub(", ","_",results$Kidney_Glomerulus)
results$Monocytes <- gsub(", ","_",results$Monocytes)
results$Whole_blood <- gsub(", ","_",results$Whole_blood)
results$Kidney_Cortex <- gsub(", ","_",results$Kidney_Cortex)
results$Peripheral_blood <- gsub(", ","_",results$Peripheral_blood)
results$Liver <- gsub(", ","_",results$Liver)

# add in annovar annotation
annovar <- read.table("NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs_annovar_annot_hg38.variant_function", header=F)
names(annovar) <- c("loc", "gene", "chr","pos", "pos2", "ref", "alt", "query_term")
annovar.results <- merge(annovar, results, by = c("query_term"), all.x=T)

annovar.results$gene <- gsub(",","__",annovar.results$gene)

write.csv(annovar.results, "NCKD_SKS_imp_CKDsubtype_only_MS_vs_all_oth_CCA_repl_943SNPs.Qtlizer.kidney_imm_liver_blood.geneInfo_only.csv", 
col.names=T, quote=F, row.names=F)
