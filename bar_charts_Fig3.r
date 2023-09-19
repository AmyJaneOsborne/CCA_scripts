# This script compares metaCCA results with published results
`%notin%` <- Negate(`%in%`)
library(ggplot2)
library(grid)
library(ggbreak)
library(stringr)
library(tidyr)

wuttke_egfr <- read.table("20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt", header=T)
wuttke_bun <- read.table("BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt", header=T)

myresults <- read.csv("supptable3_all_ckdgen_metacca.csv", header=T)
myresults$RSID <- myresults$dbSNP.RefSNP.ID

wuttke_egfr.in.myresults <- wuttke_egfr[wuttke_egfr$RSID %in% myresults$dbSNP.RefSNP.ID,]
# 982

wuttke_egfr.pBonf.in.myresults <- wuttke_egfr.in.myresults[wuttke_egfr.in.myresults$P.value < 0.05/nrow(wuttke_egfr),]
# 789

wuttke_egfr.pBonf.wuttke_bun.p0.05.in.myresults <- wuttke_bun[wuttke_bun$RSID %in% wuttke_egfr.pBonf.in.myresults$RSID & wuttke_bun$P.value < 0.05,]
# 445 snps


# # # # ========= Get all 20833 sig snps:
ckdgen_sig_all <- read.table("CKDGen_EUR_eGFR_BUN_all8346783SNPs_metaCCA_out_wchrpos.txt", header=T)
ckdgen_sig_all <- ckdgen_sig_all[ckdgen_sig_all$pBonf < 0.05,]
# 20833

# How many were sig. in Wuttke et al 2019 paper (Methods): 
#"The genome-wide significance level was set at 5 × 10-8")[for eGFR]
wuttke_egfr_pBonf.ckdgen_sig_all <- wuttke_egfr[wuttke_egfr$RSID %in% ckdgen_sig_all$RSID & wuttke_egfr$P.value < 5e-8,]
# 16545

# === BUN: " Support for kidney function relevance was categorized as ‘likely’ for all eGFR index SNPs with an inverse, significant (one-sided P < 0.05) association with BUN 
# for a given reference allele, 
#‘inconclusive’ for eGFR index SNPs whose effect on BUN was not different from 0 (P => 0.05) 
# and ‘unlikely’ for all eGFR index SNPs with a concordant (i.e. same direction), significant (one-sided P < 0.05) association with BUN for a given reference allele."

wuttke_egfr_pBonf.bun0.05.ckdgen_sig_all <- 
	wuttke_bun[wuttke_bun$RSID %in% wuttke_egfr_pBonf.ckdgen_sig_all$RSID & wuttke_bun$P.value < 0.05,]
# 8833 snps

wuttke_egfr_pBonf.ckdgen_sig_all.wuttke_bun <- merge(wuttke_egfr_pBonf.ckdgen_sig_all, wuttke_bun, by = c("RSID"), all.x=T)
wuttke_egfr_pBonf.ckdgen_sig_all.wuttke_bun.dir.1tailed0.05 <- wuttke_egfr_pBonf.ckdgen_sig_all.wuttke_bun %>%
	mutate(egfr_dir = ifelse(Effect_eGFR<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(Effect_BUN<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="+" & (P.value.y<0.05)) | (egfr_dir== "+" & bun_dir=="-" & (P.value.y<0.05)))
# 7271 SNPs

wuttke_egfr_bun <- merge(wuttke_bun, wuttke_egfr, by = c("RSID"))
wuttke_egfr_bun_wrong_dir <- wuttke_egfr_bun %>% 
	mutate(egfr_dir = ifelse(Effect_eGFR<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(Effect_BUN<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="-") | (egfr_dir== "+" & bun_dir=="+"))

# No. SNPs ID'ed by Wuttke for both eGFR (5e-8) and BUN (0.05) + effect sizes in correct directions:
wuttke_egfr_bun <- merge(wuttke_bun, wuttke_egfr, by = c("RSID"))
wuttke_egfr_bun_sig_dir <- wuttke_egfr_bun %>% 
	mutate(egfr_dir = ifelse(Effect_eGFR<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(Effect_BUN<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-")) %>%
	filter((0.5*P.value.x) < 0.05 & P.value.y < 5e-8)
length(unique(wuttke_egfr_bun_sig_dir$RSID))
# 8875

ckdgen_sig_all_dir <- ckdgen_sig_all %>% 
	mutate(egfr_dir = ifelse(Effect_eGFR<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(Effect_BUN<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-"))
# 14045

wuttke_egfr_bun_sig_dir$wuttke_egfr_bun_sig_dir <- "yes"
wuttke_egfr_bun$wuttke_egfr_bun <- "yes"

# plot comparing these with metaCCA out
wuttke_egfr_bun.wuttke_egfr_bun_sig_dir <- merge(wuttke_egfr_bun[,c(1,20)], wuttke_egfr_bun_sig_dir[,c(1,22)], by = c("RSID"),all.x=T)
ckdgen_sig_all.wuttke_egfr_bun_sig_dir <- merge(wuttke_egfr_bun.wuttke_egfr_bun_sig_dir,ckdgen_sig_all,by = c("RSID"),all.x=T)
wuttke_leadSNPs <- read.table("wuttke_kidrel1_122snps_mvp.txt", header=F)
names(wuttke_leadSNPs) <- c("RSID","chr","pos","loci")
wuttke_leadSNPs$wuttke_leadSNPs <- "yes"

ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead <- merge(ckdgen_sig_all.wuttke_egfr_bun_sig_dir, wuttke_leadSNPs[,c(1,5)],
	by = c("RSID"),all.x=T)

ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_egfr_bun_sig_dir <- 
	ifelse(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_egfr_bun_sig_dir %in% "yes", "yes","no")

ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_leadSNPs <- 
	ifelse(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_leadSNPs %in% "yes", "yes","no")

ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat <- 
	ifelse(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_leadSNPs %in% "yes", "lead",
		ifelse(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_egfr_bun_sig_dir %in% "yes", "not lead but sig.",
			"not sig."))

ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead <- ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead %>% 
	mutate(egfr_dir = ifelse(Effect_eGFR<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(Effect_BUN<0, "-", "+")) %>%
	mutate(metacca_sig_dir = ifelse((((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-")) & !is.na(pval)),"yes","no"))

forplot.ckdgen <- data.frame(wuttke_cat = character(), metacca_sig_dir = character(), count_SNPs = numeric(), stringsAsFactors=FALSE)

forplot.ckdgen[1,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[1]
forplot.ckdgen[2,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[2]
forplot.ckdgen[3,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[3]
forplot.ckdgen[4,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[1]
forplot.ckdgen[5,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[2]
forplot.ckdgen[6,1] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat)[3]

forplot.ckdgen[1,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot.ckdgen[2,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot.ckdgen[3,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot.ckdgen[4,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]
forplot.ckdgen[5,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]
forplot.ckdgen[6,2] <- unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]

forplot.ckdgen[1,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "not sig." &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot.ckdgen[2,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "not lead but sig." &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot.ckdgen[3,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "lead" &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot.ckdgen[4,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "not sig." &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))
		
forplot.ckdgen[5,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "not lead but sig." &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))

forplot.ckdgen[6,3] <- length(unique(ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[
	ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "lead" &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))
	
notleadbutsig <- ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead[ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_cat %in% "not lead but sig." &
		ckdgen_sig_all.wuttke_egfr_bun_sig_dir.lead$wuttke_egfr_bun_sig_dir %in% "yes",]


forplot.ckdgen$wuttke_cat_toshow <- ifelse(forplot.ckdgen$wuttke_cat %in% c("not sig."), "1. Not significant and/or effect sizes not compatible (n = 8,336,431)",
	ifelse(forplot.ckdgen$wuttke_cat %in% c("lead"), 
		"3. Significant and effect sizes compatible, replicated and lead (n = 122)", 
		"2. Significant and effect sizes compatible, but not lead (n = 10,230)"))
	
forplot.ckdgen$metacca_sig_dir_toshow <- ifelse(forplot.ckdgen$metacca_sig_dir %in% c("yes"), "Yes\n(n = 14,045)", "No\n(n = 8,332,738)")
		
plot.ckdgen <- ggplot(forplot.ckdgen, 
			aes(x = stringr::str_wrap(wuttke_cat_toshow,20), y = count_SNPs, fill = metacca_sig_dir_toshow)) +
		geom_bar(position="dodge", stat="identity")+
		geom_text(aes(label = count_SNPs, y = ifelse(count_SNPs<100, count_SNPs+300, count_SNPs-300)), colour = "black", size = 2.5, position = position_dodge(.9))+
		#vjust = 1.5, 
		#scale_fill_grey(start = 0.8, end = 0.6, name="SNP identified\nby metaCCA")+
		scale_fill_discrete(name="Identified by metaCCA?")+
		scale_y_break(c(8500, 8329000)) +
		scale_y_continuous(breaks = seq(0, 8331000, by = 2000), position = "left")+
		xlab("Outcome by Wuttke et al., 2019 genome-wide association study")+ 
        ylab("Number of single\nnucleotide polymorphisms")+
		labs(tag = "A")+
    theme(
    #legend.position = "none",
	#guides(fill=guide_legend(title="SNP identified by metaCCA")),
    legend.spacing = unit(0.7, 'cm'),
	legend.text=element_text(family="Arial",face="bold", size = 8,color="black"),
	legend.title=element_text(family="Arial",face="bold", size = 8,color="black"),
	panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.3),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",angle=90),
	  plot.tag = element_text(family="Arial",face="bold", size = 20,color="black")
)

jpeg("plot_CKDGen_metaCCA_vs_Wuttke_bar_scalebreak_10k_8mil_v3_colour.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot.ckdgen)
dev.off()
	

# =============
# Repeat above for the BBJ metaCCA:
bbj_egfr <- read.table("bbj-a-60.eGFR.noheader.vcf", header=F)
names(bbj_egfr) <- c("chr","pos","rsid","ref","alt","dot","pass","af","info","info2")
# in info: ##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">

#names(bbj_egfr) <- c("chr","pos","rsid","ref","alt","dot","pass","af","info","EF","SE","LP","AF","ID")
bbj_egfr <- bbj_egfr %>% 
	separate(info2, into = c("EF","SE","LP","AF","ID"), sep = ":") %>%
	select(chr, pos, rsid, ref, alt, EF, SE, LP) %>%
	mutate(pvalue_egfr = 10^(as.numeric(LP)*-1))
 
bbj_bun <- read.table("bbj-a-11.BUN.noheader.vcf", header=F)
names(bbj_bun) <- c("chr","pos","rsid","ref","alt","dot","pass","af","info","info2")

bbj_bun <- bbj_bun %>% 
	separate(info2, into = c("EF","SE","LP","AF","ID"), sep = ":") %>%
	select(chr, pos, rsid, ref, alt, EF, SE, LP) %>%
	mutate(pvalue_bun = 10^(as.numeric(LP)*-1))
	
bbj_sig_all <- read.table("metaCCA_bbj_univar_eGFR_BUN_noPrune_sigBonf5507snps_w_chrpos.txt", header=T)
bbj_sig_all <- bbj_sig_all[bbj_sig_all$pBonf < 0.05,]
# 5507
names(bbj_sig_all) <- c("rsid","r_1","log10_pval","pval","pBH","pBonf","chr","pos","allele_0","allele_1","eGFR_b","eGFR_se","BUN_b","BUN_se")

# no. SNPs sig by nmetaCCA and in comp dir for kidney function:
bbj_sig_all_dir <- bbj_sig_all %>% 
	mutate(egfr_dir = ifelse(eGFR_b<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(BUN_b<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-"))
# 3712

# No. SNPs ID'ed by Kanai for both eGFR (5e-8) and BUN (5e-8) + dir:
kanai_egfr_bun <- merge(bbj_bun, bbj_egfr, by = c("rsid"))
kanai_egfr_bun_sig_dir <- kanai_egfr_bun %>% 
	mutate(egfr_dir = ifelse(EF.y<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(EF.x<0, "-", "+")) %>%
	filter((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-")) %>%
	filter(pvalue_bun < 5e-8 & pvalue_egfr < 5e-8)
	
length(unique(kanai_egfr_bun_sig_dir$rsid))
# 1241

# How many were sig. in Kanai et al., 2018 paper? see Results: "Overall, we identified 1,407 trait-associated loci for 53 quantitative traits that satisfied a genome-wide 
#  significance threshold of P = 5.0 × 10-8 ")
bbj_egfr_pBonf.bbj_sig_all <- bbj_egfr[bbj_egfr$rsid %in% bbj_sig_all$rsid & bbj_egfr$pvalue < 5e-8,]
# 
bbj_sig_all.kanai_egfr_bun_sig_dir <- bbj_sig_all[bbj_sig_all$RSID %in% kanai_egfr_bun_sig_dir$rsid,]
# 1241

# only common AF SNPs:
commonAF <- read.table("bbj-a-60.eGFR.bbj-a-11bun.noheader.commonOnly.5973670.txt", header=T)
kanai_egfr_bun <- kanai_egfr_bun[kanai_egfr_bun$rsid %in% commonAF$rsid,]
# 5973670
kanai_egfr_bun_sig_dir <- kanai_egfr_bun_sig_dir[kanai_egfr_bun_sig_dir$rsid %in% commonAF$rsid,]
# 1241
bbj_sig_all <- bbj_sig_all[bbj_sig_all$rsid %in% commonAF$rsid,]
# 5507

kanai_egfr_bun_sig_dir$kanai_egfr_bun_sig_dir <- "yes"
kanai_egfr_bun$kanai_egfr_bun <- "yes"

# plot comparing these with metaCCA out
kanai_egfr_bun.kanai_egfr_bun_sig_dir <- merge(kanai_egfr_bun[,c(1,18)], kanai_egfr_bun_sig_dir[,c(1,20)], by = c("rsid"),all.x=T)
bbj_sig_all.kanai_egfr_bun_sig_dir <- merge(kanai_egfr_bun.kanai_egfr_bun_sig_dir, bbj_sig_all, by = c("rsid"),all.x=T)

rsid <- c("rs1533988","rs4715491","rs6026578","rs79105258","rs963837","rs12509595","rs1275609","rs549752")
kanai_leadSNPs <- data.frame(rsid)
kanai_leadSNPs$kanai_leadSNPs <- "yes"

bbj_sig_all.kanai_egfr_bun_sig_dir.lead <- merge(bbj_sig_all.kanai_egfr_bun_sig_dir, kanai_leadSNPs, by = c("rsid"),all.x=T)

bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_egfr_bun_sig_dir <- 
	ifelse(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_egfr_bun_sig_dir %in% "yes", "yes","no")

bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_leadSNPs <- 
	ifelse(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_leadSNPs %in% "yes", "yes","no")

# create one vari for lead, no-lead but sig, and ns. in Wuttke
bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat <- 
	ifelse(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_leadSNPs %in% "yes", "lead",
		ifelse(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_egfr_bun_sig_dir %in% "yes", "not lead but sig.",
			"not sig."))

bbj_sig_all.kanai_egfr_bun_sig_dir.lead <- bbj_sig_all.kanai_egfr_bun_sig_dir.lead %>% 
	mutate(egfr_dir = ifelse(eGFR_b<0, "-", "+")) %>%
	mutate(bun_dir = ifelse(BUN_b<0, "-", "+")) %>%
	mutate(metacca_sig_dir = ifelse((((egfr_dir== "-" & bun_dir=="+") | (egfr_dir== "+" & bun_dir=="-")) & !is.na(pval)),"yes","no"))

Dataset_noXY <- read.table("bbj-a-60.eGFR.bbj-a-11.BUN.noheader.noXY.txt",header=F)
bbj_sig_all.kanai_egfr_bun_sig_dir.lead <- 
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead[bbj_sig_all.kanai_egfr_bun_sig_dir.lead$rsid %in% Dataset_noXY$V3,]
# 5837593

forplot <- data.frame(kanai_cat = character(), metacca_sig_dir = character(), count_SNPs = numeric(), stringsAsFactors=FALSE)

forplot[1,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[1]
forplot[2,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[2]
forplot[3,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[3]
forplot[4,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[1]
forplot[5,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[2]
forplot[6,1] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat)[3]

forplot[1,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot[2,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot[3,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[1]
forplot[4,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]
forplot[5,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]
forplot[6,2] <- unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir)[2]

forplot[1,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "not sig." &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot[2,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "not lead but sig." &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot[3,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "lead" &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "no",1]))
		
forplot[4,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "not sig." &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))
		
forplot[5,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "not lead but sig." &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))

forplot[6,3] <- length(unique(bbj_sig_all.kanai_egfr_bun_sig_dir.lead[
	bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "lead" &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$metacca_sig_dir %in% "yes",1]))

notleadbutsig <- bbj_sig_all.kanai_egfr_bun_sig_dir.lead[bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_cat %in% "not lead but sig." &
		bbj_sig_all.kanai_egfr_bun_sig_dir.lead$kanai_egfr_bun_sig_dir %in% "yes",]

# rs7079481C>A leading to MYPN p.Pro1135Thr
notleadbutsig[notleadbutsig$rsid %in% c("rs7079481"),]
# none
bbj_sig_all.kanai_egfr_bun_sig_dir.lead[bbj_sig_all.kanai_egfr_bun_sig_dir.lead$rsid == "rs7079481",]

forplot$kanai_cat_toshow <- ifelse(forplot$kanai_cat %in% c("not sig."), "1. Not significant and/or effect sizes not compatible (n = 5,836,352)",
	ifelse(forplot$kanai_cat %in% c("lead"), "3. Significant, effect sizes compatible and lead (n = 8)", "2. Significant, effect sizes compatible but not lead (n = 1,233)"))
	
forplot$metacca_sig_dir_toshow <- ifelse(forplot$metacca_sig_dir %in% c("yes"), "Yes\n(n = 3,712)", "No\n(n = 5,833,881)")
	
plot <- ggplot(forplot, 
			aes(x = stringr::str_wrap(kanai_cat_toshow,20), y = count_SNPs, fill = metacca_sig_dir_toshow)) +
		geom_bar(position="dodge", stat="identity")+
		geom_text(aes(label = count_SNPs, y = ifelse(count_SNPs<10, count_SNPs+80, count_SNPs-100)), colour = "black", size = 2.5, position = position_dodge(.9))+
		#vjust = 1.5, 
		#scale_fill_grey(start = 0.8, end = 0.6, name="SNP identified\nby metaCCA")+
		scale_fill_discrete(name="Identified by metaCCA?")+
		scale_y_break(c(2500, 5833000)) + 
		scale_y_continuous(breaks = seq(0, 5833500, by = 500), position = "left")+
		xlab("Outcome by Kanai et al., 2018 genome-wide association study")+ 
        ylab("Number of single\nnucleotide polymorphisms")+
		labs(tag="B")+
    theme(
    #legend.position = "none",
	legend.spacing = unit(0.7, 'cm'),
	legend.text=element_text(family="Arial",face="bold", size = 8,color="black"),
	legend.title=element_text(family="Arial",face="bold", size = 8,color="black"),
	panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.3),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",angle=90),
	  #margin = margin(t = 0, r = 12, b = 0, l = 0)),
    #plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    plot.tag = element_text(family="Arial",face="bold", size = 20,color="black")
)

jpeg("plot_BBJ_metaCCA_vs_Kanai_bar_scalebreak_noXY_common_only_v2_colour.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot)
dev.off()

# to display both as A and B in same plot:
jpeg("plot_A_CKDGen_metaCCA_vs_Wuttke_bar_scalebreak_10k_8mil_v3_colour_with_BBJ_as_B.jpeg", pointsize=6, width=170, height=150,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
aplot::plot_list(plot.ckdgen,plot, ncol=1, nrow=2)
dev.off()