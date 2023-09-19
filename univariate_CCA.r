source("gene-based.cca3.2.R")

library(dplyr)
library(ggplot2)
library(grid)
require(scales)
library(sqldf)
library(stringr)
library(SuperExactTest)

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
	

# load up phenD data (egfr and bun have been adjusted for age and gender)
phenD <- read.table("NCKD_egfr_bun_baseline_order_diagdate_wo_ethnicity_w_ctrls_age_sex_adj_forCCA.csv", 
	header=T)
	
# cols needed: patient_id, case_control, egfr_adj_scaled, urea_adj_scaled, SentrixPosition
phenD <- phenD[,c(1,8,25,24,12)] 

row.names(phenD) <- phenD$SentrixPosition

for (j in 1:22){
	print(paste0("Chromosome: ", j))
	inputData <- read.table(paste0("chr", j, ".GT_data.alt_counts_nomissing.txt"),header=T,skip=6)
	print(paste0("Data loaded."))
	inputData <- inputData[!duplicated(inputData$ID),]
	inputData_torun <- inputData[,-c(1:2,4:9)]
	rownames(inputData_torun) <- inputData_torun$ID
	inputData_torun <- t(inputData_torun[,-c(1)])
	inputData_torun <- inputData_torun[rownames(inputData_torun) %in% phenD$SentrixPosition,]
	phenD_torun <- phenD[phenD$SentrixPosition %in% rownames(inputData_torun),c(3,4,2)]
	snps = inputData$ID
	l.snps = length(snps)
	geneD = as.matrix(inputData_torun)
	if(length(which(apply(geneD, 2, var) == 0))!=0){
		geneD <- geneD[ -as.numeric(which(apply(geneD, 2, var) == 0))]
	}
	nind=nrow(geneD)
	print(paste0("nind: ", nind))
	nsnps=1
	ntraits=3
	
	print('Running CCA:')
	class(geneD) <- "numeric"
	all.can <- as.matrix(apply(geneD, MARGIN = 2, cca))
	p <- apply(all.can, MARGIN = 1, fapp_func)
	
	results <- merge(all.can, as.matrix(p), by=0, all=TRUE)
	names(results) <- c('RSID','CCA_r','CCA_pval')
	write.table(results, 
		paste0("CCA_out_chr",j,".txt"),col.names=T,row.names=F,quote=F)
	
}
	
	
# Load in all results and visualise:
for (j in 1:22){
	results <- read.table(paste0("CCA_out_chr",j,".txt"),header=T)
	results <- results[,c(1:3)]
	names(results) <- c("SNP","CCA_r","CCA_pval")
	if(j==1){
		all.results <- results
	} 
	else{
		all.results <- rbind(all.results, results)
	}
}

chrposrefalt <- read.table('SNP_info_cols_chr_pos_ref_alt.txt', header=T)
names(chrposrefalt) <- c("chr","pos","SNP","ref","alt")

# Manhattan plot:
data_all <- merge(all.results, chrposrefalt, by = c("SNP"), all.x=T)
data_sig <- data_all[data_all$CCA_pval<0.05,]
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
ylim <- abs(floor(log10(min((gwas.dat$CCA_pval*5e-8))))) + 2

gwas.dat$pval_minuslog10 <- -log10(gwas.dat$CCA_pval)

sig <- 0.05

max(gwas.dat$pval_minuslog10)

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
  scale_y_continuous(breaks = c(seq(0, 14, by = 2)),expand=c(0.04, 0)) +
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

jpeg("manh_plot.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()


# ======== FOR VENN:
# add in SKS CCA results:
data_sks <- read.table("SKS_CCA_out.txt", header=T)

# 1. NCKD and SKS:
nckd_5e_8 <- data_all[data_all$CCA_pval < 5e-8,]
sks_p0.05 <- data_sks[data_sks$CCA_pval < 0.05,]
nckd_p0.05 <- data_all[data_all$CCA_pval < 0.05,]
nckd_5e_8.sks_p0.05 <- nckd_5e_8[nckd_5e_8$SNP %in% sks_p0.05$SNP,]

sks_5e_8 <- data_sks[data_sks$CCA_pval < 5e-8,]
nckd_5e_8.sks_5e_8 <- nckd_5e_8[nckd_5e_8$SNP %in% sks_5e_8$SNP,]
sks_5e_8.nckd_p0.05 <- sks_5e_8[sks_5e_8$SNP %in% nckd_p0.05$SNP,]

# 2. NCKD and CKDGen:
all_ckdgen <- read.table("CKDGen_metaCCA_out.txt", header=T)

# all ckdgen in all nckd
all_ckdgen.nckd <- all_ckdgen[all_ckdgen$RSID %in% data_all$SNP,]

# all ckdgen in nckd 5e-8
all_ckdgen.nckd_5e_8 <- all_ckdgen[all_ckdgen$RSID %in% nckd_5e_8$SNP,]

# ckdgen 5e-8 in nckd 5e-8
all_ckdgen_5e_8.nckd_5e_8 <- all_ckdgen.nckd_5e_8[all_ckdgen.nckd_5e_8$pval<5e-8,]

# ckdgen 5e-8 in nckd p0.05
all_ckdgen_5e_8.nckd_p0.05 <- all_ckdgen[all_ckdgen$RSID %in% nckd_0.05$SNP,]

# ckdgen 0.05 in nckd 5e-8
all_ckdgen_0.05.nckd_5e_8 <- all_ckdgen.nckd_5e_8[all_ckdgen.nckd_5e_8$pval<0.05,]

# all ckdgen in nckd 5e-8:
ckdgen.nckd_5e_8 <- all_ckdgen[all_ckdgen$RSID %in% nckd_5e_8$SNP,]

# ckdgen.14045, sig, dir
ckdgen.14045 <- read.table("CKDGen_metaCCA_out_14045SNPs.txt", header=F)
names(ckdgen.14045) <- c("loc", "gene", "chr", "pos_from", "pos_to", "ref", "alt", "SNP")

# ckdgen 5e-8 in nckd 5e-8
ckdgen.14045.nckd_5e_8 <- ckdgen.14045[ckdgen.14045$SNP %in% nckd_5e_8$SNP,]

# ckdgen 5e-8 in nckd p0.05
ckdgen.14045.nckd_0.05 <- ckdgen.14045[ckdgen.14045$SNP %in% nckd_p0.05$SNP,]



# 3. all ckdgen in all sks
all_ckdgen.sks <- all_ckdgen[all_ckdgen$RSID %in% data_sks$SNP,]

# all ckdgen in sks 5e-8
all_ckdgen.sks_5e_8 <- all_ckdgen[all_ckdgen$RSID %in% sks_5e_8$SNP,]

# ckdgen 5e-8 in sks 5e-8
all_ckdgen_5e_8.sks_5e_8 <- all_ckdgen.sks_5e_8[all_ckdgen.sks_5e_8$pval<5e-8,]

# ckdgen 5e-8 in sks p0,05
all_ckdgen_5e_8.sks_p0.05 <- all_ckdgen[all_ckdgen$RSID %in% sks_0.05$SNP,]

# ckdgen 0.05 in nckd 5e-8
all_ckdgen_0.05.sks_5e_8 <- all_ckdgen.sks_5e_8[all_ckdgen.sks_5e_8$pval<0.05,]

# ckdgen 5e-8 in sks p0,05
ckdgen.14045.sks_0.05 <- ckdgen.14045[ckdgen.14045$SNP %in% sks_p0.05$SNP,]

# ckdgen 5e-8 in sks 5e-8
ckdgen.14045.sks_5e_8 <- ckdgen.14045[ckdgen.14045$SNP %in% sks_5e_8$SNP,]


# STATISTICAL OVERLAPS:
# ckdgen p<5e-8, nckd p<0.05
list.2 <- list('CKDGen.p5e_8' = ckdgen.14045$SNP, 'NURTuRE-CKD.p0.05' = nckd_0.05$SNP)
Result.2 = supertest(list.2, n=length(all_ckdgen.nckd$RSID))
# n.s. overlap

Results_table.2 <- rbind(Result.2$set.names, Result.2$set.sizes, Result.2$overlap.sizes, Result.2$overlap.expected, Result.2$P.value)

# ckdgen p<5e-8, sks p<5e-8
list.sks.5e_8 <- list('CKDGen.p5e_8' = ckdgen.14045$SNP, 'SKS.p0.05' = sks_0.05$SNP)
Result.sks.5e_8 = supertest(list.sks.5e_8, n=length(all_ckdgen.sks$RSID))
# n.s. overlap


# Try all three datasets:
all_ckdgen.sks.nckd <- all_ckdgen[all_ckdgen$RSID %in% data_all$SNP & all_ckdgen$RSID %in% data_sks$SNP,]
all_ckdgen5e_8.sks0.05.nckd0.05 <- ckdgen.14045.sks_0.05[ckdgen.14045.sks_0.05$SNP %in% ckdgen.14045.nckd_0.05$SNP,]

# ckdgen 14045, sks and nckd p0.05
list.ckdgen.nckd.sks <- list('CKDGen.p5e_8' = ckdgen.14045$SNP, 'NCKD.p0.05' = nckd_0.05$SNP, 'SKS.p0.05' = sks_0.05$SNP)
Result.ckdgen.nckd.sks = supertest(list.ckdgen.nckd.sks, n=length(all_ckdgen.sks.nckd$RSID))

