# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.

library(ggplot2)
library(grid)
library(sqldf)
library(dplyr)
require(scales)
library(stringr)
library(ggrepel)
library(tidyr)
`%notin%` <- Negate(`%in%`)
sig <- 0.05

# 1. UKHLS (Fig. S2)
data_all <- read.table("univariate_results_2nddataset.txt", header=T)
names(data_all) <- c("SNP","r","pval_minuslog10","pval","pBH","pBonf", "chr","pos")
data_all$chr <- as.numeric(data_all$chr)

data_sig <- data_all[data_all$pBH<0.05,]
data_nonsig <- data_all[data_all$SNP %notin% data_sig$SNP,]

data_sig <- data_sig[ order(data_sig$chr),]
data_nonsig <- data_nonsig[ order(data_nonsig$chr),]

gwas.dat2 <- rbind(data_sig,data_nonsig)
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

notsig.dat <- data_nonsig %>%
  slice(sample(nrow(.), nrow(.) / 10))
gwas.dat <- rbind(data_sig,data_nonsig)
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
ylim <- abs(floor(log10(min(gwas.dat$pBH)))) + 2

gwas.dat$pBH_minuslog10 <- -log10(gwas.dat$pBH)

# Add SNPs found to be significiant by Wuttke et al., 2019 "Nature_sig"
# also use the metaCCA p-val so that they can be plotted.
Nature_sig <- read.table("wuttke_kidrel1_121snps_mvp.txt", header=F)
names(Nature_sig) <- c("SNP","chr","pos","loci")
nature_SNPs <- as.character(Nature_sig$SNP)
pvals130 <- gwas.dat2[gwas.dat2$SNP %in% nature_SNPs,]
Nature_sig <- merge(Nature_sig, pvals130, by=c("SNP"))
Nature_sig <- Nature_sig[,c(1:8,11:12)]
names(Nature_sig) <- c("rsid","chr","pos","loci","r","pval_minuslog10","pval","pBonf","pBH","BPcum")
Nature_sig$pBH_minuslog10 <- -log10(Nature_sig$pBH)

# for legend symbols:
Nature_sig[(nrow(Nature_sig)+1),10] <- min(gwas.dat[gwas.dat$chr==1,9])+((max(gwas.dat[gwas.dat$chr==1,9]))/5)
Nature_sig[(nrow(Nature_sig)),2] <- 1
Nature_sig[(nrow(Nature_sig)),9] <- 10^-(4.2)
Nature_sig[(nrow(Nature_sig)),11] <- 5.34

# Also import data on hg19 coords for the genes also found to be significiant by Wuttke et al., 2019:
Nature_sig_genes <- read.table("10genes_ukhls_wuttke_hg19.txt", header=F)
names(Nature_sig_genes) <- c("chr","pos_from","pos_to","gene","loci_wuttke")
gwas.dat$pos <- as.numeric(gwas.dat$pos)
Nature_sig_genes$pos_to <- as.numeric(Nature_sig_genes$pos_to)
Nature_sig_genes$pos_from <- as.numeric(Nature_sig_genes$pos_from)

# for each gene range, get SNP with max(pBH_minslog10)
for (i in 1:nrow(Nature_sig_genes)){
  range1 <- Nature_sig_genes$pos_from[i]
  range2 <- Nature_sig_genes$pos_to[i]
  chr <- Nature_sig_genes$chr[i]
  max_pBHlog <- max(gwas.dat[(gwas.dat$pos>=range1 & gwas.dat$pos<=range2 & gwas.dat$chr==chr),10])
  SNP_needed <- gwas.dat[(gwas.dat$pos>=range1 & gwas.dat$pos<=range2 & gwas.dat$chr==chr & gwas.dat$pBH_minuslog10 == max_pBHlog),1][1]
  Nature_sig_genes$pBH_minuslog10[i] <- max_pBHlog
  BPcum <- gwas.dat[gwas.dat$SNP == SNP_needed,9]
  if(BPcum)
  Nature_sig_genes$BPcum[i] <- gwas.dat[gwas.dat$SNP == SNP_needed,9]
  Nature_sig_genes$SNP[i] <- SNP_needed
}

# add extra points to nature_sig and nature_sig_genes datasets that get displayed as the legend
Nature_sig_genes[(nrow(Nature_sig_genes)+1),1] <- 1
Nature_sig_genes[(nrow(Nature_sig_genes)),7] <- min(gwas.dat[gwas.dat$chr==1,9])+((max(gwas.dat[gwas.dat$chr==1,9]))/5)
Nature_sig_genes[(nrow(Nature_sig_genes)),6] <- 5.11

box_xmax <- (min(gwas.dat[gwas.dat$chr==4,9]))+((max(gwas.dat[gwas.dat$chr==4,9])-(min(gwas.dat[gwas.dat$chr==4,9])))/3)

gwas.dat.ukhls <- gwas.dat
gwas.dat.ukhls$dataset <- "ukhls"

manhplot_ukhls <- ggplot(gwas.dat, aes(x = BPcum, y = (pBH_minuslog10), color = as.factor(chr), group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "dashed") +
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==2,3]), xmax = max(gwas.dat[gwas.dat$chr==2,3]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==4,3]), xmax = max(gwas.dat[gwas.dat$chr==4,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==6,3]), xmax = max(gwas.dat[gwas.dat$chr==6,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==8,3]), xmax = max(gwas.dat[gwas.dat$chr==8,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==10,3]), xmax = max(gwas.dat[gwas.dat$chr==10,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==12,3]), xmax = max(gwas.dat[gwas.dat$chr==12,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==14,3]), xmax = max(gwas.dat[gwas.dat$chr==14,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==16,3]), xmax = max(gwas.dat[gwas.dat$chr==16,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==18,3]), xmax = max(gwas.dat[gwas.dat$chr==18,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==20,3]), xmax = max(gwas.dat[gwas.dat$chr==20,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==22,3]), xmax = max(gwas.dat[gwas.dat$chr==22,3]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==1,9]), xmax = box_xmax, ymin = 4.9,ymax =6, fill = "white", color="black")+
  geom_point(data=Nature_sig, shape = 2,size=3, color="black") +
  geom_point(data=Nature_sig_genes, size=3, color="black", shape=0) +
  annotate(geom="text", x=(min(gwas.dat[gwas.dat$chr==1,9])+((max(gwas.dat[gwas.dat$chr==1,9]))/9)), y=5.5, label="Previously reported as \nsignificant by CKDGen GWAS:\n      SNP\n      Gene",
              color="black", family="Arial",fontface="bold", size=7/(14/5), hjust = 0, fill = "white") +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 6), breaks = c(seq(0, 5.5, by = 0.5)),expand=c(0.04, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosomal location",
       y = expression(bold(-log[10](metaCCA~bolditalic(P)~value[adjusted])))) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 7,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 7,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 7,color="black")
)

tiff("manh_plot_Fig_S1.tiff", pointsize=6, width=170, height=108,units="mm",res = 294,compression = "lzw")
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_ukhls, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()


# 2. CKDGen (Fig. S1):
data_all <- read.table("univariate_results_1stdataset.txt", header=T)
names(data_all) <- c("SNP","r","pval_minuslog10","pval","pBH","chr","pos","allele1","allele2","freq1","eff_egfr","se_egfr","pval_egfr","ntotalsum")

data_sig <- data_all[data_all$pBH<0.05,]
data_nonsig <- data_all[data_all$SNP %notin% data_sig$SNP,]
data_sig <- data_sig[ order(data_sig$chr),]
data_nonsig <- data_nonsig[ order(data_nonsig$chr),]

gwas.dat2 <- rbind(data_sig,data_nonsig)
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

gwas.dat <- rbind(data_sig,data_nonsig)
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
ylim <- abs(floor(log10(min(gwas.dat$pBH)))) + 2
gwas.dat$pBH_minuslog10 <- -log10(gwas.dat$pBH)

Nature_sig_raw <- read.table("wuttke_kidrel1_121snps_mvp.txt", header=F)
names(Nature_sig_raw) <- c("SNP","chr","pos","loci")
nature_SNPs <- as.character(Nature_sig_raw$SNP)
pvals130 <- gwas.dat[gwas.dat$SNP %in% nature_SNPs,]

Nature_sig_genes <- read.table("30genes_ckdgen_wuttke_hg19.txt", header=F)
names(Nature_sig_genes) <- c("chr","pos_from","pos_to","gene","loci_wuttke")
gwas.dat$pos <- as.numeric(gwas.dat$pos)
Nature_sig_genes$pos_to <- as.numeric(Nature_sig_genes$pos_to)
Nature_sig_genes$pos_from <- as.numeric(Nature_sig_genes$pos_from)

# for each gene range, get SNP with max(pBH_minslog10)
for (i in 1:nrow(Nature_sig_genes)){
  range1 <- Nature_sig_genes$pos_from[i]
  range2 <- Nature_sig_genes$pos_to[i]
  chr <- Nature_sig_genes$chr[i]
  Nature_sig_genes$pBH_minuslog10[i] <- max(gwas.dat[(gwas.dat$pos>=range1 & gwas.dat$pos<=range2 & gwas.dat$chr==chr),16])
  max_pBHlog <- max(gwas.dat[(gwas.dat$pos>=range1 & gwas.dat$pos<=range2 & gwas.dat$chr==chr),16])
  SNP_needed <- gwas.dat[(gwas.dat$pos>=range1 & gwas.dat$pos<=range2 & gwas.dat$chr==chr & gwas.dat$pBH_minuslog10==  max_pBHlog),1]
  Nature_sig_genes$BPcum[i] <- gwas.dat[gwas.dat$SNP == SNP_needed,15]
  Nature_sig_genes$SNP[i] <- SNP_needed
}

# Legend - add points for display
x_leg <- (max(gwas.dat[gwas.dat$chr==2,15])-min(gwas.dat[gwas.dat$chr==2,15]))/10

Nature_sig <- merge(Nature_sig_raw, pvals130, by=c("SNP","chr","pos"))
Nature_sig[(nrow(Nature_sig)+1),2] <- 2
Nature_sig[(nrow(Nature_sig)),16] <- min(gwas.dat[gwas.dat$chr==2,15])+x_leg
Nature_sig[(nrow(Nature_sig)),17] <- 86

Nature_sig_genes[(nrow(Nature_sig_genes)+1),1] <- 2
Nature_sig_genes[(nrow(Nature_sig_genes)),7] <- min(gwas.dat[gwas.dat$chr==2,15])+x_leg
Nature_sig_genes[(nrow(Nature_sig_genes)),6] <- 91

manhplot_ckdgen <- ggplot(gwas.dat, aes(x = BPcum, y = (pBH_minuslog10), color = as.factor(chr),group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6,show.legend = F) +
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==2,15]), xmax = max(gwas.dat[gwas.dat$chr==2,15]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==4,15]), xmax = max(gwas.dat[gwas.dat$chr==4,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==6,15]), xmax = max(gwas.dat[gwas.dat$chr==6,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==8,15]), xmax = max(gwas.dat[gwas.dat$chr==8,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==10,15]), xmax = max(gwas.dat[gwas.dat$chr==10,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==12,15]), xmax = max(gwas.dat[gwas.dat$chr==12,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==14,15]), xmax = max(gwas.dat[gwas.dat$chr==14,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==16,15]), xmax = max(gwas.dat[gwas.dat$chr==16,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==18,15]), xmax = max(gwas.dat[gwas.dat$chr==18,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==20,15]), xmax = max(gwas.dat[gwas.dat$chr==20,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==22,15]), xmax = max(gwas.dat[gwas.dat$chr==22,15]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = (min(gwas.dat[gwas.dat$chr==2,15])-(x_leg*2)), xmax = (min(gwas.dat[gwas.dat$chr==8,15])+(x_leg*4)), ymin = 82.5,ymax =98, fill = "white", color="black")+
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "dashed") +
  geom_point(data=Nature_sig, shape = 2,size=3, color="black") +
  geom_point(data=Nature_sig_genes, size=3, color="black", shape=0) +
  annotate(geom="text", x=249026048, y=91, label="Previously reported as significant by CKDGen GWAS:\n      Gene\n      SNP",
              color="black", family="Arial",fontface="bold", size=7/(14/5), hjust = 0, fill = "white") +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 99), breaks = c(seq(0, 95, by = 5)),expand=c(0.04, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +

  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosomal location",
     y = expression(bold(-log[10](metaCCA~bolditalic(P)~value[adjusted])))) +
  theme(
  	legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 7,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 7,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 7,color="black")
  )


# 3. BBJ (Fig. S3):
data_all <- read.table("univariate_results_3rddataset.txt", header=T)
names(data_all) <- c("SNP","r","pval_minuslog10","pval","pBH","chr","pos","allele1","allele2","egfr_b","egfr_se","bun_b","bun_se","gene")
data_sig <- data_all[data_all$pBH<0.05,]
data_nonsig <- data_all[data_all$SNP %notin% data_sig$SNP,]
data_sig <- data_sig[ order(data_sig$chr),]
data_nonsig <- data_nonsig[ order(data_nonsig$chr),]

gwas.dat2 <- rbind(data_sig,data_nonsig)
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

gwas.dat <- rbind(data_sig,data_nonsig)
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
ylim <- abs(floor(log10(min(gwas.dat$pBH)))) + 2
gwas.dat$pBH_minuslog10 <- -log10(gwas.dat$pBH)

Nature_sig_raw <- read.table("wuttke_kidrel1_122snps_mvp.txt", header=F)
names(Nature_sig_raw) <- c("SNP","chr","pos","loci")
nature_SNPs <- as.character(Nature_sig_raw$SNP)
pvals130 <- gwas.dat[gwas.dat$SNP %in% nature_SNPs,]

allgenesstudy <- sort(unique(data_all$gene))
Nature_sig_genes_79 <- array(allgenesstudy[allgenesstudy %in% Nature_sig_raw$loci])
Nature_sig_genes <- Nature_sig_raw %>% filter(grepl(paste(Nature_sig_genes_79, collapse="|"), loci))

gwas.dat$pos <- as.numeric(gwas.dat$pos)
Nature_sig_genes$pos <- as.numeric(Nature_sig_genes$pos)

# for each gene, get SNP with max(pBH_minslog10)
for (i in 1:nrow(Nature_sig_genes)){
  gene <- Nature_sig_genes$loci[i]
  all_SNPs <- gwas.dat[gwas.dat$gene %in% gene,c(1,5,6,14,15,16)]
  max_pBHminuslog10 <- max(all_SNPs$pBH_minuslog10)
  SNP <-all_SNPs[all_SNPs$pBH_minuslog10 == max_pBHminuslog10,1]
  BPcum <- all_SNPs[all_SNPs$SNP == SNP,5]
  Nature_sig_genes$BPcum[i] <- BPcum
  Nature_sig_genes$SNP_bbj[i] <- SNP
  Nature_sig_genes$pBH_minuslog10[i] <- all_SNPs[all_SNPs$SNP==SNP,6]
}
x_leg <- (max(gwas.dat[gwas.dat$chr==2,15])-min(gwas.dat[gwas.dat$chr==2,15]))/10

Nature_sig <- merge(Nature_sig_raw, pvals130, by=c("SNP","chr","pos"))
Nature_sig[(nrow(Nature_sig)+1),2] <- 2
Nature_sig[(nrow(Nature_sig)),16] <- min(gwas.dat[gwas.dat$chr==2,15])+x_leg
Nature_sig[(nrow(Nature_sig)),17] <- 39.7

Nature_sig_genes[(nrow(Nature_sig_genes)+1),2] <- 2
Nature_sig_genes[(nrow(Nature_sig_genes)),5] <- min(gwas.dat[gwas.dat$chr==2,15])+x_leg
Nature_sig_genes[(nrow(Nature_sig_genes)),7] <- 42

gwas.dat.ckdgen <- gwas.dat
gwas.dat.ckdgen$dataset <- "ckdgen"
gwas.dat.ckdgen <- gwas.dat.ckdgen[,c(1:7,15:17)]
gwas.dat.ukhls <- gwas.dat.ukhls[,c(1:5,7:11)]
gwas.dat.bbj <- gwas.dat[,c(1:7,15:16)]
gwas.dat.bbj$dataset <- "bbj"
gwas.dat.both <- rbind(gwas.dat.ckdgen,gwas.dat.ukhls)
gwas.dat.all3 <- rbind(gwas.dat.ckdgen,gwas.dat.ukhls,gwas.dat.bbj)

manhplot_both <- ggplot(gwas.dat.both, aes(x = BPcum, y = (pBH_minuslog10), color = as.factor(chr),group = as.factor(chr))) +
geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6,show.legend = F) +
geom_hline(yintercept = -log10(sig), color = "black", linetype = "dashed") +
geom_point(data=Nature_sig, shape = 2,size=3, color="black") +
geom_point(data=Nature_sig_genes, size=3, color="black", shape=0) +
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==2,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==2,8]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==4,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==4,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==6,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==6,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==8,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==8,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==10,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==10,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==12,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==12,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==14,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==14,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==16,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==16,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==18,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==18,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==20,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==20,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==22,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==22,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
scale_size_continuous(range = c(0.5,3)) +
labs(x = "Chromosomal location",
     y = expression(bold(-log[10](metaCCA~bolditalic(P)~value[adjusted])))) +
facet_wrap(~ dataset, nrow=2, ncol=1, scales = "free_y") +
theme(
	legend.position = "none",
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.background = element_rect(fill="white"),
  panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
  axis.text.x = element_text(family="Arial", face="bold", size = 6.5, vjust = 0.5,color="black"),
  axis.text.y = element_text(family="Arial",face="bold", size = 7,color="black"),
  axis.title.x = element_text(family="Arial",face="bold", size = 7,color="black"),
  axis.title.y = element_text(family="Arial",face="bold", size = 7,color="black")
)

new_labels <- c("ckdgen" = "CKDGen", "bbj" = "BioBank Japan", "ukhls" = "UKHLS")

manhplot_all3 <- ggplot(gwas.dat.all3, aes(x = BPcum, y = (pBH_minuslog10), color = as.factor(chr),group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6,show.legend = F) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "dashed") +
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==2,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==2,8]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==4,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==4,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==6,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==6,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==8,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==8,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==10,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==10,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==12,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==12,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==14,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==14,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==16,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==16,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==18,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==18,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==20,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==20,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat.both[gwas.dat.both$chr==22,8]), xmax = max(gwas.dat.both[gwas.dat.both$chr==22,8]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(tag = "b", x = "Chromosomal location",
       y = expression(bold(-log[10](bolditalic(metaCCA)~bolditalic(P)~value[adjusted])))) +
  facet_wrap(~ dataset, nrow=3, ncol=1, scales = "free_y", labeller = labeller(dataset = new_labels)) +
  theme(
  	legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(face="bold", size = 6, vjust = 0.5,color="black"),
	  axis.text.y = element_text(face="bold", size = 8,color="black"),
	  axis.title.x = element_text(face="bold", size = 8,color="black"),
	  axis.title.y = element_text(face="bold", size = 8,color="black"),
    plot.tag = element_text(face="bold", size = 20,color="black")
  )

  ann_text_ckdgen <- data.frame(BPcum = 249026048, pBH_minuslog10=70,dataset = "ckdgen")
  Nature_sig_bbj <- Nature_sig
  Nature_sig_bbj$dataset <- "bbj"
  Nature_sig_genes_bbj <- Nature_sig_genes
  Nature_sig_genes_bbj$dataset <- "bbj"
  #Nature_sig_ukhls <- Nature_sig
  #Nature_sig_ukhls$dataset <- "ukhls"
  Nature_sig_genes_ukhls <- Nature_sig_genes[c(1:9),]
  Nature_sig_genes_ukhls$dataset <- "ukhls"
  Nature_sig_ckdgen <- Nature_sig
  Nature_sig_ckdgen$dataset <- "ckdgen"
  Nature_sig_genes_ckdgen <- Nature_sig_genes
  Nature_sig_genes_ckdgen$dataset <- "ckdgen"

  Nature_sig_ckdgen[2,17] <- 60
  Nature_sig_genes_ckdgen[31,6] <- 70

# All three manh plots in one:
manhplot_all3_w_annot <- manhplot_all3 + geom_rect(data = data.frame(dataset = "ckdgen"), aes(xmin = 200469147, xmax = 1488259299, ymin = 44,ymax =84), fill="white", inherit.aes = FALSE,color="black")+
            geom_text(data = ann_text_ckdgen,label = "Previously reported as significant by CKDGen or BBJ GWAS:\n      Gene\n      SNP\n",
            color="black", family="Arial",fontface="bold", size=7/(14/5), hjust=0) +
            geom_point(data=Nature_sig_bbj[c(1:4),], shape = 2,size=3, color="black")+
            geom_point(data=Nature_sig_genes_bbj[c(1:79),], size=3, color="black", shape=0)+
            geom_point(data=Nature_sig_ckdgen, shape = 2,size=3, color="black")+
            geom_point(data=Nature_sig_genes_ckdgen, size=3, color="black", shape=0)+
            geom_point(data=Nature_sig_genes_ukhls, size=3, color="black", shape=0)

jpeg("manh_plot_all3.jpeg", pointsize=6, width=170, height=180,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_all3_w_annot, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()

# Add gene-based results as square shape:
gene_based_snps <- read.csv("genebased_forplot.csv",header=T)
gene_based_snps <- separate_rows(gene_based_snps,snp, sep=", ", convert = TRUE)
names(gene_based_snps) <- c("SNP","gene","dataset","pval_adj")
gwas.dat.all3_gene_based <- merge(gene_based_snps,gwas.dat.all3,by=c("SNP","dataset"),all.x=T)
gene_based_snps_info <- gwas.dat.all3_gene_based[,c(1:4,9:11)]
gene_based_snps_info$pBH_minuslog10 <- -log10(gene_based_snps_info$pval_adj)
sql <- sqldf("SELECT MAX(BPcum), SNP, dataset, gene, pval_adj, chr, pos, BPcum, pBH_minuslog10 FROM gene_based_snps_info WHERE dataset != 'nurture-ckd' GROUP BY gene, dataset")

row <- (nrow(sql)+1)
sql[row,6] <- 2
sql[row,8] <- 273304498
sql[row,9] <- 62
sql[row,3] <- "ckdgen"

legend <- gwas.dat.all3[1,]
legend[1,1] <- "for_legend"
legend[1,6] <- 2
legend[1,8] <- 273304498
legend[1,9] <- 70
legend[1,10] <- "ckdgen"

manhplot_all3_w_annot_genebased <- manhplot_all3 + geom_rect(data = data.frame(dataset = "ckdgen"), aes(xmin = 200469147, xmax = 1061558708, ymin = 56,ymax =82), fill="white", inherit.aes = FALSE,color="black")+
            geom_text(data = ann_text_ckdgen,label = "Analysis type\n      Univariate-SNP\n      Gene-based multivariate-SNP",color="black", fontface="bold", size=7/(14/5), hjust=0) +
            geom_point(data=sql, size=2, color="black", shape=5)+
            geom_point(data=legend, alpha = 0.75,stroke = 0, size=1.6, color="black", shape=16)

jpeg("manh_plot_all3_w_genebased.jpeg", pointsize=6, width=170, height=180,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_all3_w_annot_genebased, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()

pdf("manh_plot_all3_w_genebased.pdf", width=6, height=4.5)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_all3_w_annot_genebased, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()

# Add Fig 1A too (see script: manhplots_nurture-ckd_for_Fig1A.R)
jpeg("manh_plot_nurture_plus_all3_w_genebased.jpeg", pointsize=6, width=170, height=250,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1, heights = c(0.3, 1))))
print(manhplot_all3_w_annot_genebased, vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(manhplot_nurture, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()
