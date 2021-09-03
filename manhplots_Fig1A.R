# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.

library(ggplot2)
library(grid)
library(sqldf)
library(dplyr)
require(scales)
library(stringr)
`%notin%` <- Negate(`%in%`)

data_all <- read.table("univar_results_for_manhplot.txt", header=T)
names(data_all) <- c("SNP","r","pval","gene","pBH","pBonf", "chr","pos","ref","alt","NA","Chr","gene_pos_from","gene_pos_to","gene")
data_sig <- data_all[data_all$pval<0.01,]
data_all$chr <- as.numeric(data_all$chr)
data_sig$chr <- as.numeric(data_sig$chr)
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

notsig.dat <- data_nonsig[,c(1:3,5:10)] %>%
  slice(sample(nrow(.), nrow(.) / 10))
gwas.dat <- rbind(data_sig[,c(1:3,5:10)],data_nonsig[,c(1:3,5:10)])
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

gwas.dat$pval_minuslog10 <- -log10(gwas.dat$pval)
sig <- 0.01

box_xmax <- (min(gwas.dat[gwas.dat$chr==4,9]))+((max(gwas.dat[gwas.dat$chr==4,9])-(min(gwas.dat[gwas.dat$chr==4,9])))/3)

gwas.dat.nurture <- gwas.dat
gwas.dat.nurture$dataset <- "nurture-ckd"

# Add gene-based results to plot (as square shape):
gene_based_snps <- read.csv("genebased_forplot.csv",header=T)
gene_based_snps <- separate_rows(gene_based_snps,snp, sep=", ", convert = TRUE)
names(gene_based_snps) <- c("SNP","gene","dataset","pval_adj")
gwas.dat.gene_based <- merge(gene_based_snps[gene_based_snps$dataset=="nurture-ckd",],gwas.dat,by=c("SNP"),all.x=T)
gene_based_snps_info <- gwas.dat.gene_based[,c(1:4,9,13)]
gene_based_snps_info$pval_minuslog10 <- -log10(gene_based_snps_info$pval_adj)
sql <- sqldf("SELECT MAX(BPcum), SNP, dataset, gene, pval_adj, chr, BPcum, pval_minuslog10 FROM gene_based_snps_info WHERE dataset == 'nurture-ckd' GROUP BY gene, dataset")

manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = pval_minuslog10, color = as.factor(chr), group = as.factor(chr))) +
  geom_point(alpha = 0.75,stroke = 0, shape = 16,size=1.6) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "dashed") +
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==2,10]), xmax = max(gwas.dat[gwas.dat$chr==2,10]), ymin = 0, ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==4,10]), xmax = max(gwas.dat[gwas.dat$chr==4,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==6,10]), xmax = max(gwas.dat[gwas.dat$chr==6,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==8,10]), xmax = max(gwas.dat[gwas.dat$chr==8,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==10,10]), xmax = max(gwas.dat[gwas.dat$chr==10,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==12,10]), xmax = max(gwas.dat[gwas.dat$chr==12,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==14,10]), xmax = max(gwas.dat[gwas.dat$chr==14,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==16,10]), xmax = max(gwas.dat[gwas.dat$chr==16,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==18,10]), xmax = max(gwas.dat[gwas.dat$chr==18,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==20,10]), xmax = max(gwas.dat[gwas.dat$chr==20,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==22,10]), xmax = max(gwas.dat[gwas.dat$chr==22,10]), ymin = 0,ymax = Inf, fill = "gray65", alpha=0.2)+
  annotate("rect", xmin = min(gwas.dat[gwas.dat$chr==1,10]), xmax = box_xmax, ymin = 4.9,ymax =6, fill = "white", color="black")+
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
  scale_y_continuous(limits = c(0, 5), breaks = c(seq(0, 5, by = 1)),expand=c(0.04, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_point(data=sql, size=2, color="black", shape=5)+
  labs(tag = "a", x = "Chromosomal location",
	   y = expression(bold(-log[10](CCA~bolditalic(P)~value)))) +
  ggtitle("NURTuRE-CKD")+
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

jpeg("manh_plot_1a.jpeg", pointsize=6, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(manhplot_nurture, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()
