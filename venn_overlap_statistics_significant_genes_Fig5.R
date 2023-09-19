# This R script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
# Venn diagrams and overlap statistics for significant gene sets.

library(ggvenn)
library(qdap)
library(venn)

# For Fig. 3
CKDGen_GWAS <- read.table("wuttke_122snps.txt", header=F)
CKDGen_GWAS_genes <- CKDGen_GWAS$V4
CKDGen_metaCCA <- read.table("ckdgen_159sig_genes.txt", header=T)
CKDGen_metaCCA_genes <- multigsub(c("-"), c("_"), CKDGen_metaCCA$HGNC_Approved_sym)
UKHLS_metaCCA <- read.table("ukhls_246sig_genes.txt", header=T)
UKHLS_metaCCA_genes <- multigsub(c("-"), c("_"), UKHLS_metaCCA$HGNC_Approved_sym)
BBJ_metaCCA <- read.table("bbj_180sig_genes.txt", header=T)
BBJ_metaCCA_genes <- multigsub(c("-"), c("_"),BBJ_metaCCA$Approved_symbol)
BBJ_GWAS <- read.table("BBJ_GWAS_results_eGFR_BUN_7genes.csv", header=T)
BBJ_GWAS_genes <- BBJ_GWAS$BBJ_genes_eGFR_and_BUN

# Can replace above files for SNP results and available genes or SNPs (for Figs. S1 and S2)

x_3 <- list(
  "CKDGen\nmetaCCA\n(159)" = CKDGen_metaCCA_genes,
  "UKHLS metaCCA (246)" = UKHLS_metaCCA_genes,
  "BiobankJapan\nmetaCCA\n(180)" = BBJ_metaCCA_genes
  )

x_4 <- list(
  "CKDGen\nmetaCCA\n(159)" = CKDGen_metaCCA_genes,
  "UKHLS\nmetaCCA (246)" = UKHLS_metaCCA_genes,
  "BiobankJapan\nmetaCCA (180)" = BBJ_metaCCA_genes,
  "CKDGen\nGWAS (122)" = CKDGen_GWAS_genes
  )

x_4_BBJ <- list(
  "CKDGen\nmetaCCA\n(159)" = CKDGen_metaCCA_genes,
  "UKHLS\nmetaCCA (246)" = UKHLS_metaCCA_genes,
  "BiobankJapan\nmetaCCA (180)" = BBJ_metaCCA_genes,
  "BBJ\nGWAS (7)" = BBJ_GWAS_genes
  )

x_5 <- list(
  "CKDGen\nmetaCCA (159)" = CKDGen_metaCCA_genes,
  "UKHLS metaCCA (246)" = UKHLS_metaCCA_genes,
  "BiobankJapan\nmetaCCA (180)" = BBJ_metaCCA_genes,
  "CKDGen GWAS (122)" = CKDGen_GWAS_genes,
  "BBJ\nGWAS (7)" = BBJ_GWAS_genes
  )

CKDGen_GWAS_genes[CKDGen_GWAS_genes %in% BBJ_GWAS_genes]
x_2GWAS <- list(
  "CKDGen GWAS (122)" = CKDGen_GWAS_genes,
  "BBJ GWAS (7)" = BBJ_GWAS_genes
  )

x_2_BBJ <- list(
  "BiobankJapan\nmetaCCA (180)" = BBJ_metaCCA_genes,
  "BBJ\nGWAS (7)" = BBJ_GWAS_genes
  )

x_2_CKDGen <- list(
  "CKDGen\nmetaCCA (159)" = CKDGen_metaCCA_genes,
  "CKDGen\nGWAS (122)" = CKDGen_GWAS_genes
  )

venn_GWAS <- venn(x_2GWAS, sncs=1, ilcs =1, ggplot = TRUE)
venn3 <- venn(x_3,sncs=1, ilcs =1, ggplot = TRUE)
venn4 <- venn(x_4,sncs=1, ilcs =1,ellipse = TRUE, ggplot = TRUE)
venn5 <- venn(x_5,sncs=1, ilcs =1,ellipse = TRUE, ilabels = TRUE, ggplot = TRUE)

venn4_BBJ <- venn(x_4_BBJ,sncs=1, ilcs =1,ellipse = TRUE, ggplot = TRUE)

venn2_CKDGen <- venn(x_2_CKDGen, sncs=1, ilcs =1, ggplot = TRUE)
venn2_BBJ <- venn(x_2_BBJ,sncs=1, ilcs =1,ggplot = TRUE)

# Include NURTuRE-CKD gene set:
nurture_CKD_CCA <- read.csv("nurture_CKD_CCA_sig_genes.csv", header=T)
nurture_CKD_CCA_genes <- nurture_CKD_CCA$gene
nurture_CKD_CCA_genes <- multigsub(c("-"), c("_"), nurture_CKD_CCA_genes)

x_4_new_4 <- list(
  "CKDGen\nmetaCCA\n(159)" = CKDGen_metaCCA_genes,
  "UKHLS\nmetaCCA (246)" = UKHLS_metaCCA_genes,
  "BBJ\nmetaCCA (181)" = BBJ_metaCCA_genes,
  "NURTuRE-CKD\nCCA (277)" = nurture_CKD_CCA_genes
  )
venn <- venn(x_4_new_4,sncs=1, ellipse = TRUE, counts= TRUE, ilabels= TRUE, ggplot = TRUE, ilcs =12)

vp1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.75,layout.pos.row=1,layout.pos.col=1)
vp2 <- viewport(width = 1, height = 1, x = 0.66, y = 0.75,layout.pos.row=1,layout.pos.col=2)
vp3 <- viewport(width = 1, height = 1, x = 0.33, y = 0.5,layout.pos.row=2,layout.pos.col=1)
vp4 <- viewport(width = 1, height = 1, x = 0.66, y = 0.5,layout.pos.row=2,layout.pos.col=2)

vp5 <- viewport(width = 1, height = 1, x = 0.33, y = 0.25,layout.pos.row=1,layout.pos.col=1)
vp6 <- viewport(width = 1, height = 1, x = 0.66, y = 0.25,layout.pos.row=1,layout.pos.col=2)

# Used in Figure 3:
pdf("venn_sig_genes_Fig3.pdf", width=8, height=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(venn)
dev.off()


pdf("venn_sig_genes_4.pdf", width=8, height=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
print(venn_GWAS, vp=vp1)
print(venn3, vp=vp2)
print(venn2_CKDGen, vp=vp3)
print(venn2_BBJ, vp=vp4)
dev.off()

overlap_ckdgen_nurt <- intersect(CKDGen_metaCCA_genes,nurture_CKD_CCA_genes)
overlap_japan_nurt <- intersect(BBJ_metaCCA_genes,nurture_CKD_CCA_genes)
overlap_japan_ckdgen <- intersect(BBJ_metaCCA_genes,CKDGen_metaCCA_genes)
overlap_japan_ukhls <- intersect(BBJ_metaCCA_genes,UKHLS_metaCCA_genes)
overlap_nurt_ukhls <- intersect(nurture_CKD_CCA_genes,UKHLS_metaCCA_genes)
overlap_ckdgen_ukhls <- intersect(CKDGen_metaCCA_genes,UKHLS_metaCCA_genes)

# ====== Overlap statistics (using genes in common only):
library(SuperExactTest)

CKDGen_metaCCA <- read.table("common5292_ckdgen_100genes.txt", header=F)
CKDGen_metaCCA_genes <- multigsub(c("-"), c("_"), CKDGen_metaCCA$V1)
UKHLS_metaCCA <- read.table("common5292_ukhls_192genes.txt", header=F)
UKHLS_metaCCA_genes <- multigsub(c("-"), c("_"), UKHLS_metaCCA$V1)
BBJ_metaCCA <- read.table("common5292_bbj_122genes.txt", header=F)
BBJ_metaCCA_genes <- multigsub(c("-"), c("_"),BBJ_metaCCA$V1)

# list of common genes available for analysis between these 3 datasets:
all_5292genes <- read.csv("5292commongenes_HGNCchecked.csv",header=T)
all_5292genes_ready <- unique(multigsub(c("-"), c("_"), all_5292genes$Approved_symbol))

# list of genes sig. in NURTuRE-CKD that were available for analysis in these 3 datasets:
nurture_ckd_all_in_5292common <- nurture_CKD_CCA_genes[nurture_CKD_CCA_genes %in% all_5292genes_ready]

# list of genes available for analysis in NURTuRE-CKD:
nurture_all_avail <- read.csv("nurture-ckd_CCA_all_avail_genes.csv", header=T)
nurture_all_avail_genes <- unique(nurture_all_avail$gene)
# list of common genes available for analysis between all 4 datasets (now including NURTuRE-CKD):
nurture_ckd_all_genes_inother3 <- nurture_all_avail_genes[nurture_all_avail_genes %in% all_5292genes_ready]

# Use list of common genes between all 4 datasets: "nurture_ckd_all_genes_inother3" to get sig. common genes for each:
nurture_CKD_CCA_genesinAll4 <- nurture_CKD_CCA_genes[nurture_CKD_CCA_genes %in% nurture_ckd_all_genes_inother3]
CKDGen_metaCCA_genesinAll4 <- CKDGen_metaCCA_genes[CKDGen_metaCCA_genes %in% nurture_ckd_all_genes_inother3]
UKHLS_metaCCA_genesinAll4 <- UKHLS_metaCCA_genes[UKHLS_metaCCA_genes %in% nurture_ckd_all_genes_inother3]
BBJ_metaCCA_genesinAll4 <- BBJ_metaCCA_genes[BBJ_metaCCA_genes %in% nurture_ckd_all_genes_inother3]

# Compute overlap statistics:
list <- list('CKDGen' = CKDGen_metaCCA_genesinAll4, 'UKHLS' = UKHLS_metaCCA_genesinAll4, 'BioBankJapan' = BBJ_metaCCA_genesinAll4, 'NURTuRE_CKD' = nurture_CKD_CCA_genesinAll4)
Result = supertest(list, n=length(nurture_ckd_all_genes_inother3))
plot(Result, degree=2:3, sort.by='size')
Results_table <- summary(Result)$Table
Results_table$pBonf <- p.adjust(Results_table$P.value,method="bonferroni")
Results_table$Elements <- gsub(", ","_",Results_table$Elements)
Results_table$Intersections <- gsub(" ","_",Results_table$Intersections)
Results_table[order(Results_table$pBonf),c(1,3,8)]
write.table(Results_table,"overlap_4_genesets_statistics.txt",col.names=T,row.names=F,quote=F)
