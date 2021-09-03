# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.

# This script uses the published significant kidney cell expression gene sets by Gillies et al., 2018
# to look for enrichment of kidney cell expressed genes
# Gillies et al., 2018:
# Gillies CE, Putler R, Menon R, Otto E, Yasutake K, Nair V, Hoover P, Lieb D, Li S, Eddy S, Fermin D,
# McNulty MT; Nephrotic Syndrome Study Network (NEPTUNE), Hacohen N, Kiryluk K, Kretzler M, Wen X, Sampson MG.
# An eQTL Landscape of Kidney Tissue in Human Nephrotic Syndrome. Am J Hum Genet. 2018 Aug 2;103(2):232-244.
# doi: 10.1016/j.ajhg.2018.07.004. Epub 2018 Jul 26. PMID: 30057032; PMCID: PMC6081280.
# Table S3. Genes that Were Expressed in the Single-Cell RNA-Seq Experiment and the Cell Cluster in which Its Expresssion Was Enriched
# each tab then saved as .csv file for input.

library(hypeR)

# Get all 4 gene sets:
nurture_CKD_CCA <- read.csv("nurture_CKD_CCA_all_277sig_HGNCgenes_0.01.csv", header=T)
nurture_CKD_CCA_genes <- nurture_CKD_CCA$gene

ckdgen_CCA_all_159 <- read.table("ckdgen_159sig_genes.txt", header=T)
ckdgen_CCA_all_159_genes <- ckdgen_CCA_all_159$HGNC_Approved_sym

ukhls_CCA_all_246 <- read.table("ukhls_246sig_genes.txt", header=T)
ukhls_CCA_all_246_genes <- ukhls_CCA_all_246$HGNC_Approved_sym

bbj_CCA_all_181 <- read.table("bbj_181sig_genes.txt", header=T)
bbj_CCA_all_181_genes <- bbj_CCA_all_181$Approved_symbol

# upload Table S3 form Gillies et al., 2018 as two separate .csv files:
all_files <- list.files(path = "pub_expr_datasets", pattern = "^Gillies2018_ST3_", all.files = FALSE,
           full.names = TRUE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

table <- data.frame(gene = character(), p_val=numeric(), avg_logFC=numeric(), pct.1=numeric(), pct.2=numeric(), p_val_adj = numeric(), stringsAsFactors=FALSE)

for (i in 1:length(all_files)){
  if(i==1){
    table <- read.csv(all_files[i])
    table$dataset <- gsub("Gillies2018_ST3_","",all_files[i])
  }
  else{
    print(all_files[i])
    new_csv <- read.csv(all_files[i])
    new_csv$dataset <- gsub("Gillies2018_ST3_","",all_files[i])
    table <- rbind(table,new_csv)
  }
}

results_table_nurt <- data.frame(gene = character(), dataset_sig=character(),stringsAsFactors=FALSE)
for(i in 1:length(nurture_CKD_CCA_genes)){
  results_table_nurt[i,1] <- nurture_CKD_CCA_genes[i]
  sig_datasets <- unique(table[table$Gene == nurture_CKD_CCA_genes[i],2])
  results_table_nurt[i,2] <- paste(gsub(".csv","",sig_datasets[!is.na(sig_datasets)]), collapse="; ")
}

results_table_ukhls <- data.frame(gene = character(), dataset_sig=character(),stringsAsFactors=FALSE)
for(i in 1:length(ukhls_CCA_all_246_genes)){
  results_table_ukhls[i,1] <- ukhls_CCA_all_246_genes[i]
  sig_datasets <- unique(table[table$Gene == ukhls_CCA_all_246_genes[i],2])
  results_table_ukhls[i,2] <- paste(gsub(".csv","",sig_datasets[!is.na(sig_datasets)]), collapse="; ")
}

results_table_ckdgen <- data.frame(gene = character(), dataset_sig=character(),stringsAsFactors=FALSE)
for(i in 1:length(ckdgen_CCA_all_159_genes)){
  results_table_ckdgen[i,1] <- ckdgen_CCA_all_159_genes[i]
  sig_datasets <- unique(table[table$Gene == ckdgen_CCA_all_159_genes[i],2])
  results_table_ckdgen[i,2] <- paste(gsub(".csv","",sig_datasets[!is.na(sig_datasets)]), collapse="; ")
}

results_table_bbj <- data.frame(gene = character(), dataset_sig=character(),stringsAsFactors=FALSE)
for(i in 1:length(bbj_CCA_all_181_genes)){
  results_table_bbj[i,1] <- bbj_CCA_all_181_genes[i]
  sig_datasets <- unique(table[table$Gene == bbj_CCA_all_181_genes[i],2])
  results_table_bbj[i,2] <- paste(gsub(".csv","",sig_datasets[!is.na(sig_datasets)]), collapse="; ")
}
results_table_bbj$dataset <- "bbj"
results_table_ckdgen$dataset <- "ckdgen"
results_table_ukhls$dataset <- "ukhls"
results_table_nurt$dataset <- "nurture"
results_gillies_all4 <- rbind(results_table_bbj,results_table_ckdgen,results_table_ukhls,results_table_nurt)

list_clusters <- unique(table$cluster)

genesets <- list("T cells" = sort(table[table$cluster==list_clusters[1],1]),
                  "Proximal" = sort(table[table$cluster==list_clusters[2],1]),
                  "Endothelial1" = sort(table[table$cluster==list_clusters[3],1]),
                  "Endothelial2" = sort(table[table$cluster==list_clusters[4],1]),
                  "Mast /Immune" = sort(table[table$cluster==list_clusters[5],1]),
                  "Endothelial3" = sort(table[table$cluster==list_clusters[6],1]),
                  "Mesangial" = sort(table[table$cluster==list_clusters[7],1]),
                  "dLOH" = sort(table[table$cluster==list_clusters[8],1]),
                    "aLOH" = sort(table[table$cluster==list_clusters[9],1]),
                    "Podocyte" = sort(table[table$cluster==list_clusters[10],1]),
                    "Intercalated/CNT" = sort(table[table$cluster==list_clusters[11],1]),
                    "Collecting duct/DCT/CNT" = sort(table[table$cluster==list_clusters[12],1]),
                    "CD8+/NK" = sort(table[table$cluster==list_clusters[13],1]),
                    "Myeloid" = sort(table[table$cluster==list_clusters[14],1])
                  )

# background genes: Transcriptome Affymetrix 2.1 Array, 22,893 genes, Fig S1: https://www.cell.com/cms/10.1016/j.ajhg.2018.07.004/attachment/83682a69-f9b6-447c-9a5b-f9a9581ad0a6/mmc1
hyp_obj_nurt <- hypeR(nurture_CKD_CCA_genes, genesets, background=22893)
hyp_obj_ckdgen <- hypeR(ckdgen_CCA_all_159_genes, genesets, background=22893)
hyp_obj_bbj <- hypeR(bbj_CCA_all_181_genes, genesets, background=22893)
hyp_obj_ukhls <- hypeR(ukhls_CCA_all_246_genes, genesets, background=22893)

# create dot plot to show enrichment results:
dotplot_gillies_all4_enrich <- ggplot(hyp_obj_all4_long[hyp_obj_all4_long$fdr<0.05,], aes(x=fdr, y=reorder(label, fdr), color=fdr, size=overlap)) +
    geom_point() +
    facet_wrap(~ dataset, ncol = 4,labeller = labeller(dataset = new_labels)) +
    labs(x = "FDR", y = "Significantly enriched kidney cell type (FDR < 0.05)") +
    scale_color_continuous(low="red", high="blue", name = "FDR",
              guide=guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(breaks=c(10,15,20), range = c(2, 5), name = "Gene overlap") +
    scale_x_continuous(expand=c(0.006, 0.006)) +
    theme(axis.text.x = element_text(size=8,colour = "black"),
      axis.title.x = element_text(size=10,colour = "black",face = "bold"),
      axis.text.y = element_text(size=8,colour = "black"),
      axis.title.y = element_text(size=8,colour = "black",face = "bold"),
      legend.title = element_text(size=10,colour = "black",face = "bold"),
      panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
      strip.text.x = element_text(size = 9,face = "bold")
    )


vp <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("dotplot_gillies_all4_enrich.pdf", width=7.5, height=3)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(dotplot_gillies_all4_enrich, vp=vp)
dev.off()

hyp_obj_all4$hits <- gsub(",","/",hyp_obj_all4$hits)
write.csv(hyp_obj_all4,"dotplot_gillies_all4_enrich_table.csv",quote=F,row.names=F,col.names=T)
