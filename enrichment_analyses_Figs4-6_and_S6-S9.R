# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.

library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(grid)
library(DOSE)
library(venn)
library(ggpolypath)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("ReactomePA")
library(ReactomePA)
BiocManager::install("msigdbr")
library(msigdbr)

# ---- Get gene set results (list of genes in one column with header 'gene', and one gene symbol per row (cell), saved as .csv format):
nurture_CKD_CCA <- read.csv("nurture_277sig_genes.csv", header=T)
nurture_CKD_CCA_genes <- nurture_CKD_CCA$gene
entrez_ids <- mapIds(org.Hs.eg.db, nurture_CKD_CCA_genes, 'ENTREZID', 'SYMBOL')

ckdgen_CCA_all_159 <- read.table("ckdgen_159sig_genes.txt", header=T)
ckdgen_CCA_all_159_genes <- ckdgen_CCA_all_159$HGNC_Approved_sym
ckdgen_entrez_ids <- mapIds(org.Hs.eg.db, ckdgen_CCA_all_159_genes, 'ENTREZID', 'SYMBOL')

ukhls_CCA_all_246 <- read.table("ukhls_246sig_genes.txt", header=T)
ukhls_CCA_all_246_genes <- ukhls_CCA_all_246$HGNC_Approved_sym
ukhls_entrez_ids <- mapIds(org.Hs.eg.db, ukhls_CCA_all_246_genes, 'ENTREZID', 'SYMBOL')

bbj_CCA_all_181 <- read.table("bbj_181sig_genes.txt", header=T)
bbj_CCA_all_181_genes <- bbj_CCA_all_181$Approved_symbol
bbj_entrez_ids <- mapIds(org.Hs.eg.db, bbj_CCA_all_181_genes, 'ENTREZID', 'SYMBOL')

# === Gene Ontology ===
# NURTuRE-CKD gene set - Gene Ontology enrichment - Table S10 and Fig. S6.
nurture_all3 <- enrichGO(entrez_ids, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
dotplot_all3_35_nurture <-dotplot(nurture_all3, showCategory=35, font.size = 8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free",space="free_y")
vp1 <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("dotplot_all3_35_nurture.pdf", width=8, height=11)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(dotplot_all3_35_nurture, vp=vp1)
dev.off()
# and save results in table format:
write.csv(nurture_all3,"dotplot_all3_35_nurture_table.csv",quote=F,row.names=F,col.names=T)

# CKDGen gene set - Gene Ontology enrichment - Table S10 and Fig. S7.
ckdgen_all3 <- enrichGO(ckdgen_entrez_ids, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
dotplot_all3_30_ckdgen <-dotplot(ckdgen_all3, showCategory=32, font.size = 8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free",space="free_y")
vp1 <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("dotplot_all3_30_ckdgen.pdf", width=8, height=11)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(dotplot_all3_30_ckdgen, vp=vp1)
dev.off()
# and save results in table format:
write.csv(ckdgen_all3,"dotplot_all3_30_ckdgen_table.csv",quote=F,row.names=F,col.names=T)

# UKHLS gene set - Gene Ontology enrichment - Table S10 and Fig. S8.
ukhls_all3 <- enrichGO(ukhls_entrez_ids, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
dotplot_all3_30_ukhls <-dotplot(ukhls_all3, showCategory=40, font.size = 8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free",space="free_y")
vp1 <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("dotplot_all3_30_ukhls.pdf", width=8, height=11)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(dotplot_all3_30_ukhls, vp=vp1)
dev.off()
# and save results in table format:
write.csv(ukhls_all3,"dotplot_all3_30_ukhls_table.csv",quote=F,row.names=F,col.names=T)

# BioBank Japan gene set - Gene Ontology enrichment - Table S10 and Fig. S9.
bbj_all3 <- enrichGO(bbj_entrez_ids, OrgDb = "org.Hs.eg.db", ont="all", readable=TRUE)
dotplot_all3_30_bbj <-dotplot(bbj_all3, showCategory=21, font.size = 8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free",space="free_y")
vp1 <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("dotplot_all3_30_bbj.pdf", width=8, height=11)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(dotplot_all3_30_bbj, vp=vp1)
dev.off()
write.csv(bbj_all3,"dotplot_all3_30_bbj_table.csv",quote=F,row.names=F,col.names=T)

# === Pathway databases ===
# --- NURTuRE-CKD gene set - KEGG, Reactome, WikiPathway enrichment - Table S11 and Fig. 5a
kk_nurture <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
# no significant terms for KEGG
WP_nurture <- enrichWP(entrez_ids, organism = "Homo sapiens")
# no significant terms for WikiPathway
react_nurture <- enrichPathway(gene=entrez_ids, pvalueCutoff = 0.05, readable=TRUE)

# --- CKDGen gene set - KEGG, Reactome, WikiPathway enrichment - Table S11 and Fig. 5b
kk_ckdgen <- enrichKEGG(gene = ckdgen_entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
react_ckdgen <- enrichPathway(gene=ckdgen_entrez_ids, pvalueCutoff = 0.05, readable=TRUE)
WP_ckdgen <- enrichWP(ckdgen_entrez_ids, organism = "Homo sapiens")

eval_func = function(x){
  generatio = x[3]
  eval(parse(text = generatio))
}

WP_ckdgen_toplot <- WP_ckdgen[]
WP_ckdgen_toplot$GeneRatio2 <- apply(WP_ckdgen_toplot, 1, eval_func)
WP_ckdgen_toplot$db <- "WikiPathway"
react_ckdgen_toplot <- react_ckdgen[]
react_ckdgen_toplot$GeneRatio2 <- apply(react_ckdgen_toplot, 1, eval_func)
react_ckdgen_toplot$db <- "Reactome"
kk_ckdgen_toplot <- kk_ckdgen[]
kk_ckdgen_toplot$GeneRatio2 <- apply(kk_ckdgen_toplot, 1, eval_func)
kk_ckdgen_toplot$db <- "KEGG"

ckdgen_reac_WP_KEGG <- rbind(WP_ckdgen_toplot,react_ckdgen_toplot,kk_ckdgen_toplot)
# display long terms over multiple lines by adding \n as line breaks:
ckdgen_reac_WP_KEGG[ckdgen_reac_WP_KEGG$Description=="Transport of bile salts and organic acids,\nmetal ions and amine compounds",2] <- "Transport of bile salts and\norganic acids, metal ions\nand amine compounds"
ckdgen_reac_WP_KEGG[ckdgen_reac_WP_KEGG$Description=="FOXO-mediated transcription of oxidative\nstress, metabolic and neuronal genes",2] <- "FOXO-mediated transcription\nof oxidative stress,\nmetabolic and neuronal genes"
ckdgen_reac_WP_KEGG[ckdgen_reac_WP_KEGG$Description=="Hippo signaling regulation pathways",2] <- "Hippo signaling\nregulation pathways"

dotplot_all3_pways_ckdgen <- ggplot(ckdgen_reac_WP_KEGG, aes(x=GeneRatio2, y=reorder(Description, GeneRatio2), color=p.adjust, size=Count)) +
  geom_point() +
  facet_wrap(~ db) +
  labs(tag = "b", x = "GeneRatio", y = "Term description") +
  scale_color_continuous(low="red", high="blue", name = "p.adjust",
            guide=guide_colorbar(reverse=TRUE)) +
  scale_size(range=c(3, 8),breaks=c(4,6,8))+
  scale_x_continuous(expand=c(0.006, 0.006),limits = c(0.02,0.12),breaks=c(0.02,0.06,0.1)) +
  theme(axis.text.x = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size=10,colour = "black"),
    axis.text.y = element_text(size=9,colour = "black"),
    axis.title.y = element_text(size=10,colour = "black"),
    legend.title = element_text(size=10,colour = "black"),
    panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    plot.tag = element_text(face="bold", size = 17,color="black")
  )

# --- BioBank Japan gene set - KEGG, Reactome, WikiPathway enrichment - Table S11 and Fig. 5a
kk_bbj <- enrichKEGG(gene = bbj_entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
# 1 found.
react_bbj <- enrichPathway(gene=bbj_entrez_ids, pvalueCutoff = 0.05, readable=TRUE)
# no significant terms for Reactome.
WP_bbj <- enrichWP(bbj_entrez_ids, organism = "Homo sapiens")
# no significant terms for WikiPathways.

# --- UKHLS gene set - KEGG, Reactome, WikiPathway enrichment:
kk_ukhls <- enrichKEGG(gene = ukhls_entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
# no significant terms for KEGG.
react_ukhls <- enrichPathway(gene=ukhls_entrez_ids, pvalueCutoff = 0.05, readable=TRUE)
# no significant terms for Reactome.
WP_ukhls <- enrichWP(ukhls_entrez_ids, organism = "Homo sapiens")
# no significant terms for WikiPathway.

# --- Create plot with all sig. enriched pathway results (as facet_wrap for each gene set):
react_nurture$dataset <- "NURTuRE-CKD, Reactome"
kk_bbj$dataset <- "BioBank Japan, KEGG"

kk_bbj_toplot <- kk_bbj[]
kk_bbj_toplot$GeneRatio2 <- apply(kk_bbj_toplot, 1, eval_func)
kk_bbj_toplot$dataset <- "BioBank Japan\nKEGG"
react_nurture_toplot <- react_nurture[]
react_nurture_toplot$GeneRatio2 <- apply(react_nurture_toplot, 1, eval_func)
react_nurture_toplot$dataset <- "NURTuRE-CKD\nReactome"
bbj_kegg_and_nurt_reactome <- rbind(kk_bbj_toplot,react_nurture_toplot)

dotplot_bbj_kegg_and_nurt_reactome <- ggplot(bbj_kegg_and_nurt_reactome, aes(x=GeneRatio2, y=reorder(Description, p.adjust), color=p.adjust, size=Count)) +
  geom_point() +
  facet_wrap(~ dataset) +
  labs(tag='a', x = "GeneRatio", y = "Term description") +
  scale_color_continuous(low="red", high="blue", name = "p.adjust",
            guide=guide_colorbar(reverse=TRUE)) +
  scale_size(range=c(3, 6),breaks=c(4,10,16))+
  scale_x_continuous(expand=c(0.006, 0.006),limits = c(0.02,0.12),breaks=c(0.02,0.06,0.1)) +
  theme(axis.text.x = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size=10,colour = "black"),
    axis.text.y = element_text(size=9,colour = "black"),
    axis.title.y = element_text(size=10,colour = "black"),
    legend.title = element_text(size=8,colour = "black"),
    panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    plot.tag = element_text(face="bold", size = 17,color="black")
  )

pdf("dotplot_all3_pways_nurt_bbj_ckdgen.pdf", width=7, height=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1, heights = c(0.85, 1))))
print(dotplot_bbj_kegg_and_nurt_reactome, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(dotplot_all3_pways_ckdgen, vp=viewport(layout.pos.row=2,layout.pos.col=1))
dev.off()
# save results in table format
write.csv(ckdgen_reac_WP_KEGG,"ckdgen_pathways_table.csv",quote=F,row.names=F,col.names=T)
write.csv(bbj_kegg_and_nurt_reactome,"bbj_nurtureckd_pathways_table.csv",quote=F,row.names=F,col.names=T)


# ===== MSig database - C8 (cell-types) =====
msig_C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, entrez_gene)

# --- NURTuRE-CKD gene set - Fig. 6
msig_C8_nurture <- enricher(entrez_ids, TERM2GENE=msig_C8)

# --- CKDGen gene set - Fig. 6
msig_C8_ckdgen <- enricher(ckdgen_entrez_ids, TERM2GENE=msig_C8)

# --- UKHLS gene set - Fig. 6
msig_C8_ukhls <- enricher(ukhls_entrez_ids, TERM2GENE=msig_C8)

# --- BBJ - Fig. 6
msig_C8_bbj <- enricher(bbj_entrez_ids, TERM2GENE=msig_C8)

# Create dotplot for MSigDB C8 results for all 4 gene sets:
msig_C8_ukhls@result$dataset <- "ukhls"
msig_C8_ckdgen@result$dataset <- "ckdgen"
msig_C8_nurture@result$dataset <- "nurture"
msig_C8_bbj@result$dataset <- "bbj"
msig_C8_all <- rbind(msig_C8_ukhls@result[c(1:10),], msig_C8_ckdgen@result[c(1:10),], msig_C8_nurture@result[c(1:10),], msig_C8_bbj@result[c(1:10),])

msig_C8_all$GeneRatio2 <- apply(msig_C8_all, 1, eval_func_enrichr)
new_labels <- c("ckdgen" = "CKDGen metaCCA", "bbj" = "BioBank Japan metaCCA", "ukhls" = "UKHLS metaCCA", "nurture" = "NURTuRE-CKD CCA")

dotplot_msigdb_C8_all <- ggplot(msig_C8_all, aes(x=GeneRatio2, y=reorder(Description, Count), color=p.adjust, size=Count)) +
  geom_point() +
  facet_wrap(~ dataset, ncol = 4,labeller = labeller(dataset = new_labels)) +
  labs(x = "GeneRatio", y = "Term description") +
  scale_color_continuous(low="red", high="blue", name = "p.adjust",
            guide=guide_colorbar(reverse=TRUE)) +
  scale_size(range=c(2, 5))+
  scale_x_continuous(expand=c(0.006, 0.006)) +
  theme(axis.text.x = element_text(size=6,colour = "black"),
    axis.title.x = element_text(size=8,colour = "black"),
    axis.text.y = element_text(size=6,colour = "black"),
    axis.title.y = element_text(size=8,colour = "black"),
    legend.title = element_text(size=10,colour = "black"),
    panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "lightgrey"),
    strip.text.x = element_text(size = 6)
  )
pdf("dotplot_msigdb_C8_all.pdf", width=9, height=4)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(dotplot_msigdb_C8_all, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()
write.csv(msig_C8_all,"msig_C8_all_table.csv",quote=F,row.names=F,col.names=T)


# ----- Look for any significantly enriched Gene Ontology terms in common between 4 gene sets:
# Fig. 4a
nurture_all3_terms <- nurture_all3$Description
ckdgen_all3_terms <- ckdgen_all3$Description
ukhls_all3_terms <- ukhls_all3$Description
bbj_all3_terms <- bbj_all3$Description

overlap_4 <- list(
  "CKDGen\nmetaCCA\n(64)" = ckdgen_all3_terms,
  "UKHLS\nmetaCCA (40)" = ukhls_all3_terms,
  "BBJ\nmetaCCA (21)" = bbj_all3_terms,
  "NURTuRE_CKD (138)" = nurture_all3_terms
  )

venn_GO <- venn(overlap_4,sncs=1, ellipse = TRUE, counts= TRUE, ilabels= TRUE, ggplot = TRUE, ilcs =18)
ckdgen_nurture_GOoverlap <- intersect(nurture_all3_terms,ckdgen_all3_terms)
ukhls_nurture <- intersect(ukhls_all3_terms,nurture_all3_terms)

write.table(ukhls_nurture[order(ukhls_nurture)],"venn_4datasets_intersect_ukhls_nurture.txt",quote=F,row.names=F)

vp <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("venn_4datasets_sigGOterms.pdf", width=8, height=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(venn_GO, vp=vp)
dev.off()

# ---- Use upset and cnet plots on nurture-CKD gene set GO enrichment results to visualise network genes that overlap between:
# kidney terms (also in CKDGen gene set) and neuro-based terms (also in UKHLS gene set):
# Fig. 4b
upsetplot_nurture_allGO <-upsetplot(nurture_all3,n=30)
terms_to_show <- c("synapse organization","neuron projection guidance","metal ion transmembrane transporter activity","neuronal cell body","renal system development")
cnetplot <- enrichplot::cnetplot(nurture_all3, showCategory = terms_to_show, foldChange=nurture_all3_genes)
vp <- viewport(width = 1, height = 1, layout.pos.row=1,layout.pos.col=1)
pdf("cnetplot.pdf", width=8, height=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
print(cnetplot, vp=vp)
dev.off()
