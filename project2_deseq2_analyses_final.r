#conda activate deseq2

#R
library(DESeq2)
library(EnhancedVolcano)
library(data.table)
library(curl)
library(R.utils)
library(dplyr)
library(biomaRt)
library(gtools)

'%notin%' <- Negate('%in%')

visitdates_RNA <- read.csv("RNA_visit_dates.csv", header=T)

rna_IDs <- read.csv("RNA_sample.csv",header=F)
names(rna_IDs) <- c("PPID","RNA_ID")
rna_IDs$RNA_ID <- gsub("-","_",rna_IDs$RNA_ID)

# === Add in controls:
controls <- read.csv("Control_GSA_samples.csv", header=T)
rna_IDs[rna_IDs$PPID %in% controls$PPID,]
rna_IDs.controls <- merge(rna_IDs, controls, by = c("PPID"))

# Next: do cases vs ctrls, inc case subtypes:
rnaseq <- read.csv("rnaseq_raw.csv")
rnaseq.1st <- read.csv("counts.tab.csv",header=T)
names(rnaseq.1st)[1] <- "ensembl_id"
samplestab <- read.csv("samples.tab.csv",header=T)

genes <- read.csv("RNAseq-genes_HGNC_ensembl.csv")
ppids <- read.csv("RNAseq-patientIDs.csv")

visitdates_RNA <- visitdates_RNA[,c(1,6:8)]
names(visitdates_RNA) <- c("PPID","SampleID","Visit_Clinical.Status","Visit_Visit.Date")
visitdates_RNA$SampleID <- paste0("RNA_",visitdates_RNA$SampleID)

# === Add in controls:
rna_IDs.controls <- rna_IDs.controls[,c(2,1)]
names(rna_IDs.controls) <- c("RNA_ID","PPID")
# get baseline / visit data for controls:
rna_IDs.controls$visit <- "baseline"
rna_IDs.controls$PT <- "controls"

gene_id <- rnaseq.1st[,c(1),drop=F]
rnaseq.1st <- cbind(gene_id, rnaseq.1st[,colnames(rnaseq.1st) %in% samplestab$RNA_ID, drop=F])

control.samples <- read.csv("~/samples.control.tab.csv")
control.data <- read.csv("~/counts.control.tab.csv")
gene_id <- control.data[,c(1),drop=F]
names(gene_id) <- c("ensembl_id")
control.data <- cbind(gene_id, control.data[,colnames(control.data) %in% rna_IDs.controls$RNA_ID, drop=F])

control.samples.toBind <- control.samples[,c(1,5,6)]
names(control.samples.toBind) <- c("RNA_ID","PPID","visit")
control.samples.toBind <- control.samples.toBind[control.samples.toBind$RNA_ID %in% rna_IDs.controls$RNA_ID,]
control.samples.toBind$PT <- "control"

# merge rnaseq.1st and control.data_EUR
cases_and_controls <- merge(rnaseq.1st, control.data, by = c("ensembl_id"))

# get order of participants in data matrix:
all_RNAIDs <- colnames(cases_and_controls[,c(2:ncol(cases_and_controls))])
all_RNAIDs <- as.data.frame(all_RNAIDs)
all_RNAIDs$order <- seq.int(nrow(all_RNAIDs))
names(all_RNAIDs) <- c("RNA_ID","order")

# rbind samplestab.Sam and  control.samples.toBind, and order by order in data matrix
PT_cases_ctrls <- rbind(samplestab, control.samples.toBind)
PT_cases_ctrls <- merge(PT_cases_ctrls, all_RNAIDs, by = c("RNA_ID"))
PT_cases_ctrls <- PT_cases_ctrls[order(PT_cases_ctrls$order),]

# how many of each PT?
table(PT_cases_ctrls$PT)

PT_cases_ctrls[duplicated(PT_cases_ctrls$PPID),]
# 0 duplicated PPIDs (i.e. no potential technical replicates present)

# prepare rna matrix and replace ensemblG with gene symbols:
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl_gene_ids <- cases_and_controls$ensembl_id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=ensembl_gene_ids,mart= mart)
names(G_list) <- c("ensembl_id","hgnc_symbol")
# 53,349
G_list <- G_list[G_list$hgnc_symbol %notin% "",]
# 36,915
nrow(G_list[duplicated(G_list$hgnc_symbol),])
# 7
G_list[duplicated(G_list$hgnc_symbol),2] <- paste0(G_list[duplicated(G_list$hgnc_symbol),2],"_v2")

cases_and_controls.genes <- merge(G_list, cases_and_controls, by = c("ensembl_id"))
# 36,915

rnaseq.v2 <- cases_and_controls.genes
rownames(rnaseq.v2) <- rnaseq.v2$hgnc_symbol
rnaseq.v2 <- rnaseq.v2[,-c(1:2)]
rnaseq.v2.int <- as.data.frame(lapply(rnaseq.v2,as.integer))
rownames(rnaseq.v2.int) <- rownames(rnaseq.v2)

PTs_toInc <- sort(unique(PT_cases_ctrls$PT))[sort(unique(PT_cases_ctrls$PT)) %notin% c("","")]
# Create matrix of dataset1 ,dataset2:
PT_analyses <- combinations(n=length(PTs_toInc), r=2, v=PTs_toInc, repeats.allowed = FALSE)

# for loop:
genes_select <- genes$HGNC_symbol
plot_list = list()
plot_list2 = list()
for(i in 1:nrow(PT_analyses)){
	dataset1 = PT_analyses[i,1]
	dataset2 = PT_analyses[i,2]
	print(paste0("On run: ", i, " for ", dataset1, " vs. ", dataset2))
	PTdata_torun <- PT_cases_ctrls[PT_cases_ctrls$PT %in% c(dataset1, dataset2),]
	RNAseqdata_torun <- rnaseq.v2.int[,colnames(rnaseq.v2.int) %in% PTdata_torun$RNA_ID]

	deseq2_matrix <- DESeqDataSetFromMatrix(countData = RNAseqdata_torun,
								  colData = PTdata_torun,
								  design = ~ PT)
							

	deseq2_matrix_out <- DESeq(deseq2_matrix)
	length1 <- length(unique(PTdata_torun[PTdata_torun$PT %in% dataset1,2]))
	length2 <- length(unique(PTdata_torun[PTdata_torun$PT %in% dataset2,2]))

	# with count filter:
	# filter out genes with counts < 10:
	smallestGroupSize <- min(length1, length2)
	keep <- rowSums(counts(deseq2_matrix_out) >= 10) >= smallestGroupSize
	deseq2_matrix_out_Filt <- deseq2_matrix_out[keep,]
	res <- results(deseq2_matrix_out_Filt)
	results <- as.data.frame(res)

	# check:
	head(results[order(results$padj),c(6,2)])
	gene_tolabel <- rownames(results %>% filter(padj<0.05 & abs(log2FoldChange) >= 1))
	# save all results:
	write.csv(results, 
		paste0("dseq2_cases_",length1,dataset1,"_vs_",length2,dataset2,"_Filt_baseline_Allgenes_results.csv"), quote=F, row.names=T)

	# Try >=2 logfoldchange for labels
	gene_tolabel2 <- rownames(results %>% filter(padj<0.05 & abs(log2FoldChange) >= 2))
	plot1 <- EnhancedVolcano(results, 
		title = paste0("case gene expression (baseline):\n",dataset1," (n=",length1,") vs. ",dataset2," (n=",length2,")\n(excl. genes with low counts)"), 
		lab=rownames(results), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, ylab = bquote(~-Log[10]~ 'adjusted P-value'), 
		selectLab = c(gene_tolabel), pointSize = 2.0, labSize = 3.0, labCol = 'black',drawConnectors = TRUE, widthConnectors = 0.4, colConnectors = 'black',
		max.overlaps = 30,  maxoverlapsConnectors = 30) #max number of overlaps allowed)

	plot_list[[length(plot_list)+1]] = plot1
	
	plot2 <- EnhancedVolcano(results, 
		title = paste0("case gene expression (baseline):\n",dataset1," (n=",length1,") vs. ",dataset2," (n=",length2,")\n(excl. genes with low counts)"), 
		lab=rownames(results), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, ylab = bquote(~-Log[10]~ 'adjusted P-value'), 
		selectLab = c(gene_tolabel2), pointSize = 2.0, labSize = 3.0, labCol = 'black',drawConnectors = TRUE, widthConnectors = 0.4, colConnectors = 'black',
		max.overlaps = 30,  maxoverlapsConnectors = 30) #max number of overlaps allowed)

	plot_list2[[length(plot_list2)+1]] = plot2
	
	
	# # # # # And select genes:
	results_selectgenes <- results[rownames(results) %in% c(genes_select,paste0(genes_select,"_v2")),]
		
	# add new col for padj by no. genes analysed (benj-hoch):
	results_selectgenes$padj_select_genes_BH <- p.adjust(results_selectgenes$pvalue,method="BH")
	# save select gene results:
	write.csv(results_selectgenes, 
	paste0("dseq2_cases_",length1,dataset1,"_vs_",length2,dataset2,"_Filt_baseline_Selectgenes_results.csv"), quote=F, row.names=T)
	
}

# Save plots to tiff. Makes a separate file for each plot.
for(i in 6:nrow(PT_analyses)){
	dataset1 = PT_analyses[i,1]
	dataset2 = PT_analyses[i,2]
	print(paste0("On image saving for run: ", i, " for ", dataset1, " vs. ", dataset2))
    jpeg(filename=paste0("dseq2_cases_",dataset1,"_vs_",dataset2,"_Filt_baseline_Allgenes_Volcano_FC1.jpeg"))
    print(plot_list[[i-5]])
    dev.off()
}

for(i in 6:nrow(PT_analyses)){
	dataset1 = PT_analyses[i,1]
	dataset2 = PT_analyses[i,2]
	print(paste0("On image saving for run: ", i, " for ", dataset1, " vs. ", dataset2))
    jpeg(filename=paste0("dseq2_cases_",dataset1,"_vs_",dataset2,"_Filt_baseline_Allgenes_Volcano_FC2.jpeg"))
    print(plot_list2[[i-5]])
    dev.off()
}

