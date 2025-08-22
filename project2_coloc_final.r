library("dplyr")
library("ggplot2")
library("readr")
library("coloc")
library("GenomicRanges")
library("seqminer")
library("gwasvcf")

#remotes::install_github("mrcieu/gwasvcf")
# Skipping 9 packages not available: 
# VariantAnnotation, 
# SummarizedExperiment, S4Vectors, Rsamtools, IRanges, GenomicRanges, GenomeInfoDb, Biostrings, BiocGenerics

tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Try NURTURE-CKD CCA GTs out - since can derive all variabels needed without assuming SD of PT in CKDGen ==1
NCKD_data <- read.table("data.txt",header=T)
rep943snps <- read.table("data.SNPlist",header=F)
NCKD_data.sig <- NCKD_data[NCKD_data$SNP %in% rep943snps$V1,]
nrow(NCKD_data.sig)
# 943

# ==== In Bash: add MAF data for lead snps: rs3132450	 rs79350521	 rs1064994	 rs28630740	 rs2395231
plink2 --vcf ~/GWAS_N_APRIL2021/Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.dbSNP153IDs.maf.hwe.vcf.gz \
	--snps rs3132450, rs79350521, rs1064994, rs28630740, rs2395231  \
	--keep-allele-order \
	--double-id \
	--recode vcf \
	--out ~/GWAS_N_APRIL2021/Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.maf.hwe_MN943_5leadSNPs.vcf

# min: 6:31596138:A:G 
# max:  6:32637302:A:G

# also extract for 6:31596138 -200K, and 6:32637302 +200k
plink2 --vcf ~/GWAS_N_APRIL2021/Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.dbSNP153IDs.maf.hwe.vcf.gz \
	--chr 6 \
	--from-bp 31396138 \
	--to-bp 32837302 \
	--keep-allele-order \
	--double-id \
	--recode vcf \
	--out ~/GWAS_N_APRIL2021/Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.maf.hwe_MN_minmax_200k.vcf

# ============================

NCKD_data$uniqID <- paste0(NCKD_data$chr, ":", NCKD_data$pos, ":", NCKD_data$ref, ":", NCKD_data$alt)
NCKD_data$id <- paste0(NCKD_data$chr, ":", NCKD_data$pos)
NCKD_data$rsID <- NCKD_data$SNP

# get MAF data:
snps10 <- read.table("Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.maf.hwe_MN943_5leadSNPs.vcf.vcf", comment.char = "", skip = 30, header=T)
snps10[,c(10:ncol(snps10))] <- as.data.frame(sapply(snps10[,c(10:ncol(snps10))], function(x) gsub("0/0", 0, x)))
snps10[,c(10:ncol(snps10))] <- as.data.frame(sapply(snps10[,c(10:ncol(snps10))], function(x) gsub("0/1", 1, x)))
snps10[,c(10:ncol(snps10))] <- as.data.frame(sapply(snps10[,c(10:ncol(snps10))], function(x) gsub("1/0", 1, x)))
snps10[,c(10:ncol(snps10))] <- as.data.frame(sapply(snps10[,c(10:ncol(snps10))], function(x) gsub("1/1", 2, x)))  
snps10[,c(10:ncol(snps10))] <- as.data.frame(sapply(snps10[,c(10:ncol(snps10))], as.numeric))
AN <- ncol(snps10[,c(10:ncol(snps10))])*2
snps10$AC <- rowSums(snps10[,c(10:ncol(snps10))])
snps10$AF <- snps10$AC/AN
snps10.maf <- snps10[,c(3,2564)]
names(snps10.maf) <- c("SNP","MAF")

snps1_range <- read.table("Genotyping_data.GTS.DUP.excl.keepCtrls.maf.geno.mind.hwe2.king.EUR.beagle.hg38_1000GP.R2_0.9.maf.hwe_MN_minmax_200k.vcf.vcf", comment.char = "", skip = 30, header=T)
snps1_range[,c(10:ncol(snps1_range))] <- as.data.frame(sapply(snps1_range[,c(10:ncol(snps1_range))], function(x) gsub("0/0", 0, x)))
snps1_range[,c(10:ncol(snps1_range))] <- as.data.frame(sapply(snps1_range[,c(10:ncol(snps1_range))], function(x) gsub("0/1", 1, x)))
snps1_range[,c(10:ncol(snps1_range))] <- as.data.frame(sapply(snps1_range[,c(10:ncol(snps1_range))], function(x) gsub("1/0", 1, x)))
snps1_range[,c(10:ncol(snps1_range))] <- as.data.frame(sapply(snps1_range[,c(10:ncol(snps1_range))], function(x) gsub("1/1", 2, x)))  
snps1_range[,c(10:ncol(snps1_range))] <- as.data.frame(sapply(snps1_range[,c(10:ncol(snps1_range))], as.numeric))
AN <- ncol(snps1_range[,c(10:ncol(snps1_range))])*2
snps1_range$AC <- rowSums(snps1_range[,c(10:ncol(snps1_range))], na.rm=TRUE)
snps1_range$AF <- snps1_range$AC/AN
snps1_range.maf <- snps1_range[,c(3,2564)]
names(snps1_range.maf) <- c("SNP","MAF")

all.MAF <- rbind(snps10.maf, snps1_range.maf)

NCKD_data.MAF <- merge(NCKD_data, all.MAF, by = c("SNP"))

leadsnps <- read.table("leadSNPs.txt", header=T)
NCKD_data.MAF.leadsnps <- merge(NCKD_data.MAF, leadsnps[,c(1,2,3,4,7)], by = c("rsID"))

# for sdY and N:
SDY = 1 #(since normalised eGFR and BUN)
N = 2181


# uniqID.x = hg38, uniqID.y = hg19
eqtls_fuma <- read.table("eqtl.txt", header=T)
#eqtls_fuma.relgtex.formerge <- unique(eqtls_fuma[eqtls_fuma$tissue %in% c("Whole_Blood","Kidney_Cortex","Liver","Cells_EBV-transformed_lymphocytes"),c(1,4,13)])

# fairer to include all tissues:
eqtls_fuma.tomerge <- eqtls_fuma[,c(1,4,13)]
names(eqtls_fuma.tomerge) <- c("uniqID.y","gene","symbol")
NCKD_data.MAF.leadsnps.eqtl <- merge(NCKD_data.MAF.leadsnps, eqtls_fuma.tomerge, by = c("uniqID.y"), all.x=T)


run_mycoloc_SDY1 <- function(eqtl_sumstats, gwas_sumstats){
    eQTL_dataset = list(beta = eqtl_sumstats$beta,
                        varbeta = eqtl_sumstats$se^2,
                        N = (eqtl_sumstats$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                        MAF = eqtl_sumstats$maf, 
                        type = "quant", 
                        snp = eqtl_sumstats$id)
  gwas_dataset = list(MAF = as.numeric(datachr$MAF), 
					  snp = datachr$id,
                      position = datachr$pos,
					  type = "quant", 
					  sdY=1,
                      N = 916,
					  pvalues = datachr$pval)
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}


#3/Imported eQTL datasets (Currently GTEx_v8 only)
rnaseq_df = dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
  dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
ftp_path_list = setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)

#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
      print(ftp_path)
  }
  
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  
  #Remove rsid duplicates and multi-allelic variant
  #summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%  dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
}

#Wrap the download function around purrr::safely to avoid catch erros
safe_import = purrr::safely(import_eQTLCatalogue)

NCKD_data.MAF <- NCKD_data.MAF[!duplicated(NCKD_data.MAF$SNP),]


# for each lead snp:
for(i in 1:length(sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP)))){
	print(i)
	snp = sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP))[i]
	chr = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$SNP %in% snp,6])
	pos_plus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])+200000
	pos_minus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])-200000
	
	# now get all snps (not just lead) run in GWAS for region - defined as +/- 200kB of lead snp:
	datachr <- NCKD_data.MAF %>% filter(chr %in% chr) %>% filter(pos >= pos_minus200kb & pos <= pos_plus200kb)
	region <- unique(paste0(datachr$chr, ":", min(datachr$pos), "-", max(datachr$pos)))
	
	# can't have duplciated snps (i.e. snp cant be associated with >1 gene), so if any dup. snps, then got to run each unique snp-gene separately
	#dup_snps <- datachr[duplicated(datachr$rsID),]
	#lead_datachr = nosdY_data.ckdgen.metacca.hg38.leadsnps.eqtl %>% filter(Chr %in% chr) %>% filter(hg38_start >= pos_minus200kb & hg38_start <= pos_plus200kb)
	genes = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,17]))
	symbols = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,18]))
	for(l in 1:length(unique(genes[!is.na(genes)]))){
		gene_torun = genes[!is.na(genes)][l]

		#Import summmary stats
		summary_list2 = purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = gene_torun, column_names))
		#Extract successful results
		result_list4 = purrr::map(summary_list2, ~.$result)
		result_list4 = result_list4[!unlist(purrr::map(result_list4, is.null))]
		
		#Remove rows that have NAs for standard error
		result_filtered4 = purrr::map(result_list4, ~dplyr::filter(., !is.na(se)))
		# run coloc
		
		result_filtered5 <- result_filtered4[sapply(result_filtered4, nrow)>0]
		if(length(result_filtered5)>0){
			coloc_df_imported_out = purrr::map_df(result_filtered5, ~run_mycoloc_SDY1(., datachr), .id = "qtl_id")
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- coloc_df_imported_out$PP.H4.abf/(coloc_df_imported_out$PP.H3.abf + coloc_df_imported_out$PP.H4.abf)
		} else {
			coloc_df_imported_out = data.frame(qtl_id = NA)
			coloc_df_imported_out$nsnps	= NA
			coloc_df_imported_out$PP.H0.abf = NA	
			coloc_df_imported_out$PP.H1.abf = NA
			coloc_df_imported_out$PP.H2.abf	= NA
			coloc_df_imported_out$PP.H3.abf = NA	
			coloc_df_imported_out$PP.H4.abf = NA
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- NA
		}
		if(i == 1 & l==1){
			all_results_v6 <- coloc_df_imported_out
			
		} else{
			all_results_v6 <- rbind(coloc_df_imported_out, all_results_v6)

		}
	}
}

head(all_results_v6[order(-all_results_v6$PP.H4.abf),c(1,7,8,11)],n=30)


write.table(all_results_v6, "gtex_coloc.txt")

kidney <- all_results_v6[all_results_v6$qtl_id %in% "GTEx_V8_Kidney_Cortex",]
# qtl_id                PP.H4.abf gene     snp
# GTEx_V8_Kidney_Cortex     0.944 C4A      rs2395231
# GTEx_V8_Kidney_Cortex     0.944 C4A      rs3132450
# GTEx_V8_Kidney_Cortex     0.942 C4A      rs1064994
# GTEx_V8_Kidney_Cortex     0.928 C4A      rs79350521

# 1/Microarray datasets
microarray_df = dplyr::filter(tabix_paths, quant_method == "microarray") %>%
  dplyr::mutate(qtl_id = paste(study_label, sample_group, sep = "_"))
ftp_path_list = setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id)

#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

for(i in 1:length(sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP)))){
	print(i)
	snp = sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP))[i]
	chr = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$SNP %in% snp,6])
	pos_plus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])+200000
	pos_minus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])-200000
	
	# now get all snps (not just lead) run in GWAS for region - defined as +/- 200kB of lead snp:
	datachr <- NCKD_data.MAF %>% filter(chr %in% chr) %>% filter(pos >= pos_minus200kb & pos <= pos_plus200kb)
	region <- unique(paste0(datachr$chr, ":", min(datachr$pos), "-", max(datachr$pos)))
	
	# can't have duplciated snps (i.e. snp cant be associated with >1 gene), so if any dup. snps, then got to run each unique snp-gene separately
	#dup_snps <- datachr[duplicated(datachr$rsID),]
	#lead_datachr = nosdY_data.ckdgen.metacca.hg38.leadsnps.eqtl %>% filter(Chr %in% chr) %>% filter(hg38_start >= pos_minus200kb & hg38_start <= pos_plus200kb)
	genes = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,17]))
	symbols = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,18]))
	for(l in 1:length(unique(genes[!is.na(genes)]))){
		gene_torun = genes[!is.na(genes)][l]

		#Import summmary stats
		summary_list2 = purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = gene_torun, column_names))
		#Extract successful results
		result_list4 = purrr::map(summary_list2, ~.$result)
		result_list4 = result_list4[!unlist(purrr::map(result_list4, is.null))]
		
		#Remove rows that have NAs for standard error
		result_filtered4 = purrr::map(result_list4, ~dplyr::filter(., !is.na(se)))
		# run coloc
		
		result_filtered5 <- result_filtered4[sapply(result_filtered4, nrow)>0]
		if(length(result_filtered5)>0){
			coloc_df_imported_out = purrr::map_df(result_filtered5, ~run_mycoloc_SDY1(., datachr), .id = "qtl_id")
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- coloc_df_imported_out$PP.H4.abf/(coloc_df_imported_out$PP.H3.abf + coloc_df_imported_out$PP.H4.abf)
		} else {
			coloc_df_imported_out = data.frame(qtl_id = NA)
			coloc_df_imported_out$nsnps	= NA
			coloc_df_imported_out$PP.H0.abf = NA	
			coloc_df_imported_out$PP.H1.abf = NA
			coloc_df_imported_out$PP.H2.abf	= NA
			coloc_df_imported_out$PP.H3.abf = NA	
			coloc_df_imported_out$PP.H4.abf = NA
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- NA
		}
		if(i == 1 & l==1){
			all_results_v7 <- coloc_df_imported_out
			
		} else{
			all_results_v7 <- rbind(coloc_df_imported_out, all_results_v7)

		}
	}
}

head(all_results_v7[order(-all_results_v7$PP.H4.abf),c(1,7,8,11)],n=30)

write.table(all_results_v7, "microarray_coloc.txt")


#2/Uniformly processed RNA-seq datasets
rnaseq_df = dplyr::filter(tabix_paths, quant_method == "ge") %>%
  dplyr::mutate(qtl_id = paste(study_label, sample_group, sep = "_"))
ftp_path_list = setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)

#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

for(i in 1:length(sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP)))){
	print(i)
	snp = sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP))[i]
	chr = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$SNP %in% snp,6])
	pos_plus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])+200000
	pos_minus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])-200000
	
	# now get all snps (not just lead) run in GWAS for region - defined as +/- 200kB of lead snp:
	datachr <- NCKD_data.MAF %>% filter(chr %in% chr) %>% filter(pos >= pos_minus200kb & pos <= pos_plus200kb)
	region <- unique(paste0(datachr$chr, ":", min(datachr$pos), "-", max(datachr$pos)))
	
	# can't have duplciated snps (i.e. snp cant be associated with >1 gene), so if any dup. snps, then got to run each unique snp-gene separately
	#dup_snps <- datachr[duplicated(datachr$rsID),]
	#lead_datachr = nosdY_data.ckdgen.metacca.hg38.leadsnps.eqtl %>% filter(Chr %in% chr) %>% filter(hg38_start >= pos_minus200kb & hg38_start <= pos_plus200kb)
	genes = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,17]))
	symbols = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,18]))
	for(l in 1:length(unique(genes[!is.na(genes)]))){
		gene_torun = genes[!is.na(genes)][l]

		#Import summmary stats
		summary_list2 = purrr::map(ftp_path_list, ~safe_import(., region, selected_gene_id = gene_torun, column_names))
		#Extract successful results
		result_list4 = purrr::map(summary_list2, ~.$result)
		result_list4 = result_list4[!unlist(purrr::map(result_list4, is.null))]
		
		#Remove rows that have NAs for standard error
		result_filtered4 = purrr::map(result_list4, ~dplyr::filter(., !is.na(se)))
		# run coloc
		
		result_filtered5 <- result_filtered4[sapply(result_filtered4, nrow)>0]
		if(length(result_filtered5)>0){
			coloc_df_imported_out = purrr::map_df(result_filtered5, ~run_mycoloc_SDY1(., datachr), .id = "qtl_id")
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- coloc_df_imported_out$PP.H4.abf/(coloc_df_imported_out$PP.H3.abf + coloc_df_imported_out$PP.H4.abf)
		} else {
			coloc_df_imported_out = data.frame(qtl_id = NA)
			coloc_df_imported_out$nsnps	= NA
			coloc_df_imported_out$PP.H0.abf = NA	
			coloc_df_imported_out$PP.H1.abf = NA
			coloc_df_imported_out$PP.H2.abf	= NA
			coloc_df_imported_out$PP.H3.abf = NA	
			coloc_df_imported_out$PP.H4.abf = NA
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			coloc_df_imported_out$cond_post_prob <- NA
		}
		if(i == 1 & l==1){
			all_results_v8 <- coloc_df_imported_out
			
		} else{
			all_results_v8 <- rbind(coloc_df_imported_out, all_results_v8)

		}
	}
}

head(all_results_v8[order(-all_results_v8$PP.H4.abf),c(1,7,8,11)],n=30)

write.table(all_results_v8, "rnaseq.txt")



# # # # # For nephqtl datasets:
# 0. create lists of pos and genes to query with outside of R (using unix)
write.table(unique(NCKD_data.MAF.leadsnps.eqtl$gene), "leadSNPanalyses.eqtl.hg38.noNAgenes_geneENSG.txt", col.names=F, quote=F, row.names=F)

####=================
# on Unix:
cd ~/nephqtl
grep -wFf leadSNPanalyses.eqtl.hg38.noNAgenes_geneENSG.txt eQTLs_chrAll_50_peers_cis1000kb_Tube.NephQTL2.txt \
> eQTLs_chrAll_50_peers_cis1000kb_Tube.NephQTL2_genes.nckd_MNeph_5leadSNPanalyses.txt
grep -wFf leadSNPanalyses.eqtl.hg38.noNAgenes_geneENSG.txt eQTLs_chrAll_40_peers_cis1000kb_Glom.NephQTL2.txt \
> eQTLs_chrAll_40_peers_cis1000kb_Glom.NephQTL2_genes.nckd_MNeph_5leadSNPanalyses.txt

#=====================

nephqtl_tube <- read.table("eQTLs_chrAll_50_peers_cis1000kb_Tube.NephQTL2_genes.nckd_MNeph_5leadSNPanalyses.txt", header=F)
names(nephqtl_tube) <- c("gene","SNP_ID","CHR","POS","REF","ALT","AltFreq","beta","se","t-stat","p-value","FDR")


# 1. convert to hg38
library(rtracklayer)
hg19.gr<-GRanges(
	seqname=Rle(paste("chr",nephqtl_tube$CHR,sep="")),ranges=IRanges(start=nephqtl_tube$POS,end=nephqtl_tube$POS),snp.name=nephqtl_tube$SNP_ID)

chain.file.path<-'~/coloc/hg19ToHg38.over.chain'

c<-import.chain(chain.file.path) ## e.g. hg19ToHg38.over.chain
hg38.gr<-unlist(liftOver(hg19.gr,c))  
names(mcols(hg38.gr))<-'snp.name'

hg38.df <- data.frame(hg38.gr)
names(hg38.df) <- c("chr", "hg38_start","hg38_end","width","strand","SNP_ID")
hg38.df$chr <- gsub("chr", "", hg38.df$chr)

# 2. define id with hg38 coords
nephqtl_tube.hg38 <- merge(nephqtl_tube, hg38.df, by = c("SNP_ID"), all.x=T)
nephqtl_tube.hg38$uniqID <- paste0(nephqtl_tube.hg38$chr, ":", nephqtl_tube.hg38$hg38_start,":", toupper(nephqtl_tube.hg38$REF),":", toupper(nephqtl_tube.hg38$ALT))
nephqtl_tube.hg38$id <- paste0(nephqtl_tube.hg38$chr, ":", nephqtl_tube.hg38$hg38_start)

### repeat for glom:
nephqtl_glom <- read.table("eQTLs_chrAll_40_peers_cis1000kb_Glom.NephQTL2_genes.nckd_MNeph_5leadSNPanalyses.txt", header=F)
names(nephqtl_glom) <- c("gene","SNP_ID","CHR","POS","REF","ALT","AltFreq","beta","se","t-stat","p-value","FDR")
hg19.gr<-GRanges(
	seqname=Rle(paste("chr",nephqtl_glom$CHR,sep="")),ranges=IRanges(start=nephqtl_glom$POS,end=nephqtl_glom$POS),snp.name=nephqtl_glom$SNP_ID)
c<-import.chain(chain.file.path) ## e.g. hg19ToHg38.over.chain
hg38.gr<-unlist(liftOver(hg19.gr,c))  
names(mcols(hg38.gr))<-'snp.name'
hg38.df <- data.frame(hg38.gr)
names(hg38.df) <- c("chr", "hg38_start","hg38_end","width","strand","SNP_ID")
hg38.df$chr <- gsub("chr", "", hg38.df$chr)

nephqtl_glom.hg38 <- merge(nephqtl_glom, hg38.df, by = c("SNP_ID"), all.x=T)
nephqtl_glom.hg38$uniqID <- paste0(nephqtl_glom.hg38$chr, ":", nephqtl_glom.hg38$hg38_start,":", toupper(nephqtl_glom.hg38$REF),":", toupper(nephqtl_glom.hg38$ALT))
nephqtl_glom.hg38$id <- paste0(nephqtl_glom.hg38$chr, ":", nephqtl_glom.hg38$hg38_start)

'%notin%' <- Negate('%in%')

run_mycoloc_nephqtlglom <- function(eqtl_sumstats, gwas_sumstats){
    eQTL_dataset = list(beta = nephqtl_glom.hg38.torun$beta,
                        varbeta = nephqtl_glom.hg38.torun$se^2,
                        N = 136, # Samples size is allele number (AN) dvided by 2
                        MAF = nephqtl_glom.hg38.torun$AltFreq, 
                        type = "quant",
						snp = nephqtl_glom.hg38.torun$uniqID)
  gwas_dataset = list(MAF = as.numeric(datachr$MAF), 
					  snp = datachr$uniqID,
                      position = datachr$pos,
					  type = "quant", 
					  sdY=1,
                      N = 916,
					  pvalues = datachr$pval)
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}


run_mycoloc_nephqtltub <- function(eqtl_sumstats, gwas_sumstats){
    eQTL_dataset = list(beta = nephqtl_tube.hg38.torun$beta,
                        varbeta = nephqtl_tube.hg38.torun$se^2,
                        N = 166, # Samples size is allele number (AN) dvided by 2
                        MAF = nephqtl_tube.hg38.torun$AltFreq, 
                        type = "quant",
						snp = nephqtl_tube.hg38.torun$uniqID)
  gwas_dataset = list(MAF = as.numeric(datachr$MAF), 
					  snp = datachr$uniqID,
                      position = datachr$pos,
					  type = "quant", 
					  sdY=1,
                      N = 916,
					  pvalues = datachr$pval)
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}


for(i in 1:length(sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP)))){
	print(i)
	snp = sort(unique(NCKD_data.MAF.leadsnps.eqtl$SNP))[i]
	chr = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$SNP %in% snp,6])
	pos_plus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])+200000
	pos_minus200kb = unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,7])-200000
	
	# now get all snps (not just lead) run in GWAS for region - defined as +/- 200kB of lead snp:
	datachr <- unique(NCKD_data.MAF %>% filter(chr %in% chr) %>% filter(pos >= pos_minus200kb & pos <= pos_plus200kb))
	region <- unique(paste0(datachr$chr, ":", min(datachr$pos), "-", max(datachr$pos)))
	
	# can't have duplciated snps (i.e. snp cant be associated with >1 gene), so if any dup. snps, then got to run each unique snp-gene separately
	#dup_snps <- datachr[duplicated(datachr$rsID),]
	#lead_datachr = nosdY_data.ckdgen.metacca.hg38.leadsnps.eqtl %>% filter(Chr %in% chr) %>% filter(hg38_start >= pos_minus200kb & hg38_start <= pos_plus200kb)
	genes = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,17]))
	symbols = c(unique(NCKD_data.MAF.leadsnps.eqtl[NCKD_data.MAF.leadsnps.eqtl$rsID %in% snp,18]))
	for(l in 1:length(unique(genes[genes %notin% NA]))){
		gene_torun = genes[genes %notin% NA][l]
		print(gene_torun)

		# get all snps in region:
		# tub:
		#nephqtl_tube.hg38.torun <- unique(nephqtl_tube.hg38[(nephqtl_tube.hg38$CHR %in% chr & 
		#	nephqtl_tube.hg38$hg38_start >= min(datachr$pos) & 
		#		nephqtl_tube.hg38$hg38_end <= max(datachr$pos) & nephqtl_tube.hg38$gene %in% gene_torun),])
		#nephqtl_tube.hg38.torun <- na.omit(nephqtl_tube.hg38.torun)
		#if(nrow(nephqtl_tube.hg38.torun)>0){
		
		# glom:
		nephqtl_glom.hg38.torun <- unique(nephqtl_glom.hg38[(nephqtl_glom.hg38$chr %in% chr & nephqtl_glom.hg38$hg38_start >= min(datachr$pos) & 
			nephqtl_glom.hg38$hg38_end <= max(datachr$pos) & nephqtl_glom.hg38$gene %in% gene_torun),])
		nephqtl_glom.hg38.torun <- na.omit(nephqtl_glom.hg38.torun)
		if(nrow(nephqtl_glom.hg38.torun)>0){
			#coloc_df_imported_out = purrr::map_df(nephqtl_tube.hg38.torun, ~run_mycoloc_nephqtltub(., datachr), .id = "uniqID")
			coloc_df_imported_out = purrr::map_df(nephqtl_glom.hg38.torun, ~run_mycoloc_nephqtlglom(., datachr), .id = "uniqID")
			#coloc_df_imported_out = run_mycoloc_nephqtltub(nephqtl_tube.hg38.torun, datachr)
			coloc_df_imported_out = run_mycoloc_nephqtlglom(nephqtl_glom.hg38.torun, datachr)
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			#coloc_df_imported_out$cond_post_prob <- coloc_df_imported_out$PP.H4.abf/(coloc_df_imported_out$PP.H3.abf + coloc_df_imported_out$PP.H4.abf)
		} else {
			#coloc_df_imported_out$qtl_id = NA
			coloc_df_imported_out$nsnps	= NA
			coloc_df_imported_out$PP.H0.abf = NA	
			coloc_df_imported_out$PP.H1.abf = NA
			coloc_df_imported_out$PP.H2.abf	= NA
			coloc_df_imported_out$PP.H3.abf = NA	
			coloc_df_imported_out$PP.H4.abf = NA
			coloc_df_imported_out$gene = symbols[l]
			coloc_df_imported_out$gene_en = gene_torun
			coloc_df_imported_out$region = region
			coloc_df_imported_out$snp = snp
			#coloc_df_imported_out$cond_post_prob <- NA
		}
		if(i == 1 & l==1){
				all_results_v19 <- coloc_df_imported_out
						
		} else{
			all_results_v19 <- rbind(coloc_df_imported_out, all_results_v19)

		}		
	}
	
}

print(all_results_v19[order(-all_results_v19$PP.H4.abf),],n=30)
#write.table(all_results_v19, "coloc_MN_943_5leadsnps_usingFUMA_nephqtltub_datasets.txt", quote=F, row.names=F, col.names=T)
write.table(all_results_v19, "coloc_MN_943_5leadsnps_usingFUMA_nephqtlglom_datasets.txt", quote=F, row.names=F, col.names=T)

##### Plots:
# # # # # # # #  Create plot of PP.H4 (y-axis) vs. chr. loc. for novel lead SNPs only 
ckdgen.lead.coloc <- read.csv("leadSNPs_eqtl_coloc_alltissues.csv")
names(ckdgen.lead.coloc) <- c("qtl_id","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","gene","gene_en","region","SNP", "category")
ckdgen_leadsnps <- read.csv("leadsnps.info.forpaper.csv")
ckdgen_leadsnps$SNP <- ckdgen_leadsnps$rsID
ckdgen_leadsnps.coloc <- merge(ckdgen_leadsnps, ckdgen.lead.coloc, by = c("SNP"), all.x=T)
# 472 snps

gwas.dat2.lead.coloc <- merge(gwas.dat2, ckdgen_leadsnps.coloc , by = c("SNP","chr","pos"), all.y=T)
length(unique(gwas.dat2.lead.coloc[gwas.dat2.lead.coloc$also_lead_in_Wuttke %in% NA,1]))
# 394
length(unique(gwas.dat2.lead.coloc[gwas.dat2.lead.coloc$sig.inWuttke %in% "no" & gwas.dat2.lead.coloc$egfr_bun_eff_oppdir %notin% NA,1]))
# 157

gwas.dat2.lead.coloc_toplot <- gwas.dat2.lead.coloc[gwas.dat2.lead.coloc$sig.inWuttke %in% "no" & gwas.dat2.lead.coloc$egfr_bun_eff_oppdir %notin% NA,]
gwas.dat2.lead.coloc_toplot$genelabel <- ifelse(gwas.dat2.lead.coloc_toplot$PP.H4.abf>=0.8, gwas.dat2.lead.coloc_toplot$gene, "")
gwas.dat2.lead.coloc_toplot <- gwas.dat2.lead.coloc_toplot %>%
	group_by(SNP,gene) %>%
	mutate(genelabel_toshow_notissue = ifelse((PP.H4.abf == max(PP.H4.abf) & PP.H4.abf>=0.8), genelabel, ""))

plot_eqtl <- ggplot(gwas.dat2.lead.coloc_toplot, aes(x = BPcum, y = PP.H4.abf, color = as.factor(chr), group = as.factor(chr))) +
  geom_point(alpha = 0.5,stroke = 0, shape = 16,size=1,show.legend = F) +
  geom_hline(yintercept = 0.8, color = "black", linetype = "dotted") +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand=c(0.02, 0)) +
 # scale_y_continuous(limits = c(0, 1.4), breaks = c(seq(0, 1.4, by = 0.2)),expand=c(0.04, 0)) +
 # scale_y_continuous(limits = c(0, 1), breaks = c(seq(0, 1, by = 0.2)),expand=c(0.04, 0)) +
  scale_y_continuous(breaks = c(seq(0, 1, by = 0.1)),expand=c(0.04, 0)) +
  scale_color_manual(values = rep(c("gray65", "gray20"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  #geom_point(data=gwas.dat.genes.lead[gwas.dat.genes.lead$leadsnp_dir %notin% "",], size=3, color="black", shape=22)+
  #geom_label(label=gwas.dat.genes.lead.wuttke$novelgene, nudge_y = 4, check_overlap = T, size=1.5, colour = "black")+
  geom_label_repel(label=gwas.dat2.lead.coloc_toplot$genelabel_toshow, size=1.4, colour = "black", max.overlaps = Inf)+
  # , segment.color = NA
  labs(x = "Chromosomal location",
		y = "eQTL colocalisation PP.H4")+
	   #y = expression(bold(-log[10](CCA~bolditalic(P)~value)))) +
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
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",margin = margin(t = 0, r = 12, b = 0, l = 0))
    #plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    #plot.tag = element_text(face="bold", size = 20,color="black")
)

jpeg("plot_eqtl.jpeg", pointsize=4, width=170, height=70,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot_eqtl, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()


# Try tissue type (y-axis) by gene (x-axis)
gwas.dat2.lead.coloc_toplot_sig <- gwas.dat2.lead.coloc_toplot[gwas.dat2.lead.coloc_toplot$PP.H4.abf >=0.8 & (gwas.dat2.lead.coloc_toplot$category %notin% NA),]
gwas.dat2.lead.coloc_toplot_sig <- gwas.dat2.lead.coloc_toplot_sig %>%
	group_by(gene, category) %>%
	mutate(count_snp  = n_distinct(SNP))

plot_eqtl2 <- ggplot(gwas.dat2.lead.coloc_toplot_sig, aes(x = gene, y = category, colour = PP.H4.abf, size = as.factor(count_snp))) +
  geom_point(alpha = 0.5, stroke = 0, shape = 16, show.legend = T) +
  #scale_color_distiller(type = "seq", direction = 1,  palette = "Greys") +
  #scale_color_grey(start=0.3, end=0.9)+
  scale_size_manual(values = c("1" = 3, "2"=4), name = "Number of\ngenetic variants")+
  scale_color_gradient(low = "grey84", high = "black", name = "PP of\ncolocalisation")+
  labs(x = "Gene",	y = "Tissue or cell type category")+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme(
    legend.text=element_text(size=10),
	legend.key.size = unit(0.5, "cm"),
	legend.spacing.y = unit(0.5, 'cm'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 7, angle=45,vjust = 1, hjust=1,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",margin = margin(t = 0, r = 12, b = 0, l = 0))
    #plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    #plot.tag = element_text(face="bold", size = 20,color="black")
)

jpeg("plot_eqtl2.jpeg", pointsize=4, width=250, height=150,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot_eqtl2, vp=viewport(layout.pos.row=1,layout.pos.col=1))

dev.off()



# # # # # # # #  Create plot of PP.H4 (y-axis) vs. chr. loc. for novel lead SNPs only
tub <- read.table("coloc_leadsnps_usingFUMA_nephqtltub_datasets_v2.txt",header=T)
glom <- read.table("coloc_leadsnps_usingFUMA_nephqtlglom_datasets_v2.txt",header=T)
gtex <- read.table("gtex_coloc_leadsnps_usingFUMA_v2.txt",header=T)
marray <- read.table("microarray_leadsnps_usingFUMA_v2-all.txt",header=T)
rnaseq <- read.table("rnaseq_coloc_leadsnps_usingFUMA_v2.txt",header=T)
tub$qtl_id <- "nephQTL_tubular"
glom$qtl_id <- "nephQTL_glomerular"
glom <- glom[,c(11,1:10)]
tub <- tub[,c(1:11)]
gtex <- gtex[,c(1:11)]
marray <- marray[,c(1:11)]
rnaseq <- rnaseq[,c(1:11)]

lead.coloc <- rbind(glom,tub,gtex,marray,rnaseq)
names(lead.coloc) <- c("qtl_id","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","gene","gene_en","region","SNP")
leadsnps <- read.table("leadSNPs.txt",header=T)
leadsnps$SNP <- leadsnps$rsID
leadsnps.coloc <- merge(leadsnps, lead.coloc, by = c("SNP"), all.x=T)
# 5 snps are all novel for MN

leadsnps.coloc$genelabel <- ifelse(leadsnps.coloc$PP.H4.abf>=0.8, leadsnps.coloc$gene, "")
leadsnps.coloc_toplot <- leadsnps.coloc %>%
	group_by(SNP,gene) %>%
	mutate(genelabel_toshow_notissue = ifelse((PP.H4.abf == max(PP.H4.abf) & PP.H4.abf>=0.8), genelabel, ""))

# get category from ckdgen ss:
ckdgencats <- read.csv("leadSNPs_eqtl_coloc_alltissues.csv", header=T)
leadsnps.coloc_toplot <- merge(leadsnps.coloc_toplot, unique(ckdgencats[,c(1,12)]), by = c("qtl_id"), all.x=T)
leadsnps.coloc_toplot$category <- ifelse(leadsnps.coloc_toplot$qtl_id %in% "NephQLT_glomerular", "NephQTL_glomerular",
	ifelse(leadsnps.coloc_toplot$qtl_id %in% "NephQTL_tubular", "NephQTL_tubular", leadsnps.coloc_toplot$category))
	
# Try tissue type (y-axis) by gene (x-axis)
leadsnps.coloc_toplot_sig <- leadsnps.coloc_toplot[leadsnps.coloc_toplot$PP.H4.abf >=0.8 & (leadsnps.coloc_toplot$category %notin% c(NA, "NA")) & !is.na(leadsnps.coloc_toplot$category) & !is.na(leadsnps.coloc_toplot$gene) & (leadsnps.coloc_toplot$gene %notin% c(NA, "NA")) ,]
	
leadsnps.coloc_toplot_sig <- leadsnps.coloc_toplot_sig %>%
	group_by(gene, category) %>%
	mutate(count_snp  = n_distinct(SNP))

leadsnps.coloc_toplot_sig.cats <- leadsnps.coloc_toplot_sig %>%
	group_by(SNP, gene) %>%
	mutate(ckd_rel = ifelse(category %in% c("nephQTL_tubular","blood","kidney","immune - B-cell","immune - monocyte","immune - macrophage","immune - other","immune - T-cell", 
		"nephQTL_glomerular","liver"), 1, 0)) %>%
	filter(gene %notin% NA)
		
unique(leadsnps.coloc_toplot_sig.cats$count_snp)
#[1] 1 4 2 3

#leadsnps.coloc_toplot_sig.cats.novelSNP <- leadsnps.coloc_toplot_sig.cats %>%
	#filter(SNP %in% c("rs79350521", "rs1064994", "rs28630740", "rs2395231"))
	
library(grid)
plot3 <- ggplot(leadsnps.coloc_toplot_sig.cats, aes(x = gene, y = category, colour = PP.H4.abf, shape = as.factor(ckd_rel))) +
  geom_point(alpha = 0.5, stroke = 0, show.legend = T, size = 5) +
  #scale_color_distiller(type = "seq", direction = 1,  palette = "Greys") +
  #scale_color_grey(start=0.3, end=0.9)+
  #scale_size_manual(values = c("1" = 3, "2"=4), name = "Number of\ngenetic variants")+
  scale_shape_manual(values = c("0" = 16, "1"= 17), labels = c('No', 'Yes'), name = "CKD relevant\ntissue")+
  scale_color_gradient(low = "grey60", high = "black", name = "PP of\ncolocalisation")+
  labs(x = "Gene",	y = "Tissue or cell type category")+
  guides(shape = guide_legend(override.aes = list(size = 2)))+
  theme(
    legend.text=element_text(size=10),
	legend.key.size = unit(0.5, "cm"),
	legend.spacing.y = unit(0.5, 'cm'),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
	panel.grid = element_line(color = "lightgrey",  size = 0.5,  linetype = 1),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 8, angle=45,vjust = 1, hjust=1,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",margin = margin(t = 0, r = 12, b = 0, l = 0))
    #plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    #plot.tag = element_text(face="bold", size = 20,color="black")
)

jpeg("coloc_leadSNP_FUMA_novel_eqtl_tissueCat_gene.cats.novelSNP_v2.jpeg", pointsize=4, width=250, height=150,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot3, vp=viewport(layout.pos.row=1,layout.pos.col=1))
dev.off()


# anyway, generally shows systemic genes involved across many tissues... but not as many immune/kidney as hoped.
