# This R script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
# This script is based on the following:
# -gaCCA: Seoane JA, Campbell C, Day IN, Casas JP, Gaunt TR. Canonical correlation analysis for gene-based pleiotropy discovery.
# PLoS Comput Biol. 2014 Oct 16;10(10):e1003876. doi: 10.1371/journal.pcbi.1003876. PMID: 25329069; PMCID: PMC4199483.
# -gene-based CCA: Ferreira MA, Purcell SM. A multivariate test of association. Bioinformatics. 2009 Jan 1;25(1):132-3.
# doi: 10.1093/bioinformatics/btn563. Epub 2008 Nov 19. PMID: 19019849.
# Tang CS, Ferreira MA. A gene-based test of association using canonical correlation analysis. Bioinformatics. 2012 Mar 15;28(6):845-50.
# doi: 10.1093/bioinformatics/bts051. Epub 2012 Jan 31. PMID: 22296789.
# Original source can be found in Manuel Ferreira homepage: http://genepi.qimr.edu.au/staff/manuelF/gene/main.html

library(sqldf)
source("gene-based.cca3.2.R")

phenD_all <- read.csv("phen_data_file.csv",header=T)
row.names(phenD) <- phenD$sampleID

inputData <- read.table("GT_data_file.txt",header=T,skip=1)
inputData_forLookup <- inputData[,c(1:5)]
hg38_coords = read.table("~/hg38/glist-hg38",header=F)
names(hg38_coords) <- c("Chr","Pos_from","Pos_to","Gene")

inputData_forLookup$Gene_hg38 <- NA
Data_gene <- sqldf("select * from inputData_forLookup left join hg38_coords
                on (inputData_forLookup.POS >= hg38_coords.Pos_from and inputData_forLookup.POS <= hg38_coords.Pos_to and
                    inputData_forLookup.CHROM == hg38_coords.Chr ) ")

# remove SNPs not associated with any gene
idnil = which(is.na(Data_gene[,10]))
snps_locusID = Data_gene[-idnil,]

genes = unique(snps_locusID[,10])
l.genes = length(genes)

# Transpose so that rownames=samples and RSID=colnames
inputData_torun <- inputData[,-c(1:2,4:9)]
rownames(inputData_torun) <- inputData_torun$ID
inputData_torun <- t(inputData_torun[,-c(1)])
inputData_torun <- inputData_torun[rownames(inputData_torun) %in% phenD$sampleID,]
phenD_torun <- phenD[phenD$sampleID %in% rownames(inputData_torun),c(4:5)]

results <- data.frame(no_SNPs = numeric(), SNPs = character(), r_1 = numeric(), pval=numeric(), gene=character(), stringsAsFactors=FALSE)

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

# iterate for each gene
for(i in 1:l.genes){
	cat(genes[i])
	geneName = genes[i]
	cat("\n")

	idSnp = as.character(snps_locusID[snps_locusID$Gene==geneName,3])
	geneD = as.matrix(inputData_torun[,colnames(inputData_torun) %in% idSnp]) #retrieve genetic data of each gene

	geneD.l = dim(geneD)[2]
	phenD.l = dim(phenD_torun)[2]

	verbose = 0
	nind=nrow(geneD)
	nsnps=ncol(geneD)

	### Gene-based test
  all.can = as.numeric(cancor(dna.pruned,phenD_torun)$cor)[1]
  #all.can = as.numeric(kcca(dna.pruned,phenD_torun)$cor)
  all.can[all.can>0.99]=0.99
  ntraits=max(ncol(phenD_torun),1,na.rm=T);
  p = fapp(nsnps,ntraits,nind,all.can)

	# save results in table:
	results[i,1] = length(idSnp)
	results[i,2] = paste(idSnp, collapse="_")
	results[i,3] = all.can
	results[i,4] = p
	results[i,5] = geneName

}

# Apply multiple testing correction:
results$p_BH <- p.adjust(results$pval,method="BH")
results$p_Bonf <- p.adjust(results$pval,method="bonferroni")
results[results$p_BH < 0.05,]

write.table(results,"results_gene.txt",col.names=T,row.names=F,quote=F)


# ======= Univariate-SNP =======
results_singleSNP <- data.frame(SNP = character(), r_1 = numeric(), pval=numeric(), gene=character(), stringsAsFactors=FALSE)

for(i in 1:l.snps){
	SNP = snps[i]
	cat("\n")

	idSnp = as.character(snps_locusID[snps_locusID$ID==SNP,3])
	gene_for_snp = as.character(snps_locusID[snps_locusID$ID==SNP,10])
	geneD = as.matrix(inputData_torun[,colnames(inputData_torun) %in% idSnp]) #retrieve genetic data of each SNP
	if(sd(geneD)!=0){
		geneD.l = dim(geneD)[2]
		phenD.l = dim(phenD_torun)[2]

		verbose = 0
		nind=nrow(geneD)
		nsnps=ncol(geneD)

		keep.dna=check_mc(geneD,verbose) #remove colinearity
		dna.pruned=geneD[,keep.dna]

		### Gene-based test
		  all.can = as.numeric(cancor(dna.pruned,phenD_torun)$cor)[1]
		  all.can[all.can>0.99]=0.99
		  ntraits=max(ncol(phenD_torun),1,na.rm=T);

		  p = fapp(nsnps,ntraits,nind,all.can)

		results_singleSNP[i,1] = paste(idSnp, collapse="_")
		results_singleSNP[i,2] = all.can
		results_singleSNP[i,3] = p
		results_singleSNP[i,4] = paste(gene_for_snp, collapse="_")
	}
}

results_singleSNP$p_BH <- p.adjust(results_singleSNP$pval,method="BH")
results_singleSNP$p_Bonf <- p.adjust(results_singleSNP$pval,method="bonferroni")
results_singleSNP[results_singleSNP$p_Bonf<0.05,]
results_singleSNP[results_singleSNP$p_BH<0.05,]
write.table(results_singleSNP,"results_SNP.txt",col.names=T,row.names=F,quote=F)
