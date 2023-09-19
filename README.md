# CCA_scripts
The scripts here were developed by Amy J. Osborne and based on the following: 

-gaCCA: Seoane JA, Campbell C, Day IN, Casas JP, Gaunt TR. Canonical correlation analysis for gene-based pleiotropy discovery. PLoS Comput Biol. 2014 Oct 16;10(10):e1003876. doi: 10.1371/journal.pcbi.1003876. PMID: 25329069; PMCID: PMC4199483.

-gene-based CCA: Ferreira MA, Purcell SM. A multivariate test of association. Bioinformatics. 2009 Jan 1;25(1):132-3. doi: 10.1093/bioinformatics/btn563. Epub 2008 Nov 19. PMID: 19019849.

Tang CS, Ferreira MA. A gene-based test of association using canonical correlation analysis. Bioinformatics. 2012 Mar 15;28(6):845-50. doi: 10.1093/bioinformatics/bts051. Epub 2012 Jan 31. PMID: 22296789.
Original source can be found in Manuel Ferreira homepage: http://genepi.qimr.edu.au/staff/manuelF/gene/main.html

This software uses the following function:

-metaCCA: Cichonska A, Rousu J, Marttinen P, Kangas AJ, Soininen P, Lehtimäki T, Raitakari OT, Järvelin MR, Salomaa V, Ala-Korpela M, Ripatti S, Pirinen M. metaCCA: summary statistics-based multivariate meta-analysis of genome-wide association studies using canonical correlation analysis. Bioinformatics. 2016 Jul 1;32(13):1981-9. doi: 10.1093/bioinformatics/btw052. Epub 2016 Feb 19. PMID: 27153689; PMCID: PMC4920109.



==== Contents of scripts ====

1. Data pre-processing:

  (a) gwas_qc_filters.sh
  
  (b) checking_batch_effects_usingPCA.R

  (c) data_preparation_summary_statistics.sh 

  (d) with example files: GT_data_file.txt and phen_data_file.csv

2. Data analysis:

  (a) univariate_CCA.R

  (b) univariate_metaCCA.R

  (c) gene_based_multivariateSNP_metaCCA.R

3. Enrichment analyses and visualisation of results:

  (a) manhplots_Figs2_6.R

  (b) bar_charts_Fig3.R
  
  (c) venn_overlap_statistics_significant_genes_Fig5.R
  
  (d) enrichment_analyses_Fig4.R
