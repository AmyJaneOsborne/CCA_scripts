# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
### Dataset preparation for applying metaCCA to GWAS summary statistics datasets

# ===== 1. For BioBank Japan GWAS Summary statistics datasets (and using reference dataset: East Asian 1000GP) ====
# References:
# Kanai, M. et al. Genetic analysis of quantitative traits in the Japanese population links cell types to complex human diseases. Nat Genet 50, 390-400, doi:10.1038/s41588-018-0047-6 (2018).
# Hemani, G. et al. The MR-Base platform supports systematic causal inference across the human phenome. Elife 7, doi:10.7554/eLife.34408 (2018).

# Download published and publicly accessible GWAS summary statistics files using wget, then unzip.
wget https://gwas.mrcieu.ac.uk/files/bbj-a-60/bbj-a-60.vcf.gz
wget https://gwas.mrcieu.ac.uk/files/bbj-a-11/bbj-a-11.vcf.gz
# Count SNPs:
zcat bbj-a-60.vcf.gz | wc -l
# 6109063
zcat bbj-a-11.vcf.gz | wc -l
# 6109063

# Get same SNPs for eGFR and BUN:
zcat bbj-a-60.vcf.gz | grep -v '##' > bbj-a-60.eGFR.noheader.vcf
zcat bbj-a-11.vcf.gz | grep -v '##' > bbj-a-11.BUN.noheader.vcf
# Format ES:SE:LP:AF:ID
awk 'NR==FNR{a[$1":"$2":"$4">"$5]=" "$0; next}{ print $0 (a[$1":"$2":"$4">"$5]?a[$1":"$2":"$4">"$5]:" missing")}' OFS='\t' \
bbj-a-11.BUN.noheader.vcf \
bbj-a-60.eGFR.noheader.vcf \
| grep -v missing | awk '{print $1, $2, $3, $4, $5, $10, $20}' \
| sed 's/:/ /g' | awk '{print $1, $2, $3, $4, $5, $6, $7, $11, $12}' \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.vcf

# Download EAS 1000GP file
mkdir reference_1000GP_EAS
cd reference_1000GP_EAS
wget https://vegas2.qimrberghofer.edu.au/g1000p3_EAS.tar.gz
tar -xzf g1000p3_EAS.tar.gz

# ==== 2. For EUR 1000GP (European for CKDGen, UKHLS) =====
# Download EUR 1000GP from VEGAS2 website: https://vegas2.qimrberghofer.edu.au/
# Reference: Mishra, A. & Macgregor, S. VEGAS2: Software for more flexible gene-based testing. Twin Res Hum Genet. 18, 86-91, doi: 10.1017/thg.2014.79 (2015).
# Convert to vcf using PLINK command:  plink --bfile /path/to/yourfile --recode vcf --out /path/to/yourfile
# Reference for PLINK: Chang, C. C. et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience 4, 7, doi:10.1186/s13742-015-0047-8 (2015).

tar -xzf g1000p3_EUR.tar.gz
plink --bfile g1000p3_EUR --recode vcf --out g1000p3_EUR.vcf

# Download CKDGen data files for EUR, eGFR and BUN:
# Reference: Wuttke, M. et al. A catalog of genetic loci associated with kidney function from analyses of a million individuals. Nat Genet 51, 957-972, doi:10.1038/s41588-019-0407-x (2019).
wget https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz
wget https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt.gz
# Get SNPs for eGFR and BUN, as above -> data.txt (for example)

# Download UKHLS data files for eGFR and BUN:
# Reference for UKHLS GWAS: Prins, B. P. et al. Genome-wide analysis of health-related biomarkers in the UK Household Longitudinal Study reveals novel associations. Sci Rep 7, 11008, doi:10.1038/s41598-017-10812-1 (2017).
# Reference for EBI database: Buniello, A. et al. The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic Acids Res 47, D1005-D1012, doi:10.1093/nar/gky1120 (2019).
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005063/harmonised/28887542-GCST005063-EFO_0005208-build37.f.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005070/harmonised/28887542-GCST005070-EFO_0004741-build37.f.tsv.gz
# Get SNPs for eGFR and BUN, as above -> data.txt (for example)

# Use BCFTools to extract only SNPs needed in dataset from 1000GP:
awk '{print $1}' data.txt > data_RSID.txt
bcftools view --include ID==@data_RSID.txt g1000p3_EUR.vcf > g1000p3_EUR_data.vcf

# Count how many SNPs:
grep -v '#' g1000p3_EUR_data.vcf | wc -l
# e.g. 8,252,696

# Count how many columns:
awk '{print NF}' g1000p3_EUR_data.vcf | sort -nu | tail -n 1
# e.g. 512 cols.

# For each SNP, check that the 1000GP reference and alternative alleles match those in the data.
# If not, check whether switching the 1000GP reference and alternative alleles round then matches,
#  -> if so, switch them, and also switch the GT data to match.

# Firstly, get all 1000GP matching SNPs in dataset:
awk 'NR==FNR{a[$1]=" "$0; next}{ print $0 (a[$3]?a[$3]:" missing")}' OFS='\t' \
data.txt g1000p3_EUR_data.vcf > data_in_1000GP_EU.vcf

# Count how many SNPs:
wc -l data_in_1000GP_EU.vcf
# e.g. 8,252,727



# === THE FOLLOWING STEPS SHOWN ARE SPECIFIC FOR THE BIOBANK JAPAN DATASET, BUT WERE APPLIED TO ALL DATASETS.
# Filter to only include common vars MAF >=1% based on gnomAD genomes (European for CKDGen, UKHLS, and East Asian for BioBank Japan):
#  A. Make input file for VEP: #CHROM POS ID REF ALT QUAL FILTER INFO .vcf
awk '{print $1, $2, $3, $4, $5, "NA", "NA", "NA"}' \
bbj-a-60.eGFR.bbj-a-11.BUN.noheader.vcf \
| grep -v '#' | grep -v 'X' | grep -v 'Y'\
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.noXY.forVEP.vcf

cat header.vcf bbj-a-60.eGFR.bbj-a-11.BUN.noheader.noXY.forVEP.vcf \
> bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.vcf

# Get SNPs with gnomad genomes EUR MAF >= 1% only:
# activate conda environment for VEP
conda activate VEP_env

vep -i \
bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.vcf \
--pick_allele \
--stats_file stats_variant_effect_output_pickallele.txt \
--output_file variant_effect_output_pickallele.txt \
--offline \
--cache \
--custom ~/.vep/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_EAS

# Match back to original file -> if alt. allele in gnomAD matches data alt. (coded/effect) allele -> then use.
awk 'NR==FNR{a[$1":"$2":"$5]=" "$0; next}{ print $0 (a[$2":"$3]?a[$2":"$3]:" missing")}' \
OFS='\t' \
bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.vcf \
variant_effect_output_pickallele.txt \
| grep -v missing \
> bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.txt

# Edit file to keep only data needed for filtering step (the MAF):
awk '{print $1,$2,$3,$4,$5,$6,$7,$14}' \
bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.txt \
| grep "gnomADg_AF_EAS=" \
> bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.txt

sed -E 's/IMPACT.+gnomADg_AF_EAS/gnomADg_AF_EAS/g' \
bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.txt \
| sed -E 's/;gnomADg_FILTER=PASS//g' \
| awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' \
| sed 's/gnomADg_AF_EAS=//g' \
> bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.reformatted.txt

# Create list of SNP RSIDs with gnomad EAS genomes MAF >=0.01 only
awk '{if($8>=0.01) print $1}' \
bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.reformatted.txt \
> bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.reformatted.RSID.txt

# Get 1000GP SNPs in MAF>=1% SNPs (only based on RSID - so ref and alt may not match - this is addressed later.)
plink2 --bfile ~/reference_1000GP_EAS/g1000p3_EAS \
	--extract bbj-a-60.eGFR.bbj-a-11.BUN.noXY.forVEP.allele_matched.gnomADg_AF_EAS.reformatted.RSID.txt \
	--keep-allele-order \
	--double-id \
	--make-bed \
	--out ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.gnomADg_AF_EAS.0.01

# 6. Recode as VCF:
plink2 --bfile ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.gnomADg_AF_EAS.0.01 \
	--keep-allele-order \
	--double-id \
	--recode vcf \
	--out ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.gnomADg_AF_EAS.0.01.vcf

# ==== For each SNP matched on RSID between 1000GP-data, keep only those that have matching alt. and ref. alleles
# ==== For each step, the first command is for EUR and the second is for EAS.
# If the 1000GP ref allele == dataset alt. allele, then switch 1000GP GTs, as we want the GTs for the dataset alt. allele.
awk '{if ($4==toupper($514)) { gsub("0\\/0","2"); gsub("[12345]\\/[12345]","0"); gsub("0\\/[12345]","1"); gsub("[12345]\\/0","1"); print } }' \
data_in_1000GP_EU.vcf \
> data_in_1000GP_EU_refaltfixed.vcf

awk '{if ($4==toupper($518)) { gsub("0\\/0","2"); gsub("[12345]\\/[12345]","0"); gsub("0\\/[12345]","1"); gsub("[12345]\\/0","1"); print } }' \
reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01.vcf \
> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_refaltfixed.vcf

# If the 1000GP alt allele first letter (first character) == dataset alt. allele, then use GT=1 (first allele) 1000GP data.
# e.g. if got 2x GT1 (1/1) -> this equals two alt alleles -> substitute for a 2.
awk '{if (substr($5,1,1)==toupper($514)) {gsub("[02345]\\/[02345]","0"); gsub("1\\/1","2"); gsub("[02345]\\/1","1"); gsub("1\\/[02345]","1"); print } }' \
data_in_1000GP_EU.vcf \
> data_in_1000GP_EU_refaltOK_1.txt

awk '{if (substr($5,1,1)==toupper($518)) {gsub("[02345]\\/[02345]","0"); gsub("1\\/1","2"); gsub("[02345]\\/1","1"); gsub("1\\/[02345]","1"); print } }' \
reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01.vcf \
> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_refaltOK_1.txt


# If the 1000GP alt allele second letter (third character) == dataset alt. allele, then use GT=2 (second allele) 1000GP data.
awk '{if (substr($5,3,1)==toupper($514)) {gsub("[01345]\\/[01345]","0"); gsub("2\\/2","2"); gsub("2\\/[01345]","1"); gsub("[01345]\\/2","1"); print } }' \
data_in_1000GP_EU.vcf \
> data_in_1000GP_EU_refaltOK_2.txt

awk '{if (substr($5,3,1)==toupper($518)) {gsub("[01345]\\/[01345]","0"); gsub("2\\/2","2"); gsub("2\\/[01345]","1"); gsub("[01345]\\/2","1"); print } }' \
~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01.vcf \
> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_refaltOK_2.txt


# If the 1000GP alt allele third letter (fifth character) == dataset alt. allele, then use GT=3 (thid allele) 1000GP data.
awk '{if (substr($5,5,1)==toupper($514)) {gsub("[01245]\\/[01245]","0"); gsub("3\\/3","2"); gsub("3\\/[01245]","1"); gsub("[01245]\\/3","1"); print } }' \
data_in_1000GP_EU.vcf \
> data_in_1000GP_EU_refaltOK_3.txt

awk '{if (substr($5,5,1)==toupper($518)) {gsub("[01245]\\/[01245]","0"); gsub("3\\/3","2"); gsub("3\\/[01245]","1"); gsub("[01245]\\/3","1"); print } }' \
~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01.vcf \
> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_refaltOK_3.txt


# Concatenate all four files (or only include those which aren't empty)
# Get rid of header first:
tail -4565703 data_in_1000GP_EU_refaltOK_1.txt > data_in_1000GP_EU_refaltOK_1_nohead.txt
cat data_in_1000GP_EU_refaltfixed.vcf data_in_1000GP_EU_refaltOK_1_nohead.txt >> data_in_1000GP_EU_GTs_sorted.txt
# Check how many lines have a number 0-5/ a number 0-5 (i.e. have genotype data as expected):
grep "[0-5]/[0-5]" data_in_1000GP_EU_GTs_sorted.txt | wc -l
# e.g. 8,245,794

# Check for lost SNPs
awk 'NR==FNR{a[$3]=" "$0; next}{ print $0 (a[$3]?a[$3]:" absent")}' OFS='\t' \
data_in_1000GP_EU_GTs_sorted.txt data_in_1000GP_EU.vcf | grep absent > missing_snps.txt
# e.g. 18 1000GP snps (in missing_snps.txt) that do not have an allele that matches the effect/alternative allele for our data.

# Create version without header, to use as input for PLINK (for SNP pruning)
grep -v '#' data_in_1000GP_EU_GTs_sorted.txt > data_in_1000GP_EU_GTs_sorted_nohead.txt


# ======= Apply LD-based SNP pruning using PLINK and an r cut-off of 0.2
# 1. Split by chromosome (1 per file) and run on different servers (to speed up):
CHROMLIST=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
for j in "${CHROMLIST[@]}"
do
	echo $j
	awk -v j="$j" '{if ($1==j) print $0}' \
	data_in_1000GP_EU_GTs_sorted.txt \
	> data_in_1000GP_EU_GTs_sorted_chr"$j".txt &
done

for j in "${CHROMLIST[@]}"
do
	echo $j
	awk -v j="$j" '{if ($1==j) print $0}' \
	~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted.txt \
	> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j".txt &
done

# Save header
head -28 data_in_1000GP_EU_GTs_sorted.txt \
 > data_in_1000GP_EU_GTs_sorted.header.txt

head -28 ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.gnomADg_AF_EAS.0.01.vcf.vcf \
> ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.head.vcf

CHROMLIST1=("1" "2" "3")
CHROMLIST2=("6" "7" "8" "9" "10")
CHROMLIST3=("11" "12" "13" "14" "15")
CHROMLIST4=("16" "17" "18" "19" "20" "21" "22")
CHROMLIST5=("4" "5")

cd ~/reference_1000GP_EAS
mkdir per_chrom

# CHROMLIST1 on server1, CHROMLIST2 on server2, etc...
for j in "${CHROMLIST1[@]}"
do
	echo $j
	cat ~/reference_1000GP_EAS/g1000p3_EAS_in_bbj.head.vcf \
	~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j".txt \
	> ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j"_head.vcf
	plink --vcf ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j"_head.vcf \
	--keep-allele-order \
	--double-id \
	--make-bed \
	--out ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j"_head_makebed
done

for j in "${CHROMLIST1[@]}"
do
	echo $j
	plink --bfile ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j"_head_makebed \
	--indep-pairwise 148000 5 0.2 \
	--out ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr"$j"_head_makebed_148k_LD02
done

# Check got all prune.in files for chrs 1 - 22:
cd ~/reference_1000GP_EAS/per_chrom
ls -lt | grep "prune.in" | sort

# Concatenate all prune.in:
cat ~/reference_1000GP_EAS/per_chrom/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chr*_head_makebed_148k_LD02.prune.in \
>> ~/reference_1000GP_EAS/g1000p3_EAS_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chrALL_ld02.vcf.prune.in
# e.g. 196,253 SNPs

# Create reference genotype file (for metaCCA XX matrix) keeping pruned SNPs only:
grep -wFf \
~/reference_1000GP_EAS/g1000p3_EAS_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chrALL_ld02.vcf.prune.in \
~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted.txt \
> ~/reference_1000GP_EAS/g1000p3_EAS_WITH_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_ld02.txt

# Create genotype-phenotype file (for metaCCA XY matrix) keeping pruned SNPs only:
grep -wFf \
~/reference_1000GP_EAS/g1000p3_EAS_bbj.gnomADg_AF_EAS.0.01_GTs_sorted_chrALL_ld02.vcf.prune.in \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.vcf \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02.txt
# e.g. 196,253 SNPs


# ==== Annotate SNPs with gene names
# Download glist-hg19 from PLINK website:
cd ~
mkdir hg19_files
cd hg19_files
wget https://www.cog-genomics.org/static/bin/plink/glist-hg19

# Order hg19 gene list file by cols 1 (chr), 2 (pos) first:
sort -k1,1 -k2,2n ~/hg19_files/glist-hg19 > ~/hg19_files/glist-hg19_sort

# For each gene range, for each SNP, if SNP pos is within this gene pos. range, then add gene name to SNP line.
# (ensure autocomplete setting is off)
awk '
	NR==FNR{ range[$1,$2,$3,$4]; next }
	FNR==1
	{
		for(x in range) {
			split(x, check, SUBSEP);
			if($1==check[1] && $2>=check[2] && $2<=check[3]) print $0, check[1], check[2], check[3], check[4]
		}
	}' ~/hg19_files/glist-hg19_sort bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02.txt \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_hg19genes.txt


# Count how many unique genes
awk '{print $13}' bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_hg19genes.txt \
| sort | uniq | wc -l
# e.g. 12,827 genes

# Filter out SNPs that have eGFR and BUN effects in same direction (as only interested in kidney function relevant SNPs where
# eGFR and BUN in opposite directions - see Wuttke et al., 2019 paper)
awk '{if(($6<0 && $8>0) || ($6>0 && $8<0)) print $0}' \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02.txt \
| wc -l
# e.g. 121,163
# Save new file with above filter.
awk '{if(($6<0 && $8>0) || ($6>0 && $8<0)) print $0}' \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02.txt \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt.txt

# How many SNPs with hg19 genes?:
awk 'NR==FNR{a[$3":"$4":"$5]=" "$0; next}{ print $0 (a[$3":"$4":"$5]?a[$3":"$4":"$5]:" absent")}' OFS='\t' \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt.txt \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_hg19genes.txt \
| grep -v absent | wc -l

# Save with hg19 gene names - only one gene per SNP allowed:
awk 'NR==FNR{a[$3":"$4":"$5]=" "$0; next}{ print $0 (a[$3":"$4":"$5]?a[$3":"$4":"$5]:" absent")}' OFS='\t' \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_hg19genes.txt \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt.txt \
| grep -v absent | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$22}' \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes.txt
# Count how many genes:
awk '{print $10}' bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes.txt \
| sort | uniq | wc -l

# Save with hg19 gene names - multiple genes per SNP allowed:
awk 'NR==FNR{a[$3":"$4":"$5]=" "$0; next}{ print $0 (a[$3":"$4":"$5]?a[$3":"$4":"$5]:" absent")}' OFS='\t' \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt.txt \
~/bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_hg19genes.txt \
| grep -v absent | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13}' \
> bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes_multi.txt

# Count how many genes:
awk '{print $10}' bbj-a-60.eGFR.bbj-a-11.BUN.noheader.gnomADg_AF_EAS.0.01_GTs_sorted_ld02_dirfilt_hg19genes_multi.txt \
| sort | uniq | wc -l
