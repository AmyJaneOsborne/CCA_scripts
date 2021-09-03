# This script was written by Amy J. Osborne, amy.osborne@bristol.ac.uk.
# This script can be used to apply plink-based filtering for the quality control (QC) of genotype data
# and is based on the following sources:
# UK Biobank data analyses by MRC IEU:
#  https://data.bris.ac.uk/datasets/pnoat8cxo0u52p6ynfaekeigi/MRC%20IEU%20UK%20Biobank%20GWAS%20pipeline%20version%202.pdf
#  https://data.bris.ac.uk/datasets/1ovaau5sxunp2cv8rcy88688v/UK%20Biobank%20Genetic%20Data_MRC%20IEU%20Quality%20Control%20version%202.pdf
# Data management and summary statistics using PLINK:
# https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_3

# Using conda for plink, plink2

# define directories:
VCFdir=~/GWAS/files
mkdir -p ~/GWAS/plink_files
plinkdir=~/GWAS/plink_files
# QC variant excl. files (excluded based on GenTrain Score and Duplicate_probe - in VCF files)
mkdir -p ~/GWAS/var_excl
varexcldir=~/GWAS/var_excl
# QC sample excl. files (withdrawn samples, etc.)
mkdir -p ~/GWAS/sample_excl_files
samplexcldir=~/GWAS/sample_excl_files

mkdir -p ~/GWAS/files/count_stats
statsdir=~/GWAS/files/count_stats

# Calculate beginning stats:
zcat $VCFdir/Genotyping_file.vcf.gz | bcftools +counts > $statsdir/Genotyping_file.count.stats
# Outputs numbers of samples, SNPs, INDELs, MNPs, others, sites.

# 1. Filter out variants with GEN TRAIN SCORE < 0.75 using BCFtools (because RSIDs nor CHROM-POS-REF-ALT are not unique)
bcftools query -i 'FILTER="GEN_TRAIN_SCORE"' -f '%ID\n' $VCFdir/Genotyping_file.vcf.gz | wc -l
# Outputs xx variants with GEN TRAIN SCORE < 0.75 to be excluded.

bcftools view -e 'FILTER="GEN_TRAIN_SCORE"' $VCFdir/Genotyping_file.vcf.gz -Oz -o $varexcldir/Genotyping_file.GTS.vcf.gz
zcat $varexcldir/Genotyping_file.GTS.vcf.gz | bcftools +counts > $statsdir/Genotyping_file.GTS.count.stats.testing
# Outputs numbers of samples, SNPs, INDELs, MNPs, others, sites.

# 2. Filter out variants with DUPLICATED_PROBE by using BCFtools (because RSIDs nor CHROM-POS-REF-ALT are not unique)
zcat $varexcldir/Genotyping_file.GTS.vcf.gz | grep -v "^#" | grep "DUPLICATED_PROBE" | wc -l
# Outputs xx variants (of those remaining after step 1) with DUPLICATED PROBE, to exclude.

bcftools view -e 'FILTER="DUPLICATED_PROBE"' $varexcldir/Genotyping_file.GTS.vcf.gz -Oz -o $varexcldir/Genotyping_file.GTS.DUP.vcf.gz
zcat $varexcldir/Genotyping_file.GTS.DUP.vcf.gz | bcftools +counts > $statsdir/Genotyping_file.GTS.DUP.count.stats

# 3. Convert VCF to plink format, keep allele order.
plink --vcf $varexcldir/Genotyping_file.GTS.DUP.vcf.gz --keep-allele-order --double-id --make-bed --out $plinkdir/Genotyping_file.GTS.DUP

# 4. List duplicate variants: save output to plink.dupvar
plink --bfile $plinkdir/Genotyping_file.GTS.DUP --list-duplicate-vars --out Genotyping_file.GTS.DUP.dupvar
wc -l Genotyping_file.GTS.DUP.dupvar.dupvar
# Do rest of QC steps first, then run this step again to exclude any duplicates (to try and retain the non-filtered out sites)

# 5. Exclude samples and variants from manual sample QC check:
plink --bfile $plinkdir/Genotyping_file.GTS.DUP \
	--set-missing-var-ids @:# \
	--remove $samplexcldir/allexclusions_xxsamples.txt \
	--keep-allele-order \
	--double-id \
	--make-bed \
  --out $plinkdir/Genotyping_file.GTS.DUP.excl

# 6. Exclude variants with MAF<1%, missing GT call rate >=1.5%
plink --bfile $plinkdir/Genotyping_file.GTS.DUP.excl \
	--maf 0.01 \
	--geno 0.015 \
	--keep-allele-order \
	--double-id \
	--make-bed \
  --out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno
# xx variants removed due to missing genotype data (--geno).
# xx variants removed due to minor allele threshold(s)
# xx variants and xx people pass filters and QC.

# 7. Optional - create statistical report on missing rates:
#plink --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno --missing --out $plinkdir/CKD_GTS.DUP.excl_224.excl.maf.genomissing-stat

# 8. Remove samples with >10% low call rate SNPs (should not be any or very few, as already filtered out low-call rate SNPs):
plink --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno --mind 0.1 --keep-allele-order --double-id --make-bed --out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind

# 9. Using plink2: Hardy-Weinberg test - remove variants that violate assumptions with p-value < 0.001%:
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind --hwe 0.0001 keep-fewhet --keep-allele-order --double-id --make-bed --out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2
# --hwe: xx variants removed due to Hardy-Weinberg exact test

# 10. Using plink2: LD-based SNP pruning using default parameter settings
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2 --indep-pairwise 50 5 0.2 --out $plinkdir/CKD
# Pruning complete.  xx/xx variants removed.
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2 \
	--exclude $plinkdir/CKD.prune.out \
	--keep-allele-order \
	--double-id \
	--make-bed \
	--out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr
# xx variants remaining.

# 11. KING to remove first degree relatives (and closer relations):
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr \
	--king-cutoff 0.177 \
	--out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr
# Excluding xx variants on non-autosomes from KING-robust calculation.
# xx variants handled by initial scan (xx remaining).

plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr \
	--remove $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king.cutoff.out.id \
	--keep-allele-order \
	--double-id \
	--make-bed \
	--out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king

# was 7. Optional - create statistical report on missing rates:
plink --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr --missing --out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.missing-stat
# Total genotyping rate is xx.

# 12. PCA to look for any clustering which suggests population structure
# (may be some due to ancestry which hasn't been looked at yet)
# Then open output in R and plot with batch no. to check for batch effects (see R script: checking_missing_batch_eff.R)
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king \
     --pca 5 \
     --out $plinkdir/CKD_pca_results_king
# Excluding xx variants on non-autosomes from PCA.
# -> e.g. of output: no batch effects but can see some ancestry effects

# 13. Test for duplicates:
plink --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr --list-duplicate-vars --out $plinkdir/CKD

# 14. Optional - exclude samples and recode as VCF:
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king \
	--remove xxsamples_toexcl.txt \
	--keep-allele-order \
	--double-id \
	--recode vcf \
	--out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king.excl
# xx samples remaining

# Phasing - not required for CCA because only interested in allele counts and not whether >2 alleles are on the same chromosome.

# INCREASE LD SNP PRUNING (SAME AS FOR METACCA ANALYSES:)
# 15. Using plink2: LD-based SNP pruning using default parameter settings
plink2 --bfile $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king \
	--indep-pairwise 148000 5 0.2 \
	--keep-allele-order \
	--double-id \
	--recode vcf \
	--out $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king.prunemax
# xx/xx variants removed.

# Further preparation for CCA:
# Remove header
grep -v '##' $plinkdir/Genotyping_file.GTS.DUP.excl.maf.geno.mind.hwe2.ldpr.king.prunemax.vcf > out.prunemax.vcf
# Exclude MT and X chr variants:
grep -v 'MT' out.prunemax.vcf > out.prunemax.noMT.vcf
grep -v 'X' out.prunemax.noMT.vcf > out.prunemax.noMT.norX.vcf
# Save header line:
head -1 out.prunemax.noMT.norX.vcf > out.prunemax.noMT.norX_header.vcf
# For each SNP, get count of alt. allele only:
awk '{gsub("0\\/0","0"); gsub("1\\/0","1"); gsub("0\\/1","1"); gsub("1\\/1","2"); print }' \
out.prunemax.noMT.norX.vcf > out.prunemax.noMT.norX_alt_counts.vcf

cat out.prunemax.noMT.norX_header.vcf out.prunemax.noMT.norX_alt_counts.vcf \
> out.prunemax.noMT.norX_alt_counts_head.vcf

cp out.prunemax.noMT.norX_alt_counts_head.vcf > GT_data_file.txt

# How many variants have missing samples?
grep -v ".\\/." out.prunemax.noMT.norX_alt_counts_head.vcf | sed 's/#//g' \
> out.prunemax.noMT.norX_alt_counts_head_nomissing.txt
