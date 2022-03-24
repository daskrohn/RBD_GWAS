#!/bin/bash

FILENAME=$1
THREADS=$2

### Actual cleaning pre-imputation code 
### start like this sh STEP1_pre_imputation_step.sh FILENAME
## part 1 of 2 
## part to is written in STEP2_pre_imputation_step.sh
## Cornelis
## Last update August 3 2018
## Run using plink 1.9 and on HG19

### Sample cleaning inclusion:

# NOTE did you already filter for gentrain scores??
# if not perhaps you want to do that....

# NOTE do all samples have a gender??
# if not those samples that do not have a gender will be removed

# NOTE do all your samples have an affection status??
# if not this will cause trouble at the end of this script

## Genotyping call rates per sample
# open file call_rates.imiss and F_MISS - 1 = callrate
# all call rates saved in CALL_RATES_ALL_SAMPLES.txt

module load gcta
module spider gcc/7.3.0
module spider r/3.5.2
module spider r-bundle-bioconductor/3.9


plink --bfile $FILENAME --missing --out call_rates

rm call_rates.lmiss
rm call_rates.log
rm call_rates.hh
rm call_rates.nosex
mv call_rates.imiss log/CALL_RATES_ALL_SAMPLES.txt

## No heterozygosity outliers  
# --het from LD pruned data > use F cut-off of -0.15 and <- 0.15 for inclusion
# outliers stored here -> all_outliers.txt
# all heterozygosity is stored here -> HETEROZYGOSITY_DATA.txt

plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet

awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt

cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt

mv prunedHet.het HETEROZYGOSITY_DATA.txt

plink --bfile $FILENAME --remove all_outliers.txt --make-bed --out after_heterozyg

mv HETEROZYGOSITY_DATA.txt log

rm pruning.hh
rm pruning.nosex
rm pruned_data.bed
rm pruned_data.bim
rm pruned_data.fam
rm pruned_data.hh
rm pruned_data.log
rm pruned_data.nosex
rm prunedHet.hh
rm prunedHet.log
rm prunedHet.nosex
rm pruning.log
rm pruning.prune.in
rm pruning.prune.out
rm outliers1.txt
rm outliers2.txt
rm all_outliers.txt

## No call rate outliers for samples
# all call rates outliers are in CALL_RATE_OUTLIERS.txt

plink --bfile after_heterozyg --mind 0.05 --make-bed --out after_heterozyg_call_rate

mv after_heterozyg_call_rate.irem log/CALL_RATE_OUTLIERS.txt

rm after_heterozyg.bed
rm after_heterozyg.bim
rm after_heterozyg.fam
rm after_heterozyg.hh
rm after_heterozyg.log
rm after_heterozyg.nosex
rm after_heterozyg_call_rate.hh
rm after_heterozyg_call_rate.nosex

## No sex check fails -> check with original sample sheet 
# (Note for neuroX data that does not have GWAS back bone, use F cut-off of 0.50 instead of 0.25/0.75 and only use the PAR region’s common variants always)
# PAR = --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 
# gender failures are stored in GENDER_FAILURES.txt
# gender checks are stored in GENDER_CHECK1.txt and GENDER_CHECK2.txt

plink --bfile after_heterozyg_call_rate --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
plink --bfile after_heterozyg_call_rate --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 

touch samples_to_remove.txt
grep "PROBLEM" gender_check1.sexcheck > problems1.txt
grep "PROBLEM" gender_check2.sexcheck > problems2.txt
cat problems1.txt problems2.txt > GENDER_FAILURES.txt

cut -f 1,2 GENDER_FAILURES.txt > samples_to_remove.txt

plink --bfile after_heterozyg_call_rate --remove samples_to_remove.txt --make-bed --out after_gender

rm gender_check2.hh
rm gender_check2.log
rm gender_check2.nosex
rm gender_check1.hh
rm gender_check1.log
rm gender_check1.nosex
rm after_heterozyg_call_rate.bed
rm after_heterozyg_call_rate.bim
rm after_heterozyg_call_rate.fam
rm after_heterozyg_call_rate.log

mv gender_check1.sexcheck GENDER_CHECK1.txt
mv gender_check2.sexcheck GENDER_CHECK2.txt

mv GENDER_FAILURES.txt log

## No ancestry outliers -> based on Hapmap3 PCA plot, should be near combined CEU/TSI

# downloaded hapmap3 plink format see other script -> hapmap_prep.sh
# NOTE use other script when using NeuroX data
# neuroX_snps_for_hapmap.txt # extract these snps from neuroX
# neuroX_snps_for_hapmap_conversion.txt # conversion of IDs --update-map
# Keep in mind that this comparison with hapmap is based on the number of SNPs that overlap between your input dataset and hapmap

plink --bfile after_gender --bmerge HAPMAP_hg19_new_pos --out hapmap3_bin_snplis --make-bed

plink --bfile after_gender --flip hapmap3_bin_snplis-merge.missnp --make-bed --out after_gender3




plink --bfile after_gender3 --bmerge HAPMAP_hg19_new_pos --out hapmap3_bin_snplis --make-bed

plink --bfile after_gender3 --exclude hapmap3_bin_snplis-merge.missnp --out after_gender4 --make-bed

plink --bfile after_gender4 --bmerge HAPMAP_hg19_new_pos --out hapmap3_bin_snplis --make-bed

plink --bfile hapmap3_bin_snplis --geno 0.05 --out pca --make-bed --pca --threads $THREADS


grep "EUROPE" pca.eigenvec > eur.txt
grep "ASIA" pca.eigenvec > asia.txt
grep "AFRICA" pca.eigenvec > afri.txt
grep -v -f eur.txt pca.eigenvec | grep -v -f asia.txt | grep -v -f afri.txt > new_samples.txt
cut -d " " -f 3 after_gender.fam > new_samples_add.txt
paste new_samples_add.txt new_samples.txt > new_samples2.txt
paste eur_add.txt eur.txt > euro.txt
paste asia_add.txt asia.txt > asiao.txt
paste afri_add.txt afri.txt > afrio.txt

cat new_samples2.txt euro.txt asiao.txt afrio.txt > pca.eigenvec2

# R script for PCA plotting and filtering
R < PCA_in_R.R --no-save  
plink --bfile after_gender --keep PCA_filtered_europeans.txt --make-bed --out after_gender_heterozyg_hapmap
cat PCA_filtered_asians.txt PCA_filtered_africans.txt PCA_filtered_mixed_race.txt > hapmap_outliers33.txt
mv hapmap_outliers33.txt log


# this creates several plots and lists based on genetic ancestry

## No relatedness closer than cousin -> Pihat threshold of 0.125 
# this is a optional step we usually do it but its upto you
# also can multithread when running in linux or with better version of GCTA --threads ...

gcta64 --bfile after_gender_heterozyg_hapmap --make-grm --out GRM_matrix --autosome --maf 0.05 --thread-num $THREADS
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
awk '{print $1,$2,$6}' after_gender_heterozyg_hapmap.fam > pheno.txt
gcta64 --bfile after_gender_heterozyg_hapmap --keep GRM_matrix_0125.grm.id --make-bed --out after_gender_heterozyg_hapmap_pihat



cut -f 1,2 after_gender_heterozyg_hapmap.fam > log/IDs_before_relatedness_filter.txt
cut -f 1,2 after_gender_heterozyg_hapmap_pihat.fam > log/IDs_after_relatedness_filter.txt

## No missingness per variant geno 0.05
# high call rate (95% for GWAS, 90% for NeuroX on exome background)

plink --bfile after_gender_heterozyg_hapmap_pihat --pheno pheno.txt --make-bed --out after_gender_heterozyg_pihat_mind --geno 0.05
grep "(--geno)" after_gender_heterozyg_pihat_mind.log > log/MISSINGNESS_SNPS.txt


### Variant inclusion criteria:
# missingness by case control P > 1E-4 # needs case control status

plink --bfile after_gender_heterozyg_pihat_mind --test-missing --out missing_snps 
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
plink --bfile after_gender_heterozyg_pihat_mind --exclude missing_snps_1E4.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing1
sort -u missing_snps_1E4.txt > log/VARIANT_TEST_MISSING_SNPS.txt


# missing by haplotype P > 1E-4

plink --bfile after_gender_heterozyg_pihat_mind_missing1 --test-mishap --out missing_hap 
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\
/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt

sort -u missing_haps_1E4_final.txt > log/HAPLOTYPE_TEST_MISSING_SNPS.txt
plink --bfile after_gender_heterozyg_pihat_mind_missing1 --exclude missing_haps_1E4_final.txt --make-bed --out after_gender_heterozyg_pihat_hapmap_mind_missing12

# Optional hardy Weinberg SNP from controls

plink --bfile after_gender_heterozyg_pihat_hapmap_mind_missing12 --filter-controls --hwe 1E-4 --write-snplist
plink --bfile after_gender_heterozyg_pihat_hapmap_mind_missing12 --extract plink.snplist --make-bed --out after_gender_heterozyg_pihat_hapmap_mind_missing123

mv after_gender_heterozyg_pihat_hapmap_mind_missing123.bim data_here/FILTERED.bim
mv after_gender_heterozyg_pihat_hapmap_mind_missing123.bed data_here/FILTERED.bed
mv after_gender_heterozyg_pihat_hapmap_mind_missing123.fam data_here/FILTERED.fam

####### REMOVE A LOT TO CLEAN FOLDER...

rm pheno.txt
rm hapmap3_bin_snplis.bed
rm hapmap3_bin_snplis.bim
rm hapmap3_bin_snplis.fam
rm hapmap3_bin_snplis.log
rm hapmap3_bin_snplis.hh
rm after_gender4.bed
rm after_gender4.bim
rm after_gender4.fam
rm after_gender4.hh
rm after_gender4.log
rm hapmap3_bin_snplis-merge.missnp
rm after_gender3.bed
rm after_gender3.bim
rm after_gender3.fam
rm after_gender3.hh
rm after_gender3.log
rm after_gender.bed
rm after_gender.bim
rm after_gender.fam
rm after_gender.hh
rm after_gender.log
rm after_gender.nosex
rm afri.txt
rm asia.txt
rm eur.txt
rm pca.bed
rm pca.bim
rm pca.fam
rm pca.log
rm pca.eigenval
rm afrio.txt
rm asiao.txt
rm euro.txt
rm new_samples_add.txt
rm new_samples.txt
rm new_samples2.txt
rm after_gender_heterozyg_pihat_mind.hh
rm after_gender_heterozyg_hapmap_pihat.bed
rm after_gender_heterozyg_hapmap_pihat.bim
rm after_gender_heterozyg_hapmap_pihat.fam
rm GRM_matrix_0125.grm.id
rm GRM_matrix.grm.id
rm after_gender_heterozyg_hapmap.bim
rm after_gender_heterozyg_hapmap.log
rm after_gender_heterozyg_pihat_mind_missing1.bed
rm after_gender_heterozyg_pihat_mind_missing1.fam
rm after_gender_heterozyg_pihat_mind_missing1.hh
rm missing_snps_1E4.txt
rm after_gender_heterozyg_pihat_mind.bed
rm after_gender_heterozyg_pihat_mind.bim
rm after_gender_heterozyg_pihat_mind.fam
rm after_gender_heterozyg_pihat_mind.log
rm missing_snps.hh
rm missing_snps.log
rm missing_snps.missing
rm missing_haps_1E4_final.txt
rm missing_haps_1E4.txt
rm after_gender_heterozyg_pihat_mind_missing1.bim
rm after_gender_heterozyg_pihat_mind_missing1.log
rm missing_hap.hh
rm missing_hap.log
rm missing_hap.missing.hap
rm after_gender_heterozyg_hapmap.bed
rm after_gender_heterozyg_hapmap.fam
rm after_gender_heterozyg_hapmap.hh
rm plink.snplist
rm *.hh
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.bed
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.bim
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.fam
rm after_gender_heterozyg_pihat_hapmap_mind_missing12.log



