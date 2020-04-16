#!/bin/bash

FILENAME=$1 # name your results

# to run all at once (modify for your system):
# sbatch --cpus-per-task=26 --mem=64g --mail-type=BEGIN,TIME_LIMIT_50,END --time=48:00:00 STEP6_GWAS_launch_per_dataset.sh {FILENAME}

# this includes 6 PCs as covariate. I recommend doing a scree plot to determine the best number of PCs in your data. 


YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr1.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr1_caseControl --single wald --inVcf chr1.dose.vcf.gz  --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr2.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr2_caseControl --single wald --inVcf chr2.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr3.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr3_caseControl --single wald --inVcf chr3.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr4.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr4_caseControl --single wald --inVcf chr4.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr5.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr5_caseControl --single wald --inVcf chr5.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr6.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr6_caseControl --single wald --inVcf chr6.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr7.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr7_caseControl --single wald --inVcf chr7.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr8.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr8_caseControl --single wald --inVcf chr8.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr9.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr9_caseControl --single wald --inVcf chr9.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr10.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr10_caseControl --single wald --inVcf chr10.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 

YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr11.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr11_caseControl --single wald --inVcf chr11.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr12.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr12_caseControl --single wald --inVcf chr12.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr13.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr13_caseControl --single wald --inVcf chr13.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab  
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr14.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr14_caseControl --single wald --inVcf chr14.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab  
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr15.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr15_caseControl --single wald --inVcf chr15.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab  
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr16.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr16_caseControl --single wald --inVcf chr16.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr17.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr17_caseControl --single wald --inVcf chr17.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab  
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr18.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr18_caseControl --single wald --inVcf chr18.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr19.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr19_caseControl --single wald --inVcf chr19.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr20.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr20_caseControl --single wald --inVcf chr20.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 

YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr21.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr21_caseControl --single wald --inVcf chr21.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab  
YOUR_PATH/rvtest --noweb --hide-covar --rangeFile maf001rsq08minimums_chr22.txt --out YOUR_RESULTS_PATH/${FILENAME}.chr22_caseControl --single wald --inVcf chr22.dose.vcf.gz --dosage DS --pheno sampleInfo.tab --pheno-name PHENO --covar sampleInfo.tab --covar-name SEX_COV,AGE,PC1,PC2,PC3,PC4,PC5,PC6 --peopleIncludeFile samplesToKeep.tab 