#!/bin/bash

#CHANGE file path and password to suit current job

cd /PATH/

## add your wget codes here from Michigan Imputation Server. 


unzip -P ADD_PASSWORD '*.zip'

rm *.zip

module load plink 

mkdir plink_files

gunzip *.info.gz

mkdir plink_files/soft_calls
mkdir plink_files/hard_calls

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink --vcf chr${chnum}.dose.vcf.gz --make-bed --out s1 --double-id
  plink --bfile s1 --bmerge s1 --merge-mode 6 --out s8
  plink --bfile s1 --exclude s8.diff --make-bed --out s2
  plink --bfile s2 --list-duplicate-vars --out s4
  plink --bfile s2 --exclude s4.dupvar --make-bed --out s3
  plink --bfile s3 --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.3 --make-bed --out plink_files/soft_calls/chr${chnum}
done


for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink --vcf chr${chnum}.dose.vcf.gz --make-bed --out s1 --double-id
  plink --bfile s1 --bmerge s1 --merge-mode 6 --out s8
  plink --bfile s1 --exclude s8.diff --make-bed --out s2
  plink --bfile s2 --list-duplicate-vars --out s4
  plink --bfile s2 --exclude s4.dupvar --make-bed --out s3
  plink --bfile s3 --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.8 --make-bed --out plink_files/hard_calls/chr${chnum}
done


