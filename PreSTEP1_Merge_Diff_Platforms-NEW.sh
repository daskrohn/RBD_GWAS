#!/bin/bash

module load plink

cohort1=$1
cohort2=$2

# make a column for chr:bp instead of SNP names. 
## ONE WAY
awk '{print $1,$1,$4,$3,$4,$5,$6}' ${cohort1}.bim > andmore.bim
sed -i -e 's/\s/':'/2' andmore.bim 

rm ${cohort1}.bim
mv andmore.bim ${cohort1}.bim

awk '{print $1,$1,$4,$3,$4,$5,$6}' ${cohort2}.bim > andmore.bim
sed -i -e 's/\s/':'/2' andmore.bim 

rm ${cohort2}.bim
mv andmore.bim ${cohort2}.bim

## ANOTHER WAY
#cut -f 1,4 ${cohort1}.bim > temp1.txt
#    sed -i -e 's/\t/:/g' temp1.txt
#    cut -f 2 ${cohort1}.bim > temp2.txt
#    paste temp2.txt temp1.txt > temp3.txt
    
#paste temp3.txt ${cohort1}.bim  > rename_${cohort1}.bim



##### REMOVE palindromes. 

awk '{print $2"\t"$5"-"$6}' ${cohort1}.bim > check_palindrome.txt

awk '{if ($2=="A-T") print $1}' check_palindrome.txt > AT.txt
awk '{if ($2=="T-A") print $1}' check_palindrome.txt > TA.txt
awk '{if ($2=="G-C") print $1}' check_palindrome.txt > GC.txt
awk '{if ($2=="C-G") print $1}' check_palindrome.txt > CG.txt

cat AT.txt TA.txt GC.txt CG.txt > exclude_palindromes.txt

plink --bfile $cohort1 --exclude exclude_palindromes.txt --make-bed --out ${cohort1}_pal

rm AT.txt
rm TA.txt
rm GC.txt
rm CG.txt
rm check_palindrome.txt
mv exclude_palindromes.txt ${cohort1}_palindromes.txt

## second set

awk '{print $2"\t"$5"-"$6}' ${cohort2}.bim > check_palindrome.txt

awk '{if ($2=="A-T") print $1}' check_palindrome.txt > AT.txt
awk '{if ($2=="T-A") print $1}' check_palindrome.txt > TA.txt
awk '{if ($2=="G-C") print $1}' check_palindrome.txt > GC.txt
awk '{if ($2=="C-G") print $1}' check_palindrome.txt > CG.txt

cat AT.txt TA.txt GC.txt CG.txt > exclude_palindromes.txt

plink --bfile $cohort2 --exclude exclude_palindromes.txt --make-bed --out ${cohort2}_pal

rm AT.txt
rm TA.txt
rm GC.txt
rm CG.txt
rm check_palindrome.txt
mv exclude_palindromes.txt ${cohort2}_palindromes.txt

### MERGE



# remove duplicates 
plink --bfile ${cohort1}_pal --list-duplicate-vars --out dup_${cohort1} # for reference. 
plink --bfile ${cohort1}_pal --write-snplist --out all_snps1
cat all_snps1.snplist | sort | uniq -d > duplicated_snps1.snplist

plink --bfile ${cohort1}_pal --exclude duplicated_snps1.snplist --make-bed --out unique_${cohort1}

# next set
plink --bfile ${cohort2}_pal --list-duplicate-vars --out dup_${cohort2} # for reference. 
plink --bfile ${cohort2}_pal --write-snplist --out all_snps2
cat all_snps2.snplist | sort | uniq -d > duplicated_snps2.snplist

plink --bfile ${cohort2}_pal --exclude duplicated_snps2.snplist --make-bed --out unique_${cohort2}

# IS THERE A WAY TO KEEP ONLY ONE? Or must remove?

awk '{print $2}' unique_${cohort1}.bim > ${cohort1}_snps_to_compare.txt
awk '{print $2}' unique_${cohort2}.bim > ${cohort2}_snps_to_compare.txt


# compare common SNPs across platforms; extract and merge

comm -12 <(sort ${cohort1}_snps_to_compare.txt) <(sort ${cohort2}_snps_to_compare.txt) > comm_vars.txt

plink --bfile ${cohort1}_pal --extract comm_vars.txt --make-bed --out ${cohort1}_to_merge
plink --bfile ${cohort2}_pal --extract comm_vars.txt --make-bed --out ${cohort2}_to_merge
plink --bfile ${cohort1}_to_merge --bmerge ${cohort2}_to_merge --out MERGED_${cohort1}_${cohort2}

# may need to flip-snps. pay attention to output message
plink --bfile ${cohort1}_to_merge --flip MERGED_${cohort1}_${cohort2}.missnp --make-bed --out flipped
plink --bfile flipped --bmerge ${cohort2}_to_merge --out MERGED_${cohort1}_${cohort2}

### if it persists, and there are only a few variants with issues, I remove them. 

# plink --bfile flipped --exclude MERGED_RBD-FC_OMNI.missnp --make-bed --out more_flipped
# plink --bfile ${cohort2}_to_merge --exclude MERGED_RBD-FC_OMNI.missnp --make-bed --out flipped2
# plink --bfile more_flipped --bmerge flipped2 --out MERGED_${cohort1}_${cohort2}