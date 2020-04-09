#!/bin/bash

cohort1=$1
cohort2=$2

for cohort in $cohort1 $cohort2;
    # make a column for chr:bp instead of SNP names.
do awk '{print $1,$1,$4,$3,$4,$5,$6}' ${cohort}.bim > andmore.bim
    sed -e 's/\s/':'/2' andmore.bim > ${cohort}.bim

    ##### REMOVE palindromes.

    awk '{print $2"\t"$5"-"$6}' ${cohort}.bim > check_palindrome.txt

    awk '{if ($2=="A-T") print $1}' check_palindrome.txt > AT.txt
    awk '{if ($2=="T-A") print $1}' check_palindrome.txt > TA.txt
    awk '{if ($2=="G-C") print $1}' check_palindrome.txt > GC.txt
    awk '{if ($2=="C-G") print $1}' check_palindrome.txt > CG.txt

    cat AT.txt TA.txt GC.txt CG.txt > exclude_palindromes.txt

    plink --bfile $cohort --exclude exclude_palindromes.txt --make-bed --out ${cohort}_pal

    rm AT.txt
    rm TA.txt
    rm GC.txt
    rm CG.txt
    rm check_palindrome.txt
    mv exclude_palindromes.txt ${cohort}_palindromes.txt

    # REMOVE DUPLICATES

    echo "DUP" > dup.txt
  
    awk 'BEGIN{FS=OFS="\t";}{if(snp[$2]>0){$2="DUP";} else{snp[$2]=1} print $0}' ${cohort}_pal.bim > uniq.bim
    
    mv uniq.bim ${cohort}_pal.bim
    
    plink --bfile ${cohort}_pal --exclude dup.txt --make-bed --out unique_${cohort}
    
    rm dup.txt

    awk '{print $2}' unique_${cohort}.bim > ${cohort}_snps_to_compare.txt

done

comm -12 <(sort ${cohort1}_snps_to_compare.txt) <(sort ${cohort2}_snps_to_compare.txt) > comm_vars.txt

plink --bfile ${cohort1}_pal --extract comm_vars.txt --make-bed --out ${cohort1}_to_merge
plink --bfile ${cohort2}_pal --extract comm_vars.txt --make-bed --out ${cohort2}_to_merge
plink --bfile ${cohort1}_to_merge --bmerge ${cohort2}_to_merge --out MERGED_${cohort1}_${cohort2}

# may need to flip-snps. pay attention to output message
#plink --bfile ${cohort1}_to_merge --flip MERGED_${cohort1}_${cohort2}.missnp --make-bed --out flipped
#plink --bfile flipped --bmerge ${cohort2}_to_merge --out MERGED_${cohort1}_${cohort2}

### if it persists, and there are only a few variants with issues, I remove them.
#plink --bfile flipped --exclude MERGED_${cohort1}_${cohort2}.missnp --make-bed --out more_flipped
#plink --bfile ${cohort2}_to_merge --exclude MERGED_${cohort1}_${cohort2}.missnp --make-bed --out flipped2
#plink --bfile more_flipped --bmerge flipped2 --out MERGED_${cohort1}_${cohort2}
