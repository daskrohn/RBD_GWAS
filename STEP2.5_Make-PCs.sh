#module load flashpca 

plink --bfile FILTERED --maf 0.01 --geno 0.01 --hwe 5e-6 --make-bed --out filter_better -exclude range ~/runs/krohn/krohn/Main_0_Data/exclusion_regions_hg19.txt --memory 230000 --threads 19 
plink --bfile filter_better --indep-pairwise 1000 10 0.02 --out pruned_better_and_more --memory 230000 --threads 19
plink --bfile filter_better --extract pruned_better_and_more.prune.in --make-bed --out input_pca --memory 230000 --threads 19
plink --bfile input_pca --pca 20 --out FILENAME_pre-impute

#flashpca --bfile input_pca --suffix _$line.txt --numthreads 19
# need to learn how to use flashpca, way faster. 