

vcftools --gzvcf ../all_genes.all_samples.MIP76.vcf.gz --keep keep_RBD.txt --chr 4 --from-bp 926062 --to-bp 952643 --recode --recode-INFO-all --minDP 30 --out tmem175.rbd.30x
