### combine all info files
cat YOUR_PATH/maf001rsq03minimums_chr*.info | grep -v 'Rsq' > allChrs.Info
### combine all assoc
cat *chr*_caseControl.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs.assoc
 
####Then in R:
R
 
require(data.table)

infos <- fread("allChrs.Info")
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
assoc <- fread("allChrs.assoc")
colnames(assoc) <- c("CHROM","POS","A2","A1","N","Test","Beta","SE","Pvalue")
data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = T)
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
#final columns below
dat$SNP <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat$CHR <- dat$CHROM
dat$BP <- dat$POS
dat$OR <- exp(dat$beta)
dat0 <- dat[,c("CHR","BP","SNP","minorAllele","majorAllele","beta", "OR","se","maf","P")]
write.table(dat0, "GWAS.YOUR_NAME.tab", quote = F, sep = "\t", row.names = F)
q()
n
 
### OPTIONAL
sort -gk 10 GWAS.YOUR_NAME.tab > SORTED_gwas_RBD_Nx.tab