#!/usr/bin/env Rscript

# load in covariate file
data <- read.table("RBD_covs_header.tab", header = T) 
data$fatid <- 0
data$matid <- 0
data$SEX_COV <- data$sex-1
data$PHENO <- data$PHENO_PLINK-1
colnames(data)[1] <- "fid"
colnames(data)[2] <- "iid"
attach(data)
dat <- data[,c("fid","iid","fatid","matid","sex","SEX_COV","AGE","PHENO","PHENO_PLINK","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")]
write.table(dat,file=paste("sampleInfo_RBD.tab"), quote = F, sep = "\t", row.names = F)
write.table(dat$iid,"samplesToKeep_RBD.tab", quote = F, sep = "\t", row.names = F, col.names = F)

q()
n