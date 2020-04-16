#!/usr/bin/env Rscript

#PARSE INTO FILES
# run in the same folder as imputation results. 

R

library(plyr)

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  da <- dat[,c("SNP","ALT_Frq","Rsq")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".info",sep = ""), row.names = F, quote = F, sep = "\t")
}

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
  da <- dat[,c("range")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}

q()
n