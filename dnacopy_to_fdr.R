# take varscan copy number data to FDR values
# this is a work in progress
# this is not working as of 2014-5-25
  
library(DNAcopy)
  # loads DNAcopy R program
cn <- read.table("your.cn.file",header=F)
  
CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], 
                    data.type = 'logratio')
  # CNA.object creates a 'copy number array' data object
CNA.smoothed <- smooth.CNA(CNA.object)
  #smooths outliers reads 
segs <- segment(CNA.smoothed, verbose=1, min.width=2)
  
segs2 = segs$output
segs3 <- segments.p(segs)
# write pvalues to txt file 
# check setwd()
write.table(segs3[,1:8], file='seg_pval.txt',row.names=F,col.names=T,quote=F,sep='\t')
  
# how to read txt file and set it to variable pval
pval <- read.table('seg_pval.txt',header=T,sep='\t')
# check data visually
pval
# pvalues are kept in column 8 'pval' with NA for uncalcualted (I think these equal 1)
  
  # make vector of p values column
pvalues <- pval$pval
  # check data
pvalues
  
  # for now - remove 'NA' values
pvalues2 <- na.omit(pvalues)
  # make pvalues 2 a vector
pVal <- as.vector(pvalues2)
  
  # load fdrtool suite for FDR calculation
library(fdrtool)
FDR <- fdrtool(pVal,statistic='pvalue')
  
  # tail area-based FDR (p value corrected for multiplicity)
  # FDR of q < 0.1 is significant
FDR$qval 
  
  # calculate negative log10 qvalue 
qval.neglog10 <- -log10(FDR$qval)
  
plot(qval.neglog10)
  # the higher the -log10 value = more significant  