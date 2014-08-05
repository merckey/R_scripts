# take varscan copy number data to copynumber and DNAcopy packages
# this is a work in progress - I think this works though - will need to look at varscan ouput file
# so gets some aspects of DNAcopy similar - but segmentation isn't right
# http://varscan.sourceforge.net/copy-number-calling.html#copy-number-output

# copynumber has mm9 ref
# 2014-5-24

#### I think the problem with running tutorials is that hte size of the example data is to small from DNAcopy


setwd("/path/to/data/files")

# can't hvae DNAcopy and copynumber loaded simultaneously for plotGenome function
library(copynumber)
cn <- read.table("out.file",header=F,colClasses=c("character","numeric","numeric","numeric",
                                                  "numeric","numeric")
  # creates cn variable 'cn' for copy number file from varscan - also making sure data is numeric
cn.cp <- cbind( cn[,1], cn[,2], cn[,7])
  # take chromsome position (cn[,1]), read position (cn[,2]), 
  # and log-base-2 ratio of tumor/normal depth corrected for GC
  # or try data.frame(cn[,1],cn[,2],cn[,7])
colnames(cn.cp) <- c('Chrom','Median.bp','SAMPLE_ID')   
  # add names to data colums


cn.wins <- winsorze(data=cn.cp)
# removes outliers 

cn.seg <- pcf(data=cn.wins,gamma=12)
  # fits segmentation curves to data

head(cn.seg)
  # look at data

plotGenome(data=cn.wins,segments=cp.seg,sample=1,cex=3,assembly='mm9')
  # plots whole genome in one figure for one sample with chrom numbers on x-axis
  # do not go above cex = 3
  # can set assembly to mm9



