### raw output from varscan copynumber should be smoothed and segmented
### by circular binary segmentation (CBS)
### http://varscan.sourceforge.net/copy-number-calling.html#copy-number-segmentation


library(DNAcopy)
# loads DNAcopy R program
cn <- read.table("your.cn.file",header=F)

# defines a variable 'cn' as 'your.cn.file'
CNA.object <- CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], 
                  data.type = 'logratio',sampleid="NAMEHERE")

# CNA.object creates a 'copy number array' data object
CNA.smoothed <- smooth.CNA(CNA.object)

# apply CBS analysis 
segs <- segment(CNA.smoothed, verbose=0, min.width=2)

# visualize plot
plot(segs, plot.type='s')

# and statistical analysis
segs.p.cn <- segments.p(segs)

# A little more background statistical analysis :
# from ?segments.p
"""The p-values are obtained by applying Siegmund's approximation 
for the maximal statistic from binary segmenting consecutive 
segments within a chromosome. This p-value is only to give the 
relative importance of the change-points as the CBS is different 
from the algorithm used here."""

# correct for multiple hypothesis with FDR calculation
fdr <- p.adjust(p, p.adjust.methods, n = length(p))

segs.fdr.cn <- data.frame(segs.p.cn, fdr)

# optional sorting by most significantly altered locations
# segs.sorted.fdr.cn <- segs.fdr.cn[order(segs.fdr.cn$fdr),]

# write individual segments to file 
write.table(segs.fdr.cn, file="out.file", row.names=F, col.names=T,
            quote=F, sep="\t")

# optional writing of segments sorted by FDR significance 
# write.table(segs.sorted.fdr.cn, file='sorted.out.file', row.names=F, quote=F, sep='\t')

# should be able to use this output 'out.file' to move over in to copynumber package, etc