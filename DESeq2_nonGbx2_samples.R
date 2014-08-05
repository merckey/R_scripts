################################
# 2014-6-20
# comparision of Gbx2 High expressing tumors between Ecad + and Ecad - (or noncanonical-EMT and canonical EMT)
#$ grep Gbx2 *.txt
# PD1849E_minus.counts.txt:Gbx2  0
# PD1849E_plus.counts.txt:Gbx2	0
# PD2204E_minus.counts.txt:Gbx2	103
# PD2204E_plus.counts.txt:Gbx2	0
# PD2329E_minus.counts.txt:Gbx2	3
# PD2329E_plus.counts.txt:Gbx2	0
# PD2342E_minus.counts.txt:Gbx2	9
# PD2342E_plus.counts.txt:Gbx2	0
# PD2412E_minus.counts.txt:Gbx2	12
# PD2412E_plus.counts.txt:Gbx2	0
# PD2523E_minus.counts.txt:Gbx2	358
# PD2523E_plus.counts.txt:Gbx2	0
# PD345E_minus.counts.txt:Gbx2	0
# PD345E_plus.counts.txt:Gbx2	0
# PD798E_minus.counts.txt:Gbx2	79
# PD798E_plus.counts.txt:Gbx2	0


# calling PD2523, PD2204, PD798: GBX2 group

library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_nonGBX2_group/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_nonGBX2_group/'

sampleFiles <- grep('counts.txt',list.files(directory),value=T)

# view sampleFiles
sampleFiles

# set sampleConditions for Nicole's RNAseq ecad neg/pos experiment
sampleCondition <- c('E_minus','E_plus')
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory = directory, design=~condition)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, levels=c('E_plus','E_minus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)

# since only comparing 1 sample for each condition - multiple hypothesis testing will not work
# instead - I am sorting hits based on lowest p-value
res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# pvalue should be ranked lowest pvalue to high (most sig to least sig)

# save data 'res' to csv!
write.csv(as.data.frame(res),file='2014-7-29-DESeq2_analysis_nonGBX2group.csv')

mcols(res, use.names=T)
write.csv(as.data.frame(mcols(res,use.name=T)),file='2014-7-28-DESeq2-test-conditions-nonGBX2group.csv')

# transform raw distrbuted counts for clustering analysis
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# views clustering on individual datasets in unbiased way
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"2014-7-23-DESeq2_variance_stabilizing.png")
# axis is square root of variance over the mean for all samples
dev.off()

# clustering analysis
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-7-23-deseq2_heatmaps.png")
dev.off()


#heatmap of data
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:3000]
se#hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(assay(vsd)[select, ], col=redblue(16),
          Rowv = F, Colv = F, scale= 'non',
          dendrogram = 'none', trace = 'none', margin = c(6,6))
dev.copy(png, '2014-7-23-DESeq2_heatmap10.png')
dev.off()


# PCA analysis
# run principal component analysis on data
# good for visualizing effect of experimental covariats and batch effect
# ideal for examining primary and matching mets
# using plotPCAWithNames.R script - can change number of ntop
print(plotPCAWithSampleNames(rld, intgroup=c('condition')))
dev.copy(png, "2014-7-28-DESeq2-PCA-nongbx2-names.png")
dev.off()

# make a usable heatmap of Gbx2 dataset
# from Dave Tang's tutorial
# http://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/

# made merged_couns.txt per Dave Wheeler
library(gplots)
library("RColorBrewer")
library(preprocessCore)
data <- read.delim('sig.counts.txt',header=T,row.names='gene')
head(data)
nrow(data)
quantile(rowSums(data))
data_subset <- data[rowSums(data),]
nrow(data_subset)
heatmap(data.matrix(data_subset))
heatmap.2(data.matrix(fit2),scale='row')
colfunc <- colorRampPalette(c("blue",'black','red'))
heatmap.2(data.matrix(data), col=colfunc(100),
          Rowv = F, Colv = F, scale= 'row',key=T,keysize=1.0,density.info='none',
          dendrogram = 'none',cexCol=0.9,labRow=NA, trace = 'none', margin = c(10,5))

data2 <- log2(data)
boxplot(data2, cex.axis=0.8, las=2, main = 'original distrubtion of data',ylab='log2')
