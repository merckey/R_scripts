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
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_GBX2_group/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_GBX2_group/'

sampleFiles <- grep('s.counts.txt',list.files(directory),value=T)

sampleName <- c("PD2204E_minus","PD2204E_plus",
                "PD2523E_minus","PD2523E_plus",
                "PD798E_minus","PD798E_plus")
# view sampleFiles
sampleFiles

# set sampleConditions for Nicole's RNAseq ecad neg/pos experiment
sampleCondition <- c('E_minus','E_plus')
sampleTable <- data.frame(sampleName = sampleName, 
                          fileName = sampleFiles, 
                          condition = sampleCondition)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory = directory, 
                                       design= ~condition)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, 
                                    levels=c('E_plus','E_minus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)

# since only comparing 1 sample for each condition - multiple hypothesis testing will not work
# instead - I am sorting hits based on lowest p-value
res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# pvalue should be ranked lowest pvalue to high (most sig to least sig)

# MA plot of log fold change over the mean of normalized counts 
# set alpha to 0.05 for genes with FDR < 0.05
plotMA(res, ylim=c(-12,12), alpha=0.05)
dev.copy(png,'2014-8-6-MAplot-gbx2groupmore.png')
dev.off()

# save data 'res' to csv!
write.csv(as.data.frame(res),file='2014-8-6-DESeq2_analysis_GBX2group.csv')

mcols(res, use.names=T)
write.csv(as.data.frame(mcols(res,use.name=T)),file='2014-8-6-DESeq2-test-conditions-GBX2group.csv')

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
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-7-23-deseq2_heatmaps.png")
dev.off()


#heatmap of data
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:3000]
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
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
print(plotPCA(rld, intgroup=c('condition')))
dev.copy(png, "2014-7-23-DESeq2-PCA-gbx2.png")
dev.off()


W <- res$stat
maxCooks <- apply(assays(dds)[['cooks']],1,max)
idx <- !is.na(W)
plot(rank(W[idx]),maxCooks[idx],xlab='rank of Wald statistic', 
     ylab = "maximum Cook's distance per gene", 
     ylim = c(0,5), cex=0.4, col=rgb(0,0,0,0.3))
m <- ncol(dds)
p <- 3
abline(h=qf(0.99, p, m - p ))

plot(res$baseMean+1, -log10(res$pvalue),
     log = 'x', xlab = ' mean of normalized counts',
     ylab = expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=0.4, col=rgb(0,0,0,.3))
dev.copy(png,'2014-7-29-DESeq2-Meancountsfilter-canonical-group.png')
dev.off()


# Diagnostics plot for multiple hypothesis testing
# Benjamini-Hockberg multiple hyothesis adjustment
resFilt <- res[use & !is.na(res$pvalue),]
orderinPlot <- order(resFilt$pvalue)
showinPlot <- (resFilt$pvalue[orderinPlot] <= 0.08)
alpha <- 0.05
plot(seq(along=which(showinPlot)), resFilt$pvalue[orderinPlot][showinPlot],
     pch='.',xlab= expression(rank(p[i])),ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue),col='red3',lwd=2)

plot(1-resFilt$pvalue[orderinPlot],
     (length(resFilt$pvalue)-1):0,pch='.',
     xlab=expression(1-p[i]),ylab=expression(N(p[i])))
abline(a=0, b=slope, col='red3',lwd=2)


ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file='2014-8-12-orgGBX2_replacedoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file='2014-8-12-orgGBX2_replacedoutliers.csv')

# use for dispersion plot 
plotDispEsts(dds)
dev.copy(png, "DESeq_RNAseq_ecad_minusplus_dispersionestimates.png")
dev.off()