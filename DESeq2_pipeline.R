### bioinformatic analysis of RNA-seq data using DESeq2
### recipe is from Dave Wheeler's blog at Massey University
### http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###
### major differece from DESeq1 is that you don't need to merge individual HT-seq-count files
# will be comparing 2412, 345, 2204 Ecad +/- as group Canonical-EMT
# will be comparing 2329, 2523, 2342, 798 Ecad +/- as group non-Canonical-EMT
#
# analysis per 2014-7-23 PM on htseq count -intersections-nonempty

library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTseq_DESeq2_analysis/Counts'

# can merge individual sample fiels (i.e. control 1, control 2, etc.)
sampleFiles <- grep('PD',list.files(directory),value=T)

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

colData(ddsHTseq)$condition<-factor(colData(ddsHTseq)$condition, levels=c('E_minus','E_plus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)
# res <- results(dds)
res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# padj should be ranked lowest adj pval to high (most sig to least sig)
# save data 'res' to csv!
write.csv(as.data.frame(res),file='2014-7-23-DESeq2_pooled_withoutPD1849.csv')



attr(res, 'filterThreshold')

plot(attr(res,'filterNumRej'),type='b',ylab='number of rejections')

# visualize data wihtout independent filtering
resNoFilt <- results(dds, independentFiltering=F)
filter_table <- table(filtering=(res$padj <.1), noFiltering=(resNoFilt$padj < .1))
resNoFilt <- resNoFilt[order(resNoFilt$padj),]
head(resNoFilt)
write.table(as.data.frame(filter_table),file='filter_table')
write.csv(as.data.frame(resNoFilt),file='2014-7-23-DESeq2_results_NOFILTERING.CSV')


# plot MAplot http://en.wikipedia.org/wiki/MA_plot
plotMA(dds, ylim=c(-8,8),main = "Ecad+/- RNAseq")
dev.copy(png, "2014-7-23-Deseq2_Ecadplusminus_MAplot.png")
dev.off()

mcols(res, use.names=T)
# produces DataFrame of results of statistical tests

# opitonal annotation of res data - don't need to do this if gene names are specified in HTseq-count
#library(org.Mm.eg.db)

#keytypes(org.Mm.eg.db)

#fbids <- rownames(res)
#cols <- 'SYMBOL'
#annots <- select(org.Mm.eg.db,keys=fbids,cols=cols,keytype='ENTREZID')
#res <- cbind(entrezid=rownames(as.data.frame(res),as.data.frame(res)))
#new_res <- merge(res,annots,by.x='entrezid',by.y='ENTREZID')

#write.csv(new_res,file='DESEQ_results_experimentalname_annotated.csv')


# transform raw distrbuted counts for clustering analysis
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# scatter plot of rlog transformations between Sample conditions
head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,3:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,7:8],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,9:10],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,11:12],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,13:14],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# did not have this sample - plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

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

library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,4))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,4))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,4))

#heatmap of data
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(dds, normalized=T)[select,],col=hmcol,
          Rowv = F, Colv = F, scale='none',
          dendrogram='none',trace='none',margin=c(10,6))
dev.copy(png, "2014-7-23-DESeq2_heatmap1.png")
dev.off()
# using the rlog transformed data 'rld'
heatmap.2(assay(rld)[select, ], col=hmcol,
          Rowv=F, Colv = F, scale = 'none',
          dendrogram ='none', trace='none',margin=c(10,6))
dev.copy(png, "2014-7-23-DESeq2_heatmap2.png")
dev.off()
heatmap.2(assay(vsd)[select, ], col=hmcol,
          Rowv = F, Colv = F, scale= 'non',
          dendrogram = 'none', trace = 'none', margin = c(10,6))
dev.copy(png, '2014-7-23-DESeq2_heatmap3.png')
dev.off()

# heatmap 1 = raw counts
# heatmap 2 = regularized log transformation
# heatmap 3 = variance stabilizing transformation
# transformation shrinks the variance 

# clustering analysis
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-7-23-deseq2_heatmaps.png")
dev.off()

# PCA analysis
# run principal component analysis on data
# good for visualizing effect of experimental covariats and batch effect
# ideal for examining primary and matching mets
print(plotPCA(rld, intgroup=c('condition')))
dev.copy(png, "DESeq2_Ecad_RNAseq_PCA.png")
dev.off()


#### replacing outlier value with estimated value as predicted
#### by distrubution using "trimmed mean" approach
#### recommended if you have several replicates per treatment
#### DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file='2014-6-10-replaceoutliers-results1.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file='2014-6-10-replaceoutliers-results2.csv')

# use for dispersion plot 
plotDispEsts(dds)
dev.copy(png, "DESeq_RNAseq_ecad_minusplus_dispersionestimates.png")
dev.off()

### independent filtering to remove any tests
### that probably won't pass to reduce False discover
### by reducing total number of tests to perform

# filtering threshold
attr(res, 'filterThreshold')
plot(attr(res,'filterNumRej'), type = 'b', ylab='number of rejection')
dev.copy(png,'deseq2_filtering_threshold.png')
dev.off()

W <- mcols(dds)$WaldStatistic_condition_treated_vs_untreated
maxCooks <- apply(assays(dds)[['cooks']],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab='rank of Wald statistic',
     ylab = "maximum Cook's distance per gene",
     ylim = c(0,5), cex = 0.4, col=rgb(0,0,0,0.3))
m <- ncol(dds)
p <- 3
abline(h=qf(0.99,p,m-p))
dev.copy(png, 'deseq2_cooksdist.png')
dev.off()


# visualize pvales discarded by filtering
use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=F)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=F)
colori <- c('do not pass'="khaki", 'pass'="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = F,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

sessionInfo()
