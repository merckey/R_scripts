### bioinformatic analysis of RNA-seq data using DESeq2
### recipe is from Dave Wheeler's blog at Massey University
### http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###
### analysis per 2014-8-12 PM on htseq count -intersections-nonempty
### controlling for batch
### removing PD8640 as that sample is realy weird 
### adding PD1849 back into analysis

library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_batch2/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/Counts/'

# can merge individual sample fiels (i.e. control 1, control 2, etc.)
sampleFiles <- grep('PD',list.files(directory),value=T)

# view sampleFiles
sampleFiles

sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch1","Batch1","Batch2","Batch2")

# set sampleConditions for Nicole's RNAseq ecad neg/pos experiment
sampleCondition <- c('E_minus','E_plus')
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition,
                          Batch = sampleBatch)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                       directory = directory, 
                                       design= ~condition) # + batch)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition, 
                                      levels=c('E_plus','E_minus'))


# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)

res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# padj should be ranked lowest adj pval to high (most sig to least sig)
# save data 'res' to csv!
write.csv(as.data.frame(res),file='2014-8-8-DESeq2_allsamples_withoutPD8640_PD1849.csv')

# multifactor designs
# can analysis with more than one factor influencing the counts 
# from manual section 1.5

ddsMF <- dds
design(ddsMF) <- formula(~ Batch + condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
plotMA(ddsMF, ylim=c(-8,8),main = "Ecad+/- RNAseq batch affect controlled")
dev.copy(png, "2014-8-12-DESeq2_MAplot_batchcontrolled_allsamples2.png")
dev.off()

rld2 <- rlogTransformation(ddsMF, blind=T)
print(plotPCAWithSampleNames(rld2, intgroup=c('condition')))

resMF <- results(ddsMF)

resMF <- resMF[order(resMF$padj),]

head(resMF)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# padj should be ranked lowest adj pval to high (most sig to least sig)
# save data 'res' to csv!
write.csv(as.data.frame(resMF),file='2014-8-12-DESeq2_allsamples_withBATCHcontrol.csv')

#order by downregulated genes
resMF.log2 <- resMF[order(resMF$log2FoldChange),]
resMF.down <- resMF.log2
head(resMF.down)
#order by upregulated genes
resMF.up <- resMF[order(-resMF$log2FoldChange),]
head(resMF.up)

resMF.up.550 <- resMF.up[1:550,]
write.csv(as.data.frame(resMF.up.550),file="2014-8-31-DESeq2_lfc_2fold_up_batch2.csv")

resMF.down.116 <- resMF.down[1:116,]
write.csv(as.data.frame(resMF.down.116),file="2014-8-31-DESeq2_lfc_2fold_down_batch2.csv")

# alternative way to order by log fold change per manual
# put padj values are all kind of weird 
resL <- results(dds, lfcThreshold=2.0, altHypothesis="less")
resL <- resL[order(resL$log2FoldChange),]
head(resL)
write.csv(as.data.frame(resL),file='2014-8-31-DESeq2_lfc_2fold_down_batch2.csv')

resH <- results(dds, lfcThreshold = 2.0, altHypothesis='greater')
resH <- resH[order(-resH$log2FoldChange),]
head(resH, n = 300)



attr(resMF, 'filterThreshold')

plot(attr(resMF,'filterNumRej'),type='b',ylab='number of rejections')

# visualize data wihtout independent filtering
resNoFilt <- results(dds, independentFiltering=F)
filter_table <- table(filtering=(res$padj <.1), noFiltering=(resNoFilt$padj < .1))
resNoFilt <- resNoFilt[order(resNoFilt$padj),]
head(resNoFilt)
write.table(as.data.frame(filter_table),file='filter_table')
write.csv(as.data.frame(resNoFilt),file='2014-8-8-DESeq2_results_NOFILTERING.CSV')


# plot MAplot http://en.wikipedia.org/wiki/MA_plot
plotMA(ddsMF, ylim=c(-8,8),main = "Ecad+/- RNAseq")
dev.copy(png, "2014-9-5-Deseq2_batch2_allsamples_MAplot.png")
dev.off()

mcols(res, use.names=T)
write.csv(as.data.frame(mcols(res,use.name=T)),file='2014-8-8-DESeq2-test-conditions-allsamples.csv')
# produces DataFrame of results of statistical tests
# transform raw distrbuted counts for clustering analysis
rldMF <- rlogTransformation(ddsMF, blind=T)
vsdMF <- varianceStabilizingTransformation(ddsMF, blind=T)

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
plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
#plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")



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
# 1000 top fold change genes
select <- order(rowMeans(counts(ddsMF,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsdMF)[select,], col=my_palette,
          scale="row", key=T, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F,
          main="Transcriptome of EMT populations")
dev.copy(png, "2014-9-5-DESeq2_heatmap_top1000.png")
dev.off()

# to specify individual genes from heatmap
sample <- rowMeans(counts(ddsMF))
select[1:10]
# will give you vector of numbers
[1]  4858 19818  2035 11958 17871 11417  6058 17649 21136
[10]  5752
sample[4858] # will give you gene name and normalized value

# 550 top fold change genes
select <- order(rowMeans(counts(ddsMF,normalized=T)),decreasing=T)[1:550]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsdMF)[select,], col=my_palette,
          scale="row", key=T, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F,
          main="Transcriptome of EMT populations")
dev.copy(png, "2014-9-5-DESeq2_heatmap_top550.png")
dev.off()

select.padj = order(res$padj,decreasing=FALSE)[1:400]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
par(cex.main=0.8)
p3 <-heatmap.2(assay(vsd)[select.padj,], col=my_palette,
          scale="row", key=T, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F, margins=c(8,6),
          main="Transcriptome of EMT populations")
dev.copy(png, "2014-9-5-DESeq2_heatmap_padj_top4002.png")
dev.off()


# clustering analysis
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-8-8-deseq2_heatmaps_withoutPD8640_PD1849.png")
dev.off()

# PCA analysis
# run principal component analysis on data
# good for visualizing effect of experimental covariats and batch effect
# ideal for examining primary and matching mets

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=200)  # can change ntop from 200-500
{
  library("genefilter")
  library("lattice")
  
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))
  
  # extract sample names
  names = colnames(x)
  
  fac = factor(apply( as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  if( nlevels(fac) >= 3 )
    colours = brewer.pal(nlevels(fac), "Dark2")
  else  
    colours = c( "red", "blue" )
  
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    panel=function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
    },
    aspect = "iso", col=colours,
    main = draw.key(key = list(
      rect = list(col = colours),
      text = list(levels(fac)),
      rep = FALSE)))
}

print(plotPCAWithSampleNames(rld, intgroup=c('condition')))
dev.copy(png, "2014-9-10--PCA-batch-names.png")
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
plotDispEsts(ddsMF)
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

