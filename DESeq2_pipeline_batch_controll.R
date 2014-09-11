### bioinformatic analysis of RNA-seq data using DESeq2
### recipe is from Dave Wheeler's blog at Massey University
### http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###
### analysis per 2014-8-12 PM on htseq count -intersections-nonempty
### controlling for batch
### removing PD8640 and 1849 as that sample is realy weird 


library("DESeq2")
setwd("~/Desktop/RNAseq_Nicole_Ecad/HTSeq_DESeq2_analysis_batch2/")
directory <- '/Users/dballi/Desktop/RNAseq_Nicole_Ecad/HTseq_DESeq2_analysis/Counts'

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

# multifactor designs
# can analysis with more than one factor influencing the counts 
# from manual section 1.5

ddsMF <- dds

design(ddsMF) <- formula(~ Batch + condition)

# gut of DESeq2 analysis controlling for Batch
ddsMF <- DESeq(ddsMF)

plotMA(ddsMF, ylim=c(-8,8),main = "Ecad+/- RNAseq batch effect controlled")
dev.copy(png, "2014-8-12-DESeq2_MAplot_batchcontrolled_final.png")
dev.off()

rld2 <- rlogTransformation(ddsMF, blind=T)
print(plotPCAWithSampleNames(rld2, intgroup=c('condition')))
dev.copy(png, "2014-8-12-DESeq2_MAplot_PCA_batchcontrol_final.png")
dev.off()

resMF <- results(ddsMF)

resMF <- resMF[order(resMF$padj),]

head(resMF)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# padj should be ranked lowest adj pval to high (most sig to least sig)
# save data 'res' to csv!
write.csv(as.data.frame(resMF),file='2014-8-12-DESeq2_Batchcontrol_final.csv')

# clustering analysis
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld2)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsMF),paste(condition, type, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"2014-8-12-DESeq2_batchcontrol_clustering_final.png")
dev.off()


attr(resMF, 'filterThreshold')

plot(attr(resMF,'filterNumRej'),type='b',ylab='number of rejections')

# visualize data wihtout independent filtering
resNoFilt <- results(ddsMF, independentFiltering=F)
filter_table <- table(filtering=(resMF$padj <.1), noFiltering=(resNoFilt$padj < .1))
resNoFilt <- resNoFilt[order(resNoFilt$padj),]
head(resNoFilt)
write.table(as.data.frame(filter_table),file='filter_table')
write.csv(as.data.frame(resNoFilt),file='2014-8-12-DESeq2_batch_results_NOFILTERING.CSV')


# plot MAplot http://en.wikipedia.org/wiki/MA_plot
plotMA(dds, ylim=c(-8,8),main = "Ecad+/- RNAseq")
dev.copy(png, "2014-8-8-Deseq2_allsamples_MAplot.png")
dev.off()

mcols(res, use.names=T)
write.csv(as.data.frame(mcols(res,use.name=T)),file='2014-8-8-DESeq2-test-conditions-allsamples.csv')
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
vsd <- varianceStabilizingTransformation(ddsMF, blind=T)

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
dev.copy(png, "2014-8-12-DESeq2-PCA-batchcontrolled-all2.png")
dev.off()

sessionInfo()

