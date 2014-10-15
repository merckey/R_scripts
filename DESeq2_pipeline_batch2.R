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
load("2014-10-7-batch2.RData")
load("2014-10-7-2D-PCA_ggplot2.RData")

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

resdata <- merge(as.data.frame(resMF), as.data.frame(assay(ddsMF,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="2014-10-13-results-DESeq2-with-normalized.csv")
write.csv(as.data.frame(assay(ddsMF,normalized=T)), file="2014-10-13-DESeq2-ddsMF-normalized.csv")

# add value 8 to each cell for log fc calc
ddsMF_8 <- (assay(ddsMF) + 8)
ddsMF_8
write.csv(as.data.frame(ddsMF_8), file="2014-10-13-DESeq2-ddsMF_8.csv")


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

# log2 transformation 
logcounts <- log2(counts(ddsMF, normalized=T) + 1)

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
px <- counts(ddsMF)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsdMF)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
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
select <- order(rowMeans(counts(vsdMF,normalized=T)),decreasing=T)[1:2000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(log(counts(ddsMF)), col=my_palette,
          scale="row", key=T, keysize=1,symkey=T,density.info="none", 
          trace="none",cexCol=0.6, labRow=F,
          main="Transcriptome of EMT populations")
dev.copy(png, "2014-9-5-DESeq2_heatmap_top1000.png")
dev.off()

# to specify individual genes from heatmap
sample <- rowMeans(counts(ddsMF))
select[1:50]
# will give you vector of numbers

sample[select[1:50]] # will give you gene name and normalized value for top 32 genes

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
distsRL <- dist(t(assay(vsdMF)))
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
# 3D PCA plot and 2D PCA
library(grDevices)
library(rgl)
library("genefilter")
library("ggplot2")

rv <- rowVars(assay(vsdMF))
select <- order(rv, decreasing=T)[seq_len(min(1000,length(rv)))] # as of 10-7-14 - 500 top variable genes
pc <- prcomp(t(assay(vsdMF)[select,]))

# 2D PCA plot with ggplot2
scores <- data.frame(sampleFiles, pc$x, condition, grouping, fit$classification, stringsAsFactors=T)
scores

scores$grouping <- as.factor(scores$grouping)
pca <- ggplot(scores, aes(x=PC1, y=PC2, col=scores$grouping))  # or factor(grouping)
(pca <- pca + geom_point(size = 5))
(pca <- pca + ggtitle("Principal Components"))
(pca <- pca + scale_color_brewer(name = "", palette="Set1"))
(pca <- pca + theme(
  plot.title = element_text(face='bold'),
  legend.justification=c(0,0),
  legend.position=c(0,0), 
  legend.key = element_rect(fill='NA'),
  legend.text = element_text(size=10, face="bold"),
  axis.text.y=element_text(colour="Black"), 
  axis.text.x=element_text(colour="Black"),
  axis.title.x=element_text(face="bold"), 
  axis.title.y=element_text(face='bold'),
  panel.grid.major.x=element_blank(), 
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x=element_blank(), 
  panel.grid.minor.y=element_blank(), 
  panel.background=element_rect(color='black',fill=NA)
  ))

save(pca, file="2014-10-10-2D-PCA_ggplot2_3group.RData")
ggsave("2014-10-10-2D-PCA_ggplot2_3group.png")

# using grouping to ID sub classes
grouping <- factor(c("Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Canonical","Canonical",
                     "Non-Canonical","Non-Canonical",
                     "Non-Canonical","Non-Canonical"))

# set condition
condition <- factor(c("Ecad-","Ecad+", "Ecad-","Ecad+", 
               "Ecad-","Ecad+", "Ecad-","Ecad+", 
               "Ecad-","Ecad+", "Ecad-","Ecad+",
               "Ecad-","Ecad+", "Ecad-","Ecad+"))
condition <- as.integer(condition)



# model-based clustering alorithm (mclust)
# selecting model - multivariat mixture = VII = diagonal, varying volume and shape
fit <- Mclust(pc$x, G=3, modelNames = "VII")
fit$classification
# 3D pca with rgl plot3d using same pc 
plot3d(pc$x,size=15,lit=T, box=F, alpha=0.75, axes=F, col=grouping)
text3d(pc$x, text=sampleFiles, adj=-0.2, font=1, cex=0.8, main="Principal Components")
grid3d(side='z',at=list(z=0))
# plot3d(pc$x, col=pca$cluster)
# save as postscript pdf file then edit in gimp
rgl.snapshot("2014-10-9-3D-pca_view2.png",'png')

# or use the pca3d package 
# I prefere this 
library("pca3d")
threedpca <- pca3d(pc, show.plane=F, col=scores$grouping, radius=2, axes.color='black',
      show.axe.titles=F)
rgl.snapshot("2014-10-10-pca3d.grouping.png",'png')
save(threedpca, file="2014-10-10-3d-PCA_3group.RData")

# optional K-means based clustering
km <- kmeans(assay(vsdMF),3) # not sure how to integrate into plot
pc$cluster <- as.integer(km$cluster)
plot3d(pc$x,size=10,lit=T, axes=F,box=F)
text3d(pc$x, text=sampleFiles, adj=-0.2, font=0.5, col=colorRampPalette(c("red",'red'))(2))
plot3d(pc$x, col=pca$cluster)

save.image("~/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_batch2/2014-10-10-batch2.RData")

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


volcanoplot <- function (x, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topleft", labelsig=TRUE, textcx=1, ...) {
  with(x, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(x, padj < sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(x, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(x, padj < sigthresh & abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(x, padj < sigthresh & abs(log2FoldChange) > lfcthresh), textxy(log2FoldChange, -log10(pvalue), cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
plot(resMF$log2FoldChange,-log10(resMF$padj),pch="+",cex=0.1,
     xlab="Log Fold-Change from non-filtered eBayes",
     ylab="-log10 P-Value from non-filtered eBayes",
     main='Volcano Plot')
dev.copy("file",png)
sessionInfo()


with(resMF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resMF, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resMF, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(resMF, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(resMF, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=condition, cex=.8))
